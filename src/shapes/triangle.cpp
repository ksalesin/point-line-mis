/*
 pbrt source code is Copyright(c) 1998-2016
 Matt Pharr, Greg Humphreys, and Wenzel Jakob.

 This file is part of pbrt.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 - Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 - Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

/*
   OnColoredLine (to visualize line direction on light source), VertexColor added by Kate (12/2018)
*/

// shapes/triangle.cpp*
#include "shapes/triangle.h"
#include "texture.h"
#include "textures/constant.h"
#include "paramset.h"
#include "sampling.h"
#include "efloat.h"
#include "reflection.h"
#include "ext/rply.h"
#include <array>
#include "point-line-mis/triangle.h"
#include "point-line-mis/intersections.h"

using linesampling::Interval;
using linesampling::LineVisibilityQuery3f;
using linesampling::IntersectTriangleTriangle;

namespace pbrt {

STAT_PERCENT("Intersections/Ray-triangle Intersect() tests", nHits, nTests);
STAT_PERCENT("Intersections/Ray-triangle IntersectP() tests", nShadowHits,
		nShadowTests);
STAT_PERCENT("Intersections/Ray-triangle IntersectL() tests", nLineHits,
		nLineTests);

// Triangle Local Definitions
static void PlyErrorCallback(p_ply, const char *message) {
	Error("PLY writing error: %s", message);
}

// Triangle Method Definitions
STAT_RATIO("Scene/Triangles per triangle mesh", nTris, nMeshes);
TriangleMesh::TriangleMesh(const Transform &ObjectToWorld, int nTriangles,
		const int *vertexIndices, int nVertices, const Point3f *P,
		const Vector3f *S, const Normal3f *N, const Point2f *UV,
		const std::shared_ptr<Texture<Float>> &alphaMask,
		const std::shared_ptr<Texture<Float>> &shadowAlphaMask) :
		nTriangles(nTriangles), nVertices(nVertices), vertexIndices(
				vertexIndices, vertexIndices + 3 * nTriangles), alphaMask(
				alphaMask), shadowAlphaMask(shadowAlphaMask) {
	++nMeshes;
	nTris += nTriangles;
	triMeshBytes += sizeof(*this) + (3 * nTriangles * sizeof(int))
			+ nVertices
					* (sizeof(*P) + (N ? sizeof(*N) : 0) + (S ? sizeof(*S) : 0)
							+ (UV ? sizeof(*UV) : 0));

	// Transform mesh vertices to world space
	p.reset(new Point3f[nVertices]);
	for (int i = 0; i < nVertices; ++i)
		p[i] = ObjectToWorld(P[i]);

	// Copy _UV_, _N_, and _S_ vertex data, if present
	if (UV) {
		uv.reset(new Point2f[nVertices]);
		memcpy(uv.get(), UV, nVertices * sizeof(Point2f));
	}
	if (N) {
		n.reset(new Normal3f[nVertices]);
		for (int i = 0; i < nVertices; ++i)
			n[i] = ObjectToWorld(N[i]);
	}
	if (S) {
		s.reset(new Vector3f[nVertices]);
		for (int i = 0; i < nVertices; ++i)
			s[i] = ObjectToWorld(S[i]);
	}
}

Float SolidAngleGirard(const Point3f&p, const Point3f& p0, const Point3f& p1,
		const Point3f& p2) {
	// Project the vertices into the unit sphere around p.

	const Vector3f& v0 = Normalize(p0 - p);
	const Vector3f& v1 = Normalize(p1 - p);
	const Vector3f& v2 = Normalize(p2 - p);

	// http://math.stackexchange.com/questions/9819/area-of-a-spherical-triangle
	// Girard's theorem: surface area of a spherical triangle on a unit
	// sphere is the 'excess angle' alpha+beta+gamma-pi, where
	// alpha/beta/gamma are the interior angles at the vertices.
	//
	// Given three vertices on the sphere, a, b, c, then we can compute,
	// for example, the angle c->a->b by
	//
	// cos theta =  Dot(Cross(c, a), Cross(b, a)) /
	//              (Length(Cross(c, a)) * Length(Cross(b, a))).
	//
	Vector3f cross01 = Cross(v0, v1);
	Vector3f cross12 = Cross(v1, v2);
	Vector3f cross20 = Cross(v2, v0);

	// Some of these vectors may be degenerate. In this case, we don't want
	// to normalize them so that we don't hit an assert. This is fine,
	// since the corresponding dot products below will be zero.
	if (cross01.LengthSquared() > 0)
		cross01 = Normalize(cross01);
	if (cross12.LengthSquared() > 0)
		cross12 = Normalize(cross12);
	if (cross20.LengthSquared() > 0)
		cross20 = Normalize(cross20);

	// We only need to do three cross products to evaluate the angles at
	// all three vertices, though, since we can take advantage of the fact
	// that Cross(a, b) = -Cross(b, a).
	return std::abs(
			std::acos(Clamp(Dot(cross01, -cross12), -1, 1))
					+ std::acos(Clamp(Dot(cross12, -cross20), -1, 1))
					+ std::acos(Clamp(Dot(cross20, -cross01), -1, 1)) - Pi);
}
std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
		const Transform *ObjectToWorld, const Transform *WorldToObject,
		bool reverseOrientation, int nTriangles, const int *vertexIndices,
		int nVertices, const Point3f *p, const Vector3f *s, const Normal3f *n,
		const Point2f *uv, const std::shared_ptr<Texture<Float>> &alphaMask,
		const std::shared_ptr<Texture<Float>> &shadowAlphaMask) {
	std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>(
			*ObjectToWorld, nTriangles, vertexIndices, nVertices, p, s, n, uv,
			alphaMask, shadowAlphaMask);
	std::vector<std::shared_ptr<Shape>> tris;
	tris.reserve(nTriangles);
	for (int i = 0; i < nTriangles; ++i) {
		std::shared_ptr<Triangle> triangle = std::make_shared<Triangle>(
				ObjectToWorld, WorldToObject, reverseOrientation, mesh, i);
		if (triangle->Area() < 1e-9)
			Warning("found a degenerate triangle with area %.6f\n",
					triangle->Area());
		else
			tris.push_back(triangle);
	}
	return tris;
}

bool WritePlyFile(const std::string &filename, int nTriangles,
		const int *vertexIndices, int nVertices, const Point3f *P,
		const Vector3f *S, const Normal3f *N, const Point2f *UV) {
	p_ply plyFile = ply_create(filename.c_str(), PLY_DEFAULT, PlyErrorCallback,
			0, nullptr);
	if (plyFile != nullptr) {
		ply_add_element(plyFile, "vertex", nVertices);
		ply_add_scalar_property(plyFile, "x", PLY_FLOAT);
		ply_add_scalar_property(plyFile, "y", PLY_FLOAT);
		ply_add_scalar_property(plyFile, "z", PLY_FLOAT);
		if (N != nullptr) {
			ply_add_scalar_property(plyFile, "nx", PLY_FLOAT);
			ply_add_scalar_property(plyFile, "ny", PLY_FLOAT);
			ply_add_scalar_property(plyFile, "nz", PLY_FLOAT);
		}
		if (UV != nullptr) {
			ply_add_scalar_property(plyFile, "s", PLY_FLOAT);
			ply_add_scalar_property(plyFile, "t", PLY_FLOAT);
		}
		if (S != nullptr)
			Warning("PLY mesh in \"%s\" will be missing tangent vectors \"S\".",
					filename.c_str());

		ply_add_element(plyFile, "face", nTriangles);
		ply_add_list_property(plyFile, "vertex_indices", PLY_UINT8, PLY_INT);
		ply_write_header(plyFile);

		for (int i = 0; i < nVertices; ++i) {
			ply_write(plyFile, P[i].x);
			ply_write(plyFile, P[i].y);
			ply_write(plyFile, P[i].z);
			if (N) {
				ply_write(plyFile, N[i].x);
				ply_write(plyFile, N[i].y);
				ply_write(plyFile, N[i].z);
			}
			if (UV) {
				ply_write(plyFile, UV[i].x);
				ply_write(plyFile, UV[i].y);
			}
		}

		for (int i = 0; i < nTriangles; ++i) {
			ply_write(plyFile, 3);
			ply_write(plyFile, vertexIndices[3 * i]);
			ply_write(plyFile, vertexIndices[3 * i + 1]);
			ply_write(plyFile, vertexIndices[3 * i + 2]);
		}
		ply_close(plyFile);
		return true;
	}
	return false;
}

Bounds3f Triangle::ObjectBound() const {
	// Get triangle vertices in _p0_, _p1_, and _p2_
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];
	return Union(Bounds3f((*WorldToObject)(p0), (*WorldToObject)(p1)),
			(*WorldToObject)(p2));
}

Bounds3f Triangle::WorldBound() const {
	// Get triangle vertices in _p0_, _p1_, and _p2_
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];
	return Union(Bounds3f(p0, p1), p2);
}

bool Triangle::Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
		bool testAlphaTexture) const {
	ProfilePhase p(Prof::TriIntersect);
	++nTests;
	// Get triangle vertices in _p0_, _p1_, and _p2_
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];

	// Perform ray--triangle intersection test

	// Transform triangle vertices to ray coordinate space

	// Translate vertices based on ray origin
	Point3f p0t = p0 - Vector3f(ray.o);
	Point3f p1t = p1 - Vector3f(ray.o);
	Point3f p2t = p2 - Vector3f(ray.o);

	// Permute components of triangle vertices and ray direction
	int kz = MaxDimension(Abs(ray.d));
	int kx = kz + 1;
	if (kx == 3)
		kx = 0;
	int ky = kx + 1;
	if (ky == 3)
		ky = 0;
	Vector3f d = Permute(ray.d, kx, ky, kz);
	p0t = Permute(p0t, kx, ky, kz);
	p1t = Permute(p1t, kx, ky, kz);
	p2t = Permute(p2t, kx, ky, kz);

	// Apply shear transformation to translated vertex positions
	Float Sx = -d.x / d.z;
	Float Sy = -d.y / d.z;
	Float Sz = 1.f / d.z;
	p0t.x += Sx * p0t.z;
	p0t.y += Sy * p0t.z;
	p1t.x += Sx * p1t.z;
	p1t.y += Sy * p1t.z;
	p2t.x += Sx * p2t.z;
	p2t.y += Sy * p2t.z;

	// Compute edge function coefficients _e0_, _e1_, and _e2_
	Float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
	Float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
	Float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

	// Fall back to double precision test at triangle edges
	if (sizeof(Float) == sizeof(float)
			&& (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
		double p2txp1ty = (double) p2t.x * (double) p1t.y;
		double p2typ1tx = (double) p2t.y * (double) p1t.x;
		e0 = (float) (p2typ1tx - p2txp1ty);
		double p0txp2ty = (double) p0t.x * (double) p2t.y;
		double p0typ2tx = (double) p0t.y * (double) p2t.x;
		e1 = (float) (p0typ2tx - p0txp2ty);
		double p1txp0ty = (double) p1t.x * (double) p0t.y;
		double p1typ0tx = (double) p1t.y * (double) p0t.x;
		e2 = (float) (p1typ0tx - p1txp0ty);
	}

	// Perform triangle edge and determinant tests
	if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
		return false;
	Float det = e0 + e1 + e2;
	if (det == 0)
		return false;

	// Compute scaled hit distance to triangle and test against ray $t$ range
	p0t.z *= Sz;
	p1t.z *= Sz;
	p2t.z *= Sz;
	Float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
	if (det < 0 && (tScaled >= 0 || tScaled < ray.tMax * det))
		return false;
	else if (det > 0 && (tScaled <= 0 || tScaled > ray.tMax * det))
		return false;

	// Compute barycentric coordinates and $t$ value for triangle intersection
	Float invDet = 1 / det;
	Float b0 = e0 * invDet;
	Float b1 = e1 * invDet;
	Float b2 = e2 * invDet;
	Float t = tScaled * invDet;

	// Ensure that computed triangle $t$ is conservatively greater than zero

	// Compute $\delta_z$ term for triangle $t$ error bounds
	Float maxZt = MaxComponent(Abs(Vector3f(p0t.z, p1t.z, p2t.z)));
	Float deltaZ = gamma(3) * maxZt;

	// Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
	Float maxXt = MaxComponent(Abs(Vector3f(p0t.x, p1t.x, p2t.x)));
	Float maxYt = MaxComponent(Abs(Vector3f(p0t.y, p1t.y, p2t.y)));
	Float deltaX = gamma(5) * (maxXt + maxZt);
	Float deltaY = gamma(5) * (maxYt + maxZt);

	// Compute $\delta_e$ term for triangle $t$ error bounds
	Float deltaE = 2
			* (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

	// Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
	Float maxE = MaxComponent(Abs(Vector3f(e0, e1, e2)));
	Float deltaT = 3
			* (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE)
			* std::abs(invDet);
	if (t <= deltaT)
		return false;

	// Compute triangle partial derivatives
	Vector3f dpdu, dpdv;
	Point2f uv[3];
	GetUVs(uv);

	// Compute deltas for triangle partial derivatives
	Vector2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
	Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
	Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
	bool degenerateUV = std::abs(determinant) < 1e-8;
	if (!degenerateUV) {
		Float invdet = 1 / determinant;
		dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
		dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
	}
	if (degenerateUV || Cross(dpdu, dpdv).LengthSquared() == 0)
		// Handle zero determinant for triangle partial derivative matrix
		CoordinateSystem(Normalize(Cross(p2 - p0, p1 - p0)), &dpdu, &dpdv);

	// Compute error bounds for triangle intersection
	Float xAbsSum = (std::abs(b0 * p0.x) + std::abs(b1 * p1.x)
			+ std::abs(b2 * p2.x));
	Float yAbsSum = (std::abs(b0 * p0.y) + std::abs(b1 * p1.y)
			+ std::abs(b2 * p2.y));
	Float zAbsSum = (std::abs(b0 * p0.z) + std::abs(b1 * p1.z)
			+ std::abs(b2 * p2.z));
	Vector3f pError = gamma(7) * Vector3f(xAbsSum, yAbsSum, zAbsSum);

	// Interpolate $(u,v)$ parametric coordinates and hit point
	Point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
	Point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

	// Test intersection against alpha texture, if present
	if (testAlphaTexture && mesh->alphaMask) {
		SurfaceInteraction isectLocal(pHit, Vector3f(0, 0, 0), uvHit, -ray.d,
				dpdu, dpdv, Normal3f(0, 0, 0), Normal3f(0, 0, 0), ray.time,
				this);
		if (mesh->alphaMask->Evaluate(isectLocal) == 0)
			return false;
	}

	// Fill in _SurfaceInteraction_ from triangle hit
	*isect = SurfaceInteraction(pHit, pError, uvHit, -ray.d, dpdu, dpdv,
			Normal3f(0, 0, 0), Normal3f(0, 0, 0), ray.time, this);

	// Override surface normal in _isect_ for triangle
	isect->n = isect->shading.n = Normal3f(Normalize(Cross(dp02, dp12)));
	if (mesh->n || mesh->s) {
		// Initialize _Triangle_ shading geometry

		// Compute shading normal _ns_ for triangle
		Normal3f ns;
		if (mesh->n) {
			ns = (b0 * mesh->n[v[0]] + b1 * mesh->n[v[1]] + b2 * mesh->n[v[2]]);
			if (ns.LengthSquared() > 0)
				ns = Normalize(ns);
			else
				ns = isect->n;
		} else
			ns = isect->n;

		// Compute shading tangent _ss_ for triangle
		Vector3f ss;
		if (mesh->s) {
			ss = (b0 * mesh->s[v[0]] + b1 * mesh->s[v[1]] + b2 * mesh->s[v[2]]);
			if (ss.LengthSquared() > 0)
				ss = Normalize(ss);
			else
				ss = Normalize(isect->dpdu);
		} else
			ss = Normalize(isect->dpdu);

		// Compute shading bitangent _ts_ for triangle and adjust _ss_
		Vector3f ts = Cross(ss, ns);
		if (ts.LengthSquared() > 0.f) {
			ts = Normalize(ts);
			ss = Cross(ts, ns);
		} else
			CoordinateSystem((Vector3f) ns, &ss, &ts);

		// Compute $\dndu$ and $\dndv$ for triangle shading geometry
		Normal3f dndu, dndv;
		if (mesh->n) {
			// Compute deltas for triangle partial derivatives of normal
			Vector2f duv02 = uv[0] - uv[2];
			Vector2f duv12 = uv[1] - uv[2];
			Normal3f dn1 = mesh->n[v[0]] - mesh->n[v[2]];
			Normal3f dn2 = mesh->n[v[1]] - mesh->n[v[2]];
			Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
			bool degenerateUV = std::abs(determinant) < 1e-8;
			if (degenerateUV)
				dndu = dndv = Normal3f(0, 0, 0);
			else {
				Float invDet = 1 / determinant;
				dndu = (duv12[1] * dn1 - duv02[1] * dn2) * invDet;
				dndv = (-duv12[0] * dn1 + duv02[0] * dn2) * invDet;
			}
		} else
			dndu = dndv = Normal3f(0, 0, 0);
		isect->SetShadingGeometry(ss, ts, dndu, dndv, true);
	}

	// Ensure correct orientation of the normal (this breaks refraction!)
	if (Dot(ray.d, isect->n) > 0)
		isect->n = -isect->n;
	if (Dot(isect->n, isect->shading.n) < 0)
		isect->shading.n = -isect->shading.n;

//	// Ensure correct orientation of the geometric normal
//	if (mesh->n)
//		isect->n = Faceforward(isect->n, isect->shading.n);
//	else if (reverseOrientation ^ transformSwapsHandedness)
//		isect->n = isect->shading.n = -isect->n;

	*tHit = t;
	++nHits;
	return true;
}

bool Triangle::IntersectP(const Ray &ray, bool testAlphaTexture) const {
	ProfilePhase p(Prof::TriIntersectP);
	++nShadowTests;
	++ray.triangleIsectPTests;

	// Get triangle vertices in _p0_, _p1_, and _p2_
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];

	// Perform ray--triangle intersection test

	// Transform triangle vertices to ray coordinate space

	// Translate vertices based on ray origin
	Point3f p0t = p0 - Vector3f(ray.o);
	Point3f p1t = p1 - Vector3f(ray.o);
	Point3f p2t = p2 - Vector3f(ray.o);

	// Permute components of triangle vertices and ray direction
	int kz = MaxDimension(Abs(ray.d));
	int kx = kz + 1;
	if (kx == 3)
		kx = 0;
	int ky = kx + 1;
	if (ky == 3)
		ky = 0;
	Vector3f d = Permute(ray.d, kx, ky, kz);
	p0t = Permute(p0t, kx, ky, kz);
	p1t = Permute(p1t, kx, ky, kz);
	p2t = Permute(p2t, kx, ky, kz);

	// Apply shear transformation to translated vertex positions
	Float Sx = -d.x / d.z;
	Float Sy = -d.y / d.z;
	Float Sz = 1.f / d.z;
	p0t.x += Sx * p0t.z;
	p0t.y += Sy * p0t.z;
	p1t.x += Sx * p1t.z;
	p1t.y += Sy * p1t.z;
	p2t.x += Sx * p2t.z;
	p2t.y += Sy * p2t.z;

	// Compute edge function coefficients _e0_, _e1_, and _e2_
	Float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
	Float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
	Float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

	// Fall back to double precision test at triangle edges
	if (sizeof(Float) == sizeof(float)
			&& (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
		double p2txp1ty = (double) p2t.x * (double) p1t.y;
		double p2typ1tx = (double) p2t.y * (double) p1t.x;
		e0 = (float) (p2typ1tx - p2txp1ty);
		double p0txp2ty = (double) p0t.x * (double) p2t.y;
		double p0typ2tx = (double) p0t.y * (double) p2t.x;
		e1 = (float) (p0typ2tx - p0txp2ty);
		double p1txp0ty = (double) p1t.x * (double) p0t.y;
		double p1typ0tx = (double) p1t.y * (double) p0t.x;
		e2 = (float) (p1typ0tx - p1txp0ty);
	}

	// Perform triangle edge and determinant tests
	if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
		return false;
	Float det = e0 + e1 + e2;
	if (det == 0)
		return false;

	// Compute scaled hit distance to triangle and test against ray $t$ range
	p0t.z *= Sz;
	p1t.z *= Sz;
	p2t.z *= Sz;
	Float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
	if (det < 0 && (tScaled >= 0 || tScaled < ray.tMax * det))
		return false;
	else if (det > 0 && (tScaled <= 0 || tScaled > ray.tMax * det))
		return false;

	// Compute barycentric coordinates and $t$ value for triangle intersection
	Float invDet = 1 / det;
	Float b0 = e0 * invDet;
	Float b1 = e1 * invDet;
	Float b2 = e2 * invDet;
	Float t = tScaled * invDet;

	// Ensure that computed triangle $t$ is conservatively greater than zero

	// Compute $\delta_z$ term for triangle $t$ error bounds
	Float maxZt = MaxComponent(Abs(Vector3f(p0t.z, p1t.z, p2t.z)));
	Float deltaZ = gamma(3) * maxZt;

	// Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
	Float maxXt = MaxComponent(Abs(Vector3f(p0t.x, p1t.x, p2t.x)));
	Float maxYt = MaxComponent(Abs(Vector3f(p0t.y, p1t.y, p2t.y)));
	Float deltaX = gamma(5) * (maxXt + maxZt);
	Float deltaY = gamma(5) * (maxYt + maxZt);

	// Compute $\delta_e$ term for triangle $t$ error bounds
	Float deltaE = 2
			* (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

	// Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
	Float maxE = MaxComponent(Abs(Vector3f(e0, e1, e2)));
	Float deltaT = 3
			* (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE)
			* std::abs(invDet);
	if (t <= deltaT)
		return false;

	// Test shadow ray intersection against alpha texture, if present
	if (testAlphaTexture && (mesh->alphaMask || mesh->shadowAlphaMask)) {
		// Compute triangle partial derivatives
		Vector3f dpdu, dpdv;
		Point2f uv[3];
		GetUVs(uv);

		// Compute deltas for triangle partial derivatives
		Vector2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
		Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
		Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
		bool degenerateUV = std::abs(determinant) < 1e-8;
		if (!degenerateUV) {
			Float invdet = 1 / determinant;
			dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
			dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
		}
		if (degenerateUV || Cross(dpdu, dpdv).LengthSquared() == 0)
			// Handle zero determinant for triangle partial derivative matrix
			CoordinateSystem(Normalize(Cross(p2 - p0, p1 - p0)), &dpdu, &dpdv);

		// Interpolate $(u,v)$ parametric coordinates and hit point
		Point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
		Point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];
		SurfaceInteraction isectLocal(pHit, Vector3f(0, 0, 0), uvHit, -ray.d,
				dpdu, dpdv, Normal3f(0, 0, 0), Normal3f(0, 0, 0), ray.time,
				this);
		if (mesh->alphaMask && mesh->alphaMask->Evaluate(isectLocal) == 0)
			return false;
		if (mesh->shadowAlphaMask
				&& mesh->shadowAlphaMask->Evaluate(isectLocal) == 0)
			return false;
	}
	++nShadowHits;
	return true;
}

bool Triangle::IntersectL(const LineVisibilityQuery3f &query,
		const Normal3f &normal, Float time,
		linesampling::IntervalCollection &intersection) const {
	ProfilePhase profile(Prof::TriIntersectL);
	++nLineTests;

	// Get triangle vertices in _p0_, _p1_, and _p2_
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];

	linesampling::Line3f isect;
	linesampling::Triangle3f triangle(p0, p1, p2);

	if (!linesampling::IntersectTriangleTriangle(query, triangle, isect))
		return false;

	const Interval& isectInterval = FindShadowedInterval(query, isect);

	if (intersection.subtract(isectInterval)) {
		++nLineHits;
		return true;
	} else
		return false;
}

Float Triangle::Area() const {
	// Get triangle vertices in _p0_, _p1_, and _p2_
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];
	return 0.5 * Cross(p1 - p0, p2 - p0).Length();
}

Interaction Triangle::Sample(const Point2f &u, Float *pdf) const {
	Point2f b = UniformSampleTriangle(u);
	// Get triangle vertices in _p0_, _p1_, and _p2_
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];
	Interaction it;
	it.p = b[0] * p0 + b[1] * p1 + (1 - b[0] - b[1]) * p2;
	// Compute surface normal for sampled point on triangle
	it.n = Normalize(Normal3f(Cross(p1 - p0, p2 - p0)));
	// Ensure correct orientation of the geometric normal; follow the same
	// approach as was used in Triangle::Intersect().
	if (mesh->n) {
		Normal3f ns(
				b[0] * mesh->n[v[0]] + b[1] * mesh->n[v[1]]
						+ (1 - b[0] - b[1]) * mesh->n[v[2]]);
		it.n = Faceforward(it.n, ns);
	} else if (reverseOrientation ^ transformSwapsHandedness)
		it.n *= -1;

	// Compute error bounds for sampled point on triangle
	Point3f pAbsSum = Abs(b[0] * p0) + Abs(b[1] * p1)
			+ Abs((1 - b[0] - b[1]) * p2);
	it.pError = gamma(6) * Vector3f(pAbsSum.x, pAbsSum.y, pAbsSum.z);
	*pdf = 1 / Area();
	return it;
}

Float Triangle::SolidAngle(const Point3f &p, int nSamples) const {
	// Project the vertices into the unit sphere around p.
	return SolidAngleGirard(p, mesh->p[v[0]], mesh->p[v[1]], mesh->p[v[2]]);
}

bool Triangle::CanEvaluateAnalytically(const SurfaceInteraction& isect,
		Spectrum * spectrum) const {
	// Quickly terminate when analytical integration is impossible to avoid clipping
	if (!isect.bsdf->CanEvaluateAnalytically())
		return false;

	if (spectrum) {
		const niels::Triangle t(mesh->p[v[0]], mesh->p[v[1]], mesh->p[v[2]]);

		// Clip against the surface normal
		std::vector<niels::Triangle> geometryNormalClip;
		t.Clip(isect.p, isect.n, geometryNormalClip);

		Spectrum sum(0.f);
		for (const niels::Triangle& triangle : geometryNormalClip) {
			Spectrum analytical;
			if (!isect.bsdf->CanEvaluateAnalytically(isect, triangle,
					&analytical))
				return false;
			sum += analytical;
		}

		*spectrum = sum;
	}

	return true;
}

bool Triangle::CanIlluminate(const Interaction& it) const {
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];

	const Vector3f& u = p1 - p0;
	const Vector3f& v = p2 - p0;
	const Float dot = Dot(p0 - it.p, Cross(u, v));

	if (reverseOrientation ^ transformSwapsHandedness) {
		if (dot <= 0)
			return false;
	} else {
		if (dot >= 0)
			return false;
	}

	const Float d0 = Dot(p0 - it.p, it.n);
	const Float d1 = Dot(p1 - it.p, it.n);
	const Float d2 = Dot(p2 - it.p, it.n);

	if (d0 <= 0 && d1 <= 0 && d2 <= 0)
		return false;

	return true;
}

void MinMidMax(Float x, Float y, Float z, int &i_min, int &i_mid, int &i_max) {
	if (x < y && x < z) {
		i_min = 0;

		if (y < z) {
			i_mid = 1;
			i_max = 2;
		} else {
			i_mid = 2;
			i_max = 1;
		}
	} else if (y < z) {
		i_min = 1;

		if (x < z) {
			i_mid = 0;
			i_max = 2;
		} else {
			i_mid = 2;
			i_max = 0;
		}
	} else {
		i_min = 2;

		if (x < y) {
			i_mid = 0;
			i_max = 1;
		} else {
			i_mid = 1;
			i_max = 0;
		}
	}
}

LineInteraction Triangle::LineSample(const Point2f &u,
		bool importanceSample) const {
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];

	const Point3f p[3] = { p0, p1, p2 };

	const Vector3f& e1 = p1 - p0;
	const Vector3f& e2 = p2 - p0;
	const Vector3f& cross = Cross(e1, e2);
	const Float area = cross.Length() * 0.5;

	const Vector3f& n1 = Normalize(e1);
	const Vector3f& n2 = Normalize(e2);
	const Normal3f& n = Normal3f(Normalize(cross));
	const Normal3f& normal =
			reverseOrientation ^ transformSwapsHandedness ? -n : n;
	// Create rotated orthonormal system in triangle plane
	const Vector3f& w = Normalize(Cross(cross, e1));

	const Float angle = Pi * u.x;
	const Float x1 = std::cos(angle);
	const Float x2 = std::sin(angle);

	const Vector3f& x = Normalize(x1 * n1 + x2 * w);
	const Vector3f& y = Normalize(Cross(x, normal));

	// perform sanity checks
	if (AbsDot(x, normal) > EPSILON)
		Error("x not in plane of triangle!");
	if (AbsDot(y, normal) > EPSILON)
		Error("y not in plane of triangle!");
	if (AbsDot(x, y) > EPSILON)
		Error("x and y not orthonormal!");
	if (std::abs(x.Length() - 1.0) > EPSILON)
		Error("x is not unit length!");
	if (std::abs(y.Length() - 1.0) > EPSILON)
		Error("y is not unit length!");
	if (std::abs(acos(Dot(x, n1)) - angle) > 1e-5)
		Error("Triangle::LineSample incorrect angle! %.16f != %.16f", angle,
				acos(Dot(x, n1)));

	// project vertices on x/y with origin p0
	const Point2f proj_p0(0.0, 0.0);
	const Point2f proj_p1(Dot(e1, x), Dot(e1, y));
	const Point2f proj_p2(Dot(e2, x), Dot(e2, y));
	Point2f proj_p[] = { proj_p0, proj_p1, proj_p2 };

	int pMinI, pMidI, pMaxI;

	MinMidMax(proj_p0.x, proj_p1.x, proj_p2.x, pMinI, pMidI, pMaxI);

	const Point2f& pMin = proj_p[pMinI];
	const Point2f& pMid = proj_p[pMidI];
	const Point2f& pMax = proj_p[pMaxI];

	Float xOffSet, y1, y2, pdf;

	if (importanceSample) {
		Float invRange = 1.0 / (pMax.x - pMin.x);
		Float a_l = (pMid.x - pMin.x) * invRange;
		Float a_r = (pMax.x - pMid.x) * invRange;

		Float uu;
		if (u.y < a_l) {
			uu = std::sqrt(u.y * a_l);
		} else {
			uu = 1.0 - std::sqrt(a_r * (1.0 - u.y));
		}

		xOffSet = pMin.x * (1.0 - uu) + pMax.x * uu;
		y1 = pMin.y * (1.0 - uu) + pMax.y * uu;

		if (xOffSet <= pMid.x) {
			Float s = (xOffSet - pMin.x) / (pMid.x - pMin.x);
			y2 = pMin.y * (1.0 - s) + s * pMid.y;
		} else {
			Float s = (xOffSet - pMid.x) / (pMax.x - pMid.x);
			y2 = pMid.y * (1.0 - s) + s * pMax.y;
		}
		// PDF proportional to length of line, integrates to 1 over whole triangle
		pdf = std::abs(y1 - y2) / area;
	} else {
		// Linearly interpolate along [a,b] using offset u as parameter
		xOffSet = pMin.x * (1.0 - u.y) + pMax.x * u.y;
		y1 = pMin.y * (1.0 - u.y) + pMax.y * u.y;

		if (xOffSet <= pMid.x) {
			Float s = (xOffSet - pMin.x) / (pMid.x - pMin.x);
			y2 = pMin.y * (1.0 - s) + s * pMid.y;
		} else {
			Float s = (xOffSet - pMid.x) / (pMax.x - pMid.x);
			y2 = pMid.y * (1.0 - s) + s * pMax.y;
		}
		pdf = 1.0 / (pMax.x - pMin.x);
	}

	LineInteraction result;
	// Calculate final line segment endpoints in triangle space
	const Point3f& o = p0 + x * xOffSet;
	result.line.p1 = o + y1 * y;
	result.line.p2 = o + y2 * y;
	result.pdf = pdf;
	result.length = std::abs(y1 - y2);
	result.n = normal;

	return result;
}

bool Triangle::OnColoredLine(const Interaction &intr, float direction) const {
	// Below determination of bases is copied exactly from Triangle::LineSample
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];

	const Point3f p[3] = { p0, p1, p2 };

	const Vector3f& e1 = p1 - p0;
	const Vector3f& e2 = p2 - p0;
	const Vector3f& cross = Cross(e1, e2);
	const Float area = cross.Length() * 0.5;

	const Vector3f& n1 = Normalize(e1);
	const Vector3f& n2 = Normalize(e2);
	const Normal3f& n = Normal3f(Normalize(cross));
	const Normal3f& normal =
			reverseOrientation ^ transformSwapsHandedness ? -n : n;
	// Create rotated orthonormal system in triangle plane
	const Vector3f& w = Normalize(Cross(cross, e1));

	const Float angle = Pi * direction;
	const Float x1 = std::cos(angle);
	const Float x2 = std::sin(angle);

	const Vector3f& x = Normalize(x1 * n1 + x2 * w);
	const Vector3f& y = Normalize(Cross(x, normal));

	// project vertices on x/y with origin p0
	const Point2f proj_p0(0.0, 0.0);
	const Point2f proj_p1(Dot(e1, x), Dot(e1, y));
	const Point2f proj_p2(Dot(e2, x), Dot(e2, y));
	Point2f proj_p[] = { proj_p0, proj_p1, proj_p2 };

	// Since triangle vertices are already stored in world space, we don't need to worry about converting the intersection point from world -> triangle space
	const Vector3f& e3 = intr.p - p0;
	float itsOffset = Dot(e3, x);

	int pMinI, pMidI, pMaxI;

	MinMidMax(proj_p0.x, proj_p1.x, proj_p2.x, pMinI, pMidI, pMaxI);

	const Point2f& pMin = proj_p[pMinI];
	const Point2f& pMax = proj_p[pMaxI];

	// Below this line is all Kate
	// Usually an offset u in [a,b] is given for finding the line in the triangle
	// Here we define a set of lines of width .06 between a and b

	int numLines = 5;
	float a = pMin.x;
	float b = pMax.x;
	float step = (b - a) / numLines;
	float lineWidth = (b - a) * 0.03;

	for (int i = 0; i < numLines; i++) {
		float offset = i * step + a;

		if ((std::abs(itsOffset - offset)) < lineWidth) {
			return true;
		}
	}

	return false;
}

Float Triangle::PdfLine(const Interaction &ref, const Vector3f &wi, Point3f &pLight,
        LineInteraction &lineInteraction, const Float& lineAngle, bool importanceSample) const {
    
    // Get point on quad light in wi direction
    Ray ray = ref.SpawnRay(wi);
    Float tHit;
    SurfaceInteraction isectLight;
    if (!Intersect(ray, &tHit, &isectLight, false))
        return 0.f;
    
    pLight = isectLight.p;

    // Get vertices of triangle, create a triangle-aligned orthonormal basis
    const Point3f &p0 = mesh->p[v[0]];
    const Point3f &p1 = mesh->p[v[1]];
    const Point3f &p2 = mesh->p[v[2]];

    const Point3f p[3] = { p0, p1, p2 };

    const Vector3f& e1 = p1 - p0;
    const Vector3f& e2 = p2 - p0;
    const Vector3f& cross = Cross(e1, e2);
    const Float area = cross.Length() * 0.5;

    const Vector3f& n1 = Normalize(e1);
    const Vector3f& n2 = Normalize(e2);
    const Normal3f& n = Normal3f(Normalize(cross));
    const Normal3f& normal =
            reverseOrientation ^ transformSwapsHandedness ? -n : n;
    const Vector3f& w = Normalize(Cross(cross, e1));

    // Get line direction angle (offset from n1)
    const Float angle = Pi * lineAngle;
    const Float x1 = std::cos(angle);
    const Float x2 = std::sin(angle);

    // Create an angle-aligned orthonormal basis (x, y, normal)
    const Vector3f& x = Normalize(x1 * n1 + x2 * w);
    const Vector3f& y = Normalize(Cross(x, normal));

    // perform sanity checks
    if (AbsDot(x, normal) > EPSILON)
        Error("x not in plane of triangle!");
    if (AbsDot(y, normal) > EPSILON)
        Error("y not in plane of triangle!");
    if (AbsDot(x, y) > EPSILON)
        Error("x and y not orthonormal!");
    if (std::abs(x.Length() - 1.0) > EPSILON)
        Error("x is not unit length!");
    if (std::abs(y.Length() - 1.0) > EPSILON)
        Error("y is not unit length!");
    if (std::abs(acos(Dot(x, n1)) - angle) > 1e-5)
        Error("Triangle::LineSample incorrect angle! %.16f != %.16f", angle,
                acos(Dot(x, n1)));

    // project vertices on x/y with origin p0
    const Point2f proj_p0(0.0, 0.0);
    const Point2f proj_p1(Dot(e1, x), Dot(e1, y));
    const Point2f proj_p2(Dot(e2, x), Dot(e2, y));

    const Vector3f vp = pLight - p0;
    const Point2f proj_pLight(Dot(vp, x), Dot(vp, y));

    std::vector<Point2f> proj_p{proj_p0, proj_p1, proj_p2};

    // Put vertices in order along x-axis
    std::sort(proj_p.begin(), proj_p.end(), [](const Point2f& lhs, const Point2f& rhs)
    {
        return lhs.x < rhs.x;
    });

    const Point2f pMin = proj_p[0];
    const Point2f pMid = proj_p[1];
    const Point2f pMax = proj_p[2];

    // Get line offset parameter
    double u = (proj_pLight.x - pMin.x) / (pMax.x - pMin.x);

    // Find line segment given random offset u.y along x-axis
    Float lx, ly1, ly2, pdf;

    lx = pMin.x * (1.0 - u) + pMax.x * u;
    ly1 = pMin.y * (1.0 - u) + pMax.y * u;

    if (lx <= pMid.x) {
        Float s = (lx - pMin.x) / (pMid.x - pMin.x);
        ly2 = pMin.y * (1.0 - s) + s * pMid.y;
    } else {
        Float s = (lx - pMid.x) / (pMax.x - pMid.x);
        ly2 = pMid.y * (1.0 - s) + s * pMax.y;
    }

    if (importanceSample) {
        // pdf proportional to length of line
        pdf = std::abs(ly1 - ly2) / area;
    } else {
        // linearly interpolate uniformly along offset axis
        pdf = 1.0 / (pMax.x - pMin.x);
    }

    // Calculate final line segment endpoints in triangle space
    const Point3f& lo = p0 + x * lx;
    lineInteraction.line.p1 = lo + ly1 * y;
    lineInteraction.line.p2 = lo + ly2 * y;
    lineInteraction.pdf = pdf;
    lineInteraction.length = std::abs(ly1 - ly2);
    lineInteraction.n = normal;

    return pdf;
}

std::vector<std::shared_ptr<Shape>> CreateTriangleMeshShape(
		const Transform *o2w, const Transform *w2o, bool reverseOrientation,
		const ParamSet &params,
		std::map<std::string, std::shared_ptr<Texture<Float>>>*floatTextures) {
			int nvi, npi, nuvi, nsi, nni;
			const int *vi = params.FindInt("indices", &nvi);
			const Point3f *P = params.FindPoint3f("P", &npi);
			const Point2f *uvs = params.FindPoint2f("uv", &nuvi);
			if (!uvs) uvs = params.FindPoint2f("st", &nuvi);
			std::vector<Point2f> tempUVs;
			if (!uvs) {
				const Float *fuv = params.FindFloat("uv", &nuvi);
				if (!fuv) fuv = params.FindFloat("st", &nuvi);
				if (fuv) {
					nuvi /= 2;
					tempUVs.reserve(nuvi);
					for (int i = 0; i < nuvi; ++i)
					tempUVs.push_back(Point2f(fuv[2 * i], fuv[2 * i + 1]));
					uvs = &tempUVs[0];
				}
			}
			if (uvs) {
				if (nuvi < npi) {
					Error(
					"Not enough of \"uv\"s for triangle mesh.  Expected %d, "
					"found %d.  Discarding.",
					npi, nuvi);
					uvs = nullptr;
				} else if (nuvi > npi)
				Warning(
				"More \"uv\"s provided than will be used for triangle "
				"mesh.  (%d expcted, %d found)",
				npi, nuvi);
			}
			if (!vi) {
				Error(
				"Vertex indices \"indices\" not provided with triangle mesh shape");
				return std::vector<std::shared_ptr<Shape>>();
			}
			if (!P) {
				Error("Vertex positions \"P\" not provided with triangle mesh shape");
				return std::vector<std::shared_ptr<Shape>>();
			}
			const Vector3f *S = params.FindVector3f("S", &nsi);
			if (S && nsi != npi) {
				Error("Number of \"S\"s for triangle mesh must match \"P\"s");
				S = nullptr;
			}
			const Normal3f *N = params.FindNormal3f("N", &nni);
			if (N && nni != npi) {
				Error("Number of \"N\"s for triangle mesh must match \"P\"s");
				N = nullptr;
			}
			for (int i = 0; i < nvi; ++i)
			if (vi[i] >= npi) {
				Error(
				"trianglemesh has out of-bounds vertex index %d (%d \"P\" "
				"values were given",
				vi[i], npi);
				return std::vector<std::shared_ptr<Shape>>();
			}

			std::shared_ptr<Texture<Float>> alphaTex;
			std::string alphaTexName = params.FindTexture("alpha");
			if (alphaTexName != "") {
				if (floatTextures->find(alphaTexName) != floatTextures->end())
				alphaTex = (*floatTextures)[alphaTexName];
				else
				Error("Couldn't find float texture \"%s\" for \"alpha\" parameter",
				alphaTexName.c_str());
			} else if (params.FindOneFloat("alpha", 1.f) == 0.f)
			alphaTex.reset(new ConstantTexture<Float>(0.f));

			std::shared_ptr<Texture<Float>> shadowAlphaTex;
			std::string shadowAlphaTexName = params.FindTexture("shadowalpha");
			if (shadowAlphaTexName != "") {
				if (floatTextures->find(shadowAlphaTexName) != floatTextures->end())
				shadowAlphaTex = (*floatTextures)[shadowAlphaTexName];
				else
				Error(
				"Couldn't find float texture \"%s\" for \"shadowalpha\" "
				"parameter",
				shadowAlphaTexName.c_str());
			} else if (params.FindOneFloat("shadowalpha", 1.f) == 0.f)
			shadowAlphaTex.reset(new ConstantTexture<Float>(0.f));

			return CreateTriangleMesh(o2w, w2o, reverseOrientation, nvi / 3, vi, npi, P,
			S, N, uvs, alphaTex, shadowAlphaTex);
		}

	}
	// namespace pbrt
