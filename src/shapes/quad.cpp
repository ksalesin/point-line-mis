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

// shapes/quad.cpp*
#include "shapes/quad.h"
#include "shapes/triangle.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"
#include "efloat.h"
#include "reflection.h"
#include "point-line-mis/triangle.h"

namespace pbrt {

STAT_PERCENT("Intersections/Ray-quad Intersect() tests", nQuadHits, nQuadTests);
STAT_PERCENT("Intersections/Ray-quad IntersectP() tests", nQuadShadowHits,
		nQuadShadowTests);

Bounds3f Quad::ObjectBound() const {
	// Get quad vertices in _p0_, _p1_, _p2_, and _p3_
	const Point3f &p0 = (*WorldToObject)(GetPoint(0));
	const Point3f &p1 = (*WorldToObject)(GetPoint(1));
	const Point3f &p2 = (*WorldToObject)(GetPoint(2));
	const Point3f &p3 = (*WorldToObject)(GetPoint(3));
	return Union(Bounds3f(p0, p1), Bounds3f(p2, p3));
}

Bounds3f Quad::WorldBound() const {
	// Get quad vertices in _p0_, _p1_, _p2_, and _p3_
	const Point3f &p0 = GetPoint(0);
	const Point3f &p1 = GetPoint(1);
	const Point3f &p2 = GetPoint(2);
	const Point3f &p3 = GetPoint(3);
	return Union(Bounds3f(p0, p1), Bounds3f(p2, p3));
}

bool Quad::Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
		bool testAlphaTexture) const {
	// Compute $\VEC{s}_1$
	++nQuadTests;

	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point3f& p0 = GetPoint(0);
	const Point3f& p1 = GetPoint(1);
	const Point3f& p2 = GetPoint(3);
	const Vector3f& e1 = p1 - p0;
	const Vector3f& e2 = p2 - p0;
	const Vector3f& s1 = Cross(ray.d, e2);
	const Float divisor = Dot(s1, e1);

	if (divisor == 0.)
		return false;
	const Float invDivisor = 1.f / divisor;

	// Compute first barycentric coordinate
	const Vector3f& s = ray.o - p0;
	Float b1 = Dot(s, s1) * invDivisor;
	if (b1 < 0. || b1 > 1.)
		return false;

	// Compute second barycentric coordinate
	const Vector3f& s2 = Cross(s, e1);
	const Float b2 = Dot(ray.d, s2) * invDivisor;
	if (b2 < 0. || b2 > 1.)
		return false;

	// Compute _t_ to intersection point
	const Float t = Dot(e2, s2) * invDivisor;
	if (t < 0 || t > ray.tMax)
		return false;

	// Compute triangle partial derivatives
	const Float b0 = 1.0 - b1 - b2;

	// Compute error bounds for triangle intersection
	const Float xAbsSum = (std::abs(b0 * p0.x) + std::abs(b1 * p1.x)
			+ std::abs(b2 * p2.x));
	const Float yAbsSum = (std::abs(b0 * p0.y) + std::abs(b1 * p1.y)
			+ std::abs(b2 * p2.y));
	const Float zAbsSum = (std::abs(b0 * p0.z) + std::abs(b1 * p1.z)
			+ std::abs(b2 * p2.z));
	const Vector3f& pError = gamma(7) * Vector3f(xAbsSum, yAbsSum, zAbsSum);

	// Interpolate $(u,v)$ parametric coordinates and hit point
	const Point3f& pHit = b0 * p0 + b1 * p1 + b2 * p2;
	const Point2f uvHit(b1 + b2, b2);

	// Fill in _SurfaceInteraction_ from triangle hit
	*isect = SurfaceInteraction(pHit, pError, uvHit, -ray.d, u, v,
			Normal3f(0, 0, 0), Normal3f(0, 0, 0), ray.time, this);

	// Override surface normal in _isect_ for triangle
	isect->n = isect->shading.n = Normal3f(Normalize(Cross(u, v)));

	// Ensure correct orientation of the normal (this breaks refraction!)
	if (Dot(ray.d, isect->n) > 0)
		isect->n = -isect->n;
	if (Dot(isect->n, isect->shading.n) < 0)
		isect->shading.n = -isect->shading.n;

//	// Ensure correct orientation of the geometric normal
//	if (reverseOrientation ^ transformSwapsHandedness)
//		isect->n = isect->shading.n = -isect->n;

	*tHit = t;
	++nQuadHits;
	return true;
}

bool Quad::IntersectP(const Ray &ray, bool testAlphaTexture) const {
	// Compute $\VEC{s}_1$
	++nQuadShadowTests;

	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point3f& p0 = GetPoint(0);
	const Point3f& p1 = GetPoint(1);
	const Point3f& p2 = GetPoint(3);
	const Vector3f& e1 = p1 - p0;
	const Vector3f& e2 = p2 - p0;
	const Vector3f& s1 = Cross(ray.d, e2);
	const Float divisor = Dot(s1, e1);

	if (divisor == 0.)
		return false;
	const Float invDivisor = 1.f / divisor;

	// Compute first barycentric coordinate
	const Vector3f& s = ray.o - p0;
	Float b1 = Dot(s, s1) * invDivisor;
	if (b1 < 0. || b1 > 1.)
		return false;

	// Compute second barycentric coordinate
	const Vector3f& s2 = Cross(s, e1);
	const Float b2 = Dot(ray.d, s2) * invDivisor;
	if (b2 < 0. || b2 > 1.)
		return false;

	// Compute _t_ to intersection point
	const Float t = Dot(e2, s2) * invDivisor;
	if (t < 0 || t > ray.tMax)
		return false;

	++nQuadShadowHits;
	return true;
}

Float Quad::Area() const {
	return Cross(u, v).Length();
}

Point2f Quad::GetUV(const Point3f& p) const {
	const Vector3f op = p - this->o;
	
	Float du = Dot(op, this->un);
	Float dv = Dot(op, this->vn);

	Float fu = du / this->eu;
	Float fv = dv / this->ev;

	return Point2f(fu,fv);
}

Interaction Quad::SampleUV(const Point2f& uv) const {
	// could use the actual uv coordinates to do something smarter here
	Interaction it;

	const Point3f& p0 = GetPoint(0);
	const Point3f& p1 = GetPoint(1);
	const Point3f& p2 = GetPoint(2);

	it.p = this->o + (uv.x * this->u) + (uv.y * this->v);
	// Compute surface normal for sampled point on triangle
	const Vector3f& n = Normalize(Cross(this->u, this->v));

	it.n = Normal3f(n);
	// Ensure correct orientation of the geometric normal; follow the same
	// approach as was used in Triangle::Intersect().
	if (reverseOrientation ^ transformSwapsHandedness)
		it.n *= -1;

	Point3f pAbsSum = Abs(u[0] * p0) + Abs(u[1] * p1)
			+ Abs((1 - u[0] - u[1]) * p2);
	it.pError = gamma(6) * Vector3f(pAbsSum.x, pAbsSum.y, pAbsSum.z);

	return it;
}


Interaction Quad::Sample(const Point2f &u, Float *pdf) const {
	Interaction it;

	const Point3f& p0 = GetPoint(0);
	const Point3f& p1 = GetPoint(1);
	const Point3f& p2 = GetPoint(2);

	it.p = this->o + (u.x * this->u) + (u.y * this->v);
	// Compute surface normal for sampled point on triangle
	const Vector3f& n = Normalize(Cross(this->u, this->v));

	it.n = Normal3f(n);
	// Ensure correct orientation of the geometric normal; follow the same
	// approach as was used in Triangle::Intersect().
	if (reverseOrientation ^ transformSwapsHandedness)
		it.n *= -1;

	Point3f pAbsSum = Abs(u[0] * p0) + Abs(u[1] * p1)
			+ Abs((1 - u[0] - u[1]) * p2);
	it.pError = gamma(6) * Vector3f(pAbsSum.x, pAbsSum.y, pAbsSum.z);

	*pdf = 1 / Area();
	return it;
}

Float Quad::SolidAngle(const Point3f &p, int nSamples) const {
	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point3f& p0 = GetPoint(0);
	const Point3f& p1 = GetPoint(1);
	const Point3f& p2 = GetPoint(2);
	const Point3f& p3 = GetPoint(3);

	return SolidAngleGirard(p, p0, p1, p2) + SolidAngleGirard(p, p0, p2, p3);
}

bool Quad::CanEvaluateAnalytically(const SurfaceInteraction& isect,
		Spectrum * spectrum) const {
	// Quickly terminate when analytical integration is impossible to avoid clipping
	if (!isect.bsdf->CanEvaluateAnalytically())
		return false;

	if (spectrum) {
		// Get triangle vertices in _p1_, _p2_, and _p3_
		const Point3f& p0 = GetPoint(0);
		const Point3f& p1 = GetPoint(1);
		const Point3f& p2 = GetPoint(2);
		const Point3f& p3 = GetPoint(3);

		niels::Triangle t0(p0, p1, p3);
		niels::Triangle t1(p3, p1, p2);

		// Clip against the surface normal
		std::vector<niels::Triangle> geometryNormalClip;
		t0.Clip(isect.p, isect.n, geometryNormalClip);
		t1.Clip(isect.p, isect.n, geometryNormalClip);

		Spectrum sum(0.f);
		for (const niels::Triangle& triangle : geometryNormalClip) {
			Spectrum analytical(0.f);
			if (!isect.bsdf->CanEvaluateAnalytically(isect, triangle,
					&analytical))
				return false;
			sum += analytical;
		}
		*spectrum = sum;
	}

	return true;
}

bool Quad::CanIlluminate(const Interaction& it) const {
	// check if point is below this quad
	const Float dot = Dot(GetPoint(0) - it.p, Cross(u, v));
	if (reverseOrientation ^ transformSwapsHandedness) {
		if (dot <= 0)
			return false;
	} else {
		if (dot >= 0)
			return false;
	}

	for (int i = 0; i < 4; ++i)
		if (Dot(GetPoint(i) - it.p, it.n) > 0)
			return true;
	return false;
}

LineInteraction Quad::LineSample(const Point2f &u, bool importanceSample) const {

	// Get edges of quad, create a quad-aligned orthonormal basis (n1, w, normal)
	const Vector3f e1 = this->u;
	const Vector3f e2 = this->v;
	const Vector3f e3 = this->u + this->v;
	const Vector3f cross = Cross(e1, e2);
	const Float area = cross.Length();

	const Vector3f n1 = Normalize(e1);
	const Vector3f n2 = Normalize(e2);
	Normal3f normal = Normal3f(Normalize(cross));
	normal = reverseOrientation ^ transformSwapsHandedness ? -normal : normal;

	const Vector3f w = Normalize(Cross(cross, e1));

	// Get line direction angle (offset from n1)
	const Float angle = Pi * u.x;
	const Float x1 = std::cos(angle);
	const Float x2 = std::sin(angle);

	// Create an angle-aligned orthonormal basis (x, y, normal)
	const Vector3f& x = Normalize(x1 * n1 + x2 * w);
	const Vector3f& y = Normalize(Cross(x, normal));

	// Sanity checks
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

	// Project vertices onto x/y basis with origin o
	const Point2f proj_p0(0.0, 0.0);
	const Point2f proj_p1(Dot(e1, x), Dot(e1, y));
	const Point2f proj_p2(Dot(e2, x), Dot(e2, y));
	const Point2f proj_p3(Dot(e3, x), Dot(e3, y));

	std::vector<Point2f> proj_p{proj_p0, proj_p1, proj_p2, proj_p3};

	// Put vertices in order along x-axis
	std::sort(proj_p.begin(), proj_p.end(), [](const Point2f& lhs, const Point2f& rhs)
	{
		return lhs.x < rhs.x;
	});

	const Point2f pMin = proj_p[0];
	const Point2f pMid1 = proj_p[1];
	const Point2f pMid2 = proj_p[2];
	const Point2f pMax = proj_p[3];

	// Find line segment given random offset u.y along x-axis
	Float lx, ly1, ly2, pdf;

	// importance sample line segments such that pdf is uniform over area
	if (importanceSample) {
		Float invRange = 1.0 / (pMax.x - pMin.x);
		Float alpha = (pMid1.x - pMin.x) * invRange;
		Float beta = (pMid2.x - pMin.x) * invRange;
		Float gamma = 1 + beta - alpha;

		// transform random offset to distribution proportional to length of line
		Float uu;
		if (u.y < alpha) {
			uu = std::sqrt(u.y * alpha * gamma);
		} else if (u.y < beta) {
			uu = (u.y * gamma + alpha) / 2;
		} else {
			uu = 1 - std::sqrt((1-beta) * gamma * (1-u.y));
		}

		lx = pMin.x * (1.0 - uu) + pMax.x * uu;

		if (lx < pMid1.x) {
			Float s1 = (lx - pMin.x) / (pMid1.x - pMin.x);
			ly1 = pMin.y * (1-s1) + pMid1.y * s1;
			Float s2 = (lx - pMin.x) / (pMid2.x - pMin.x);
			ly2 = pMin.y * (1-s2) + pMid2.y * s2;

		} else if (lx < pMid2.x) {
			Float s1 = (lx - pMid1.x) / (pMax.x - pMid1.x);
			ly1 = pMid1.y * (1-s1) + pMax.y * s1;
			Float s2 = (lx - pMin.x) / (pMid2.x - pMin.x);
			ly2 = pMin.y * (1-s2) + pMid2.y * s2;

		} else {
			Float s1 = (lx - pMid1.x) / (pMax.x - pMid1.x);
			ly1 = pMid1.y * (1-s1) + pMax.y * s1;
			Float s2 = (lx - pMid2.x) / (pMax.x - pMid2.x);
			ly2 = pMid2.y * (1-s2) + pMax.y * s2;
		}

		// pdf proportional to length of line
		pdf = std::abs(ly1 - ly2) / area;

	} else {
		// linearly interpolate uniformly along offset axis
		lx = pMin.x * (1.0 - u.y) + pMax.x * u.y;

		if (lx < pMid1.x) {
			Float s1 = (lx - pMin.x) / (pMid1.x - pMin.x);
			ly1 = pMin.y * (1-s1) + pMid1.y * s1;
			Float s2 = (lx - pMin.x) / (pMid2.x - pMin.x);
			ly2 = pMin.y * (1-s2) + pMid2.y * s2;

		} else if (lx < pMid2.x) {
			Float s1 = (lx - pMid1.x) / (pMax.x - pMid1.x);
			ly1 = pMid1.y * (1-s1) + pMax.y * s1;
			Float s2 = (lx - pMin.x) / (pMid2.x - pMin.x);
			ly2 = pMin.y * (1-s2) + pMid2.y * s2;

		} else {
			Float s1 = (lx - pMid1.x) / (pMax.x - pMid1.x);
			ly1 = pMid1.y * (1-s1) + pMax.y * s1;
			Float s2 = (lx - pMid2.x) / (pMax.x - pMid2.x);
			ly2 = pMid2.y * (1-s2) + pMax.y * s2;
		}

		pdf = 1.0 / (pMax.x - pMin.x);
	}

	LineInteraction result;
	
	// Calculate final line segment endpoints in quad space
	const Point3f& lo = this->o + x * lx;
	result.line.p1 = lo + ly1 * y;
	result.line.p2 = lo + ly2 * y;
	result.pdf = pdf;
	result.length = std::abs(ly1 - ly2);
	result.n = normal;

	return result;
}

bool Quad::OnColoredLine(const Interaction &intr, float direction) const {

	// Get edges of quad, create a quad-aligned orthonormal basis (n1, w, normal)
	const Vector3f e1 = this->u;
	const Vector3f e2 = this->v;
	const Vector3f e3 = this->u + this->v;
	const Vector3f cross = Cross(e1, e2);
	const Float area = cross.Length();

	const Vector3f n1 = Normalize(e1);
	const Vector3f n2 = Normalize(e2);
	Normal3f normal = Normal3f(Normalize(cross));
	normal = reverseOrientation ^ transformSwapsHandedness ? -normal : normal;
	const Vector3f w = Normalize(Cross(cross, e1));

	// Get line direction angle (offset from n1)
	const Float angle = Pi * direction;
	const Float x1 = std::cos(angle);
	const Float x2 = std::sin(angle);

	// Create an angle-aligned orthonormal basis (x, y, normal)
	const Vector3f x = Normalize(x1 * n1 + x2 * w);
	const Vector3f y = Normalize(Cross(x, normal));

	// Sanity checks
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

	// Project vertices onto x/y basis with origin o
	const Point2f proj_p0(0.0, 0.0);
	const Point2f proj_p1(Dot(e1, x), Dot(e1, y));
	const Point2f proj_p2(Dot(e2, x), Dot(e2, y));
	const Point2f proj_p3(Dot(e3, x), Dot(e3, y));

	std::vector<Point2f> proj_p{proj_p0, proj_p1, proj_p2, proj_p3};

	// Put vertices in order along x-axis
	std::sort(proj_p.begin(), proj_p.end(), [](const Point2f& lhs, const Point2f& rhs)
	{
		return lhs.x < rhs.x;
	});

	const Point2f pMin = proj_p[0];
	const Point2f pMax = proj_p[3];

	// Project intersection point onto x/y basis
	const Vector3f e4 = intr.p - this->o;
	float itsOffset = Dot(e4, x);

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

Float Quad::PdfLine(const Interaction &ref, const Vector3f &wi, Point3f &pLight,
			LineInteraction &lineInteraction, const Float& lineAngle, bool importanceSample) const {
	
	// Get point on quad light in wi direction
	Ray ray = ref.SpawnRay(wi);
	Float tHit;
	SurfaceInteraction isectLight;
	if (!Intersect(ray, &tHit, &isectLight, false))
		return 0.f;

	pLight = isectLight.p;

	// Get edges of quad, create a quad-aligned orthonormal basis (n1, w, normal)
	const Vector3f e1 = this->u;
	const Vector3f e2 = this->v;
	const Vector3f e3 = this->u + this->v;
	const Vector3f cross = Cross(e1, e2);
	const Float area = cross.Length();

	const Vector3f n1 = Normalize(e1);
	const Vector3f n2 = Normalize(e2);
	Normal3f normal = Normal3f(Normalize(cross));
	normal = reverseOrientation ^ transformSwapsHandedness ? -normal : normal;
	const Vector3f w = Normalize(Cross(cross, e1));

	// Get line direction angle (offset from n1)
	const Float angle = Pi * lineAngle;
	const Float x1 = std::cos(angle);
	const Float x2 = std::sin(angle);

	// Create an angle-aligned orthonormal basis (x, y, normal)
	const Vector3f x = Normalize(x1 * n1 + x2 * w);
	const Vector3f y = Normalize(Cross(x, normal));

	// Sanity checks
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

	// Project vertices and sample point onto x/y basis with origin o
	const Point2f proj_p0(0.0, 0.0);
	const Point2f proj_p1(Dot(e1, x), Dot(e1, y));
	const Point2f proj_p2(Dot(e2, x), Dot(e2, y));
	const Point2f proj_p3(Dot(e3, x), Dot(e3, y));

	const Vector3f vp = pLight - (*this).o;
	const Point2f proj_pLight(Dot(vp, x), Dot(vp, y));

	std::vector<Point2f> proj_p{proj_p0, proj_p1, proj_p2, proj_p3};

	// Put vertices in order along x-axis
	std::sort(proj_p.begin(), proj_p.end(), [](const Point2f& lhs, const Point2f& rhs)
	{
		return lhs.x < rhs.x;
	});

	const Point2f pMin = proj_p[0];
	const Point2f pMid1 = proj_p[1];
	const Point2f pMid2 = proj_p[2];
	const Point2f pMax = proj_p[3];

	// Get line offset parameter
	double u = (proj_pLight.x - pMin.x) / (pMax.x - pMin.x);

	// Find line segment given random offset u.y along x-axis
	Float lx, ly1, ly2, pdf;

	lx = pMin.x * (1.0 - u) + pMax.x * u;

	if (lx < pMid1.x) {
		Float s1 = (lx - pMin.x) / (pMid1.x - pMin.x);
		ly1 = pMin.y * (1-s1) + pMid1.y * s1;
		Float s2 = (lx - pMin.x) / (pMid2.x - pMin.x);
		ly2 = pMin.y * (1-s2) + pMid2.y * s2;

	} else if (lx < pMid2.x) {
		Float s1 = (lx - pMid1.x) / (pMax.x - pMid1.x);
		ly1 = pMid1.y * (1-s1) + pMax.y * s1;
		Float s2 = (lx - pMin.x) / (pMid2.x - pMin.x);
		ly2 = pMin.y * (1-s2) + pMid2.y * s2;

	} else {
		Float s1 = (lx - pMid1.x) / (pMax.x - pMid1.x);
		ly1 = pMid1.y * (1-s1) + pMax.y * s1;
		Float s2 = (lx - pMid2.x) / (pMax.x - pMid2.x);
		ly2 = pMid2.y * (1-s2) + pMax.y * s2;
	}

	if (importanceSample) {
		// pdf proportional to length of line
		pdf = std::abs(ly1 - ly2) / area;
	} else {
		// linearly interpolate uniformly along offset axis
		pdf = 1.0 / (pMax.x - pMin.x);
	}

	// Calculate final line segment endpoints in triangle space
	const Point3f& lo = (*this).o + x * lx;
	lineInteraction.line.p1 = lo + ly1 * y;
	lineInteraction.line.p2 = lo + ly2 * y;
	lineInteraction.pdf = pdf;
	lineInteraction.length = std::abs(ly1 - ly2);
	lineInteraction.n = normal;

	return pdf;
}

// Adapted from "An Area-Preserving Parameterization for Spherical Rectangles" 
// by Carlos Ureña, Marcos Fajardo, and Alan King
// https://www.arnoldrenderer.com/research/egsr2013_spherical_rectangle.pdf
Interaction Quad::SampleSolidAngle(const Point3f& p, const Point2f &sample, Float* pdf) const {
	Interaction it;

	// Local coordinate system
	const Point3f o = p;
	const Vector3f x = Normalize(this->u);
	const Vector3f y = Normalize(this->v);
	      Vector3f z = Cross(x, y);

	// Constants
	Vector3f d = this->o - p;
	Float z0 = Dot(d, z);
	if (z0 > 0) {
		z *= -1;
		z0 *= -1;
	}

	Float z0sq = z0 * z0;
	Float x0 = Dot(d, x);
	Float y0 = Dot(d, y);
	Float x1 = x0 + this->u.Length();
	Float y1 = y0 + this->v.Length();
	Float y0sq = y0 * y0;
	Float y1sq = y1 * y1;

	// Vectors to four vertices of quad
	const Vector3f v00(x0, y0, z0);
	const Vector3f v01(x0, y1, z0);
	const Vector3f v10(x1, y0, z0);
	const Vector3f v11(x1, y1, z0);

	// Normals of edges of sph quad
	const Vector3f n0 = Normalize(Cross(v00, v10));
	const Vector3f n1 = Normalize(Cross(v10, v11));
	const Vector3f n2 = Normalize(Cross(v11, v01));
	const Vector3f n3 = Normalize(Cross(v01, v00));

	// Internal angles of sph quad
	Float g0 = std::acos(-Dot(n0, n1));
	Float g1 = std::acos(-Dot(n1, n2));
	Float g2 = std::acos(-Dot(n2, n3));
	Float g3 = std::acos(-Dot(n3, n0));

	// More constants
	Float b0 = n0.z;
	Float b1 = n2.z;
	Float b0sq = b0 * b0;
	Float k = 2 * M_PI - g2 - g3;

	// Solid angle of sph squad (pdf uniform over SA)
	Float SA = g0 + g1 - k;
	*pdf = 1 / SA;

	// Compute 'cu'
	Float u = sample.x;
	Float v = sample.y;

	Float au = u * SA + k;
	Float fu = (std::cos(au) * b0 - b1) / std::sin(au);
	Float sign = (fu > 0) ? 1 : -1;
	Float cu = sign * 1 / std::sqrt(fu*fu + b0sq);
	cu = std::max((Float)-1., std::min(cu, (Float)1.));

	// Compute 'xu'
	Float xu = -(cu * z0) / std::sqrt(1 - cu*cu);
	xu = std::max(x0, std::min(xu, x1));

	// Compute 'yv'
	Float r = std::sqrt(xu*xu + z0sq);
	Float h0 = y0 / std::sqrt(r*r + y0sq);
	Float h1 = y1 / std::sqrt(r*r + y1sq);
	Float hv = h0 + v * (h1 - h0);
	Float hv2 = hv * hv;
	Float yv = (hv2 < 1-EPSILON) ? (hv * r) / std::sqrt(1-hv2) : y1;

	// Transform (xu,yv,z0) to world coords
	Point3f soln = o + xu * x + yv * y + z0 * z;

	it.p = soln;
	it.n = Normal3f(z.x, z.y, z.z);
	return it;
}

// Adapted from "An Area-Preserving Parameterization for Spherical Rectangles" 
// by Carlos Ureña, Marcos Fajardo, and Alan King
// https://www.arnoldrenderer.com/research/egsr2013_spherical_rectangle.pdf
Float Quad::PdfSolidAngle(const Point3f& p) const {
	// Local coordinate system
	const Point3f o = p;
	const Vector3f x = Normalize(this->u);
	const Vector3f y = Normalize(this->v);
	      Vector3f z = Cross(x, y);

	// Constants
	Vector3f d = this->o - p;
	Float z0 = Dot(d, z);
	if (z0 > 0) {
		z *= -1;
		z0 *= -1;
	}

	Float z0sq = z0 * z0;
	Float x0 = Dot(d, x);
	Float y0 = Dot(d, y);
	Float x1 = x0 + this->u.Length();
	Float y1 = y0 + this->v.Length();
	Float y0sq = y0 * y0;
	Float y1sq = y1 * y1;

	// Vectors to four vertices of quad
	const Vector3f v00(x0, y0, z0);
	const Vector3f v01(x0, y1, z0);
	const Vector3f v10(x1, y0, z0);
	const Vector3f v11(x1, y1, z0);

	// Normals of edges of sph quad
	const Vector3f n0 = Normalize(Cross(v00, v10));
	const Vector3f n1 = Normalize(Cross(v10, v11));
	const Vector3f n2 = Normalize(Cross(v11, v01));
	const Vector3f n3 = Normalize(Cross(v01, v00));

	// Internal angles of sph quad
	Float g0 = std::acos(-Dot(n0, n1));
	Float g1 = std::acos(-Dot(n1, n2));
	Float g2 = std::acos(-Dot(n2, n3));
	Float g3 = std::acos(-Dot(n3, n0));

	// More constants
	Float b0 = n0.z;
	Float b1 = n2.z;
	Float b0sq = b0 * b0;
	Float k = 2 * M_PI - g2 - g3;

	// Solid angle of sph squad (pdf uniform over SA)
	Float SA = g0 + g1 - k;

	return 1 / SA;
}

// Adapted from "An Area-Preserving Parameterization for Spherical Rectangles" 
// by Carlos Ureña, Marcos Fajardo, and Alan King
// https://www.arnoldrenderer.com/research/egsr2013_spherical_rectangle.pdf
LineInteraction Quad::SampleSQLine(const Point3f& p, const Float u, int strategy) const {
	// Local coordinate system
	const Point3f o = p;
	Vector3f x, y, z;
	Float ex, ey;

	// which edge of quad light to align lines with
	if (strategy == 1) {
		x = this->un;
		y = this->vn;
		z = this->wn;

		ex = this->eu;
		ey = this->ev;
	} else {
		x = this->vn;
		y = this->un;
		z = this->wn;

		ex = this->ev;
		ey = this->eu;
	}

	// Constants
	Vector3f d = this->o - o;
	Float z0 = Dot(d, z);
	if (z0 > 0) {
		z *= -1;
		z0 *= -1;
	}

	Float z0sq = z0 * z0;
	Float x0 = Dot(d, x);
	Float y0 = Dot(d, y);
	Float x1 = x0 + ex;
	Float y1 = y0 + ey;
	Float y0sq = y0 * y0;
	Float y1sq = y1 * y1;

	// Vectors to four vertices of quad
	const Vector3f v00(x0, y0, z0);
	const Vector3f v01(x0, y1, z0);
	const Vector3f v10(x1, y0, z0);
	const Vector3f v11(x1, y1, z0);

	// Normals of edges of sph quad
	const Vector3f n0 = Normalize(Cross(v00, v10));
	const Vector3f n1 = Normalize(Cross(v10, v11));
	const Vector3f n2 = Normalize(Cross(v11, v01));
	const Vector3f n3 = Normalize(Cross(v01, v00));

	// Internal angles of sph quad
	Float g0 = std::acos(-Dot(n0, n1));
	Float g1 = std::acos(-Dot(n1, n2));
	Float g2 = std::acos(-Dot(n2, n3));
	Float g3 = std::acos(-Dot(n3, n0));

	// More constants
	Float b0 = n0.z;
	Float b1 = n2.z;
	Float b0sq = b0 * b0;
	Float k = 2 * M_PI - g2 - g3;

	// Solid angle of sph squad (pdf uniform over SA)
	Float SA = g0 + g1 - k;

	// Compute 'cu'
	Float au = u * SA + k;
	Float fu = (std::cos(au) * b0 - b1) / std::sin(au);
	Float sign = (fu > 0) ? 1 : -1;
	Float cu = sign * 1 / std::sqrt(fu*fu + b0sq);
	cu = std::max((Float)-1., std::min(cu, (Float)1.));

	// Compute 'xu'
	Float xu = -(cu * z0) / std::sqrt(1 - cu*cu);
	xu = std::max(x0, std::min(xu, x1));

	// Transform endpoints of xu segment to world space
	Point3f p1 = o + xu * x + y0 * y + z0 * z;
	Point3f p2 = o + xu * x + y1 * y + z0 * z;

	// This is joint pdf(xu, yv), which is later used to find marginal pdf(xu)
	Float pdf = 1 / SA;

	LineInteraction result;
	result.line.p1 = p1;
	result.line.p2 = p2;
	result.pdf = pdf;
	result.length = y1-y0;
	result.n = Normal3f(z);
	result.xu = xu;

	return result;
}

// Adapted from "An Area-Preserving Parameterization for Spherical Rectangles" 
// by Carlos Ureña, Marcos Fajardo, and Alan King
// https://www.arnoldrenderer.com/research/egsr2013_spherical_rectangle.pdf
Point3f Quad::SampleSQPoint(const IntervalCollection intervals, const LineInteraction &interaction, 
		const Point3f& p, const Float v, Float *pointPdf, int strategy) const {

	// Local coordinate system
	const Point3f o = p;
	Vector3f x, y, z;
	Float ex, ey;

	// which edge of quad light to align lines with
	if (strategy == 1) {
		x = this->un;
		y = this->vn;
		z = this->wn;

		ex = this->eu;
		ey = this->ev;
	} else {
		x = this->vn;
		y = this->un;
		z = this->wn;

		ex = this->ev;
		ey = this->eu;
	}

	// Constants
	Vector3f op = this->o - o;
	Float z0 = Dot(op, z);
	if (z0 > 0) {
		z *= -1;
		z0 *= -1;
	}

	Float z0sq = z0 * z0;
	Float x0 = Dot(op, x);
	Float y0 = Dot(op, y);
	Float y0sq = y0*y0;
	Float x1 = x0 + ex;
	Float y1 = y0 + ey;
	Float y1sq = y1*y1;

	// Retrieve 'xu'
	Float xu = interaction.xu;

	Float dsq = xu*xu + z0sq;
	Float d = std::sqrt(dsq);

	// find visible region heights [h0,h1] along y-axis in 'sph quad space'
	std::vector<Point2f> visHeights;
	Float sum_heights = 0.;

	Point3f start = interaction.line.p1;
	Point3f end = interaction.line.p2;
	Vector3f ld = end - start;

	Float h0 = y0 / std::sqrt(dsq + y0*y0);
	Float h1 = y1 / std::sqrt(dsq + y1*y1);

	IntervalConstIterator it = intervals.begin();
	while (it != intervals.end()) {
		// find visible segment endpoints in world space
		Point3f v1 = start + it->min * ld;
		Point3f v2 = start + it->max * ld;

		// find y-coordinates of visible segment endpoints in 'sph quad space'
		Float yv1 = Dot(v1-o, y);
		Float yv2 = Dot(v2-o, y);

		Float hv0 = yv1 / std::sqrt(dsq + yv1*yv1);
		Float hv1 = yv2 / std::sqrt(dsq + yv2*yv2);

		// add heights to list
		visHeights.push_back(Point2f(hv0, hv1));
		sum_heights += std::abs(hv1 - hv0);

		it++;
	}

	// lerp between all visible heights
	Float hv;
	Float sum = 0.;

	for (int i = 0; i < visHeights.size(); i++) {
		Float h1 = visHeights[i].x;
		Float h2 = visHeights[i].y;

		Float dh = std::abs(h2-h1) / sum_heights;
		Float newsum = sum + dh;

		if (v <= newsum) {
			// Stretch v to [0,1] within region
			Float v_lerp = (v - sum) / dh;

			// Interpolate along region
			hv = h1 + v_lerp * (h2 - h1);
			break;
		}
		sum = newsum;
	}

	// sample point (xu, yv, z0) on quad light
	Float hv2 = hv*hv;
	Float yv = (hv*d) / std::sqrt(1-hv2);
	yv = std::min(y1, std::max(yv, y0)); // clamp to [y0,y1]

	// condensed form of second part of marginal line pdf * conditional point pdf
	*pointPdf = (sum_heights > 1e-6) ? (h1-h0) / sum_heights : 0.;

	// Transform (xu, yv, z0) to world coordinates
	return o + xu*x + yv*y + z0*z;
}

Float Quad::PdfSQLine(const Interaction &ref, const Vector3f& wi, Point3f& pLight, LineInteraction &interaction, int strategy) const {
	
	// Get point on quad light in wi direction
	Ray ray = ref.SpawnRay(wi);
	Float tHit;
	SurfaceInteraction isectLight;
	if (!Intersect(ray, &tHit, &isectLight, false))
		return 0.f;

	pLight = isectLight.p;
	
	// Local coordinate system
	const Point3f o = ref.p;
	Vector3f x, y, z;
	Float ex, ey;

	if (strategy == 1) {
		x = this->un;
		y = this->vn;
		z = this->wn;

		ex = this->eu;
		ey = this->ev;
	} else {
		x = this->vn;
		y = this->un;
		z = this->wn;

		ex = this->ev;
		ey = this->eu;
	}

	// Constants
	Vector3f d = this->o - o;
	Float z0 = Dot(d, z);
	if (z0 > 0) {
		z *= -1;
		z0 *= -1;
	}

	Float z0sq = z0 * z0;
	Float x0 = Dot(d, x);
	Float y0 = Dot(d, y);
	Float x1 = x0 + ex;
	Float y1 = y0 + ey;
	Float y0sq = y0 * y0;
	Float y1sq = y1 * y1;

	// Vectors to four vertices of quad
	const Vector3f v00(x0, y0, z0);
	const Vector3f v01(x0, y1, z0);
	const Vector3f v10(x1, y0, z0);
	const Vector3f v11(x1, y1, z0);

	// Normals of edges of sph quad
	const Vector3f n0 = Normalize(Cross(v00, v10));
	const Vector3f n1 = Normalize(Cross(v10, v11));
	const Vector3f n2 = Normalize(Cross(v11, v01));
	const Vector3f n3 = Normalize(Cross(v01, v00));

	// Internal angles of sph quad
	Float g0 = std::acos(-Dot(n0, n1));
	Float g1 = std::acos(-Dot(n1, n2));
	Float g2 = std::acos(-Dot(n2, n3));
	Float g3 = std::acos(-Dot(n3, n0));

	// More constants
	Float b0 = n0.z;
	Float b1 = n2.z;
	Float b0sq = b0 * b0;
	Float k = 2 * M_PI - g2 - g3;

	// Solid angle of sph squad (pdf uniform over SA)
	Float SA = g0 + g1 - k;

	// Recover xu in 'sph quad space'
	Float px = pLight.x * x.x + pLight.y * x.y + pLight.z * x.z;
	Float ox = o.x * x.x + o.y * x.y + o.z * x.z;
	Float xu = px - ox;
	xu = std::max(x0, std::min(xu, x1));

	// Transform endpoints of xu segment to world space
	Point3f p1 = o + xu * x + y0 * y + z0 * z;
	Point3f p2 = o + xu * x + y1 * y + z0 * z;

	// This is joint pdf(xu, yv), which is later used to find marginal pdf(xu)
	Float pdf = 1 / SA;

	interaction.line.p1 = p1;
	interaction.line.p2 = p2;
	interaction.pdf = pdf;
	interaction.length = y1-y0;
	interaction.n = Normal3f(z);
	interaction.xu = xu;

	return pdf;
}

Float Quad::PdfSQPoint(const IntervalCollection intervals, const LineInteraction &interaction, 
		const Point3f& p, int strategy) const {

	// No need to check for visibility, would have exited early in integrator

	// Local coordinate system
	const Point3f o = p;
	Vector3f x, y, z;
	Float ex, ey;

	if (strategy == 1) {
		x = this->un;
		y = this->vn;
		z = this->wn;

		ex = this->eu;
		ey = this->ev;
	} else {
		x = this->vn;
		y = this->un;
		z = this->wn;

		ex = this->ev;
		ey = this->eu;
	}

	// Constants
	Vector3f op = this->o - o;
	Float z0 = Dot(op, z);
	if (z0 > 0) {
		z *= -1;
		z0 *= -1;
	}

	Float z0sq = z0 * z0;
	Float x0 = Dot(op, x);
	Float y0 = Dot(op, y);
	Float y0sq = y0*y0;
	Float x1 = x0 + ex;
	Float y1 = y0 + ey;
	Float y1sq = y1*y1;

	// Retrieve 'xu'
	Float xu = interaction.xu;

	Float dsq = xu*xu + z0sq;
	Float d = std::sqrt(dsq);

	// find visible region heights [h0,h1] along y-axis in 'sph quad space'
	Float sum_heights = 0.;

	Point3f start = interaction.line.p1;
	Point3f end = interaction.line.p2;
	Vector3f ld = end - start;

	Float h0 = y0 / std::sqrt(dsq + y0*y0);
	Float h1 = y1 / std::sqrt(dsq + y1*y1);

	IntervalConstIterator it = intervals.begin();
	while (it != intervals.end()) {
		// find visible segment endpoints in world space
		Point3f v1 = start + it->min * ld;
		Point3f v2 = start + it->max * ld;

		// find y-coordinates of visible segment endpoints in 'sph quad space'
		Float yv1 = Dot(v1-o, y);
		Float yv2 = Dot(v2-o, y);

		Float hv0 = yv1 / std::sqrt(dsq + yv1*yv1);
		Float hv1 = yv2 / std::sqrt(dsq + yv2*yv2);

		// add heights to list
		sum_heights += std::abs(hv1 - hv0);

		it++;
	}

	// condensed form of second part of marginal line pdf * conditional point pdf
	return (sum_heights > 1e-6) ? (h1-h0) / sum_heights : 0.;
}

Float smoothstepEaseIn(Float x) {
    if (x <= 0)
        return 0.;

    if (x >= 1)
        return 1.;

    return -(x-2) * x;
}

Float smoothstepCubic(Float x) {
    if (x <= 0)
        return 0.;

    if (x >= 1)
        return 1.;

    Float x2 = x * x;
    return 3*x2 - 2*x*x2;
}

Float Quad::smoothstep(const Point3f& pLight, const Float smoothstepWidth) const {
	// Get smoothstep distance (normalized over half distance from start to end)
	Float ulen = this->u.Length();
	Float vlen = this->v.Length();

	Float uhalf = ulen/2.;
	Float vhalf = vlen/2.;

	Vector3f po = pLight - this->o;

	Float udot = Dot(po,this->u);
	Float udist = std::min(udot, ulen-udot);
	
	Float vdot = Dot(po,this->v);
	Float vdist = std::min(vdot, vlen-vdot);

	// Float uratio = udist / uhalf; // normalized distance
	// Float vratio = vdist / vhalf; // normalized distance
	Float uratio = udist; // global distance
	Float vratio = vdist; // global distance

	if (uhalf < EPSILON || uhalf < EPSILON || vdist < EPSILON || udist < EPSILON)
		return 0.; // avoid nans

    Float xu = uratio / smoothstepWidth;
    Float xv = vratio / smoothstepWidth;

	if (isnan(xu) || isnan(xv))
		return 0.;

	Float ssu = smoothstepEaseIn(xu);
	Float ssv = smoothstepEaseIn(xv);

	return ssu * ssv; // multiply
    // return std::min(ssu, ssv); // min
}

std::shared_ptr<Shape> CreateQuadShape(const Transform *o2w,
		const Transform *w2o, bool reverseOrientation, const ParamSet &params) {
	const Point3f& o = params.FindOnePoint3f("o", Point3f(0, 0, 0));
	const Vector3f& u = params.FindOneVector3f("u", Vector3f(1, 0, 0));
	const Vector3f& v = params.FindOneVector3f("v", Vector3f(0, 0, 1));

	const Point3f& ow = (*o2w)(o);
	const Vector3f& uw = (*o2w)(u);
	const Vector3f& vw = (*o2w)(v);

	// Init spherical quad constants
	const Vector3f& un = Normalize(uw);
	const Vector3f& vn = Normalize(vw);
	const Vector3f& wn = Normalize(Cross(uw,vw));

	const double eu = uw.Length();
	const double ev = vw.Length();

	return std::make_shared<Quad>(o2w, w2o, reverseOrientation, ow, uw, vw, un, vn, wn, eu, ev);
}

}
// namespace pbrt
