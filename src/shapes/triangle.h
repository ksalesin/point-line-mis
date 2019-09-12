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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SHAPES_TRIANGLE_H
#define PBRT_SHAPES_TRIANGLE_H

// shapes/triangle.h*
#include "shape.h"
#include "stats.h"
#include <map>
#include "point-line-mis/lineinteraction.h"

using linesampling::LineInteraction;

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/Triangle meshes", triMeshBytes);

// Triangle Declarations
struct TriangleMesh {
	// TriangleMesh Public Methods
	TriangleMesh(const Transform &ObjectToWorld, int nTriangles,
			const int *vertexIndices, int nVertices, const Point3f *P,
			const Vector3f *S, const Normal3f *N, const Point2f *uv,
			const std::shared_ptr<Texture<Float>> &alphaMask,
			const std::shared_ptr<Texture<Float>> &shadowAlphaMask);

	// TriangleMesh Data
	const int nTriangles, nVertices;
	std::vector<int> vertexIndices;
	std::unique_ptr<Point3f[]> p;
	std::unique_ptr<Normal3f[]> n;
	std::unique_ptr<Vector3f[]> s;
	std::unique_ptr<Point2f[]> uv;
	std::shared_ptr<Texture<Float>> alphaMask, shadowAlphaMask;
};

class Triangle: public Shape {
public:
	// Triangle Public Methods
	Triangle(const Transform *ObjectToWorld, const Transform *WorldToObject,
			bool reverseOrientation, const std::shared_ptr<TriangleMesh> &mesh,
			int triNumber) :
			Shape(ObjectToWorld, WorldToObject, reverseOrientation), mesh(mesh) {
		v = &mesh->vertexIndices[3 * triNumber];
		triMeshBytes += sizeof(*this);
	}

	inline const Point3f& operator[](int index) const {
		return mesh->p[v[index]];
	}

	Bounds3f ObjectBound() const;
	Bounds3f WorldBound() const;
	bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
			bool testAlphaTexture = true) const;
	bool IntersectP(const Ray &ray, bool testAlphaTexture = true) const;
	bool IntersectL(const LineVisibilityQuery3f &query, const Normal3f &normal,
			Float time, IntervalCollection &intersection) const;
	Float Area() const;

	using Shape::Sample;  // Bring in the other Sample() overload.
	Interaction Sample(const Point2f &u, Float *pdf) const;

	// Returns the solid angle subtended by the triangle w.r.t. the given
	// reference point p.
	Float SolidAngle(const Point3f &p, int nSamples = 0) const;

	bool CanEvaluateAnalytically(const SurfaceInteraction& isect,
			Spectrum * spectrum = nullptr) const override;
	bool CanIlluminate(const Interaction& isect) const;

	LineInteraction LineSample(const Point2f &u,
			bool importanceSample = true) const override;
	
	Float PdfLine(const Interaction &ref, const Vector3f &wi, Point3f &pLight,
			LineInteraction &lineInteraction, const Float& lineAngle, bool importanceSample) const override;

	// Visualizes lines on triangle in given direction (same as how line lighting would be evaluated)
	bool OnColoredLine(const Interaction &intr, float direction) const override;

private:
	// Triangle Private Methods
	void GetUVs(Point2f uv[3]) const {
		if (mesh->uv) {
			uv[0] = mesh->uv[v[0]];
			uv[1] = mesh->uv[v[1]];
			uv[2] = mesh->uv[v[2]];
		} else {
			uv[0] = Point2f(0, 0);
			uv[1] = Point2f(1, 0);
			uv[2] = Point2f(1, 1);
		}
	}

	// Triangle Private Data
	std::shared_ptr<TriangleMesh> mesh;
	const int *v;
};

Float SolidAngleGirard(const Point3f&p, const Point3f& p0, const Point3f& p1,
		const Point3f& p2);

std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(const Transform *o2w,
		const Transform *w2o, bool reverseOrientation, int nTriangles,
		const int *vertexIndices, int nVertices, const Point3f *p,
		const Vector3f *s, const Normal3f *n, const Point2f *uv,
		const std::shared_ptr<Texture<Float>> &alphaTexture,
		const std::shared_ptr<Texture<Float>> &shadowAlphaTexture);
std::vector<std::shared_ptr<Shape>> CreateTriangleMeshShape(
		const Transform *o2w, const Transform *w2o, bool reverseOrientation,
		const ParamSet &params,
		std::map<std::string, std::shared_ptr<Texture<Float>>>*floatTextures =
nullptr);

bool WritePlyFile(const std::string &filename, int nTriangles,
		const int *vertexIndices, int nVertices, const Point3f *P,
		const Vector3f *S, const Normal3f *N, const Point2f *UV);

}  // namespace pbrt

#endif  // PBRT_SHAPES_TRIANGLE_H
