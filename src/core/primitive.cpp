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

// core/primitive.cpp*
#include "primitive.h"
#include "light.h"
#include "interaction.h"
#include "stats.h"

namespace pbrt {

STAT_PERCENT("Primitive/GeometricPrimitive IntersectP() tests", nPHits,
		nPTests);
STAT_PERCENT("Primitive/GeometricPrimitive IntersectL() tests", nLHits,
		nLTests);

using linesampling::LineVisibilityQuery3f;
using linesampling::IntervalCollection;

std::mutex TransformedPrimitiveMutex;
std::map<uint32_t, uint32_t> TransformedPrimitiveRemapping;
std::atomic<uint32_t> Primitive::nextPrimitiveIdentifier { 0 };

// Primitive Method Definitions
Primitive::~Primitive() {
}
const AreaLight *Aggregate::GetAreaLight() const {
	LOG(FATAL)<<
	"Aggregate::GetAreaLight() method"
	"called; should have gone to GeometricPrimitive";
	return nullptr;
}

const Material *Aggregate::GetMaterial() const {
	LOG(FATAL)<<
	"Aggregate::GetMaterial() method"
	"called; should have gone to GeometricPrimitive";
	return nullptr;
}

void Aggregate::ComputeScatteringFunctions(SurfaceInteraction *isect,
		MemoryArena &arena, TransportMode mode, bool allowMultipleLobes) const {
	LOG(FATAL)<<
	"Aggregate::ComputeScatteringFunctions() method"
	"called; should have gone to GeometricPrimitive";
}

// TransformedPrimitive Method Definitions
bool TransformedPrimitive::Intersect(const Ray &r,
		SurfaceInteraction *isect) const {
	// Compute _ray_ after transformation by _PrimitiveToWorld_
	Transform InterpolatedPrimToWorld;
	PrimitiveToWorld.Interpolate(r.time, &InterpolatedPrimToWorld);
	const Ray& ray = Inverse(InterpolatedPrimToWorld)(r);
	if (!primitive->Intersect(ray, isect))
		return false;
	r.tMax = ray.tMax;

	// Transform instance's intersection data to world space
	if (!InterpolatedPrimToWorld.IsIdentity())
		*isect = InterpolatedPrimToWorld(*isect);

	isect->instance = this;

	CHECK_GE(Dot(isect->n, isect->shading.n), 0);
	return true;
}

bool TransformedPrimitive::IntersectP(const Ray &r) const {
	Transform InterpolatedPrimToWorld;
	PrimitiveToWorld.Interpolate(r.time, &InterpolatedPrimToWorld);
	Transform InterpolatedWorldToPrim = Inverse(InterpolatedPrimToWorld);

	const Ray& ray = InterpolatedWorldToPrim(r);
	const bool result = primitive->IntersectP(ray);

	r.boxIsectPTests += ray.boxIsectPTests;
	r.triangleIsectPTests += ray.triangleIsectPTests;
	r.geometricPrimitiveIsectPTests += ray.geometricPrimitiveIsectPTests;

	return result;
}

bool TransformedPrimitive::IntersectL(const LineVisibilityQuery3f& query,
		const Normal3f &n, Float time, IntervalCollection &intersection) const {
	Transform InterpolatedPrimToWorld;
	PrimitiveToWorld.Interpolate(time, &InterpolatedPrimToWorld);
	Transform InterpolatedWorldToPrim = Inverse(InterpolatedPrimToWorld);

	const Point3f& tp0 = InterpolatedWorldToPrim(query.p0);
	const Point3f& tp1 = InterpolatedWorldToPrim(query.p1);
	const Point3f& tp2 = InterpolatedWorldToPrim(query.p2);
	const Normal3f& tn = InterpolatedWorldToPrim(n);
	const LineVisibilityQuery3f tquery(tp0, tp1, tp2);
	return primitive->IntersectL(tquery, tn, time, intersection);
}

bool TransformedPrimitive::IntersectP(const Ray &r, uint32_t ignoreId) const {
	if (ignoreId == 0)
		return IntersectP(r);
	else if (id == ignoreId)
		return false;

	Transform InterpolatedPrimToWorld;
	PrimitiveToWorld.Interpolate(r.time, &InterpolatedPrimToWorld);
	Transform InterpolatedWorldToPrim = Inverse(InterpolatedPrimToWorld);

	const Ray& ray = InterpolatedWorldToPrim(r);
	const bool result = primitive->IntersectP(ray, ignoreId);

	r.boxIsectPTests += ray.boxIsectPTests;
	r.triangleIsectPTests += ray.triangleIsectPTests;
	r.geometricPrimitiveIsectPTests += ray.geometricPrimitiveIsectPTests;
	return result;
}

Float TransformedPrimitive::Area() const {
	if (PrimitiveToWorld.HasScale()) {
		Float s = PrimitiveToWorld.GetScale();

		if (s < 0) {
			return 0;
		} else {
			return primitive->Area() * s * s;
		}
	} else {
		return primitive->Area();
	}
}

bool TransformedPrimitive::IntersectAll(const Ray& r,
		GroupIsects& intersections, int groupSize) const {
	// Compute _ray_ after transformation by _PrimitiveToWorld_
	Transform InterpolatedPrimToWorld;
	PrimitiveToWorld.Interpolate(r.time, &InterpolatedPrimToWorld);
	const Ray& ray = Inverse(InterpolatedPrimToWorld)(r);

	GroupIsects localIntersections;
	if (!primitive->IntersectAll(ray, localIntersections, groupSize))
		return false;

	const bool canReplaceGroup = Count() <= groupSize;

	for (const GroupIsect& localIntersection : localIntersections) {
		SurfaceInteraction isect(localIntersection.first);
		if (!InterpolatedPrimToWorld.IsIdentity())
			isect = InterpolatedPrimToWorld(isect);
		isect.instance = this;

		if (canReplaceGroup) {
			intersections.push_back(std::make_pair(isect, shared_from_this()));
			continue;
		}
		/*
		 * To avoid having excessive doubles in our occlusion map,
		 * we create a unique identifier for the transformed group.
		 *
		 * To do this, we use the __TransformedPrimitiveRemapping__ map
		 * which is can be accessed thread-safely via the
		 * __TransformedPrimitiveMutex__.
		 *
		 * When a new Transformed primitive needs to be returned, we
		 * perform a lookup in the map to see whether an instance has
		 * already been returned.
		 *
		 * When no instance has been returned, we generate a new unique
		 * identifier using __nextPrimitiveIdentifier__
		 *
		 * Otherwise, we return the id.
		 */

		mutex.lock();
		auto it = remapped.find(localIntersection.second->id);

		std::shared_ptr<const Primitive> primitive;
		if (it == remapped.end()) {
			// Create a transformed primitive
			const std::shared_ptr<const TransformedPrimitive> transformed =
					std::make_shared<const TransformedPrimitive>(
							localIntersection.second, PrimitiveToWorld);

			primitive = std::dynamic_pointer_cast<const Primitive>(transformed);
			std::pair<uint32_t, std::shared_ptr<const Primitive>> pair =
					std::make_pair(localIntersection.second->id, primitive);

			remapped.insert(pair);
		} else
			primitive = it->second;

		mutex.unlock();

		intersections.push_back(std::make_pair(isect, primitive));
	}

	return true;
}

// GeometricPrimitive Method Definitions
Bounds3f GeometricPrimitive::WorldBound() const {
	return shape->WorldBound();
}

bool GeometricPrimitive::IntersectP(const Ray &r) const {
	const bool result = shape->IntersectP(r);
	++r.geometricPrimitiveIsectPTests;
	++nPTests;
	if (result)
		++nPHits;
	return result;
}

bool GeometricPrimitive::IntersectL(const LineVisibilityQuery3f &query,
		const Normal3f &normal, Float time,
		IntervalCollection &intersection) const {
	const bool result = shape->IntersectL(query, normal, time, intersection);
	++nLTests;
	if (result)
		++nLHits;
	return result;
}

bool GeometricPrimitive::Intersect(const Ray &r,
		SurfaceInteraction *isect) const {
	Float tHit;
	if (!shape->Intersect(r, &tHit, isect))
		return false;
	r.tMax = tHit;
	isect->primitive = this;
	isect->instance = nullptr;
	isect->t = tHit;

	//	CHECK_GE(Dot(isect->n, isect->shading.n), 0.);
	// Initialize _SurfaceInteraction::mediumInterface_ after _Shape_
	// intersection
	if (mediumInterface.IsMediumTransition())
		isect->mediumInterface = mediumInterface;
	else
		isect->mediumInterface = MediumInterface(r.medium);
	return true;
}

bool GeometricPrimitive::IntersectAll(const Ray& r, GroupIsects& intersections,
		int groupSize) const {
	Float tHit;
	SurfaceInteraction isect;
	if (!shape->Intersect(r, &tHit, &isect))
		return false;

	isect.primitive = this;
	isect.instance = nullptr;
	isect.t = tHit;
//	CHECK_GE(Dot(isect.n, isect.shading.n), 0.);

	// Initialize _SurfaceInteraction::mediumInterface_ after _Shape_
	// intersection
	if (mediumInterface.IsMediumTransition())
		isect.mediumInterface = mediumInterface;
	else
		isect.mediumInterface = MediumInterface(r.medium);

	// Add the surface interaction to the intersections
	// Note: a geometric primitive is always equal to its group
	const GroupIsect& groupIsect = std::make_pair(isect, shared_from_this());
	intersections.push_back(groupIsect);

	return true;
}

const AreaLight *GeometricPrimitive::GetAreaLight() const {
	return areaLight.get();
}

const Material *GeometricPrimitive::GetMaterial() const {
	return material.get();
}

void GeometricPrimitive::ComputeScatteringFunctions(SurfaceInteraction *isect,
		MemoryArena &arena, TransportMode mode, bool allowMultipleLobes) const {
	ProfilePhase p(Prof::ComputeScatteringFuncs);
	if (material)
		material->ComputeScatteringFunctions(isect, arena, mode,
				allowMultipleLobes);
//	CHECK_GE(Dot(isect->n, isect->shading.n), 0.);
}

}  // namespace pbrt
