/*
 * list.cpp
 *
 *  Created on: Jul 16, 2017
 *      Author: niels
 */

#include "list.h"
#include "paramset.h"

namespace pbrt {

ListAggregate::ListAggregate(
		const std::vector<std::shared_ptr<const Primitive>> &p) :
		primitives(p) {
	for (const std::shared_ptr<const Primitive>& primitive : p) {
		count += primitive->Count();
		area += primitive->Area();
		bounds = Union(bounds, primitive->WorldBound());
	}
}

inline bool ListAggregate::Intersect(const Ray &ray,
		SurfaceInteraction *isect) const {
	if (!bounds.IntersectP(ray))
		return false;

	bool hit = false;
	for (const std::shared_ptr<const Primitive>& p : primitives)
		if (p->Intersect(ray, isect))
			hit = true;
	return hit;
}

bool ListAggregate::IntersectP(const Ray &ray) const {
	if (!bounds.IntersectP(ray))
		return false;

	for (const std::shared_ptr<const Primitive>& p : primitives)
		if (p->IntersectP(ray))
			return true;
	return false;
}

bool ListAggregate::IntersectP(const Ray &ray, uint32_t ignoreId) const {
	if (ignoreId == 0)
		return IntersectP(ray);
	else if (id == ignoreId)
		return false;

	if (!bounds.IntersectP(ray))
		return false;

	for (const std::shared_ptr<const Primitive>& p : primitives)
		if (p->IntersectP(ray, ignoreId))
			return true;
	return false;
}

bool ListAggregate::IntersectAll(const Ray& ray, GroupIsects& intersections,
		int groupSize) const {
	if (!bounds.IntersectP(ray))
		return false;

	bool hit = false;
	for (const std::shared_ptr<const Primitive>& p : primitives)
		if (p->IntersectAll(ray, intersections, groupSize))
			hit = true;
	return hit;
}

Bounds3f ListAggregate::WorldBound() const {
	return bounds;
}

std::shared_ptr<const ListAggregate> CreateListAccelerator(
		const std::vector<std::shared_ptr<const Primitive>> &prims,
		const ParamSet &ps) {
	return std::make_shared<const ListAggregate>(prims);
}

}
