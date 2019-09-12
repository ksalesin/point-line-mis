/*
 * list.h
 *
 *  Created on: Jul 16, 2017
 *      Author: niels
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef ACCELERATORS_LIST_H
#define ACCELERATORS_LIST_H

// accelerators/lists.h*
#include "pbrt.h"
#include "primitive.h"
#include "interaction.h"

namespace pbrt {

// ListAggregate Declarations
class ListAggregate: public Aggregate {
public:
	ListAggregate(const std::vector<std::shared_ptr<const Primitive>> &p);

	Bounds3f WorldBound() const;
	bool Intersect(const Ray &ray, SurfaceInteraction *isect) const override;
	bool IntersectP(const Ray &ray) const override;
	bool IntersectP(const Ray &ray, uint32_t ignoreId) const override;
	bool IntersectAll(const Ray& ray, GroupIsects& intersections,
			int groupSize = 1) const override;

	inline uint32_t Count() const override {
		return count;
	}
	inline Float Area() const override {
		return area;
	}
private:
	const std::vector<std::shared_ptr<const Primitive>> primitives;
	Bounds3f bounds;
	uint32_t count = 0;
	Float area = 0;
};

std::shared_ptr<const ListAggregate> CreateListAccelerator(
		const std::vector<std::shared_ptr<const Primitive>> &prims,
		const ParamSet &ps);

}

#endif /* ACCELERATORS_LIST_H_ */
