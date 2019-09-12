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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_ACCELERATORS_BVH_H
#define PBRT_ACCELERATORS_BVH_H

// accelerators/bvh.h*
#include "pbrt.h"
#include "primitive.h"
#include <atomic>

namespace pbrt {

// BVHAccel Forward Declarations
struct BVHBuildNode;
struct BVHPrimitiveInfo;
struct MortonPrimitive;
struct LinearBVHNode;
class BVHSubtree;

// BVHAccel Declarations
class BVHAccel: public Aggregate {
public:
	// BVHAccel Public Types
	enum class SplitMethod {
		SAH, HLBVH, Middle, EqualCounts
	};

	// BVHAccel Public Methods
	BVHAccel(const std::vector<std::shared_ptr<const Primitive>> &p,
			int maxPrimsInNode = 1, SplitMethod splitMethod = SplitMethod::SAH);
	Bounds3f WorldBound() const;
	~BVHAccel();

	inline bool Intersect(const Ray &ray, SurfaceInteraction *isect) const
			override {
		return IntersectFrom(ray, isect, 0);
	}
	inline bool IntersectP(const Ray &ray) const override {
		return IntersectPFrom(ray, 0);
	}
	inline bool IntersectP(const Ray &ray, uint32_t ignoreId) const override {
		return IntersectPFrom(ray, ignoreId, 0);
	}
	inline bool IntersectL(const LineVisibilityQuery3f &query,
			const Normal3f &normal, Float time,
			IntervalCollection& collection) const {
		return IntersectLFrom(query, normal, time, collection, 0);
	}
	bool IntersectAll(const Ray& ray, GroupIsects& intersections,
			int groupSize = 1) const {
		return IntersectAllFrom(ray, intersections, groupSize, 0);
	}

	inline uint32_t Count() const override {
		return count;
	}

	inline Float Area() const override {
		return area;
	}
private:
// BVHAccel Private Methods
	BVHBuildNode *recursiveBuild(MemoryArena &arena,
			std::vector<BVHPrimitiveInfo> &primitiveInfo, int start, int end,
			int *totalNodes,
			std::vector<std::shared_ptr<const Primitive>> &orderedPrims);
	BVHBuildNode *HLBVHBuild(MemoryArena &arena,
			const std::vector<BVHPrimitiveInfo> &primitiveInfo, int *totalNodes,
			std::vector<std::shared_ptr<const Primitive>> &orderedPrims) const;
	BVHBuildNode *emitLBVH(BVHBuildNode *&buildNodes,
			const std::vector<BVHPrimitiveInfo> &primitiveInfo,
			MortonPrimitive *mortonPrims, int nPrimitives, int *totalNodes,
			std::vector<std::shared_ptr<const Primitive>> &orderedPrims,
			std::atomic<int> *orderedPrimsOffset, int bitIndex) const;
	BVHBuildNode *buildUpperSAH(MemoryArena &arena,
			std::vector<BVHBuildNode *> &treeletRoots, int start, int end,
			int *totalNodes) const;

	bool IntersectFrom(const Ray &ray, SurfaceInteraction *isect,
			uint32_t from) const;
	bool IntersectPFrom(const Ray &ray, uint32_t from) const;
	bool IntersectPFrom(const Ray &ray, uint32_t ignoreId, uint32_t from) const;
	bool IntersectLFrom(const LineVisibilityQuery3f &query,
			const Normal3f &normal, Float time, IntervalCollection& collection,
			uint32_t from) const;
	bool IntersectAllFrom(const Ray& ray, GroupIsects& intersections,
			int groupSize = 1, uint32_t from = 0) const;
	bool IntersectAllNoGroupFrom(const Ray& ray, GroupIsects& intersections,
			uint32_t from = 0) const;

	int flattenBVHTree(BVHBuildNode *node, int *offset);

// BVHAccel Private Data
	const int maxPrimsInNode;
	const SplitMethod splitMethod;
	std::vector<std::shared_ptr<const Primitive>> primitives;
	std::vector<std::shared_ptr<const BVHSubtree>> subtrees;
	LinearBVHNode *nodes = nullptr;

	uint32_t count = 0;
	Float area = 0;
	friend class BVHSubtree;
};

class BVHSubtree: public Aggregate {
public:
	BVHSubtree(const BVHAccel* bvh, uint32_t offset) :
			Aggregate(), bvh(bvh), offset(offset) {
	}

	Bounds3f WorldBound() const override;
	uint32_t Count() const override;
	Float Area() const override;

	bool Intersect(const Ray& ray, SurfaceInteraction *isect) const;
	bool IntersectP(const Ray& ray) const override;
	bool IntersectP(const Ray& ray, uint32_t ignoreId) const override;
private:
	const BVHAccel* bvh;
	const uint32_t offset;

};

std::shared_ptr<const BVHAccel> CreateBVHAccelerator(
		const std::vector<std::shared_ptr<const Primitive>> &prims,
		const ParamSet &ps);

}  // namespace pbrt

#endif  // PBRT_ACCELERATORS_BVH_H
