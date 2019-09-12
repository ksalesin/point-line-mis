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

// accelerators/bvh.cpp*
#include "bvh.h"

#include "accelerators/bvh.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "parallel.h"
#include <algorithm>

#include "point-line-mis/trianglebox.h"

using linesampling::IntersectTriangleBox;

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/BVH tree", treeBytes);
STAT_PERCENT("BVH/Ray-box Intersect() tests", nHits, nTests);
STAT_PERCENT("BVH/Ray-box IntersectP() tests", nShadowHits, nShadowTests);
STAT_PERCENT("BVH/Triangle-box IntersectP() tests", nTriangleBoxHits,
		nTriangleBoxTests);
STAT_RATIO("BVH/Primitives per leaf node", totalPrimitives, totalLeafNodes);
STAT_COUNTER("BVH/Interior nodes", interiorNodes);
STAT_COUNTER("BVH/Leaf nodes", leafNodes);

// BVHAccel Local Declarations
struct BVHPrimitiveInfo {
	BVHPrimitiveInfo() {
	}
	BVHPrimitiveInfo(size_t primitiveNumber, const Bounds3f &bounds) :
			primitiveNumber(primitiveNumber), bounds(bounds), centroid(
					.5f * bounds.pMin + .5f * bounds.pMax) {
	}
	size_t primitiveNumber;
	Bounds3f bounds;
	Point3f centroid;
};

struct BVHBuildNode {
	// BVHBuildNode Public Methods
	void InitLeaf(int first, int nPrim, uint32_t n, const Bounds3f &b,
			Float a) {
		firstPrimOffset = first;
		nPrimitives = nPrim;
		bounds = b;
		children[0] = children[1] = nullptr;
		count = n;
		area = a;
		++leafNodes;
		++totalLeafNodes;
		totalPrimitives += n;
	}
	void InitInterior(int axis, BVHBuildNode *c0, BVHBuildNode *c1) {
		children[0] = c0;
		children[1] = c1;
		if (c0) {
			count += c0->count;
			area += c0->area;
		}
		if (c1) {
			count += c1->count;
			area += c1->area;
		}
		bounds = Union(c0->bounds, c1->bounds);
		splitAxis = axis;
		nPrimitives = 0;
		++interiorNodes;
	}

	int Depth() {
		if (nPrimitives > 0)
			return 0;
		else {
			int leftDepth = children[0] ? children[0]->Depth() : 0;
			int rightDepth = children[1] ? children[1]->Depth() : 0;
			return std::max(leftDepth, rightDepth) + 1;
		}
	}
	Bounds3f bounds;
	BVHBuildNode *children[2];
	int splitAxis, firstPrimOffset, nPrimitives;
	uint32_t count = 0;
	Float area = 0;
};

struct MortonPrimitive {
	int primitiveIndex;
	uint32_t mortonCode;
};

struct LBVHTreelet {
	int startIndex, nPrimitives;
	BVHBuildNode *buildNodes;
};

struct LinearBVHNode {
	Bounds3f bounds;
	union {
		int primitivesOffset;   // leaf
		int secondChildOffset;  // interior
	};
	uint16_t nPrimitives;  // 0 -> interior node
	uint8_t axis;          // interior node: xyz
	uint32_t count;
	Float area;
	//uint8_t pad[1];        // ensure 32 byte total size
	uint8_t depth; // maximum depth of 255
	// size:
	// bounds: 48 bytes
	// union:   4 bytes
	// nprim:   2 bytes
	// axis     1 byte
	// count    4 bytes
	// area     4 bytes
	// pad      1 bytes
};

// BVHAccel Utility Functions
inline uint32_t LeftShift3(uint32_t x) {
	CHECK_LE(x, (1 << 10));
	if (x == (1 << 10))
		--x;
#ifdef PBRT_HAVE_BINARY_CONSTANTS
	x = (x | (x << 16)) & 0b00000011000000000000000011111111;
	// x = ---- --98 ---- ---- ---- ---- 7654 3210
	x = (x | (x << 8)) & 0b00000011000000001111000000001111;
	// x = ---- --98 ---- ---- 7654 ---- ---- 3210
	x = (x | (x << 4)) & 0b00000011000011000011000011000011;
	// x = ---- --98 ---- 76-- --54 ---- 32-- --10
	x = (x | (x << 2)) & 0b00001001001001001001001001001001;
	// x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
#else
	x = (x | (x << 16)) & 0x30000ff;
	// x = ---- --98 ---- ---- ---- ---- 7654 3210
	x = (x | (x << 8)) & 0x300f00f;
	// x = ---- --98 ---- ---- 7654 ---- ---- 3210
	x = (x | (x << 4)) & 0x30c30c3;
	// x = ---- --98 ---- 76-- --54 ---- 32-- --10
	x = (x | (x << 2)) & 0x9249249;
	// x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
#endif // PBRT_HAVE_BINARY_CONSTANTS
	return x;
}

inline uint32_t EncodeMorton3(const Vector3f &v) {
	CHECK_GE(v.x, 0);
	CHECK_GE(v.y, 0);
	CHECK_GE(v.z, 0);
	return (LeftShift3(v.z) << 2) | (LeftShift3(v.y) << 1) | LeftShift3(v.x);
}

static void RadixSort(std::vector<MortonPrimitive> *v) {
	std::vector<MortonPrimitive> tempVector(v->size());
	PBRT_CONSTEXPR int bitsPerPass = 6;
	PBRT_CONSTEXPR int nBits = 30;
	static_assert((nBits % bitsPerPass) == 0,
			"Radix sort bitsPerPass must evenly divide nBits");
	PBRT_CONSTEXPR int nPasses = nBits / bitsPerPass;

	for (int pass = 0; pass < nPasses; ++pass) {
		// Perform one pass of radix sort, sorting _bitsPerPass_ bits
		int lowBit = pass * bitsPerPass;

		// Set in and out vector pointers for radix sort pass
		std::vector<MortonPrimitive> &in = (pass & 1) ? tempVector : *v;
		std::vector<MortonPrimitive> &out = (pass & 1) ? *v : tempVector;

		// Count number of zero bits in array for current radix sort bit
		PBRT_CONSTEXPR int nBuckets = 1 << bitsPerPass;
		int bucketCount[nBuckets] = { 0 };
		PBRT_CONSTEXPR int bitMask = (1 << bitsPerPass) - 1;
		for (const MortonPrimitive &mp : in) {
			int bucket = (mp.mortonCode >> lowBit) & bitMask;
			CHECK_GE(bucket, 0);
			CHECK_LT(bucket, nBuckets);
			++bucketCount[bucket];
		}

		// Compute starting index in output array for each bucket
		int outIndex[nBuckets];
		outIndex[0] = 0;
		for (int i = 1; i < nBuckets; ++i)
			outIndex[i] = outIndex[i - 1] + bucketCount[i - 1];

		// Store sorted values in output array
		for (const MortonPrimitive &mp : in) {
			int bucket = (mp.mortonCode >> lowBit) & bitMask;
			out[outIndex[bucket]++] = mp;
		}
	}
	// Copy final result from _tempVector_, if needed
	if (nPasses & 1)
		std::swap(*v, tempVector);
}

// BVHAccel Method Definitions
BVHAccel::BVHAccel(const std::vector<std::shared_ptr<const Primitive>> &p,
		int maxPrimsInNode, SplitMethod splitMethod) :
		maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod), primitives(
				p) {

	ProfilePhase _(Prof::AccelConstruction);
	if (primitives.empty())
		return;
	// Build BVH from _primitives_

	// Initialize _primitiveInfo_ array for primitives
	std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
	for (size_t i = 0; i < primitives.size(); ++i) {
		primitiveInfo[i] = {i, primitives[i]->WorldBound()};
		count += primitives[i]->Count();
		area += primitives[i]->Area();
	}

	// Build BVH tree for primitives using _primitiveInfo_
	MemoryArena arena(1024 * 1024);
	int totalNodes = 0;
	std::vector<std::shared_ptr<const Primitive>> orderedPrims;
	orderedPrims.reserve(primitives.size());
	BVHBuildNode *root;
	if (splitMethod == SplitMethod::HLBVH)
		root = HLBVHBuild(arena, primitiveInfo, &totalNodes, orderedPrims);
	else
		root = recursiveBuild(arena, primitiveInfo, 0, primitives.size(),
				&totalNodes, orderedPrims);
	primitives.swap(orderedPrims);
	LOG(INFO)<< StringPrintf("BVH created with %d nodes for %d "
			"primitives (%.2f MB)", totalNodes,
			(int)primitives.size(),
			float(totalNodes * sizeof(LinearBVHNode)) /
			(1024.f * 1024.f));

	// Compute representation of depth-first traversal of BVH tree
	treeBytes += totalNodes * sizeof(LinearBVHNode) + sizeof(*this)
			+ primitives.size() * sizeof(primitives[0]);
	nodes = AllocAligned<LinearBVHNode>(totalNodes);
	int offset = 0;
	flattenBVHTree(root, &offset);
	CHECK_EQ(totalNodes, offset);

	// create the subtrees
	for (int i = 0; i < totalNodes; ++i)
		subtrees.push_back(std::make_shared<BVHSubtree>(this, i));
}

Bounds3f BVHAccel::WorldBound() const {
	return nodes ? nodes[0].bounds : Bounds3f();
}

struct BucketInfo {
	int count = 0;
	Bounds3f bounds;
};

BVHBuildNode *BVHAccel::recursiveBuild(MemoryArena &arena,
		std::vector<BVHPrimitiveInfo> &primitiveInfo, int start, int end,
		int *totalNodes,
		std::vector<std::shared_ptr<const Primitive>> &orderedPrims) {
	CHECK_NE(start, end);
	BVHBuildNode *node = arena.Alloc<BVHBuildNode>();
	(*totalNodes)++;
	// Compute bounds of all primitives in BVH node
	Bounds3f bounds;
	for (int i = start; i < end; ++i)
		bounds = Union(bounds, primitiveInfo[i].bounds);
	int nPrimitives = end - start;
	if (nPrimitives == 1) {
		// Create leaf _BVHBuildNode_
		int firstPrimOffset = orderedPrims.size();
		Float area = 0;
		uint32_t count = 0;
		for (int i = start; i < end; ++i) {
			int primNum = primitiveInfo[i].primitiveNumber;
			orderedPrims.push_back(primitives[primNum]);
			area += primitives[primNum]->Area();
			count += primitives[primNum]->Count();
		}
		node->InitLeaf(firstPrimOffset, nPrimitives, count, bounds, area);
		return node;
	} else {
		// Compute bound of primitive centroids, choose split dimension _dim_
		Bounds3f centroidBounds;
		for (int i = start; i < end; ++i)
			centroidBounds = Union(centroidBounds, primitiveInfo[i].centroid);
		int dim = centroidBounds.MaximumExtent();

		// Partition primitives into two sets and build children
		int mid = (start + end) / 2;
		if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
			// Create leaf _BVHBuildNode_
			int firstPrimOffset = orderedPrims.size();
			Float area = 0;
			uint32_t count = 0;
			for (int i = start; i < end; ++i) {
				int primNum = primitiveInfo[i].primitiveNumber;
				orderedPrims.push_back(primitives[primNum]);
				area += primitives[primNum]->Area();
				count += primitives[primNum]->Count();
			}
			node->InitLeaf(firstPrimOffset, nPrimitives, count, bounds, area);
			return node;
		} else {
			// Partition primitives based on _splitMethod_
			switch (splitMethod) {
			case SplitMethod::Middle: {
				// Partition primitives through node's midpoint
				Float pmid = (centroidBounds.pMin[dim]
						+ centroidBounds.pMax[dim]) / 2;
				BVHPrimitiveInfo *midPtr = std::partition(&primitiveInfo[start],
						&primitiveInfo[end - 1] + 1,
						[dim, pmid](const BVHPrimitiveInfo &pi) {
							return pi.centroid[dim] < pmid;
						});
				mid = midPtr - &primitiveInfo[0];
				// For lots of prims with large overlapping bounding boxes, this
				// may fail to partition; in that case don't break and fall
				// through
				// to EqualCounts.
				if (mid != start && mid != end)
					break;
			}
			case SplitMethod::EqualCounts: {
				// Partition primitives into equally-sized subsets
				mid = (start + end) / 2;
				std::nth_element(&primitiveInfo[start], &primitiveInfo[mid],
						&primitiveInfo[end - 1] + 1,
						[dim](const BVHPrimitiveInfo &a,
								const BVHPrimitiveInfo &b) {
							return a.centroid[dim] < b.centroid[dim];
						});
				break;
			}
			case SplitMethod::SAH:
			default: {
				// Partition primitives using approximate SAH
				if (nPrimitives <= 2) {
					// Partition primitives into equally-sized subsets
					mid = (start + end) / 2;
					std::nth_element(&primitiveInfo[start], &primitiveInfo[mid],
							&primitiveInfo[end - 1] + 1,
							[dim](const BVHPrimitiveInfo &a,
									const BVHPrimitiveInfo &b) {
								return a.centroid[dim] <
								b.centroid[dim];
							});
				} else {
					// Allocate _BucketInfo_ for SAH partition buckets
					const int nBuckets = std::max(24,
							std::min(128, 4 * nPrimitives));
					BucketInfo buckets[nBuckets];

					// Initialize _BucketInfo_ for SAH partition buckets
					for (int i = start; i < end; ++i) {
						int b = nBuckets
								* centroidBounds.Offset(
										primitiveInfo[i].centroid)[dim];
						if (b == nBuckets)
							b = nBuckets - 1;
						CHECK_GE(b, 0);
						CHECK_LT(b, nBuckets);
						buckets[b].count++;
						buckets[b].bounds = Union(buckets[b].bounds,
								primitiveInfo[i].bounds);
					}

					// Compute costs for splitting after each bucket
					Float cost[nBuckets - 1];
					for (int i = 0; i < nBuckets - 1; ++i) {
						Bounds3f b0, b1;
						int count0 = 0, count1 = 0;
						for (int j = 0; j <= i; ++j) {
							b0 = Union(b0, buckets[j].bounds);
							count0 += buckets[j].count;
						}
						for (int j = i + 1; j < nBuckets; ++j) {
							b1 = Union(b1, buckets[j].bounds);
							count1 += buckets[j].count;
						}
						cost[i] = 1
								+ (count0 * b0.SurfaceArea()
										+ count1 * b1.SurfaceArea())
										/ bounds.SurfaceArea();
					}

					// Find bucket to split at that minimizes SAH metric
					Float minCost = cost[0];
					int minCostSplitBucket = 0;
					for (int i = 1; i < nBuckets - 1; ++i) {
						if (cost[i] < minCost) {
							minCost = cost[i];
							minCostSplitBucket = i;
						}
					}

					// Either create leaf or split primitives at selected SAH
					// bucket
					Float leafCost = nPrimitives;
					if (nPrimitives > maxPrimsInNode || minCost < leafCost) {
						BVHPrimitiveInfo *pmid = std::partition(
								&primitiveInfo[start],
								&primitiveInfo[end - 1] + 1,
								[=](const BVHPrimitiveInfo &pi) {
									int b = nBuckets *
									centroidBounds.Offset(pi.centroid)[dim];
									if (b == nBuckets) b = nBuckets - 1;
									CHECK_GE(b, 0);
									CHECK_LT(b, nBuckets);
									return b <= minCostSplitBucket;
								});
						mid = pmid - &primitiveInfo[0];
					} else {
						// Create leaf _BVHBuildNode_
						int firstPrimOffset = orderedPrims.size();
						Float area = 0;
						uint32_t count;
						for (int i = start; i < end; ++i) {
							int primNum = primitiveInfo[i].primitiveNumber;
							orderedPrims.push_back(primitives[primNum]);
							area += primitives[primNum]->Area();
							count += primitives[primNum]->Count();
						}
						node->InitLeaf(firstPrimOffset, nPrimitives, count,
								bounds, area);
						return node;
					}
				}
				break;
			}
			}
			node->InitInterior(dim,
					recursiveBuild(arena, primitiveInfo, start, mid, totalNodes,
							orderedPrims),
					recursiveBuild(arena, primitiveInfo, mid, end, totalNodes,
							orderedPrims));
		}
	}
	return node;
}

BVHBuildNode *BVHAccel::HLBVHBuild(MemoryArena &arena,
		const std::vector<BVHPrimitiveInfo> &primitiveInfo, int *totalNodes,
		std::vector<std::shared_ptr<const Primitive>> &orderedPrims) const {
	// Compute bounding box of all primitive centroids
	Bounds3f bounds;
	for (const BVHPrimitiveInfo &pi : primitiveInfo)
		bounds = Union(bounds, pi.centroid);

	// Compute Morton indices of primitives
	std::vector<MortonPrimitive> mortonPrims(primitiveInfo.size());
	ParallelFor([&](int i) {
		// Initialize _mortonPrims[i]_ for _i_th primitive
			PBRT_CONSTEXPR int mortonBits = 10;
			PBRT_CONSTEXPR int mortonScale = 1 << mortonBits;
			mortonPrims[i].primitiveIndex = primitiveInfo[i].primitiveNumber;
			Vector3f centroidOffset = bounds.Offset(primitiveInfo[i].centroid);
			mortonPrims[i].mortonCode = EncodeMorton3(centroidOffset * mortonScale);
		}, primitiveInfo.size(), 512);

	// Radix sort primitive Morton indices
	RadixSort(&mortonPrims);

	// Create LBVH treelets at bottom of BVH

	// Find intervals of primitives for each treelet
	std::vector<LBVHTreelet> treeletsToBuild;
	for (int start = 0, end = 1; end <= (int) mortonPrims.size(); ++end) {
#ifdef PBRT_HAVE_BINARY_CONSTANTS
		uint32_t mask = 0b00111111111111000000000000000000;
#else
		uint32_t mask = 0x3ffc0000;
#endif
		if (end == (int) mortonPrims.size()
				|| ((mortonPrims[start].mortonCode & mask)
						!= (mortonPrims[end].mortonCode & mask))) {
			// Add entry to _treeletsToBuild_ for this treelet
			int nPrimitives = end - start;
			int maxBVHNodes = 2 * nPrimitives;
			BVHBuildNode *nodes = arena.Alloc < BVHBuildNode
					> (maxBVHNodes, false);
			treeletsToBuild.push_back( { start, nPrimitives, nodes });
			start = end;
		}
	}

	// Create LBVHs for treelets in parallel
	std::atomic<int> atomicTotal(0), orderedPrimsOffset(0);
	orderedPrims.resize(primitives.size());
	ParallelFor([&](int i) {
		// Generate _i_th LBVH treelet
			int nodesCreated = 0;
			const int firstBitIndex = 29 - 12;
			LBVHTreelet &tr = treeletsToBuild[i];
			tr.buildNodes =
			emitLBVH(tr.buildNodes, primitiveInfo, &mortonPrims[tr.startIndex],
					tr.nPrimitives, &nodesCreated, orderedPrims,
					&orderedPrimsOffset, firstBitIndex);
			atomicTotal += nodesCreated;
		}, treeletsToBuild.size());
	*totalNodes = atomicTotal;

	// Create and return SAH BVH from LBVH treelets
	std::vector<BVHBuildNode *> finishedTreelets;
	finishedTreelets.reserve(treeletsToBuild.size());
	for (LBVHTreelet &treelet : treeletsToBuild)
		finishedTreelets.push_back(treelet.buildNodes);
	return buildUpperSAH(arena, finishedTreelets, 0, finishedTreelets.size(),
			totalNodes);
}

BVHBuildNode *BVHAccel::emitLBVH(BVHBuildNode *&buildNodes,
		const std::vector<BVHPrimitiveInfo> &primitiveInfo,
		MortonPrimitive *mortonPrims, int nPrimitives, int *totalNodes,
		std::vector<std::shared_ptr<const Primitive>> &orderedPrims,
		std::atomic<int> *orderedPrimsOffset, int bitIndex) const {
	CHECK_GT(nPrimitives, 0);
	if (bitIndex == -1 || nPrimitives < maxPrimsInNode) {
		// Create and return leaf node of LBVH treelet
		(*totalNodes)++;
		BVHBuildNode *node = buildNodes++;
		Bounds3f bounds;
		int firstPrimOffset = orderedPrimsOffset->fetch_add(nPrimitives);
		Float area = 0;
		uint32_t count = 0;
		for (int i = 0; i < nPrimitives; ++i) {
			int primitiveIndex = mortonPrims[i].primitiveIndex;
			orderedPrims[firstPrimOffset + i] = primitives[primitiveIndex];
			bounds = Union(bounds, primitiveInfo[primitiveIndex].bounds);
			area += primitives[primitiveIndex]->Area();
			count += primitives[primitiveIndex]->Count();
		}
		node->InitLeaf(firstPrimOffset, nPrimitives, count, bounds, area);
		return node;
	} else {
		int mask = 1 << bitIndex;
		// Advance to next subtree level if there's no LBVH split for this bit
		if ((mortonPrims[0].mortonCode & mask)
				== (mortonPrims[nPrimitives - 1].mortonCode & mask))
			return emitLBVH(buildNodes, primitiveInfo, mortonPrims, nPrimitives,
					totalNodes, orderedPrims, orderedPrimsOffset, bitIndex - 1);

		// Find LBVH split point for this dimension
		int searchStart = 0, searchEnd = nPrimitives - 1;
		while (searchStart + 1 != searchEnd) {
			CHECK_NE(searchStart, searchEnd);
			int mid = (searchStart + searchEnd) / 2;
			if ((mortonPrims[searchStart].mortonCode & mask)
					== (mortonPrims[mid].mortonCode & mask))
				searchStart = mid;
			else {
				CHECK_EQ(mortonPrims[mid].mortonCode & mask,
						mortonPrims[searchEnd].mortonCode & mask);
				searchEnd = mid;
			}
		}
		int splitOffset = searchEnd;
		CHECK_LE(splitOffset, nPrimitives - 1);
		CHECK_NE(mortonPrims[splitOffset - 1].mortonCode & mask,
				mortonPrims[splitOffset].mortonCode & mask);

		// Create and return interior LBVH node
		(*totalNodes)++;
		BVHBuildNode *node = buildNodes++;
		BVHBuildNode *lbvh[2] = { emitLBVH(buildNodes, primitiveInfo,
				mortonPrims, splitOffset, totalNodes, orderedPrims,
				orderedPrimsOffset, bitIndex - 1), emitLBVH(buildNodes,
				primitiveInfo, &mortonPrims[splitOffset],
				nPrimitives - splitOffset, totalNodes, orderedPrims,
				orderedPrimsOffset, bitIndex - 1) };
		int axis = bitIndex % 3;
		node->InitInterior(axis, lbvh[0], lbvh[1]);
		return node;
	}
}

BVHBuildNode *BVHAccel::buildUpperSAH(MemoryArena &arena,
		std::vector<BVHBuildNode *> &treeletRoots, int start, int end,
		int *totalNodes) const {
	CHECK_LT(start, end);
	int nNodes = end - start;
	if (nNodes == 1)
		return treeletRoots[start];
	(*totalNodes)++;
	BVHBuildNode *node = arena.Alloc<BVHBuildNode>();

	// Compute bounds of all nodes under this HLBVH node
	Bounds3f bounds;
	for (int i = start; i < end; ++i)
		bounds = Union(bounds, treeletRoots[i]->bounds);

	// Compute bound of HLBVH node centroids, choose split dimension _dim_
	Bounds3f centroidBounds;
	for (int i = start; i < end; ++i) {
		Point3f centroid = (treeletRoots[i]->bounds.pMin
				+ treeletRoots[i]->bounds.pMax) * 0.5f;
		centroidBounds = Union(centroidBounds, centroid);
	}
	int dim = centroidBounds.MaximumExtent();
	// FIXME: if this hits, what do we need to do?
	// Make sure the SAH split below does something... ?
	CHECK_NE(centroidBounds.pMax[dim], centroidBounds.pMin[dim]);

	// Allocate _BucketInfo_ for SAH partition buckets
	PBRT_CONSTEXPR int nBuckets = 12;
	struct BucketInfo {
		int count = 0;
		Bounds3f bounds;
	};
	BucketInfo buckets[nBuckets];

	// Initialize _BucketInfo_ for HLBVH SAH partition buckets
	for (int i = start; i < end; ++i) {
		Float centroid = (treeletRoots[i]->bounds.pMin[dim]
				+ treeletRoots[i]->bounds.pMax[dim]) * 0.5f;
		int b =
				nBuckets
						* ((centroid - centroidBounds.pMin[dim])
								/ (centroidBounds.pMax[dim]
										- centroidBounds.pMin[dim]));
		if (b == nBuckets)
			b = nBuckets - 1;
		CHECK_GE(b, 0);
		CHECK_LT(b, nBuckets);
		buckets[b].count++;
		buckets[b].bounds = Union(buckets[b].bounds, treeletRoots[i]->bounds);
	}

	// Compute costs for splitting after each bucket
	Float cost[nBuckets - 1];
	for (int i = 0; i < nBuckets - 1; ++i) {
		Bounds3f b0, b1;
		int count0 = 0, count1 = 0;
		for (int j = 0; j <= i; ++j) {
			b0 = Union(b0, buckets[j].bounds);
			count0 += buckets[j].count;
		}
		for (int j = i + 1; j < nBuckets; ++j) {
			b1 = Union(b1, buckets[j].bounds);
			count1 += buckets[j].count;
		}
		cost[i] = .125f
				+ (count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea())
						/ bounds.SurfaceArea();
	}

	// Find bucket to split at that minimizes SAH metric
	Float minCost = cost[0];
	int minCostSplitBucket = 0;
	for (int i = 1; i < nBuckets - 1; ++i) {
		if (cost[i] < minCost) {
			minCost = cost[i];
			minCostSplitBucket = i;
		}
	}

	// Split nodes and create interior HLBVH SAH node
	BVHBuildNode **pmid =
			std::partition(&treeletRoots[start], &treeletRoots[end - 1] + 1,
					[=](const BVHBuildNode *node) {
						Float centroid =
						(node->bounds.pMin[dim] + node->bounds.pMax[dim]) * 0.5f;
						int b = nBuckets *
						((centroid - centroidBounds.pMin[dim]) /
								(centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
						if (b == nBuckets) b = nBuckets - 1;
						CHECK_GE(b, 0);
						CHECK_LT(b, nBuckets);
						return b <= minCostSplitBucket;
					});
	int mid = pmid - &treeletRoots[0];
	CHECK_GT(mid, start);
	CHECK_LT(mid, end);
	node->InitInterior(dim,
			this->buildUpperSAH(arena, treeletRoots, start, mid, totalNodes),
			this->buildUpperSAH(arena, treeletRoots, mid, end, totalNodes));
	return node;
}

int BVHAccel::flattenBVHTree(BVHBuildNode *node, int *offset) {
	LinearBVHNode *linearNode = &nodes[*offset];
	linearNode->bounds = node->bounds;
	linearNode->count = node->count;
	linearNode->area = node->area;
	linearNode->depth = node->Depth();

	int myOffset = (*offset)++;
	if (node->nPrimitives > 0) {
		CHECK(!node->children[0] && !node->children[1]);
		CHECK_LT(node->nPrimitives, 65536);
		linearNode->primitivesOffset = node->firstPrimOffset;
		linearNode->nPrimitives = node->nPrimitives;
	} else {
		// Create interior flattened BVH node
		linearNode->axis = node->splitAxis;
		linearNode->nPrimitives = 0;
		flattenBVHTree(node->children[0], offset);
		linearNode->secondChildOffset = flattenBVHTree(node->children[1],
				offset);
	}
	return myOffset;
}

BVHAccel::~BVHAccel() {
	FreeAligned(nodes);
}

bool BVHAccel::IntersectFrom(const Ray &ray, SurfaceInteraction *isect,
		uint32_t from) const {
	if (!nodes)
		return false;
	ProfilePhase p(Prof::AccelIntersect);
	bool hit = false;
	Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
	int dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
	// Follow ray through BVH nodes to find primitive intersections
	int toVisitOffset = 0, currentNodeIndex = from;
	int nodesToVisit[64];

	int64_t localHits = 0, localTests = 0;
	while (true) {
		const LinearBVHNode *node = &nodes[currentNodeIndex];
		++localTests;

		// Check ray against BVH node
		if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
			++localHits;

			if (node->nPrimitives > 0) {
				// Intersect ray with primitives in leaf BVH node
				for (int i = 0; i < node->nPrimitives; ++i)
					if (primitives[node->primitivesOffset + i]->Intersect(ray,
							isect))
						hit = true;
				if (toVisitOffset == 0)
					break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			} else {
				// Put far BVH node on _nodesToVisit_ stack, advance to near
				// node
				if (dirIsNeg[node->axis]) {
					nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
					currentNodeIndex = node->secondChildOffset;
				} else {
					nodesToVisit[toVisitOffset++] = node->secondChildOffset;
					currentNodeIndex = currentNodeIndex + 1;
				}
			}
		} else {
			if (toVisitOffset == 0)
				break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
		}
	}

	nTests += localTests;
	nHits += localHits;
	return hit;
}

bool BVHAccel::IntersectPFrom(const Ray &ray, uint32_t from) const {
	if (!nodes)
		return false;
	ProfilePhase p(Prof::AccelIntersectP);
	Vector3f invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
	int dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
	int nodesToVisit[64];
	int toVisitOffset = 0, currentNodeIndex = from;
	int64_t localHits = 0, localTests = 0;

	while (true) {
		const LinearBVHNode *node = &nodes[currentNodeIndex];
		++localTests;
		++ray.boxIsectPTests;
		if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
			++localHits;

			// Process BVH node _node_ for traversal
			if (node->nPrimitives > 0) {
				for (int i = 0; i < node->nPrimitives; ++i) {
					if (primitives[node->primitivesOffset + i]->IntersectP(
							ray)) {
						nShadowTests += localTests;
						nShadowHits += localHits;
						return true;
					}
				}
				if (toVisitOffset == 0)
					break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			} else {
				if (dirIsNeg[node->axis]) {
					/// second child first
					nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
					currentNodeIndex = node->secondChildOffset;
				} else {
					nodesToVisit[toVisitOffset++] = node->secondChildOffset;
					currentNodeIndex = currentNodeIndex + 1;
				}
			}
		} else {
			if (toVisitOffset == 0)
				break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
		}
	}

	nShadowTests += localTests;
	nShadowHits += localHits;

	return false;
}

bool BVHAccel::IntersectPFrom(const Ray &ray, uint32_t ignoreId,
		uint32_t from) const {
	if (!nodes)
		return false;
	else if (ignoreId == 0)
		return IntersectPFrom(ray, from);
	else if (id == ignoreId)
		return false;

	ProfilePhase p(Prof::AccelIntersectP);
	Vector3f invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
	int dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
	int nodesToVisit[64];
	int toVisitOffset = 0, currentNodeIndex = from;
	int64_t localHits = 0, localTests = 0;

	while (true) {
		const LinearBVHNode *node = &nodes[currentNodeIndex];

		// ignore subtrees!
		if (subtrees[currentNodeIndex]->id == ignoreId) {
			if (toVisitOffset == 0)
				break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
			continue;
		}

		++localTests;
		if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
			++localHits;
			++ray.boxIsectPTests;

			// Process BVH node _node_ for traversal
			if (node->nPrimitives > 0) {
				for (int i = 0; i < node->nPrimitives; ++i) {
					if (primitives[node->primitivesOffset + i]->IntersectP(ray,
							ignoreId)) {
						nShadowTests += localTests;
						nShadowHits += localHits;
						return true;
					}
				}
				if (toVisitOffset == 0)
					break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			} else {
				if (dirIsNeg[node->axis]) {
					/// second child first
					nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
					currentNodeIndex = node->secondChildOffset;
				} else {
					nodesToVisit[toVisitOffset++] = node->secondChildOffset;
					currentNodeIndex = currentNodeIndex + 1;
				}
			}
		} else {
			if (toVisitOffset == 0)
				break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
		}
	}

	nShadowTests += localTests;
	nShadowHits += localHits;

	return false;
}

bool BVHAccel::IntersectLFrom(const LineVisibilityQuery3f &query,
		const Normal3f &normal, Float time, IntervalCollection& collection,
		uint32_t from) const {
	if (!nodes)
		return false;
	ProfilePhase p(Prof::AccelIntersectL);

	const Vector3f& d = 0.5 * (query.p1 + query.p2) - query.p0;
	const int dirIsNeg[3] = { d.x < 0, d.y < 0, d.z < 0 };

	int nodesToVisit[64];
	int toVisitOffset = 0, currentNodeIndex = from;
	bool hit = false;

	uint32_t hits = 0;
	uint32_t tests = 0;
	while (!collection.empty()) {
		const LinearBVHNode *node = &nodes[currentNodeIndex];

		++tests;
		if (IntersectTriangleBox(node->bounds, query)) {
			++hits;
			if (node->nPrimitives > 0) {
				// Intersect ray with primitives in leaf BVH node
				for (int i = 0; i < node->nPrimitives; ++i) {
					const std::shared_ptr<const Primitive>& primitive =
							primitives[node->primitivesOffset + i];

					if (primitive->IntersectL(query, normal, time, collection))
						hit = true;
				}
				if (toVisitOffset == 0)
					break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			} else {
				if (dirIsNeg[node->axis]) {
					// second child first
					nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
					currentNodeIndex = node->secondChildOffset;
				} else {
					nodesToVisit[toVisitOffset++] = node->secondChildOffset;
					currentNodeIndex = currentNodeIndex + 1;
				}
			}
		} else {
			if (toVisitOffset == 0)
				break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
		}
	}

	nTriangleBoxHits += hits;
	nTriangleBoxTests += tests;

	return hit;

}

bool BVHAccel::IntersectAllFrom(const Ray& ray, GroupIsects& intersections,
		int groupSize, uint32_t from) const {
	// if the groupsize is 1, use the version where we do not keep checking
	// the groupsize in each iteration!

	if (groupSize == 1)
		return IntersectAllNoGroupFrom(ray, intersections, from);

	if (!nodes)
		return false;
	ProfilePhase p(Prof::AccelIntersect);
	bool hit = false;
	Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
	int dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
	// Follow ray through BVH nodes to find primitive intersections
	int toVisitOffset = 0, currentNodeIndex = from;
	int nodesToVisit[64];
	while (true) {
		const LinearBVHNode *node = &nodes[currentNodeIndex];

		if (node->count <= groupSize) {
			GroupIsects localIsects;
			if (IntersectAllNoGroupFrom(ray, localIsects, currentNodeIndex)) {
				hit = true;

				for (GroupIsect& localIsect : localIsects)
					localIsect.second = subtrees[currentNodeIndex];

				intersections.insert(intersections.end(), localIsects.begin(),
						localIsects.end());
			}

			if (toVisitOffset == 0)
				break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
		} else {
			// Check ray against BVH node
			if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
				if (node->nPrimitives > 0) {
					// Intersect ray with primitives in leaf BVH node
					for (int i = 0; i < node->nPrimitives; ++i)
						if (primitives[node->primitivesOffset + i]->IntersectAll(
								ray, intersections, groupSize))
							hit = true;
					if (toVisitOffset == 0)
						break;
					currentNodeIndex = nodesToVisit[--toVisitOffset];
				} else {
					// Put far BVH node on _nodesToVisit_ stack, advance to near
					// node
					if (dirIsNeg[node->axis]) {
						nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
						currentNodeIndex = node->secondChildOffset;
					} else {
						nodesToVisit[toVisitOffset++] = node->secondChildOffset;
						currentNodeIndex = currentNodeIndex + 1;
					}
				}
			} else {
				if (toVisitOffset == 0)
					break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			}
		}
	}
	return hit;
}

bool BVHAccel::IntersectAllNoGroupFrom(const Ray& ray,
		GroupIsects& intersections, uint32_t from) const {
	if (!nodes)
		return false;
	ProfilePhase p(Prof::AccelIntersect);
	bool hit = false;
	Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
	int dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
	// Follow ray through BVH nodes to find primitive intersections
	int toVisitOffset = 0, currentNodeIndex = from;
	int nodesToVisit[64];
	while (true) {
		const LinearBVHNode *node = &nodes[currentNodeIndex];
		// Check ray against BVH node
		if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
			if (node->nPrimitives > 0) {
				// Intersect ray with primitives in leaf BVH node
				for (int i = 0; i < node->nPrimitives; ++i)
					if (primitives[node->primitivesOffset + i]->IntersectAll(
							ray, intersections, 1))
						hit = true;
				if (toVisitOffset == 0)
					break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			} else {
				// Put far BVH node on _nodesToVisit_ stack, advance to near
				// node
				if (dirIsNeg[node->axis]) {
					nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
					currentNodeIndex = node->secondChildOffset;
				} else {
					nodesToVisit[toVisitOffset++] = node->secondChildOffset;
					currentNodeIndex = currentNodeIndex + 1;
				}
			}
		} else {
			if (toVisitOffset == 0)
				break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
		}
	}
	return hit;
}

Bounds3f BVHSubtree::WorldBound() const {
	return bvh->nodes[offset].bounds;
}

uint32_t BVHSubtree::Count() const {
	return bvh->nodes[offset].count;
}

Float BVHSubtree::Area() const {
	return bvh->nodes[offset].area;
}

bool BVHSubtree::Intersect(const Ray& ray, SurfaceInteraction *isect) const {
	return bvh->IntersectFrom(ray, isect, offset);
}

bool BVHSubtree::IntersectP(const Ray& ray) const {
	return bvh->IntersectPFrom(ray, offset);
}

bool BVHSubtree::IntersectP(const Ray& ray, uint32_t ignoreId) const {
	if (ignoreId == 0)
		return IntersectP(ray);
	else if (id == ignoreId)
		return false;
	return bvh->IntersectPFrom(ray, ignoreId, offset);
}

std::shared_ptr<const BVHAccel> CreateBVHAccelerator(
		const std::vector<std::shared_ptr<const Primitive>> &prims,
		const ParamSet &ps) {
	std::string splitMethodName = ps.FindOneString("splitmethod", "sah");
	BVHAccel::SplitMethod splitMethod;
	if (splitMethodName == "sah")
		splitMethod = BVHAccel::SplitMethod::SAH;
	else if (splitMethodName == "hlbvh")
		splitMethod = BVHAccel::SplitMethod::HLBVH;
	else if (splitMethodName == "middle")
		splitMethod = BVHAccel::SplitMethod::Middle;
	else if (splitMethodName == "equal")
		splitMethod = BVHAccel::SplitMethod::EqualCounts;
	else {
		Warning("BVH split method \"%s\" unknown.  Using \"sah\".",
				splitMethodName.c_str());
		splitMethod = BVHAccel::SplitMethod::SAH;
	}

	int maxPrimsInNode = ps.FindOneInt("maxnodeprims", 1);
	return std::make_shared<const BVHAccel>(prims, maxPrimsInNode, splitMethod);
}

}  // namespace pbrt
