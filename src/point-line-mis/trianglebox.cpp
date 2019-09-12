/*
 * intersections.cpp
 *
 *  Created on: Jan 26, 2016
 *      Author: niels
 */

#include "trianglebox.h"

#include "core/geometry.h"
#include "core/pbrt.h"
#include "point-line-mis/intersections.h"
#include "point-line-mis/interval.h"
#include "point-line-mis/line3.h"
#include "point-line-mis/triangle3.h"

using pbrt::Float;

namespace linesampling {

inline bool IntersectTrianglePlane(int axis, Float position,
		const Bounds3f &box, const LineVisibilityQuery3f &query,
		Interval &interval) {
	// offset of plane u to origin;
	Float dv0 = query.p0[axis] - position; // distance of v0 to plane of triangle u
	Float dv1 = query.p1[axis] - position; // distance of v1 to plane of triangle u
	Float dv2 = query.p2[axis] - position; // distance of v2 to plane of triangle u

	// test whether they have the same sign and different from 0.0
	const Float dv1dv0 = dv1 * dv0;
	const Float dv1dv2 = dv1 * dv2;
	if (dv1dv0 > 0 && dv1dv2 > 0)
		return false;

	// find first intersection point
	Float p0_a2;
	Float p0_a3;
	Float p1_a2;
	Float p1_a3;
	int a1 = axis;
	int a2 = (axis + 1) % 3;
	int a3 = (axis + 2) % 3;
	if (dv1dv2 > 0) {
		// v1 v2 same side
		const Float s = dv1 / (dv1 - dv0);
		const Float t = dv2 / (dv2 - dv0);
		p0_a2 = query.p1[a2] * (1.0 - s) + s * query.p0[a2];
		p0_a3 = query.p1[a3] * (1.0 - s) + s * query.p0[a3];
		p1_a2 = query.p2[a2] * (1.0 - t) + t * query.p0[a2];
		p1_a3 = query.p2[a3] * (1.0 - t) + t * query.p0[a3];
	} else if (dv1dv0 > 0) {
		// v0 v1 same side
		const Float s = dv0 / (dv0 - dv2);
		const Float t = dv1 / (dv1 - dv2);
		p0_a2 = query.p0[a2] * (1.0 - s) + s * query.p2[a2];
		p0_a3 = query.p0[a3] * (1.0 - s) + s * query.p2[a3];
		p1_a2 = query.p1[a2] * (1.0 - t) + t * query.p2[a2];
		p1_a3 = query.p1[a3] * (1.0 - t) + t * query.p2[a3];

	} else {
		// v0 v2 same side
		const Float s = dv0 / (dv0 - dv1);
		const Float t = dv2 / (dv2 - dv1);
		p0_a2 = query.p0[a2] * (1.0 - s) + s * query.p1[a2];
		p0_a3 = query.p0[a3] * (1.0 - s) + s * query.p1[a3];
		p1_a2 = query.p2[a2] * (1.0 - t) + t * query.p1[a2];
		p1_a3 = query.p2[a3] * (1.0 - t) + t * query.p1[a3];
	}

	const Float dy = p1_a2 - p0_a2;
	const Float dz = p1_a3 - p0_a3;
	const Float invDy = 1.0 / dy;
	const Float invDz = 1.0 / dz;
	const int yIsNeg = dy < 0 ? 1 : 0;
	const int zIsNeg = dz < 0 ? 1 : 0;

	const Float tMinY = (box[yIsNeg][a2] - p0_a2) * invDy;
	const Float tMaxY = (box[1 - yIsNeg][a2] - p0_a2) * invDy;
	const Float tMinZ = (box[zIsNeg][a3] - p0_a3) * invDz;
	const Float tMaxZ = (box[1 - zIsNeg][a3] - p0_a3) * invDz;

	if (tMinY > tMaxZ || tMinZ > tMaxY)
		return false;

	Float min = std::max(tMinY, tMinZ);
	Float max = std::min(tMaxY, tMaxZ);

//	if (dv1dv2 > 0) {
//		if (min < 0 && max > 1.0) {
//			interval.min = 0.0;
//			interval.max = 1.0;
//			return true;
//		}
//	}

	if (min > 1.0 || max < 0.0)
		return false;

	min = std::max(0.0, min);
	max = std::min(1.0, max);

	Point3f isect0, isect1;
	isect0[a1] = position;
	isect0[a2] = p0_a2 + dy * min;
	isect0[a3] = p0_a3 + dz * min;
	isect1[a1] = position;
	isect1[a2] = p0_a2 + dy * max;
	isect1[a3] = p0_a3 + dz * max;

	interval = FindShadowedInterval(query.p0, query.p1, query.p2, isect0,
			isect1);

	return true;
}

inline bool IntersectTrianglePlaneX(Float position, const Bounds3f &box,
		const LineVisibilityQuery3f &query, Interval &interval) {
	return IntersectTrianglePlane(0, position, box, query, interval);
//	// offset of plane u to origin;
//	const Float dv0 = query.p0.x - position; // distance of v0 to plane of triangle u
//	const Float dv1 = query.p1.x - position; // distance of v1 to plane of triangle u
//	const Float dv2 = query.p2.x - position; // distance of v2 to plane of triangle u
//
//	// test whether they have the same sign and different from 0.0
//	const Float dv1dv0 = dv1 * dv0;
//	const Float dv1dv2 = dv1 * dv2;
//	if (dv1dv0 > 0 && dv1dv2 > 0)
//		return false;
//
//	// find first intersection point
//	Float p0y;
//	Float p0z;
//	Float p1y;
//	Float p1z;
//	if (dv1dv2 > 0) {
//		// v1 v2 same side
//		const Float s = dv1 / (dv1 - dv0);
//		const Float t = dv2 / (dv2 - dv0);
//		p0y = query.p1.y * (1.0 - s) + s * query.p0.y;
//		p0z = query.p1.z * (1.0 - s) + s * query.p0.z;
//		p1y = query.p2.y * (1.0 - t) + t * query.p0.y;
//		p1z = query.p2.z * (1.0 - t) + t * query.p0.z;
//	} else if (dv1dv0 > 0) {
//		// v0 v1 same side
//		const Float s = dv0 / (dv0 - dv2);
//		const Float t = dv1 / (dv1 - dv2);
//		p0y = query.p0.y * (1.0 - s) + s * query.p2.y;
//		p0z = query.p0.z * (1.0 - s) + s * query.p2.z;
//		p1y = query.p1.y * (1.0 - t) + t * query.p2.y;
//		p1z = query.p1.z * (1.0 - t) + t * query.p2.z;
//
//	} else {
//		// v0 v2 same side
//		const Float s = dv0 / (dv0 - dv1);
//		const Float t = dv2 / (dv2 - dv1);
//		p0y = query.p0.y * (1.0 - s) + s * query.p1.y;
//		p0z = query.p0.z * (1.0 - s) + s * query.p1.z;
//		p1y = query.p2.y * (1.0 - t) + t * query.p1.y;
//		p1z = query.p2.z * (1.0 - t) + t * query.p1.z;
//	}
//
//	const Float dy = p1y - p0y;
//	const Float dz = p1z - p0z;
//	const Float invDy = 1.0 / dy;
//	const Float invDz = 1.0 / dz;
//	const int yIsNeg = dy < 0 ? 1 : 0;
//	const int zIsNeg = dz < 0 ? 1 : 0;
//
//	const Float tMinY = (box[yIsNeg].y - p0y) * invDy;
//	const Float tMaxY = (box[1 - yIsNeg].y - p0y) * invDy;
//	const Float tMinZ = (box[zIsNeg].z - p0z) * invDz;
//	const Float tMaxZ = (box[1 - zIsNeg].z - p0z) * invDz;
//
//	if (tMinY > tMaxZ || tMinZ > tMaxY)
//		return false;
//
//	Float min = std::max(tMinY, tMinZ);
//	Float max = std::min(tMaxY, tMaxZ);
//
//	if (dv1dv2 > 0) {
//		if (min < 0 && max > 1) {
//			interval.min = 0;
//			interval.max = 1;
//			return true;
//		}
//	}
//
//	if (min > 1.0 || max < 0.0)
//		return false;
//
//	min = std::max(0.0, min);
//	max = std::min(1.0, max);
//
//	const Point3f isect0(position, p0y + dy * min, p0z + dz * min);
//	const Point3f isect1(position, p0y + dy * max, p0z + dz * max);
//
//	interval = FindShadowedInterval(query.p0, query.p1, query.p2, isect0,
//			isect1);
//
//	return true;
}

inline bool IntersectTrianglePlaneY(Float position, const Bounds3f &box,
		const LineVisibilityQuery3f &query, Interval &interval) {
	return IntersectTrianglePlane(1, position, box, query, interval);

//	// offset of plane u to origin;
//	const Float dv0 = query.p0.y - position; // distance of v0 to plane of triangle u
//	const Float dv1 = query.p1.y - position; // distance of v1 to plane of triangle u
//	const Float dv2 = query.p2.y - position; // distance of v2 to plane of triangle u
//
//	// test whether they have the same sign and different from 0.0
//	const Float dv1dv0 = dv1 * dv0;
//	const Float dv1dv2 = dv1 * dv2;
//	if (dv1dv0 > 0 && dv1dv2 > 0)
//		return false;
//
//	// find first intersection point
//	Float p0x;
//	Float p0z;
//	Float p1x;
//	Float p1z;
//	if (dv1dv2 > 0) {
//		// v1 v2 same side
//		const Float s = (dv1) / (dv1 - dv0);
//		const Float t = (dv2) / (dv2 - dv0);
//		p0x = query.p1.x * (1.0 - s) + s * query.p0.x;
//		p0z = query.p1.z * (1.0 - s) + s * query.p0.z;
//		p1x = query.p2.x * (1.0 - t) + t * query.p0.x;
//		p1z = query.p2.z * (1.0 - t) + t * query.p0.z;
//	} else if (dv1dv0 > 0) {
//		// v0 v1 same side
//		const Float s = (dv0) / (dv0 - dv2);
//		const Float t = (dv1) / (dv1 - dv2);
//		p0x = query.p0.x * (1.0 - s) + s * query.p2.x;
//		p0z = query.p0.z * (1.0 - s) + s * query.p2.z;
//		p1x = query.p1.x * (1.0 - t) + t * query.p2.x;
//		p1z = query.p1.z * (1.0 - t) + t * query.p2.z;
//
//	} else {
//		// v0 v2 same side
//		const Float s = (dv0) / (dv0 - dv1);
//		const Float t = (dv2) / (dv2 - dv1);
//		p0x = query.p0.x * (1.0 - s) + s * query.p1.x;
//		p0z = query.p0.z * (1.0 - s) + s * query.p1.z;
//		p1x = query.p2.x * (1.0 - t) + t * query.p1.x;
//		p1z = query.p2.z * (1.0 - t) + t * query.p1.z;
//	}
//
//	const Float dx = p1x - p0x;
//	const Float dz = p1z - p0z;
//	const Float invDx = 1.0 / dx;
//	const Float invDz = 1.0 / dz;
//	const int xIsNeg = dx < 0 ? 1 : 0;
//	const int zIsNeg = dz < 0 ? 1 : 0;
//
//	const Float tMinX = (box[xIsNeg].x - p0x) * invDx;
//	const Float tMaxX = (box[1 - xIsNeg].x - p0x) * invDx;
//	const Float tMinZ = (box[zIsNeg].z - p0z) * invDz;
//	const Float tMaxZ = (box[1 - zIsNeg].z - p0z) * invDz;
//
//	if (tMinX > tMaxZ || tMinZ > tMaxX)
//		return false;
//
//	Float min = std::max(tMinX, tMinZ);
//	Float max = std::min(tMaxX, tMaxZ);
//
//	if (dv1dv2 > 0) {
//		if (min < 0 && max > 1) {
//			interval.min = 0;
//			interval.max = 1;
//			return true;
//		}
//	}
//
//	if (min > 1.0 || max < 0.0)
//		return false;
//
//	min = std::max(0.0, min);
//	max = std::min(1.0, max);
//
//	const Point3f isect0(p0x + dx * min, position, p0z + dz * min);
//	const Point3f isect1(p0x + dx * max, position, p0z + dz * max);
//
//	interval = FindShadowedInterval(query.p0, query.p1, query.p2, isect0,
//			isect1);
//
//	return true;
}

inline bool IntersectTrianglePlaneZ(Float position, const Bounds3f &box,
		const LineVisibilityQuery3f &query, Interval &interval) {
//	// offset of plane u to origin;
//	const Float dv0 = query.p0.z - position; // distance of v0 to plane of triangle u
//	const Float dv1 = query.p1.z - position; // distance of v1 to plane of triangle u
//	const Float dv2 = query.p2.z - position; // distance of v2 to plane of triangle u
//
//	// test whether they have the same sign and different from 0.0
//	const Float dv1dv0 = dv1 * dv0;
//	const Float dv1dv2 = dv1 * dv2;
//	if (dv1dv0 > 0 && dv1dv2 > 0)
//		return false;
//
//	// find first intersection point
//	Float p0x;
//	Float p0y;
//	Float p1x;
//	Float p1y;
//	if (dv1dv2 > 0) {
//		// v1 v2 same side
//		const Float s = (dv1) / (dv1 - dv0);
//		const Float t = (dv2) / (dv2 - dv0);
//		p0x = query.p1.x * (1.0 - s) + s * query.p0.x;
//		p0y = query.p1.y * (1.0 - s) + s * query.p0.y;
//		p1x = query.p2.x * (1.0 - t) + t * query.p0.x;
//		p1y = query.p2.y * (1.0 - t) + t * query.p0.y;
//	} else if (dv1dv0 > 0) {
//		// v0 v1 same side
//		const Float s = (dv0) / (dv0 - dv2);
//		const Float t = (dv1) / (dv1 - dv2);
//		p0x = query.p0.x * (1.0 - s) + s * query.p2.x;
//		p0y = query.p0.y * (1.0 - s) + s * query.p2.y;
//		p1x = query.p1.x * (1.0 - t) + t * query.p2.x;
//		p1y = query.p1.y * (1.0 - t) + t * query.p2.y;
//
//	} else {
//		// v0 v2 same side
//		const Float s = (dv0) / (dv0 - dv1);
//		const Float t = (dv2) / (dv2 - dv1);
//		p0x = query.p0.x * (1.0 - s) + s * query.p1.x;
//		p0y = query.p0.y * (1.0 - s) + s * query.p1.y;
//		p1x = query.p2.x * (1.0 - t) + t * query.p1.x;
//		p1y = query.p2.y * (1.0 - t) + t * query.p1.y;
//	}
//
//	const Float dx = p1x - p0x;
//	const Float dy = p1y - p0y;
//	const Float invDx = 1.0 / dx;
//	const Float invDy = 1.0 / dy;
//	const int xIsNeg = dx < 0 ? 1 : 0;
//	const int yIsNeg = dy < 0 ? 1 : 0;
//
//	const Float tMinX = (box[xIsNeg].x - p0x) * invDx;
//	const Float tMaxX = (box[1 - xIsNeg].x - p0x) * invDx;
//	const Float tMinY = (box[yIsNeg].y - p0y) * invDy;
//	const Float tMaxY = (box[1 - yIsNeg].y - p0y) * invDy;
//
//	if (tMinX > tMaxY || tMinY > tMaxX)
//		return false;
//
//	Float min = std::max(tMinX, tMinY);
//	Float max = std::min(tMaxX, tMaxY);
//
//	if (dv1dv2 > 0) {
//		if (min < 0 && max > 1) {
//			interval.min = 0;
//			interval.max = 1;
//			return true;
//		}
//	}
//
//	if (min > 1.0 || max < 0.0)
//		return false;
//
//	min = std::max(0.0, min);
//	max = std::min(1.0, max);
//
//	const Point3f isect0(p0x + dx * min, p0y + dy * min, position);
//	const Point3f isect1(p0x + dx * max, p0y + dy * max, position);
//
//	interval = FindShadowedInterval(query.p0, query.p1, query.p2, isect0,
//			isect1);
//
//	return true;
	return IntersectTrianglePlane(2, position, box, query, interval);

}

bool IntersectTriangleBoxBounds(const Bounds3f& bounds,
		const LineVisibilityQuery3f& triangle, Interval &interval) {
	static const Vector3f e[3] = { Vector3f(1, 0, 0), Vector3f(0, 1, 0),
			Vector3f(0, 0, 1) };

	// quick accept
	if (bounds.isInside(triangle.p0)) {
		interval.min = 0;
		interval.max = 1;
		return true;
	}

	if (bounds.isInside(triangle.p1) && bounds.isInside(triangle.p2)) {
		interval.min = 0;
		interval.max = 1;
		return true;
	}

	// false when the bounding boxes do not overlap
	if (!Overlaps(bounds, triangle.bounds))
		return false;

	// false when the plane of the triangle does not intersect the box
	if (!IntersectPlaneBox(triangle.n, triangle.p0, bounds))
		return false;

	/*************************************************
	 * Check if the x axis test is correct
	 ************************************************/

	bool result = false;

	Interval t;
	if (IntersectTrianglePlaneX(bounds.pMin.x, bounds, triangle, t)) {
		interval.min = std::min(interval.min, t.min);
		interval.max = std::max(interval.max, t.max);
		result = true;
	}
	if (IntersectTrianglePlaneX(bounds.pMax.x, bounds, triangle, t)) {
		interval.min = std::min(interval.min, t.min);
		interval.max = std::max(interval.max, t.max);
		result = true;
	}

	if (IntersectTrianglePlaneY(bounds.pMin.y, bounds, triangle, t)) {
		interval.min = std::min(interval.min, t.min);
		interval.max = std::max(interval.max, t.max);
		result = true;
	}
	if (IntersectTrianglePlaneY(bounds.pMax.y, bounds, triangle, t)) {
		interval.min = std::min(interval.min, t.min);
		interval.max = std::max(interval.max, t.max);
		result = true;
	}

	if (IntersectTrianglePlaneZ(bounds.pMin.z, bounds, triangle, t)) {
		interval.min = std::min(interval.min, t.min);
		interval.max = std::max(interval.max, t.max);
		result = true;
	}
	if (IntersectTrianglePlaneZ(bounds.pMax.z, bounds, triangle, t)) {
		interval.min = std::min(interval.min, t.min);
		interval.max = std::max(interval.max, t.max);
		result = true;
	}

	return result;
}

bool IntersectTriangleBoxBoundsOptimized(const Bounds3f& bounds,
		const LineVisibilityQuery3f& triangle, Interval &interval) {
	static const Vector3f e[3] = { Vector3f(1, 0, 0), Vector3f(0, 1, 0),
			Vector3f(0, 0, 1) };

	// quick accept
	if (bounds.isInside(triangle.p0)) {
		interval.min = 0;
		interval.max = 1;
		return true;
	}

	if (bounds.isInside(triangle.p1) && bounds.isInside(triangle.p2)) {
		interval.min = 0;
		interval.max = 1;
		return true;
	}

	// false when the bounding boxes do not overlap
	if (!Overlaps(bounds, triangle.bounds))
		return false;

	// false when the plane of the triangle does not intersect the box
	if (!IntersectPlaneBox(triangle.n, triangle.p0, bounds))
		return false;

	/******************************************************************************
	 * Check the position of the origin to determine the planes to be intersected *
	 ******************************************************************************/

	bool result = false;
	Interval t;

	if (triangle.p0.x < bounds.pMin.x) {
		// left side of the box
		if (triangle.p0.y < bounds.pMin.y) {
			// below the box
			if (triangle.p0.z < bounds.pMin.z) {
				// behind the box

				/* CAN ONLY SEE:
				 * MIN_X_PLANE
				 * MIN_Y PLANE
				 * MIN_Z_PLANE
				 */
				if (IntersectTrianglePlaneX(bounds.pMin.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMin.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMin.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else if (triangle.p0.z > bounds.pMax.z) {
				// in front of the box

				/* CAN ONLY SEE:
				 *
				 * MIN_X_PLANE
				 * MIN_Y_PLANE
				 * MAX_Z_PLANE
				 */
				if (IntersectTrianglePlaneX(bounds.pMin.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMin.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMax.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else {
				// middle of box

				/* CAN ONLY SEE:
				 *
				 * MIN_X_PLANE
				 * MIN_Y_PLANE
				 */
				if (IntersectTrianglePlaneX(bounds.pMin.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMin.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			}
		} else if (triangle.p0.y > bounds.pMax.y) {
			// on top of the box
			if (triangle.p0.z < bounds.pMin.z) {
				// behind the box

				/*
				 * CAN ONLY SEE
				 *
				 * MIN_X_PLANE
				 * MAX_Y_PLANE
				 * MIN_Z_PLANE
				 */
				if (IntersectTrianglePlaneX(bounds.pMin.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMax.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMin.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else if (triangle.p0.z > bounds.pMax.z) {
				// in front of box

				/*
				 * CAN ONLY SEE
				 *
				 * MIN_X_PLANE
				 * MAX_Y_PLANE
				 * MAX_Z_PLANE
				 */

				if (IntersectTrianglePlaneX(bounds.pMin.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMax.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMax.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else {
				// in center of box

				/*
				 * CAN ONLY SEE
				 *
				 * MIN_X_PLANE
				 * MAX_Y_PLANE
				 */
				if (IntersectTrianglePlaneX(bounds.pMin.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMax.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			}
		} else {
			// in the middle of y
			if (triangle.p0.z < bounds.pMin.z) {
				// behind the box

				/*
				 * CAN ONLY SEE:
				 *
				 * MIN_X_PLANE
				 * MIN_Z_PLANE
				 */
				if (IntersectTrianglePlaneX(bounds.pMin.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMin.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else if (triangle.p0.z > bounds.pMax.z) {
				// in front of the box

				/*
				 * CAN ONLY SEE:
				 *
				 * MIN_X_PLANE
				 * MAX_Z_PLANE
				 */
				if (IntersectTrianglePlaneX(bounds.pMin.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMax.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else {
				/*
				 * CAN ONLY SEE:
				 *
				 * MIN_X_PLANE
				 *
				 * NEEDS FURTHER INSPECTION
				 */
				// TODO
				if (IntersectTrianglePlaneX(bounds.pMin.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			}
		}
	} else if (triangle.p0.x > bounds.pMax.x) {
		// right side of the box
		if (triangle.p0.y < bounds.pMin.y) {
			// below the box
			if (triangle.p0.z < bounds.pMin.z) {
				// behind the cube

				/*
				 * CAN ONLY SEE
				 *
				 * MAX_X_PLANE
				 * MIN_Y_PLANE
				 * MIN_Z_PLANE
				 */
				if (IntersectTrianglePlaneX(bounds.pMax.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMin.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMin.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else if (triangle.p0.z > bounds.pMax.z) {
				// in front of the cube

				/*
				 * CAN ONLY SEE
				 *
				 * MAX_X_PLANE
				 * MIN_Y_PLANE
				 * MAX_Z_PLANE
				 */
				if (IntersectTrianglePlaneX(bounds.pMax.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMin.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMax.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else {
				// in middle along z

				/*
				 * can only see
				 *
				 * max_x_plane
				 * min_y_plane
				 */
				if (IntersectTrianglePlaneX(bounds.pMax.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMin.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			}
		} else if (triangle.p0.y > bounds.pMax.y) {
			// above the box
			if (triangle.p0.z < bounds.pMin.z) {
				// behind the cube

				/*
				 * can only see
				 *
				 * max_x_plane
				 * max_y_plane
				 * min_z_plane
				 */
				if (IntersectTrianglePlaneX(bounds.pMax.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMax.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMin.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else if (triangle.p0.z > bounds.pMax.z) {
				// behind the cube

				/*
				 * can only see
				 *
				 * max_x_plane
				 * max_y_plane
				 * min_z_plane
				 */
				if (IntersectTrianglePlaneX(bounds.pMax.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMax.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMax.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else {
				/*
				 * can only see
				 *
				 * max_x_plane
				 * max_y_plane
				 */
				if (IntersectTrianglePlaneX(bounds.pMax.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMax.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			}
		} else {
			// in middle along y
			if (triangle.p0.z < bounds.pMin.z) {
				// behind the cube

				/*
				 * can only see
				 *
				 * max_x_plane
				 * min_z_plane
				 */
				if (IntersectTrianglePlaneX(bounds.pMax.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMin.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else if (triangle.p0.z > bounds.pMax.z) {
				/*
				 * can only see
				 *
				 * max_x_plane
				 * min_z_plane
				 */
				if (IntersectTrianglePlaneX(bounds.pMax.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMax.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else {
				/*
				 * can only see
				 *
				 * max_x_plane
				 */
				if (IntersectTrianglePlaneX(bounds.pMax.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
				// TODO
				// TODO
				// TODO
			}
		}
	} else {
		// middle along x
		if (triangle.p0.y < bounds.pMin.y) {
			// below the box
			if (triangle.p0.z < bounds.pMin.z) {
				// behind of the box

				/*
				 * can only see
				 *
				 * min_y
				 * min_z
				 */
				if (IntersectTrianglePlaneY(bounds.pMin.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMin.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else if (triangle.p0.z > bounds.pMax.z) {
				// in front of the box
				/*
				 * can only see
				 *
				 * min_y
				 * max_z
				 */
				if (IntersectTrianglePlaneY(bounds.pMin.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMax.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else {
				// in between z

				/*
				 * can only see
				 *
				 * min_y
				 */
				if (IntersectTrianglePlaneY(bounds.pMin.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			}
		} else if (triangle.p0.y > bounds.pMax.y) {
			// below the box
			if (triangle.p0.z < bounds.pMin.z) {
				// behind of the box

				/*
				 * can only see
				 *
				 * max_y
				 * min_z
				 */
				if (IntersectTrianglePlaneY(bounds.pMax.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMin.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else if (triangle.p0.z > bounds.pMax.z) {
				// in front of the box
				/*
				 * can only see
				 *
				 * min_y
				 * max_z
				 */
				if (IntersectTrianglePlaneY(bounds.pMax.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMax.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else {
				// in between z

				/*
				 * can only see
				 *
				 * min_y
				 */
				if (IntersectTrianglePlaneY(bounds.pMax.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}

//				if (IntersectTrianglePlaneY(bounds.pMin.y, bounds, triangle,
//						t)) {
//					if (t.min < interval.min - 1e-5)
//						Severe("should not happen; t.min! %d", result);
//					if (t.max > interval.max + 1e-5)
//						Severe("should not happen; t.max! %d", result);
////					interval.min = std::min(interval.min, t.min);
////					interval.max = std::max(interval.max, t.max);
////					result = true;
//				}

				return result;
				// TODO
				// TODO
//				return result;
			}
		} else {
			// y along the middle
			if (triangle.p0.z < bounds.pMin.z) {
				// behind of the box
				if (IntersectTrianglePlaneZ(bounds.pMin.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else if (triangle.p0.z > bounds.pMax.z) {
				// behind of the box
				if (IntersectTrianglePlaneZ(bounds.pMax.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			} else {
				if (IntersectTrianglePlaneX(bounds.pMin.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneX(bounds.pMax.x, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMin.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneY(bounds.pMax.y, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMin.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				if (IntersectTrianglePlaneZ(bounds.pMax.z, bounds, triangle,
						t)) {
					interval.min = std::min(interval.min, t.min);
					interval.max = std::max(interval.max, t.max);
					result = true;
				}
				return result;
			}
		}
	}
}

}
