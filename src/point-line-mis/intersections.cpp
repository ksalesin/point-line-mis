/*
 * intersections.cpp
 *
 *  Created on: Jan 26, 2016
 *      Author: niels
 */

#include "point-line-mis/intersections.h"

namespace linesampling {

Float ShortestLineParameter(const Point3f &origin, const Point3f& destination,
		const Point3f &x, const Point3f &y) {
	return ShortestLineParameter(origin, destination - origin, x, y);
}

Float ShortestLineParameter(const Point3f &origin, const Vector3f &direction,
		const Point3f &x, const Point3f &y) {
	const Vector3f xy = y - x;
	const Vector3f w = origin - x;

	const Float a = Dot(direction, direction);
	const Float b = Dot(direction, xy);
	const Float c = Dot(xy, xy);
	const Float d = Dot(direction, w);
	const Float e = Dot(xy, w);

	const Float D = a * c - b * b;
	if (D < EPSILON)
		return 0;
	else
		return (b * e - c * d) / D;
}

/*
 * Calculates the intersection between a plane with normal
 * <n> and distance <d> to the origin with the line segment
 * spanned by the two points <p1> and <p2>.
 *
 * @param n			the normal of the plane
 * @param d			the distance of the plane to the origin
 * @param origin 	the origin of the line segment
 * @param direction	the direction of the line segment
 */

bool IntersectPlaneLine(const Normal3f &n, Float d, const Point3f &origin,
		const Vector3f &direction, Point3f *isect, Float *isectParameter) {

	Float numerator = Dot(n, origin) + d;
	Float denominator = Dot(n, direction);

	// line is parallel with the plane.
	if (std::abs(denominator) < EPSILON)
		return false;

	Float t = numerator / denominator;
	if (isectParameter)
		*isectParameter = t;
	if (isect)
		*isect = origin + t * direction;
	return true;
}

bool IntersectTriangleTriangle(const Point3f& u0, const Point3f& u1,
		const Point3f& u2, const Normal3f &un, const Point3f& v0,
		const Point3f& v1, const Point3f& v2, const Normal3f &vn,
		Line3f &isect) {
	const Float du = -Dot(un, u0); // offset of plane u to origin
	Float dv0 = Dot(un, v0) + du; // distance of v0 to plane of triangle u
	Float dv1 = Dot(un, v1) + du; // distance of v1 to plane of triangle u
	Float dv2 = Dot(un, v2) + du; // distance of v2 to plane of triangle u

	// test whether they have the same sign and different from 0.0
	const Float dv0dv1 = dv0 * dv1;
	const Float dv0dv2 = dv0 * dv2;
	if (dv0dv1 > 0 && dv0dv2 > 0)
		return false;

	const Float dv = -Dot(vn, v0); // offset of plane u to origin
	Float du0 = Dot(vn, u0) + dv; // distance of u0 to plane of triangle v
	Float du1 = Dot(vn, u1) + dv; // distance of u1 to plane of triangle v
	Float du2 = Dot(vn, u2) + dv; // distance of u2 to plane of triangle v

	// test whether they have the same sign and different from 0.0
	const Float du0du1 = du0 * du1;
	const Float du0du2 = du0 * du2;
	if (du0du1 > 0 && du0du2 > 0)
		return false;

	Point3f O;
	Vector3f D;
	IntersectPlanePlane(un, du, vn, dv, O, D);

	// get the parametric intersection between the edges of the triangles and the intersection line
	Float ut0, ut1;
	if (du0du1 > 0) {
		// u0 and u1 same side
		ut0 = ShortestLineParameter(O, D, u0, u2);
		ut1 = ShortestLineParameter(O, D, u1, u2);
	} else if (du0du2 > 0) {
		// u0 and u2 same side
		ut0 = ShortestLineParameter(O, D, u0, u1);
		ut1 = ShortestLineParameter(O, D, u2, u1);
	} else {
		// u1 and u2 same side
		ut0 = ShortestLineParameter(O, D, u1, u0);
		ut1 = ShortestLineParameter(O, D, u2, u0);
	}

	Float vt0, vt1;
	if (dv0dv1 > 0) {
		// v0 and v1 same side
		vt0 = ShortestLineParameter(O, D, v0, v2);
		vt1 = ShortestLineParameter(O, D, v1, v2);
	} else if (dv0dv2 > 0) {
		// v0 and v2 same side
		vt0 = ShortestLineParameter(O, D, v0, v1);
		vt1 = ShortestLineParameter(O, D, v2, v1);
	} else {
		// v1 and v2 same side
		vt0 = ShortestLineParameter(O, D, v1, v0);
		vt1 = ShortestLineParameter(O, D, v2, v0);
	}

	// find the intersection
	Float minu = std::min(ut0, ut1);
	Float maxu = std::max(ut0, ut1);
	Float minv = std::min(vt0, vt1);
	Float maxv = std::max(vt0, vt1);

	if (minu > maxv || maxu < minv)
		return false;

	Float min = std::max(minu, minv);
	Float max = std::min(maxu, maxv);

	if (std::abs(max - min) < EPSILON)
		return false;

	isect.p1 = O + D * min;
	isect.p2 = O + D * max;

	return true;
}

bool IntersectTriangleTriangle(const Triangle3f &u, const Triangle3f &v,
		Line3f &isect) {
	return IntersectTriangleTriangle(u.p0, u.p1, u.p2, u.n, v.p0, v.p1, v.p2,
			v.n, isect);
}

Interval FindShadowedInterval(const Point3f& origin, const Point3f& lineOrigin,
		const Point3f& lineDestination, const Point3f& p1, const Point3f&p2) {
	Float t0 = ShortestLineParameter(lineOrigin, lineDestination, origin, p1);
	Float t1 = ShortestLineParameter(lineOrigin, lineDestination, origin, p2);
	if (t0 > t1)
		std::swap(t0, t1);
	return Interval(t0, t1);
}

Interval FindShadowedInterval(const LineVisibilityQuery3f& query,
		const Line3f &isect) {
	Float t0 = ShortestLineParameter(query.p1, query.p2, query.p0, isect.p1);
	Float t1 = ShortestLineParameter(query.p1, query.p2, query.p0, isect.p2);

	if (t0 > t1)
		std::swap(t0, t1);
	return Interval(t0, t1);
}

void IntersectPlanePlane(const Normal3f &n1, Float d1, const Normal3f &n2,
		Float d2, Point3f &origin, Vector3f &direction) {
	direction = Vector3f(Cross(n1, n2));

	// find the most numerically stable axis
	Float adx = std::abs(direction.x);
	Float ady = std::abs(direction.y);
	Float adz = std::abs(direction.z);

	// find any origin for the line intersection
	if (adx > ady && adx > adz) {
		Float invDenominator = 1.0 / direction.x;
		origin.x = 0;
		origin.y = (d2 * n1.z - d1 * n2.z) * invDenominator;
		origin.z = (d1 * n2.y - d2 * n1.y) * invDenominator;
	} else if (ady > adz) {
		Float invDenominator = -1.0 / direction.y;
		origin.x = (d2 * n1.z - d1 * n2.z) * invDenominator;
		origin.y = 0;
		origin.z = (d1 * n2.x - d2 * n1.x) * invDenominator;
	} else {
		Float invDenominator = 1.0 / direction.z;
		origin.x = (d2 * n1.y - d1 * n2.y) * invDenominator;
		origin.y = (d1 * n2.x - d2 * n1.x) * invDenominator;
		origin.z = 0;
	}
}


/*
 * returns whether the plane and the box overlap.
 */
bool IntersectPlaneBox(const Normal3f &n, const Point3f& p,
		const Bounds3f& box) {
	Point3f vmin, vmax;
	for (uint32_t q = 0; q < 3; ++q) {
		const Float v = p[q];

		if (n[q] > 0.0) {
			vmin[q] = box.pMin[q] - v;
			vmax[q] = box.pMax[q] - v;
		} else {
			vmin[q] = box.pMax[q] - v;
			vmax[q] = box.pMin[q] - v;
		}
	}

	if (pbrt::Dot(n, vmin) > 0.0)
		return false;
	if (pbrt::Dot(n, vmax) < 0.0)
		return false;
	return true;
}

bool IntersectTriangleBox(const Bounds3f& bounds,
		const Triangle3f& triangle) {
	// quick accept
	if (bounds.isInside(triangle.p0))
		return true;
	if (bounds.isInside(triangle.p1))
		return true;
	if (bounds.isInside(triangle.p2))
		return true;

	// false when the plane of the triangle does not intersect the box
	if (!IntersectPlaneBox(triangle.n, triangle.p0, bounds))
		return false;

	// false when the bounding boxes do not overlap
	if (!Overlaps(bounds, triangle.bounds))
		return false;

	// center the box at the origin
	const Point3f& boxCenter = 0.5 * bounds.pMin + 0.5 * bounds.pMax;
	const Vector3f& boxHalfWidth = bounds.pMax - boxCenter;

	// translate the triangle center
	const Vector3f& p0 = triangle.p0 - boxCenter;
	const Vector3f& p1 = triangle.p1 - boxCenter;
	const Vector3f& p2 = triangle.p2 - boxCenter;
	const Vector3f& f0 = p1 - p0;
	const Vector3f& f1 = p2 - p1;
	const Vector3f& f2 = p0 - p2;

	const Float f0y_abs = std::abs(f0.y);
	const Float f0z_abs = std::abs(f0.z);

	/**************************************************************************
	 * a00
	 **************************************************************************/
	{
		const Float ap0 = f0.y * p0.z - f0.z * p0.y;
		const Float ap2 = f0.y * p2.z - f0.z * p2.y;
		const Float min = std::min(ap0, ap2);
		const Float max = std::max(ap0, ap2);
		const Float r = boxHalfWidth.y * f0z_abs + boxHalfWidth.z * f0y_abs;

		if (min > r || max < -r)
			return false;
	}

	const Float f1y_abs = std::abs(f1.y);
	const Float f1z_abs = std::abs(f1.z);

	/**************************************************************************
	 * a01
	 **************************************************************************/

	{
		const Float ap0 = f1.y * p0.z - f1.z * p0.y;
		const Float ap2 = f1.y * p2.z - f1.z * p2.y;

		const Float min = std::min(ap0, ap2);
		const Float max = std::max(ap0, ap2);
		const Float r = boxHalfWidth.y * f1z_abs + boxHalfWidth.z * f1y_abs;

		if (min > r || max < -r)
			return false;
	}

	const Float f2y_abs = std::abs(f2.y);
	const Float f2z_abs = std::abs(f2.z);

	/**************************************************************************
	 * a02
	 **************************************************************************/

	{
		const Float ap0 = f2.y * p0.z - f2.z * p0.y;
		const Float ap1 = f2.y * p1.z - f2.z * p1.y;

		const Float min = std::min(ap0, ap1);
		const Float max = std::max(ap0, ap1);
		const Float r = boxHalfWidth.y * f2z_abs + boxHalfWidth.z * f2y_abs;

		if (min > r || max < -r)
			return false;
	}

	const Float f0x_abs = std::abs(f0.x);

	/**************************************************************************
	 * a10
	 **************************************************************************/
	{
		const Float ap0 = f0.z * p0.x - f0.x * p0.z;
		const Float ap2 = f0.z * p2.x - f0.x * p2.z;

		const Float min = std::min(ap0, ap2);
		const Float max = std::max(ap0, ap2);
		const Float r = boxHalfWidth.x * f0z_abs + boxHalfWidth.z * f0x_abs;

		if (min > r || max < -r)
			return false;
	}

	const Float f1x_abs = std::abs(f1.x);

	/**************************************************************************
	 * a11
	 **************************************************************************/
	{
		const Float ap0 = f1.z * p0.x - f1.x * p0.z;
		const Float ap2 = f1.z * p2.x - f1.x * p2.z;

		const Float min = std::min(ap0, ap2);
		const Float max = std::max(ap0, ap2);
		const Float r = boxHalfWidth.x * f1z_abs + boxHalfWidth.z * f1x_abs;

		if (min > r || max < -r)
			return false;
	}

	const Float f2x_abs = std::abs(f2.x);

	/**************************************************************************
	 * a12
	 **************************************************************************/
	{
		const Float ap0 = f2.z * p0.x - f2.x * p0.z;
		const Float ap1 = f2.z * p1.x - f2.x * p1.z;

		const Float min = std::min(ap0, ap1);
		const Float max = std::max(ap0, ap1);
		const Float r = boxHalfWidth.x * f2z_abs + boxHalfWidth.z * f2x_abs;

		if (min > r || max < -r)
			return false;
	}

	/**************************************************************************
	 * a20
	 **************************************************************************/
	{
		const Float ap0 = f0.x * p0.y - f0.y * p0.x;
		const Float ap2 = f0.x * p2.y - f0.y * p2.x;

		const Float min = std::min(ap0, ap2);
		const Float max = std::max(ap0, ap2);
		const Float r = boxHalfWidth.x * f0y_abs + boxHalfWidth.y * f0x_abs;

		if (min > r || max < -r)
			return false;
	}

	/**************************************************************************
	 * a21
	 **************************************************************************/
	{
		const Float ap0 = f1.x * p0.y - f1.y * p0.x;
		const Float ap2 = f1.x * p2.y - f1.y * p2.x;

		const Float min = std::min(ap0, ap2);
		const Float max = std::max(ap0, ap2);
		const Float r = boxHalfWidth.x * f1y_abs + boxHalfWidth.y * f1x_abs;

		if (min > r || max < -r)
			return false;
	}

	/**************************************************************************
	 * a22
	 **************************************************************************/
	{
		const Float ap0 = f2.x * p0.y - f2.y * p0.x;
		const Float ap1 = f2.x * p1.y - f2.y * p1.x;

		const Float min = std::min(ap0, ap1);
		const Float max = std::max(ap0, ap1);
		const Float r = boxHalfWidth.x * f2y_abs + boxHalfWidth.y * f2x_abs;

		if (min > r || max < -r)
			return false;
	}

	return true;
}

}

//	Float min, max;
//	/**************************************************************************
//	 * a00
//	 **************************************************************************/
//	{
//		min = f0.y * p0.z - f0.z * p0.y;
//		max = f0.y * p2.z - f0.z * p2.y;
//		if (min > max)
//			std::swap(min, max);
//		const Float r = boxHalfWidth.y * f0z_abs + boxHalfWidth.z * f0y_abs;
//
//		if (min > r || max < -r)
//			return false;
//	}
//
//	const Float f1y_abs = std::abs(f1.y);
//	const Float f1z_abs = std::abs(f1.z);
//
//	/**************************************************************************
//	 * a01
//	 **************************************************************************/
//
//	{
//		min = f1.y * p0.z - f1.z * p0.y;
//		max = f1.y * p2.z - f1.z * p2.y;
//		if (min > max)
//			std::swap(min, max);
//		const Float r = boxHalfWidth.y * f1z_abs + boxHalfWidth.z * f1y_abs;
//
//		if (min > r || max < -r)
//			return false;
//	}
//
//	const Float f2y_abs = std::abs(f2.y);
//	const Float f2z_abs = std::abs(f2.z);
//
//	/**************************************************************************
//	 * a02
//	 **************************************************************************/
//
//	{
//		min = f2.y * p0.z - f2.z * p0.y;
//		max = f2.y * p1.z - f2.z * p1.y;
//		if (min > max)
//			std::swap(min, max);
//		const Float r = boxHalfWidth.y * f2z_abs + boxHalfWidth.z * f2y_abs;
//
//		if (min > r || max < -r)
//			return false;
//	}
//
//	const Float f0x_abs = std::abs(f0.x);
//
//	/**************************************************************************
//	 * a10
//	 **************************************************************************/
//	{
//		min = f0.z * p0.x - f0.x * p0.z;
//		max = f0.z * p2.x - f0.x * p2.z;
//		if (min > max)
//			std::swap(min, max);
//		const Float r = boxHalfWidth.x * f0z_abs + boxHalfWidth.z * f0x_abs;
//
//		if (min > r || max < -r)
//			return false;
//	}
//
//	const Float f1x_abs = std::abs(f1.x);
//
//	/**************************************************************************
//	 * a11
//	 **************************************************************************/
//	{
//		min = f1.z * p0.x - f1.x * p0.z;
//		max = f1.z * p2.x - f1.x * p2.z;
//		if (min > max)
//			std::swap(min, max);
//		const Float r = boxHalfWidth.x * f1z_abs + boxHalfWidth.z * f1x_abs;
//
//		if (min > r || max < -r)
//			return false;
//	}
//
//	const Float f2x_abs = std::abs(f2.x);
//
//	/**************************************************************************
//	 * a12
//	 **************************************************************************/
//	{
//		min = f2.z * p0.x - f2.x * p0.z;
//		max = f2.z * p1.x - f2.x * p1.z;
//		if (min > max)
//			std::swap(min, max);
//		const Float r = boxHalfWidth.x * f2z_abs + boxHalfWidth.z * f2x_abs;
//
//		if (min > r || max < -r)
//			return false;
//	}
//
//	/**************************************************************************
//	 * a20
//	 **************************************************************************/
//	{
//		min = f0.x * p0.y - f0.y * p0.x;
//		max = f0.x * p2.y - f0.y * p2.x;
//		if (min > max)
//			std::swap(min, max);
//		const Float r = boxHalfWidth.x * f0y_abs + boxHalfWidth.y * f0x_abs;
//
//		if (min > r || max < -r)
//			return false;
//	}
//
//	/**************************************************************************
//	 * a21
//	 **************************************************************************/
//	{
//		min = f1.x * p0.y - f1.y * p0.x;
//		max = f1.x * p2.y - f1.y * p2.x;
//		if (min > max)
//			std::swap(min, max);
//		const Float r = boxHalfWidth.x * f1y_abs + boxHalfWidth.y * f1x_abs;
//
//		if (min > r || max < -r)
//			return false;
//	}
//
//	/**************************************************************************
//	 * a22
//	 **************************************************************************/
//	{
//		min = f2.x * p0.y - f2.y * p0.x;
//		max = f2.x * p1.y - f2.y * p1.x;
//		if (min > max)
//			std::swap(min, max);
//		const Float r = boxHalfWidth.x * f2y_abs + boxHalfWidth.y * f2x_abs;
//
//		if (min > r || max < -r)
//			return false;
//	}

