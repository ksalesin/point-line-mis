/*
 * intersections.h
 *
 *  Created on: Jan 25, 2016
 *      Author: niels
 */

#ifndef PBRT_LINESAMPLING_INTERSECTIONS_H_
#define PBRT_LINESAMPLING_INTERSECTIONS_H_

#include "geometry.h"
#include "pbrt.h"
#include "point-line-mis/line3.h"
#include "point-line-mis/triangle3.h"
#include "point-line-mis/interval.h"
#include "point-line-mis/lineinteraction.h"

using pbrt::Float;
using pbrt::Point3f;
using pbrt::Vector3f;
using pbrt::Normal3f;
using pbrt::Bounds3f;
using pbrt::Dot;
using pbrt::Cross;

namespace linesampling {

Interval FindShadowedInterval(const Point3f& origin, const Point3f& lineOrigin,
		const Point3f& lineDestination, const Point3f& p1, const Point3f&p2);

Interval FindShadowedInterval(const LineVisibilityQuery3f& query,
		const Line3f &isect);

Float ShortestLineParameter(const Point3f&, const Point3f&, const Point3f&,
		const Point3f&);

Float ShortestLineParameter(const Point3f&, const Vector3f&, const Point3f&,
		const Point3f&);

void IntersectPlanePlane(const Normal3f &n1, Float d1, const Normal3f &n2,
		Float d2, Point3f &origin, Vector3f &direction);

bool IntersectTriangleTriangle(const Triangle3f &u, const Triangle3f &v,
		Line3f &isect);

bool IntersectPlaneBox(const Normal3f &n, const Point3f& p,
		const Bounds3f& box);

bool IntersectTriangleBox(const Bounds3f& bounds, const Triangle3f& triangle);

}
#endif /* LINESAMPLING_INTERSECTIONS_H_ */
