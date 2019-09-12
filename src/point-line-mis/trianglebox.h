/*
 * intersections.h
 *
 *  Created on: Jan 25, 2016
 *      Author: niels
 */

#ifndef PBRT_LINESAMPLING_TRIANGLEBOX_H_
#define PBRT_LINESAMPLING_TRIANGLEBOX_H_

#include "geometry.h"
#include "pbrt.h"
#include "point-line-mis/line3.h"
#include "point-line-mis/triangle3.h"
#include "point-line-mis/interval.h"
#include "point-line-mis/intersections.h"
#include "point-line-mis/lineinteraction.h"

namespace linesampling {

//bool IntersectTrianglePlaneX(Float position, const Bounds3f &box,
//		const LineVisibilityQuery3f &query, Interval &interval);
//
//bool IntersectTrianglePlaneY(Float position, const Bounds3f &box,
//		const LineVisibilityQuery3f &query, Interval &interval);
//
//bool IntersectTrianglePlaneZ(Float position, const Bounds3f &box,
//		const LineVisibilityQuery3f &query, Interval &interval);

bool IntersectTriangleBoxBounds(const Bounds3f& bounds,
		const LineVisibilityQuery3f& triangle, Interval& interval);

bool IntersectTriangleBoxBoundsOptimized(const Bounds3f& bounds,
		const LineVisibilityQuery3f& triangle, Interval& interval);

}
#endif /* PBRT_LINESAMPLING_TRIANGLEBOX_H_ */
