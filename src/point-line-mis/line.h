/*
 * line.h
 *
 *  Created on: Jul 17, 2017
 *      Author: niels
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PHD_STOCHASTIC_LINE_H
#define PHD_STOCHASTIC_LINE_H

#include "pbrt.h"
#include "geometry.h"

using pbrt::Float;
using pbrt::Point3f;
using pbrt::Normal3f;
using pbrt::Vector3f;

namespace niels {

class Line {
public:
	Line(const Point3f& p0, const Point3f& p1);

	Point3f operator()(Float t) const;
	bool IntersectPlane(const Point3f& p, const Normal3f& n, Float *tHit) const;

	bool AbovePlane(const Point3f& p, const Normal3f& n) const;
private:
	const Point3f& p0;
	const Point3f& p1;
};

bool IntersectLinePlane(const Point3f& p0, const Point3f& p1, const Point3f &p,
		const Normal3f& n, Point3f *isect);

}

#endif /* PHD_STOCHASTIC_LINE_H_ */
