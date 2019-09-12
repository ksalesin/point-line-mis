/*
 * line.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: niels
 */

#include "point-line-mis/line.h"

namespace niels {

Line::Line(const Point3f& p0, const Point3f &p1) :
		p0(p0), p1(p1) {
}

Point3f Line::operator()(Float t) const {
	return Lerp(t, p0, p1);
}

bool Line::IntersectPlane(const Point3f& p, const Normal3f& n,
		Float *tHit) const {
	const Vector3f& d = p0 - p1;

	Float numerator = Dot(p - p0, n);
	Float denominator = Dot(n, d);

	if (denominator == 0) {
		// plane and line are parallel
		if (numerator != 0)
			return false;
		if (tHit)
			*tHit = 0;
		return true;
	} else {
		const Float t = numerator / denominator;
		if (t < 0 || t > 1)
			return false;

		if (tHit)
			*tHit = t;
		return true;
	}
}

bool Line::AbovePlane(const Point3f& p, const Normal3f& n) const {
	return Dot(p0 - p, n) >= 0 && Dot(p1 - p, n) >= 0;
}

bool IntersectLinePlane(const Point3f& p0, const Point3f& p1, const Point3f &p,
		const Normal3f& n, Point3f *isect) {
	const Vector3f& d = p1 - p0;

	const Float numerator = Dot(p - p0, n);
	const Float denominator = Dot(n, d);

	if (denominator == 0) {
		// plane and line are parallel
		if (numerator != 0)
			return false;
		if (isect)
			*isect = p0;
		return true;
	} else {
		const Float t = numerator / denominator;
		if (isect)
			*isect = Lerp(t, p0, p1);
		return true;
	}
}
}

