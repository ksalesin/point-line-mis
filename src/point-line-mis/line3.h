/*
 * line3.h
 *
 *  Created on: Jan 4, 2016
 *      Author: niels
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_LINESAMPLING_LINE_H
#define PBRT_LINESAMPLING_LINE_H

#include "geometry.h"

using pbrt::Float;
using pbrt::Point3;
using pbrt::Vector3;

// linesampling namespace
namespace linesampling {

// Line3 Public Declaration
template<typename T>
class Line3 {
public:
	// Line3 Public Methods
	Line3() {
	}

	Line3(const Point3<T>& from, const Point3<T>& to) :
			p1(from), p2(to) {
	}

	Line3(const Point3<T>& origin, const Vector3<T> &direction) :
			p1(origin), p2(origin + direction) {
	}

	T length() const {
		return (p2 - p1).Length();
	}

	Point3<T> p1;
	Point3<T> p2;
};

typedef Line3<Float> Line3f;

}
#endif /* LINESAMPLING_LINE_H_ */
