/*
 * triangle.h
 *
 *  Created on: Jan 13, 2016
 *      Author: niels
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_LINESAMPLING_TRIANGLE_H
#define PBRT_LINESAMPLING_TRIANGLE_H

// linesampling/triangle.h*
#include "pbrt.h"
#include "geometry.h"

using pbrt::Point3;
using pbrt::Vector3;
using pbrt::Normal3;
using pbrt::Bounds3;
using pbrt::Cross;
using pbrt::Float;

namespace linesampling {

// Triangle3 Declaration
template<typename T>
class Triangle3 {
public:
	// Triangle3 Public Methods
	Triangle3(const Point3<T>& p0, const Point3<T>& p1, const Point3<T>& p2) :
			p0(p0), p1(p1), p2(p2) {
		const Vector3<T> &e1 = p1 - p0;
		const Vector3<T> &e2 = p2 - p0;
		const Vector3<T> &cross = Cross(e1, e2);

		T length = cross.Length();

		bounds = Union(Bounds3<T>(p0, p1), p2);
		area = 0.5 * length;
		if (length == 0) {
			n = Normal3<T>(cross);
		} else {
			n = Normal3<T>(cross / length);
		}
		
		Float ax = std::abs(n.x);
		Float ay = std::abs(n.y);
		Float az = std::abs(n.z);
		if (ax < ay && ax < az) {
			nMinIndex = 0;
			nMidIndex = ay < az ? 1 : 2;
			nMaxIndex = ay < az ? 2 : 1;
		} else if (n.y < n.z) {
			nMinIndex = 1;
			nMidIndex = ax < az ? 0 : 2;
			nMaxIndex = ax < az ? 2 : 0;
		} else {
			nMinIndex = 2;
			nMidIndex = ax < ay ? 0 : 1;
			nMaxIndex = ax < ay ? 1 : 0;
		}
	}

	T Area() const {
		return 0.5 * Cross(p1 - p0, p2 - p0).Length();
	}

	Point3<T> p0;
	Point3<T> p1;
	Point3<T> p2;
	Normal3<T> n;
	Bounds3<T> bounds;
	T area;
	uint32_t nMinIndex;
	uint32_t nMidIndex;
	uint32_t nMaxIndex;
};

typedef Triangle3<Float> Triangle3f;

}

#endif /* LINESAMPLING_TRIANGLE_H */
