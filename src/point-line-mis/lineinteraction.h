/*
 * lineinteraction.h
 *
 *  Created on: Jan 4, 2016
 *      Author: niels
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_LINESAMPLING_LINEINTERACTION_H_
#define PBRT_LINESAMPLING_LINEINTERACTION_H_

// linesampling/lineinteraction.h*
#include "pbrt.h"
#include "geometry.h"
#include "medium.h"
#include "point-line-mis/line3.h"
#include "point-line-mis/triangle3.h"

using pbrt::Float;
using pbrt::Point3f;
using pbrt::Vector3f;
using pbrt::Normal3f;
using pbrt::MediumInterface;

// linesampling namespace
namespace linesampling {

struct LineEvaluation {
	LineEvaluation(const Normal3f& n, const Normal3f& sn, const Normal3f& ln,
			const Vector3f& wi, const Vector3f& p1, const Vector3f& p2) :
			n(n), sn(sn), ln(ln), wi(wi), p1(p1), p2(p2) {
	}

	const Normal3f n;
	const Normal3f sn;
	const Normal3f ln;
	const Vector3f wi;
	const Vector3f p1;
	const Vector3f p2;
};

// LineInteraction Public Declaration
class LineInteraction {
public:
//	// LineInteraction Public Methods
//	LineInteraction() :
//			time(0), pdf(1.0), length(0) {
//	}
//
//	LineInteraction(const Line3f& line, const Normal3f& n, Float time,
//			Float pdf, Float length, const MediumInterface& interface) :
//			line(line), n(n), time(time), pdf(pdf), length(length), mediumInterface(
//					interface) {
//	}
//
//	std::string to_string() const {
//		std::ostringstream stream;
//		stream.precision(16);
//		stream << line.p1.() << " -> " + line.p2.to_string()
//				<< "; pdf = " << pdf << "; length = " << length;
//		return stream.str();
//	}

	Line3f line; 	// the line which is sampled
	Normal3f n;		// the normal at the light source
	Float time;		// time when the sample was generated
	Float pdf;		// probability that the line sample was generated
	Float length;	// the length of the line segment
	Float xu;       // for spherical quad sampling, the offset xu that was sampled
	MediumInterface mediumInterface;
};

template<typename T>
class LineVisibilityQuery3: public Triangle3<T> {
public:
	// LineVisibilityQuery Public Methods
	LineVisibilityQuery3(const Point3<T> &queryPoint, const Line3<T> &line) :
			Triangle3<T>(queryPoint, line.p1, line.p2) {
		u = line.p2 - line.p1;
		w = line.p1 - line.p0;
	}

	LineVisibilityQuery3(const Point3<T> &queryPoint, const Point3<T> &p1,
			const Point3<T> &p2) :
			Triangle3<T>(queryPoint, p1, p2) {
		u = p2 - p1;
		w = p1 - queryPoint;
	}

	Vector3<T> u;
	Vector3<T> w;
};

typedef LineVisibilityQuery3<Float> LineVisibilityQuery3f;

}
#endif /* LINESAMPLING_LINEINTERACTION_H_ */
