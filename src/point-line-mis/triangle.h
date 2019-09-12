/*
 * triangle.h
 *
 *  Created on: Jul 17, 2017
 *      Author: niels
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PHD_STOCHASTIC_TRIANGLE_H
#define PHD_STOCHASTIC_TRIANGLE_H

#include "pbrt.h"
#include "geometry.h"
#include "interaction.h"

using pbrt::Float;
using pbrt::Point3f;
using pbrt::Vector3f;
using pbrt::Normal3f;
using pbrt::Normalize;
using pbrt::Cross;
using pbrt::AbsDot;
using pbrt::Dot;
using pbrt::Inv2Pi;
using pbrt::SurfaceInteraction;
using pbrt::IntPow;
using pbrt::Pi;
using pbrt::Clamp;

namespace niels {

class Triangle {
public:
	Triangle(const Point3f& p0, const Point3f& p1, const Point3f& p2);

	void Clip(const Point3f& p, const Normal3f& n,
			std::vector<Triangle> &clipped) const;

	void Orient(const Triangle& triangle);
	Normal3f Normal() const;

	bool HasVertex(const Point3f& p, Float epsilon = 0) const;
	Float Area() const;
	bool Above(const Point3f& p, const Normal3f& n) const;

	Float Diffuse(const Point3f& p, const Normal3f& nn) const;
	Float SolidAngle(const Point3f& p) const;
	Float SignedSolidAngle(const Point3f& p) const;
	Float Phong(const SurfaceInteraction& isect, int n) const;

	void Print(FILE* file) const;
	const Point3f& operator[](int index) const;
private:
	Point3f vertices[3];
};

}

#endif /* PHD_STOCHASTIC_TRIANGLE_H_ */
