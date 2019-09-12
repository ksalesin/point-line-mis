/*
 pbrt source code is Copyright(c) 1998-2016
 Matt Pharr, Greg Humphreys, and Wenzel Jakob.

 This file is part of pbrt.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 - Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 - Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_SHAPE_H
#define PBRT_CORE_SHAPE_H

// core/shape.h*
#include "pbrt.h"
#include "geometry.h"
#include "interaction.h"
#include "memory.h"
#include "transform.h"
#include "point-line-mis/interval.h"
#include "point-line-mis/lineinteraction.h"

using linesampling::LineInteraction;
using linesampling::LineVisibilityQuery3f;
using linesampling::IntervalCollection;

namespace pbrt {

// Shape Declarations
class Shape {
public:
	// Shape Interface
	Shape(const Transform *ObjectToWorld, const Transform *WorldToObject,
			bool reverseOrientation);
	virtual ~Shape();
	virtual Bounds3f ObjectBound() const = 0;
	virtual Bounds3f WorldBound() const;
	virtual bool Intersect(const Ray &ray, Float *tHit,
			SurfaceInteraction *isect, bool testAlphaTexture = true) const = 0;
	virtual bool IntersectP(const Ray &ray,
			bool testAlphaTexture = true) const {
		return Intersect(ray, nullptr, nullptr, testAlphaTexture);
	}
	virtual bool IntersectL(const LineVisibilityQuery3f &query,
			const Normal3f &normal, Float time,
			IntervalCollection &intersection) const {
		// Error("Unimplemented Shape::IntersectL method called!");
		return false;
	}

	virtual Float Area() const = 0;
	// Sample a point on the surface of the shape and return the PDF with
	// respect to area on the surface.
	virtual Interaction Sample(const Point2f &u, Float *pdf) const = 0;
	virtual Float Pdf(const Interaction &) const {
		return 1 / Area();
	}

	// Sample a point on the shape given a reference point |ref| and
	// return the PDF with respect to solid angle from |ref|.
	virtual Interaction Sample(const Interaction &ref, const Point2f &u,
			Float *pdf) const;
	virtual Float Pdf(const Interaction &ref, const Vector3f &wi) const;
	virtual Float Pdf(SurfaceInteraction& shapeRef, const Interaction &ref, const Vector3f &wi) const;

	// Returns the solid angle subtended by the shape w.r.t. the reference
	// point p, given in world space. Some shapes compute this value in
	// closed-form, while the default implementation uses Monte Carlo
	// integration; the nSamples parameter determines how many samples are
	// used in this case.
	virtual Float SolidAngle(const Point3f &p, int nSamples = 512) const;

	virtual bool CanEvaluateAnalytically(const SurfaceInteraction& isect,
			Spectrum * spectrum) const {
		return false;
	}
	virtual bool CanIlluminate(const Interaction& isect) const {
		return true;
	}

	// Sample a line on shape by surface area method
	virtual LineInteraction LineSample(const Point2f &u, bool importanceSample =
			true) const {
		Error("Unimplemented Shape::LineSample method called!");
		return linesampling::LineInteraction();
	}
	virtual Float PdfLine(const Interaction &ref, const Vector3f &wi, Point3f &pLight,
			LineInteraction &lineInteraction, const Float& lineAngle, bool importanceSample) const {
		Error("Unimplemented Shape::PdfLine method called!");
		return 0.f;
	}

	// Sample a line on shape by solid angle method
	virtual LineInteraction SampleSQLine(const Point3f& p, const Float u, int strategy) const {
		Error("Unimplemented Shape::SampleSQLine method called!");
		return linesampling::LineInteraction();
	}
	virtual Float PdfSQLine(const Interaction &ref, const Vector3f& wi, Point3f& pLight, LineInteraction &interaction, int strategy) const {
		Error("Unimplemented Shape::PdfSQLine method called!");
		return 0.;
	}

	// Sample a point within given line by solid angle method
	virtual Point3f SampleSQPoint(const IntervalCollection intervals, const LineInteraction &interaction, 
			const Point3f& p, const Float v, Float *pointPdf, int strategy) const {
		Error("Unimplemented Shape::SampleSQPoint method called!");
		return Point3f();
	}
	virtual Float PdfSQPoint(const IntervalCollection intervals, const LineInteraction &interaction, 
			const Point3f& p, int strategy) const {
		Error("Unimplemented Shape::PdfSQPoint method called!");
		return 0.;
	}

	// Get UV coordinates of given point on quad
	virtual Point2f GetUV(const Point3f& p) const {
		Error("Unimplemented Shape::GetUV method called!");
		return Point2f();
	}

	// Sample point on quad given UV coordinates
	virtual Interaction SampleUV(const Point2f& uv) const {
		Error("Unimplemented Shape::SampleUV method called!");
		return Interaction();
	}

	virtual Float smoothstep(const Point3f& pLight, const Float smoothstepWidth) const {
		Error("Unimplemented Shape::smoothstep method called!");
		return 0.;
	}

	// Sample a point uniformly on shape by solid angle method
	virtual Interaction SampleSolidAngle(const Point3f& p, const Point2f &sample, Float* pdf) const {
		Error("Unimplemented Shape::SampleSolidAngle (in Shape) method called!");
		return Interaction();
	}
	virtual Float PdfSolidAngle(const Point3f& p) const {
		Error("Unimplemented Shape::PdfSolidAngle method called!");
		return 0.f;
	}

	// Visualizes lines on quad in given direction (same as how line lighting would be evaluated)
	virtual bool OnColoredLine(const Interaction &intr, float direction) const {
		Error("Unimplemented Shape::OnColoredLine method called!");
		return false;
	}

	// Sample a point on the shape given a reference point |ref| and
	// return the PDF with respect to solid angle from |ref|.
	virtual Interaction SampleSolidAngle(const Interaction &ref, const Point2f &u, Float *pdf) const;
	virtual Float PdfSolidAngle(const Interaction &ref, const Vector3f &wi) const;

	// Shape Public Data
	const Transform *ObjectToWorld, *WorldToObject;
	const bool reverseOrientation;
	const bool transformSwapsHandedness;
};

}  // namespace pbrt

#endif  // PBRT_CORE_SHAPE_H
