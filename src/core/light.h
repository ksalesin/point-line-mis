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

/* 
   LLines (to visualize line direction on light source) added by Kate (12/2018)
*/

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_LIGHT_H
#define PBRT_CORE_LIGHT_H

// core/light.h*
#include "pbrt.h"
#include "memory.h"
#include "interaction.h"
#include "point-line-mis/lineinteraction.h"
#include "point-line-mis/interval.h"

using linesampling::LineInteraction;
using linesampling::IntervalCollection;

namespace pbrt {

// LightFlags Declarations
enum class LightFlags
	: int {
		DeltaPosition = 1, DeltaDirection = 2, Area = 4, Infinite = 8
};

inline bool IsDeltaLight(int flags) {
	return flags & (int) LightFlags::DeltaPosition
			|| flags & (int) LightFlags::DeltaDirection;
}

// Light Declarations
class Light {
public:
	// Light Interface
	virtual ~Light();
	Light(int flags, const Transform &LightToWorld,
			const MediumInterface &mediumInterface, int xSamples = 1,
			int ySamples = 1);
	virtual Spectrum Sample_Li(const Interaction &ref, const Point2f &u,
			Vector3f *wi, Float *pdf, VisibilityTester *vis) const = 0;
	virtual Spectrum Power() const = 0;
	virtual void Preprocess(const Scene &scene) {
	}
	virtual Spectrum Le(const RayDifferential &r) const;
	virtual Float Pdf_Li(const Interaction &ref, const Vector3f &wi) const = 0;
	virtual Spectrum Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
			Ray *ray, Normal3f *nLight, Float *pdfPos, Float *pdfDir) const = 0;
	virtual void Pdf_Le(const Ray &ray, const Normal3f &nLight, Float *pdfPos,
			Float *pdfDir) const = 0;

	virtual Spectrum L(const Interaction &intr, const Vector3f &w) const {
		Error("Unimplemented Light::L called!");
		return 0.f;
	}

	virtual bool CanEvaluateAnalytically(const SurfaceInteraction& isect,
			Spectrum * spectrum = nullptr) const {
		return false;
	}

	virtual bool CanIlluminate(const Interaction& isect) const {
		return true;
	}
	
	// Sample line on light by surface area method
	virtual Spectrum Sample_Line_Li(const Interaction &ref, const Point2f &u,
			LineInteraction &lineInteraction, bool importanceSample = true) const {
		Error("Unimplemented Light::Sample_Line_Li called!");
		return Spectrum(0.f);
	}
	virtual Float PdfLine(const Interaction &ref, const Vector3f &wi, Point3f &pLight,
			LineInteraction &lineInteraction, const Float& lineAngle, bool importanceSample) const {
		Error("Unimplemented Light::PdfLine called!");
		return 0.f;
	}

	// Sample point uniformly on light by solid angle method
	virtual Spectrum Sample_Li_SA(const Interaction &ref, const Point2f &u,
			Vector3f *wi, Float *pdf, VisibilityTester *vis) const {
		Error("Unimplemented Light::Sample_Li_SA called!");
		return 0.f;
	}
	virtual Float Pdf_SA(const Interaction &ref, const Vector3f &wi) const {
		Error("Unimplemented Light::PdfSA called!");
		return 0.f;
	}

	// Sample line on light by solid angle method
	virtual Spectrum Sample_SphQuadLine_Li(const Point3f& p, const Float u, int strategy,
			LineInteraction &interaction) const {
		Error("Unimplemented Light::Sample_SphQuadLine_Li called!");
		return 0.f;
	}
	virtual Float PdfSQLine(const Interaction &ref, const Vector3f& wi, Point3f& pLight, LineInteraction &interaction, int strategy) const {
		Error("Unimplemented Light::PdfSQLine method called!");
		return 0.;
	}

	// Sample a point within given line by solid angle method
	virtual Point3f Sample_SphQuadPointInLine(const IntervalCollection collection, const LineInteraction &interaction, 
			const Point3f& p, const Float v, Float *pointPdf, int strategy) const {
		Error("Unimplemented Light::Sample_SphQuadPointInLine called!");
		return Point3f();
	}
	virtual Float PdfSQPoint(const IntervalCollection intervals, const LineInteraction &interaction, 
			const Point3f& p, int strategy) const {
		Error("Unimplemented Light::PdfSQPoint method called!");
		return 0.;
	}

	virtual Float smoothstep(const Point3f& pLight, const Float smoothstepWidth) const {
		Error("Unimplemented Light::smoothstep method called!");
		return 0.;
	}

	// Light Public Data
	const int flags;
	const int xSamples;
	const int ySamples;
	const MediumInterface mediumInterface;

protected:
	// Light Protected Data
	const Transform LightToWorld, WorldToLight;
};

class VisibilityTester {
public:
	VisibilityTester() {
	}
	// VisibilityTester Public Methods
	VisibilityTester(const Interaction &p0, const Interaction &p1) :
			p0(p0), p1(p1) {
	}
	const Interaction &P0() const {
		return p0;
	}
	const Interaction &P1() const {
		return p1;
	}
	bool Unoccluded(const Scene &scene) const;
	Spectrum Tr(const Scene &scene, Sampler &sampler) const;

private:
	Interaction p0, p1;
};

class AreaLight: public Light {
public:
	// AreaLight Interface
	AreaLight(const Transform &LightToWorld, const MediumInterface &medium,
			int xSamples, int ySamples);
	// virtual Spectrum L(const Interaction &intr, const Vector3f &w);
	virtual Spectrum LScaled(const Interaction &intr, const Vector3f &w) const = 0;
	virtual Spectrum LLines(const Interaction &intr, const Vector3f &w, float direction) const {
		Error("Unimplemented Light::LLines method called!");
		return Spectrum(0.f);
	}
};

}  // namespace pbrt

#endif  // PBRT_CORE_LIGHT_H
