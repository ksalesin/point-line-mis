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

#ifndef PBRT_LIGHTS_DIFFUSE_H
#define PBRT_LIGHTS_DIFFUSE_H

// lights/diffuse.h*
#include "pbrt.h"
#include "light.h"
#include "primitive.h"
#include "texture.h"
#include "mipmap.h"

namespace pbrt {

// DiffuseAreaLight Declarations
class DiffuseAreaLight: public AreaLight {
public:
	// DiffuseAreaLight Public Methods
	DiffuseAreaLight(const Transform &LightToWorld,
			const MediumInterface &mediumInterface, const Spectrum &Le,
			int xSamples, int ySamples, const std::shared_ptr<Shape> &shape, 
			const std::string &texmap, bool twoSided = false);
		
	Spectrum L(const Interaction &intr, const Vector3f &w) const override;
	Spectrum L(const Interaction &intr, const Vector3f &w, const Point2f& uv) const;
	Spectrum LScaled(const Interaction &intr, const Vector3f &w) const;

	Spectrum Power() const;
	Spectrum Sample_Li(const Interaction &ref, const Point2f &u, Vector3f *wo,
			Float *pdf, VisibilityTester *vis) const override;
	Spectrum Sample_Texture_Li(const Interaction &ref, const Point2f &u,
		Vector3f *wi, Float *pdf, VisibilityTester *vis) const;
	Float Pdf_Li(const Interaction &, const Vector3f &) const override;
	Float Pdf_Li_Texture(const Interaction &, const Vector3f &) const;

	Spectrum Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
			Ray *ray, Normal3f *nLight, Float *pdfPos, Float *pdfDir) const override;
	void Pdf_Le(const Ray &, const Normal3f &, Float *pdfPos,
			Float *pdfDir) const override;

	virtual bool CanEvaluateAnalytically(const SurfaceInteraction& isect,
			Spectrum * spectrum = nullptr) const override;
	virtual bool CanIlluminate(const Interaction& isect) const override;

	// Sample line on light by surface area method
	Spectrum Sample_Line_Li(const Interaction &ref, const Point2f &u,
			linesampling::LineInteraction &interaction, bool importanceSample =
					true) const override;
	virtual Float PdfLine(const Interaction &ref, const Vector3f &wi, Point3f &pLight,
			LineInteraction &lineInteraction, const Float& lineAngle, bool importanceSample) const override;
	
	// Sample point uniformly on light by solid angle method
	virtual Spectrum Sample_Li_SA(const Interaction &ref, const Point2f &u,
			Vector3f *wi, Float *pdf, VisibilityTester *vis) const override;
	virtual Float Pdf_SA(const Interaction &ref, const Vector3f &wi) const override;

	// Sample line on light by solid angle method
	virtual Spectrum Sample_SphQuadLine_Li(const Point3f& p, const Float u, int strategy,
			LineInteraction &interaction) const override;
	Float PdfSQLine(const Interaction &ref, const Vector3f& wi, Point3f& pLight, LineInteraction &interaction, int strategy) const override {
		return shape->PdfSQLine(ref, wi, pLight, interaction, strategy);
	}

	// Sample a point within given line by solid angle method
	virtual Point3f Sample_SphQuadPointInLine(const IntervalCollection collection, const LineInteraction &interaction, 
			const Point3f& p, const Float v, Float *pointPdf, int strategy) const override {
		return shape->SampleSQPoint(collection, interaction, p, v, pointPdf, strategy);
	}
	Float PdfSQPoint(const IntervalCollection intervals, const LineInteraction &interaction, 
			const Point3f& p, int strategy) const override {
		return shape->PdfSQPoint(intervals, interaction, p, strategy);
	}

	// Visualizes lines on light in given direction (same as how line lighting would be evaluated)
	Spectrum LLines(const Interaction &intr, const Vector3f &w, float direction) const override {
		if (Lemit.IsBlack())
			fprintf(stderr, "Lemit is black!");

		if (!twoSided && Dot(intr.n, w) < 0)
			return Spectrum(0.f);

		return (shape->OnColoredLine(intr, direction)) ? Spectrum::FromRGB(0.f, 0.f, 1.f) : Lemit;
	}

	virtual Float smoothstep(const Point3f& pLight, const Float smoothstepWidth) const override {
		return shape->smoothstep(pLight, smoothstepWidth);
	}

protected:
	// DiffuseAreaLight Protected Data
	const Spectrum Lemit;
	std::shared_ptr<Shape> shape;
	// Added after book publication: by default, DiffuseAreaLights still
	// only emit in the hemimsphere around the surface normal.  However,
	// this behavior can now be overridden to give emission on both sides.
	const bool twoSided;
	const Float area;

	// Texture
	bool hasTexture;
	std::unique_ptr<MIPMap<RGBSpectrum>> Lmap;
	std::unique_ptr<Distribution2D> distribution;
};

std::shared_ptr<AreaLight> CreateDiffuseAreaLight(const Transform &light2world,
		const Medium *medium, const ParamSet &paramSet,
		const std::shared_ptr<Shape> &shape);

}  // namespace pbrt

#endif  // PBRT_LIGHTS_DIFFUSE_H
