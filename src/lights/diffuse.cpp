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

// lights/diffuse.cpp*
#include "lights/diffuse.h"
#include "paramset.h"
#include "sampling.h"
#include "shapes/quad.h"
#include "shapes/triangle.h"
#include "imageio.h"
#include "stats.h"

namespace pbrt {

// DiffuseAreaLight Method Definitions
DiffuseAreaLight::DiffuseAreaLight(const Transform &LightToWorld,
		const MediumInterface &mediumInterface, const Spectrum &Lemit,
		int xSamples, int ySamples, const std::shared_ptr<Shape> &shape, 
		const std::string &texmap, bool twoSided) :
		AreaLight(LightToWorld, mediumInterface, xSamples, ySamples), Lemit(Lemit), 
		shape(shape), twoSided(twoSided), area(shape->Area()) {
	// Warn if light has transformation with non-uniform scale, though not
	// for Triangles, since this doesn't matter for them.
	if (WorldToLight.HasScale()
			&& !(dynamic_cast<const Triangle *>(shape.get()) != nullptr
					|| dynamic_cast<const Quad *>(shape.get()) != nullptr))
		Warning("Scaling detected in world to light transformation! "
				"The system has numerous assumptions, implicit and explicit, "
				"that this transform will have no scale factors in it. "
				"Proceed at your own risk; your image may have errors.");

	hasTexture = texmap != "";
	
	// Read texel data from _texmap_ and initialize _Lmap_
	Point2i resolution;
	std::unique_ptr<RGBSpectrum[]> texels(nullptr);
	if (hasTexture) {
		texels = ReadImage(texmap, &resolution);
		if (texels)
			for (int i = 0; i < resolution.x * resolution.y; ++i)
				texels[i] *= Lemit.ToRGBSpectrum();
	}
	if (!texels) {
		resolution.x = resolution.y = 1;
		texels = std::unique_ptr<RGBSpectrum[]>(new RGBSpectrum[1]);
		texels[0] = Lemit.ToRGBSpectrum();
	}
	Lmap.reset(new MIPMap<RGBSpectrum>(resolution, texels.get()));

	// Initialize sampling PDFs

	// Compute scalar-valued image _img_ from environment map
	int width = 2 * Lmap->Width(), height = 2 * Lmap->Height();
	std::unique_ptr<Float[]> img(new Float[width * height]);
	float fwidth = 0.5f / std::min(width, height);
	ParallelFor([&](int64_t v) {
		Float vp = (v + .5f) / (Float)height;
		// Float sinTheta = std::sin(Pi * (v + .5f) / height);
		for (int u = 0; u < width; ++u) {
			Float up = (u + .5f) / (Float)width;
			img[u + v * width] = Lmap->Lookup(Point2f(up, vp), fwidth).y(); // y = luminance
			// img[u + v * width] *= sinTheta;
		}
	}, height, 32);

	// Compute sampling distributions for rows and columns of image
	distribution.reset(new Distribution2D(img.get(), width, height));
}

Spectrum DiffuseAreaLight::Power() const {
	return (twoSided ? 2 : 1) * Lemit * area * Pi;
}

Spectrum DiffuseAreaLight::L(const Interaction &intr, const Vector3f &w) const {
	if (Lemit.IsBlack())
		fprintf(stderr, "Lemit is black!");
	
	if (!twoSided && Dot(intr.n, w) <= 0) 
		return Spectrum(0.f);

	if (hasTexture) {
		const Point2f uv = shape->GetUV(intr.p);
		return Spectrum(Lmap->Lookup(uv), SpectrumType::Illuminant);
	}

	return Lemit;
}

Spectrum DiffuseAreaLight::L(const Interaction &intr, const Vector3f &w, const Point2f& uv) const {
	if (Lemit.IsBlack())
		fprintf(stderr, "Lemit is black!");

	return (twoSided || Dot(intr.n, w) > 0) ? Spectrum(Lmap->Lookup(uv), SpectrumType::Illuminant) : Spectrum(0.f) ;
}

Spectrum DiffuseAreaLight::LScaled(const Interaction &intr, const Vector3f &w) const {
	if (Lemit.IsBlack())
		fprintf(stderr, "Lemit is black!");

	Spectrum L;
	if (hasTexture) {
		const Point2f uv = shape->GetUV(intr.p);
		L = Spectrum(Lmap->Lookup(uv), SpectrumType::Illuminant) / Lemit;
	} else {
		L = Lemit;
	}

	// Float Lmax = L.MaxComponentValue();
	// Spectrum LScaled = L / Lmax;
	Spectrum LScaled = L;

	return (twoSided || Dot(intr.n, w) > 0) ? LScaled : Spectrum(0.f);
}

Spectrum DiffuseAreaLight::Sample_Li(const Interaction &ref, const Point2f &u,
		Vector3f *wi, Float *pdf, VisibilityTester *vis) const {
	
	if (hasTexture)
		return Sample_Texture_Li(ref, u, wi, pdf, vis);
		
	ProfilePhase _(Prof::LightSample);
	Interaction pShape = shape->Sample(ref, u, pdf);
	pShape.mediumInterface = mediumInterface;

	if (*pdf == 0 || (pShape.p - ref.p).LengthSquared() == 0) {
		*pdf = 0;
		return Spectrum(0.f);
	}
	*wi = Normalize(pShape.p - ref.p);
	*vis = VisibilityTester(ref, pShape);
	return L(pShape, -*wi);
}

Spectrum DiffuseAreaLight::Sample_Texture_Li(const Interaction &ref, const Point2f &u,
		Vector3f *wi, Float *pdf, VisibilityTester *vis) const {
	ProfilePhase _(Prof::LightSample);

	// Find (u,v) sample coordinates (importance sample texels)
	Float uvPdf;
	Point2f uv = distribution->SampleContinuous(u, &uvPdf);

	// Naive for comparison, not importance sampled
	// Point2f uv = u;
	// Float uvPdf = 1.; // choice of (u,v) uniform over [0,1]x[0,1]

	// Convert sample point to direction
	Interaction pShape = shape->SampleUV(uv);
	pShape.mediumInterface = mediumInterface;

	*wi = Normalize(pShape.p - ref.p);
	*vis = VisibilityTester(ref, pShape);

	// Convert from area measure, as returned by the Sample() call
	// above, to solid angle measure.
	*pdf = uvPdf * 1./area * DistanceSquared(ref.p, pShape.p) / AbsDot(pShape.n, -*wi);

	if (*pdf == 0 || std::isinf(*pdf) || wi->LengthSquared() == 0) {
		*pdf = 0;
		return Spectrum(0.f);
	}

	return L(pShape, -*wi, uv);
}

Float DiffuseAreaLight::Pdf_Li(const Interaction &ref,
		const Vector3f &wi) const {

	if (hasTexture) 
		return Pdf_Li_Texture(ref, wi);

	ProfilePhase _(Prof::LightPdf);
	return shape->Pdf(ref, wi);
}

Float DiffuseAreaLight::Pdf_Li_Texture(const Interaction &ref,
		const Vector3f &wi) const {
	ProfilePhase _(Prof::LightPdf);
	SurfaceInteraction shapeRef;
	Float shapePdf = shape->Pdf(shapeRef, ref, wi);

	Point2f uv = shape->GetUV(shapeRef.p);
	Float texPdf = distribution->Pdf(uv);

	return texPdf * shapePdf;
}

Spectrum DiffuseAreaLight::Sample_Li_SA(const Interaction &ref, const Point2f &u,
		Vector3f *wi, Float *pdf, VisibilityTester *vis) const {
	ProfilePhase _(Prof::LightSample);
	Interaction pShape = shape->SampleSolidAngle(ref, u, pdf);
	pShape.mediumInterface = mediumInterface;

	if (*pdf == 0 || (pShape.p - ref.p).LengthSquared() == 0) {
		*pdf = 0;
		return 0.f;
	}
	*wi = Normalize(pShape.p - ref.p);
	*vis = VisibilityTester(ref, pShape);
	return L(pShape, -*wi);
}

Float DiffuseAreaLight::Pdf_SA(const Interaction &ref, const Vector3f &wi) const {
	ProfilePhase _(Prof::LightPdf);
	return shape->PdfSolidAngle(ref, wi);
}

Spectrum DiffuseAreaLight::Sample_Line_Li(const Interaction &ref,
		const Point2f &u, linesampling::LineInteraction &interaction,
		bool importanceSample) const {
	interaction = shape->LineSample(u, importanceSample);
	interaction.mediumInterface = mediumInterface;
	return Lemit;
}
Float DiffuseAreaLight::PdfLine(const Interaction &ref, const Vector3f &wi, Point3f &pLight,
			LineInteraction &lineInteraction, const Float& lineAngle, bool importanceSample) const {
	return shape->PdfLine(ref, wi, pLight, lineInteraction, lineAngle, importanceSample);
}

Spectrum DiffuseAreaLight::Sample_SphQuadLine_Li(const Point3f& p, const Float u, int strategy,
		LineInteraction &interaction) const {
	interaction = shape->SampleSQLine(p, u, strategy);
	interaction.mediumInterface = mediumInterface;
	return Lemit;
}

Spectrum DiffuseAreaLight::Sample_Le(const Point2f &u1, const Point2f &u2,
		Float time, Ray *ray, Normal3f *nLight, Float *pdfPos,
		Float *pdfDir) const {
	ProfilePhase _(Prof::LightSample);
// Sample a point on the area light's _Shape_, _pShape_
	Interaction pShape = shape->Sample(u1, pdfPos);
	pShape.mediumInterface = mediumInterface;
	*nLight = pShape.n;

// Sample a cosine-weighted outgoing direction _w_ for area light
	Vector3f w;
	if (twoSided) {
		Point2f u = u2;
		// Choose a side to sample and then remap u[0] to [0,1] before
		// applying cosine-weighted hemisphere sampling for the chosen side.
		if (u[0] < .5) {
			u[0] = std::min(u[0] * 2, OneMinusEpsilon);
			w = CosineSampleHemisphere(u);
		} else {
			u[0] = std::min((u[0] - .5f) * 2, OneMinusEpsilon);
			w = CosineSampleHemisphere(u);
			w.z *= -1;
		}
		*pdfDir = 0.5f * CosineHemispherePdf(std::abs(w.z));
	} else {
		w = CosineSampleHemisphere(u2);
		*pdfDir = CosineHemispherePdf(w.z);
	}

	Vector3f v1, v2, n(pShape.n);
	CoordinateSystem(n, &v1, &v2);
	w = w.x * v1 + w.y * v2 + w.z * n;
	*ray = pShape.SpawnRay(w);
	return L(pShape, w);
}

void DiffuseAreaLight::Pdf_Le(const Ray &ray, const Normal3f &n, Float *pdfPos,
		Float *pdfDir) const {
	ProfilePhase _(Prof::LightPdf);
	Interaction it(ray.o, n, Vector3f(), Vector3f(n), ray.time,
			mediumInterface);
	*pdfPos = shape->Pdf(it);
	*pdfDir =
			twoSided ?
					(.5 * CosineHemispherePdf(AbsDot(n, ray.d))) :
					CosineHemispherePdf(Dot(n, ray.d));
}

bool DiffuseAreaLight::CanEvaluateAnalytically(const SurfaceInteraction& isect,
		Spectrum * spectrum) const {
	bool result = shape->CanEvaluateAnalytically(isect, spectrum);
	if (spectrum)
		*spectrum *= Lemit;
	return result;
}

bool DiffuseAreaLight::CanIlluminate(const Interaction& isect) const {
	return shape->CanIlluminate(isect);
}

std::shared_ptr<AreaLight> CreateDiffuseAreaLight(const Transform &light2world,
		const Medium *medium, const ParamSet &paramSet,
		const std::shared_ptr<Shape> &shape) {
	Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
	Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
	int xSamples = paramSet.FindOneInt("xsamples", 1);
	int ySamples = paramSet.FindOneInt("ysamples", 1);
	bool twoSided = paramSet.FindOneBool("twosided", false);
	std::string texmap = paramSet.FindOneFilename("mapname", "");

	if (PbrtOptions.quickRender) {
		xSamples = std::max(1, xSamples / 4);
		ySamples = std::max(1, ySamples / 4);
	}
	return std::make_shared<DiffuseAreaLight>(light2world, medium, L * sc,
			xSamples, ySamples, shape, texmap, twoSided);
}

}  // namespace pbrt
