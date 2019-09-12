/*
 pbrt source code is Copyright(c) 1998-2015
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

#include "directlinelighting.h"

#include "stats.h"
#include "paramset.h"
#include "interval.h"
#include "lineinteraction.h"
#include "camera.h"

using linesampling::Line3f;
using linesampling::LineInteraction;
using linesampling::IntervalCollection;
using linesampling::LineEvaluation;

namespace pbrt {

// DirectLightingIntegrator Method Definitions
void DirectLineLightingIntegrator::Preprocess(const Scene &scene,
		Sampler &sampler) {
	// Compute number of samples to use for each light
	for (const auto &light : scene.lights) {
		nLightSamples.push_back(
				std::make_pair(light->xSamples, light->ySamples));
	}

	// Request samples for sampling all lights
	for (int i = 0; i < maxDepth; ++i) {
		for (size_t j = 0; j < scene.lights.size(); ++j) {
			const std::pair<int, int>& lightSample = nLightSamples[j];

			if (strategy == LineSampleStrategy::UNIFORM) {
				// KN: for each 2d sample, one number is random line direction, other number is random offset
				sampler.Request2DArray(lightSample);
			} else if (strategy == LineSampleStrategy::FIXED) {
				// KN: each number is random offset
				sampler.Request1DArray(lightSample.first * lightSample.second);
			} else if (strategy == LineSampleStrategy::PARALLEL) {
				// KN: each number is random offset
				// sampler.Request2DArray(lightSample.first * lightSample.second);
				sampler.Request2DPixelArray(lightSample);
			}
		}
	}

	// Generate line directions
	if (strategy == LineSampleStrategy::FIXED
			|| strategy == LineSampleStrategy::FIXED_RANDOM) {
		for (size_t j = 0; j < scene.lights.size(); ++j) {
			if (strategy == LineSampleStrategy::FIXED) {
				lightDirections.push_back(fixedDirection);
			} else {
				lightDirections.push_back(sampler.UniformFloat());
			}
		}
	}
}

// Integrator Utility Functions
Spectrum DirectLineLightingIntegrator::UniformSampleAllLights(
		const Interaction &it, const Scene &scene, MemoryArena &arena,
		Sampler &sampler, const std::vector<std::pair<int, int>> &nLightSamples,
		bool handleMedia) const {
	ProfilePhase p(Prof::DirectLighting);
	Spectrum L(0.f);

	for (size_t j = 0; j < scene.lights.size(); ++j) {
		// Accumulate contribution of j'th light to L
		const std::shared_ptr<Light> &light = scene.lights[j];
		const std::pair<int, int>& lightSample = nLightSamples[j];
		const int nSamples = lightSample.first * lightSample.second;
		const Float invSamples = 1.0 / Float(nSamples);

		if (strategy == LineSampleStrategy::FIXED
				|| strategy == LineSampleStrategy::FIXED_RANDOM) {
			const Float lineDirection = lightDirections[j];
			const Float *lineOffsets = sampler.Get1DArray(nSamples);
			// KN: running with strategy "FIXED_RANDOM" results in the following error
			if (!lineOffsets)
				Error("DirectLineLightingIntegrator::UniformSampleAllLights::"
						"lineOffsets is nullptr!");
			if (!light->CanIlluminate(it))
				continue;

			// KN: for each random offset u, estimate direct lighting contribution
			for (int k = 0; k < nSamples; ++k)
				L += this->EstimateDirect(it, *light, lineDirection,
						lineOffsets[k], scene, sampler, arena, handleMedia)
						* invSamples;
		} else if (strategy == LineSampleStrategy::UNIFORM) {
			const Point2f* lineSamples = sampler.Get2DArray(lightSample.first,
					lightSample.second);
			if (!lineSamples)
				Error("DirectLineLightingIntegrator::UniformSampleAllLights::"
						"lineSamplers is nullptr!");
			if (!light->CanIlluminate(it))
				continue;

			for (int k = 0; k < nSamples; ++k)
				// KN: first entry in lineSamples acts as line direction, second entry acts as random offset u in [a,b]
				L += this->EstimateDirect(it, *light, lineSamples[k].x,
						lineSamples[k].y, scene, sampler, arena, handleMedia)
						* invSamples;
		} else if (strategy == LineSampleStrategy::PARALLEL) {
			// KN: get a random line direction that will be constant for this shading point only
			// const Float lineDirection = sampler.Get1DPixelSample();

			// KN: changed to use line below for line direction since line above always returned zero
			const Float lineDirection = sampler.UniformFloat();
			// const Float *lineOffsets = sampler.Get1DArray(nSamples);
			const Point2f *lineOffsets = sampler.Get2DPixelArray(lightSample); // only using first random number of pair for line offset

			if (!lineOffsets)
				Error("DirectLineLightingIntegrator::UniformSampleAllLights::"
						"lineOffsets is nullptr!");
			if (!light->CanIlluminate(it))
				continue;
			for (int k = 0; k < nSamples; ++k)
				L += this->EstimateDirect(it, *light, lineDirection,
						lineOffsets[k].x, scene, sampler, arena, handleMedia)
						* invSamples;
		} else {
			Error("Unimplemented strategy!");
		}
	}
	return L;
}

Spectrum DirectLineLightingIntegrator::Li(const RayDifferential &ray,
		const Scene &scene, Sampler &sampler, MemoryArena &arena,
		int depth) const {
	ProfilePhase p(Prof::SamplerIntegratorLi);
	Spectrum L(0.f);

// Find closest ray intersection or return background radiance
	SurfaceInteraction isect;
	if (!scene.Intersect(ray, &isect)) {
		if (method != RenderMethod::LINESAMPLECOUNT)
			for (const auto &light : scene.lights)
				L += light->Le(ray);
		return L;
	}

// Compute scattering functions for surface interaction
	isect.ComputeScatteringFunctions(ray, arena);

	if (!isect.bsdf) {
		if (method == RenderMethod::LINESAMPLECOUNT)
			return L;
		else
			return Li(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
	}

	const Vector3f& wo = isect.wo;

	// Compute emitted light if ray hit an area light source
	if (method == RenderMethod::RENDER) {
		if (showLineDirection) {
			L += isect.LeLines(isect.wo, fixedDirection);
		} else {
			L += isect.Le(isect.wo);
		}
	}

	if (scene.lights.size() > 0) {
		// Compute direct lighting for _DirectLightingIntegrator_ integrator
		L += UniformSampleAllLights(isect, scene, arena, sampler,
				nLightSamples);
	}
	if (depth + 1 < maxDepth) {
		Vector3f wi;
		// Trace rays for specular reflection and rethisfraction
		L += SpecularReflect(ray, isect, scene, sampler, arena, depth);
		L += SpecularTransmit(ray, isect, scene, sampler, arena, depth);
	}
	return L;
}

Spectrum DirectLineLightingIntegrator::EstimateDirect(const Interaction &it,
		const Light &light, Float lightDirection, Float lightPos,
		const Scene &scene, Sampler &sampler, MemoryArena &arena,
		bool handleMedia, bool specular) const {

	/**************************************************************************
	 * Sample the light source
	 **************************************************************************/

	linesampling::LineInteraction lineInteraction;
	const Point2f pLight(lightDirection, lightPos);
	const Spectrum& Li = light.Sample_Line_Li(it, pLight, lineInteraction,
			importanceSample);

	// exit if no contribution, or when the line sample is too small.
	if (Li.IsBlack() || lineInteraction.length < EPSILON)
		return Spectrum(0.f);

	/**************************************************************************
	 * Retrieve the intersection information
	 **************************************************************************/

	// retrieve the interaction
	const SurfaceInteraction &isect = (const SurfaceInteraction &) it;

	// check whether the normals have to be inverted.
	const Point3f &p = isect.p;
	const Normal3f &n = isect.n;
	const Normal3f &sn = isect.shading.n;
	const Vector3f &reflection = Reflect(isect.wo, sn);

	/**********************************************************************
	 * Clip the line against the plane spanned by the real surface normal *
	 **********************************************************************/

	const Line3f& sampleLine = lineInteraction.line;
	const Normal3f &ln = lineInteraction.n;
	const Vector3f &wo = isect.wo;

	Point3f p1 = sampleLine.p1;
	Point3f p2 = sampleLine.p2;
	Vector3f e1 = p1 - p;
	Vector3f e2 = p2 - p;
	Vector3f d = p2 - p1;

	// check whether the point is completely beneath the line sample
	if (Dot(p - p1, ln) <= 0)
		return Spectrum(0.f);

	// exit when the line is under the plane of the real surface normal
	if (Dot(e1, n) <= 0 && Dot(e2, n) <= 0)
		return Spectrum(0.f);

	// intersect the line with the plane spanned by the normal and the intersection point
	const Float s = -Dot(n, e1) / Dot(n, d);
	if (s >= -EPSILON && s <= 1.0 + EPSILON) {
		if (Dot(n, e1) > EPSILON) {
			// line.p1 is on the correct side, line.p2 is under the plane
			p2 = p1 + d * (s - EPSILON);
			e2 = p2 - p;

			// if line.p1 is on the correct side and s is smaller than
			// epsilon, nothing of the line segment remains
			if (s < EPSILON)
				return Spectrum(0.f);
		} else if (Dot(n, e2) > 0) {
			// line.p2 is on the correct side, line.p1 is under the plane
			p1 = p1 + d * (s + EPSILON);
			e1 = p1 - p;

			// if line.p2 is on the correct side and s is larger than
			// 1 - epsilon, nothing of the line segment remains
			if (s > 1.0 - EPSILON)
				return Spectrum(0.f);
		} else
			return Spectrum(0.f);
		d = p2 - p1;
	}

	// perform the intersection test
	IntervalCollection collection(0, 1);
	{
		// add a small epsilon to the line sample to avoid self intersections
		const Float offset = EPSILON;
		const Vector3f& isectLineOffset = Vector3f(ln * offset);
		const Point3f& isectPoint = p + Vector3f(offset * n);
		const Point3f& isectLineOrigin = p1 + isectLineOffset;
		const Point3f& isectLineDestination = p2 + isectLineOffset;

		// construct the query
		linesampling::LineVisibilityQuery3f query(isectPoint, isectLineOrigin,
				isectLineDestination);

		// perform the intersection test
		// If line segment is completely obscured by primitives (collection of line intervals is empty), then exit
		if (scene.IntersectL(query, n, it.time, collection)
				&& collection.empty())
			return Spectrum(0.f);
	}

	if (method == RenderMethod::LINESAMPLECOUNT)
		return Spectrum(collection.size());

	// split the interval containing the shading normal intersection
	if (Dot(n, sn) < OneMinusEpsilon) {
		const Float t = -Dot(sn, e1) / Dot(sn, d);
		collection.split(t);
	}

	/**************************************************************************
	 * Transform the vectors to a two dimensional plane for easier evaluation
	 **************************************************************************/

	const LineEvaluation e(n, sn, ln, isect.wo, e1, e2);
	const Spectrum& f = isect.bsdf->f(e, collection);

	return f * Li / lineInteraction.pdf;
}

DirectLineLightingIntegrator *CreateDirectLineLightingIntegrator(
		const ParamSet &params, std::shared_ptr<Sampler> sampler,
		std::shared_ptr<const Camera> camera) {
	int maxDepth = params.FindOneInt("maxdepth", 1);

	int np;
	const int *pb = params.FindInt("pixelbounds", &np);
	Bounds2i pixelBounds = camera->film->GetSampleBounds();
	if (pb) {
		if (np != 4)
			Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
					np);
		else {
			pixelBounds = Intersect(pixelBounds, Bounds2i { { pb[0], pb[2] }, {
					pb[1], pb[3] } });
			if (pixelBounds.Area() == 0)
				Error("Degenerate \"pixelbounds\" specified.");
		}
	}

	// max seconds
	float maxSeconds = params.FindOneFloat("max-seconds", -1.0);

	bool importanceSample = params.FindOneBool("importance", true);
	bool showLineDirection = params.FindOneBool("showdirection", false);

	/***************************************************************************
	 * Find the render method
	 **************************************************************************/

	const std::string& methodString = params.FindOneString("method", "render");
	DirectLineLightingIntegrator::RenderMethod method;
	if (methodString == "lines" || methodString == "linecount") {
		method = DirectLineLightingIntegrator::RenderMethod::LINESAMPLECOUNT;
	} else {
		if (methodString != "render") {
			Error("Unknown render method string '%s' for "
					"DirectLineLightingIntegrator specified!\nUsing "
					"'render' instead!", methodString.c_str());
		}
		method = DirectLineLightingIntegrator::RenderMethod::RENDER;
	}

	/***************************************************************************
	 * Find the line sampling strategy
	 **************************************************************************/

	const std::string& lineStrategyString = params.FindOneString("linestrategy",
			params.FindOneString("strategy", "parallel"));
	DirectLineLightingIntegrator::LineSampleStrategy strategy;

	// std::cout << "line strategy: " << lineStrategyString << "\n";

	Float fixedDirection = Infinity;

	if (lineStrategyString == "uniform") {
		strategy = DirectLineLightingIntegrator::LineSampleStrategy::UNIFORM;
	} else if (lineStrategyString == "fixed") {
		fixedDirection = params.FindOneFloat("fixeddirection", 0);
		strategy = DirectLineLightingIntegrator::LineSampleStrategy::FIXED;
	} else if (lineStrategyString == "fixedrandom") {
		strategy =
				DirectLineLightingIntegrator::LineSampleStrategy::FIXED_RANDOM;
	} else {
		if (lineStrategyString != "parallel") {
			Error("Unknown render method string '%s' for "
					"DirectLineLightingIntegrator specified!\nUsing "
					"'parallel' instead!", methodString.c_str());
		}
		strategy = DirectLineLightingIntegrator::LineSampleStrategy::PARALLEL;
	}

	return new DirectLineLightingIntegrator(maxDepth, camera, sampler,
			pixelBounds, method, strategy, importanceSample, showLineDirection, fixedDirection, maxSeconds);
}

}
