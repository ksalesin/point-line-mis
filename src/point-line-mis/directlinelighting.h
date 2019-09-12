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

/* 
   Option to visualize line direction on light source (for fixed strategy) added by Kate
*/

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_LINESAMPLING_DIRECTLINELIGHTING_H
#define PBRT_LINESAMPLING_DIRECTLINELIGHTING_H

// linesampling/directlinelighting.h*
#include "pbrt.h"
#include "integrator.h"
#include "scene.h"
#include "sampling.h"

namespace pbrt {

// DirectLightingIntegrator Declarations
class DirectLineLightingIntegrator: public SamplerIntegrator {
public:
	enum class RenderMethod {
		RENDER, LINESAMPLECOUNT
	};

	enum class LineSampleStrategy {
		FIXED, /* fixed direction for all light sources, same for all shading points (input in scene file with "fixeddirection") */
		FIXED_RANDOM, /* fixed direction but randomly sampled for each light (instead of input in scene file) */
		UNIFORM, /* random direction for every light sample */
		PARALLEL /* parallel direction for each shading point, randomly initialized $x$ */
	};

	// DirectLightingIntegrator Public Methods
	DirectLineLightingIntegrator(int maxDepth,
			std::shared_ptr<const Camera> camera,
			std::shared_ptr<Sampler> sampler, 
			const Bounds2i &pixelBounds,
			const RenderMethod& method, 
			const LineSampleStrategy& strategy,
			bool importanceSample, 
			bool showLineDirection, 
			Float fixedDirection,
			Float maxSeconds) :
			SamplerIntegrator(camera, sampler, pixelBounds, maxSeconds), maxDepth(maxDepth), method(
					method), strategy(strategy), importanceSample(importanceSample), 
					showLineDirection(showLineDirection), fixedDirection(fixedDirection) {
	}

	void Preprocess(const Scene &scene, Sampler &sampler) override;

	Spectrum Li(const RayDifferential &ray, const Scene &scene,
			Sampler &sampler, MemoryArena &arena, int depth) const override;

	Spectrum UniformSampleAllLights(const Interaction &it, const Scene &scene,
			MemoryArena &arena, Sampler &sampler,
			const std::vector<std::pair<int, int>> &nLightSamples,
			bool handleMedia = false) const override;

	Spectrum EstimateDirect(const Interaction &it, const Light &light,
			Float lightDir, Float LightOffset, const Scene &scene,
			Sampler &sampler, MemoryArena &arena, bool handleMedia = false,
			bool specular = false) const;
private:
	// DirectLightingIntegrator Private Data
	const int maxDepth;
	const bool importanceSample;
	const bool showLineDirection;
	const RenderMethod method;
	const LineSampleStrategy strategy;
	std::vector<std::pair<int, int>> nLightSamples;

	/*
	 * this vector is only used when the LineSampleStrategy is fixed!
	 * during the preprocess a fixed random direction is chosen for each light
	 * source, which is reused for each shading point.
	 */
	std::vector<Float> lightDirections;
	Float fixedDirection; // only used with FIXED strategy
};

DirectLineLightingIntegrator *CreateDirectLineLightingIntegrator(
		const ParamSet &params, std::shared_ptr<Sampler> sampler,
		std::shared_ptr<const Camera> camera);

}

#endif  // PBRT_LINESAMPLING_DIRECTLINELIGHTING_H
