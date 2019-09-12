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
   File created by Kate on 3/27/2019.
*/

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_LINESAMPLING_POINTLINEDIRECT_H
#define PBRT_LINESAMPLING_POINTLINEDIRECT_H

#include "pbrt.h"
#include "integrator.h"
#include "scene.h"
#include "sampling.h"

namespace pbrt {

// PointLineDirectIntegratpr Declarations
class PointLineDirectIntegrator: public SamplerIntegrator {
public:
	enum class RenderMethod {
		RENDER,           // normal
        LINESAMPLECOUNT,  // render number of lines into which line segment is cut after visibility estimation 
        DRAWLINES,        // render normally, but also draw fixed line direction on light source
        INTVIZ         	  // visualize integrand for given sampling scheme at given scene intersection point
	};

	enum class SampleStrategy {
		LINE_SURF_AREA,     // line sampling, uniform over surface area of light
		LINE_SOLID_ANGLE,   // line sampling, uniform over solid angle of light
		PT_SURF_AREA,     	// point sampling, uniform over surface area of light
		PT_SOLID_ANGLE,    	// point sampling, uniform over solid angle of light
		PT_BRDF,       		// point sampling, uses PBRT's scheme for sampling a shade point's BRDF
		NONE 				// placeholder for no strategy
	};

	// for progress reporter, used for research images
	// static long nSceneTraversals;

	// PointLineDirectIntegrator Public Methods
	PointLineDirectIntegrator(int maxDepth,
			std::shared_ptr<const Camera> camera,
			std::shared_ptr<Sampler> sampler, 
			const Bounds2i &pixelBounds,
			const RenderMethod& method, 
			const SampleStrategy& strategy1, const SampleStrategy& strategy2, const SampleStrategy& strategy3, 
			int nStrategies, bool stratified, bool mis,
			Float lineAngle1, Float lineAngle2, Float lineAngle3,
			bool lineIS1, bool lineIS2, bool lineIS3, Float smoothstepWidth, Float ratio,
			float maxSeconds, long maxSceneTraversals) :
			SamplerIntegrator(camera, sampler, pixelBounds, maxSeconds, maxSceneTraversals), 
				_maxDepth(maxDepth), _method(method), 
				_strategy1(strategy1), _strategy2(strategy2), _strategy3(strategy3), _nStrategies(nStrategies),
				_lineAngle1(lineAngle1), _lineAngle2(lineAngle2), _lineAngle3(lineAngle3), 
				_lineIS1(lineIS1), _lineIS2(lineIS2), _lineIS3(lineIS3), 
				_smoothstepWidth(smoothstepWidth), _ratio(ratio), _stratified(stratified), _mis(mis) {}

	void Preprocess(const Scene &scene, Sampler &sampler) override;

	Spectrum Li(const RayDifferential &ray, const Scene &scene,
			Sampler &sampler, MemoryArena &arena, int depth) const override;

	Spectrum UniformSampleAllLights(const Interaction &it, const Scene &scene,
			MemoryArena &arena, Sampler &sampler,
			const std::vector<std::pair<int, int>> &nLightSamples,
			bool handleMedia = false) const override;

	Spectrum UniformSampleOneLight(const Interaction &it,
			const Scene &scene, MemoryArena &arena, Sampler &sampler,
			const std::vector<std::pair<int, int>> &nLightSamples,
			bool handleMedia = false, const Distribution1D *lightDistrib =
					nullptr) const;

	// Overrides base class, does nothing
	Spectrum EstimateDirect(const Interaction &it, const Light &light,
			const Point2f &uLight, const Scene &scene,
			Sampler &sampler, MemoryArena &arena, bool handleMedia = false,
			bool specular = false) const override { return Spectrum(0.f); }

	// Estimates direct lighting for all sample strategies, MIS if >1 strategy
	Spectrum EstimateDirect(const Interaction &it, const Light &light,
			Point2f &u, int strategy, int nSamples1, int nSamples2, int nSamples3, 
			const Scene &scene, Sampler &sampler, MemoryArena &arena, 
			bool handleMedia = false, bool specular = false) const;

	// Routes traffic to correct estimate method given sampling strategies
	Spectrum EstimateIntermediate(const Interaction &it, const Light &light,
			const Point2f &uLight, const Scene &scene, Sampler &sampler, MemoryArena &arena, const int strategy,
			Vector3f *wi, Float *pdf, Float *smoothstep, bool handleMedia = false, bool specular = false) const;

	// Estimates direct lighting by sampling a line on light source, then sampling point within unoccluded regions of line
	Spectrum EstimateDirectLine(const Interaction &it, const Light &light,
			const Point2f& uLight, const int strategy, const Scene &scene, Sampler &sampler, MemoryArena &arena, 
			Vector3f *wi, Float *pdf, Float *smoothstep, bool handleMedia = false, bool specular = false) const;
	
	// Routes traffic to correct pdf method given sampling strategies
	Float PdfIntermediate(const Interaction &it, const Light &light, const Scene &scene,
			Sampler &sampler, MemoryArena &arena, const int strategy, Vector3f& wi, Float *smoothstep,
			bool handleMedia = false, bool specular = false) const;
	
	Float PdfLine(const Interaction &it, const Light &light, const Scene &scene,
			Sampler &sampler, MemoryArena &arena, const int strategy, Vector3f& wi, Float *smoothstep,
			bool handleMedia = false, bool specular = false) const;

private:
	// Smoothstep a point p that lies between two other points, start and end
	Float smoothstepCubic(Float x) const;
	Float smoothstepEaseIn(Float x) const;
	Float smoothstepPoints(Point3f& p, Point3f& start, Point3f& end) const;
	Float smoothstepPoints(Point3f& p, Point3f& start, Point3f& end, Float occDistWorld1, Float occDistWorld2) const;

	std::vector<std::pair<int, int>> _nLightSamples;

	const int _maxDepth;
	const RenderMethod _method;

	const Float _smoothstepWidth; // width of smoothstep radius of influence
	const bool _stratified; 	  // whether samples should be split along domain and remapped or used as is
	const bool _mis; 			  // whether to MIS (true) or average (false) for multiple strategies
	const Float _ratio; 		  // ratio of point samples to line samples for MIS/avg

	int _nStrategies;
	const SampleStrategy _strategy1;
	const SampleStrategy _strategy2;
	const SampleStrategy _strategy3;

	// line directions on light source
	Float _lineAngle1;
	Float _lineAngle2;
	Float _lineAngle3;

	// whether to importance sample lines by length
	const bool _lineIS1;
	const bool _lineIS2;
	const bool _lineIS3;
};

PointLineDirectIntegrator *CreatePointLineDirectIntegrator(
		const ParamSet &params, std::shared_ptr<Sampler> sampler,
		std::shared_ptr<const Camera> camera);

}

#endif  // PBRT_LINESAMPLING_POINTLINEDIRECT_H
