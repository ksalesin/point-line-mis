/*
 * group.h
 *
 *  Created on: Aug 14, 2017
 *      Author: niels
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef INTEGRATORS_GROUP_H_
#define INTEGRATORS_GROUP_H_

// integrators/directlighting.h*
#include "pbrt.h"
#include "integrator.h"
#include "scene.h"

namespace pbrt {

// GroupIntegrator Declarations
class GroupIntegrator: public SamplerIntegrator {
public:
	// GroupIntegrator Public Methods
	GroupIntegrator(std::shared_ptr<const Camera> camera,
			std::shared_ptr<Sampler> sampler, const Bounds2i &pixelBounds,
			int groupSize = 1) :
			SamplerIntegrator(camera, sampler, pixelBounds), groupSize(
					groupSize) {
	}

	Spectrum Li(const RayDifferential &ray, const Scene &scene,
			Sampler &sampler, MemoryArena &arena, int depth) const;
	void Preprocess(const Scene &scene, Sampler &sampler) override {
	}

protected:
	const int groupSize;
};

GroupIntegrator *CreateGroupIntegrator(const ParamSet &params,
		std::shared_ptr<Sampler> sampler, std::shared_ptr<const Camera> camera);

} // namespace pbrt

#endif /* INTEGRATORS_GROUP_H_ */
