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

#ifndef PBRT_CORE_SCENE_H
#define PBRT_CORE_SCENE_H

// core/scene.h*
#include "pbrt.h"
#include "primitive.h"
#include "integrator.h"

namespace pbrt {

// Scene Declarations
class Scene {
public:
	// Scene Public Methods
	Scene(std::shared_ptr<const Primitive> aggregate,
			const std::vector<std::shared_ptr<Light>> &lights,
			std::shared_ptr<const Primitive> inside = nullptr,
			std::shared_ptr<const Primitive> outside = nullptr,
			std::shared_ptr<const Primitive> mixed = nullptr) :
			lights(lights), aggregate(aggregate) {
		// Scene Constructor Implementation
		worldBound = aggregate->WorldBound();
		for (const auto &light : lights) {
			light->Preprocess(*this);
			if (light->flags & (int) LightFlags::Infinite)
				infiniteLights.push_back(light);
		}

		proxies[ProxyType::INSIDE] = inside;
		proxies[ProxyType::OUTSIDE] = outside;
		proxies[ProxyType::MIXED] = mixed;
		proxyBound = worldBound;

		for (int i = 0; i < 3; ++i)
			if (proxies[i])
				proxyBound = Union(proxyBound, proxies[i]->WorldBound());

	}
	const Bounds3f &WorldBound() const {
		return worldBound;
	}
	const Bounds3f &ProxyBound() const {
		return proxyBound;
	}

	bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
	bool IntersectP(const Ray &ray) const;
	bool IntersectP(const Ray &ray, uint32_t ignoreId) const;
	bool IntersectP(const Ray &ray, const ProxyType& type,
			uint32_t ignoreId) const;
	bool IntersectL(const LineVisibilityQuery3f &query, const Normal3f &normal,
			Float time, IntervalCollection&) const;

	bool IntersectAll(const Ray& r, GroupIsects&, int groupSize = 1) const;
	bool IntersectTr(Ray ray, Sampler &sampler, SurfaceInteraction *isect,
			Spectrum *transmittance) const;
	Float Area() const;

	// Scene Public Data
	std::vector<std::shared_ptr<Light>> lights;
	// Store infinite light sources separately for cases where we only want
	// to loop over them.
	std::vector<std::shared_ptr<Light>> infiniteLights;

private:
	// Scene Private Data
	std::shared_ptr<const Primitive> aggregate;
	std::shared_ptr<const Primitive> proxies[3];
	Bounds3f worldBound;
	Bounds3f proxyBound;
};

}  // namespace pbrt

#endif  // PBRT_CORE_SCENE_H
