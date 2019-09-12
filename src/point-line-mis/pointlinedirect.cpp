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
   Option to visualize line direction on light source (for fixed strategy) added by Kate (12/2018)
*/

#include "pointlinedirect.h"

#include "stats.h"
#include "paramset.h"
#include "interval.h"
#include "lineinteraction.h"
#include "camera.h"

using linesampling::Line3f;
using linesampling::LineInteraction;
using linesampling::IntervalCollection;
using linesampling::Interval;
using linesampling::IntervalConstIterator;
using linesampling::LineEvaluation;

namespace pbrt {

STAT_COUNTER("Scene traversals", nSceneTraversals);
// long PointLineDirectIntegrator::nSceneTraversals = 0;

Float PointLineDirectIntegrator::smoothstepEaseIn(Float x) const {
    if (x <= 0)
        return 0.;

    if (x >= 1)
        return 1.;

    return -(x-2) * x;
}

Float PointLineDirectIntegrator::smoothstepCubic(Float x) const {
    if (x <= 0)
        return 0.;

    if (x >= 1)
        return 1.;

    Float x2 = x * x;
    return 3*x2 - 2*x*x2;
}

Float PointLineDirectIntegrator::smoothstepPoints(Point3f& p, Point3f& start, Point3f& end) const {
	// Get smoothstep distance (normalized over half distance from start to end)
	Float dmax = (start-end).Length()/2.;
    Float d = std::min((p-start).Length(), (p-end).Length());
	// Float ratio = d / dmax; // normalized distance
	Float ratio = d; // global distance

	if (dmax < 1e-10 || ratio < 1e-10)
		return 0.; // avoid nans

    Float x = ratio/_smoothstepWidth;

	if (isnan(x))
		return 0.;
		
	// use either cubic or ease-in
    return smoothstepEaseIn(x);
}


Float PointLineDirectIntegrator::smoothstepPoints(Point3f& p, Point3f& start, Point3f& end, Float occDistWorld1, Float occDistWorld2) const {
	
	// global distance
	Float d1 = (p-start).Length();
	Float d2 = (p-end).Length();

	Float o1 = occDistWorld1/_smoothstepWidth;
	Float o2 = occDistWorld2/_smoothstepWidth;

	Float sso1 = smoothstepEaseIn(o1);
	Float sso2 = smoothstepEaseIn(o2);
	
	// smoothstep width adjusted by closeness to occluded zone
	Float ssw_new1 = _smoothstepWidth * sso1;
	Float ssw_new2 = _smoothstepWidth * sso2;
	
	Float x1 = d1/ssw_new1;
	Float x2 = d2/ssw_new2;
	
	return smoothstepEaseIn(x1) * smoothstepEaseIn(x2); // multiply
	// return std::min(smoothstepEaseIn(x1), smoothstepEaseIn(x2)); // min
}

// PointLineDirectIntegrator Method Definitions
void PointLineDirectIntegrator::Preprocess(const Scene &scene,
		Sampler &sampler) {

	// Compute number of samples to use for each light
	for (const auto &light : scene.lights) {
		if (light->xSamples != light->ySamples)
			Error("Sample count X must equal sample count Y.");
		
		_nLightSamples.push_back(std::make_pair(light->xSamples, light->ySamples));
	}

	// Request samples for sampling all lights
	for (int i = 0; i < _maxDepth; ++i) {
		for (size_t j = 0; j < scene.lights.size(); ++j) {
			std::pair<int,int> lightSamples = _nLightSamples[j];

			// If just one sample total, take one sample from each strategy
			if (lightSamples.first == 1 && lightSamples.second == 1) {
				for (int k = 0; k < _nStrategies; k++) {
					sampler.Request2DArray(lightSamples);
				}
				
			// Otherwise split up the samples and allocate an equal amount to each strategy
			} else {
				sampler.Request2DArray(lightSamples);
			}
		}
	}
}

// Integrator Utility Functions
Spectrum PointLineDirectIntegrator::Li(const RayDifferential &ray,
		const Scene &scene, Sampler &sampler, MemoryArena &arena,
		int depth) const {
	ProfilePhase p(Prof::SamplerIntegratorLi);
	Spectrum L(0.f);

    // Find closest ray intersection or return background radiance
	SurfaceInteraction isect;
	if (!scene.Intersect(ray, &isect)) {
		if (_method != RenderMethod::LINESAMPLECOUNT)
			for (const auto &light : scene.lights)
				L += light->Le(ray);
		return L;
	}

    // Compute scattering functions for surface interaction
	isect.ComputeScatteringFunctions(ray, arena);

	if (!isect.bsdf) {
		if (_method == RenderMethod::LINESAMPLECOUNT)
			return L;
		else
			return Li(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
	}

	const Vector3f& wo = isect.wo;

	// Compute emitted light if ray hit an area light source
	if (_method == RenderMethod::RENDER) {
        L += isect.LeScaled(isect.wo);
    } else if (_method == RenderMethod::DRAWLINES) {
        L += isect.LeLines(isect.wo, _lineAngle1);
    }

	if (scene.lights.size() > 0) {
		// Compute direct lighting for _PointLineDirectIntegrator_ integrator
		L += UniformSampleAllLights(isect, scene, arena, sampler, _nLightSamples);
	}

	if (depth + 1 < _maxDepth) {
		Vector3f wi;
		// Trace rays for specular reflection and refraction
		L += SpecularReflect(ray, isect, scene, sampler, arena, depth);
		L += SpecularTransmit(ray, isect, scene, sampler, arena, depth);
	}

	return L;
}

Spectrum PointLineDirectIntegrator::UniformSampleAllLights(
		const Interaction &it, const Scene &scene, MemoryArena &arena,
		Sampler &sampler, const std::vector<std::pair<int, int>> &nLightSamples,
		bool handleMedia) const {
	
	ProfilePhase p(Prof::DirectLighting);
	Spectrum L(0.f);

	for (size_t j = 0; j < scene.lights.size(); ++j) {
		// Accumulate contribution of j'th light to L
		const std::shared_ptr<Light> &light = scene.lights[j];
		const std::pair<int, int>& lightSample = _nLightSamples[j];
		const int nSamples = lightSample.first * lightSample.second;

        if (!light->CanIlluminate(it))
            continue;

		// Use Shirley remapping to split random samples between MIS strategies
		const Point2f *allSamples = sampler.Get2DArray(lightSample.first, lightSample.second);

		std::vector<Point2f> samples1, samples2, samples3;

		if (nSamples == 1) {
			samples1.push_back(allSamples[0]);

			if (_nStrategies > 1) {
				const Point2f *moreSamples = sampler.Get2DArray(lightSample.first, lightSample.second);
				samples2.push_back(moreSamples[0]);
			}

			if (_nStrategies > 2) {
				const Point2f *evenMoreSamples = sampler.Get2DArray(lightSample.first, lightSample.second);
				samples3.push_back(evenMoreSamples[0]);
			}

		} else {
			if (_nStrategies == 1) {
				for (int i = 0; i < nSamples; i++) {
					samples1.push_back(allSamples[i]);
				}

			} else if (_nStrategies == 2) {
				int halfSamples = nSamples/2;
				double sqrtSamples = std::sqrt(nSamples);

				if ((int)sqrtSamples % 2 != 0)
					Error("For remapping equally between 2 MIS strategies, must use some sample count N x N where N is divisible by 2.");
				
				if (_stratified) {
					// split samples along first dimension
					for (int i = 0; i < nSamples; i++) {
						Point2f s = allSamples[i];
						if (s.x < 0.5) {
							s.x = s.x * 2;
							samples1.push_back(s);
						} else {
							s.x = s.x * 2 - 1;
							samples2.push_back(s);
						}
					}
				} else {
					int nSamples1 = floor(_ratio * nSamples);
    				int nSamples2 = nSamples - nSamples1;

					for (int i = 0; i < nSamples1; i++) {
						samples1.push_back(allSamples[i]);
					}
					for (int i = nSamples1; i < nSamples1+nSamples2; i++)
						samples2.push_back(allSamples[i]);
				}

			} else {
				int thirdSamples = nSamples / 3;
				double sqrtSamples = std::sqrt(nSamples);

				if ((int)sqrtSamples % 3 != 0)
					Error("For remapping equally between 3 MIS strategies, must use some sample count N x N where N is divisible by 3.");
				
				if (_stratified) {
					// split samples along first dimension
					double oneThird = 1./3.;
					for (int i = 0; i < nSamples; i++) {
						Point2f s = allSamples[i];
						if (s.x < oneThird) {
							s.x = s.x * 3;
							samples1.push_back(s);
						} else if (s.x < 1-oneThird) {
							s.x = s.x * 3 - 1;
							samples2.push_back(s);
						} else {
							s.x = s.x * 3 - 2;
							samples3.push_back(s);
						}
					}
				} else {
					int nSamples1 = std::ceil(_ratio * nSamples);
    				int nSamples2 = std::ceil(_ratio * nSamples);
    				int nSamples3 = nSamples - nSamples1 - nSamples2;

					for (int i = 0; i < nSamples1; i++)
						samples1.push_back(allSamples[i]);
					for (int i = nSamples1; i < nSamples1+nSamples2; i++)
						samples2.push_back(allSamples[i]);
					for (int i = nSamples1+nSamples2; i < nSamples1+nSamples2+nSamples3; i++)
						samples3.push_back(allSamples[i]);
				}
			}
		}

		int nSamples1 = samples1.size(), nSamples2 = samples2.size(), nSamples3 = samples3.size();
		Float invSamples1 = 1./nSamples1, invSamples2 = 1./nSamples2, invSamples3 = 1./nSamples3;

		for (Point2f sample : samples1) {
			L += this->EstimateDirect(it, *light, sample, 1, nSamples1, nSamples2, nSamples3, scene, sampler, arena)
							* invSamples1;
		}
		for (Point2f sample : samples2) {
			L += this->EstimateDirect(it, *light, sample, 2, nSamples1, nSamples2, nSamples3, scene, sampler, arena)
							* invSamples2;
		}
		for (Point2f sample : samples3) {
			L += this->EstimateDirect(it, *light, sample, 3, nSamples1, nSamples2, nSamples3, scene, sampler, arena)
							* invSamples3;
		}

		// if averaging, need to average total samples to get correct result (assumes equal sample split between strategies)
		if (nSamples3 > 0 && !_mis) {
			L /= 3.;
		} else if (nSamples2 > 0 && !_mis) {
			L /= 2.;
		}
	}

	return L;
}

Spectrum PointLineDirectIntegrator::UniformSampleOneLight(
		const Interaction &it,
			const Scene &scene, MemoryArena &arena, Sampler &sampler,
			const std::vector<std::pair<int, int>> &nLightSamples,
			bool handleMedia, const Distribution1D *lightDistrib) const {
	
	ProfilePhase p(Prof::DirectLighting);
	Spectrum L(0.f);

	// Randomly choose a single light to sample, _light_
	int nLights = int(scene.lights.size());
	if (nLights == 0)
		return Spectrum(0.f);

	int lightNum;
	Float lightPdf;
	if (lightDistrib) {
		lightNum = lightDistrib->SampleDiscrete(sampler.Get1D(), &lightPdf);
		if (lightPdf == 0)
			return Spectrum(0.f);
	} else {
		lightNum = std::min((int) (sampler.Get1D() * nLights), nLights - 1);
		lightPdf = Float(1) / nLights;
	}

	const std::shared_ptr<Light> &light = scene.lights[lightNum];

	if (!light->CanIlluminate(it))
		return L;

	// Only sample that one light
	for (size_t j = lightNum; j < lightNum+1; ++j) {
		// Accumulate contribution of j'th light to L
		const std::shared_ptr<Light> &light = scene.lights[j];
		const std::pair<int, int>& lightSample = _nLightSamples[j];
		const int nSamples = lightSample.first * lightSample.second;
		double invSamples = 1./nSamples;

        if (!light->CanIlluminate(it))
            continue;

		// Use Shirley remapping to split random samples between MIS strategies
		const Point2f *allSamples = sampler.Get2DArray(lightSample.first, lightSample.second);

		std::vector<Point2f> samples1, samples2, samples3;

		if (nSamples == 1) {
			samples1.push_back(allSamples[0]);

			if (_nStrategies > 1) {
				const Point2f *moreSamples = sampler.Get2DArray(lightSample.first, lightSample.second);
				samples2.push_back(moreSamples[0]);
			}

			if (_nStrategies > 2) {
				const Point2f *evenMoreSamples = sampler.Get2DArray(lightSample.first, lightSample.second);
				samples3.push_back(evenMoreSamples[0]);
			}

		} else {
			if (_nStrategies == 1) {
				for (int i = 0; i < nSamples; i++) {
					samples1.push_back(allSamples[i]);
				}

			} else if (_nStrategies == 2) {
				int halfSamples = nSamples/2;
				double sqrtSamples = std::sqrt(nSamples);

				if ((int)sqrtSamples % 2 != 0)
					Error("For remapping equally between 2 MIS strategies, must use some sample count N x N where N is divisible by 2.");
				
				if (_stratified) {
					// split samples along first dimension
					for (int i = 0; i < nSamples; i++) {
						Point2f s = allSamples[i];
						if (s.x < 0.5) {
							s.x = s.x * 2;
							samples1.push_back(s);
						} else {
							s.x = s.x * 2 - 1;
							samples2.push_back(s);
						}
					}
				} else {
					int nSamples1 = floor(_ratio * nSamples);
    				int nSamples2 = nSamples - nSamples1;

					for (int i = 0; i < nSamples1; i++) {
						samples1.push_back(allSamples[i]);
					}
					for (int i = nSamples1; i < nSamples1+nSamples2; i++)
						samples2.push_back(allSamples[i]);
				}

			} else {
				int thirdSamples = nSamples / 3;
				double sqrtSamples = std::sqrt(nSamples);

				if ((int)sqrtSamples % 3 != 0)
					Error("For remapping equally between 3 MIS strategies, must use some sample count N x N where N is divisible by 3.");
				
				if (_stratified) {
					// split samples along first dimension
					double oneThird = 1./3.;
					for (int i = 0; i < nSamples; i++) {
						Point2f s = allSamples[i];
						if (s.x < oneThird) {
							s.x = s.x * 3;
							samples1.push_back(s);
						} else if (s.x < 1-oneThird) {
							s.x = s.x * 3 - 1;
							samples2.push_back(s);
						} else {
							s.x = s.x * 3 - 2;
							samples3.push_back(s);
						}
					}
				} else {
					int nSamples1 = std::ceil(_ratio * nSamples);
    				int nSamples2 = std::ceil(_ratio * nSamples);
    				int nSamples3 = nSamples - nSamples1 - nSamples2;

					for (int i = 0; i < nSamples1; i++)
						samples1.push_back(allSamples[i]);
					for (int i = nSamples1; i < nSamples1+nSamples2; i++)
						samples2.push_back(allSamples[i]);
					for (int i = nSamples1+nSamples2; i < nSamples1+nSamples2+nSamples3; i++)
						samples3.push_back(allSamples[i]);
				}
			}
		}

		int nSamples1 = samples1.size(), nSamples2 = samples2.size(), nSamples3 = samples3.size();
		Float invSamples1 = 1./nSamples1, invSamples2 = 1./nSamples2, invSamples3 = 1./nSamples3;

		for (Point2f sample : samples1) {
			L += this->EstimateDirect(it, *light, sample, 1, nSamples1, nSamples2, nSamples3, scene, sampler, arena)
							* invSamples1;
		}
		for (Point2f sample : samples2) {
			L += this->EstimateDirect(it, *light, sample, 2, nSamples1, nSamples2, nSamples3, scene, sampler, arena)
							* invSamples2;
		}
		for (Point2f sample : samples3) {
			L += this->EstimateDirect(it, *light, sample, 3, nSamples1, nSamples2, nSamples3, scene, sampler, arena)
							* invSamples3;
		}
	}

	return L;
}

Spectrum PointLineDirectIntegrator::EstimateDirect(const Interaction &it, 
		const Light &light, Point2f &u, int strategy, int nSamples1, int nSamples2, int nSamples3, 
		const Scene &scene, Sampler &sampler, MemoryArena &arena, bool handleMedia, bool specular) const {
	
	Spectrum Ld(0.f);

	/**************************************************************************
	 * Only one strategy
	 **************************************************************************/

	if (_nStrategies == 1) {
		Vector3f wi;
		Float pdf, smoothstep;

		Ld += this->EstimateIntermediate(it, light, u, scene, sampler, arena, strategy, &wi, &pdf, &smoothstep);

	/**************************************************************************
	 * MIS/Avg between 2 strategies
	 **************************************************************************/

	} else if (_nStrategies == 2) {
		// MIS between two strategies
		Vector3f wi;
		Float pdf_f, smoothstep_f = 1.f;

		// Sample one strategy
		Spectrum Ld1 = this->EstimateIntermediate(it, light, u, scene, sampler, arena, strategy, &wi, &pdf_f, &smoothstep_f);

		int samples_f = (strategy == 1) ? nSamples1 : nSamples2;
		int samples_g = (strategy == 1) ? nSamples2 : nSamples1;

		if (samples_g == 0 || !_mis) {
			Ld += Ld1;
		} else if (!Ld1.IsBlack()) {
			// Get pdf of other strategy
			int strategy_g = (strategy == 1) ? 2 : 1;
			Float smoothstep_g = 1.f;
			float pdf_g = this->PdfIntermediate(it, light, scene, sampler, arena, strategy_g, wi, &smoothstep_g);

			// Get MIS weight
			Float weight;
			if (_smoothstepWidth != 0) {
				// !!! IMPORTANT: when smoothing, assumes strategy 1 is a strategy that needs smoothing and strategy 2 is a line strategy
				weight = SmoothBalanceHeuristicF(samples_f, pdf_f, samples_g, pdf_g, smoothstep_g);
			} else {
				weight = BalanceHeuristic(samples_f, pdf_f, samples_g, pdf_g);
			}

			Ld += Ld1 * weight;
		}

	/**************************************************************************
	 * MIS/Avg between 3 strategies
	 **************************************************************************/

	} else {
		// MIS between three strategies
		Vector3f wi;
		Float pdf_f, smoothstep_f = 1.f;

		// Sample strategy 1
		Spectrum Ld1 = this->EstimateIntermediate(it, light, u, scene, sampler, arena, strategy, &wi, &pdf_f, &smoothstep_f);
		
		int samples_f, samples_g, samples_h;
		if (strategy == 1) {
			samples_f = nSamples1;
			samples_g = nSamples2;
			samples_h = nSamples3;
		} else if (strategy == 2) {
			samples_f = nSamples2;
			samples_g = nSamples1;
			samples_h = nSamples3;
		} else {
			samples_f = nSamples3;
			samples_g = nSamples1;
			samples_h = nSamples2;
		}

		if ((samples_g == 0 && samples_h == 0) || !_mis) {
			Ld += Ld1;
		} else if (!Ld1.IsBlack()) {
			int strategy_g, strategy_h;
			if (strategy == 1) {
				strategy_g = 2;
				strategy_h = 3;
			} else if (strategy == 2) {
				strategy_g = 1;
				strategy_h = 3;
			} else {
				strategy_g = 1;
				strategy_h = 2;
			}

			// Get pdf of other strategies
			Float smoothstep_g = 1.f, smoothstep_h = 1.f;
			float pdf_g = this->PdfIntermediate(it, light, scene, sampler, arena, strategy_g, wi, &smoothstep_g);
			float pdf_h = this->PdfIntermediate(it, light, scene, sampler, arena, strategy_h, wi, &smoothstep_h);

			// Get MIS weight
			Float weight;
			if (_smoothstepWidth != 0) {
				// !!! IMPORTANT: when smoothing, assumes strategy 1 is a strategy that needs smoothing and strategies 2 and 3 are line strategies
				if (strategy == 1) {
					weight = SmoothBalanceHeuristicF(samples_f, pdf_f, samples_g, pdf_g, samples_h, pdf_h, smoothstep_g, smoothstep_h);
				} else {
					weight = SmoothBalanceHeuristicG(samples_f, pdf_f, samples_g, pdf_g, samples_h, pdf_h, smoothstep_f, smoothstep_h);
				}
			} else {
				weight = BalanceHeuristic(samples_f, pdf_f, samples_g, pdf_g, samples_h, pdf_h);
			}

			Ld += Ld1 * weight;
		}
	}

	return Ld;
}

Spectrum PointLineDirectIntegrator::EstimateIntermediate(const Interaction &it, 
		const Light &light, const Point2f &uLight, const Scene &scene,
		Sampler &sampler, MemoryArena &arena, const int strategy, 
		Vector3f *wi, Float *pdf, Float *smoothstep, bool handleMedia, bool specular) const {

	Spectrum Ld(0.f);

	SampleStrategy strategyType;
	if (strategy == 1) {
		strategyType = _strategy1;
	} else if (strategy == 2) {
		strategyType = _strategy2;
	} else {
		strategyType = _strategy3;
	}

	/**************************************************************************
	 * PT_SURF_AREA: sample uniformly over surface area of light source
	 **************************************************************************/

	if (strategyType == SampleStrategy::PT_SURF_AREA) {
		VisibilityTester visibility;
		Spectrum Li = light.Sample_Li(it, uLight, wi, pdf, &visibility);

		VLOG(2) << "EstimateDirect uLight:" << uLight << " -> Li: " << Li
							<< ", wi: " << wi << ", pdf: " << *pdf;

		if (!Li.IsBlack() && *pdf > 0) {
			// Evaluate BSDF for light sampling strategy
			const SurfaceInteraction &isect = (const SurfaceInteraction &) it;

			if (Dot(*wi, isect.n) < 0)
				return Ld;

			Spectrum f = isect.bsdf->f(isect.wo, *wi, BSDF_NOSPECULAR)
								* AbsDot(*wi, isect.shading.n);

			if (!f.IsBlack()) {
				// Compute effect of visibility for light source sample
				// StartTimer("SamplerIntegrator/VisibilityTime");
				if (!visibility.Unoccluded(scene)) {
					VLOG(2) << "  shadow ray blocked";
					Li = Spectrum(0.f);
				} else {
					VLOG(2) << "  shadow ray unoccluded";
				}
				// StopTimer("SamplerIntegrator/VisibilityTime");
				nSceneTraversals++;

				// Add light's contribution to reflected radiance
				if (!Li.IsBlack()) {
					Ld += f * Li / *pdf;
				}
			}
		}

	/**************************************************************************
	 * PT_BRDF: importance sample BRDF (picking wi over shading hemisphere)
	 **************************************************************************/

	} else if (strategyType == SampleStrategy::PT_BRDF) {
		BxDFType bsdfFlags = specular ? BSDF_ALL : BxDFType(BSDF_ALL & ~BSDF_SPECULAR);

		// Sample BSDF
		if (!IsDeltaLight(light.flags)) {
			BxDFType sampledType;
			bool sampledSpecular = (sampledType & BSDF_SPECULAR) != 0;
			const SurfaceInteraction &isect = (const SurfaceInteraction &)it;

			Spectrum f = isect.bsdf->Sample_f(isect.wo, wi, uLight, pdf, bsdfFlags, &sampledType)
									* AbsDot(*wi, isect.shading.n);

			VLOG(2) << "  BSDF / phase sampling f: " << f << ", scatteringPdf: " << *pdf;
				
			if (!f.IsBlack() && *pdf > 0) {
				// Find intersection and compute transmittance
				SurfaceInteraction lightIsect;
				Ray ray = it.SpawnRay(*wi);
				// StartTimer("SamplerIntegrator/VisibilityTime");
				bool foundSurfaceInteraction = scene.Intersect(ray, &lightIsect);
				// StopTimer("SamplerIntegrator/VisibilityTime");
				// look for intersection of direction and light source
				nSceneTraversals++;

				// Add light contribution from material sampling
				Spectrum Li(0.f);
				if (foundSurfaceInteraction) {
					if (lightIsect.primitive->GetAreaLight() == &light)
						Li = lightIsect.Le(-*wi);
				} else  {
					Li = light.Le(ray);
				}

				if (!Li.IsBlack()) {
					Ld += f * Li / *pdf;
				} 
			}
		}

	/********************************************************************************************
	 * PT_SOLID_ANGLE: point sample uniformly over solid angle subtended by light over shading point
	 ********************************************************************************************/

	} else if (strategyType == SampleStrategy::PT_SOLID_ANGLE) {

		VisibilityTester visibility;
		const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
		Spectrum Li = light.Sample_Li_SA(it, uLight, wi, pdf, &visibility);

		VLOG(2) << "EstimateDirect uLight:" << uLight << " -> Li: " << Li
							<< ", wi: " << wi << ", pdf: " << *pdf;

		if (!Li.IsBlack() && *pdf > 0) {
			// Evaluate BSDF for light sampling strategy
			const SurfaceInteraction &isect = (const SurfaceInteraction &) it;

			if (Dot(*wi, isect.n) < 0)
				return Ld;

			Spectrum f = isect.bsdf->f(isect.wo, *wi, BSDF_NOSPECULAR)
								* AbsDot(*wi, isect.shading.n);

			if (!f.IsBlack()) {
				// Compute effect of visibility for light source sample
				// StartTimer("SamplerIntegrator/VisibilityTime");
				if (!visibility.Unoccluded(scene)) {
					VLOG(2) << "  shadow ray blocked";
					Li = Spectrum(0.f);
				} else {
					VLOG(2) << "  shadow ray unoccluded";
				}
				// StopTimer("SamplerIntegrator/VisibilityTime");
				nSceneTraversals++;

				// Add light's contribution to reflected radiance
				if (!Li.IsBlack()) {
					Ld += f * Li / *pdf;
				}
			}
		}

	/**************************************************************************
	 * LINE: line sampling
	 **************************************************************************/

	} else {
		Ld += this->EstimateDirectLine(it, light, uLight, strategy, scene, sampler, arena, wi, pdf, smoothstep);
	}

	return Ld;
}

Spectrum PointLineDirectIntegrator::EstimateDirectLine(const Interaction &it,
		const Light &light, const Point2f& uLight, const int strategy, 
		const Scene &scene, Sampler &sampler, MemoryArena &arena,
		Vector3f *wi, Float *pdf, Float *smoothstep, bool handleMedia, bool specular) const {

	Spectrum Ld(0.f);

	SampleStrategy strategyType;
	Float lineAngle;
	bool lineIS;

	if (strategy == 1) {
		strategyType = _strategy1;
		lineAngle = _lineAngle1;
		lineIS = _lineIS1;
	} else if (strategy == 2) {
		strategyType = _strategy2;
		lineAngle = _lineAngle2;
		lineIS = _lineIS2;
	} else {
		strategyType = _strategy3;
		lineAngle = _lineAngle3;
		lineIS = _lineIS3;
	}

	/**************************************************************************
	 * Sample a line on the light source
	 **************************************************************************/

	// retrieve the interaction
	const SurfaceInteraction &isect = (const SurfaceInteraction &) it;

	Spectrum Li;
	linesampling::LineInteraction lineInteraction;
	if (strategyType == SampleStrategy::LINE_SOLID_ANGLE) {
		Li = light.Sample_SphQuadLine_Li(isect.p, uLight.x, lineAngle, lineInteraction);
	} else {
		Li = light.Sample_Line_Li(it, Point2f(lineAngle, uLight.x), lineInteraction, lineIS);
	}

	// exit if no contribution, or when the line sample is too small
	if (Li.IsBlack() || lineInteraction.length < EPSILON)
		return Ld;

	/**************************************************************************
	 * Retrieve the intersection information
	 **************************************************************************/

	// check whether the normals have to be inverted
	const Point3f &p = isect.p;
	const Normal3f &n = isect.n;
	const Normal3f &sn = isect.shading.n;
	const Vector3f &reflection = Reflect(isect.wo, sn);

	/**********************************************************************
	 * Clip the line against the plane spanned by the real surface normal
	 **********************************************************************/

	const Line3f& sampleLine = lineInteraction.line;
	const Normal3f &ln = lineInteraction.n;
	const Vector3f &wo = isect.wo;

	// make shadow triangle
	Point3f p1 = sampleLine.p1;
	Point3f p2 = sampleLine.p2;
	Vector3f e1 = p1 - p;
	Vector3f e2 = p2 - p;
	Vector3f d = p2 - p1;

	// check whether shadow triangle is degenerate
	Vector3f cross = Cross(e1, e2);
	if (cross.Length() < EPSILON)
		return Ld;

	// check whether the point is completely beneath the line sample
	if (Dot(p - p1, ln) <= 0)
		return Ld;

	// check if line is under the plane of the real surface normal
	if (Dot(e1, n) <= 0 && Dot(e2, n) <= 0)
		return Ld;

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
				return Ld;
		} else if (Dot(n, e2) > 0) {
			// line.p2 is on the correct side, line.p1 is under the plane
			p1 = p1 + d * (s + EPSILON);
			e1 = p1 - p;

			// if line.p2 is on the correct side and s is larger than
			// 1 - epsilon, nothing of the line segment remains
			if (s > 1.0 - EPSILON)
				return Ld;
		} else {
			return Ld;
		}
		// reassign since p1 or p2 may have changed
		d = p2 - p1;
	}

	/********************************************************************************
	 * Intersect shadow triangle with the scene, subtract shadowed regions of line
	 ********************************************************************************/

	IntervalCollection collection(0, 1);

    // add a small epsilon to the line sample to avoid self intersections
    const Vector3f& isectLineOffset = Vector3f(ln * EPSILON);
    const Point3f& isectPoint = p + Vector3f(EPSILON * n);
    const Point3f& isectLineOrigin = p1 + isectLineOffset;
    const Point3f& isectLineDestination = p2 + isectLineOffset;

    // construct the shadow triangle query
    linesampling::LineVisibilityQuery3f query(isectPoint, isectLineOrigin,
            isectLineDestination);

    // perform the intersection test
    // check if line segment is completely obscured
	nSceneTraversals++;
	// StartTimer("SamplerIntegrator/VisibilityTime");
	bool intersectL = scene.IntersectL(query, n, it.time, collection);
	// StopTimer("SamplerIntegrator/VisibilityTime");
    if (intersectL && collection.empty())
		return Spectrum(0.f);

	if (_method == RenderMethod::LINESAMPLECOUNT)
		return Spectrum(collection.size());

	/**************************************************************************
	 * Sample point within line
	 **************************************************************************/

	bool foundEndpts = false;
	Point3f pLight, v1, v2;
	Float pointPdf, occDist1, occDist2;
	linesampling::Interval pointInterval;

	if (strategyType == SampleStrategy::LINE_SOLID_ANGLE) {
		pLight = light.Sample_SphQuadPointInLine(collection, lineInteraction,
									isect.p, uLight.y, &pointPdf, lineAngle);

	} else {
		Float visPdf;
		Float t = collection.SampleInterval(pointInterval, &visPdf, &occDist1, &occDist2, uLight.y);

		// Get endpoints of visible region in world space
		foundEndpts = true;
		v1 = p1 + pointInterval.min * d;
		v2 = p1 + pointInterval.max * d;

		pLight = v1 + t * (v2 - v1);
		pointPdf = (1 / d.Length()) * visPdf;
	}

	// visible intervals are degenerate
	if (pointPdf < EPSILON)
		return Ld;

	*wi = Normalize(pLight-p);

	// Smoothstep
	if (_smoothstepWidth != 0) {
		if (!foundEndpts) {
			Float u = (pLight-p1).Length() / d.Length();
			collection.GetPointInterval(pointInterval, &occDist1, &occDist2, u);

			// Get endpoints of visible region in world space
			v1 = p1 + pointInterval.min * d;
			v2 = p1 + pointInterval.max * d;
		}

		// transform occDist to world space
		occDist1 *= d.Length();
		occDist2 *= d.Length();

		Float smoothstepLight = light.smoothstep(pLight, _smoothstepWidth);
		// Float smoothstepLine = this->smoothstepPoints(pLight, v1, v2); // old method
		Float smoothstepLine = this->smoothstepPoints(pLight, v1, v2, occDist1, occDist2); // new method (decreasing ss width)
		*smoothstep = smoothstepLight * smoothstepLine; // multiply
		// *smoothstep = std::min(smoothstepLight, smoothstepLine); // min
	}

	/**************************************************************************
	 * Shading
	 **************************************************************************/

	Spectrum f = isect.bsdf->f(isect.wo, *wi, BSDF_NOSPECULAR)
					* AbsDot(*wi, isect.shading.n);
 
	// No need to check for visibility, etc. as we would have exited early after line-scene intersection if occluded
	Interaction iLight;
	iLight.p = pLight;
	iLight.n = lineInteraction.n;
	Li = light.L(iLight, -*wi); // handles textures

	if (strategyType == SampleStrategy::LINE_SOLID_ANGLE) {
		*pdf = lineInteraction.pdf * pointPdf;
	} else {
		// transform pdf to solid angle space
		double saPdf = DistanceSquared(isect.p, pLight) / AbsDot(lineInteraction.n, -*wi);
		*pdf = lineInteraction.pdf * pointPdf * saPdf;	
	}

	return f * Li / *pdf;
}

Float PointLineDirectIntegrator::PdfIntermediate(const Interaction &it, const Light &light, const Scene &scene,
		Sampler &sampler, MemoryArena &arena, const int strategy, Vector3f &wi, Float *smoothstep,
		bool handleMedia, bool specular) const {

	SampleStrategy strategyType;
	if (strategy == 1) {
		strategyType = _strategy1;
	} else if (strategy == 2) {
		strategyType = _strategy2;
	} else {
		strategyType = _strategy3;
	}

	if (strategyType == SampleStrategy::PT_SURF_AREA) {
		return light.Pdf_Li(it, wi);
		
	} else if (strategyType == SampleStrategy::PT_BRDF) {
		BxDFType bsdfFlags = BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
		const SurfaceInteraction &isect = (const SurfaceInteraction &)it;

		return isect.bsdf->Pdf(isect.wo, wi, bsdfFlags);

	} else if (strategyType == SampleStrategy::PT_SOLID_ANGLE) {
		const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
		return light.Pdf_SA(isect, wi);

	} else {
		return this->PdfLine(it, light, scene, sampler, arena, strategy, wi, smoothstep);
	}
}

Float PointLineDirectIntegrator::PdfLine(const Interaction &it, const Light &light, const Scene &scene,
		Sampler &sampler, MemoryArena &arena, const int strategy, Vector3f& wi, Float *smoothstep,
		bool handleMedia, bool specular) const {

	Spectrum Ld(0.f);

	SampleStrategy strategyType;
	Float lineAngle;
	bool lineIS;

	if (strategy == 1) {
		strategyType = _strategy1;
		lineAngle = _lineAngle1;
		lineIS = _lineIS1;
	} else if (strategy == 2) {
		strategyType = _strategy2;
		lineAngle = _lineAngle2;
		lineIS = _lineIS2;
	} else {
		strategyType = _strategy3;
		lineAngle = _lineAngle3;
		lineIS = _lineIS3;
	}

	/*********************************************************************************
	 * Retrieve line and point on light source given this direction wi and line angle
	 *********************************************************************************/

	// retrieve the interaction
	const SurfaceInteraction &isect = (const SurfaceInteraction &) it;

	Point3f pLight;
	LineInteraction lineInteraction;
	Float linePdf;

	if (strategyType == SampleStrategy::LINE_SOLID_ANGLE) {
		linePdf = light.PdfSQLine(it, wi, pLight, lineInteraction, lineAngle);
	} else {
		linePdf = light.PdfLine(it, wi, pLight, lineInteraction, lineAngle, lineIS);
	}

	/**************************************************************************
	 * Retrieve the intersection information
	 **************************************************************************/

	// transform pdf to solid angle space
	double saPdf = DistanceSquared(isect.p, pLight) / AbsDot(lineInteraction.n, -wi);

	// check whether the normals have to be inverted
	const Point3f &p = isect.p;
	const Normal3f &n = isect.n;
	const Normal3f &sn = isect.shading.n;
	const Vector3f &reflection = Reflect(isect.wo, sn);

	/**********************************************************************
	 * Clip the line against the plane spanned by the real surface normal
	 **********************************************************************/

	const Line3f& sampleLine = lineInteraction.line;
	const Normal3f &ln = lineInteraction.n;
	const Vector3f &wo = isect.wo;

	// make shadow triangle
	Point3f p1 = sampleLine.p1;
	Point3f p2 = sampleLine.p2;
	Vector3f e1 = p1 - p;
	Vector3f e2 = p2 - p;
	Vector3f d = p2 - p1;

	// check whether the point is completely beneath the line sample
	if (Dot(p - p1, ln) <= 0)
		return 0.f;

	// check if line is under the plane of the real surface normal
	if (Dot(e1, n) <= 0 && Dot(e2, n) <= 0)
		return 0.f;

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
				return 0.f;
		} else if (Dot(n, e2) > 0) {
			// line.p2 is on the correct side, line.p1 is under the plane
			p1 = p1 + d * (s + EPSILON);
			e1 = p1 - p;

			// if line.p2 is on the correct side and s is larger than
			// 1 - epsilon, nothing of the line segment remains
			if (s > 1.0 - EPSILON)
				return 0.f;
		} else {
			return 0.f;
		}
		// reassign since p1 or p2 may have changed
		d = p2 - p1;
	}

	/********************************************************************************
	 * Intersect shadow triangle with the scene, subtract shadowed regions of line
	 ********************************************************************************/

	IntervalCollection collection(0, 1);

    // add a small epsilon to the line sample to avoid self intersections
    const Vector3f& isectLineOffset = Vector3f(ln * EPSILON);
    const Point3f& isectPoint = p + Vector3f(EPSILON * n);
    const Point3f& isectLineOrigin = p1 + isectLineOffset;
    const Point3f& isectLineDestination = p2 + isectLineOffset;

    // construct the shadow triangle query
    linesampling::LineVisibilityQuery3f query(isectPoint, isectLineOrigin,
            isectLineDestination);

    // perform the intersection test
	// StartTimer("SamplerIntegrator/VisibilityTime");
    scene.IntersectL(query, n, it.time, collection);
	// StopTimer("SamplerIntegrator/VisibilityTime");
	nSceneTraversals++;

	/**************************************************************************
	 * Get pdf of point along line
	 **************************************************************************/

	Float pointPdf, occDist1, occDist2;
	linesampling::Interval pointInterval;
	
	if (strategyType == SampleStrategy::LINE_SOLID_ANGLE) {
		// Check if point is shadowed
		Float u = (pLight-p1).Length() / d.Length();
		if (!collection.GetPointInterval(pointInterval, &occDist1, &occDist2, u))
			return 0.f;

		pointPdf = light.PdfSQPoint(collection, lineInteraction, isect.p, lineAngle);

	} else {
		Float u = (pLight-p1).Length() / d.Length();
		Float visPdf = collection.VisPdf(pointInterval, &occDist1, &occDist2, u);
		if (visPdf == 0)
			return 0.;
		
		pointPdf = visPdf * (1 / d.Length());
	}

	// Smoothstep
	if (_smoothstepWidth != 0) {
		// Get endpoints of visible region in world space
		Point3f v1 = p1 + pointInterval.min * d;
		Point3f v2 = p1 + pointInterval.max * d;

		// Transform occDist to world space
		occDist1 *= d.Length();
		occDist2 *= d.Length();

		Float smoothstepLight = light.smoothstep(pLight, _smoothstepWidth);
		// Float smoothstepLine = this->smoothstepPoints(pLight, v1, v2); // old method
		Float smoothstepLine = this->smoothstepPoints(pLight, v1, v2, occDist1, occDist2); // new method (decreasing ss width)
		*smoothstep = smoothstepLight * smoothstepLine; // multiply
		// *smoothstep = std::min(smoothstepLight, smoothstepLine); // min
	}

	return (strategyType == SampleStrategy::LINE_SOLID_ANGLE) ? linePdf * pointPdf : linePdf * pointPdf * saPdf;
}

void AssignStrategy(PointLineDirectIntegrator::SampleStrategy& strategy, const std::string& strategyString) {
	if (strategyString == "pt-surf-area") {
		strategy = PointLineDirectIntegrator::SampleStrategy::PT_SURF_AREA;
	} else if (strategyString == "pt-solid-angle") {
		strategy = PointLineDirectIntegrator::SampleStrategy::PT_SOLID_ANGLE;
	} else if (strategyString == "brdf") {
		strategy = PointLineDirectIntegrator::SampleStrategy::PT_BRDF;
	} else if (strategyString == "line-surf-area") {
		strategy = PointLineDirectIntegrator::SampleStrategy::LINE_SURF_AREA;
	} else if (strategyString == "line-solid-angle") {
		strategy = PointLineDirectIntegrator::SampleStrategy::LINE_SOLID_ANGLE;
	} else {
		if (strategyString != "none") {
			Error("Unknown render method string for "
					"PointLineDirectIntegrator specified!\n"
					"Choices are: pt-surf-area | pt-solid-angle | brdf | "
					"line-surf-area | line-solid-angle");
		}
		strategy = PointLineDirectIntegrator::SampleStrategy::NONE;
	}
}

PointLineDirectIntegrator *CreatePointLineDirectIntegrator(
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

	/***************************************************************************
	 * Find the render method
	 **************************************************************************/

	const std::string& methodString = params.FindOneString("method", "render");

	PointLineDirectIntegrator::RenderMethod method;
	if (methodString == "lines" || methodString == "linecount") {
		method = PointLineDirectIntegrator::RenderMethod::LINESAMPLECOUNT;
	} else if (methodString == "drawlines") {
        method = PointLineDirectIntegrator::RenderMethod::DRAWLINES;
	} else if (methodString == "intviz") {
		method = PointLineDirectIntegrator::RenderMethod::INTVIZ;
    } else {
		if (methodString != "render") {
			Error("Unknown render method string '%s' for "
					"DirectLineLightingIntegrator specified!\nUsing "
					"'render' instead!", methodString.c_str());
		}
		method = PointLineDirectIntegrator::RenderMethod::RENDER;
	}

	/***************************************************************************
	 * Find sampling strategies
	 **************************************************************************/

	const std::string& strategyString1 = params.FindOneString("strategy1", "none");
	const std::string& strategyString2 = params.FindOneString("strategy2", "none");
	const std::string& strategyString3 = params.FindOneString("strategy3", "none");

	PointLineDirectIntegrator::SampleStrategy strategy1;
	PointLineDirectIntegrator::SampleStrategy strategy2;
	PointLineDirectIntegrator::SampleStrategy strategy3;

	AssignStrategy(strategy1, strategyString1);
	AssignStrategy(strategy2, strategyString2);
	AssignStrategy(strategy3, strategyString3);

	if (strategy1 == PointLineDirectIntegrator::SampleStrategy::NONE)
		Error("Sample strategy 1 cannot be none!");
	
	int nStrategies = 1;
	if (strategy2 != PointLineDirectIntegrator::SampleStrategy::NONE) 
		nStrategies++;
	if (strategy3 != PointLineDirectIntegrator::SampleStrategy::NONE) 
		nStrategies++;

	/***************************************************************************
	 * Find other sampler params
	 **************************************************************************/

	// importance sample lines by length
	bool lineIS1 = params.FindOneBool("importance1", true);
	bool lineIS2 = params.FindOneBool("importance2", true);
	bool lineIS3 = params.FindOneBool("importance3", true);

	// line directions
	Float lineAngle1 = params.FindOneFloat("angle1", 0);
	Float lineAngle2 = params.FindOneFloat("angle2", 0);
	Float lineAngle3 = params.FindOneFloat("angle3", 0);

	// distance to smoothstep from edges of light and occluders (in world space)
	Float smoothstepWidth = params.FindOneFloat("smoothstep", 0.0);

	// max seconds
	float maxSeconds = params.FindOneFloat("max-seconds", -1.0);

	// max scene traversals
	long maxSceneTraversals = params.FindOneInt("max-traversals", -1);

	// for MIS, whether to split and remap one sample array or create separate sample arrays
	bool stratified = params.FindOneBool("stratified", false);

	// whether to MIS or average for multiple strategies (MIS if true, average if false)
	bool mis = params.FindOneBool("mis", true);

	// ratio of point samples to line samples for MIS (assume points is strategy 1)
	Float ratio = params.FindOneFloat("ratio", 0.5);

	return new PointLineDirectIntegrator(maxDepth, camera, sampler, pixelBounds, method, 
			strategy1, strategy2, strategy3, nStrategies, stratified, mis, lineAngle1, lineAngle2, lineAngle3,
			lineIS1, lineIS2, lineIS3, smoothstepWidth, ratio, maxSeconds, maxSceneTraversals);
}

}
