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

// samplers/stratified.cpp*
#include "samplers/stratified.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {

// StratifiedSampler Method Definitions
void StratifiedSampler::StartPixel(const Point2i &p) {
	ProfilePhase _(Prof::StartPixel);

	// Generete single array for pixel
	for (size_t i = 0; i < pixel2DArraySizes.size(); ++i) {
		const std::pair<int, int> &arraySizes = pixel2DArraySizes[i];
		const int n = arraySizes.first * arraySizes.second;

		// Change to multijitter for research
		if (multijitter) {
			MultijitterSample2D(&pixelArray2D[i][0], arraySizes.first,
				arraySizes.second, rng);
		} else {
			StratifiedSample2D(&pixelArray2D[i][0], arraySizes.first,
					arraySizes.second, rng, true); // plain jitter (not uniform)
			Shuffle(&pixelArray2D[i][0], n, 1, rng); // latin hypercube
		}
	}

	// Generate single stratified samples for the pixel
	for (size_t i = 0; i < samples1D.size(); ++i) {
		StratifiedSample1D(&samples1D[i][0], xPixelSamples * yPixelSamples, rng,
				true);
		Shuffle(&samples1D[i][0], xPixelSamples * yPixelSamples, 1, rng);
	}
	for (size_t i = 0; i < samples2D.size(); ++i) {
		StratifiedSample2D(&samples2D[i][0], xPixelSamples, yPixelSamples, rng,
				true);
		Shuffle(&samples2D[i][0], xPixelSamples * yPixelSamples, 1, rng);
	}

	// Generate arrays of stratified samples for the pixel
	for (size_t i = 0; i < samples1DArraySizes.size(); ++i)
		for (int64_t j = 0; j < samplesPerPixel; ++j) {
			int count = samples1DArraySizes[i];
			StratifiedSample1D(&sampleArray1D[i][j * count], count, rng,
					true);
			Shuffle(&sampleArray1D[i][j * count], count, 1, rng);
		}

	for (size_t i = 0; i < samples2DArraySizes.size(); ++i)
		for (int64_t j = 0; j < samplesPerPixel; ++j) {
			const std::pair<int, int> &arraySizes = samples2DArraySizes[i];
			const int n = arraySizes.first * arraySizes.second;

			if (multijitter) {
				MultijitterSample2D(&sampleArray2D[i][j * n], arraySizes.first,
					arraySizes.second, rng);
			} else {
				StratifiedSample2D(&sampleArray2D[i][j * n], arraySizes.first,
						arraySizes.second, rng, true); // plain jitter (not uniform)
				Shuffle(&sampleArray2D[i][j * n], n, 1, rng); // latin hypercube
			}
		}
	PixelSampler::StartPixel(p);
}

std::unique_ptr<Sampler> StratifiedSampler::Clone(int seed) {
	StratifiedSampler *ss = new StratifiedSampler(*this);
	ss->rng.SetSequence(seed);
	return std::unique_ptr<Sampler>(ss);
}

StratifiedSampler *CreateStratifiedSampler(const ParamSet &params) {
	bool fixed = params.FindOneBool("fixedpixels", true);
	bool multijitter = params.FindOneBool("multijitter", true);
	int xsamp = params.FindOneInt("xsamples", 4);
	int ysamp = params.FindOneInt("ysamples", 4);
	int sd = params.FindOneInt("dimensions", 1);
	if (PbrtOptions.quickRender)
		xsamp = ysamp = 1;
	return new StratifiedSampler(fixed, xsamp, ysamp, multijitter, sd);
}

}  // namespace pbrt
