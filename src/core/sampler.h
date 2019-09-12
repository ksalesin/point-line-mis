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

#ifndef PBRT_CORE_SAMPLER_H
#define PBRT_CORE_SAMPLER_H

// core/sampler.h*
#include "pbrt.h"
#include "geometry.h"
#include "rng.h"
#include <inttypes.h>

namespace pbrt {

// Sampler Declarations
class Sampler {
public:
	// Sampler Interface
	virtual ~Sampler();
	Sampler(bool fixedPixelSamples, int64_t samplesPerPixel);
	virtual void StartPixel(const Point2i &p);
	virtual Float Get1D() = 0;
	virtual Point2f Get2D() = 0;
	virtual uint32_t UniformUInt32() {
		return rng.UniformUInt32();
	}
	virtual Float UniformFloat() {
		return rng.UniformFloat();
	}
	CameraSample GetCameraSample(const Point2i &pRaster);

	void Request1DArray(int n);
	void Request2DArray(const std::pair<int, int>& n);
	void Request2DArray(int x, int y) {
		Request2DArray(std::make_pair(x, y));
	}
	void Request2DPixelArray(const std::pair<int, int>& n);
	void Request2DPixelArray(int x, int y) {
		Request2DPixelArray(std::make_pair(x, y));
	}
	virtual int RoundCount(int n) const {
		return n;
	}
	const Float *Get1DArray(int n);
	const Point2f *Get2DArray(const std::pair<int, int> pair) {
		return Get2DArray(pair.first, pair.second);
	}
	const Point2f *Get2DArray(int xSamples, int ySamples);
	const Point2f *Get2DPixelArray(const std::pair<int, int> pair) {
		return Get2DPixelArray(pair.first, pair.second);
	}
	const Point2f *Get2DPixelArray(int xSamples, int ySamples);
	virtual bool StartNextSample();
	virtual std::unique_ptr<Sampler> Clone(int seed) = 0;
	virtual bool SetSampleNumber(int64_t sampleNum);
	std::string StateString() const {
		return StringPrintf("(%d,%d), sample %" PRId64, currentPixel.x,
				currentPixel.y, currentPixelSampleIndex);
	}
	int64_t CurrentSampleNumber() const {
		return currentPixelSampleIndex;
	}

	// Sampler Public Data
	const int64_t samplesPerPixel;
	const bool fixedPixelSamples;
protected:
	// Sampler Protected Data
	Point2i currentPixel;
	int64_t currentPixelSampleIndex;

	std::vector<int> samples1DArraySizes;
	std::vector<std::pair<int, int>> samples2DArraySizes;
	std::vector<std::pair<int, int>> pixel2DArraySizes;
	std::vector<std::vector<Float>> sampleArray1D;
	std::vector<std::vector<Point2f>> sampleArray2D;
	std::vector<std::vector<Point2f>> pixelArray2D;
	RNG rng;

private:
	// Sampler Private Data
	size_t array1DOffset, array2DOffset, pixelArray2DOffset;
};

class PixelSampler: public Sampler {
public:
	// PixelSampler Public Methods
	PixelSampler(bool fixedPixelSamples, int64_t samplesPerPixel,
			int nSampledDimensions);
	bool StartNextSample();
	bool SetSampleNumber(int64_t);
	Float Get1D();
	Point2f Get2D();

protected:
	// PixelSampler Protected Data
	std::vector<std::vector<Float>> samples1D;
	std::vector<std::vector<Point2f>> samples2D;
	int current1DDimension = 0, current2DDimension = 0;
};

class GlobalSampler: public Sampler {
public:
	// GlobalSampler Public Methods
	bool StartNextSample();
	void StartPixel(const Point2i &);
	bool SetSampleNumber(int64_t sampleNum);
	Float Get1D();
	Point2f Get2D();
	GlobalSampler(bool fixedPixelSamples, int64_t samplesPerPixel) :
			Sampler(fixedPixelSamples, samplesPerPixel) {
	}
	virtual int64_t GetIndexForSample(int64_t sampleNum) const = 0;
	virtual Float SampleDimension(int64_t index, int dimension) const = 0;

private:
	// GlobalSampler Private Data
	int dimension;
	int64_t intervalSampleIndex;
	static const int arrayStartDim = 5;
	int arrayEndDim;
};

}  // namespace pbrt

#endif  // PBRT_CORE_SAMPLER_H
