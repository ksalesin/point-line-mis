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

#ifndef PBRT_CORE_STATS_H
#define PBRT_CORE_STATS_H

// core/stats.h*
#include "pbrt.h"
#include <map>
#include <chrono>
#include <string>
#include <functional>

namespace pbrt {

// Statistics Declarations
class StatsAccumulator;
class StatRegisterer {
public:
	// StatRegisterer Public Methods
	StatRegisterer(std::function<void(StatsAccumulator &)> func) {
		if (!funcs)
			funcs = new std::vector<std::function<void(StatsAccumulator &)>>;
		funcs->push_back(func);
	}
	static void CallCallbacks(StatsAccumulator &accum);

private:
	// StatRegisterer Private Data
	static std::vector<std::function<void(StatsAccumulator &)>> *funcs;
};

void AddStatistic(const std::string &name, const std::string& value);
void AddStatistic(const std::string &name, Float value);
void AddStatistic(const std::string &name, int64_t value);

void StartTimer(const std::string &name);
void StopTimer(const std::string &name);
void PrintStats(FILE *dest);
void PrintEEAStats(FILE *timeDest);
void PrintEEAStats(FILE *timeDest, FILE *travDest);
void ClearStats();
void ReportThreadStats();

class StatsAccumulator {
public:
	// StatsAccumulator Public Methods
	void ReportCounter(const std::string &name, int64_t val) {
		counters[name] += val;
	}

	void StartTimer(const std::string &name) {
		auto cpuIterator = cpuStartTimes.find(name);
		if (cpuIterator != cpuStartTimes.end())
			Warning("CPU timer '%s' is already running!", name.c_str());

		auto wallIterator = wallStartTimes.find(name);
		if (wallIterator != wallStartTimes.end())
			Warning("Wall clock timer '%s' is already running!", name.c_str());

		cpuStartTimes[name] = std::clock();
		wallStartTimes[name] = std::chrono::high_resolution_clock::now();
	}

	void StopTimer(const std::string &name) {
		auto cpuIterator = cpuStartTimes.find(name);
		if (cpuIterator == cpuStartTimes.end())
			Warning("CPU timer '%s' was not started!", name.c_str());
		auto wallIterator = wallStartTimes.find(name);
		if (wallIterator == wallStartTimes.end())
			Warning("Wall timer '%s' was not started!", name.c_str());

		// Get the current times
		auto nowCPU = std::clock();
		auto nowWall = std::chrono::high_resolution_clock::now();

		// Get the difference
		Float cpuTime = Float(nowCPU - cpuIterator->second)
				/ Float(CLOCKS_PER_SEC);
		std::chrono::duration<double> dur = nowWall - wallIterator->second;
		Float wallTime = dur.count();
		// Float wallTime = Float(
		// 		std::chrono::duration_cast<std::chrono::milliseconds>(
		// 				nowWall - wallIterator->second).count()) / 1000.0;

		cpuTimerValues[name] += cpuTime;
		wallTimerValues[name] += wallTime;

		cpuStartTimes.erase(cpuIterator);
		wallStartTimes.erase(wallIterator);
	}

	void AddStatistic(const std::string &name, const std::string& value) {
		if (stringStats.find(name) != stringStats.end())
			Error("Statistics already contains string statistic with name '%s'",
					name.c_str());
		stringStats[name] = value;
	}

	void AddStatistic(const std::string &name, Float value) {
		if (floatStats.find(name) != floatStats.end())
			Error("Statistics already contains float statistic with name '%s'",
					name.c_str());
		floatStats[name] = value;
	}

	void AddStatistic(const std::string &name, int64_t value) {
		if (intStats.find(name) != intStats.end())
			Error("Statistics already contains int statistic with name '%s'",
					name.c_str());
		intStats[name] = value;
	}

	void ReportMemoryCounter(const std::string &name, int64_t val) {
		memoryCounters[name] += val;
	}
	void ReportIntDistribution(const std::string &name, int64_t sum,
			int64_t count, int64_t min, int64_t max) {
		intDistributionSums[name] += sum;
		intDistributionCounts[name] += count;
		if (intDistributionMins.find(name) == intDistributionMins.end())
			intDistributionMins[name] = min;
		else
			intDistributionMins[name] = std::min(intDistributionMins[name],
					min);
		if (intDistributionMaxs.find(name) == intDistributionMaxs.end())
			intDistributionMaxs[name] = max;
		else
			intDistributionMaxs[name] = std::max(intDistributionMaxs[name],
					max);
	}
	void ReportFloatDistribution(const std::string &name, double sum,
			int64_t count, double min, double max) {
		floatDistributionSums[name] += sum;
		floatDistributionCounts[name] += count;
		if (floatDistributionMins.find(name) == floatDistributionMins.end())
			floatDistributionMins[name] = min;
		else
			floatDistributionMins[name] = std::min(floatDistributionMins[name],
					min);
		if (floatDistributionMaxs.find(name) == floatDistributionMaxs.end())
			floatDistributionMaxs[name] = max;
		else
			floatDistributionMaxs[name] = std::max(floatDistributionMaxs[name],
					max);
	}
	void ReportPercentage(const std::string &name, int64_t num, int64_t denom) {
		percentages[name].first += num;
		percentages[name].second += denom;
	}
	void ReportRatio(const std::string &name, int64_t num, int64_t denom) {
		ratios[name].first += num;
		ratios[name].second += denom;
	}

	void Print(FILE *file);
	void PrintEEA(FILE *timeFile);
	void PrintEEA(FILE *timeFile, FILE *travFile);
	void Clear();

private:
	// StatsAccumulator Private Data
	// std::map<std::string, std::chrono::system_clock::time_point> wallStartTimes;
	std::map<std::string, std::chrono::high_resolution_clock::time_point> wallStartTimes;
	std::map<std::string, clock_t> cpuStartTimes;

	std::map<std::string, Float> wallTimerValues;
	std::map<std::string, Float> cpuTimerValues;

	std::map<std::string, std::string> stringStats;
	std::map<std::string, int64_t> intStats;
	std::map<std::string, Float> floatStats;

	std::map<std::string, int64_t> counters;
	std::map<std::string, int64_t> memoryCounters;
	std::map<std::string, int64_t> intDistributionSums;
	std::map<std::string, int64_t> intDistributionCounts;
	std::map<std::string, int64_t> intDistributionMins;
	std::map<std::string, int64_t> intDistributionMaxs;
	std::map<std::string, double> floatDistributionSums;
	std::map<std::string, int64_t> floatDistributionCounts;
	std::map<std::string, double> floatDistributionMins;
	std::map<std::string, double> floatDistributionMaxs;
	std::map<std::string, std::pair<int64_t, int64_t>> percentages;
	std::map<std::string, std::pair<int64_t, int64_t>> ratios;
};

enum class Prof {
	SceneConstruction,
	AccelConstruction,
	TextureLoading,
	MIPMapCreation,

	IntegratorRender,
	SamplerIntegratorLi,
	SPPMCameraPass,
	SPPMGridConstruction,
	SPPMPhotonPass,
	SPPMStatsUpdate,
	BDPTGenerateSubpath,
	BDPTConnectSubpaths,
	LightDistribLookup,
	LightDistribSpinWait,
	LightDistribCreation,
	DirectLighting,
	BSDFEvaluation,
	BSDFSampling,
	BSDFPdf,
	BSSRDFEvaluation,
	BSSRDFSampling,
	PhaseFuncEvaluation,
	PhaseFuncSampling,
	AccelIntersect,
	AccelIntersectP,
	LightSample,
	LightPdf,
	MediumSample,
	MediumTr,
	TriIntersect,
	TriIntersectP,
	CurveIntersect,
	CurveIntersectP,
	ShapeIntersect,
	ShapeIntersectP,
	ComputeScatteringFuncs,
	GenerateCameraRay,
	MergeFilmTile,
	SplatFilm,
	AddFilmSample,
	StartPixel,
	GetSample,
	TexFiltTrilerp,
	TexFiltEWA,
	PhotonLookup,
	StochasticProxies,
	StochasticProxiesPreprocess,
	AccelIntersectL,
	TriIntersectL,
	NumProfCategories
};

static_assert((int)Prof::NumProfCategories <= 64,
		"No more than 64 profiling categories may be defined.");

inline uint64_t ProfToBits(Prof p) {
	return 1ull << (int) p;
}

static const char *ProfNames[] = { "Scene parsing and creation",
		"Acceleration structure creation", "Texture loading",
		"MIP map generation",

		"Integrator::Render()", "SamplerIntegrator::Li()", "SPPM camera pass",
		"SPPM grid construction", "SPPM photon pass",
		"SPPM photon statistics update", "BDPT subpath generation",
		"BDPT subpath connections", "SpatialLightDistribution lookup",
		"SpatialLightDistribution spin wait",
		"SpatialLightDistribution creation", "Direct lighting", "BSDF::f()",
		"BSDF::Sample_f()", "BSDF::PDF()", "BSSRDF::f()", "BSSRDF::Sample_f()",
		"PhaseFunction::p()", "PhaseFunction::Sample_p()",
		"Accelerator::Intersect()", "Accelerator::IntersectP()",
		"Light::Sample_*()", "Light::Pdf()", "Medium::Sample()", "Medium::Tr()",
		"Triangle::Intersect()", "Triangle::IntersectP()", "Curve::Intersect()",
		"Curve::IntersectP()", "Other Shape::Intersect()",
		"Other Shape::IntersectP()", "Material::ComputeScatteringFunctions()",
		"Camera::GenerateRay[Differential]()", "Film::MergeTile()",
		"Film::AddSplat()", "Film::AddSample()", "Sampler::StartPixelSample()",
		"Sampler::GetSample[12]D()", "MIPMap::Lookup() (trilinear)",
		"MIPMap::Lookup() (EWA)", "OcclusionMap::Lookup()",
		"Proxies::Direct Lighting", "Proxies::Preprocess",
		"Accelerator::IntersectL()", "Triangle::IntersectL()" };

static_assert((int)Prof::NumProfCategories ==
		sizeof(ProfNames) / sizeof(ProfNames[0]),
		"ProfNames[] array and Prof enumerant have different "
		"numbers of entries!");

extern PBRT_THREAD_LOCAL uint64_t ProfilerState;
inline uint64_t CurrentProfilerState() {
	return ProfilerState;
}

class ProfilePhase {
public:
// ProfilePhase Public Methods
	ProfilePhase(Prof p) {
		categoryBit = ProfToBits(p);
		reset = (ProfilerState & categoryBit) == 0;
		ProfilerState |= categoryBit;
	}
	~ProfilePhase() {
		if (reset)
			ProfilerState &= ~categoryBit;
	}
	ProfilePhase(const ProfilePhase &) = delete;
	ProfilePhase &operator=(const ProfilePhase &) = delete;

private:
// ProfilePhase Private Data
	bool reset;
	uint64_t categoryBit;
};

void InitProfiler();
void SuspendProfiler();
void ResumeProfiler();
void ProfilerWorkerThreadInit();
void ReportProfilerResults(FILE *dest);
void ClearProfiler();
void CleanupProfiler();

// Statistics Macros
#define STAT_COUNTER(title, var)                           \
    static PBRT_THREAD_LOCAL int64_t var;                  \
    static void STATS_FUNC##var(StatsAccumulator &accum) { \
        accum.ReportCounter(title, var);                   \
        var = 0;                                           \
    }                                                      \
    static StatRegisterer STATS_REG##var(STATS_FUNC##var)
#define STAT_MEMORY_COUNTER(title, var)                    \
    static PBRT_THREAD_LOCAL int64_t var;                  \
    static void STATS_FUNC##var(StatsAccumulator &accum) { \
        accum.ReportMemoryCounter(title, var);             \
        var = 0;                                           \
    }                                                      \
    static StatRegisterer STATS_REG##var(STATS_FUNC##var)

// Work around lack of support for constexpr in VS2013.
#ifdef PBRT_IS_MSVC2013
#define STATS_INT64_T_MIN LLONG_MAX
#define STATS_INT64_T_MAX _I64_MIN
#define STATS_DBL_T_MIN DBL_MAX
#define STATS_DBL_T_MAX -DBL_MAX
#else
#define STATS_INT64_T_MIN std::numeric_limits<int64_t>::max()
#define STATS_INT64_T_MAX std::numeric_limits<int64_t>::lowest()
#define STATS_DBL_T_MIN std::numeric_limits<double>::max()
#define STATS_DBL_T_MAX std::numeric_limits<double>::lowest()
#endif

#define STAT_INT_DISTRIBUTION(title, var)                                  \
    static PBRT_THREAD_LOCAL int64_t var##sum;                             \
    static PBRT_THREAD_LOCAL int64_t var##count;                           \
    static PBRT_THREAD_LOCAL int64_t var##min = (STATS_INT64_T_MIN);       \
    static PBRT_THREAD_LOCAL int64_t var##max = (STATS_INT64_T_MAX);       \
    static void STATS_FUNC##var(StatsAccumulator &accum) {                 \
        accum.ReportIntDistribution(title, var##sum, var##count, var##min, \
                                    var##max);                             \
        var##sum = 0;                                                      \
        var##count = 0;                                                    \
        var##min = std::numeric_limits<int64_t>::max();                    \
        var##max = std::numeric_limits<int64_t>::lowest();                 \
    }                                                                      \
    static StatRegisterer STATS_REG##var(STATS_FUNC##var)

#define STAT_FLOAT_DISTRIBUTION(title, var)                                  \
    static PBRT_THREAD_LOCAL double var##sum;                                \
    static PBRT_THREAD_LOCAL int64_t var##count;                             \
    static PBRT_THREAD_LOCAL double var##min = (STATS_DBL_T_MIN);            \
    static PBRT_THREAD_LOCAL double var##max = (STATS_DBL_T_MAX);            \
    static void STATS_FUNC##var(StatsAccumulator &accum) {                   \
        accum.ReportFloatDistribution(title, var##sum, var##count, var##min, \
                                      var##max);                             \
        var##sum = 0;                                                        \
        var##count = 0;                                                      \
        var##min = std::numeric_limits<double>::max();                       \
        var##max = std::numeric_limits<double>::lowest();                    \
    }                                                                        \
    static StatRegisterer STATS_REG##var(STATS_FUNC##var)

#define ReportValue(var, value)                                   \
    do {                                                          \
        var##sum += value;                                        \
        var##count += 1;                                          \
        var##min = std::min(var##min, decltype(var##min)(value)); \
        var##max = std::max(var##max, decltype(var##min)(value)); \
    } while (0)

#define STAT_PERCENT(title, numVar, denomVar)                 \
    static PBRT_THREAD_LOCAL int64_t numVar, denomVar;        \
    static void STATS_FUNC##numVar(StatsAccumulator &accum) { \
        accum.ReportPercentage(title, numVar, denomVar);      \
        numVar = denomVar = 0;                                \
    }                                                         \
    static StatRegisterer STATS_REG##numVar(STATS_FUNC##numVar)

#define STAT_RATIO(title, numVar, denomVar)                   \
    static PBRT_THREAD_LOCAL int64_t numVar, denomVar;        \
    static void STATS_FUNC##numVar(StatsAccumulator &accum) { \
        accum.ReportRatio(title, numVar, denomVar);           \
        numVar = denomVar = 0;                                \
    }                                                         \
    static StatRegisterer STATS_REG##numVar(STATS_FUNC##numVar)

}  // namespace pbrt

#endif  // PBRT_CORE_STATS_H
