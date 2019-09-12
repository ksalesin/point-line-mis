/*
 * interval.h
 *
 *  Created on: Jan 13, 2016
 *      Author: niels
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_LINESAMPLING_LINESAMPLING_INTERVAL_H
#define PBRT_LINESAMPLING_LINESAMPLING_INTERVAL_H

// linesampling/interval.h*
#include "pbrt.h"
#include <list>

using pbrt::Float;
using pbrt::Infinity;
using pbrt::Error;
using pbrt::EPSILON;

namespace linesampling {

class Interval {
public:
	Interval() :
			min(Infinity), max(-Infinity) {
	}

	Interval(Float value) :
			min(value), max(value) {
	}

	Interval(Float min, Float max) :
			min(min), max(max) {
		if (min > max) {
			Error("Interval::error in constructor! %.16f %.16f\n", min, max);
		}
	}

	inline bool contains(const Float& value) const {
		return value >= min && value <= max;
	}

	inline Float length() const {
		return max - min;
	}

	inline bool overlaps(const Interval& interval) const {
		if (min > interval.max)
			return false;
		if (max < interval.min)
			return false;
		return true;
	}

	inline bool operator<(const Interval& interval) const {
		return min < interval.min;
	}

	inline bool operator<(const Float value) const {
		return min < value;
	}

	inline bool empty() const {
		return max - min <= EPSILON;
	}

	inline void Print() {
		fprintf(stdout, "[%.10f, %.10f]\n", min, max);
	}

	Float min;
	Float max;
};

typedef std::list<Interval>::iterator IntervalIterator;
typedef std::list<Interval>::reverse_iterator IntervalReverseIterator;
typedef std::list<Interval>::const_iterator IntervalConstIterator;
typedef std::list<Interval>::const_reverse_iterator IntervalConstReverseIterator;
class IntervalCollection {
public:
	IntervalCollection(Float min, Float max);

	Float minimum() const;
	Float maximum() const;
	void split(Float value);
	void clipToMax(Float max);
	void clipToMin(Float min);
	void clip(Float min, Float max);
	bool empty() const;
	bool contains(Float value) const;
	bool overlaps(const Interval &interval) const;
	bool overlaps(const Interval &interval, bool& complete) const;

	bool completeOverlap(const Interval &interval) const;

	uint32_t size() const;

	IntervalIterator begin();
	IntervalConstIterator begin() const;
	IntervalReverseIterator rbegin();
	IntervalConstReverseIterator rbegin() const;

	IntervalIterator end();
	IntervalConstIterator end() const;
	IntervalReverseIterator rend();
	IntervalConstReverseIterator rend() const;

	bool subtract(const Interval& interval);

	bool GetPointInterval(Interval& interval, Float* occDist1, Float* occDist2, const float u);
	Float SampleInterval(Interval& interval, Float* visPdf, const float u);
	Float SampleInterval(Interval& interval, Float* visPdf, Float* occDist1, Float* occDist2, const float u);
	
	Float VisPdf(Interval& interval, float u);
	Float VisPdf(Interval& interval, Float* occDist1, Float* occDist2, float u);

private:
//	Float length;
	std::list<Interval> intervals;
};
}

#endif /* LINESAMPLING_LINESAMPLING_INTERVAL_H */
