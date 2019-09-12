/*
 * interval.cpp
 *
 *  Created on: Jan 27, 2016
 *      Author: niels
 */

#include "point-line-mis/interval.h"
#include <algorithm>
#include "geometry.h"

namespace linesampling {

class IntervalComparator {
public:
	bool operator()(const Interval &interval, Float value) {
		return interval.min < value;
	}

	bool operator()(Float value, const Interval &interval) {
		return value < interval.min;
	}
};

IntervalCollection::IntervalCollection(Float min, Float max) {
	intervals.push_back(Interval(min, max));
}

Float IntervalCollection::minimum() const {
	if (empty())
		return Infinity;
	return intervals.front().min;
}
Float IntervalCollection::maximum() const {
	if (empty())
		return -Infinity;
	return intervals.back().max;
}

bool IntervalCollection::empty() const {
	return intervals.empty();
}

void IntervalCollection::clipToMax(Float max) {
	if (maximum() > max) {
		IntervalIterator it = end();
		while (it != begin()) {
			--it;

			if (it->min >= max) {
				it = intervals.erase(it);
			} else if (it->max >= max) {
				it->max = max;
				break;
			}

		}
	}

	if (maximum() > max)
		fprintf(stderr, "IntervalCollection::Clip:: max is incorrect!\n");
}

void IntervalCollection::clipToMin(Float min) {
	if (minimum() < min) {
		IntervalIterator it = begin();
		while (it != end()) {
			if (it->max <= min)
				it = intervals.erase(it);
			else if (it->min <= min) {
				it->min = min;
				break;
			} else
				++it;
		}
	}

	if (minimum() < min)
		fprintf(stderr,
				"IntervalCollection::Clip:: minimum is incorrect! %f vs %f\n",
				minimum(), min);
}
void IntervalCollection::clip(Float min, Float max) {
	clipToMin(min);
	clipToMax(max);
}

void IntervalCollection::split(Float value) {
	if (value >= maximum() || value <= minimum())
		return;

	IntervalIterator it = begin();

	while (it != end()) {
		if (it->min <= value && value <= it->max) {
			const Float max = it->max;
			it->max = value;

			++it;
			intervals.insert(it, Interval(value, max));

			return;
		}
		++it;
	}
}

bool IntervalCollection::contains(Float value) const {
	IntervalConstIterator it = begin();

//	IntervalConstIterator it = std::upper_bound(begin(), end(), value,
//			IntervalComparator());

	while (it != end()) {
		if (it->min > value)
			return false;
		if (it->contains(value))
			return true;

		++it;
	}
	return false;
}

uint32_t IntervalCollection::size() const {
	return intervals.size();
}

IntervalIterator IntervalCollection::begin() {
	return intervals.begin();
}

IntervalConstIterator IntervalCollection::begin() const {
	return intervals.begin();
}

IntervalReverseIterator IntervalCollection::rbegin() {
	return intervals.rbegin();
}

IntervalConstReverseIterator IntervalCollection::rbegin() const {
	return intervals.rbegin();
}

IntervalIterator IntervalCollection::end() {
	return intervals.end();
}

IntervalConstIterator IntervalCollection::end() const {
	return intervals.end();
}

IntervalReverseIterator IntervalCollection::rend() {
	return intervals.rend();
}

IntervalConstReverseIterator IntervalCollection::rend() const {
	return intervals.rend();
}

bool IntervalCollection::overlaps(const Interval &interval) const {
	if (empty() || interval.empty())
		return false;

	{
		const Float min = intervals.front().min;
		const Float max = intervals.back().max;

		if (interval.max < min || interval.min > max)
			return false;
	}

	IntervalConstIterator start = std::upper_bound(begin(), end(),
			interval.min - (interval.max - interval.min) - EPSILON,
			IntervalComparator());
	if (start != begin())
		--start;

	while (start != end()) {
		// exit condition
		if (start->min > interval.max)
			break;
		else if (start->max < interval.min)
			++start;
		else
			return true;
	}
	return false;
}

bool IntervalCollection::overlaps(const Interval &interval,
		bool& complete) const {
	if (empty() || interval.empty())
		return false;
	complete = false;

	{
		const Float min = intervals.front().min;
		const Float max = intervals.back().max;

		if (interval.max < min || interval.min > max)
			return false;
	}

	IntervalConstIterator start = std::upper_bound(begin(), end(),
			interval.min - (interval.max - interval.min) - EPSILON,
			IntervalComparator());
	if (start != begin())
		--start;

	while (start != end()) {
		// exit condition
		if (start->min > interval.max)
			break;
		else if (start->max < interval.min)
			++start;
		else {
			complete = start->min <= interval.min && start->max >= interval.max;
			return true;
		}
	}
	return false;
}

bool IntervalCollection::completeOverlap(const Interval &interval) const {
	if (empty() || interval.empty())
		return false;

	IntervalConstIterator start = std::upper_bound(begin(), end(),
			interval.min - (interval.max - interval.min) - EPSILON,
			IntervalComparator());
	if (start != begin())
		--start;

	while (start != end()) {
		// exit condition
		if (start->min > interval.max)
			break;
		else if (start->max < interval.min)
			++start;
		else if (start->min <= interval.min && start->max >= interval.max)
			return true;
		return false;
	}
	return false;
}

bool IntervalCollection::subtract(const Interval& interval) {
	if (interval.empty())
		return false;

	{
		const Float min = intervals.front().min;
		const Float max = intervals.back().max;

		if (interval.max < min || interval.min > max)
			return false;
	}

	const static IntervalComparator comp;
	const Float minimum = 2.0 * interval.min - interval.max;

	IntervalIterator start = std::lower_bound(begin(), end(), minimum, comp);
	if (start != begin())
		--start;

	bool result = false;

	while (start != end()) {
		// test overlap
		if (start->min > interval.max)
			break;
		if (start->max < interval.min) {
			++start;
			continue;
		}

		if (interval.min <= start->min && interval.max >= start->max) {
			// the interval is larger, we just need to remove i
			start = intervals.erase(start);
		} else if (interval.min <= start->min) {
			// the interval starts on the left
			if (std::abs(start->max - interval.max) <= EPSILON) {
				start = intervals.erase(start);
			} else {
				*start = Interval(interval.max, start->max);
				++start;
			}
		} else if (interval.max >= start->max) {
			// the interval starts on the right
			if (std::abs(interval.min - start->min) <= EPSILON) {
				start = intervals.erase(start);
			} else {
				*start = Interval(start->min, interval.min);
				++start;
			}
		} else {
			Interval i0(start->min, interval.min);
			Interval i1(interval.max, start->max);

			Float i0l = std::abs(interval.min - start->min);
			Float i1l = std::abs(start->max - interval.max);

			if (i0l <= EPSILON && i1l <= EPSILON) {
				start = intervals.erase(start);
			} else if (i0l <= EPSILON) {
				*start = i1;
				++start;
			} else if (i1l <= EPSILON) {
				*start = i0;
				++start;
			} else {
				*start = i0;
				++start;
				start = intervals.insert(start, i1);
				++start;
			}
		}

		result = true;
	}

	return result;
}

bool IntervalCollection::GetPointInterval(Interval& interval, Float* occDist1, Float* occDist2, const float u) {
	if (empty())
		return false;

	IntervalConstIterator it = begin();
	while (it != end()) {
		if (it->min <= u && it->max >= u) {
			interval = (*it);
			
			// get occlusion distance of nearest occluder(s)
			auto next = std::next(it,1);
			auto prev = std::prev(it,1);
			*occDist1 = (it == begin()) ? it->min : std::abs(it->min - prev->max);
			*occDist2 = (next == end()) ? 1 - it->max : std::abs(next->min - it->max);

			return true;
		}
		it++;
	}

	return false;
}

Float IntervalCollection::SampleInterval(Interval& interval, Float* visPdf, const float u) {

	if (empty()) {
		// sample uniformly along occluded line (would never reach because early exits if empty)
		*visPdf = 0;
		interval = Interval(0,1);
		return u;
	}

	Float visLen = 0;

	IntervalConstIterator it = begin();
	while (it != end()) {
		visLen += it->length();
		it++;
	}

	*visPdf = 1 / visLen;
	Float sum = 0;

	it = begin();
	while (it != end()) {
		Float regionLen = it->length()/visLen;
		Float newsum = sum + regionLen;

		if (u <= newsum) {
			// Stretch u to [0,1] within visible region
			Float v = (u - sum) / regionLen;
			interval = (*it);
			return v;
		}

		sum = newsum;
		it++;
	}

	return 0.;
}

Float IntervalCollection::SampleInterval(Interval& interval, Float* visPdf, Float* occDist1, Float* occDist2, const float u) {

	if (empty()) {
		// sample uniformly along occluded line (would never reach because early exits if empty)
		*visPdf = 0; *occDist1 = 0; *occDist2 = 0;
		interval = Interval(0,1);
		return u;
	}

	Float visLen = 0;

	IntervalConstIterator it = begin();
	while (it != end()) {
		visLen += it->length();
		it++;
	}

	*visPdf = 1 / visLen;
	Float sum = 0;

	it = begin();
	while (it != end()) {
		Float regionLen = it->length()/visLen;
		Float newsum = sum + regionLen;

		if (u <= newsum) {
			// Stretch u to [0,1] within visible region
			Float v = (u - sum) / regionLen;
			interval = (*it);

			// get occlusion distance of nearest occluder(s)
			auto next = std::next(it,1);
			auto prev = std::prev(it,1);
			*occDist1 = (it == begin()) ? it->min : std::abs(it->min - prev->max);
			*occDist2 = (next == end()) ? 1 - it->max : std::abs(next->min - it->max);

			return v;
		}

		sum = newsum;
		it++;
	}

	return 0.;
}

Float IntervalCollection::VisPdf(Interval& interval, float u) {
	
	if (empty())
		return 0.;

	bool isShadowed = true;
	float visLen = 0;
	IntervalConstIterator it = begin();
	while (it != end()) {
		if (it->min <= u && it->max >= u) {
			isShadowed = false;
			interval = (*it);
		}

		visLen += it->length();
		++it;
	}

	// If point is shadowed, would never be created by line sampling strategy
	if (isShadowed) {
		interval = Interval(0,1);
		return 0.f;
	}	

	return 1 / visLen;
}

Float IntervalCollection::VisPdf(Interval& interval, Float* occDist1, Float* occDist2, float u) {
	
	if (empty())
		return 0.;

	bool isShadowed = true;
	float visLen = 0;
	IntervalConstIterator it = begin();
	while (it != end()) {
		if (it->min <= u && it->max >= u) {
			isShadowed = false;
			interval = (*it);

			// get occlusion distance of nearest occluder
			auto next = std::next(it,1);
			auto prev = std::prev(it,1);	
			*occDist1 = (it == begin()) ? it->min : std::abs(it->min - prev->max);
			*occDist2 = (next == end()) ? 1 - it->max : std::abs(next->min - it->max);
		}

		visLen += it->length();
		it++;
	}

	// If point is shadowed, would never be created by line sampling strategy
	if (isShadowed) {
		interval = Interval(0,1);
		return 0.f;
	}	

	return 1 / visLen;
}

} // namespace linesampling
