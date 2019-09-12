/*
 * phongcache.cpp
 *
 *  Created on: Sep 14, 2017
 *      Author: niels
 */

#include "phongcache.h"

using pbrt::Pi;
using pbrt::PiOver2;
using pbrt::TwoPi;
using pbrt::TwoOverPi;
using pbrt::IntegrateCosinus;

namespace linesampling {

PhongCache::PhongCache(int n, int precision) :
		n(n), precision(precision), points(new Float[precision + 1]), even(
				n % 2 == 0) {

	Float invPoints = 1.0 / Float(precision);

	for (int i = 0; i <= precision; ++i) {
		Float x = PiOver2 * i * invPoints;   // between 0 and Pi/2
		points[i] = IntegrateCosinus(x, n);
	}
	offset = 2.0 * IntegrateCosinus(PiOver2, n);
}

Float PhongCache::Get(Float u) const {
	if (u < 0)
		return -Get(-u);
	else if (u <= PiOver2) {
		const Float uOffset = u * TwoOverPi * precision;
		const int minimumIndex = static_cast<int>(std::floor(uOffset));
		const int maximumIndex = std::min(precision, minimumIndex + 1);
		const Float mu = uOffset - minimumIndex;
		return (1.0 - mu) * points[minimumIndex] + (mu) * points[maximumIndex];

	} else if (u <= Pi) {
		if (even)
			return offset - Get(Pi - u);
		else
			return Get(Pi - u); // mirror

	} else if (u <= TwoPi) {
		if (even)
			return offset + Get(u - Pi);
		else
			return -Get(u - Pi);
	} else
		return Get(fmod(u, TwoPi));
	return 0;
}

}
