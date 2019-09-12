/*
 * phongcache.h
 *
 *  Created on: Sep 14, 2017
 *      Author: niels
 */

#ifndef PHD_LINESAMPLING_PHONGCACHE_H_
#define PHD_LINESAMPLING_PHONGCACHE_H_

#include "pbrt.h"

using pbrt::Float;

namespace linesampling {

class PhongCache {
public:
	PhongCache(int n, int precision);
	~PhongCache() {
		delete[] points;
	}
	Float Get(Float u) const;

private:

	const int n;
	const int precision;
	const bool even;
	Float offset;
	Float * points;

};

}

#endif /* PHD_LINESAMPLING_PHONGCACHE_H_ */
