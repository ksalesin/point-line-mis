/*
 * box.cpp
 *
 *  Created on: Aug 26, 2017
 *      Author: niels
 */

#include "box.h"
#include "paramset.h"

namespace pbrt {

std::shared_ptr<Shape> CreateBoxShape(const Transform *o2w,
		const Transform *w2o, bool reverseOrientation, const ParamSet &params) {
	Point3f p0 = params.FindOnePoint3f("p0", Point3f(0, 0, 0));
	Point3f p1 = params.FindOnePoint3f("p1", Point3f(1, 1, 1));

	return std::make_shared<Box>(o2w, w2o, reverseOrientation, p0, p1);
}
}
