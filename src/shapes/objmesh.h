/*
 * objmesh.h
 *
 *  Created on: Jul 20, 2017
 *      Author: niels
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef SHAPES_OBJMESH_H_
#define SHAPES_OBJMESH_H_

// shapes/objmesh.h*
#include "pbrt.h"
#include "shapes/triangle.h"

namespace pbrt {
std::vector<std::shared_ptr<Shape>> CreateOBJMesh(const Transform *o2w,
		const Transform *w2o, bool reverseOrientation, const ParamSet &params,
		std::map<std::string, std::shared_ptr<Texture<Float>>>*floatTextures =
		nullptr);

	}

#endif /* SHAPES_OBJMESH_H_ */
