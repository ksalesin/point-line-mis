/*
 * box.h
 *
 *  Created on: Aug 26, 2017
 *      Author: niels
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef SHAPES_BOX_H_
#define SHAPES_BOX_H_

// shapes/box.h*
#include "shape.h"
#include "stats.h"

namespace pbrt {

// Box declaration
class Box: public Shape {
public:
	// Box Public Methods
	Box(const Transform *ObjectToWorld, const Transform *WorldToObject,
			bool reverseOrientation, const Point3f& p0, const Point3f& p1) :
			Shape(ObjectToWorld, WorldToObject, reverseOrientation), bounds(p0,
					p1) {

	}

	Bounds3f ObjectBound() const {
		return bounds;
	}

	bool Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
			bool testAlphaTexture) const {
		ProfilePhase p(Prof::ShapeIntersect);

		// transform the ray
		Ray ray = (*WorldToObject)(r);

		// parametric indices
		Float t0 = 0, t1 = ray.tMax;
		int hitDimension; // dimension  of the plane which is hit
		bool hitSwapped; // whether a swap occured after the hit

		for (int i = 0; i < 3; ++i) {
			// Update interval for _i_th bounding box slab
			Float invRayDir = 1 / ray.d[i];
			Float tNear = (bounds.pMin[i] - ray.o[i]) * invRayDir;
			Float tFar = (bounds.pMax[i] - ray.o[i]) * invRayDir;

			// Update parametric interval from slab intersection $t$ values
			bool swapped = false;
			if (tNear > tFar) {
				std::swap(tNear, tFar);
				swapped = true;
			}

			// Update _tFar_ to ensure robust ray--bounds intersection
			tFar *= 1 + 2 * gamma(3);

			if (tNear > t0) {
				isect->n = Normal3f(0, 0, 0);
				isect->n[i] = swapped ? 1 : -1;
				hitDimension = i;
				hitSwapped = swapped;
				t0 = tNear;
			}

			t1 = tFar < t1 ? tFar : t1;
			if (t0 > t1)
				return false;
		}
		if (tHit)
			*tHit = t0;
		int nextIndex = (hitDimension + 1) % 3;
		int previousIndex = (hitDimension + 2) % 3;

		Point3f pHit = ray(t0);
		Vector3f diag = bounds.Diagonal();
		Float u = (pHit[nextIndex] - bounds.pMin[nextIndex])
				/ (diag[nextIndex]);
		Float v = (pHit[previousIndex] - bounds.pMin[previousIndex])
				/ (diag[previousIndex]);
		Vector3f pError = gamma(5) * Abs((Vector3f) pHit);
		Vector3f dpdu, dpdv;
		Normal3f dndu, dndv;

		dpdu[previousIndex] = (hitSwapped ? -1.0 : 1.0) / diag[previousIndex];
		dpdv[nextIndex] = 1.0 / diag[nextIndex];

		*isect = (*ObjectToWorld)(
				SurfaceInteraction(pHit, pError, Point2f(u, v), -ray.d, dpdu,
						dpdv, dndu, dndv, ray.time, this));
		return true;
	}

	bool IntersectP(const Ray &r, bool testAlphaTexture) const {
		//ProfilePhase p(Prof::ShapeIntersectP);
		const Ray& ray = (*WorldToObject)(r);
//		Float t0, t1;
		return bounds.IntersectP(ray);
//		if (!bounds.IntersectP(ray, &t0, &t1))
//			return false;
//		return (t0 < ray.tMax) && (t1 > 0);
	}

	Float Area() const override {
		return bounds.SurfaceArea();
	}

	Interaction Sample(const Point2f &u, Float *pdf) const override {
		Error("Unimplemented Box::Sample() called!");
		return Interaction();
	}

	Float Pdf(const Interaction &ref, const Vector3f &wi) const override {
		Error("Unimplemented Box::Pdf called!");
		return 0.;
	}
private:
	const Bounds3f bounds;
};

std::shared_ptr<Shape> CreateBoxShape(const Transform *o2w,
		const Transform *w2o, bool reverseOrientation, const ParamSet &params);

} //namespace pbrt

#endif /* SHAPES_BOX_H_ */
