/*
 * quad.h
 *
 *  Created on: Jul 28, 2017
 *      Author: niels
 */

#ifndef SHAPES_QUAD_H_
#define SHAPES_QUAD_H_

#include "shape.h"
#include "stats.h"

namespace pbrt {

using linesampling::IntervalConstIterator;

// Quad Public Declaration
class Quad: public Shape {
public:
	// Quad Public Methods
	Quad(const Transform *ObjectToWorld, const Transform *WorldToObject,
		bool reverseOrientation, const Point3f& o, const Vector3f& u,
		const Vector3f& v, const Vector3f& un, const Vector3f& vn, const Vector3f& wn,
		const double eu, const double ev) :
		Shape(ObjectToWorld, WorldToObject, reverseOrientation), 
				o(o), u(u), v(v), un(un), vn(vn), wn(wn), eu(eu), ev(ev) { }

	inline Point3f GetPoint(int index) const {
		switch (index) {
		case 0:
			return o;
		case 1:
			return o + u;
		case 2:
			return o + u + v;
		case 3:
			return o + v;
		default:
			return o;
		}
	}
	inline Point3f operator[](int index) const {
		return GetPoint(index);
	}

	Bounds3f ObjectBound() const;
	Bounds3f WorldBound() const;
	bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
			bool testAlphaTexture = true) const;
	bool IntersectP(const Ray &ray, bool testAlphaTexture = true) const;
	Float Area() const;

	using Shape::Sample;  // Bring in the other Sample() overload.
	Interaction Sample(const Point2f &u, Float *pdf) const;

	// Returns the solid angle subtended by the triangle w.r.t. the given
	// reference point p.
	Float SolidAngle(const Point3f &p, int nSamples = 0) const override;

	bool CanEvaluateAnalytically(const SurfaceInteraction& isect,
			Spectrum * spectrum = nullptr) const override;

	bool CanIlluminate(const Interaction& isect) const override;

	// Sample a line on quad by surface area method
	LineInteraction LineSample(const Point2f &u,
			bool importanceSample = true) const override;
	Float PdfLine(const Interaction &ref, const Vector3f &wi, Point3f &pLight,
			LineInteraction &lineInteraction, const Float& lineAngle, bool importanceSample) const override;

	// Sample a line on quad by solid angle method
	LineInteraction SampleSQLine(const Point3f& p, const Float u, int strategy) const override;
	Float PdfSQLine(const Interaction &ref, const Vector3f& wi, Point3f& pLight, LineInteraction &interaction, int strategy) const override;

	// Sample a point within given line by solid angle method
	Point3f SampleSQPoint(const IntervalCollection intervals, const LineInteraction &interaction, 
			const Point3f& p, const Float v, Float *pointPdf, int strategy) const override;
	Float PdfSQPoint(const IntervalCollection intervals, const LineInteraction &interaction, 
			const Point3f& p, int strategy) const override;

	// Sample a point uniformly on quad by solid angle method (Ureña et al. implementation)
	Interaction SampleSolidAngle(const Point3f& p, const Point2f &sample, Float* pdf) const override;
	Float PdfSolidAngle(const Point3f& p) const override;

	// Get UV coordinates of given point on quad
	virtual Point2f GetUV(const Point3f& p) const override;

	// Sample point on quad given UV coordinates
	virtual Interaction SampleUV(const Point2f& uv) const override; 

	// Visualizes lines on quad in given direction (same as how line lighting would be evaluated)
	bool OnColoredLine(const Interaction &intr, float direction) const override;

	virtual Float smoothstep(const Point3f& pLight, const Float smoothstepWidth) const override;

private:
	// Triangle Private Methods
	void GetUVs(Point2f uv[4]) const {
		uv[0] = Point2f(0, 0);
		uv[1] = Point2f(1, 0);
		uv[2] = Point2f(0, 1);
		uv[3] = Point2f(1, 1);
	}

	// Quad Private Data
	const Point3f o;
	const Vector3f u;
	const Vector3f v;

	const Vector3f un;
	const Vector3f vn;
	const Vector3f wn;

	const double eu;
	const double ev;
};

std::shared_ptr<Shape> CreateQuadShape(const Transform *o2w,
		const Transform *w2o, bool reverseOrientation, const ParamSet &params);

}

#endif /* SHAPES_QUAD_H_ */
