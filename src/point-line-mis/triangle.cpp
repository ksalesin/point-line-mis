/*
 * triangle.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: niels
 */

#include "triangle.h"
#include "reflection.h"
#include "point-line-mis/line.h"

using pbrt::Reflect;

namespace niels {

Triangle::Triangle(const Point3f &p0, const Point3f &p1, const Point3f &p2) {
	vertices[0] = p0;
	vertices[1] = p1;
	vertices[2] = p2;
}

Float Triangle::Area() const {
	const Vector3f& e01 = vertices[1] - vertices[0];
	const Vector3f& e02 = vertices[2] - vertices[0];

	const Vector3f& cross = Cross(e01, e02);
	const Float length = cross.Length();

	return 0.5 * length;
}

void Triangle::Orient(const Triangle& triangle) {
	const Vector3f& e01 = vertices[1] - vertices[0];
	const Vector3f& e02 = vertices[2] - vertices[0];
	const Vector3f& e = Cross(e01, e02);

	const Vector3f& t01 = triangle.vertices[1] - triangle.vertices[0];
	const Vector3f& t02 = triangle.vertices[2] - triangle.vertices[0];
	const Vector3f& t = Cross(t01, t02);

	if (Dot(t, e) < 0) {
		std::swap(vertices[1], vertices[2]);
	}
}

Float Triangle::Diffuse(const Point3f& p, const Normal3f& n) const {
	Vector3f sum;

	Vector3f v[3];
	for (int i = 0; i < 3; ++i)
		v[i] = Normalize(vertices[i] - p);

	for (int i = 0; i < 3; ++i) {
		const Vector3f& vi = v[i];
		const Vector3f& vk = v[(i + 1) % 3];

		const Float dot = Dot(vi, vk);
		const Float angle = std::acos(dot);
		const Vector3f& cross = Normalize(Cross(vi, vk));

		sum += cross * angle;
	}

	return AbsDot(sum, n) * Inv2Pi;
}

Normal3f Triangle::Normal() const {
	const Vector3f& e01 = vertices[1] - vertices[0];
	const Vector3f& e02 = vertices[2] - vertices[0];

	return Normal3f(Normalize(Cross(e01, e02)));
}

bool Triangle::HasVertex(const Point3f& p, Float epsilon) const {
	for (const Point3f& point : vertices)
		if ((point - p).Length() <= epsilon)
			return true;
	return false;
}

const Point3f& Triangle::operator[](int index) const {
	return vertices[index];
}

Float Triangle::SolidAngle(const Point3f& p) const {
	return std::abs(SignedSolidAngle(p));
}

Float Triangle::SignedSolidAngle(const Point3f& p) const {
	// Project the vertices into the unit sphere around p.
	Vector3f pSphere[3] = { Normalize(vertices[0] - p), Normalize(
			vertices[1] - p), Normalize(vertices[2] - p) };

	// http://math.stackexchange.com/questions/9819/area-of-a-spherical-triangle
	// Girard's theorem: surface area of a spherical triangle on a unit
	// sphere is the 'excess angle' alpha+beta+gamma-pi, where
	// alpha/beta/gamma are the interior angles at the vertices.
	//
	// Given three vertices on the sphere, a, b, c, then we can compute,
	// for example, the angle c->a->b by
	//
	// cos theta =  Dot(Cross(c, a), Cross(b, a)) /
	//              (Length(Cross(c, a)) * Length(Cross(b, a))).
	//
	Vector3f cross01 = (Cross(pSphere[0], pSphere[1]));
	Vector3f cross12 = (Cross(pSphere[1], pSphere[2]));
	Vector3f cross20 = (Cross(pSphere[2], pSphere[0]));

	// Some of these vectors may be degenerate. In this case, we don't want
	// to normalize them so that we don't hit an assert. This is fine,
	// since the corresponding dot products below will be zero.
	if (cross01.LengthSquared() > 0)
		cross01 = Normalize(cross01);
	if (cross12.LengthSquared() > 0)
		cross12 = Normalize(cross12);
	if (cross20.LengthSquared() > 0)
		cross20 = Normalize(cross20);

	// We only need to do three cross products to evaluate the angles at
	// all three vertices, though, since we can take advantage of the fact
	// that Cross(a, b) = -Cross(b, a).
	return Pi - std::acos(Clamp(Dot(cross01, -cross12), -1, 1))
			- std::acos(Clamp(Dot(cross12, -cross20), -1, 1))
			- std::acos(Clamp(Dot(cross20, -cross01), -1, 1));
}

Float CosSumIntegral(Float x, Float c, int n) {
	int i = n % 2;

	const Float sinx = std::sin(x);
	const Float cosx = std::cos(x);

	Float F = (i == 0) ? x : sinx;
	Float S = 0;

	while (i <= n) {
		S += IntPow(c, i) * F;

		F = (IntPow(cosx, i + i) * sinx + (i + 1) * F) / Float(i + 2);
		i += 2;
	}
	return S;
}

/*
 * Implementation of CosSumIntegral (as defined by Arvo)
 *
 * Jim Arvo: Applications of Irradiance Tensors to the Simulation of
 * Non-Lambertian Phenomena, p. 5
 */
Float CosSumIntegral(Float x, Float y, Float c, int n) {
	// LINE 1: i <= if even(n) then 0 else 1;
	int i = n % 2;

	// initialize constants
	const Float cosx = std::cos(x);
	const Float cosy = std::cos(y);
	const Float sinx = std::sin(x);
	const Float siny = std::sin(y);
	const Float cosx2 = cosx * cosx;
	const Float cosy2 = cosy * cosy;
	const Float c2 = c * c;

	// LINE 2: F <= if even(n) then y-x else sin(y) - sin(x);
	Float F;

	// LINE 3: S <= 0;
	Float S = 0;

	// loop variables
	Float ci, cosx_loop, cosy_loop;

	if (i == 0) {
		// F <= if even(n) then y - x;
		F = y - x;
		ci = 1;
		cosx_loop = cosx;
		cosy_loop = cosy;
	} else {
		// F <= if odd(n) then sin(y) - sin(x);
		F = siny - sinx;
		ci = c;
		cosx_loop = cosx2;
		cosy_loop = cosy2;

	}

	// LINE 4: while (i<=n) do;
	while (i <= n) {
		// if i >= 0 then S <= S + c^i * F
		S += F * ci;

		// LINE 5: T <= cos^(i+1)*siny - cos^{i+1}*sinx
		const Float T = cosy_loop * siny - cosx_loop * sinx;

		// LINE 6: F = [T + (i + 1) * F] / (i + 2)
		F = (T + Float(i + 1) * F) / Float(i + 2);

		// LINE 7: i <= i + 2
		i += 2;

		// update loop variables
		ci *= c2;
		cosy_loop *= cosy2;
		cosx_loop *= cosx2;
	}

	return S;
}

Float LineIntegral(const Vector3f& A, const Vector3f& B, const Vector3f& w,
		int n) {
	// LINE 1: if (n < 0)
	if (n < 0)
		return 0;

	const Vector3f& nA = Normalize(A);
	const Vector3f& nB = Normalize(B);
	const Vector3f& nn = Normalize(Cross(nA, nB));

//	if (AbsDot(w, nA) < 1e-9 && AbsDot(w, nB) < 1e-9)
//		return 0;

	// Sanity checks
	CHECK_NEAR(Dot(nA, nn), 0, 1e-8);
	CHECK_NEAR(Dot(nB, nn), 0, 1e-8);

	// LINE 2: s <= Normalize[A];
	const Vector3f& s = nA;
	// LINE 3: t <= Normalize[(I -ss^T)B]
	const Vector3f& t = Normalize(Cross(nn, s));

	CHECK_NEAR(Dot(s, nA), 1, 1e-8);
	CHECK_NEAR(Dot(s, t), 0, 1e-8);
	CHECK_NEAR(Dot(t, nn), 0, 1e-8);

	// LINE 4: a <= w * s
	const Float a = Dot(w, s);
	// LINE 5: b <= w * t
	const Float b = Dot(w, t);
	// LINE 6: c <= sqrt(a*a + b*b)
	const Float c = std::sqrt(a * a + b * b);

	CHECK_NEAR(w.Length(), 1, 1e-8);

	const Float l = std::acos(Dot(nA, nB));
	if (l == 0)
		return 0;
	const Float sigma = b == 0 ? 0 : ((b < 0 ? -1 : 1) * std::acos(a / c));

	return CosSumIntegral(-sigma, l - sigma, c, n);
}

Float BoundaryIntegral(const Point3f& p, const Triangle& triangle,
		const Vector3f& w, int n) {
	Float b = 0;

	for (int i = 0; i < 3; ++i) {
		const Vector3f& A = triangle[i] - p;
		const Vector3f& B = triangle[(i + 1) % 3] - p;
		const Vector3f& cross = Cross(A, B);

		if (cross.LengthSquared() == 0)
			continue;

		b += Dot(Normalize(cross), w) * LineIntegral(A, B, w, n);
	}

	return b;
}

Float Triangle::Phong(const SurfaceInteraction& isect, int n) const {
	// Correctly orient the normal
	const Vector3f& r = Normalize(Reflect(isect.wo, isect.shading.n));

	Float sum = 0;
	std::vector<Triangle> clipped;
	Clip(isect.p, Normal3f(r), clipped);

	for (const Triangle& triangle : clipped) {
		Float boundary = -BoundaryIntegral(isect.p, triangle, r, n - 1);
		if (n % 2 == 0) {
			const Float sa = triangle.SolidAngle(isect.p);
			if (boundary > 0)
				sum -= sa;
			else
				sum += sa;
		}

		sum += boundary;

	}

	return std::abs(sum) / Float(n + 1);
}

void Triangle::Clip(const Point3f& p, const Normal3f& n,
		std::vector<Triangle> &clipped) const {

	// Normalize the normal
	const Normal3f& nn = Normalize(n);

	// Distances from vertices to plane
	const Float d[3] = { Dot(vertices[0] - p, nn), Dot(vertices[1] - p, nn),
			Dot(vertices[2] - p, nn) };

	// find the vertex with the largest index
	int i0;
	if (d[0] > d[1] && d[0] > d[2])
		i0 = 0;
	else if (d[1] > d[2])
		i0 = 1;
	else
		i0 = 2;

	// if the largest distance is smaller than zero, we can exit
	if (d[i0] <= 0)
		return;

	// retrieve the next indices
	int i1 = (i0 + 1) % 3;
	int i2 = (i0 + 2) % 3;

	// sort in order for i1 to be the second largest distance
	if (d[i1] < d[i2])
		std::swap(i1, i2);

	// check the smallest distance
	if (d[i2] >= 0) {
		clipped.push_back(*this);
		return;
	}

	if (d[i1] <= 0) {
		// i1 and i2 are on the other side
		Point3f p01, p02;
		if (!IntersectLinePlane(vertices[i0], vertices[i1], p, n, &p01)
				|| !IntersectLinePlane(vertices[i0], vertices[i2], p, n, &p02))
			return;

		Triangle c(vertices[i0], p02, p01);
		if (c.Area() < 1e-9)
			return;

		c.Orient(*this);
		CHECK_GT(Dot(c.Normal(), Normal()), 0);
		clipped.push_back(c);
	} else {
		/***********************************************************************
		 * i1 is on the same side, i2 is on the other side
		 *
		 * i0-----------i1
		 * \            /
		 *  \          /
		 *   \        /
		 *  p02------p12-----------------
		 *     \    /
		 *      \  /
		 *       i2
		 *
		 **********************************************************************/

		// i1 is on the same side, i2 is on the other side
		Point3f p02, p12;
		if (!IntersectLinePlane(vertices[i0], vertices[i2], p, n, &p02)
				|| !IntersectLinePlane(vertices[i1], vertices[i2], p, n, &p12))
			return;

		Triangle c1(vertices[i0], p12, p02);
		Triangle c2(vertices[i0], vertices[i1], p12);
		if (c1.Area() > 1e-9) {
			c1.Orient(*this);
			CHECK_GT(Dot(c1.Normal(), Normal()), 0);
			clipped.push_back(c1);
		}

		if (c2.Area() > 1e-9) {
			c2.Orient(*this);
			CHECK_GT(Dot(c2.Normal(), Normal()), 0);
			clipped.push_back(c2);
		}
	}
}

void Triangle::Print(FILE* file) const {
	fprintf(file, "[Triangle]:\n");
	for (int i = 0; i < 3; ++i)
		fprintf(file, "%f %f %f\n", vertices[i].x, vertices[i].y,
				vertices[i].z);

}
bool Triangle::Above(const Point3f& p, const Normal3f &n) const {
	return Dot(vertices[0] - p, n) >= 0 && Dot(vertices[1] - p, n) >= 0
			&& Dot(vertices[2] - p, n) >= 0;
}

}
