#include "tests/gtest/gtest.h"

#include "pbrt.h"
#include "point-line-mis/triangle.h"

using namespace pbrt;

TEST(TriangleAbove, TriangleClipping) {

	niels::Triangle t0(Point3f(0, 0, 1), Point3f(1, 0, 1), Point3f(1, 1, 1));

	std::vector<niels::Triangle> clipped;
	t0.Clip(Point3f(0, 0, 0), Normal3f(0, 0, 1), clipped);

	EXPECT_EQ(1u, clipped.size());
	EXPECT_EQ(t0[0], clipped[0][0]);
	EXPECT_EQ(t0[1], clipped[0][1]);
	EXPECT_EQ(t0[2], clipped[0][2]);
}

TEST(TriangleOne, TriangleClipping) {
	niels::Triangle t0(Point3f(0, 0, 0), Point3f(1, 0, 1), Point3f(1, 0, -1));

	std::vector<niels::Triangle> clipped;
	t0.Clip(Point3f(0, 0, 0), Normal3f(0, 0, 1), clipped);

	EXPECT_EQ(1u, clipped.size());
	const niels::Triangle clippedTriangle = clipped[0];
	EXPECT_EQ(true, clippedTriangle.HasVertex(Point3f(0, 0, 0), 1e-6));
	EXPECT_EQ(true, clippedTriangle.HasVertex(Point3f(1, 0, 1), 1e-6));
	EXPECT_EQ(true, clippedTriangle.HasVertex(Point3f(1, 0, 0), 1e-6));
}

TEST(TriangleTwo, TriangleClipping) {
	niels::Triangle t0(Point3f(0, 0, -1), Point3f(4, 0, -1), Point3f(2, 0, 1));

	std::vector<niels::Triangle> clipped;
	t0.Clip(Point3f(0, 0, 0), Normal3f(0, 0, 1), clipped);

	EXPECT_EQ(1u, clipped.size());
	const niels::Triangle clippedTriangle = clipped[0];
	EXPECT_EQ(true, clippedTriangle.HasVertex(Point3f(1, 0, 0), 1e-6));
	EXPECT_EQ(true, clippedTriangle.HasVertex(Point3f(3, 0, 0), 1e-6));
	EXPECT_EQ(true, clippedTriangle.HasVertex(Point3f(2, 0, 1), 1e-6));
}

TEST(TriangleThree, TriangleClipping) {
	niels::Triangle t0(Point3f(0, 0, 1), Point3f(4, 0, 1), Point3f(2, 0, -1));

	std::vector<niels::Triangle> clipped;
	t0.Clip(Point3f(0, 0, 0), Normal3f(0, 0, 1), clipped);

	EXPECT_EQ(2u, clipped.size());
	const niels::Triangle& firstTriangle = clipped[0];
	const niels::Triangle& secondTriangle = clipped[1];

	firstTriangle.Print(stdout);
	secondTriangle.Print(stdout);

	EXPECT_EQ(true, firstTriangle.HasVertex(Point3f(4, 0, 1), 1e-6));
	EXPECT_EQ(true, firstTriangle.HasVertex(Point3f(1, 0, 0), 1e-6));
	EXPECT_EQ(true, firstTriangle.HasVertex(Point3f(0, 0, 1), 1e-6));

	EXPECT_EQ(true, secondTriangle.HasVertex(Point3f(4, 0, 1), 1e-6));
	EXPECT_EQ(true, secondTriangle.HasVertex(Point3f(3, 0, 0), 1e-6));
	EXPECT_EQ(true, secondTriangle.HasVertex(Point3f(1, 0, 0), 1e-6));

}
