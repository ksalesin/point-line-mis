/*
 * objmesh.cpp
 *
 *  Created on: Jul 20, 2017
 *      Author: niels
 */

#include "objmesh.h"
#include "paramset.h"
#include "textures/constant.h"

#include <fstream>

namespace pbrt {

typedef std::map<std::string, std::shared_ptr<Texture<Float>>>TextureMaps;

class OBJVertex {
public:
	OBJVertex() {
	}
	OBJVertex(const Point3f& point, const Normal3f& normal,
			const Point2f& texture) :
			point(point), normal(normal), texture(texture) {
	}

	bool operator<(const OBJVertex& vertex) const {
		if (point.x < vertex.point.x)
			return true;
		if (point.x > vertex.point.x)
			return false;

		if (point.y < vertex.point.y)
			return true;
		if (point.y > vertex.point.y)
			return false;

		if (point.z < vertex.point.z)
			return true;
		if (point.z > vertex.point.z)
			return false;

		if (normal.x < vertex.normal.x)
			return true;
		if (normal.x > vertex.normal.x)
			return false;

		if (normal.y < vertex.normal.y)
			return true;
		if (normal.y > vertex.normal.y)
			return false;

		if (normal.z < vertex.normal.z)
			return true;
		if (normal.z > vertex.normal.z)
			return false;

		if (texture.x < vertex.texture.x)
			return true;
		if (texture.x > vertex.texture.x)
			return false;

		if (texture.y < vertex.texture.y)
			return true;

		return false;
	}

	Point3f point;
	Normal3f normal;
	Point2f texture;
};

std::vector<std::shared_ptr<Shape>> CreateOBJMesh(const Transform *o2w,
		const Transform *w2o, bool reverseOrientation, const ParamSet &params,
		TextureMaps *floatTextures) {
	const std::string filename = params.FindOneFilename("filename", "");
	std::vector<Point3f> objVertices;
	std::vector<Normal3f> objNormals;
	std::vector<Point2f> objTextureCoordinates;
	std::vector<std::shared_ptr<Shape>> result;

	std::vector<OBJVertex> meshVertices;
	std::map<OBJVertex, int> meshMapping;
	std::vector<int> meshIndices;

	std::string token;
	std::ifstream stream(filename);

	while (stream >> token) {
		if (token == "v") {
			/* Read a vertex */
			std::string xs, ys, zs;
			stream >> xs >> ys >> zs;
			Point3f p;
			p.x = Float(std::atof(xs.c_str()));
			p.y = Float(std::atof(ys.c_str()));
			p.z = Float(std::atof(zs.c_str()));
			objVertices.push_back(p);
		} else if (token == "vn") {
			/* Read a vertex */
			std::string xs, ys, zs;
			stream >> xs >> ys >> zs;
			Normal3f vn;
			vn.x = Float(std::atof(xs.c_str()));
			vn.y = Float(std::atof(ys.c_str()));
			vn.z = Float(std::atof(zs.c_str()));
			objNormals.push_back(vn);
		} else if (token == "vt") {
			/* Read a vertex */
			std::string us, vs;
			stream >> us >> vs;
			Point2f vt;
			vt.x = Float(std::atof(us.c_str()));
			vt.y = Float(std::atof(vs.c_str()));
			objTextureCoordinates.push_back(vt);
		} else if (token == "f") {
			for (int i = 0; i < 3; ++i) {
				stream >> token;

				int first = token.find("/");
				int second = token.find("/", first + 1);

				int i0 = std::atoi(token.substr(0, first).c_str());
				int i1 = std::atoi(token.substr(first + 1, second).c_str());
				int i2 = std::atoi(token.substr(second + 1).c_str());

				const Point3f& point = objVertices[i0 - 1];
				const Point2f& texture = objTextureCoordinates[i1 - 1];
				const Normal3f& normal = objNormals[i2 - 1];

				OBJVertex vertex(point, normal, texture);

				const auto iterator = meshMapping.find(vertex);
				if (iterator != meshMapping.end())
					meshIndices.push_back(iterator->second);
				else {
					int index = meshVertices.size();
					meshIndices.push_back(index);
					meshMapping.insert(std::make_pair(vertex, index));
					meshVertices.push_back(vertex);
				}

//				fprintf(stdout, "%s -> %d %d %d\n", token.c_str(), i0, i1, i2);

			}
		}
	}

	const int nTriangles = meshIndices.size() / 3;
	const int vertexCount = meshVertices.size();
	Point3f* p = new Point3f[vertexCount];
	Normal3f* n = new Normal3f[vertexCount];
	Point2f* uv = new Point2f[vertexCount];

	for (int i = 0; i < vertexCount; ++i) {
		p[i] = meshVertices[i].point;
		n[i] = meshVertices[i].normal;
		uv[i] = meshVertices[i].texture;
	}

	std::shared_ptr<Texture<Float>> alphaTex;
	std::string alphaTexName = params.FindTexture("alpha");
	if (alphaTexName != "") {
		if (floatTextures->find(alphaTexName) != floatTextures->end())
			alphaTex = (*floatTextures)[alphaTexName];
		else
			Error("Couldn't find float texture \"%s\" for \"alpha\" parameter",
					alphaTexName.c_str());
	} else if (params.FindOneFloat("alpha", 1.f) == 0.f) {
		alphaTex.reset(new ConstantTexture<Float>(0.f));
	}

	std::shared_ptr<Texture<Float>> shadowAlphaTex;
	std::string shadowAlphaTexName = params.FindTexture("shadowalpha");
	if (shadowAlphaTexName != "") {
		if (floatTextures->find(shadowAlphaTexName) != floatTextures->end())
			shadowAlphaTex = (*floatTextures)[shadowAlphaTexName];
		else
			Error("Couldn't find float texture \"%s\" for \"shadowalpha\" "
					"parameter", shadowAlphaTexName.c_str());
	} else if (params.FindOneFloat("shadowalpha", 1.f) == 0.f)
		shadowAlphaTex.reset(new ConstantTexture<Float>(0.f));

	std::vector<std::shared_ptr<Shape>> results = CreateTriangleMesh(o2w, w2o,
			reverseOrientation, nTriangles, &meshIndices[0], vertexCount, p,
			nullptr, n, uv, alphaTex, shadowAlphaTex);

	delete[] p;
	delete[] n;
	delete[] uv;

	return results;

}

}
