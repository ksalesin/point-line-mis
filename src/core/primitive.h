/*
 pbrt source code is Copyright(c) 1998-2016
 Matt Pharr, Greg Humphreys, and Wenzel Jakob.

 This file is part of pbrt.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 - Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 - Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_PRIMITIVE_H
#define PBRT_CORE_PRIMITIVE_H

// core/primitive.h*
#include "pbrt.h"
#include "shape.h"
#include "material.h"
#include "medium.h"
#include "point-line-mis/interval.h"
#include "point-line-mis/lineinteraction.h"
#include <mutex>

namespace pbrt {

//
using linesampling::LineVisibilityQuery3f;
using linesampling::IntervalCollection;

// Forward Declarations
class Primitive;

/*
 * Typedefs for a PrimitiveGroup
 *
 * The first element of the pair is a regular surface interaction.
 * The second element is a Primitive which is either:
 * => the intersected primitive of the surface interaction itself as a
 * 		std::shared_ptr (i.e. isect->primitive);
 * => a group of Primitive(s) which also contains the intersected primitive
 */
typedef std::pair<SurfaceInteraction, std::shared_ptr<const Primitive>> GroupIsect;
typedef std::vector<GroupIsect> GroupIsects;

/*
 * Utility typedef for shorter code.
 */
typedef std::vector<std::shared_ptr<const Primitive>> PrimitiveList;

// Comparator for sorting GroupIsects based on their ray parameter
struct GroupIsectComparator {
	inline bool operator()(const GroupIsect& i1, const GroupIsect& i2) {
		return i1.first.t < i2.first.t;
	}
};

// Enumeration of the proxies
enum ProxyType {
	INSIDE = 0, OUTSIDE = 1, MIXED = 2
};

// Primitive Declarations
class Primitive: public std::enable_shared_from_this<Primitive> {
public:
	// Primitive Interface
	Primitive() :
			id(++nextPrimitiveIdentifier) {
	}

	virtual ~Primitive();
	virtual Bounds3f WorldBound() const = 0;
	virtual bool Intersect(const Ray &r, SurfaceInteraction *) const = 0;
	virtual bool IntersectP(const Ray &r) const = 0;
	virtual bool IntersectP(const Ray &r, uint32_t ignoreId) const = 0;
	virtual bool IntersectL(const LineVisibilityQuery3f& query,
			const Normal3f &n, Float time,
			IntervalCollection& collection) const {
		LOG(FATAL)<< "Unimplemented Primitive::IntersectL!";
		return false;
	}
	virtual bool IntersectAll(const Ray& r, GroupIsects&,
			int groupSize = 1) const {
		LOG(FATAL)<< "Unimplemented Primitive::IntersectAll!";
		return false;
	}
	virtual const AreaLight *GetAreaLight() const = 0;
	virtual const Material *GetMaterial() const = 0;
	virtual void ComputeScatteringFunctions(SurfaceInteraction *isect,
			MemoryArena &arena, TransportMode mode,
			bool allowMultipleLobes) const = 0;
	virtual uint32_t Count() const = 0;
	virtual Float Area() const = 0;

	virtual std::shared_ptr<const Primitive> GetProxy(const ProxyType& type) const {
		return nullptr;
	}

	virtual uint32_t GetProxyId(const ProxyType& type) const {
		const std::shared_ptr<const Primitive>& proxy = GetProxy(type);
		if (proxy == nullptr) {
			return 0;
		}
		else {
			return proxy->id;
		}
	}

	const uint32_t id;
	static std::atomic<uint32_t> nextPrimitiveIdentifier;
};

// GeometricPrimitive Declarations
class GeometricPrimitive: public Primitive {
public:
	// GeometricPrimitive Public Methods
	virtual Bounds3f WorldBound() const;
	virtual bool Intersect(const Ray &r, SurfaceInteraction *isect) const;
	virtual bool IntersectP(const Ray &r) const;
	virtual bool IntersectP(const Ray &r, uint32_t ignoreId) const override {
		if (id == ignoreId)
			return false;
		return IntersectP(r);
	}
	virtual bool IntersectL(const LineVisibilityQuery3f& query,
			const Normal3f &n, Float time,
			IntervalCollection& collection) const;

	virtual bool IntersectAll(const Ray& r, GroupIsects&,
			int groupSize = 1) const override;
	GeometricPrimitive(const std::shared_ptr<Shape> &shape,
			const std::shared_ptr<Material> &material,
			const std::shared_ptr<AreaLight> &areaLight,
			const MediumInterface &mediumInterface) :
			shape(shape), material(material), areaLight(areaLight), mediumInterface(
					mediumInterface) {
	}
	const AreaLight *GetAreaLight() const;
	const Material *GetMaterial() const;
	void ComputeScatteringFunctions(SurfaceInteraction *isect,
			MemoryArena &arena, TransportMode mode,
			bool allowMultipleLobes) const;
	inline uint32_t Count() const override {
		return 1;
	}
	inline Float Area() const override {
		return shape->Area();
	}
private:
	// GeometricPrimitive Private Data
	std::shared_ptr<Shape> shape;
	std::shared_ptr<Material> material;
	std::shared_ptr<AreaLight> areaLight;
	MediumInterface mediumInterface;
};

// TransformedPrimitive Declarations
class TransformedPrimitive: public Primitive {
public:
	// TransformedPrimitive Public Methods
	TransformedPrimitive(const std::shared_ptr<const Primitive> &primitive,
			const AnimatedTransform &PrimitiveToWorld) :
			primitive(primitive), PrimitiveToWorld(PrimitiveToWorld) {
		proxies[ProxyType::INSIDE] = nullptr;
		proxies[ProxyType::OUTSIDE] = nullptr;
		proxies[ProxyType::MIXED] = nullptr;
	}

	// TransformedPrimitive Public Methods
	TransformedPrimitive(const std::shared_ptr<const Primitive> &primitive,
			const std::shared_ptr<const Primitive> &inside,
			const std::shared_ptr<const Primitive> &outside,
			const std::shared_ptr<const Primitive> &mixed,
			const AnimatedTransform &PrimitiveToWorld) :
			primitive(primitive), PrimitiveToWorld(PrimitiveToWorld) {
		proxies[ProxyType::INSIDE] = inside;
		proxies[ProxyType::OUTSIDE] = outside;
		proxies[ProxyType::MIXED] = mixed;
	}

	bool Intersect(const Ray &r, SurfaceInteraction *in) const;
	bool IntersectP(const Ray &r) const;
	bool IntersectP(const Ray &r, uint32_t ignoreId) const override;
	bool IntersectL(const LineVisibilityQuery3f& query, const Normal3f &n,
			Float time, IntervalCollection &intersection) const override;

	bool IntersectAll(const Ray& r, GroupIsects&, int groupSize = 1) const
			override;
	const AreaLight *GetAreaLight() const {
		return nullptr;
	}
	const Material *GetMaterial() const {
		return nullptr;
	}
	void ComputeScatteringFunctions(SurfaceInteraction *isect,
			MemoryArena &arena, TransportMode mode,
			bool allowMultipleLobes) const {
		LOG(FATAL)<<
		"TransformedPrimitive::ComputeScatteringFunctions() shouldn't be "
		"called";
	}
	Bounds3f WorldBound() const {
		return PrimitiveToWorld.MotionBounds(primitive->WorldBound());
	}

	virtual uint32_t Count() const override {
		return primitive->Count();
	}

	virtual std::shared_ptr<const Primitive> GetProxy(const ProxyType& type) const override {
		return proxies[type];
	}

	inline Float Area() const override;
private:
	// TransformedPrimitive Private Data
	std::shared_ptr<const Primitive> primitive;
	const AnimatedTransform PrimitiveToWorld;

	// List with the proxies
	std::shared_ptr<const Primitive> proxies[3];

	// A cache for the returned stuff
	mutable std::mutex mutex;
	mutable std::map<uint32_t, std::shared_ptr<const Primitive>> remapped;
};

// Aggregate Declarations
class Aggregate: public Primitive {
public:
	// Aggregate Public Methods
	Aggregate() :
			Primitive() {
	}

	const AreaLight *GetAreaLight() const;
	const Material *GetMaterial() const;
	void ComputeScatteringFunctions(SurfaceInteraction *isect,
			MemoryArena &arena, TransportMode mode,
			bool allowMultipleLobes) const;
};

}  // namespace pbrt

#endif  // PBRT_CORE_PRIMITIVE_H
