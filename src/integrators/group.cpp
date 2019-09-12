#include "group.h"
#include "paramset.h"
#include "camera.h"
#include "film.h"

namespace pbrt {

static uint32_t hash(uint32_t key) {
	uint32_t hash = key;
	hash += (hash << 10);
	hash ^= (hash >> 6);
	hash += (hash << 3);
	hash ^= (hash >> 11);
	hash += (hash << 15);
	return hash;
}

Spectrum GroupIntegrator::Li(const RayDifferential &ray, const Scene &scene,
		Sampler &sampler, MemoryArena &arena, int depth) const {
	ProfilePhase p(Prof::SamplerIntegratorLi);

	// Find closest ray intersection or return background radiance
	SurfaceInteraction isect;
	GroupIsects isects;
	if (!scene.IntersectAll(ray, isects, groupSize)) {
		return Spectrum(0.f);
	}

	/* Sort the intersections */
	std::sort(isects.begin(), isects.end(),
			[]( const GroupIsect& lhs, const GroupIsect& rhs) {
				return lhs.first.t < rhs.first.t;
			});

	uint32_t id = isects[0].second->id;
	uint32_t h = hash(id);

	Float r = Float(h & 0xff) / 255.0;
	Float g = Float((h >> 8) & 0xff) / 255.0;
	Float b = Float((h >> 16) & 0xff) / 255.0;

	return Spectrum::FromRGB(r, g, b);
}

GroupIntegrator *CreateGroupIntegrator(const ParamSet &params,
		std::shared_ptr<Sampler> sampler,
		std::shared_ptr<const Camera> camera) {

	int np;
	const int *pb = params.FindInt("pixelbounds", &np);
	Bounds2i pixelBounds = camera->film->GetSampleBounds();
	if (pb) {
		if (np != 4)
			Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
					np);
		else {
			pixelBounds = Intersect(pixelBounds, Bounds2i { { pb[0], pb[2] }, {
					pb[1], pb[3] } });
			if (pixelBounds.Area() == 0)
				Error("Degenerate \"pixelbounds\" specified.");
		}
	}

	int groupSize = std::max(1, params.FindOneInt("groupsize", 1));
	return new GroupIntegrator(camera, sampler, pixelBounds, groupSize);
}
}
