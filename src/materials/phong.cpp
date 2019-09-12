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

// materials/matte.cpp*
#include "materials/phong.h"
#include "paramset.h"
#include "reflection.h"
#include "interaction.h"
#include "texture.h"
#include "interaction.h"

namespace pbrt {

// MatteMaterial Method Definitions
void PhongMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
		MemoryArena &arena, TransportMode mode, bool allowMultipleLobes) const {
	// Perform bump mapping with _bumpMap_, if present
	if (bumpMap)
		Bump(bumpMap, si);

	// Evaluate textures for _MatteMaterial_ material and allocate BRDF
	si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
	Spectrum diffuse = Kd->Evaluate(*si).Clamp();
	Spectrum phong = Ks->Evaluate(*si).Clamp();

	if (!diffuse.IsBlack())
		si->bsdf->Add(ARENA_ALLOC(arena, LambertianReflection)(diffuse));
	if (!phong.IsBlack())
		si->bsdf->Add(
		ARENA_ALLOC(arena, PhongReflection)(phong, n, cache_n, cache_n2));
}

PhongMaterial *CreatePhongMaterial(const TextureParams &mp) {
	std::shared_ptr<Texture<Spectrum>> Kd = mp.GetSpectrumTexture("Kd",
			Spectrum(0.5f));
	std::shared_ptr<Texture<Spectrum>> Ks = mp.GetSpectrumTexture("Ks",
			Spectrum(0.25f));
	std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull(
			"bumpmap");
	int n = mp.FindInt("n", mp.FindInt("exponent", 10));
	int cacheSize = mp.FindInt("cachesize", mp.FindInt("cache", 0));
	return new PhongMaterial(Kd, Ks, bumpMap, n, cacheSize);
}

}  // namespace pbrt
