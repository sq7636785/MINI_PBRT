
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
#include "me/matte.h"
#include "paramset.h"
#include "reflection.h"
#include "interaction.h"
#include "texture.h"
#include "interaction.h"

namespace pbrt {

// MatteMaterial Method Definitions
	void MatteMaterial::ComputeScatteringFunctions(SurfaceInteraction *si, 
		                                           MemoryArena &arena, 
		                                           TransportMode mode, 
		                                           bool allowMultipleLobes) const {
		//bump mapping if present
		if (bumpMap) { Bump(bumpMap, si); }

		//evaluate parameters and allocate bsdf
		si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
		Spectrum r= Kd->Evaluate(*si).Clamp();
		//here is roughness, if not equal to zero, then lambertianReflect, else OrenNayar
		Float sig = Clamp(sigma->Evaluate(*si), 0, 90);
		if (!r.IsBlack()) {
			if (sig == 0) {
				si->bsdf->Add(ARENA_ALLOC(arena, LambertianReflection)(r));
			} else {
				// not implementation, should be diffuse with roughness.
				si->bsdf->Add(ARENA_ALLOC(arena, LambertianReflection)(r));
			}
		}
	}

	MatteMaterial* CreateMatteMaterial(const TextureParams &mp) {
		std::shared_ptr<Texture<Spectrum>> Kd = mp.GetSpectrumTexture("Kd", Spectrum(0.5f));
		std::shared_ptr<Texture<Float>> sigma = mp.GetFloatTexture("sigma", 0.f);
		std::shared_ptr<Texture<Float>> bump = mp.GetFloatTextureOrNull("bumpmap");
		return new MatteMaterial(Kd, sigma, bump);
	}


}  // namespace pbrt
