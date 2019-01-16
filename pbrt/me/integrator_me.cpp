
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

// core/integrator.cpp*
#include "integrator.h"
#include "scene.h"
#include "interaction.h"
#include "sampling.h"
#include "parallel.h"
#include "film.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "stats.h"

namespace pbrt {

STAT_COUNTER("Integrator/Camera rays traced", nCameraRays);


// Integrator Utility Functions
Spectrum UniformSampleAllLights(const Interaction &it, const Scene &scene,
	MemoryArena &arena, Sampler &sampler,
	const std::vector<int> &nLightSamples,
	bool handleMedia) {
	ProfilePhase p(Prof::DirectLighting);
	Spectrum L(0.f);
	for (size_t j = 0; j < scene.lights.size(); ++j) {
		const std::shared_ptr<Light>& light = scene.lights[j];
		int nSamples = nLightSamples[j];
		const Point2f *uLightArray = sampler.Get2DArray(nSamples);
		const Point2f *uScatteringArray = sampler.Get2DArray(nSamples);
		if (!uLightArray || !uScatteringArray) {
			//if sampler array use done, then get single sample point.
			Point2f uLight = sampler.Get2D();
			Point2f uScatter = sampler.Get2D();
			L += EstimateDirect(it, uScatter, *light, uLight, scene, sampler, arena, handleMedia);
		} else {
			Spectrum Ld(0.f);
			for (size_t k = 0; k < nSamples; ++k) {
				Ld += EstimateDirect(it, uScatteringArray[k], *light, uLightArray[k], scene, sampler, arena, handleMedia);
			}
			L += Ld / nSamples;
		}
	}
	return L;
}

Spectrum UniformSampleOneLight(const Interaction &it, const Scene &scene,
	MemoryArena &arena, Sampler &sampler,
	bool handleMedia, const Distribution1D *lightDistrib) {
	ProfilePhase p(Prof::DirectLighting);
	// Randomly choose a single light to sample, _light_
	int nLights = int(scene.lights.size());
	if (nLights == 0) return Spectrum(0.f);
	int lightNum;
	Float lightPdf;
	if (lightDistrib) {
		lightNum = lightDistrib->SampleDiscrete(sampler.Get1D(), &lightPdf);
		if (lightPdf == 0) return Spectrum(0.f);
	} else {
		lightNum = std::min((int)(sampler.Get1D() * nLights), nLights - 1);
		lightPdf = Float(1) / nLights;
	}
	const std::shared_ptr<Light>& light = scene.lights[lightNum];
	Point2f uScattering = sampler.Get2D();
	Point2f uLight = sampler.Get2D();
	return EstimateDirect(it, uScattering, *light, uLight,
		scene, sampler, arena, handleMedia) / lightPdf;
}
Spectrum EstimateDirect(const Interaction &it, const Point2f &uScattering,
                        const Light &light, const Point2f &uLight,
                        const Scene &scene, Sampler &sampler,
                        MemoryArena &arena, bool handleMedia, bool specular) {
	BxDFType bsdfFlags =
		specular ? BSDF_ALL : BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
	
	Vector3f wi;
	Float lightPdf = 0;
	Float scatteringPdf = 0;
	VisibilityTester vis;

	int nSampleLight = 1;
	Spectrum Ld(0.f);
	for (int lightSample = 0; lightSample < nSampleLight; ++lightSample) {
		Spectrum Li = light.Sample_Li(it, uLight, &wi, &lightPdf, &vis);



		if (lightPdf > 0.f && !Li.IsBlack()) {
			Spectrum f;
			if (it.IsSurfaceInteraction()) {
				const SurfaceInteraction& isect = (const SurfaceInteraction &)(it);
				f = isect.bsdf->f(isect.wo, wi, bsdfFlags) * AbsDot(wi, isect.shading.n);
				scatteringPdf = isect.bsdf->Pdf(isect.wo, wi, bsdfFlags);
			} else {
				const MediumInteraction &mi = (const MediumInteraction &)it;
				Float p = mi.phase->p(mi.wo, wi);
				f = Spectrum(p);
				scatteringPdf = p;
			}
			if (!f.IsBlack()) {
				//visibility
				if (handleMedia) {
					Li *= vis.Tr(scene, sampler);
				} else {
					if (!vis.Unoccluded(scene)) {
						Li = Spectrum(0.f);
						//加入头发的透射.
						//Li *= vis.Tr(scene);
					} else {

					}
				}

				if (!Li.IsBlack()) {
					if (IsDeltaLight(light.flags)) {
						Ld += f * Li / lightPdf;
					} else {
						//mis first step , sample light
						Float weight = PowerHeuristic(1, lightPdf, 1, scatteringPdf);
						Ld += f * Li * weight / lightPdf;
					}
				}
			}
		}

		//sample bsdf for mis
		if (!IsDeltaLight(light.flags)) {
			Spectrum f;
			bool sampledSpecular = false;
			if (it.IsSurfaceInteraction()) {
				//sample bsdf direction
				BxDFType sampledType;
				const SurfaceInteraction& isect = (const SurfaceInteraction&)(it);
				f = isect.bsdf->Sample_f(isect.wo, &wi, uScattering, &scatteringPdf, bsdfFlags, &sampledType);
				f *= AbsDot(wi, isect.shading.n);
				sampledSpecular = (sampledType & BSDF_SPECULAR) != 0;
			} else {
				const MediumInteraction &mi = (const MediumInteraction &)it;
				Float p = mi.phase->Sample_p(mi.wo, &wi, uScattering);
				f = Spectrum(p);
				scatteringPdf = p;
			}

			if (!f.IsBlack() && scatteringPdf > 0.f) {
				Float weight = 1;
				if (!sampledSpecular) {
					lightPdf = light.Pdf_Li(it, wi);
					if (lightPdf == 0) {
						return Ld;
					}
					weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
				}

				SurfaceInteraction lightIsect;
				Ray ray = it.SpawnRay(wi);
				Spectrum Tr(1.f);
				bool foundSurfaceInteracton = handleMedia ?
					scene.IntersectTr(ray, sampler, &lightIsect, &Tr) :
					scene.Intersect(ray, &lightIsect);

				Spectrum Li(0.f);
				if (foundSurfaceInteracton) {
					if (lightIsect.primitive->GetAreaLight() == &light) {
						Li = lightIsect.Le(-wi);
					}
				} else {
					Li = light.Le(ray);
				}
				if (!f.IsBlack()) {
					Ld += Li * f * weight * Tr / scatteringPdf;
				}
			}
		}
	}
	Ld /= static_cast<Float>(nSampleLight);

    return Ld;
}


Spectrum SamplerIntegrator::SpecularReflect(
    const RayDifferential &ray, const SurfaceInteraction &isect,
    const Scene &scene, Sampler &sampler, MemoryArena &arena, int depth) const {
    // Compute specular reflection direction _wi_ and BSDF value
	Vector3f wo = isect.wo;
	Vector3f wi;
	Float pdf;
	BxDFType type = BxDFType(BSDF_REFLECTION | BSDF_SPECULAR);
	Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, type);
	const Normal3f &ns = isect.shading.n;
	if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
		RayDifferential rd = isect.SpawnRay(wi);
		if (ray.hasDifferentials) {

		}
		return f * Li(rd, scene, sampler, arena, depth + 1) * AbsDot(wi, ns) / pdf;
	} else {
		return Spectrum(0.f);
	}
}

Spectrum SamplerIntegrator::SpecularTransmit(
    const RayDifferential &ray, const SurfaceInteraction &isect,
    const Scene &scene, Sampler &sampler, MemoryArena &arena, int depth) const {
    Vector3f wo = isect.wo;
	Vector3f wi;
	const Normal3f ns = isect.shading.n;
	BxDFType type = BxDFType(BSDF_SPECULAR | BSDF_TRANSMISSION);
	Float pdf;
	Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, type);
	Spectrum L(0.f);
	if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
		RayDifferential rd = isect.SpawnRay(wi);
		if (ray.hasDifferentials) {

		}
		L = f * Li(rd, scene, sampler, arena, depth + 1) * AbsDot(wi, ns) / pdf;
	}

	return L;
}


void SamplerIntegrator::Render(const Scene &scene) {
	Preprocess(scene, *sampler);
	// Render image tiles in parallel

	// Compute number of tiles, _nTiles_, to use for parallel rendering
	Bounds2i sampleBounds = camera->film->GetSampleBounds();
	Vector2i sampleExtent = sampleBounds.Diagonal();
	const int tileSize = 16;
	Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
		(sampleExtent.y + tileSize - 1) / tileSize);
	ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
	{
		ParallelFor2D([&](Point2i tile) {
			// Render section of image corresponding to _tile_

			// Allocate _MemoryArena_ for tile
			MemoryArena arena;

			// Get sampler instance for tile
			int seed = tile.y * nTiles.x + tile.x;
			std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);

			// Compute sample bounds for tile
			int x0 = sampleBounds.pMin.x + tile.x * tileSize;
			int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
			int y0 = sampleBounds.pMin.y + tile.y * tileSize;
			int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
			Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
			LOG(INFO) << "Starting image tile " << tileBounds;

			// Get _FilmTile_ for tile
			std::unique_ptr<FilmTile> filmTile =
				camera->film->GetFilmTile(tileBounds);

			// Loop over pixels in tile to render them
			for (Point2i pixel : tileBounds) {
				{
					ProfilePhase pp(Prof::StartPixel);
					tileSampler->StartPixel(pixel);
				}

				// Do this check after the StartPixel() call; this keeps
				// the usage of RNG values from (most) Samplers that use
				// RNGs consistent, which improves reproducability /
				// debugging.
				if (!InsideExclusive(pixel, pixelBounds))
					continue;

				do {
					// Initialize _CameraSample_ for current sample
					CameraSample cameraSample =
						tileSampler->GetCameraSample(pixel);

					// Generate camera ray for current sample
					//                     RayDifferential ray;
					//                     Float rayWeight =
					//                         camera->GenerateRayDifferential(cameraSample, &ray);
					//                     ray.ScaleDifferentials(
					//                         1 / std::sqrt((Float)tileSampler->samplesPerPixel));
					Ray ray;
					Float rayWeight =
						camera->GenerateRay(cameraSample, &ray);

					++nCameraRays;

					// Evaluate radiance along camera ray
					Spectrum L(0.f);


					if (rayWeight > 0) L = Li(ray, scene, *tileSampler, arena);
					/************************************************************************/
					/* test volume light                                                    */
					/************************************************************************/
					//if (rayWeight > 0) L = VolumeLiRadiance(ray, scene, *tileSampler, arena);
					//if (rayWeight > 0) L = VolumeLiRadiance(ray, scene, *tileSampler, arena);

					// Issue warning if unexpected radiance value returned
					if (L.HasNaNs()) {
						LOG(ERROR) << StringPrintf(
							"Not-a-number radiance value returned "
							"for pixel (%d, %d), sample %d. Setting to black.",
							pixel.x, pixel.y,
							(int)tileSampler->CurrentSampleNumber());
						L = Spectrum(0.f);
					} else if (L.y() < -1e-5) {
						LOG(ERROR) << StringPrintf(
							"Negative luminance value, %f, returned "
							"for pixel (%d, %d), sample %d. Setting to black.",
							L.y(), pixel.x, pixel.y,
							(int)tileSampler->CurrentSampleNumber());
						L = Spectrum(0.f);
					} else if (std::isinf(L.y())) {
						LOG(ERROR) << StringPrintf(
							"Infinite luminance value returned "
							"for pixel (%d, %d), sample %d. Setting to black.",
							pixel.x, pixel.y,
							(int)tileSampler->CurrentSampleNumber());
						L = Spectrum(0.f);
					}
					VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
						ray << " -> L = " << L;

					// Add camera ray's contribution to image
					filmTile->AddSample(cameraSample.pFilm, L, rayWeight);

					// Free _MemoryArena_ memory from computing image sample
					// value
					arena.Reset();
				} while (tileSampler->StartNextSample());
			}
			LOG(INFO) << "Finished image tile " << tileBounds;

			// Merge image tile into _Film_
			camera->film->MergeFilmTile(std::move(filmTile));
			reporter.Update();
		}, nTiles);
		reporter.Done();
	}
	LOG(INFO) << "Rendering finished";

	// Save final image after rendering
	camera->film->WriteImage();
}


}  // namespace pbrt
