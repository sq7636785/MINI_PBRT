
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


// integrators/directlighting.cpp*
#include "me/directlighting.h"
#include "interaction.h"
#include "paramset.h"
#include "camera.h"
#include "film.h"
#include "stats.h"
#include "sampling.h"

namespace pbrt {

	// DirectLightingIntegrator Method Definitions
	void DirectLightingIntegrator::Preprocess(const Scene &scene,
		Sampler &sampler) {
		if (strategy == LightStrategy::UniformSampleAll) {
			// Compute number of samples to use for each light
			for (const auto &light : scene.lights)
				nLightSamples.push_back(sampler.RoundCount(light->nSamples));

			// Request samples for sampling all lights
			for (int i = 0; i < maxDepth; ++i) {
				for (size_t j = 0; j < scene.lights.size(); ++j) {
					sampler.Request2DArray(nLightSamples[j]);
					sampler.Request2DArray(nLightSamples[j]);
				}
			}
		}
	}

	Spectrum DirectLightingIntegrator::Li(const RayDifferential &ray,
		const Scene &scene, Sampler &sampler,
		MemoryArena &arena, int depth) const {
		ProfilePhase p(Prof::SamplerIntegratorLi);
		Spectrum L(0.f);

		//no intersection, calculate 
		SurfaceInteraction isect;
		if (!scene.Intersect(ray, &isect)) {
			for (const auto &light : scene.lights) {
				L += light->Le(ray);
			}
			return L;
		}

		isect.ComputeScatteringFunctions(ray, arena);
		if (!isect.bsdf) {
			return Li(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
		}

		Vector3f wo = isect.wo;
		//compute  emitted light if hit an area light
		L += isect.Le(wo);
		if (scene.lights.size() > 0) {
			if (strategy == LightStrategy::UniformSampleAll) {
				L += UniformSampleAllLights(isect, scene, arena, sampler, nLightSamples);
			} else {
				L += UniformSampleOneLight(isect, scene, arena, sampler);
			}
		}
		// 	if (depth + 1 < maxDepth) {
		// 		L += SpecularReflect(ray, isect, scene, sampler, arena, depth);
		// 		L += SpecularTransmit(ray, isect, scene, sampler, arena, depth);
		// 	}

		return L;
	}


	//volume Li
	Spectrum DirectLightingIntegrator::VolumeLi(const RayDifferential& ray, const Scene& scene, Sampler& sampler, MemoryArena &arena, int depth) const {
		ProfilePhase p(Prof::SamplerIntegratorLi);
		Spectrum L(0.f);

		//no intersection, calculate 
		SurfaceInteraction isect;
		if (!scene.Intersect(ray, &isect)) {
			for (const auto &light : scene.lights) {
				L += light->Le(ray);
			}
			return L;
		}
		BxDFType bsdfFlags = BxDFType(BSDF_ALL & ~BSDF_SPECULAR);

		isect.ComputeScatteringFunctions(ray, arena);
		if (!isect.bsdf) {
			return VolumeLi(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
		}
		L += isect.Le(isect.wo);
		//目前仅支持单光源
		for (size_t i = 0; i < scene.lights.size(); ++i) {
			const std::shared_ptr<Light> light = scene.lights[i];
			Vector3f wi;
			int idx = scene.volume->GetIdxFromPoint(isect.p.x, isect.p.y, isect.p.z);
			const Voxel& v = scene.volume->voxel[idx];

			Point3f p = (v.bound.pMin + v.bound.pMax) / 2.0;
			light->SampleWi(p, &wi);
			Spectrum f = isect.bsdf->f(isect.wo, wi, bsdfFlags) * AbsDot(wi, isect.shading.n);;
			if (!f.IsBlack()) {
				Spectrum irrandiance = RGBSpectrum::FromRGB(v.rgb);
				L += irrandiance * f;
				// 			if (L.y() > 1.0) {
				// 				std::cout << std::endl;
				// 				std::cout << "irrandiance: "<< irrandiance << std::endl;
				// 				std::cout << "f: " <<  f << std::endl;
				// 				std::cout <<  "L: " << L << std::endl;
				// 				std::cout << "Wi: " << wi << std::endl;
				// 				std::cout << "Wo: " << isect.wo << std::endl;
				// 				std::cout << idx << std::endl;
				// 			}
			}

		}
		return L;
	}

	Spectrum DirectLightingIntegrator::VolumeLiRadiance(const RayDifferential& ray, const Scene& scene, Sampler& sampler, MemoryArena &arena, int depth) const {
		ProfilePhase p(Prof::SamplerIntegratorLi);
		Spectrum L(0.f);

		//no intersection, calculate 
		SurfaceInteraction isect;
		if (!scene.Intersect(ray, &isect)) {
			for (const auto &light : scene.lights) {
				L += light->Le(ray);
			}
			return L;
		}
		BxDFType bsdfFlags = BxDFType(BSDF_ALL & ~BSDF_SPECULAR);

		isect.ComputeScatteringFunctions(ray, arena);
		if (!isect.bsdf) {
			return VolumeLiRadiance(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
		}
		L += isect.Le(isect.wo);
		
		//reconstruct radiance spherical function
		int idx = scene.volume->GetIdxFromPoint(isect.p.x, isect.p.y, isect.p.z);
		const Voxel& v = scene.volume->voxel[idx];




//#define UNIFORM_SAMPLE
#define HAIR_SAMPLE
//#define LIGHT_SAMPLE
//#define BSDF_MATRIX
#ifdef BSDF_MATRIX
#pragma region _BSDF_Matrix
		std::vector<Spectrum> cLightOri = v.shC;
		std::vector<Spectrum> cLihgtRot(SHTerms(scene.volume->shL), Spectrum(0.f));
		std::vector<Spectrum> dOutR(SHTerms(scene.volume->shL), Spectrum(0.f));

		Vector3f globalN(isect.n.x, isect.n.y, isect.n.z);
		Vector3f localN(0.f, 0.f, 1.f);
		scene.volume->RotateSH(globalN, localN, cLightOri.data(), cLihgtRot.data());
		scene.volume->SHMatrixTransV(cLihgtRot, dOutR.data());
		//reconstruct
		Spectrum outRadiance(0.f);
		std::vector<Float> ylm(SHTerms(scene.volume->shL));
		SHEvaluate(isect.wo, scene.volume->shL, ylm.data());

		for (int i = 0; i < SHTerms(scene.volume->shL); ++i) {
			outRadiance += dOutR[i] * ylm[i];
		}
		if (outRadiance.y() > 0.f) {
			L += outRadiance;
		}
#endif
#pragma endregion

#ifdef UNIFORM_SAMPLE
#pragma region _unifomSample
		bool isReconstuct = false;
		for (const auto& c : v.shC) {
			if (c != 0.f) {
				isReconstuct = true;
				break;
			}
		}
		if (isReconstuct) {
			//reconstruct radiance
			//std::cout << "reconstruct" << std::endl;

			Spectrum mcEstimator(0.f);
			for (int j = 0; j < scene.volume->nSHSample; ++j) {
				Vector3f wi = scene.volume->shSample[j].w;
				Spectrum radiance(0.f);
				for (int i = 0; i < SHTerms(scene.volume->shL); ++i) {
					radiance += v.shC[i] * scene.volume->shSample[j].y[i];
				}
				if (radiance.y() > 0.f) {
					Spectrum f = isect.bsdf->f(isect.wo, wi, bsdfFlags) * AbsDot(wi, isect.shading.n);
					if (!f.IsBlack()) {
						mcEstimator += radiance * f;
					}
				}
			}
			Float scale = 2 * Pi / static_cast<Float>(scene.volume->nSHSample);
			mcEstimator *= scale;
			L += mcEstimator;
			//std::cout << L << std::endl;
		}
#pragma endregion
#endif

#ifdef LIGHT_SAMPLE
#pragma region _sampleLight
		int nSampleSphere = 2500;
		RNG rng;
		Spectrum mcEestimate(0.f);
		for (int i = 0; i < nSampleSphere; ++i) {
			Point2f u(rng.UniformFloat(), rng.UniformFloat());
			Vector3f wi = UniformSampleHemisphere(u);
			Float pdf = UniformHemispherePdf();
			std::vector<Float> ylm(SHTerms(scene.volume->shL));
			SHEvaluate(wi, scene.volume->shL, ylm.data());
			Spectrum radiance(0.f);
			for (int k = 0; k < SHTerms(scene.volume->shL); ++k) {
				radiance += v.shC[k] * ylm[k];
			}
			if (radiance.y() > 0.f) {
				Spectrum f = isect.bsdf->f(isect.wo, wi, bsdfFlags) * AbsDot(wi, isect.shading.n);
				mcEestimate += radiance * f / (pdf * nSampleSphere);
			}
		}
		L += mcEestimate;
#pragma endregion
#endif

#ifdef HAIR_SAMPLE		
#pragma region _hariSample
		int nSampleHair = 50;
		RNG rng;
		Spectrum mcEestimate(0.f);
		for (int i = 0; i < nSampleHair; ++i) {
			Point2f u(rng.UniformFloat(), rng.UniformFloat());
			Vector3f wi;
			Float pdf;
			Spectrum f = isect.bsdf->Sample_f(isect.wo, &wi, u, &pdf, bsdfFlags);//* AbsDot(wi, isect.shading.n);
			std::vector<Float> ylm(SHTerms(scene.volume->shL));
			SHEvaluate(wi, scene.volume->shL, ylm.data());
			Spectrum radiance(0.f);
			for (int k = 0; k < SHTerms(scene.volume->shL); ++k) {
				radiance += v.shC[k] * ylm[k];
			}
			if (radiance.y() > 0.f) {

				mcEestimate += radiance * f * AbsDot(wi, isect.shading.n) / pdf;
			}
		}
		L += (mcEestimate / nSampleHair);
#pragma endregion
#endif
		return L;
	}



	Spectrum DirectLightingIntegrator::VolumeSHMatrixWo(const RayDifferential& ray, const Scene& scene, Sampler& sampler, MemoryArena &arena, int depth) const {
		return Spectrum(0.f);
	}


	DirectLightingIntegrator *CreateDirectLightingIntegrator(
		const ParamSet &params, std::shared_ptr<Sampler> sampler,
		std::shared_ptr<const Camera> camera) {
		int maxDepth = params.FindOneInt("maxdepth", 5);
		LightStrategy strategy;
		std::string st = params.FindOneString("strategy", "all");
		if (st == "one")
			strategy = LightStrategy::UniformSampleOne;
		else if (st == "all")
			strategy = LightStrategy::UniformSampleAll;
		else {
			Warning(
				"Strategy \"%s\" for direct lighting unknown. "
				"Using \"all\".",
				st.c_str());
			strategy = LightStrategy::UniformSampleAll;
		}
		int np;
		const int *pb = params.FindInt("pixelbounds", &np);
		Bounds2i pixelBounds = camera->film->GetSampleBounds();
		if (pb) {
			if (np != 4)
				Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
					np);
			else {
				pixelBounds = Intersect(pixelBounds,
					Bounds2i{ { pb[0], pb[2] },{ pb[1], pb[3] } });
				if (pixelBounds.Area() == 0)
					Error("Degenerate \"pixelbounds\" specified.");
			}
		}
		return new DirectLightingIntegrator(strategy, maxDepth, camera, sampler,
			pixelBounds);
	}

}  // namespace pbrt
