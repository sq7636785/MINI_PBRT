#include <stdlib.h>

#include "scene.h"
#include "interaction.h"
#include "sampling.h"
#include "hair.h"


namespace pbrt {
	void Scene::VolumeIrrandiance() {

		//对所有光源采样， 计算辐照度
		for (size_t j = 0; j < lights.size(); ++j) {
			const std::shared_ptr<Light>& light = lights[j];
			
			for (auto &v : volume->voxel) {
				VisibilityTester vis;
				Interaction it;
				Vector3f wi;
				Float pdf;

				//注意， 这里没有加入随机数
				Point2f uLight;

				it.p = (v.bound.pMin + v.bound.pMax) / 2;
				Spectrum Tr(1.0);
				Spectrum Li = light->Sample_Li(it, uLight, &wi, &pdf, &vis); 

				if (pdf > 0.f && !Li.IsBlack()) {
					//	std::cout << "before Tr: " << Li << std::endl; 
					auto tr = vis.Tr(*this);
					//std::cout << tr << std::endl;
					Li *= tr;
					//	std::cout << "after Tr: " << Li << std::endl;
						//irrandiance
//					Float cosTheta = 1.0;

					//平均方向不为0. 说明有头发， 那么取平均头发方向的切向为法向量.
// 					if (v.avgDirection.LengthSquared() != 0) {
// 						Float sinTheta = AbsDot(v.avgDirection, wi);
// 						cosTheta = std::sqrt(1.0 - sinTheta * sinTheta);
// 					}
				} else {
					Li = Spectrum(0.f);
				}
		//		std::cout << "after Dot: " << Li << std::endl << std::endl;
				Float tmpRGB[3];
				Li.ToRGB(tmpRGB);
				v.rgb[0] += tmpRGB[0];
				v.rgb[1] += tmpRGB[1];
				v.rgb[2] += tmpRGB[2];
			}
		}
	}


	

	void Scene::VolumeSHRadiance() {
		//for each direction

		//for each voxel
		for (int voxelId = 0; voxelId < volume->voxel.size(); ++voxelId) {
			Voxel& v = volume->voxel[voxelId];
			VisibilityTester vis;
			Interaction it;
			it.p = (v.bound.pMin + v.bound.pMax) / 2;
			//for each light
			int validSample = 0;
			for (size_t j = 0; j < lights.size(); ++j) {
				const std::shared_ptr<Light>& light = lights[j];
				//for each sample
#define UNIFORM_SAMPLE
//#define SAMPLE_LIGHT
#ifdef UNIFORM_SAMPLE
				for (int i = 0; i < volume->nSHSample; ++i) {
					Vector3f wi = volume->shSample[i].w;
					Spectrum Tr(1.0);
					Float pdf;
					//Spectrum Li = light->Li(it, wi, &pdf, &vis);
					Spectrum Li = SphericalLightFunc(wi);
					pdf = UniformHemispherePdf();
					if (!Li.IsBlack()) {
						++validSample;
					//	 						auto tr = vis.Tr(*this);
					//	 						Li *= tr;
												//std::cout << Li.y() << std::endl;
						for (int k = 0; k < SHTerms(volume->shL); ++k) {
							v.shC[k] += volume->shSample[i].y[k] * Li / (pdf * volume->nSHSample);
						}
					}
				}
			}
#endif

#ifdef SAMPLE_LIGHT
			int nDim = 25;
			int nSampleLight = nDim * nDim;
			RNG rng;
			std::vector<Point2f> u(nSampleLight);
			StratifiedSample2D(u.data(), nDim, nDim, rng);
			for (int i = 0; i < nSampleLight; ++i) {
				VisibilityTester vis;
				Interaction it;
				Vector3f wi;
				Float pdf;

				//注意， 这里没有加入随机数

				it.p = (v.bound.pMin + v.bound.pMax) / 2;
				Spectrum Li = light->Sample_Li(it, u[i], &wi, &pdf, &vis);

				if (pdf > 0.f && !Li.IsBlack()) {
					//pdf = UniformHemispherePdf();
					auto tr = vis.Tr(*this);
					Li *= tr;
					std::vector<Float> ylm(SHTerms(volume->shL));
					SHEvaluate(wi, volume->shL, ylm.data());
					for (int k = 0; k < SHTerms(volume->shL); ++k) {
						v.shC[k] += ylm[k] * Li / (pdf * nSampleLight);
					}
				}
			}
		}
#endif
			

				//std::cout << validSample << std::endl;
				//Float scale = 1.f / static_cast<Float>(validSample);
	// 			for (auto &c : v.shC) {
	// 				std::cout << c << ' ' << std::endl;
	// 				//std::cout << scale << ' ' << c << std::endl;
	// 			}
			//std::cout << voxelId << std::endl;
		}
	}

	void Scene::VolumeBSDFMatrix() {
		//HairBSDF bsdf()
		Float h = 0.5;
		Float e = 1.55;
		Float bm = 0.25;
		Float bn = 0.3;
		Float a = 2;
		Float s[3] = { 0.5447, 0.9601, 1.781 };
		Spectrum sig_a = Spectrum::FromRGB(s);
		HairBSDF bsdf(h, e, sig_a, bm, bn, a);
		volume->ComputeBSDFMatrix(bsdf);
	}
}