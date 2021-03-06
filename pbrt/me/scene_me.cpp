#include <stdlib.h>

#include "scene.h"
#include "interaction.h"
#include "sampling.h"
#include "hair.h"


namespace pbrt {


	inline bool ValidCheck(Spectrum& f, std::string info) {
		if (std::isinf(f.y()) || std::isnan(f.y())) {
			return false;
		}
		if (f.y() < 0.0) {
			return false;
		}
		if (f.y() > 1.5 ) {
			std::cout << info << std::endl <<  f.y() << ' ' << f << std::endl;
			return false;
		}
		return true;
	}

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
#pragma omp parallel for schedule(static,1)
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
//#define UNIFORM_SAMPLE
#define SAMPLE_LIGHT
#ifdef UNIFORM_SAMPLE
				for (int i = 0; i < volume->nSHSample; ++i) {
					Vector3f wi = volume->shSample[i].w;
					Spectrum Tr(1.0);
					Float pdf;
					Spectrum Li = light->Li(it, wi, &pdf, &vis);
					//Spectrum Li = SphericalLightFunc(wi);
					pdf = UniformSpherePdf();
					if (!Li.IsBlack()) {
						++validSample;
						auto tr = vis.Tr(*this);
						Li *= tr;
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

	void Scene::InitHairBSDF(bool isLight) {
		bool lightHair = true;
		Float h = 0.5;
		Float e = 1.55;
		Float bm = 0.3;
		Float bn = 0.3;
		Float a = 2;
		Float s[3] = { 0.1257, 0.2091, 0.411 };
		Spectrum sig_a = Spectrum::FromRGB(s);
		if (!lightHair) {
			h = 0.5;
			e = 1.55;
			bm = 0.25;
			bn = 0.3;
			a = 2;
			Float s[3] = { 0.5447, 0.9601, 1.781 };
			sig_a = Spectrum::FromRGB(s);
		}
		hairBSDF = std::make_shared<HairBSDF>(h, e, sig_a, bm, bn, a);
	}

	void Scene::VolumeBSDFMatrix(int oSamples, int iSamples) {
		//HairBSDF bsdf()
		volume->ComputeBSDFMatrix(*hairBSDF, oSamples, iSamples);
	}




#define UPDATE_PAPER
	//calculate sh parameters. 
	//Args:
	//    startPoint: Segment startPoint
	//    endPoint: Segment endPoint
	//	  power: the light power
	//Returns:
	//    int, endPoint idx in volume
	//breif:
	//     accept two point, respect a light segment, the two point can in one voxel or two, update the parameters,
	//     if in two, each voxel should calculate the contribution and update
	//     multiple the attention when light leave the voxel.  ps: only in leave voxel.
	int Scene::UpdateILFromTwoPoint(const Point3f& startPoint, const Point3f& endPoint, Spectrum* power, int& moveLen) {
		int nextIdx = volume->GetIdxFromPoint(endPoint.x, endPoint.y, endPoint.z);
		int curIdx = volume->GetIdxFromPoint(startPoint.x, startPoint.y, startPoint.z);
		Vector3f d = Normalize(endPoint - startPoint);

		//in one voxel
		if (nextIdx == curIdx) {
			Float length = (endPoint - startPoint).Length();
			volume->voxel[curIdx].lightLength += length;
			volume->UpdateSHFromLightSegment(d, length, curIdx, *power);
		} else {
			Ray photonRay(startPoint, d);
			Float t1, t2;
			Float curVoxelLength = 0.f, nextVoxelLenth = 0.f;
			if (volume->voxel[curIdx].bound.IntersectP(photonRay, &t1, &t2)) {
				Point3f boundIntersectLoc = photonRay.o + photonRay.d * t2;
				curVoxelLength = (boundIntersectLoc - startPoint).Length();
				nextVoxelLenth = (endPoint - boundIntersectLoc).Length();
			}
			volume->UpdateSHFromLightSegment(d, curVoxelLength, curIdx, *power);
			volume->voxel[curIdx].lightLength += curVoxelLength;
			
			//这里有点问题， 先用总距离乘垂直衰减， 不行再尝试其他的.
			//判断当前体素内有没有头发
			if (volume->voxel[curIdx].sigma != 0) {
#ifndef UPDATE_PAPER
				*power *= std::exp(-volume->voxel[curIdx].lightLength * volume->voxel[curIdx].sigma);
#else
				*power *= std::exp(-volume->voxel[curIdx].lightLength * volume->voxel[curIdx].sigma * std::sin(std::acos(Dot(volume->voxel[curIdx].avgDirection, d))));
				//*power *= std::exp(curVoxelLength * volume->voxel[curIdx].sigma * std::sin(std::acos(Dot(volume->voxel[curIdx].avgDirection, d))));
				if (!ValidCheck(*power, "UpdateILFromTwoPoint   attention power!")) {
					return -1;
				}
#endif
				moveLen = 0;
				//attention according sigma(theta)
			} 
			else {
				if (++moveLen == 15) {
					//std::cout << moveLen << " moveEnd" <<std::endl;
					return -1;
				}
			}
			if (nextIdx != -1) {
				volume->UpdateSHFromLightSegment(d, nextVoxelLenth, nextIdx, *power);
				volume->voxel[nextIdx].lightLength += nextVoxelLenth;
			}
		}

		return nextIdx;
	}


//#define ScatterInfo
#define sampleHairDirection


	void Scene::VolumeIndirectLight(int sampleNum) {
		Float marchSize = volume->GetMaxDimDelta() * 3.2;
		std::vector<Point2f> u1(sampleNum), u2(sampleNum), bsdfU(sampleNum);
		RNG rng;
		MemoryArena arena;
		for (int i = 0; i < sampleNum; ++i) {
			u1[i] = { rng.UniformFloat(), rng.UniformFloat() };
			u2[i] = { rng.UniformFloat(), rng.UniformFloat() };
			bsdfU[i] = { rng.UniformFloat(), rng.UniformFloat() };
		}

		int intersectNum = 0;
		int scatterNum = 0;
		//for each light
		for (size_t j = 0; j < lights.size(); ++j) {
			const std::shared_ptr<Light>& light = lights[j];
			//sample light ray
			const int nDim = static_cast<int>(std::floor(std::sqrt(sampleNum)));
			//StratifiedSample2D(u1.data(), nDim, nDim, rng);
			//StratifiedSample2D(u2.data(), nDim, nDim, rng);
			//StratifiedSample2D(bsdfU.data(), nDim, nDim, rng);

			Float scale = 1.0 / static_cast<Float>(sampleNum);
			bool inverseN = false;
#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic,1) private(arena) //234942  // 439885
#endif
			for (int i = 0; i < sampleNum; ++i) {
				Ray photonRay;
				Normal3f nLight;
				Float pdfPos, pdfDir;
				Spectrum Le = light->Sample_Le(u1[i], u2[i], inverseN ? -1.0 : Infinity_, &photonRay, &nLight, &pdfPos, &pdfDir);
				if (pdfPos == 0 || pdfDir == 0 || Le.IsBlack()) {
					continue;
				}
				//Spectrum beta = (AbsDot(nLight, photonRay.d) * Le) / (pdfPos * pdfDir);
				Spectrum lightPower = light->Power();
				if (lightPower.IsBlack()) {
					continue;
				}

// 				Float envScale = 512 * 512 * scale;
// 				Le *= envScale;
				//the photon energy
				//Le *= (beta * scale);
				if (light->flags == 8) {
					Le = Le * lightPower * scale * 4.0 * 30.0 / 7.0;
				} else {
					Le = lightPower * scale;
				}

				if (!ValidCheck(Le, "sample Le and scale")) {
					continue;
				}

				//ray scene intersect test
				SurfaceInteraction isect;
				if (!Intersect(photonRay, &isect)) {
					continue;
				}


				isect.ComputeScatteringFunctions(photonRay, arena, true, TransportMode::Importance);
				if (!isect.bsdf) {
					photonRay = isect.SpawnRay(photonRay.d);
					continue;
				}



				++intersectNum;


				const BSDF &photonBSDF = *isect.bsdf;
				Vector3f wi, wo = -photonRay.d;
				Float pdf;
				BxDFType flags;
				Spectrum fr = photonBSDF.Sample_f(wo, &wi, bsdfU[i], &pdf, BSDF_ALL, &flags);
				if (pdf <= 0) {
					continue;
				}
				fr = fr * AbsDot(wi, isect.shading.n) / pdf;
				if (!ValidCheck(fr, "first intersect hair bsdf")) {
					continue;
				}

				Le *= fr;

				if (!ValidCheck(Le, "first intersect hair * Le")) {
					continue;
				}

				//indirect lighting
				//first scatter in voxel
			
				int curIdx = volume->GetIdxFromPoint(isect.p.x, isect.p.y, isect.p.z);
				Voxel& v = volume->voxel[curIdx];

				Point3f curPoint = isect.p;
				Point3f nextPoint = isect.p + marchSize * wi;
				int moveLen = 1;

				//update current voxel
				//curIdx = UpdateILFromTwoPoint(curPoint, nextPoint, &Le, moveLen);

				// 第一次到第二次之间的能量要不要记录。
				curPoint = nextPoint;
				curIdx = volume->GetIdxFromPoint(nextPoint.x, nextPoint.y, nextPoint.z);



				//std::cout << std::endl;
				//light traverse, everytime move unit marchsize
				while (!Le.IsBlack() && curIdx > 0) {
					//scattering event according to voxel parameters
					Voxel& curV = volume->voxel[curIdx];
					Float attentionSintheta = std::sqrt(1.f - std::pow(Dot(curV.avgDirection, wi), 2));
					Float p = 1 - std::exp(-curV.sigma * attentionSintheta); //* marchSize);
					//std::cout << p << std::endl;

// 					int x, y, z;
// 					volume->GetXYZ(curIdx, &x, &y, &z);
// 					std::cout << curIdx << ' ' << x << ' ' << y << ' ' << z << std::endl;

					bool isScatter = rng.UniformFloat() < p;
					if (isScatter && curV.sigma != 0) {
#pragma omp atomic
						++scatterNum;
						//update ray direction
						//这里少了一个选择n个过程， 因为采样方向都是假定n是z轴正方向的。
						//所以要有一个世界-》局部， 局部-》世界的过程， 就想bsdf里面做的那样。

						//根据论文， 先根据头发方向分布的标准差采样一个方向v， 然后推算出法向n。
						//先用平均头发方向代替以下

						//先用随机方向代替以下
						wo = -wi;
						Point2f scatterU = { rng.UniformFloat(), rng.UniformFloat() };
						fr = ScatterEvent(wo, &wi, curV.directionV, curV.avgDirection, scatterU, &pdf, &flags);


						Le *= fr;
						

#ifdef ScatterInfo
						std::cout << "scatter! scatter p: " << p << std::endl;
						Float oPhi = Degrees(SphericalPhi(wo));
						Float oTheta = Degrees(SphericalTheta(wo));
						Float iPhi = Degrees(SphericalPhi(wi));
						Float iTheta = Degrees(SphericalTheta(wi));
						std::cout << "wo phi: " << oPhi << " wo theta: " << oTheta << std::endl;
						std::cout << "wi phi: " << iPhi << " wi theta: " << iTheta << std::endl;
						std::cout << "attention: " << (fr * AbsCosTheta(wi) / pdf) << std::endl;
#endif
					}

					nextPoint = curPoint + marchSize * wi;
					curIdx = UpdateILFromTwoPoint(curPoint, nextPoint, &Le, moveLen);
					curPoint = nextPoint;
				}
				//std::cout << i << std::endl;
			}
		}
		std::cout << "intersectNum: " << intersectNum << " scatterNum: " << scatterNum << std::endl;
	}


	Spectrum Scene::ScatterEvent(const Vector3f& wo, Vector3f* wi, Float v,
		                     const Vector3f& hairDirection, const Point2f& u, Float* pdf, BxDFType* type) const {
#ifdef sampleHairDirection
		Vector3f ss, ts, ns;
		ss = hairDirection;

		//sampling a hair direction according to v.
#ifdef UPDATE_PAPER
		CoordinateSystem(ss, &ts, &ns);
		Float theta, phi;
		theta = std::acos(1.0 - u[0] * v);
		phi = u[1] * Pi * 2;
		Vector3f tV = { std::sin(theta) * std::cos(phi),  std::sin(theta) * std::sin(phi), std::cos(theta) };
		ss = tV.x * ts + tV.y * ns + tV.z * ss;
#endif
		

		//ts = Normalize(Cross(hairDirection, wo));

		//ss is hair direction, according it generate N
		ts = Normalize(Cross(ss, wo));
		ns = Normalize(Cross(ss, ts));
		//world to loacal
		Vector3f localWo = Vector3f(Dot(wo, ss), Dot(wo, ts), Dot(wo, ns));
#else
		Vector3f localWo = wo;
#endif
		Spectrum fr = hairBSDF->Sample_f(localWo, wi, u, pdf, type);
		if (*pdf <= 0) {
			fr = Spectrum(0.0);
		}
#ifdef sampleHairDirection
		//local to world
		*wi = Normalize(Vector3f(
			ss.x * wi->x + ts.x * wi->y + ns.x * wi->z,
			ss.y * wi->x + ts.y * wi->y + ns.y * wi->z,
			ss.z * wi->x + ts.z * wi->y + ns.z * wi->z
		));
#endif

		fr = fr * AbsDot(*wi, ns) / *pdf;
		if (!ValidCheck(fr, "scatter event fr")) {
			fr = Spectrum(0.0);
		}


		return fr;
	}

}