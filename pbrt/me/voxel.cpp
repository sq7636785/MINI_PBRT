#include "voxel.h"
#include "rng.h"
#include "sampling.h"
#include <fstream>
#include <omp.h>
#include "spherical_harmonics.h"
namespace pbrt {
	inline bool ValidCheck2(Spectrum& f) {
		if (std::isinf(f.y()) || std::isnan(f.y())) {
			return false;
		}
		if (f.y() < 0.0) {
			return false;
		}
		if (f.y() > 1.0) {
			std::cout << f.y() << "render" << std::endl;
		}
		return true;
	}


//#define UNIFORM_REFLECT
#define BSDF_SAMPLE
	//#define BRUTEFORCE

	Volume::Volume(const Bounds3f& bound, Float partitionNum /* = 100.0 */, int shL /* = 4 */, int nSHSample /* = 10000 */)
		: worldBound(bound), partitionNum(partitionNum), shL(shL), nSHSample(nSHSample) {
		Vector3f axisLen = worldBound.Diagonal();
		Float invPartitionNum = 1.0 / partitionNum;
		xDelta = axisLen.x * invPartitionNum;
		yDelta = axisLen.y * invPartitionNum;
		zDelta = axisLen.z * invPartitionNum;
		vDelta = xDelta * yDelta * zDelta;
	}

	bool Volume::ConstructVolume() {

		Point3f basePoint = worldBound.pMin;

//		std::cout << xDelta << ' ' << yDelta << ' ' << zDelta << std::endl;

		//construct voxel
		for (int z = 1; z <= partitionNum; ++z) {
			for (int y = 1; y <= partitionNum; ++y) {
				for (int x = 1; x <= partitionNum; ++x) {
					Voxel tmp;
					Point3f curPointMin = basePoint + Point3f((x - 1) * xDelta, (y - 1) * yDelta, (z - 1) * zDelta);
					Point3f curPointMax = basePoint + Point3f(x * xDelta, y * yDelta, z * zDelta);
					tmp.bound = Bounds3f(curPointMin, curPointMax);
					tmp.shC.assign(SHTerms(shL), Spectrum(0.f));
					
					voxel.push_back(tmp);
				}
			}
		}
		bool constructVoxel =  (voxel[0].bound.pMin - worldBound.pMin).Length() < 0.001 && 
							   (voxel[voxel.size() - 1].bound.pMax - worldBound.pMax).Length() < 0.001;

#ifdef UNIFORM_REFLECT
		shSample.assign(nSHSample, SHSample());
		RNG rng;
		std::vector<Point2f> u(nSHSample);
		int sqrtNSample = static_cast<int>(std::sqrt(nSHSample));
		StratifiedSample2D(u.data(), sqrtNSample, sqrtNSample, rng);

		for (int i = 0; i < nSHSample; ++i) {
			shSample[i].w = UniformSampleSphere(u[i]);
			shSample[i].y.assign(SHTerms(shL), 0.f);
			SHEvaluate(shSample[i].w, shL, shSample[i].y.data());
		}
#endif
		bsdfMatrix.assign(SHTerms(shL) * SHTerms(shL), Spectrum(0.f));
		return constructVoxel;
	}

	
	bool Volume::CalculateVoxel(const std::vector<std::shared_ptr<Primitive>>& curves, bool mode){
		int validVoxel = 0;
		
		Float validLength = std::min(xDelta, std::min(yDelta, zDelta)) / 2.0;
		Float invD = 1.0 / validLength;
		


		std::vector<CalUtil> innerData(voxel.size());

#pragma region __CurveByVoxel
		//对每一根头发去找他能覆盖到的体素的范围，再来计算。
		if (mode) {
			int curvesSize = curves.size();
#pragma omp parallel for schedule(static,1)
			for (int curveId = 0; curveId < curvesSize; ++curveId) {
				auto boundSet = GetVoxelSet(*curves[curveId]);
				for (auto voxelId : boundSet) {

					Float distance;
					Vector3f direction;
					Float width;

					if (curves[curveId]->GetShape()->DistanceToPoint(voxel[voxelId].bound, &width, &distance, &direction)) {
						if (distance < validLength) {
#pragma omp critical
							{
								innerData[voxelId].distances.push_back(distance);
								innerData[voxelId].directions.push_back(direction);
								innerData[voxelId].width.push_back(width);
							}
						}
					}
				}
				//std::cout << curveId << std::endl;
			}


			int voxelSize = voxel.size();
//#pragma omp parallel for schedule(static,1)
			for (int idx = 0; idx < voxelSize; ++idx) {
				Float Num = static_cast<Float>(innerData[idx].distances.size());
				if (Num != 0.0) {
					Float avgDiameter = 0.0;
					Vector3f avgDirection;
					for (size_t i = 0; i < innerData[idx].distances.size(); ++i) {
						Float weight = 1 - innerData[idx].distances[i] * invD;
						avgDirection += innerData[idx].directions[i] * weight;
						avgDiameter += innerData[idx].width[i];
					}
					avgDiameter /= Num;

					voxel[idx].sigma = avgDiameter * Num * InvPi * invD * invD;
					voxel[idx].avgDirection = Normalize(avgDirection);
					++validVoxel;
					
					//calculate v
					Float avgDot = 0.f;
					for (size_t i = 0; i < innerData[idx].directions.size(); ++i) {
						avgDot += Dot(Normalize(innerData[idx].directions[i]), voxel[idx].avgDirection);
					}
					avgDot /= Num;
					Float directionV = 0.f;
					for (size_t i = 0; i < innerData[idx].directions.size(); ++i) {
						directionV += std::pow(Dot(Normalize(innerData[idx].directions[i]), voxel[idx].avgDirection) - avgDot, 2);
					}
					voxel[idx].directionV = std::sqrt(directionV / Num);
					voxel[idx].hairNum = Num;
					//std::cout << voxel[idx].directionV << std::endl;
				}
			}
		}
#pragma endregion

		
#pragma region __VoxelByCurve
		if (!mode) {
			for (size_t idx = 0; idx < voxel.size(); ++idx) {
				Vector3f diagnal = voxel[idx].bound.Diagonal();
				for (size_t curveId = 0; curveId < curves.size(); ++curveId) {
					Float distance;
					Vector3f direction;
					Float width;
					if (Overlaps(curves[curveId]->WorldBound(), voxel[idx].bound)) {
						if (curves[curveId]->GetShape()->DistanceToPoint(voxel[idx].bound, &width, &distance, &direction)) {
							if (distance < validLength) {
								innerData[idx].distances.push_back(distance);
								innerData[idx].directions.push_back(direction);
								innerData[idx].width.push_back(width);
							}
						}
					}
				}
				Float Num = static_cast<Float>(innerData[idx].distances.size());
				if (Num != 0.0) {
					Float avgDiameter = 0.0;
					Vector3f avgDirection;
					for (size_t i = 0; i < innerData[idx].distances.size(); ++i) {
						Float weight = 1 - innerData[idx].distances[i] * invD;
						avgDirection += innerData[idx].directions[i] * weight;
						avgDiameter += innerData[idx].width[i];
					}
					avgDiameter /= Num;
					voxel[idx].sigma = 2.0 * avgDiameter * Num * InvPi * invD * invD;
					voxel[idx].avgDirection = Normalize(avgDirection);
					++validVoxel;
				}
				//calculate v
				Float avgDot = 0.f;
				for (size_t i = 0; i < innerData[idx].directions.size(); ++i) {
					avgDot += Dot(Normalize(innerData[idx].directions[i]), voxel[idx].avgDirection);
				}
				avgDot /= Num;
				Float directionV = 0.f;
				for (size_t i = 0; i < innerData[idx].directions.size(); ++i) {
					directionV += std::pow(Dot(Normalize(innerData[idx].directions[i]), voxel[idx].avgDirection) - avgDot, 2);
				}
				voxel[idx].directionV = std::sqrt(directionV / Num);
				std::cout << voxel[idx].directionV << std::endl;
			}
		}
#pragma endregion


		return validVoxel > 0;
	}

	Spectrum Volume::Tr(const Point3f& p0, const Point3f& p1) const {
		Vector3f rayD = p1 - p0;


		Float t1, t2;
		if (!worldBound.IntersectP(Ray(p0, rayD, 1.0), &t1, &t2)) {
			return Spectrum(1.0);
		}

		Point3f startPos = p0 + rayD * (t1 - 0.0001);
		Point3f endPos = p0 + rayD * (t2 - 0.0001);
		rayD = Normalize(rayD);
		Ray ray(startPos, rayD);

		Spectrum tr(1.0);
		
		int xMove = rayD.x > 0.0 ? 1 : -1;
		int yMove = rayD.y > 0.0 ? 1 : -1;
		int zMove = rayD.z > 0.0 ? 1 : -1;

		//start idx
		Vector3f deltaPosStart = startPos - worldBound.pMin;
		int xIdx = static_cast<int>(deltaPosStart.x / xDelta);
		int yIdx = static_cast<int>(deltaPosStart.y / yDelta);
		int zIdx = static_cast<int>(deltaPosStart.z / zDelta);
		
		//end idx
		Vector3f deltaPosEnd = endPos - worldBound.pMin;
		int endX = static_cast<int>(deltaPosEnd.x / xDelta);
		int endY = static_cast<int>(deltaPosEnd.y / yDelta);
		int endZ = static_cast<int>(deltaPosEnd.z / zDelta);

		int curIdx = GetIdx(xIdx, yIdx, zIdx);
		int endIdx = GetIdx(endX, endY, endZ);

		while (curIdx != endIdx) {
			
			if (voxel[curIdx].bound.IntersectP(ray, &t1, &t2)) {
				Point3f startLoc = ray.o + ray.d * t1;
				Point3f endLoc = ray.o + ray.d * t2;
				Float length = (endLoc - startLoc).Length();
				//tr *= (1.0 - voxel[curIdx].sigma * length);
				tr *= std::exp(-voxel[curIdx].sigma * length);
				ray.o = endLoc + ray.d * 0.0001f;
				xIdx += xMove;
				yIdx += yMove;
				zIdx += zMove;
				curIdx = GetIdx(xIdx, yIdx, zIdx);
			} else {
				break;
			}
		}

		return tr;
	}


	void Volume::SaveData(std::string fileName) const {
		std::ofstream out(fileName);
		if (!out) {
			std::cout << "can open file" << std::endl;
			return;
		}

		for (size_t i = 0; i < voxel.size(); ++i) {
			int x, y, z;
			GetXYZ(i, &x, &y, &z);
			out << x << ' ' << y << ' ' << z << ' ';
			out << voxel[i].sigma << ' ' << voxel[i].avgDirection.x << ' ' << voxel[i].avgDirection.y << ' ' << voxel[i].avgDirection.z << ' ';
			for (auto &c : voxel[i].shC) {
				out << c.y() << ' ';
			}
			out << std::endl;
//			out << voxel[i].rgb[0] << ' ' << voxel[i].rgb[1] << ' ' << voxel[i].rgb[2] << std::endl;
		}
		out.close();
	}

	void Volume::LoadData(std::string fileName) {
		std::ifstream in(fileName);
		if (!in) {
			std::cout << "can open file" << std::endl;
			return;
		}
		for (size_t i = 0; i < voxel.size(); ++i) {
			int x, y, z;
			Float sigma, avgDX, avgDY, avgDZ;
			//in >> voxel[i].sigma >> voxel[i].avgDirection.x >> voxel[i].avgDirection.y >> voxel[i].avgDirection.z;
			in >> x >> y >> z;
			if (i != GetIdx(x, y, z)) {
				std::cout << "idx error" << std::endl;
			}

			in >> sigma >> avgDX >> avgDY >> avgDZ;
			voxel[i].sigma = sigma;
			voxel[i].avgDirection.x = avgDX;
			voxel[i].avgDirection.y = avgDY;
			voxel[i].avgDirection.z = avgDZ;
			//std::cout << x1 << ' ' << x2 << ' ' << x3 << ' ' << x4 << std::endl;
			voxel[i].directionV = 0.0;
		}
		in.close();
	}

	std::vector<int> Volume::GetVoxelSet(const Primitive& curve) const {
		Bounds3f curveBound = curve.WorldBound();
		std::vector<int> result;
		int total = partitionNum * partitionNum * partitionNum;
		for (Float z = curveBound.pMin.z; z < curveBound.pMax.z; z += zDelta) {
			for (Float y = curveBound.pMin.y; y < curveBound.pMax.y; y += yDelta) {
				for (Float x = curveBound.pMin.x; x < curveBound.pMax.x; x += xDelta) {
					int idx = GetIdxFromPoint(x, y, z);
					if (idx > 0 && idx < total) {
						result.push_back(idx);
					}
				}
			}
		}
		return result;
	}



	void Volume::ComputeBSDFMatrix(const BxDF& bsdf, int oSamples, int iSamples /*= 50*50*/) {

		// Precompute directions $\w{}$ and SH values for directions
		omp_lock_t ompLock;
		omp_init_lock(&ompLock);


#ifdef UNIFORM_REFLECT
// 		std::vector<Float> Ylm(SHTerms(shL) * nSamples);
// 		std::vector<Vector3f> w(nSamples);
// 		std::vector<Point2f> u(nSamples);
// 		RNG rng;
// 		int nDim = static_cast<int>(std::sqrt(nSamples));
// 		StratifiedSample2D(u.data(), nDim, nDim, rng);
// 
// 		for (int i = 0; i < nSamples; ++i) {
// 			w[i] = UniformSampleHemisphere(u[i]);
// 			SHEvaluate(w[i], shL, &Ylm[SHTerms(shL)*i]);
// 		}

		oSamples = nSHSample;
		iSamples = nSHSample;
//#pragma omp parallel for schedule(static,1)
		for (int osamp = 0; osamp < oSamples; ++osamp) {
			const Vector3f &wo = shSample[osamp].w;
			for (int isamp = 0; isamp < iSamples; ++isamp) {
				const Vector3f &wi = shSample[isamp].w;
				// Update BSDF matrix elements for sampled directions
				Spectrum f = bsdf.f(wo, wi);
				if (!f.IsBlack()) {
					Float pdf = UniformSpherePdf() * UniformSpherePdf();
					//n -> z+
					f *= (AbsCosTheta(wi)) / (pdf * oSamples * iSamples);
					for (int i = 0; i < SHTerms(shL); ++i) {
						for (int j = 0; j < SHTerms(shL); ++j) {
							//omp_set_lock(&ompLock);
							bsdfMatrix[i*SHTerms(shL) + j] += f * shSample[isamp].y[j] *
								shSample[osamp].y[i];
							//omp_unset_lock(&ompLock);
						}
					}
				}
			}
			std::cout << osamp << std::endl;
		}
#endif
#ifdef BSDF_SAMPLE
		//int oSamples = 10000;//40000;
		//int iSamples = 64;//256;

		std::vector<Float> Ylm(SHTerms(shL) * oSamples);
		std::vector<Vector3f> w(oSamples);
		std::vector<Point2f> u(oSamples);
		RNG rng;
		int nDim = static_cast<int>(std::sqrt(oSamples));
		StratifiedSample2D(u.data(), nDim, nDim, rng);

		for (int i = 0; i < oSamples; ++i) {
			w[i] = UniformSampleSphere(u[i]);
			SHEvaluate(w[i], shL, &Ylm[SHTerms(shL)*i]);
		}

		
		
		// Compute double spherical integral for BSDF matrix
//#pragma omp parallel for schedule(dynamic, 1)
		for (int osamp = 0; osamp < oSamples; ++osamp) {
			const Vector3f &wo = w[osamp];
			
			std::vector<Point2f> uWi(iSamples);
			nDim = static_cast<int>(std::sqrt(iSamples));
			StratifiedSample2D(uWi.data(), nDim, nDim, rng);

			for (int isamp = 0; isamp < iSamples; ++isamp) {
				//sample bsdf
				Vector3f wi;
				Float wiPdf;
				// Update BSDF matrix elements for sampled directions
				Spectrum f = bsdf.Sample_f(wo, &wi, uWi[isamp], &wiPdf);
				std::vector<Float> ylmWi(SHTerms(shL));
				SHEvaluate(wi, shL, ylmWi.data());

				if (!f.IsBlack()) {
					Float pdf = UniformSpherePdf() * wiPdf;
					//n -> z+
					f *= (AbsCosTheta(wi)) / (pdf * oSamples * iSamples);
					for (int i = 0; i < SHTerms(shL); ++i) {
						for (int j = 0; j < SHTerms(shL); ++j) {
//#pragma omp critical
						//	omp_set_lock(&ompLock);
							bsdfMatrix[i*SHTerms(shL) + j] += f * ylmWi[j] * Ylm[osamp*SHTerms(shL) + i];
						//	omp_unset_lock(&ompLock);
						}
					}
				}
			}
			//std::cout << osamp << std::endl;
		}
		
#endif
#ifdef BRUTEFORCE
	//	int oSamples = 160;
	//	int iSamples = 160;
		std::vector<Float> Ylm(SHTerms(shL) * oSamples);
		std::vector<Vector3f> wo(oSamples);
		std::vector<Point2f> u(oSamples);
		RNG rng;
		int oDim = static_cast<int>(std::sqrt(oSamples));
		StratifiedSample2D(u.data(), oDim, oDim, rng);

		for (int i = 0; i < oSamples; ++i) {
			wo[i] = UniformSampleSphere(u[i]);
			SHEvaluate(wo[i], shL, &Ylm[SHTerms(shL)*i]);
		}

#pragma omp parallel for schedule(static,1)
		for (int i = 0; i < SHTerms(shL); ++i) {
			for (int j = 0; j < SHTerms(shL); ++j) {
				Spectrum oEstimat(0.f);
				for (int oD = 0; oD < oSamples; ++oD) {

					Spectrum iEstimat(0.f);
					std::vector<Point2f> uWi(iSamples);
					int iDim = static_cast<int>(std::sqrt(iSamples));
					StratifiedSample2D(uWi.data(), iDim, iDim, rng);
					for (int iD = 0; iD < iSamples; ++iD) {
						Vector3f wi;
						Float wiPdf;
						// Update BSDF matrix elements for sampled directions
						Spectrum f = bsdf.Sample_f(wo[oD], &wi, uWi[iD], &wiPdf);
						f *= AbsCosTheta(wi);
						std::vector<Float> ylmWi(SHTerms(shL));
						SHEvaluate(wi, shL, ylmWi.data());
						iEstimat += f * Ylm[SHTerms(shL) * oD + i] * ylmWi[j] / wiPdf;
					}
					iEstimat /= static_cast<Float>(iSamples);
					oEstimat += (iEstimat / UniformSpherePdf());
				}
				oEstimat /= static_cast<Float>(oSamples);
				bsdfMatrix[i * SHTerms(shL) + j] = oEstimat;
			}
		}
#endif
		omp_destroy_lock(&ompLock);

		std::ofstream outBsdfMatrix;
		std::string folderPath = "./bsdf_matrix/";
		folderPath += std::to_string(shL) + ".txt";
		outBsdfMatrix.open(folderPath);
		for (int i = 0; i < SHTerms(shL); ++i) {
			for (int j = 0; j < SHTerms(shL); ++j) {
				Float tmpRGB[3];
				bsdfMatrix[i * SHTerms(shL) + j].ToRGB(tmpRGB);
				outBsdfMatrix << tmpRGB[0] << ' ' << tmpRGB[1] << ' ' << tmpRGB[2] << std::endl;
			}
		}
		outBsdfMatrix.close();
	}


	void Volume::LoadBsdfMatrix(std::string fileName) {
		std::ifstream in(fileName);
		if (!in) {
			std::cout << "can open file" << std::endl;
			return;
		}
		for (int i = 0; i < SHTerms(shL); ++i) {
			for (int j = 0; j < SHTerms(shL); ++j) {
				Float tmpRGB[3];
				in >> tmpRGB[0] >> tmpRGB[1] >> tmpRGB[2];
				bsdfMatrix[i * SHTerms(shL) + j] = Spectrum::FromRGB(tmpRGB);
			}
		}
		in.close();
	}



	static
		Transform Rotate(Float cosTheta, Float sinTheta, const Vector3f &axis) {
		Vector3f a = Normalize(axis);
		// 		Float sinTheta = std::sin(Radians(theta));
		// 		Float cosTheta = std::cos(Radians(theta));
		Matrix4x4 m;
		// Compute rotation of first basis vector
		m.m[0][0] = a.x * a.x + (1 - a.x * a.x) * cosTheta;
		m.m[0][1] = a.x * a.y * (1 - cosTheta) - a.z * sinTheta;
		m.m[0][2] = a.x * a.z * (1 - cosTheta) + a.y * sinTheta;
		m.m[0][3] = 0;

		// Compute rotations of second and third basis vectors
		m.m[1][0] = a.x * a.y * (1 - cosTheta) + a.z * sinTheta;
		m.m[1][1] = a.y * a.y + (1 - a.y * a.y) * cosTheta;
		m.m[1][2] = a.y * a.z * (1 - cosTheta) - a.x * sinTheta;
		m.m[1][3] = 0;

		m.m[2][0] = a.x * a.z * (1 - cosTheta) - a.y * sinTheta;
		m.m[2][1] = a.y * a.z * (1 - cosTheta) + a.x * sinTheta;
		m.m[2][2] = a.z * a.z + (1 - a.z * a.z) * cosTheta;
		m.m[2][3] = 0;
		return Transform(m, Transpose(m));
	}



	void Volume::RotateSH(const Vector3f& v1, const Vector3f& v2, Spectrum* cIn, Spectrum* cOut) const {
		Vector3f vAxis = Cross(v1, v2);
 		Float theta = std::acos(Dot(v1, v2));
		Eigen::Quaterniond r(Eigen::AngleAxisd(
			theta, Eigen::Vector3d(vAxis.x, vAxis.y, vAxis.z).normalized()));
		std::unique_ptr<sh::Rotation> rz_sh(sh::Rotation::Create(shL, r));

		std::vector<Float> c1(SHTerms(shL));
		std::vector<Float> c2(SHTerms(shL));
		std::vector<Float> c3(SHTerms(shL));
		bool isSameLight = true;
		for (int i = 0; i < SHTerms(shL); ++i) {
			Float tmpRGB[3];
			cIn[i].ToRGB(tmpRGB);
			c1[i] = tmpRGB[0];
			c2[i] = tmpRGB[1];
			c3[i] = tmpRGB[2];
			if (c1[i] != c2[i] || c1[i] != c3[i] || c2[i] != c3[i]) {
				isSameLight = false;
			}
		}
		if (isSameLight) {
			rz_sh->Apply(c1, &c1);
			c2 = c1;
			c3 = c1;
		} else {
			rz_sh->Apply(c1, &c1);
			rz_sh->Apply(c2, &c2);
			rz_sh->Apply(c3, &c3);
		}
		for (int i = 0; i < SHTerms(shL); ++i) {
			Float tmpRGB[3];
			tmpRGB[0] = c1[i];
			tmpRGB[1] = c2[i];
			tmpRGB[2] = c3[i];
			cOut[i] = Spectrum::FromRGB(tmpRGB);
		}
		

// 		Float cosTheta = Dot(v1, v2);
// 		Float sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
// 		Transform trans = Rotate(cosTheta, sinTheta, vAxis);
//  		Transform trans = Rotate(theta * 180.f / Pi, vAxis);
// 		SHRotate(cIn, cOut, trans.GetMatrix(), shL);
// 		for (int i = 0; i < SHTerms(shL); ++i) {
// 			std::cout << cOut[i].y() << ' ' << oriOut[i].y() << std::endl;
// 		}

	}

	void Volume::SHMatrixTransV(const std::vector<Spectrum>& c, Spectrum* c_out) const {
		for (int i = 0; i < SHTerms(shL); ++i) {
			c_out[i] = Spectrum(0.f);
			for (int j = 0; j < SHTerms(shL); ++j)
				c_out[i] += bsdfMatrix[SHTerms(shL) * i + j] * c[j];
		}
	}

	void Volume::UpdateSHFromLightSegment(const Vector3f& wi, Float length, int idx, const Spectrum& power) {
		Voxel& v = voxel[idx];
		std::vector<Float> Ylm(SHTerms(shL));
		SHEvaluate(wi, shL, Ylm.data());
		//std::cout << "idx: " << idx << "  power: " << power << std::endl;
#pragma omp critical
		{
			for (int i = 0; i < SHTerms(shL); ++i) {
				v.shC[i] += power * length * Ylm[i] / vDelta;
			}
		}
	}


	//这里插值有一点bug，因为有的体素内是没有光照信息的， 这样可能会有黑的地方。
	std::vector<Spectrum> Volume::InterpolateVoxelSHC(const Point3f& p) const {
		Float xOffeset = (p.x - worldBound.pMin.x) / xDelta;
		Float yOffeset = (p.y - worldBound.pMin.y) / yDelta;
		Float zOffeset = (p.z - worldBound.pMin.z) / zDelta;
		int xL = static_cast<int>(xOffeset);
		int yL = static_cast<int>(yOffeset);
		int zL = static_cast<int>(zOffeset);

		if (xL >= partitionNum - 1|| yL >= partitionNum - 1|| zL >= partitionNum - 1) {
			return voxel[GetIdx(xL, yL, zL)].shC;
		}

		xOffeset = xOffeset - xL;
		yOffeset = yOffeset - yL;
		zOffeset = zOffeset - zL;

		//x interpolate
		std::vector<Spectrum> xyLzL = LinerInterpolateSHC(voxel[GetIdx(xL, yL,     zL)].shC,     voxel[GetIdx(xL + 1, yL,     zL)].shC, xOffeset);
		std::vector<Spectrum> xyHzL = LinerInterpolateSHC(voxel[GetIdx(xL, yL + 1, zL)].shC,     voxel[GetIdx(xL + 1, yL + 1, zL)].shC, xOffeset);
		std::vector<Spectrum> xyLZH = LinerInterpolateSHC(voxel[GetIdx(xL, yL,     zL + 1)].shC, voxel[GetIdx(xL + 1, yL,     zL + 1)].shC, xOffeset);
		std::vector<Spectrum> xyHzH = LinerInterpolateSHC(voxel[GetIdx(xL, yL + 1, zL + 1)].shC, voxel[GetIdx(xL + 1, yL + 1, zL + 1)].shC, xOffeset);
		
		//y interpolate
		std::vector<Spectrum> yzL = LinerInterpolateSHC(xyLzL, xyHzL, yOffeset);
		std::vector<Spectrum> yzH = LinerInterpolateSHC(xyLZH, xyHzH, yOffeset);

		//zinterpolate
		std::vector<Spectrum> result = LinerInterpolateSHC(yzL, yzH, zOffeset);
		return result;
	}

	std::vector<Spectrum> Volume::LinerInterpolateSHC(const std::vector<Spectrum>& c1, const std::vector<Spectrum>& c2, Float offeset) const {
		std::vector<Spectrum> res(SHTerms(shL));
		for (int i = 0; i < SHTerms(shL); ++i) {
			res[i] = c1[i] * (1.0 - offeset) + c2[i] * offeset;
		}
		return res;
	}

	void Volume::BoxFilterSHC(int radius) {
		std::vector<std::vector<Spectrum>> filterSHC;
		filterSHC.assign(partitionNum * partitionNum * partitionNum, std::vector<Spectrum>(SHTerms(shL), Spectrum(0.f)));
		int dx[8] = { -1, 1, 1, -1, -1, 1, 1, -1 };
		int dy[8] = { -1, -1, 1, 1, -1, -1, 1, 1 };
		int dz[8] = { -1, -1, -1, -1, 1, 1, 1, 1 };

		Float scale = 1.0 / 9.0;
		int endIdx = static_cast<int>(partitionNum - radius + 1);
#pragma omp parallel for schedule(static,1)
		for (int x = radius; x < endIdx; ++x) {
			for (int y = radius; y < endIdx; ++y) {
				for (int z = radius; z < endIdx; ++z) {
					int curIdx = GetIdx(x, y, z);
					//add center
					for (int i = 0; i < SHTerms(shL); ++i) {
						filterSHC[curIdx][i] = voxel[curIdx].shC[i] * scale;
					}
					//add edge
					for (int i = 0; i < 8; ++i) {
						std::vector<Spectrum>& curVoxelC = voxel[GetIdx(x + dx[i], y + dy[i], z + dz[i])].shC;
						for (int j = 0; j < SHTerms(shL); ++j) {
							filterSHC[curIdx][j] +=curVoxelC[j] * scale;
						}
					}

				}
			}
		}
		for (int i = 0; i < filterSHC.size(); ++i) {
			voxel[i].shC = std::move(filterSHC[i]);
		}
 	}

}