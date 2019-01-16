#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif


#ifndef PBRT_VOXEL_H
#define PBRT_VOXEL_H

#include "geometry.h"
#include "primitive.h"
#include "light.h"
#include "sh.h"
#include "reflection.h"

namespace pbrt {

	struct CalUtil {
		CalUtil() {}
		std::vector<Float>		distances;
		std::vector<Float>		width;
		std::vector<Vector3f>	directions;
	};

	struct Voxel {
		Voxel() : sigma(0.0), directionV(0.0), lightLength(0.0) {
			rgb[0] = rgb[1] = rgb[2] = 0.0;
		}
		Bounds3f				bound;
		Float					sigma;
		Vector3f				avgDirection;
		Float					directionV;
		Float					rgb[3];
		Float					lightLength;
		std::vector<Spectrum>	shC;
	};

	struct SHSample {
		SHSample(){}
		Vector3f			w;
		std::vector<Float>	y;
	};

	class Volume {
	  public:
		Volume(const Bounds3f& bound, Float partitionNum = 100.0, int shL = 4, int nSHSample = 10000);

		bool ConstructVolume();
		bool CalculateVoxel(const std::vector<std::shared_ptr<Primitive>>& curves, bool mode = false);
		Spectrum Tr(const Point3f& p0, const Point3f& p1) const ;
		void SaveData(std::string fileName) const;
		void LoadData(std::string fileName);
		void LoadBsdfMatrix(std::string fileName);

		std::vector<int> GetVoxelSet(const Primitive& curve) const;
		int GetIdx(int x, int y, int z) const;
		void GetXYZ(int idx, int* x, int* y, int* z) const;
		int GetIdxFromPoint(Float x, Float y, Float z) const;
		void SetSHPara(int lmax, int nSample);

		void UpdateSHFromLightSegment(const Vector3f& wi, Float length, int idx, const Spectrum& power);

		void ComputeBSDFMatrix(BxDF& bsdf, int nSamples = 50*50);
		void RotateSH(const Vector3f& v1, const Vector3f& v2, Spectrum* cIn, Spectrum* cOut) const;
		void SHMatrixTransV(const std::vector<Spectrum>& c, Spectrum* c_out) const;

		Float GetMinDimDelta() const;
		Float GetMaxDimDelta() const;

		std::vector<SHSample>	shSample;
		std::vector<Spectrum>	bsdfMatrix;
		std::vector<Voxel>		voxel;
		int						shL;
		int						nSHSample;


	  private:
		

		Float					partitionNum;
		Bounds3f				worldBound;

		Float					xDelta;
		Float					yDelta;
		Float					zDelta;
		Float					vDelta;
	};


	inline
		int Volume::GetIdx(int x, int y, int z) const {
		int idx = x + y * partitionNum + z * partitionNum * partitionNum;
		if (idx < 0 || idx >= partitionNum * partitionNum *partitionNum) {
			return -1;
		}
		return x + y * partitionNum + z * partitionNum * partitionNum;
	}
	inline
		void Volume::GetXYZ(int idx, int* x, int* y, int* z) const {
		if (idx < 0 || idx >= partitionNum * partitionNum *partitionNum) {
			*x = 0;
			*y = 0;
			*z = 0;
			return;
		}
		int pN = static_cast<int>(partitionNum);
		*z = idx / (pN * pN);
		*y = (idx % (pN * pN)) / pN;
		*x = idx % pN;
	}

	inline
		int Volume::GetIdxFromPoint(Float x, Float y, Float z) const {
		int xIdx = static_cast<int>((x - worldBound.pMin.x) / xDelta);
		int yIdx = static_cast<int>((y - worldBound.pMin.y) / yDelta);
		int zIdx = static_cast<int>((z - worldBound.pMin.z) / zDelta);
		return GetIdx(xIdx, yIdx, zIdx);
	}


	inline
		void Volume::SetSHPara(int lmax, int nSample) {

	}

	static Spectrum SphericalLightFunc(Vector3f& w) {
		Float theta = SphericalTheta(w);

		Spectrum l(1.f);
		l *= exp(-std::abs(theta));
		//l *= (std::max(0.0, 5.0 * std::cos(theta) - 4.0) + std::max(0.0, -4.0 * std::sin(theta - Pi) * std::cos(phi - 2.5) - 3.0));
		return l;
	}

	inline 
		Float Volume::GetMinDimDelta() const {
		Float marchSize = xDelta > yDelta ? (yDelta > zDelta ? zDelta : yDelta) :
											(xDelta > zDelta ? zDelta : xDelta);
		return marchSize;
	}

	inline
		Float Volume::GetMaxDimDelta() const {
		Float marchSize = xDelta < yDelta ? (yDelta < zDelta ? zDelta : yDelta) :
											(xDelta < zDelta ? zDelta : xDelta);
		return marchSize;
	}
}





#endif