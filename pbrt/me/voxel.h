#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif


#ifndef PBRT_VOXEL_H
#define PBRT_VOXEL_H

#include "geometry.h"
#include "primitive.h"
#include "light.h"

namespace pbrt {

	struct CalUtil {
		CalUtil() {}
		std::vector<Float>		distances;
		std::vector<Float>		width;
		std::vector<Vector3f>	directions;
	};

	struct Voxel {
		Voxel() : sigma(0.0),directionV(0.0) {
			rgb[0] = rgb[1] = rgb[2] = 0.0;
		}
		Bounds3f	bound;
		Float		sigma;
		Vector3f	avgDirection;
		Float		directionV;
		Float       rgb[3];
	};

	class Volume {
	  public:
		Volume(const Bounds3f& bound, Float partitionNum = 100.0);

		bool ConstructVolume();
		bool CalculateVoxel(const std::vector<std::shared_ptr<Primitive>>& curves, bool mode = false);
		Spectrum Tr(const Point3f& p0, const Point3f& p1) const ;
		void SaveData(std::string fileName) const;
		void LoadData(std::string fileName);
		std::vector<int> GetVoxelSet(const Primitive& curve) const;


		std::vector<Voxel>	voxel;

	  private:
		int GetIdx(int x, int y, int z) const;
		void GetXYZ(const int idx, int* x, int* y, int* z) const;
		int GetIdxFromPoint(Float x, Float y, Float z) const;

		Float				partitionNum;
		Bounds3f			worldBound;

		Float				xDelta;
		Float				yDelta;
		Float				zDelta;
	};

}





#endif