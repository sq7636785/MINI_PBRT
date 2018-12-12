#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif


#ifndef PBRT_VOXEL_H
#define PBRT_VOXEL_H

#include "geometry.h"
#include "primitive.h"

namespace pbrt {

	struct CalUtil {
		std::vector<Float>		distances;
		std::vector<Float>		width;
		std::vector<Vector3f>	directions;
	};

	struct Voxel {
		Bounds3f	bound;
		Float		sigma;
		Vector3f	avgDirection;
		Float		directionV;
	};

	class Volume {
	  public:
		Volume(const Bounds3f& bound, const Float partitionNum = 100.0) 
			: worldBound(bound), partitionNum(partitionNum)
		{ }

		bool ConstructVolume();
		bool CalculateVoxel(const std::vector<std::shared_ptr<Primitive>>& curves);
		Spectrum Tr(const Point3f& p0, const Point3f& p1);
		 
	  private:
		std::vector<Voxel>	voxel;
		Float				partitionNum;
		Bounds3f			worldBound;
	};

}





#endif