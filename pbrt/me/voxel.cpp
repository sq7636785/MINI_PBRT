#include "voxel.h"


namespace pbrt {

	bool Volume::ConstructVolume() {
		Vector3f axisLen = worldBound.Diagonal();
		Float invPartitionNum = 1.0 / partitionNum;
		Float xDelta = axisLen.x * invPartitionNum;
		Float yDelta = axisLen.y * invPartitionNum;
		Float zDelta = axisLen.z * invPartitionNum;

		Point3f basePoint = worldBound.pMin;

		for (int z = 1; z <= partitionNum; ++z) {
			for (int y = 1; y <= partitionNum; ++y) {
				for (int x = 1; x <= partitionNum; ++x) {
					Voxel tmp;
					Point3f curPointMin = basePoint + Point3f((x - 1) * xDelta, (y - 1) * yDelta, (z - 1) * zDelta);
					Point3f curPointMax = basePoint + Point3f(x * xDelta, y * yDelta, z * zDelta);
					tmp.bound = Bounds3f(curPointMin, curPointMax);
					voxel.push_back(tmp);
				}
			}
		}
		return (voxel[0].bound.pMin - worldBound.pMin).Length() < 0.001 && 
			   (voxel[voxel.size() - 1].bound.pMax - worldBound.pMax).Length() < 0.001;
		
	}

	
	bool Volume::CalculateVoxel(const std::vector<std::shared_ptr<Primitive>>& curves){
		int validVoxel = 0;
		std::vector<CalUtil> innerData(voxel.size());
		for (size_t idx = 0; idx < voxel.size(); ++idx) {
			Vector3f diagnal = voxel[idx].bound.Diagonal();
			Float validLength = std::min(diagnal.x, std::min(diagnal.y, diagnal.z)) / 2.0;
			Float invD = 1.0 / validLength;
			for (size_t curveId = 0; curveId < curves.size(); ++curveId) {
				Float distance;
				Vector3f derection;
				Float width;
				if (curves[curveId]->GetShape()->DistanceToPoint(voxel[idx].bound, &width, &distance, &derection)) {
					if (distance < validLength) {
					innerData[idx].distances.push_back(distance);
					innerData[idx].directions.push_back(derection);
					innerData[idx].width.push_back(width);
					}
				}
			}
			Float Num = static_cast<Float>(innerData[idx].distances.size());
			if (Num != 0.0) {
				Float avgDiameter = 0.0;
				for (size_t i = 0; i < innerData[idx].distances.size(); ++i) {
					Float weight = 1 - innerData[idx].distances[i] * invD;
					voxel[idx].avgDirection += innerData[idx].directions[i] * weight;
					avgDiameter += innerData[idx].width[i];
				}
				avgDiameter /= Num;
				voxel[idx].sigma = 2.0 * avgDiameter * Num * InvPi * invD * invD;
				++validVoxel;
			} else {
				voxel[idx].sigma = 0.0;
			}
		}

		return validVoxel > 0;
	}

	pbrt::Spectrum Volume::Tr(const Point3f& p0, const Point3f& p1) {
		return Spectrum(0.0);
	}

}