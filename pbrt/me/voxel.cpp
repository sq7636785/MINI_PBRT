#include "voxel.h"
#include <fstream>

namespace pbrt {

	Volume::Volume(const Bounds3f& bound, Float partitionNum /*= 100.0*/) 
		: worldBound(bound), partitionNum(partitionNum) {
		Vector3f axisLen = worldBound.Diagonal();
		Float invPartitionNum = 1.0 / partitionNum;
		xDelta = axisLen.x * invPartitionNum;
		yDelta = axisLen.y * invPartitionNum;
		zDelta = axisLen.z * invPartitionNum;
	}

	bool Volume::ConstructVolume() {

		Point3f basePoint = worldBound.pMin;

//		std::cout << xDelta << ' ' << yDelta << ' ' << zDelta << std::endl;

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

	
	bool Volume::CalculateVoxel(const std::vector<std::shared_ptr<Primitive>>& curves, bool mode){
		int validVoxel = 0;
		
		Float validLength = std::min(xDelta, std::min(yDelta, zDelta)) / 2.0;
		Float invD = 1.0 / validLength;
		


		std::vector<CalUtil> innerData(voxel.size());

#pragma region __CurveByVoxel
		//对每一根头发去找他能覆盖到的体素的范围，再来计算。
		if (mode) {
			for (size_t curveId = 0; curveId < curves.size(); ++curveId) {
				auto boundSet = GetVoxelSet(*curves[curveId]);
				for (auto voxelId : boundSet) {

					Float distance;
					Vector3f direction;
					Float width;

					if (curves[curveId]->GetShape()->DistanceToPoint(voxel[voxelId].bound, &width, &distance, &direction)) {
						if (distance < validLength) {
							innerData[voxelId].distances.push_back(distance);
							innerData[voxelId].directions.push_back(direction);
							innerData[voxelId].width.push_back(width);
						}
					}
				}
				//std::cout << curveId << std::endl;
			}

			for (size_t idx = 0; idx < voxel.size(); ++idx) {
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
				//std::cout << idx << std::endl;
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
				//std::cout << idx << std::endl;
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
			out << voxel[i].sigma << ' ' << voxel[i].avgDirection.x << ' ' << voxel[i].avgDirection.y << ' ' << voxel[i].avgDirection.z << std::endl;
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
					if (idx != 0 && idx < total) {
						result.push_back(idx);
					}
				}
			}
		}
		return result;
	}



	inline
	int Volume::GetIdx(int x, int y, int z) const {
		return x + y * partitionNum + z * partitionNum * partitionNum;
	}
	inline
		void Volume::GetXYZ(int idx, int* x, int* y, int* z) const {
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

}