#include <stdlib.h>

#include "scene.h"
#include "interaction.h"


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
				
			//	std::cout << "before Tr: " << Li << std::endl; 
				auto tr = vis.Tr(*this);
				//std::cout << tr << std::endl;
				Li *= tr;
			//	std::cout << "after Tr: " << Li << std::endl;
				//irrandiance
				Float cosTheta = 1.0;
				
				//平均方向不为0. 说明有头发， 那么取平均头发方向的切向为法向量.
				if (v.avgDirection.LengthSquared() != 0) {
					Float sinTheta = AbsDot(v.avgDirection, wi);
					cosTheta = std::sqrt(1.0 - sinTheta * sinTheta);
				}
				Li *= cosTheta;
		//		std::cout << "after Dot: " << Li << std::endl << std::endl;
				Float tmpRGB[3];
				Li.ToRGB(tmpRGB);
				v.rgb[0] += tmpRGB[0];
				v.rgb[1] += tmpRGB[1];
				v.rgb[2] += tmpRGB[2];
			}
		}
	}
}