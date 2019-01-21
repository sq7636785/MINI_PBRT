#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif


#ifndef PBRT_SHAPES_DISK_H
#define PBRT_SHAPES_DISK_H

#include "shape.h"

namespace pbrt {
	class Disk : public Shape {
	  public:
		Disk(const Transform *ObjectToWorld, const Transform *WorldToObject, bool reverseOrientation, Float height, Float radius, Float innerRadius, Float phiMax)
			:Shape(ObjectToWorld, WorldToObject, reverseOrientation),
			 height(height),
			 radius(radius),
			 innerRadius(innerRadius),
			 phiMax(Radians(Clamp(phiMax,0, 360))) { }
		Bounds3f ObjectBound() const;
		bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
		bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect, bool testAlphaTexture) const;
		Float Area() const;
		Interaction Sample(const Point2f &u, Float *pdf) const;

	  private:
		const Float height;
		const Float radius;
		const Float innerRadius;
		const Float phiMax;
	};

	std::shared_ptr<Shape> CreateDiskShape(const Transform *o2w,
		const Transform *w2o,
		bool reverseOrientation,
		const ParamSet &params);

}





#endif