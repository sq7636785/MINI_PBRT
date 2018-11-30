#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif


#ifndef PBRT_SHAPES_SPHERE_H
#define PBRT_SHAPES_SPHERE_H

#include "shape.h"

namespace pbrt {
	class Sphere : public Shape {
	  public:
		  Sphere(const Transform *ObjectToWorld, const Transform *WorldToObject, bool reverseOrientation, Float radius, Float zMin, Float zMax, Float phiMax)
			  :Shape(ObjectToWorld, WorldToObject, reverseOrientation),
			  radius(radius),
			  zMin(Clamp(std::min(zMin, zMax), -radius, radius)),
			  zMax(Clamp(std::max(zMin, zMax), -radius, radius)),
			  thetaMin(std::acos(Clamp(std::min(zMin, zMax) / radius, -1, 1))),
			  thetaMax(std::acos(Clamp(std::max(zMin, zMax) / radius, -1, 1))),
			  phiMax(Radians(Clamp(phiMax, 0, 360))) { }

		  Bounds3f ObjectBound() const;
		  bool Intersect(const Ray& ray, Float* tHit, SurfaceInteraction* isect, bool testAlphaTexture) const;
		  bool IntersectP(const Ray& ray, bool testAphaTexture) const;
		  Float Area() const;
		  Interaction Sample(const Point2f &u, Float *pdf) const;
		  



	  private:
		const Float radius;
		const Float zMin;
		const Float zMax;
		const Float thetaMin;
		const Float thetaMax;
		const Float phiMax;

	};


	std::shared_ptr<Shape> CreateSphereShape(const Transform *o2w,
		const Transform *w2o,
		bool reverseOrientation,
		const ParamSet &params);

}





#endif