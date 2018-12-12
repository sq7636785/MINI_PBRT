
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SHAPES_CURVE_H
#define PBRT_SHAPES_CURVE_H

// shapes/curve.h*
#include "shape.h"

namespace pbrt {

	enum class CurveType {
		Flat,
		Cylinder,
		Ribbon
	};

	struct CurveCommon {
		CurveCommon(const Point3f c[4], Float w0, Float w1, CurveType type, const Normal3f* norm);
		const CurveType type;
		Point3f cpObj[4];
		Float width[2];
		Normal3f n[2];
		Float normalAngle;
		Float invSinNormalAngle;
	};


	class Curve : public Shape {
	  public:
		Curve(const Transform* ObjectToWorld, const Transform* WorldToObject,
			  bool reverseOrientation, const std::shared_ptr<CurveCommon>& common,
			  Float uMin, Float uMax)
			:Shape(ObjectToWorld, WorldToObject, reverseOrientation),
			 common(common),
		     uMin(uMin),
			 uMax(uMax) { }
		
		Bounds3f ObjectBound() const;
		bool Intersect(const Ray& ray, Float* tHit, SurfaceInteraction* isect, bool testAlphaTexture) const;
		Float Area() const;
		Interaction Sample(const Point2f &u, Float *pdf) const;
		bool DistanceToPoint(Bounds3f &worldVoxel, Float* width, Float* distance, Vector3f* direction) const;
	  private:
		  bool recursiveIntersect(const Ray& r, Float* tHit, SurfaceInteraction* isect, 
			                      const Point3f cp[4], const Transform& rayToObject,
			                      Float u0, Float u1, int depth) const;
		

	private:
		const std::shared_ptr<CurveCommon> common;
		const Float uMin, uMax;
	};



std::vector<std::shared_ptr<Shape>> CreateCurveShape(const Transform *o2w,
                                                     const Transform *w2o,
                                                     bool reverseOrientation,
                                                     const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SHAPES_CURVE_H
