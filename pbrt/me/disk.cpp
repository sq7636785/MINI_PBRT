#include "me/disk.h"
#include "sampling.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"


namespace pbrt {

	Bounds3f Disk::ObjectBound() const {
		return Bounds3f(Point3f(-radius, -radius, height), Point3f(radius, radius, height));
	}

	Float Disk::Area() const {
		return phiMax * 0.5 * (radius * radius - innerRadius * innerRadius);
	}

	Interaction Disk::Sample(const Point2f &u, Float *pdf) const{
		Point2f pd = ConcentricSampleDisk(u);
		Point3f pObj = Point3f(pd.x * radius, pd.y * radius, height);
		Interaction it;
		it.n = Normalize((*ObjectToWorld)(Normal3f(0, 0, 1)));
		if (reverseOrientation) { it.n *= -1; }
		it.p = (*ObjectToWorld)(pObj, Vector3f(0, 0, 0), &it.pError);
		*pdf = 1 / Area();
		return it;
	}

	bool Disk::IntersectP(const Ray &r, bool testAlphaTexture) const {
		ProfilePhase p(Prof::ShapeIntersectP);

		Vector3f oErr, dErr;
		Ray ray = (*WorldToObject)(r, &oErr, &dErr);

		if (ray.d.z == 0) { return false; }
		Float tShapeHit = (height - ray.o.z) / ray.d.z;
		if (tShapeHit <=0 || tShapeHit >= ray.tMax) { return false; }

		Point3f pHit = ray(tShapeHit);
		Float dist2 = pHit.x * pHit.x + pHit.y * pHit.y;
		if (dist2 > radius * radius || dist2 < innerRadius * innerRadius) { return false; }

		Float phi = std::atan2(pHit.y, pHit.x);
		if (phi < 0) { phi += 2 * Pi; }
		if (phi > phiMax) { return false; }
		return true;
	}


	bool Disk::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect, bool testAlphaTexture) const {
		ProfilePhase p(Prof::ShapeIntersect);

		Vector3f dErr, oErr;
		Ray ray = (*WorldToObject)(r, &oErr, &dErr);

		if (ray.d.z == 0) { return false; }
		Float tShapeHit = (height - ray.o.z) / ray.d.z;
		if (tShapeHit <= 0 || tShapeHit >= ray.tMax) { return false; }
		
		Point3f pHit = ray(tShapeHit);
		Float dist2 = pHit.x * pHit.x + pHit.y * pHit.y;
		if (dist2 > radius * radius || dist2 < innerRadius * innerRadius) { return false; }

		Float phi = std::atan2(pHit.y, pHit.x);
		if (phi < 0) { phi += 2 * Pi; }
		if (phi > phiMax) { return false; }

		//dpdu dpdv
		Float u = phi / phiMax;
		Float rHit = std::sqrt(dist2);
		Float v = 1 - ((rHit - innerRadius) / (radius - innerRadius));
		Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
		Vector3f dpdv = Vector3f(pHit.x, pHit.y, 0.0) * (radius - innerRadius) / rHit;
		
		Normal3f dndu(0, 0, 0);
		Normal3f dndv(0, 0, 0);

		pHit.z = height;
		Vector3f pError(0, 0, 0);

		*isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v), -ray.d, dpdu, dpdv, dndu, dndv, ray.time, this));
		*tHit = (Float)tShapeHit;
		return true;
	}

	std::shared_ptr<Shape> CreateDiskShape(const Transform *o2w,
	            	                       const Transform *w2o,
		                                   bool reverseOrientation,
		                                   const ParamSet &params) {
		Float height = params.FindOneFloat("height", 0.0);
		Float radius = params.FindOneFloat("radius", 1.0);
		Float innerRadius = params.FindOneFloat("innerradius", 0);
		Float phimax = params.FindOneFloat("phimax", 360.0);
		return std::make_shared<Disk>(o2w, w2o, reverseOrientation, height, radius,innerRadius, phimax);
	}


}
