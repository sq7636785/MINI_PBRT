#include "me/sphere.h"
#include "sampling.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"


namespace pbrt {
	Bounds3f Sphere::ObjectBound() const {
		return Bounds3f(Point3f(-radius, -radius, zMin), Point3f(radius, radius, zMax));
	}

	Float Sphere::Area() const {
		return phiMax * radius * (zMax - zMin);
	}

	Interaction Sphere::Sample(const Point2f &u, Float *pdf)  const {
		Point3f pObj = Point3f(0, 0, 0) + radius * UniformSampleHemisphere(u);
		Interaction it;
		it.n = Normalize((*ObjectToWorld)(Normal3f(pObj.x, pObj.y, pObj.z)));
		if (reverseOrientation) { it.n *= -1; }

		//pError
		pObj *= radius / Distance(pObj, Point3f(0, 0, 0));
		Vector3f pObjError = gamma(5) * Abs(Vector3f(pObj));
		it.p = (*ObjectToWorld)(pObj, pObjError, &it.pError);
		*pdf = 1 / Area();
		return it;
	}

	bool Sphere::Intersect(const Ray& r, Float* tHit, SurfaceInteraction* isect, bool testAlphaTexture) const {
		ProfilePhase p(Prof::ShapeIntersect);
		Float phi;
		Point3f pHit;
		Vector3f oErr, dErr;
		Ray ray = (*WorldToObject)(r, &oErr, &dErr);

		EFloat ox(ray.o.x, oErr.x);
		EFloat oy(ray.o.y, oErr.y);
		EFloat oz(ray.o.z, oErr.z);

		EFloat dx(ray.d.x, dErr.x);
		EFloat dy(ray.d.y, dErr.y);
		EFloat dz(ray.d.z, dErr.z);

		//coefficient
		EFloat a = dx * dx + dy * dy + dz * dz;
		EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
		EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

		//solve
		EFloat t0, t1;
		if (!Quadratic(a, b, c, &t0, &t1)) { return false; }

		// Check quadric shape _t0_ and _t1_ for nearest intersection
		if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
		EFloat tShapeHit = t0;
		if (tShapeHit.LowerBound() <= 0) {
			tShapeHit = t1;
			if (tShapeHit.UpperBound() > ray.tMax) return false;
		}


		pHit = ray((Float)tShapeHit);

		//refine
		pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
		if (pHit.x == 0 && pHit.y == 0) {
			pHit.x = 1e-5f * radius;
		}
		phi = std::atan2(pHit.y, pHit.x);
		if (phi < 0) {
			phi += 2 * Pi;
		}

		// Test sphere intersection against clipping parameters
		if ((zMin > -radius && pHit.z < zMin) || (zMax < radius && pHit.z > zMax) ||
			phi > phiMax) {
			if (tShapeHit == t1) { return false; }
			if (t1.UpperBound() > ray.tMax) { return false; }

			tShapeHit = t1;
			// Compute sphere hit position and $\phi$
			pHit = ray((Float)tShapeHit);

			// Refine sphere intersection point
			pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
			if (pHit.x == 0 && pHit.y == 0) { pHit.x = 1e-5f * radius; }
			phi = std::atan2(pHit.y, pHit.x);
			if (phi < 0) phi += 2 * Pi;
			if ((zMin > -radius && pHit.z < zMin) ||
				(zMax < radius && pHit.z > zMax) || phi > phiMax)
				return false;
		}

		//dpdu dpdv
		Float u = phi / phiMax;
		Float theta = std::acos(Clamp(pHit.z / radius, -1, 1));
		Float v = (theta - thetaMin) / (thetaMax - thetaMin);

		Float zRadius = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
		Float invZRadius = 1 / zRadius;
		Float cosPhi = pHit.x * invZRadius;
		Float sinPhi = pHit.y * invZRadius;
		
		Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
		Vector3f dpdv = (thetaMax - thetaMin) * Vector3f(pHit.z * cosPhi, pHit.z * sinPhi, -radius * std::sin(theta));

		Vector3f pError = gamma(5) * Abs((Vector3f)pHit);

		*isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v), -ray.d, dpdu, dpdv, Normal3f(), Normal3f(), ray.time, this));
		*tHit = (Float)tShapeHit;
		return true;
	}

	bool Sphere::IntersectP(const Ray& r, bool testAphaTexture) const {
		ProfilePhase p(Prof::ShapeIntersectP);
		Float phi;
		Point3f pHit;
		Vector3f oErr, dErr;
		Ray ray = (*WorldToObject)(r, &oErr, &dErr);

		EFloat ox(ray.o.x, oErr.x);
		EFloat oy(ray.o.y, oErr.y);
		EFloat oz(ray.o.z, oErr.z);

		EFloat dx(ray.d.x, dErr.x);
		EFloat dy(ray.d.y, dErr.y);
		EFloat dz(ray.d.z, dErr.z);

		//coefficient
		EFloat a = dx * dx + dy * dy + dz * dz;
		EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
		EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

		//solve
		EFloat t0, t1;
		if (!Quadratic(a, b, c, &t0, &t1)) { return false; }

		// Check quadric shape _t0_ and _t1_ for nearest intersection
		if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
		EFloat tShapeHit = t0;
		if (tShapeHit.LowerBound() <= 0) {
			tShapeHit = t1;
			if (tShapeHit.UpperBound() > ray.tMax) return false;
		}


		pHit = ray((Float)tShapeHit);

		//refine
		pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
		if (pHit.x == 0 && pHit.y == 0) {
			pHit.x = 1e-5f * radius;
		}
		phi = std::atan2(pHit.y, pHit.x);
		if (phi < 0) {
			phi += 2 * Pi;
		}

		// Test sphere intersection against clipping parameters
		if ((zMin > -radius && pHit.z < zMin) || (zMax < radius && pHit.z > zMax) ||
			phi > phiMax) {
			if (tShapeHit == t1) { return false; }
			if (t1.UpperBound() > ray.tMax) { return false; }

			tShapeHit = t1;
			// Compute sphere hit position and $\phi$
			pHit = ray((Float)tShapeHit);

			// Refine sphere intersection point
			pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
			if (pHit.x == 0 && pHit.y == 0) { pHit.x = 1e-5f * radius; }
			phi = std::atan2(pHit.y, pHit.x);
			if (phi < 0) phi += 2 * Pi;
			if ((zMin > -radius && pHit.z < zMin) ||
				(zMax < radius && pHit.z > zMax) || phi > phiMax)
				return false;
		}
		return true;
	}

// 	std::shared_ptr<Shape> CreatSphereShape(const Transform* o2w,
// 		const Transform* w2o, bool reverseOrientation, const ParamSet& params) {
// 		Float radius = params.FindOneFloat("radius", 1.0f);
// 		Float zmin = params.FindOneFloat("zmin", -radius);
// 		Float zmax = params.FindOneFloat("zmax", radius);
// 		Float phimax = params.FindOneFloat("phimax", 360.0f);
// 
// 		return std::make_shared<Sphere>(o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax);
// 	}
	std::shared_ptr<Shape> CreateSphereShape(const Transform *o2w,
		const Transform *w2o,
		bool reverseOrientation,
		const ParamSet &params) {
		Float radius = params.FindOneFloat("radius", 1.f);
		Float zmin = params.FindOneFloat("zmin", -radius);
		Float zmax = params.FindOneFloat("zmax", radius);
		Float phimax = params.FindOneFloat("phimax", 360.f);
		return std::make_shared<Sphere>(o2w, w2o, reverseOrientation, radius, zmin,
			zmax, phimax);
	}
}
