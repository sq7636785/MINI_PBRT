#include "shape.h"
#include "stats.h"

namespace pbrt {
	Interaction Shape::Sample(const Interaction &ref, const Point2f &u,
		Float *pdf) const {
		Interaction intr = Sample(u, pdf);
		Vector3f wi = intr.p - ref.p;
		if (wi.LengthSquared() == 0) {
			*pdf = 0.f;
		} else {
			*pdf *= DistanceSquared(ref.p, intr.p) / AbsDot(intr.n, -wi);
			if (std::isnan(*pdf)) {
				*pdf = 0.f;
			}
		}
		return intr;
	}

	Float Shape::Pdf(const Interaction &ref, const Vector3f &wi) const {
		// Intersect sample ray with area light geometry
		Ray ray = ref.SpawnRay(wi);
		Float tHit;
		SurfaceInteraction isectLight;
		if (!Intersect(ray, &tHit, &isectLight, false)) {
			return 0.f;
		}
		Float pdf =  DistanceSquared(ref.p, isectLight.p) / (AbsDot(isectLight.n, -wi) * Area());
		if (std::isnan(pdf)) {
			pdf = 0.f;
		}
		return pdf;
	}
}