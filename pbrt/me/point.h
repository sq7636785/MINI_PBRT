
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_LIGHTS_POINT_H
#define PBRT_LIGHTS_POINT_H

// lights/point.h*
#include "pbrt.h"
#include "light.h"
#include "shape.h"

namespace pbrt {

	// PointLight Declarations
	class PointLight : public Light {
	public:
		// PointLight Public Methods
		PointLight(const Transform &LightToWorld, const MediumInterface &mediumInterface, const Spectrum &I)
			:Light((int)LightFlags::DeltaPosition, LightToWorld, mediumInterface),
			pLight(LightToWorld(Point3f(0, 0, 0))),
			I(I) {
		}

		Spectrum Sample_Li(const Interaction &ref, const Point2f &u, Vector3f *wi, Float *pdf, VisibilityTester *vis) const;
		Spectrum Power() const;
		Float Pdf_Li(const Interaction& ref, const Vector3f& wi) const;
		Spectrum Sample_Le(const Point2f &u1, const Point2f &u2, Float time, Ray *ray, Normal3f *nLight, Float *pdfPos, Float *pdfDir) const;
		void Pdf_Le(const Ray &, const Normal3f &, Float *pdfPos, Float *pdfDir) const;


	private:
		const Point3f pLight;
		const Spectrum I;
	};

	std::shared_ptr<PointLight> CreatePointLight(const Transform &light2world,
		const Medium *medium,
		const ParamSet &paramSet);

}  // namespace pbrt

#endif  // PBRT_LIGHTS_POINT_H
