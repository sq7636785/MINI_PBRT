//light point
#include "me/point.h"
#include "scene.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {

// PointLight Method Definitions
	Spectrum PointLight::Sample_Li(const Interaction &ref, const Point2f &u, Vector3f *wi, Float *pdf, VisibilityTester *vis) const {
		ProfilePhase _(Prof::LightSample);
		*wi = Normalize(pLight - ref.p);
		*pdf = 1.f;
		//vistest 肯定是两个点之间的， 这里用interaction， 主要是生成光线， 需要注意的是这里的光线方向是没有被Normalize的。
		*vis = VisibilityTester(ref, Interaction(pLight, ref.time, mediumInterface));
		return I / DistanceSquared(pLight, ref.p);
	}

	Float PointLight::Pdf_Li(const Interaction& ref, const Vector3f& wi) const {
		return 0;
	}

	Spectrum PointLight::Power() const {
		return 4 * Pi * I;
	}

	Spectrum PointLight::Sample_Le(const Point2f &u1, const Point2f &u2, Float time, Ray *ray, Normal3f *nLight, Float *pdfPos, Float *pdfDir) const {
		ProfilePhase _(Prof::LightSample);
		*ray = Ray(pLight, UniformSampleSphere(u1), Infinity_, time, mediumInterface.inside);
		*nLight = Normal3f(ray->d);
		*pdfPos = 1;
		*pdfDir = UniformSpherePdf();
		return I;
	}

	void PointLight::Pdf_Le(const Ray &, const Normal3f &, Float *pdfPos, Float *pdfDir) const {
		ProfilePhase _(Prof::LightPdf);
		*pdfPos = 0;
		*pdfDir = UniformSpherePdf();
	}

	//add 
	void PointLight::SampleWi(const Point3f& it, Vector3f* wi) const {
		*wi = Normalize(pLight - it);
	}

	//add
	Spectrum PointLight::Li(const Interaction& ref, const Vector3f& w, Float* pdf, VisibilityTester* vis) const {
		Vector3f lightD = Normalize(pLight - ref.p);
		if (lightD == w) {
			*pdf = 1.f;
			*vis = VisibilityTester(ref, Interaction(pLight, ref.time, mediumInterface));
			return I;
		} else {
			*pdf = 0.f;
			return Spectrum(0.f);
		}
	}



	std::shared_ptr<PointLight> CreatePointLight(const Transform &light2world,
                                             const Medium *medium,
                                             const ParamSet &paramSet) {
    Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    Point3f P = paramSet.FindOnePoint3f("from", Point3f(0, 0, 0));
    Transform l2w = Translate(Vector3f(P.x, P.y, P.z)) * light2world;
    return std::make_shared<PointLight>(l2w, medium, I * sc);
}



}  // namespace pbrt
