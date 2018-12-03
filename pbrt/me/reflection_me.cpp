
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// core/reflection.cpp*
#include "reflection.h"
#include "spectrum.h"
#include "sampler.h"
#include "sampling.h"
#include "interpolation.h"
#include "scene.h"
#include "interaction.h"
#include "stats.h"
#include <stdarg.h>

namespace pbrt {

// BxDF Utility Functions

Spectrum LambertianReflection::f(const Vector3f &wo, const Vector3f &wi) const {
    return R * InvPi;
}


Spectrum LambertianReflection::rho(int, const Point2f *, const Point2f *) const {
	return R;
}

Spectrum LambertianReflection::rho(const Vector3f &, int, const Point2f *) const {
	return R;
}

Spectrum BxDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                        Float *pdf, BxDFType *sampledType) const {
	//cosine weight hemisphere sample
	*wi = CosineSampleHemisphere(u);
	if (wi->z < 0) { wi->z *= -1; }
	*pdf = Pdf(wo, *wi);
	return f(wo, *wi);
}

Float BxDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    //cosine weight sample pdf cosTheta / Pi
	return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi: 0;
}


//no wo, sample nSample direcion and average
Spectrum BxDF::rho(const Vector3f &w, int nSamples, const Point2f *u) const {
	Spectrum r(0.f);
	for (int i = 0; i < nSamples; ++i) {
		Vector3f wi;
		Float pdf;
		Spectrum f = Sample_f(w, &wi, u[i], &pdf);
		r += f * AbsCosTheta(wi) / pdf;
	}
	return r / nSamples;
}

Spectrum BxDF::rho(int nSamples, const Point2f *u1, const Point2f *u2) const {
    Spectrum r(0.f);
    for (int i = 0; i < nSamples; ++i) {
        // Estimate one term of $\rho_\roman{hh}$
		Vector3f wo, wi;
		Float pdfo, pdfi;
		wo = CosineSampleHemisphere(u1[i]);
		pdfo = UniformHemispherePdf();
		pdfi = 0;
		Spectrum f = Sample_f(wo, &wi, u2[i], &pdfi);
		if (pdfi > 0) {

		}
    }
    return r / (Pi * nSamples);
}



Vector3f BSDF::WorldToLocal(const Vector3f &v) const {
	return Vector3f(Dot(v, ns), Dot(v, ts), Dot(v, ss));
}
Vector3f BSDF::LocalToWorld(const Vector3f &v) const {
	return Vector3f(
		ns.x * v.x + ts.x * v.y + ss.x * v.z,
		ns.y * v.x + ts.y * v.y + ss.y * v.z,
		ns.z * v.x + ts.z * v.y + ss.z * v.z
	);
}

// BSDF Method Definitions
Spectrum BSDF::f(const Vector3f &woW, const Vector3f &wiW,
                 BxDFType flags) const {
    ProfilePhase pp(Prof::BSDFEvaluation);
	Vector3f wi = WorldToLocal(wiW);
	Vector3f wo = WorldToLocal(woW);
	Spectrum r(0.f);
	if (wo.z == 0) {return r; }
	bool reflct = Dot(ng, wiW) * Dot(ng, woW) > 0;
	for (int i = 0; i < nBxDFs; ++i) {
		if (bxdfs[i]->MatchesFlags(flags) &&
			(reflct && (bxdfs[i]->type & BSDF_REFLECTION)) ||
			(!reflct) && (bxdfs[i]->type & BSDF_TRANSMISSION)) {
			r += bxdfs[i]->f(wo, wi);
		}
	}
	return r;
}

Spectrum BSDF::rho(int nSamples, const Point2f *samples1,
                   const Point2f *samples2, BxDFType flags) const {
    Spectrum ret(0.f);
	for (int i = 0; i < nBxDFs; ++i) {
		if (bxdfs[i]->MatchesFlags(flags)) {
			ret += bxdfs[i]->rho(nSamples, samples1, samples2);
		}
	}
    return ret;
}

Spectrum BSDF::rho(const Vector3f &wo, int nSamples, const Point2f *samples,
                   BxDFType flags) const {
    Spectrum ret(0.f);
	for (int i = 0; i < nBxDFs; ++i) {
		if (bxdfs[i]->MatchesFlags(flags)) {
			ret += bxdfs[i]->rho(wo, nSamples, samples);
		}
	}
    return ret;
}

Spectrum BSDF::Sample_f(const Vector3f &woWorld, Vector3f *wiWorld,
                        const Point2f &u, Float *pdf, BxDFType type,
                        BxDFType *sampledType) const {
    ProfilePhase pp(Prof::BSDFSampling);
    // Choose which _BxDF_ to sample
	int matchingComps = NumComponents(type);
	if (matchingComps == 0) {
		*pdf = 0;
		if (sampledType) {
			*sampledType = BxDFType(0);
		}
		return Spectrum(0.f); 
	}
	
	//随机选一个bxdf来计算wi
	int comp = std::min((int)std::floor(u[0] * matchingComps), matchingComps - 1);
	BxDF *bxdf = nullptr;
	int count = comp;
	for (int i = 0; i < nBxDFs; ++i) {
		if (bxdfs[i]->MatchesFlags(type) && count-- == 0) {
			bxdf = bxdfs[i];
			break;
		}
	}

	CHECK(bxdf != nullptr);
	// Remap _BxDF_ sample _u_ to $[0,1)^2$
	Point2f uRemapped(std::min(u[0] * matchingComps - comp, OneMinusEpsilon), u[1]);
	Vector3f wo = WorldToLocal(woWorld);
	if (wo.z == 0) { return Spectrum(0.f); }
	Vector3f wi;
	*pdf = 0;
	if (sampledType) { * sampledType = bxdf->type; }
	Spectrum f = bxdf->Sample_f(wo, &wi, uRemapped, pdf, sampledType);

	if (pdf == 0) {
		if (sampledType) {*sampledType = BxDFType(0); }
	}
	*wiWorld = LocalToWorld(wi);

	//compute pdf
	for (int i = 0; i < nBxDFs; ++i) {
		if (bxdfs[i] != bxdf && bxdfs[i]->MatchesFlags(type)) {
			*pdf += bxdfs[i]->Pdf(wo, wi);
		}
	}
	if (matchingComps > 1) { *pdf /= matchingComps; }

	//compute f
	if (!(bxdf->type & BSDF_SPECULAR)) {
		f = 0;
		bool reflect = Dot(*wiWorld, ng) * Dot(woWorld, ng) > 0;
		for (int i = 0; i < nBxDFs; ++i) {
			if (bxdfs[i]->MatchesFlags(type) &&
				((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
				(!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION)))) {
				f += bxdfs[i]->f(wo, wi);
			}
		}
	}
    return f;
}

Float BSDF::Pdf(const Vector3f &woWorld, const Vector3f &wiWorld,
                BxDFType flags) const {
    ProfilePhase pp(Prof::BSDFPdf);
	if (nBxDFs == 0) { return 0.f; }
	Float v = 0.f;
	int nMatchComponents = 0;
    Vector3f wi = WorldToLocal(wiWorld);
	Vector3f wo = WorldToLocal(woWorld);
	for (int i = 0; i < nBxDFs; ++i) {
		if (bxdfs[i]->MatchesFlags(flags)) {
			v += bxdfs[i]->Pdf(wo, wi);
			++nMatchComponents;
		}
	}
	v = nMatchComponents == 0 ? 0.f : v / nMatchComponents;
    return v;
}


}  // namespace pbrt
