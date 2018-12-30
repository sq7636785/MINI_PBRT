
/*
pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.
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


// core/sh.cpp*

#include "sh.h"
#include "sampling.h"
#include "scene.h"
#include "integrator.h"
#include "imageio.h"

namespace pbrt {
	// Spherical Harmonics Local Definitions
	static void legendrep(Float x, int lmax, Float *out) {
#define P(l,m) out[SHIndex(l,m)]
		// Compute $m=0$ Legendre values using recurrence
		P(0, 0) = 1.f;
		P(1, 0) = x;
		for (int l = 2; l <= lmax; ++l) {
			P(l, 0) = ((2 * l - 1)*x*P(l - 1, 0) - (l - 1)*P(l - 2, 0)) / l;
			assert(!std::isnan(P(l, 0)));
			assert(!std::isinf(P(l, 0)));
		}

		// Compute $m=l$ edge using Legendre recurrence
		Float neg = -1.f;
		Float dfact = 1.f;
		Float xroot = std::sqrt(std::max(0.0, 1.0 - x * x));
		Float xpow = xroot;
		for (int l = 1; l <= lmax; ++l) {
			P(l, l) = neg * dfact * xpow;
			assert(!std::isnan(P(l, l)));
			assert(!std::isinf(P(l, l)));
			neg *= -1.f;      // neg = (-1)^l
			dfact *= 2 * l + 1; // dfact = (2*l-1)!!
			xpow *= xroot;    // xpow = powf(1.f - x*x, Float(l) * 0.5f);
		}

		// Compute $m=l-1$ edge using Legendre recurrence
		for (int l = 2; l <= lmax; ++l) {
			P(l, l - 1) = x * (2 * l - 1) * P(l - 1, l - 1);
			assert(!std::isnan(P(l, l - 1)));
			assert(!std::isinf(P(l, l - 1)));
		}

		// Compute $m=1, \ldots, l-2$ values using Legendre recurrence
		for (int l = 3; l <= lmax; ++l)
			for (int m = 1; m <= l - 2; ++m) {
				P(l, m) = ((2 * (l - 1) + 1) * x * P(l - 1, m) -
					(l - 1 + m) * P(l - 2, m)) / (l - m);
				assert(!std::isnan(P(l, m)));
				assert(!std::isinf(P(l, m)));
			}
#if 0
		// wrap up with the negative m ones now
		// P(l,-m)(x) = -1^m (l-m)!/(l+m)! P(l,m)(x)
		for (int l = 1; l <= lmax; ++l) {
			Float fa = 1.f, fb = fact(2 * l);
			// fa = fact(l+m), fb = fact(l-m)
			for (int m = -l; m < 0; ++m) {
				Float neg = ((-m) & 0x1) ? -1.f : 1.f;
				P(l, m) = neg * fa / fb * P(l, -m);
				fb /= l - m;
				fa *= (l + m + 1) > 1 ? (l + m + 1) : 1.;
			}
		}
#endif
#undef P
	}


	static inline Float fact(Float v);
	static inline Float divfact(int a, int b);
	static inline Float K(int l, int m) {
		return std::sqrt((2.f * l + 1.f) * Inv4Pi * divfact(l, m));
	}


	static inline Float divfact(int a, int b) {
		if (b == 0) return 1.f;
		Float fa = a, fb = fabsf(b);
		Float v = 1.f;
		for (Float x = fa - fb + 1.f; x <= fa + fb; x += 1.f)
			v *= x;
		return 1.f / v;
	}


	// n!! = 1 if n==0 or 1, otherwise n * (n-2)!!
	static Float dfact(Float v) {
		if (v <= 1.f) return 1.f;
		return v * dfact(v - 2.f);
	}


	static inline Float fact(Float v) {
		if (v <= 1.f) return 1.f;
		return v * fact(v - 1.f);
	}


	static void sinCosIndexed(Float s, Float c, int n,
		Float *sout, Float *cout) {
		Float si = 0, ci = 1;
		for (int i = 0; i < n; ++i) {
			// Compute $\sin{}i\phi$ and $\cos{}i\phi$ using recurrence
			*sout++ = si;
			*cout++ = ci;
			Float oldsi = si;
			si = si * c + ci * s;
			ci = ci * c - oldsi * s;
		}
	}


	static void toZYZ(const Matrix4x4 &m, Float *alpha, Float *beta, Float *gamma) {
#define M(a, b) (m.m[a][b])

		Float sy = std::sqrt(M(2, 1)*M(2, 1) + M(2, 0)*M(2, 0));
		if (sy > 16 * FLT_EPSILON) {
			*gamma = -atan2f(M(1, 2), -M(0, 2));
			*beta = -atan2f(sy, M(2, 2));
			*alpha = -atan2f(M(2, 1), M(2, 0));
		} else {
			*gamma = 0;
			*beta = -atan2f(sy, M(2, 2));
			*alpha = -atan2f(-M(1, 0), M(1, 1));
		}
#undef M
	}


	static inline Float lambda(Float l) {
		return std::sqrt((4.f * Pi) / (2.f * l + 1.));
	}


	void SHComputeDiffuseTransfer(const Point3f &p, const Normal3f &n,
		Float rayEpsilon, const Scene *scene, RNG &rng, int nSamples,
		int lmax, Spectrum *c_transfer) {
		for (int i = 0; i < SHTerms(lmax); ++i)
			c_transfer[i] = 0.f;

		Float *Ylm = ALLOCA(Float, SHTerms(lmax));
		for (int i = 0; i < nSamples; ++i) {
			// Sample _i_th direction and compute estimate for transfer coefficients
			Point2f u(rng.UniformFloat(), rng.UniformFloat());
			Vector3f w = UniformSampleSphere(u);
			Float pdf = UniformSpherePdf();
			if (Dot(w, n) > 0.f && !scene->IntersectP(Ray(p, w, rayEpsilon))) {
				// Accumulate contribution of direction $\w{}$ to transfer coefficients
				SHEvaluate(w, lmax, Ylm);
				for (int j = 0; j < SHTerms(lmax); ++j)
					c_transfer[j] += (Ylm[j] * AbsDot(w, n)) / (pdf * nSamples);
			}
		}
	}


	// Spherical Harmonics Definitions
	void SHEvaluate(const Vector3f &w, int lmax, Float *out) {
		if (lmax > 28) {
			Error("SHEvaluate() runs out of numerical precision for lmax > 28. "
				"If you need more bands, try recompiling using doubles.");
			exit(1);
		}

		// Compute Legendre polynomial values for $\cos\theta$
		assert(w.Length() > .995f && w.Length() < 1.005f);
		legendrep(w.z, lmax, out);

		// Compute $K_l^m$ coefficients
		Float *Klm = ALLOCA(Float, SHTerms(lmax));
		for (int l = 0; l <= lmax; ++l)
			for (int m = -l; m <= l; ++m)
				Klm[SHIndex(l, m)] = K(l, m);

		// Compute $\sin\phi$ and $\cos\phi$ values
		Float *sins = ALLOCA(Float, lmax + 1), *coss = ALLOCA(Float, lmax + 1);
		Float xyLen = std::sqrt(std::max(0.0, 1.0 - w.z*w.z));
		if (xyLen == 0.f) {
			for (int i = 0; i <= lmax; ++i) sins[i] = 0.f;
			for (int i = 0; i <= lmax; ++i) coss[i] = 1.f;
		} else
			sinCosIndexed(w.y / xyLen, w.x / xyLen, lmax + 1, sins, coss);

		// Apply SH definitions to compute final $(l,m)$ values
		static const Float sqrt2 = std::sqrt(2.f);
		for (int l = 0; l <= lmax; ++l) {
			for (int m = -l; m < 0; ++m) {
				out[SHIndex(l, m)] = sqrt2 * Klm[SHIndex(l, m)] *
					out[SHIndex(l, -m)] * sins[-m];
				assert(!std::isnan(out[SHIndex(l, m)]));
				assert(!std::isinf(out[SHIndex(l, m)]));
			}
			out[SHIndex(l, 0)] *= Klm[SHIndex(l, 0)];
			for (int m = 1; m <= l; ++m) {
				out[SHIndex(l, m)] *= sqrt2 * Klm[SHIndex(l, m)] * coss[m];
				assert(!std::isnan(out[SHIndex(l, m)]));
				assert(!std::isinf(out[SHIndex(l, m)]));
			}
		}
	}


#if 0
	// Believe this is correct, but not well tested
	void SHEvaluate(Float costheta, Float cosphi, Float sinphi, int lmax, Float *out) {
		legendrep(costheta, lmax, out);

		Float *Klm = ALLOCA(Float, SHTerms(lmax));
		klm(lmax, Klm);

		Float sqrt2 = sqrtf(2.f);
		Float sins[(lmax + 1)], coss[(lmax + 1)];
		sinCosIndexed(sinphi, cosphi, lmax + 1, sins, coss);

		for (int l = 0; l <= lmax; ++l) {
			for (int m = -l; m < 0; ++m)
				// sin(-x) = -sin(x)
				out[SHIndex(l, m)] = sqrt2 * Klm[SHIndex(l, m)] *
				out[SHIndex(l, -m)] * sins[-m];

			out[SHIndex(l, 0)] *= Klm[SHIndex(l, 0)];

			for (int m = 1; m <= l; ++m)
				out[SHIndex(l, m)] *= sqrt2 * Klm[SHIndex(l, m)] * coss[m];
		}
	}


#endif
	

	void SHProjectIncidentDirectRadiance(const Point3f &p, Float pEpsilon,
		Float time, MemoryArena &arena, const Scene *scene,
		bool computeLightVis, int lmax, RNG &rng, Spectrum *c_d) {
		
	}



	void SHReduceRinging(Spectrum *c, int lmax, Float lambda) {
		for (int l = 0; l <= lmax; ++l) {
			Float scale = 1.f / (1.f + lambda * l * l * (l + 1) * (l + 1));
			for (int m = -l; m <= l; ++m)
				c[SHIndex(l, m)] *= scale;
		}
	}


	void SHRotate(const Spectrum *c_in, Spectrum *c_out, const Matrix4x4 &m,
		int lmax) {
		Float alpha, beta, gamma;
		toZYZ(m, &alpha, &beta, &gamma);
		Spectrum *work = ALLOCA(Spectrum, SHTerms(lmax));
		SHRotateZ(c_in, c_out, gamma, lmax);
		SHRotateXPlus(c_out, work, lmax);
		SHRotateZ(work, c_out, beta, lmax);
		SHRotateXMinus(c_out, work, lmax);
		SHRotateZ(work, c_out, alpha, lmax);
	}


	void SHRotateZ(const Spectrum *c_in, Spectrum *c_out, Float alpha,
		int lmax) {
		assert(c_in != c_out);
		c_out[0] = c_in[0];
		if (lmax == 0) return;
		// Precompute sine and cosine terms for $z$-axis SH rotation
		Float *ct = ALLOCA(Float, lmax + 1);
		Float *st = ALLOCA(Float, lmax + 1);
		sinCosIndexed(sinf(alpha), cosf(alpha), lmax + 1, st, ct);
		for (int l = 1; l <= lmax; ++l) {
			// Rotate coefficients for band _l_ about $z$
			for (int m = -l; m < 0; ++m)
				c_out[SHIndex(l, m)] =
				(ct[-m] * c_in[SHIndex(l, m)] +
					-st[-m] * c_in[SHIndex(l, -m)]);
			c_out[SHIndex(l, 0)] = c_in[SHIndex(l, 0)];
			for (int m = 1; m <= l; ++m)
				c_out[SHIndex(l, m)] =
				(ct[m] * c_in[SHIndex(l, m)] +
					st[m] * c_in[SHIndex(l, -m)]);
		}
	}


	void SHConvolveCosTheta(int lmax, const Spectrum *c_in,
		Spectrum *c_out) {
		static const Float c_costheta[18] = { 0.8862268925, 1.0233267546,
			0.4954159260, 0.0000000000, -0.1107783690, 0.0000000000,
			0.0499271341, 0.0000000000, -0.0285469331, 0.0000000000,
			0.0185080823, 0.0000000000, -0.0129818395, 0.0000000000,
			0.0096125342, 0.0000000000, -0.0074057109, 0.0000000000 };
		for (int l = 0; l <= lmax; ++l)
			for (int m = -l; m <= l; ++m) {
				int o = SHIndex(l, m);
				if (l < 18) c_out[o] = lambda(l) * c_in[o] * c_costheta[l];
				else        c_out[o] = 0.f;
			}
	}


	void SHConvolvePhong(int lmax, Float n, const Spectrum *c_in,
		Spectrum *c_out) {
		for (int l = 0; l <= lmax; ++l) {
			Float c_phong = expf(-(l*l) / (2.f * n));
			for (int m = -l; m <= l; ++m) {
				int o = SHIndex(l, m);
				c_out[o] = lambda(l) * c_in[o] * c_phong;
			}
		}
	}



	void SHComputeTransferMatrix(const Point3f &p, Float rayEpsilon,
		const Scene *scene, RNG &rng, int nSamples, int lmax,
		Spectrum *T) {
		for (int i = 0; i < SHTerms(lmax)*SHTerms(lmax); ++i)
			T[i] = 0.f;
		
		Float *Ylm = ALLOCA(Float, SHTerms(lmax));
		for (int i = 0; i < nSamples; ++i) {
			// Compute Monte Carlo estimate of $i$th sample for transfer matrix
			Point2f u(rng.UniformFloat(), rng.UniformFloat());
			Vector3f w = UniformSampleSphere(u);
			Float pdf = UniformSpherePdf();
			if (!scene->IntersectP(Ray(p, w, rayEpsilon))) {
				// Update transfer matrix for unoccluded direction
				SHEvaluate(w, lmax, Ylm);
				for (int j = 0; j < SHTerms(lmax); ++j)
					for (int k = 0; k < SHTerms(lmax); ++k)
						T[j*SHTerms(lmax) + k] += (Ylm[j] * Ylm[k]) / (pdf * nSamples);
			}
		}
	}


	void SHComputeBSDFMatrix(const Spectrum &Kd, const Spectrum &Ks,
		Float roughness, RNG &rng, int nSamples, int lmax, Spectrum *B) {
		
	}


	void SHMatrixVectorMultiply(const Spectrum *M, const Spectrum *v,
		Spectrum *vout, int lmax) {
		for (int i = 0; i < SHTerms(lmax); ++i) {
			vout[i] = 0.f;
			for (int j = 0; j < SHTerms(lmax); ++j)
				vout[i] += M[SHTerms(lmax) * i + j] * v[j];
		}
	}
}