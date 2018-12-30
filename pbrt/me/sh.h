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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_SH_H
#define PBRT_CORE_SH_H

// core/sh.h*
#include "pbrt.h"
#include "geometry.h"
#include "spectrum.h"

namespace pbrt {
	// Spherical Harmonics Declarations
	inline int SHTerms(int lmax) {
		return (lmax + 1) * (lmax + 1);
	}


	inline int SHIndex(int l, int m) {
		return l * l + l + m;
	}


	void SHEvaluate(const Vector3f &v, int lmax, Float *out);
	void SHComputeDiffuseTransfer(const Point3f &p, const Normal3f &n, Float rayEpsilon,
		const Scene *scene, RNG &rng, int nSamples, int lmax, Spectrum *c_transfer);

	
	
	void SHReduceRinging(Spectrum *c, int lmax, Float lambda = .005f);
	void SHRotate(const Spectrum *c_in, Spectrum *c_out, const Matrix4x4 &m,
		int lmax);
	void SHRotateZ(const Spectrum *c_in, Spectrum *c_out, Float alpha, int lmax);
	void SHRotateXMinus(const Spectrum *c_in, Spectrum *c_out, int lmax);
	void SHRotateXPlus(const Spectrum *c_in, Spectrum *c_out, int lmax);
	//void SHSwapYZ(const Spectrum *c_in, Spectrum *c_out, int lmax);
	void SHConvolveCosTheta(int lmax, const Spectrum *c_in, Spectrum *c_out);
	void SHConvolvePhong(int lmax, Float n, const Spectrum *c_in, Spectrum *c_out);
	
	void SHComputeTransferMatrix(const Point3f &p, Float rayEpsilon,
		const Scene *scene, RNG &rng, int nSamples, int lmax, Spectrum *T);
	void SHComputeBSDFMatrix(const Spectrum &Kd, const Spectrum &Ks,
		Float roughness, RNG &rng, int nSamples, int lmax, Spectrum *B);
	void SHMatrixVectorMultiply(const Spectrum *M, const Spectrum *v,
		Spectrum *vout, int lmax);
}
#endif // PBRT_CORE_SH_H