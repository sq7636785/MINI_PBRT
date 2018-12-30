
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


// samplers/stratified.cpp*
#include "stratified.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {

// StratifiedSampler Method Definitions
	void StratifiedSampler::StartPixel(const Point2i& p) {
		//generate single stratified samples for the pixel
		for (size_t i = 0; i < samples1D.size(); ++i) {
			StratifiedSample1D(&samples1D[i][0], xPixelSamples * yPixelSamples, rng, jitterSamples);
			Shuffle(&samples1D[i][0], xPixelSamples * yPixelSamples, 1, rng);
		}
		for (size_t i = 0; i < samples2D.size(); ++i) {
			StratifiedSample2D(&samples2D[i][0], xPixelSamples, yPixelSamples, rng, jitterSamples);
			Shuffle(&samples2D[i][0], xPixelSamples * yPixelSamples, 1, rng);
		}

		//generate array stratified samples for pixel
		//这里是申请一段随机数， 在光线跟踪中用的， 在渲染每一个摄像机sample时， 可能会需要一系列个sample去追踪光线。
		for (size_t i = 0; i < samples1DArraySizes.size(); ++i) {
			for (int64_t j = 0; j < samplesPerPixel; ++j) {
				int count = samples1DArraySizes[i];
				StratifiedSample1D(&sampleArray1D[i][j * count], count, rng, jitterSamples);
				Shuffle(&sampleArray1D[i][j*count], count, 1, rng);
			}
		}

		for (size_t i = 0; i < samples2DArraySizes.size(); ++i) {
			for (int64_t j = 0; j < samplesPerPixel; ++j) {
				int count = samples2DArraySizes[i];
				LatinHypercube(&sampleArray2D[i][j * count].x, count, 2, rng);
			}
		}
		Sampler::StartPixel(p);
	}


std::unique_ptr<Sampler> StratifiedSampler::Clone(int seed) {
	StratifiedSampler* ns = new StratifiedSampler(*this);
	ns->rng.SetSequence(seed);
	return std::unique_ptr<Sampler>(ns);
}

StratifiedSampler *CreateStratifiedSampler(const ParamSet &params) {
    bool jitter = params.FindOneBool("jitter", true);
    //int xsamp = params.FindOneInt("xsamples", 4);
	//int ysamp = params.FindOneInt("ysamples", 4);
	int nsamp = params.FindOneInt("pixelsamples", 16);
	int xsamp = static_cast<int>(std::sqrt(nsamp));
	int ysamp = xsamp;
    int sd = params.FindOneInt("dimensions", 4);
    if (PbrtOptions.quickRender) xsamp = ysamp = 1;
    return new StratifiedSampler(xsamp, ysamp, jitter, sd);
}

}  // namespace pbrt
