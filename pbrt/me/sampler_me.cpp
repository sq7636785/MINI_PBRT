#include "sampler.h"
#include "sampling.h"
#include "camera.h"
#include "stats.h"

namespace pbrt {

	//对pixelSampler来说， 首先是对像素内采点， 即， 你要知道一个像素内要采多少个 点， samplerPerPixel
	//然后这里还有一个维度， 可以理解为对不同的地方的采样， 比如， 像素上光线的位置， dof中镜头上的位置，都算是维度
	//初始化是给每个维度初始化像素采样个数那么多数据， 用超了了返回随机采样
	//一个采样点用完之后要开始下一个采样点,将currenPixelIndex++， dimensionOffset置为0

	PixelSampler::PixelSampler(int64_t samplesPerPixel, int nSampledDimensions)
		: Sampler(samplesPerPixel) {
		for (int i = 0; i < nSampledDimensions; ++i) {
			samples1D.push_back(std::vector<Float>(samplesPerPixel));
			samples2D.push_back(std::vector<Point2f>(samplesPerPixel));
		}
	}

	bool PixelSampler::StartNextSample() {
		current1DDimension = current2DDimension = 0;
		return Sampler::StartNextSample();
	}

	bool PixelSampler::SetSampleNumber(int64_t sampleNum) {
		current1DDimension = current2DDimension = 0;
		return Sampler::SetSampleNumber(sampleNum);
	}

	Float PixelSampler::Get1D() {
		ProfilePhase _(Prof::GetSample);
		CHECK_LT(currentPixelSampleIndex, samplesPerPixel);
		if (current1DDimension < samples1D.size()) {
			return samples1D[current1DDimension++][currentPixelSampleIndex];
		} else {
			return rng.UniformFloat();
		}
	}

	Point2f PixelSampler::Get2D() {
		ProfilePhase _(Prof::GetSample);
		CHECK_LT(currentPixelSampleIndex, samplesPerPixel);
		if (current2DDimension < samples2D.size()) {
			return samples2D[current2DDimension++][currentPixelSampleIndex];
		} else {
			return Point2f(rng.UniformFloat(), rng.UniformFloat());
		}
	}

}