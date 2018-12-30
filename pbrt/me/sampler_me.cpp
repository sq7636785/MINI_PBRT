#include "sampler.h"
#include "sampling.h"
#include "camera.h"
#include "stats.h"

namespace pbrt {

	//��pixelSampler��˵�� �����Ƕ������ڲɵ㣬 ���� ��Ҫ֪��һ��������Ҫ�ɶ��ٸ� �㣬 samplerPerPixel
	//Ȼ�����ﻹ��һ��ά�ȣ� �������Ϊ�Բ�ͬ�ĵط��Ĳ����� ���磬 �����Ϲ��ߵ�λ�ã� dof�о�ͷ�ϵ�λ�ã�������ά��
	//��ʼ���Ǹ�ÿ��ά�ȳ�ʼ�����ز���������ô�����ݣ� �ó����˷����������
	//һ������������֮��Ҫ��ʼ��һ��������,��currenPixelIndex++�� dimensionOffset��Ϊ0

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