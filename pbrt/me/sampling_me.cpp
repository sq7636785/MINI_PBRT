#include "sampling.h"
#include "geometry.h"
#include "shape.h"


namespace pbrt {
	void StratifiedSample1D(Float *samp, int nSamples, RNG &rng, bool jitter) {
		Float invSampleNum = 1.f / static_cast<Float>(nSamples);
		for (int i = 0; i < nSamples; ++i) {
			Float delta = jitter ? rng.UniformFloat() : 0.5f;
			samp[i] = std::min((i + delta) * invSampleNum, OneMinusEpsilon);
		}
	}

	void StratifiedSample2D(Point2f *samp, int nx, int ny, RNG &rng, bool jitter) {
		Float dx = 1.f / static_cast<Float>(nx);
		Float dy = 1.f / static_cast<Float>(ny);
		for (int i = 0; i < nx; ++i) {
			for (int j = 0; j < ny; ++j) {
				Float deltaX = jitter ? rng.UniformFloat() : 0.5f;
				Float deltaY = jitter ? rng.UniformFloat() : 0.5f;
				samp->x = std::min((i + deltaX) * dx, OneMinusEpsilon);
				samp->y = std::min((j + deltaY) * dy, OneMinusEpsilon);
				++samp;
			}
		}
	}
}