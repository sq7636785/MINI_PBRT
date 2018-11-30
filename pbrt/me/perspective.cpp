
// cameras/perspective.cpp*
#include "me/perspective.h"
#include "paramset.h"
#include "sampler.h"
#include "sampling.h"
#include "light.h"
#include "stats.h"

namespace pbrt {

	// PerspectiveCamera Method Definitions
	PerspectiveCamera::PerspectiveCamera(const AnimatedTransform &CameraToWorld, 
		                                 const Bounds2f &screenWindow, 
		                                 Float shutterOpen, 
		                                 Float shutterClose, 
		                                 Float lensRadius, 
		                                 Float focalDistance, 
		                                 Float fov, 
		                                 Film *film, 
		                                 const Medium *medium) 
		: ProjectiveCamera(CameraToWorld, Perspective(fov, 1e-2f, 1000.f), 
		                  screenWindow, shutterOpen, shutterClose, lensRadius, focalDistance, film, medium){
		dxCamera = (RasterToCamera(Point3f(1, 0, 0)) - RasterToCamera(Point3f(0, 0, 0)));
		dyCamera = (RasterToCamera(Point3f(0, 1, 0)) - RasterToCamera(Point3f(0, 0, 0)));

		Point2i res = film->fullResolution;
		Point3f pMin = RasterToCamera(Point3f(0, 0, 0));
		Point3f pMax = RasterToCamera(Point3f(res.x, res.y, 0));
		pMin /= pMin.z;
		pMax /= pMax.z;
		A = std::abs((pMax.x - pMin.x) * (pMax.y - pMin.y));
	}

	Float PerspectiveCamera::GenerateRay(const CameraSample &sample, Ray *ray) const{
		ProfilePhase prof(Prof::GenerateCameraRay);
		//compute raster and camera sample positions
		Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
		Point3f pCamera = RasterToCamera(pFilm);
		*ray = Ray(Point3f(0, 0, 0), Normalize(Vector3f(pCamera)));
		//DOF
		if (lensRadius > 0) {
			Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);
			//compute point on plane of focus
			Float ft = focalDistance / ray->d.z;
			Point3f pFocus = (*ray)(ft);
			
			//update ray
			ray->o = Point3f(pLens.x, pLens.y, 0);
			ray->d = Normalize(pFocus - ray->o);
		}

		ray->time = Lerp(sample.time, shutterOpen, shutterClose);
		ray->medium = medium;
		//camera to world
		*ray = (CameraToWorld)(*ray);
		return 1;
	}


Float PerspectiveCamera::GenerateRayDifferential(const CameraSample &sample,
                                                 RayDifferential *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    // Compute raster and camera sample positions
    Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
    Point3f pCamera = RasterToCamera(pFilm);
    Vector3f dir = Normalize(Vector3f(pCamera.x, pCamera.y, pCamera.z));
    *ray = RayDifferential(Point3f(0, 0, 0), dir);
    // Modify ray for depth of field
    if (lensRadius > 0) {
        // Sample point on lens
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);

        // Compute point on plane of focus
        Float ft = focalDistance / ray->d.z;
        Point3f pFocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point3f(pLens.x, pLens.y, 0);
        ray->d = Normalize(pFocus - ray->o);
    }

    // Compute offset rays for _PerspectiveCamera_ ray differentials
    if (lensRadius > 0) {
        // Compute _PerspectiveCamera_ ray differentials accounting for lens

        // Sample point on lens
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);
        Vector3f dx = Normalize(Vector3f(pCamera + dxCamera));
        Float ft = focalDistance / dx.z;
        Point3f pFocus = Point3f(0, 0, 0) + (ft * dx);
        ray->rxOrigin = Point3f(pLens.x, pLens.y, 0);
        ray->rxDirection = Normalize(pFocus - ray->rxOrigin);

        Vector3f dy = Normalize(Vector3f(pCamera + dyCamera));
        ft = focalDistance / dy.z;
        pFocus = Point3f(0, 0, 0) + (ft * dy);
        ray->ryOrigin = Point3f(pLens.x, pLens.y, 0);
        ray->ryDirection = Normalize(pFocus - ray->ryOrigin);
    } else {
        ray->rxOrigin = ray->ryOrigin = ray->o;
        ray->rxDirection = Normalize(Vector3f(pCamera) + dxCamera);
        ray->ryDirection = Normalize(Vector3f(pCamera) + dyCamera);
    }
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->medium = medium;
    *ray = CameraToWorld(*ray);
    ray->hasDifferentials = true;
    return 1;
}

Spectrum PerspectiveCamera::We(const Ray &ray, Point2f *pRaster2) const {
    // Interpolate camera matrix and check if $\w{}$ is forward-facing
    Transform c2w;
    CameraToWorld.Interpolate(ray.time, &c2w);
    Float cosTheta = Dot(ray.d, c2w(Vector3f(0, 0, 1)));
    if (cosTheta <= 0) return 0;

    // Map ray $(\p{}, \w{})$ onto the raster grid
    Point3f pFocus = ray((lensRadius > 0 ? focalDistance : 1) / cosTheta);
    Point3f pRaster = Inverse(RasterToCamera)(Inverse(c2w)(pFocus));

    // Return raster position if requested
    if (pRaster2) *pRaster2 = Point2f(pRaster.x, pRaster.y);

    // Return zero importance for out of bounds points
    Bounds2i sampleBounds = film->GetSampleBounds();
    if (pRaster.x < sampleBounds.pMin.x || pRaster.x >= sampleBounds.pMax.x ||
        pRaster.y < sampleBounds.pMin.y || pRaster.y >= sampleBounds.pMax.y)
        return 0;

    // Compute lens area of perspective camera
    Float lensArea = lensRadius != 0 ? (Pi * lensRadius * lensRadius) : 1;

    // Return importance for point on image plane
    Float cos2Theta = cosTheta * cosTheta;
    return Spectrum(1 / (A * lensArea * cos2Theta * cos2Theta));
}

void PerspectiveCamera::Pdf_We(const Ray &ray, Float *pdfPos,
                               Float *pdfDir) const {
    // Interpolate camera matrix and fail if $\w{}$ is not forward-facing
    Transform c2w;
    CameraToWorld.Interpolate(ray.time, &c2w);
    Float cosTheta = Dot(ray.d, c2w(Vector3f(0, 0, 1)));
    if (cosTheta <= 0) {
        *pdfPos = *pdfDir = 0;
        return;
    }

    // Map ray $(\p{}, \w{})$ onto the raster grid
    Point3f pFocus = ray((lensRadius > 0 ? focalDistance : 1) / cosTheta);
    Point3f pRaster = Inverse(RasterToCamera)(Inverse(c2w)(pFocus));

    // Return zero probability for out of bounds points
    Bounds2i sampleBounds = film->GetSampleBounds();
    if (pRaster.x < sampleBounds.pMin.x || pRaster.x >= sampleBounds.pMax.x ||
        pRaster.y < sampleBounds.pMin.y || pRaster.y >= sampleBounds.pMax.y) {
        *pdfPos = *pdfDir = 0;
        return;
    }

    // Compute lens area of perspective camera
    Float lensArea = lensRadius != 0 ? (Pi * lensRadius * lensRadius) : 1;
    *pdfPos = 1 / lensArea;
    *pdfDir = 1 / (A * cosTheta * cosTheta * cosTheta);
}

Spectrum PerspectiveCamera::Sample_Wi(const Interaction &ref, const Point2f &u,
                                      Vector3f *wi, Float *pdf,
                                      Point2f *pRaster,
                                      VisibilityTester *vis) const {
    // Uniformly sample a lens interaction _lensIntr_
    Point2f pLens = lensRadius * ConcentricSampleDisk(u);
    Point3f pLensWorld = CameraToWorld(ref.time, Point3f(pLens.x, pLens.y, 0));
    Interaction lensIntr(pLensWorld, ref.time, medium);
    lensIntr.n = Normal3f(CameraToWorld(ref.time, Vector3f(0, 0, 1)));

    // Populate arguments and compute the importance value
    *vis = VisibilityTester(ref, lensIntr);
    *wi = lensIntr.p - ref.p;
    Float dist = wi->Length();
    *wi /= dist;

    // Compute PDF for importance arriving at _ref_

    // Compute lens area of perspective camera
    Float lensArea = lensRadius != 0 ? (Pi * lensRadius * lensRadius) : 1;
    *pdf = (dist * dist) / (AbsDot(lensIntr.n, *wi) * lensArea);
    return We(lensIntr.SpawnRay(-*wi), pRaster);
}

PerspectiveCamera *CreatePerspectiveCamera(const ParamSet &params,
                                           const AnimatedTransform &cam2world,
                                           Film *film, const Medium *medium) {
    // Extract common camera parameters from _ParamSet_
    Float shutterOpen = params.FindOneFloat("shutteropen", 0.0f);
	Float shutterClose = params.FindOneFloat("shutterclose", 1.0f);
	if (shutterOpen < shutterClose) {
		std::swap(shutterOpen, shutterClose);
	}

	Float lensRadius = params.FindOneFloat("lensradius", 0.0f);
	Float focalDistance = params.FindOneFloat("focaldistance", 1e6);
	Float frame = params.FindOneFloat("frameaspectratio", Float(film->fullResolution.x) / Float(film->fullResolution.y));
	Bounds2f screen;
	if (frame > 1.f) {
		screen.pMin.x = -frame;
		screen.pMax.x = frame;
		screen.pMin.y = -1.f;
		screen.pMax.y = 1.f;
	} else {
		screen.pMin.x = -1.f;
		screen.pMax.x = 1.f;
		screen.pMin.y = -1.f / frame;
		screen.pMax.y = 1.f / frame;
	}

	int swi;
	const Float *sw = params.FindFloat("screenwindow", &swi);
	if (sw) {
		if (swi == 4) {
			screen.pMin.x = sw[0];
			screen.pMax.x = sw[1];
			screen.pMin.y = sw[2];
			screen.pMax.y = sw[3];
		} else
			Error("\"screenwindow\" should have four values");
	}

	Float fov = params.FindOneFloat("fov", 90.0f);
	Float halfFov = params.FindOneFloat("halffov", -1.0f);
	if (halfFov > 0.f) {
		fov = 2 * halfFov;
	}

	return new PerspectiveCamera(cam2world, screen, shutterOpen, shutterClose, lensRadius, focalDistance, fov, film, medium);

}

}  // namespace pbrt