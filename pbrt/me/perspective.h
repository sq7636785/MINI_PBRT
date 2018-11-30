#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CAMERAS_PERSPECTIVE_H
#define PBRT_CAMERAS_PERSPECTIVE_H

// cameras/perspective.h*
#include "pbrt.h"
#include "camera.h"
#include "film.h"

namespace pbrt {

// PerspectiveCamera Declarations
class PerspectiveCamera : public ProjectiveCamera {
  public:
    // PerspectiveCamera Public Methods
    PerspectiveCamera(const AnimatedTransform &CameraToWorld,
                      const Bounds2f &screenWindow, Float shutterOpen,
                      Float shutterClose, Float lensRadius, Float focalDistance,
                      Float fov, Film *film, const Medium *medium);
    Float GenerateRay(const CameraSample &sample, Ray *) const;
    Float GenerateRayDifferential(const CameraSample &sample,
                                  RayDifferential *ray) const;
    Spectrum We(const Ray &ray, Point2f *pRaster2 = nullptr) const;
    void Pdf_We(const Ray &ray, Float *pdfPos, Float *pdfDir) const;
    Spectrum Sample_Wi(const Interaction &ref, const Point2f &sample,
                       Vector3f *wi, Float *pdf, Point2f *pRaster,
                       VisibilityTester *vis) const;

  private:
    // PerspectiveCamera Private Data
    Vector3f dxCamera, dyCamera;
    Float A;
};

PerspectiveCamera *CreatePerspectiveCamera(const ParamSet &params,
                                           const AnimatedTransform &cam2world,
                                           Film *film, const Medium *medium);

}  // namespace pbrt

#endif  // PBRT_CAMERAS_PERSPECTIVE_H
