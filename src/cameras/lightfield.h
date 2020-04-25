#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CAMERAS_LIGHTFIELD_H
#define PBRT_CAMERAS_LIGHTFIELD_H

// cameras/lightfield.h*
#include "pbrt.h"
#include "camera.h"
#include "perspective.h"
#include "film.h"

#include <deque>

namespace pbrt {

class LightfieldCamera : public Camera {
  public:
    LightfieldCamera(const AnimatedTransform &CameraToWorld,
                     Float fov,
                     const Vector2i &camerasPerDim,
                     const Bounds2f &gridBounds,
                     Film *film, const Medium *medium);

    Float GenerateRayDifferential(const CameraSample &sample,
                                  RayDifferential *) const;
    Float GenerateRay(const CameraSample &sample, Ray *) const;

  private:
    // Add member variables here for your implementation.
    // This includes any parameters you need in GenerateRayDifferential
    // as well as a deque of Transforms and a deque of PerspectiveCameras
    // if following the assignment's recommended implementation.
};

LightfieldCamera *CreateLightfieldCamera(const ParamSet &params,
                                         const AnimatedTransform &cam2world,
                                         Film *film, const Medium *medium);

}  // namespace pbrt

#endif  // PBRT_CAMERAS_LIGHTFIELD_H
