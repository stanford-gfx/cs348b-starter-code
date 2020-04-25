// cameras/lightfield.cpp*
#include "cameras/lightfield.h"
#include "paramset.h"
#include "sampler.h"
#include "sampling.h"
#include "light.h"
#include "stats.h"
#include "filters/box.h"
#include <exception>

namespace pbrt {

// A helper function to create new Films for the data cameras in the
// lightfield array. Simply pass in the size in pixels of the film
// and pass this pointer directly to the PerspectiveCamera constructor.
// Note that these films are not actually used to store images for
// the data cameras, but a Camera requires a film for inferring its
// resolution.
static Film *MakeFilm(const Point2i &resolution) {
    std::unique_ptr<Filter> filter(new BoxFilter(Vector2f(0.5f, 0.5f)));

    return new Film(resolution, Bounds2f(Point2f(0, 0), Point2f(1, 1)),
                    std::move(filter), 35, "", 1.f);
}

// Constructor for the LightfieldCamera. This function should create
// any parameters required by GenerateRayDifferential, as well as
// camerasPerDim.x * camerasPerDim.y PBRT PerspectiveCameras that will
// serve as the data cameras in the lightfield array. Each data camera
// will require a Transform as well, which should be stored within the
// LightfieldCamera as described in the assignment handout.
// Sections 2.9.3 and 6.2.2 of the PBRT book describe the
// AnimatedTransform and PerspectiveCamera, respectively, and may be
// helpful.
LightfieldCamera::LightfieldCamera(const AnimatedTransform &CameraToWorld,
                                   Float fov,
                                   const Vector2i &camerasPerDim,
                                   const Bounds2f &gridBounds,
                                   Film *film, const Medium *medium)
    : Camera(CameraToWorld, 0, 0, film, medium)
{
    throw std::runtime_error("Unimplemented");
}

// Generates a new ray corresponding to 'sample' and stores the new ray
// in 'ray'. This function should find the appropriate data camera based
// on the film coordinates of 'sample', and then call the corresponding
// PerspectiveCamera's GenerateRayDifferential function to do the actual
// ray generation.
// Section 6.1 of PBRT describes GenerateRay and GenerateRayDifferential
// in the context of the Camera parent class.
Float LightfieldCamera::GenerateRayDifferential(const CameraSample &sample,
                                                RayDifferential *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);

    throw std::runtime_error("Unimplemented");
}

// This function generates rays that don't require differentials
// simply by delegating to GenerateRayDifferential function that you
// implement and throwing away the differentials. Although this is
// inefficient, it removes the need for a near copy paste between
// GenerateRayDifferentials and GenerateRay.
Float LightfieldCamera::GenerateRay(const CameraSample &sample,
                                    Ray *ray) const {
    RayDifferential rayDiff;
    Float w = GenerateRayDifferential(sample, &rayDiff);

    *ray = rayDiff;

    return w;
}

// API entry point: reads parameters for LightfieldCamera and forwards
// them to constructor. Parameters are the per data camera FOV,
// the number of cameras in the x and y dimensions, and the world space
// bounds of the lightfield camera array.
LightfieldCamera *CreateLightfieldCamera(const ParamSet &params,
                                         const AnimatedTransform &cam2world,
                                         Film *film, const Medium *medium) {
    // Extract common camera parameters from _ParamSet_
    Float fov = params.FindOneFloat("fov", 90.);

    Vector2i camerasPerDim(16, 16);
    int ndims;
    const int *dims = params.FindInt("camerasperdim", &ndims);
    if (dims) {
        if (ndims == 2) {
            camerasPerDim = Vector2i(dims[0], dims[1]);
        } else {
            Error("\"camerasperdim\" should have two values");
        }
    }

    Bounds2f gridBounds(Point2f(-0.6, -0.6), Point2f(0.6, 0.6));
    int nbounds;
    const float *bounds = params.FindFloat("cameragridbounds", &nbounds);
    if (bounds) {
        if (nbounds == 4) {
            gridBounds.pMin = Point2f(bounds[0], bounds[2]);
            gridBounds.pMax = Point2f(bounds[1], bounds[3]);
        } else {
            Error("\"cameragridbounds\" should have four values");
        }
    }

    return new LightfieldCamera(cam2world, fov, camerasPerDim,
                                gridBounds, film, medium);
}

}  // namespace pbrt
