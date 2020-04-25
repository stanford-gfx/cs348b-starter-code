#include <ctype.h>
#include <exception>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <memory>
#include "fileutil.h"
#include "imageio.h"
#include "pbrt.h"
#include "spectrum.h"
#include "parallel.h"
#include "transform.h"
#include "cameras/perspective.h"
#include "filters/box.h"
#include "samplers/halton.h"

#include "sampling.h"

using namespace pbrt;
using namespace std;

// Configuration parameters for LightfieldRenderer.
// These default values are overwritten in main() by
// the provided command line arguments.
struct LightfieldParameters {
    Float inputFov {90.0f};
    Float outputFov {45.0f};
    Vector2i camerasPerDim {1, 1};
    Bounds2f griddim {Point2f {0, 0}};
    Point2i outputDim {256, 256};
    Point3f cameraPos {0, 0, -0.02f};
    Vector3f cameraRot {0, 0, 0};
    Float focalDistance {10.0f};
    Float lensRadius {0.0f};
    int samplesPerPixel {1};
};

// Constructs a film that should be passed in when constructing PBRT cameras
Film * MakeFilm(const Point2i &resolution, const string &filename) {
    auto filter = unique_ptr<BoxFilter>(new BoxFilter(Vector2f(0.5f, 0.5f)));

    return new Film(resolution,
                    Bounds2f(Point2f(0.f, 0.f), Point2f(1.f, 1.f)),
                    move(filter), 35, filename, 1);
}

// Helper class to read in the stored light field image and query
// the radiance of rays within the light field.
class LightfieldManager {
public:
    LightfieldManager(const LightfieldParameters &params,
                      const string &inputFilename)
        : lightfieldResolution(),
          lightfieldImage(ReadImage(inputFilename, &lightfieldResolution)),
          params(params)
          // Initialize any new member variables you add here
    {
        if (!lightfieldImage) {
            cerr << "Failed to read lightfield image" << endl;
            exit(EXIT_FAILURE);
        }

        // More complex initialization / construction of member variables
        // can happen here.
    }

    PerspectiveCamera MakeVirtualCamera(const string &outputFilename) const;

    bool ComputeRadiance(const Ray &virtualRay, RGBSpectrum *radiance) const;

private:
    // Dimensions of the input lightfield in pixels
    Point2i lightfieldResolution;
    // The input lightfield stored as a row major image
    unique_ptr<RGBSpectrum[]> lightfieldImage;
    // Copy of lightfield parameters
    const LightfieldParameters params;

    // Add any parameters or transforms here as member variables for use by
    // your MakeVirtualCamera and ComputeRadiance implementations
};

// Construct a PBRT PerspectiveCamera that will provide the functionality
// for the virtual camera that launches rays through the lightfield.
// As mentioned previously, sections 2.9.3 and 6.2.2 of the PBRT book
// describe the AnimatedTransform and PerspectiveCamera, respectively,
// and may be helpful.
PerspectiveCamera LightfieldManager::MakeVirtualCamera(const string &outputFilename) const {
    // 3 - 6 lines
    throw runtime_error("Unimplemented");
}

// Compute the radiance along 'virtualRay'. Store the result in *radiance,
// and return true if 'virtualRay' is within the bounds of the lightfield.
// Otherwise, return false.
// A short description of screen space and raster space can be found at the
// beginning of chapter 6.2 in the PBRT book which may be helpful.
bool LightfieldManager::ComputeRadiance(const Ray &virtualRay,
                                        RGBSpectrum *radiance) const {
    // 50 - 150 lines
    throw runtime_error("Unimplemented");
}

void Render(const LightfieldParameters &params,
            const string &inputFilename,
            const string &outputFilename) {

    // Construct the LightfieldManager, which holds the input lightfield
    // and all transforms and parameters necessary to sample values from it.
    LightfieldManager lightfieldManager(params, inputFilename);

    // Make a convenience bounds variable for getting the film tile and
    // iterating over pixel indices in the output image later.
    const Bounds2i imageBounds(Point2i(0, 0), params.outputDim);

    // Create the random sampler used by the virtual
    // camera when generating rays.
    HaltonSampler sampler(params.samplesPerPixel, imageBounds);

    // Construct the virtual camera (represented by a PBRT PerspectiveCamera).
    PerspectiveCamera virtualCamera =
        lightfieldManager.MakeVirtualCamera(outputFilename);

    // Output pixel values will be written into this FilmTile that is merged
    // into the film at the end of the function.
    unique_ptr<FilmTile> outputTile =
        virtualCamera.film->GetFilmTile(imageBounds);

    // Loop over pixels in the output image, generating rays from the virtual
    // camera and accumulating the resulting radiances into the output film.
    for (Point2i pixel : imageBounds) {
        sampler.StartPixel(pixel);

        do {
            CameraSample cameraSample = sampler.GetCameraSample(pixel);

            Ray virtualRay;
            virtualCamera.GenerateRay(cameraSample, &virtualRay);

            RGBSpectrum radiance;
            if (!lightfieldManager.ComputeRadiance(virtualRay, &radiance)) {
                continue;
            }

            outputTile->AddSample(cameraSample.pFilm, radiance, 1.f);
        } while (sampler.StartNextSample());
    }

    // Write the image
    virtualCamera.film->MergeFilmTile(move(outputTile));
    virtualCamera.film->WriteImage();
}

static void usage(const char *name) {
    cerr << name << " [options] input_lightfield.exr output_file.exr" << endl <<
        R"(Options:
    --camsperdim <x> <y>                    Number of data cameras in the X and
                                            Y direction, respectively
    --camerapos <x> <y> <z>                 Position of the desired camera for
                                            reconstruction (treating the data
                                            camera plane as z = 0)
    --camerarot <x> <y> <z>                 Rotation of the desired around the
                                            x, y, and z axes in degrees. Starting
                                            orientation matches the orientation of
                                            the data cameras.
    --focaldistance <fd>                    Focal distance of the virtual camera:
                                            distance along the camera's z axis to the
                                            focal plane
    --griddim <xmin> <xmax> <ymin> <ymax>   The bounding box of the data cameras
                                            positions in the x and y dimensions
                                            (in meters)
    --inputfov <f>                          FOV (in degrees) of the input
                                            lightfield data cameras
    --outputfov <f>                         FOV (in degrees) of the output image
    --outputdim <w> <h>                     Output size in pixels per dimension
    --lensradius <r>                        The radius of the thin lens in the
                                            thin lens approximation for the output
                                            camera.
    --samplesperpixel <p>                   The number of samples to take (per pixel)
                                            in Monte Carlo the integration of the
                                            final image.
)";
    exit(EXIT_FAILURE);
}

int main(int argc, const char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = 1; // Warning and above

    LightfieldParameters params;

    int i;
    for (i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') break;
        if (!strcmp(argv[i], "--inputfov")) {
            ++i;
            params.inputFov = atof(argv[i]);
        } else if (!strcmp(argv[i], "--outputfov")) {
            ++i;
            params.outputFov = atof(argv[i]);
        } else if (!strcmp(argv[i], "--camsperdim")) {
            ++i;
            params.camerasPerDim.x = atoi(argv[i]);
            ++i;
            params.camerasPerDim.y = atoi(argv[i]);
        } else if (!strcmp(argv[i], "--outputdim")) {
            ++i;
            params.outputDim.x = atoi(argv[i]);
            ++i;
            params.outputDim.y = atoi(argv[i]);
        } else if (!strcmp(argv[i], "--camerapos")) {
            ++i;
            params.cameraPos.x = atof(argv[i]);
            ++i;
            params.cameraPos.y = atof(argv[i]);
            ++i;
            params.cameraPos.z = atof(argv[i]);
        } else if (!strcmp(argv[i], "--camerarot")) {
            ++i;
            params.cameraRot.x = atof(argv[i]);
            ++i;
            params.cameraRot.y = atof(argv[i]);
            ++i;
            params.cameraRot.z = atof(argv[i]);
        } else if (!strcmp(argv[i], "--griddim")) {
            ++i;
            params.griddim.pMin.x = atof(argv[i]);
            ++i;
            params.griddim.pMax.x = atof(argv[i]);
            ++i;
            params.griddim.pMin.y = atof(argv[i]);
            ++i;
            params.griddim.pMax.y = atof(argv[i]);
        }  else if (!strcmp(argv[i], "--focaldistance")) {
            ++i;
            params.focalDistance = atof(argv[i]);
        } else if (!strcmp(argv[i], "--lensradius")) {
            ++i;
            params.lensRadius = atof(argv[i]);
        }  else if (!strcmp(argv[i], "--samplesperpixel")) {
            ++i;
            params.samplesPerPixel = atoi(argv[i]);
        } else {
            usage(argv[0]);
        }
    }

    if (i + 1 >= argc)
        usage(argv[0]);

    const char *inFilename = argv[i], *outFilename = argv[i + 1];

    Render(params, inFilename, outFilename);

    return 0;
}
