
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_CAMERA_H
#define PBRT_CORE_CAMERA_H

// core/camera.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"
#include "film.h"

namespace pbrt {

/// <summary>
/// Camera Declarations holds camera options and defines the interface that all camera
/// implementations must provide
/// </summary>
class Camera {
  public:
    // Camera Interface

	/// <summary>
	/// Initializes a new instance of the <see cref="Camera"/> class.
	/// </summary>
	/// <param name="CameraToWorld">The transformation that places the camera in the scene over time.</param>
	/// <param name="shutterOpen">The shutter open time.</param>
	/// <param name="shutterClose">The shutter close time.</param>
	/// <param name="film">Pointer to a film object, that represents the final image.</param>
	/// <param name="medium">The scattering medium that the camera lies in.</param>
	Camera(const AnimatedTransform& CameraToWorld, Float shutterOpen,
           Float shutterClose, Film *film, const Medium *medium);

	/// <summary>
	/// Finalizes an instance of the <see cref="Camera"/> class.
	/// </summary>
	virtual ~Camera();

	/// <summary>
	/// Generates a world space ray corresponding to a sample position on the film plane.
	/// </summary>
	/// <param name="sample">The sample.</param>
	/// <param name="ray">The ray.</param>
	/// <returns></returns>
	virtual Float GenerateRay(const CameraSample& sample, Ray* ray) const = 0;

	/// <summary>
	/// Geerate a world space ray corresponding to a sample position on the filme plane and 
	/// computes information about the image area (shift xy directions) that the ray is sampling. (e.g. information for anti-aliasing computations)
	/// </summary>
	/// <param name="sample">The sample.</param>
	/// <param name="rd">The rd.</param>
	/// <returns></returns>
	virtual Float GenerateRayDifferential(const CameraSample& sample,
                                          RayDifferential *rd) const;

	/// <summary>
	/// Wes the specified ray.
	/// </summary>
	/// <param name="ray">The ray.</param>
	/// <param name="pRaster2">The p raster2.</param>
	/// <returns></returns>
	virtual Spectrum We(const Ray& ray, Point2f* pRaster2 = nullptr) const;
	
	/// <summary>
	/// PDFs the we.
	/// </summary>
	/// <param name="ray">The ray.</param>
	/// <param name="pdfPos">The PDF position.</param>
	/// <param name="pdfDir">The PDF dir.</param>
	virtual void Pdf_We(const Ray& ray, Float* pdfPos, Float* pdfDir) const;

	/// <summary>
	/// Samples the wi.
	/// </summary>
	/// <param name="ref">The reference.</param>
	/// <param name="u">The u.</param>
	/// <param name="wi">The wi.</param>
	/// <param name="pdf">The PDF.</param>
	/// <param name="pRaster">The p raster.</param>
	/// <param name="vis">The vis.</param>
	/// <returns></returns>
	virtual Spectrum Sample_Wi(const Interaction& ref, const Point2f& u,
                               Vector3f *wi, Float *pdf, Point2f *pRaster,
                               VisibilityTester *vis) const;

    // Camera Public Data
    AnimatedTransform CameraToWorld;
	const Float shutterOpen, shutterClose;
    Film *film;
    const Medium *medium;
};

/// <summary>
/// Holds all of the sample values needed to specify a camera ray. 
/// </summary>
struct CameraSample {
    Point2f pFilm;	/// point on the film to which the generated ray carries radiance
    Point2f pLens;	/// point on the lens the ray passes through
    Float time;		/// the time at which the ray should sample the scene (use value to linearly interpolate whithing the shutter time range)
};

inline std::ostream &operator<<(std::ostream &os, const CameraSample &cs) {
    os << "[ pFilm: " << cs.pFilm << " , pLens: " << cs.pLens <<
        StringPrintf(", time %f ]", cs.time);
    return os;
}

class ProjectiveCamera : public Camera {
  public:
    // ProjectiveCamera Public Methods
    ProjectiveCamera(const AnimatedTransform &CameraToWorld,
                     const Transform &CameraToScreen,
                     const Bounds2f &screenWindow, Float shutterOpen,
                     Float shutterClose, Float lensr, Float focald, Film *film,
                     const Medium *medium)
        : Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
          CameraToScreen(CameraToScreen) {
        // Initialize depth of field parameters
        lensRadius = lensr;
        focalDistance = focald;

        // Compute projective camera transformations

        // Compute projective camera screen transformations
        ScreenToRaster =
            Scale(film->fullResolution.x, film->fullResolution.y, 1) *
            Scale(1 / (screenWindow.pMax.x - screenWindow.pMin.x),
                  1 / (screenWindow.pMin.y - screenWindow.pMax.y), 1) *
            Translate(Vector3f(-screenWindow.pMin.x, -screenWindow.pMax.y, 0));
        RasterToScreen = Inverse(ScreenToRaster);
        RasterToCamera = Inverse(CameraToScreen) * RasterToScreen;
    }

  protected:
    // ProjectiveCamera Protected Data
    Transform CameraToScreen, RasterToCamera;
    Transform ScreenToRaster, RasterToScreen;
    Float lensRadius, focalDistance;
};

}  // namespace pbrt

#endif  // PBRT_CORE_CAMERA_H
