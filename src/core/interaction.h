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

#ifndef PBRT_CORE_INTERACTION_H
#define PBRT_CORE_INTERACTION_H

 // core/interaction.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"
//#include "medium.h"
//#include "material.h"

namespace pbrt {
	// Interaction Declarations
	// Generic class, provides some common member variables and methods
	struct Interaction {
		// Interaction Public Methods
		Interaction() : time(0) {}
		Interaction(const Point3f &p, const Normal3f &n, const Vector3f &pError,
					const Vector3f &wo, Float time,
					const MediumInterface &mediumInterface)
			: p(p),
			time(time),
			pError(pError),
			wo(Normalize(wo)),
			n(n),
			mediumInterface(mediumInterface) {}

		/// <summary>
		/// Determines whether [is surface interaction].
		/// </summary>
		/// <returns>
		///   <c>true</c> if [is surface interaction]; otherwise, <c>false</c>.
		/// </returns>
		bool IsSurfaceInteraction() const { return n != Normal3f(); }
		/// <summary>
		/// Generate rays leaving intersection points
		/// </summary>
		/// <param name="d">The d.</param>
		/// <returns></returns>
		Ray SpawnRay(const Vector3f& d) const {
			Point3f o = OffsetRayOrigin(p, pError, n, d); 
			return Ray(o, d, Infinity, time, GetMedium(d));
		}
		/// <summary>
		/// 
		/// </summary>
		/// <param name="p2">Point</param>
		/// <returns></returns>
		Ray SpawnRayTo(const Point3f& p2) const {
			Point3f origin = OffsetRayOrigin(p, pError, n, p2 - p);
			Vector3f d = p2 - p;
			// set the tMax value of shadow rays to be just under one
			// so they stop before the surface of light sources
			return Ray(origin, d, 1 - ShadowEpsilon, time, GetMedium(d));
		}

		/// <summary>
		/// Spawns the ray to.
		/// </summary>
		/// <param name="it">Interaction</param>
		/// <returns></returns>
		Ray SpawnRayTo(const Interaction& it) const {
			Point3f origin = OffsetRayOrigin(p, pError, n, it.p - p);
			Point3f target = OffsetRayOrigin(it.p, it.pError, it.n, origin - it.p);
			Vector3f d = target - origin;
			return Ray(origin, d, 1 - ShadowEpsilon, time, GetMedium(d));
		}
		Interaction(const Point3f &p, const Vector3f &wo, Float time,
			const MediumInterface &mediumInterface)
			: p(p), time(time), wo(wo), mediumInterface(mediumInterface) {}
		Interaction(const Point3f &p, Float time,
			const MediumInterface &mediumInterface)
			: p(p), time(time), mediumInterface(mediumInterface) {}
		bool IsMediumInteraction() const { return !IsSurfaceInteraction(); }
		const Medium *GetMedium(const Vector3f &w) const {
			return Dot(w, n) > 0 ? mediumInterface.outside : mediumInterface.inside;
		}
		const Medium *GetMedium() const {
			CHECK_EQ(mediumInterface.inside, mediumInterface.outside);
			return mediumInterface.inside;
		}

		// Interaction Public Data
		// all interactions must hava a point p and time associated with them
		Point3f p;
		Float time;
		// gives an conservative bound on the floating-point error,
		// where the point p was computed by ray intersection 
		Vector3f pError;
		// stores the negative ray direction (outgoing direction when computing lighting at points)
		Vector3f wo;
		// Surface normal at the point
		Normal3f n;
		// scattering media at their point
		MediumInterface mediumInterface;
	};
	
	class MediumInteraction : public Interaction {
	public:
		// MediumInteraction Public Methods
		MediumInteraction() : phase(nullptr) {}
		MediumInteraction(const Point3f &p, const Vector3f &wo, Float time,
			const Medium *medium, const PhaseFunction *phase)
			: Interaction(p, wo, time, medium), phase(phase) {}
		bool IsValid() const { return phase != nullptr; }

		// MediumInteraction Public Data
		const PhaseFunction *phase;
	};

	// SurfaceInteraction Declarations
	// Represents local information at a point on a 2D surface.
	// Serves to cleanly isolate the geometric portion of the ray tracer from the shading and illumination portions
	class SurfaceInteraction : public Interaction {
		public:
			// SurfaceInteraction Public Methods
			SurfaceInteraction() {}
			SurfaceInteraction(	const Point3f &p, const Vector3f &pError,
								const Point2f &uv, const Vector3f &wo,
								const Vector3f &dpdu, const Vector3f &dpdv,
								const Normal3f &dndu, const Normal3f &dndv, Float time,
								const Shape *sh, int faceIndex = 0);

			/// <summary>
			/// Updates the shading geometry.
			/// </summary>
			/// <param name="dpdu">The dpdu.</param>
			/// <param name="dpdv">The DPDV.</param>
			/// <param name="dndu">The dndu.</param>
			/// <param name="dndv">The DNDV.</param>
			/// <param name="orientationIsAuthoritative">if set to <c>true</c> [orientation is authoritative].</param>
			void SetShadingGeometry(const Vector3f &dpdu, const Vector3f &dpdv,
									const Normal3f &dndu, const Normal3f &dndv,
									bool orientationIsAuthoritative);

			void ComputeScatteringFunctions(const RayDifferential &ray, MemoryArena &arena,
											bool allowMultipleLobes = false,
											TransportMode mode = TransportMode::Radiance);

			void ComputeDifferentials(const RayDifferential &r) const;

			Spectrum Le(const Vector3f &w) const;

			// SurfaceInteraction Public Data
			Point2f uv;	// coordinates from the parameterization of the surface
			Vector3f dpdu, dpdv;	// parametric derivatives of the point dp/du and dp/dv
			Normal3f dndu, dndv;	// partial derivatives of the surface normal
									// show the differential change in surface normal as we 
									// move u and v along the surface
			const Shape *shape = nullptr;	// pointer to the shape that the point lies on

			// Stores second instance of a surface normal and the various partial derivatives to
			// represent possibly perturbed values of these quantities as can be generated by bump
			// mapping or interpolated per-vertex normals with triangles.
			struct {
				Normal3f n;
				Vector3f dpdu, dpdv;
				Normal3f dndu, dndv;
			} shading; // some parts of the system use this shading geometry, 
			// while others need to work with the original quantities

			// Pointer to the Primitive that the ray hit
			const Primitive *primitive = nullptr;

			// BSDF and BSSRDF pointers for the point are stored in the SurfaceInteraction
			BSDF *bsdf = nullptr;
			BSSRDF *bssrdf = nullptr;

			mutable Vector3f dpdx, dpdy;
			mutable Float dudx = 0, dvdx = 0, dudy = 0, dvdy = 0;

			// Added after book publication. Shapes can optionally provide a face
			// index with an intersection point for use in Ptex texture lookups.
			// If Ptex isn't being used, then this value is ignored.
			int faceIndex = 0;
	};

}  // namespace pbrt

#endif // PBRT_CORE_INTERACTION_H