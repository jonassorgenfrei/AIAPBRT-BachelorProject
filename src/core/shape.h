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

#ifndef PBRT_CORE_SHAPE_H
#define PBRT_CORE_SHAPE_H

 // core/shape.h*
#include "pbrt.h"
#include "geometry.h"
#include "interaction.h"
#include "memory.h"
#include "transform.h"

namespace pbrt {

	// Shape Declarations
	class Shape {
	public:
		// Shape Interface
		// Shapes are defined in object coordinate space

		/// <summary>
		/// Initializes a new instance of the <see cref="Shape"/> class.
		/// </summary>
		/// <param name="ObjectToWorld">Transformation from Object Space to world space</param>
		/// <param name="WorldToObject">Inverse Transformation of <see cred="ObjectToWorld">from world to object space</param>
		/// <param name="reverseOrientation">if set to <c>true</c> the Orientation is reversed/ the surface normal should be reversed.</param>
		Shape(	const Transform *ObjectToWorld,
				const Transform *WorldToObject,
				bool reverseOrientation);
	
		/// <summary>
		/// Finalizes an instance of the <see cref="Shape"/> class.
		/// </summary>
		virtual ~Shape();

		/// <summary>
		/// Returns an (axis-aligned) bounding box in the shape's object space
		/// </summary>
		/// <returns>Returns an (axis-aligned) bounding box in the shape's object space</returns>
		virtual Bounds3f ObjectBound() const = 0;
	
		/// <summary>
		/// Returns an (axis-aligned) bounding box in the world space
		/// </summary>
		/// <returns></returns>
		virtual Bounds3f WorldBound() const;

		// Intersection Test
		// Method that test for ray intersections with their shape
		// To Keep in mind:
		//	=> ignore intersections after Ray::tMax
		//	=> rays are in world space 

		/// <summary>
		/// Returns geometric information about a single ray-shape intersection corresponding to the first intersection
		/// if any in the (0, tMax) parametric range along the ray.
		/// </summary>
		/// <param name="ray">The ray.</param>
		/// <param name="tHit">Pointer to store the parametric distance along the ray if intersection is found. (for multiple intersection, the closest one is returned)</param>
		/// <param name="isect">Pointer to store information about the intersection. (Completly caputres the local geometric pproperties of a surface)</param>
		/// <param name="testAlphaTexture">if set to <c>true</c> the function shoudl perform an alpha cutting operation.</param>
		/// <returns></returns>
		virtual bool Intersect(const Ray &ray,
								Float *tHit,
								SurfaceInteraction *isect,
								bool testAlphaTexture = true) const = 0;

		/// <summary>
		/// Predicate function that determines whether or not an intersection occurs, without returning any details about the intersection itself.
		/// </summary>
		/// <param name="ray">The ray.</param>
		/// <param name="testAlphaTexture">if set to <c>true</c> [test alpha texture].</param>
		/// <returns></returns>
		virtual bool IntersectP(const Ray &ray,
								bool testAlphaTexture = true) const {
			return Intersect(ray, nullptr, nullptr, testAlphaTexture); // default case: just ignores additional information
			// should be implemented more efficient
		}

		/// <summary>
		/// Returns the area of the surface. (e.g. to use as area lights)
		/// </summary>
		/// <returns></returns>
		virtual Float Area() const = 0;

		// Sample a point on the surface of the shape and return the PDF with
		// respect to area on the surface.
		virtual Interaction Sample(const Point2f &u, Float *pdf) const = 0;
		virtual Float Pdf(const Interaction &) const { return 1 / Area(); }

		// Sample a point on the shape given a reference point |ref| and
		// return the PDF with respect to solid angle from |ref|.
		virtual Interaction Sample(const Interaction &ref, const Point2f &u,
			Float *pdf) const;
		virtual Float Pdf(const Interaction &ref, const Vector3f &wi) const;

		// Returns the solid angle subtended by the shape w.r.t. the reference
		// point p, given in world space. Some shapes compute this value in
		// closed-form, while the default implementation uses Monte Carlo
		// integration; the nSamples parameter determines how many samples are
		// used in this case.
		virtual Float SolidAngle(const Point3f &p, int nSamples = 512) const;

		// Shape Public Data
		const Transform *ObjectToWorld, *WorldToObject; // pointer to transformation from obj space to world space and its inverse
		const bool reverseOrientation;	
		const bool transformSwapsHandedness; // stores the return value of the SwapsHandedness() for their object-to-world transformation
	};

}  // namespace pbrt

#endif // PBRT_CORE_SHAPE_H