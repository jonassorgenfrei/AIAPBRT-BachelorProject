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

#ifndef PBRT_SHAPES_CYLINDER_H
#define PBRT_SHAPES_CYLINDER_H

 // shapes/cylinder.h*
#include "shape.h"

namespace pbrt {

	// Cylinder Declarations
	class Cylinder : public Shape {
	public:
		// Cylinder Public Methods
		/// <summary>
		/// Initializes a new instance of the <see cref="Cylinder"/> class.
		/// </summary>
		/// <param name="ObjectToWorld">The object to world transformation.</param>
		/// <param name="WorldToObject">The world to object transformation. (Reverse object to world transformation)</param>
		/// <param name="reverseOrientation">if set to <c>true</c> [reverse orientation].</param>
		/// <param name="radius">The radius.</param>
		/// <param name="zMin">The z minimum.</param>
		/// <param name="zMax">The z maximum.</param>
		/// <param name="phiMax">The phi maximum. Max Sweep value</param>
		Cylinder(const Transform *ObjectToWorld,
					const Transform *WorldToObject,
					bool reverseOrientation,
					Float radius,
					Float zMin,
					Float zMax,
					Float phiMax)
			: Shape(ObjectToWorld, WorldToObject, reverseOrientation),
			radius(radius),
			zMin(std::min(zMin, zMax)),
			zMax(std::max(zMin, zMax)),
			phiMax(Radians(Clamp(phiMax, 0, 360))) {}

		/// <summary>
		/// Returns an (axis-aligned) bounding box in the shape's object space
		/// </summary>
		/// <returns>
		/// Returns an (axis-aligned) bounding box in the shape's object space
		/// </returns>
		Bounds3f ObjectBound() const;

		/// <summary>
		/// Returns geometric information about a single ray-shape intersection corresponding to the first intersection
		/// if any in the (0, tMax) parametric range along the ray.
		/// </summary>
		/// <param name="ray">The ray.</param>
		/// <param name="tHit">Pointer to store the parametric distance along the ray if intersection is found. (for multiple intersection, the closest one is returned)</param>
		/// <param name="isect">Pointer to store information about the intersection. (Completly caputres the local geometric pproperties of a surface)</param>
		/// <param name="testAlphaTexture">if set to <c>true</c> the function shoudl perform an alpha cutting operation.</param>
		/// <returns></returns>
		bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
			bool testAlphaTexture) const;

		/// <summary>
		/// Predicate function that determines whether or not an intersection occurs, without returning any details about the intersection itself.
		/// </summary>
		/// <param name="ray">The ray.</param>
		/// <param name="testAlphaTexture">if set to <c>true</c> [test alpha texture].</param>
		/// <returns></returns>
		bool IntersectP(const Ray &ray, bool testAlphaTexture) const;

		/// <summary>
		/// Returns the area of the surface. (e.g. to use as area lights)
		/// A cylinder is just a rolled-up rectangle. Unrolled the height is zmax-zmin and the width is rphimax
		/// </summary>
		/// <returns></returns>
		Float Area() const;

		Interaction Sample(const Point2f &u, Float *pdf) const;

	protected:
		// Cylinder Private Data
		const Float radius, zMin, zMax, phiMax;
	};

	std::shared_ptr<Cylinder> CreateCylinderShape(const Transform *o2w,
		const Transform *w2o,
		bool reverseOrientation,
		const ParamSet &params);

}  // namespace pbrt

#endif // PBRT_SHAPES_CYLINDER_H