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

#ifndef PBRT_SHAPES_CURVE_H
#define PBRT_SHAPES_CURVE_H

 // shapes/curve.h*
#include "shape.h"

namespace pbrt {
	struct CurveCommon;
	// CurveType Declarations
	enum class CurveType {	Flat,		// always oriented to face the ray being intersected with them (for fine swept cylindrical shapes like hair or fur)
							Cylinder,	// span a few pixels on the screen (like spaghetti seen from not too far away), the Curve shape can compute a shading normal that makes the curve appear actually be a cylinder
							Ribbon		// shapes that don't actually have a cylindrical cross section (such as a blade of grass)
							};

	// CurveCommon Declarations
	struct CurveCommon {
		// initializes member variables with values passed into it for the control
		// points, curve width etc.
		// Curve points should be in the curve's object space.
		CurveCommon(const Point3f c[4],
					Float w0,
					Float w1,
					CurveType type,
					const Normal3f* norm);
		const CurveType type;
		Point3f cpObj[4];
		Float width[2];
		Normal3f n[2];	// for ribbon curves, stores a surface normal to orient the curve 
		// at each point
		Float normalAngle, invSinNormalAngle;
	};

	// Curve Declarations
	class Curve : public Shape {
	public:
		// Curve Public Methods
		Curve(	const Transform* ObjectToWorld,
				const Transform* WorldToObject,
				bool reverseOrientation,
				const std::shared_ptr<CurveCommon>& common,
				Float uMin,
				Float uMax)
			: Shape(ObjectToWorld, WorldToObject, reverseOrientation),
			common(common),
			uMin(uMin),
			uMax(uMax) {}
		Bounds3f ObjectBound() const;

		/// <summary>
		/// Returns geometric information about a single ray-shape intersection 
		/// corresponding to the first intersection if any in the (0, tMax) parametric 
		/// range along the ray.
		/// </summary>
		/// <param name="ray">The ray.</param>
		/// <param name="tHit">Pointer to store the parametric distance along the ray if intersection is found. (for multiple intersection, the closest one is returned)</param>
		/// <param name="isect">Pointer to store information about the intersection. (Completly caputres the local geometric pproperties of a surface)</param>
		/// <param name="testAlphaTexture">if set to <c>true</c> the function shoudl perform an alpha cutting operation.</param>
		/// <returns></returns>
		bool Intersect(	const Ray& ray,
						Float* tHit,
						SurfaceInteraction* isect,
						bool testAlphaTexture) const;

		/// <summary>
		/// Returns the area of the surface. (e.g. to use as area lights)
		/// </summary>
		/// <returns></returns>
		Float Area() const;

		/// <summary>
		/// Samples the specified u.
		/// </summary>
		/// <param name="u">The u.</param>
		/// <param name="pdf">The PDF.</param>
		/// <returns></returns>
		Interaction Sample(const Point2f& u,
							Float* pdf) const;

	private:
		// Curve Private Methods
		/// <summary>
		/// Tests whether the given ray intersects the given curve segment over the given parameter  [u0, u1]
		/// </summary>
		/// <param name="r">The r.</param>
		/// <param name="tHit">The t hit.</param>
		/// <param name="isect">The isect.</param>
		/// <param name="cp">The cp.</param>
		/// <param name="rayToObject">The ray to object.</param>
		/// <param name="u0">The u0.</param>
		/// <param name="u1">The u1.</param>
		/// <param name="depth">The depth.</param>
		/// <returns></returns>
		bool recursiveIntersect(const Ray& r,
								Float* tHit,
								SurfaceInteraction* isect,
								const Point3f cp[4],
								const Transform& rayToObject,
								Float u0,
								Float u1,
								int depth) const;

		// Curve Private Data
		const std::shared_ptr<CurveCommon> common;	// pointer to a CurveCommon structure,
		// stores the control points and other information about the curve that is shared
		// across curve segments
		const Float uMin, uMax;	// parametric range of u values 
	};

	std::vector<std::shared_ptr<Shape>> CreateCurveShape(	const Transform* o2w,
															const Transform* w2o,
															bool reverseOrientation,
															const ParamSet& params);

}  // namespace pbrt

#endif // PBRT_SHAPES_CURVE_H
