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

#ifndef PBRT_CORE_SCENE_H
#define PBRT_CORE_SCENE_H

 // core/scene.h*
#include "pbrt.h"
#include "geometry.h"
//#include "primitive.h"
//#include "light.h"

namespace pbrt {

/// <summary>
/// Scene Declarations
/// </summary>
class Scene {
public:
	// Scene Public Methods
	Scene() {};	// TODO DELTE DEFAULT
	/// <summary>
	/// Initializes a new instance of the <see cref="Scene"/> class.
	/// </summary>
	/// <param name="aggregate">The aggregate.</param>
	/// <param name="lights">The lights.</param>
/*	Scene(std::shared_ptr<Primitive> aggregate,
		const std::vector<std::shared_ptr<Light>> &lights)
		: lights(lights), aggregate(aggregate) {
		// Scene Constructor Implementation
		worldBound = aggregate->WorldBound();
		for (const auto & light : lights) {
			light->Preprocess(*this);
			if (light->flags & (int)LightFlags::Infinite)
				infiniteLights.push_back(light);
		}
	}

	/// <summary>
	/// Worlds the bound.
	/// </summary>
	/// <returns></returns>
	const Bounds3f &WorldBound() const { return worldBound; }

	/// <summary>
	/// Traces the given ray into the scene.
	/// If Intersection it fills in the provided SurfaceInteraction structure with information
	/// about the closest intersection point along the ray.
	/// </summary>
	/// <param name="ray">The ray.</param>
	/// <param name="isect">The isect.</param>
	/// <returns>Indicating whether the ray intersected any of the primitives</returns>
	bool Intersect(const Ray &ray,
					SurfaceInteraction *isect) const;

	/// <summary>
	/// Checks for the existence of intersections along the ray.
	/// Does not return any information about those intersections.
	/// Because this routine doesn�t need to search for the closest intersection or compute
	/// any additional information about intersections, it is generally more efficient than Scene::Intersect(). 
	/// This routine is used for shadow rays. 
	/// </summary>
	/// <param name="ray">The ray.</param>
	/// <returns></returns>
	bool IntersectP(const Ray &ray) const;

	/// <summary>
	/// Intersects the tr.
	/// </summary>
	/// <param name="ray">The ray.</param>
	/// <param name="sampler">The sampler.</param>
	/// <param name="isect">The isect.</param>
	/// <param name="transmittance">The transmittance.</param>
	/// <returns></returns>
	bool IntersectTr(Ray ray,
					Sampler &sampler,
					SurfaceInteraction * isect,
					Spectrum *transmittance) const;
	

	// Scene Public Data
	std::vector<std::shared_ptr<Light>> lights;
	// Store infinite light sources separately for cases where we only want
	// to loop over them.
	std::vector<std::shared_ptr<Light>> infiniteLights;*/

/*private:
	// Scene Private Data
	std::shared_ptr<Primitive> aggregate;
	Bounds3f WorldBound;*/
};

}  // namespace pbrt

#endif // PBRT_CORE_SCENE_H