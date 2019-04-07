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

#ifndef PBRT_CORE_INTEGRATOR_H
#define PBRT_CORE_INTEGRATOR_H

 // core/integrator.h*
#include "pbrt.h"
// #include "primitive.h"
//#include "spectrum.h"
//#include "light.h"
//#include "reflection.h"
//#include "sampler.h"
//#include "material.h"


namespace pbrt {

	/// <summary>
	/// Integrator Interface
	/// </summary>
	class Integrator {
public:
	// Integrator Interface

	/// <summary>
	/// Finalizes an instance of the <see cref="Integrator"/> class.
	/// </summary>
	virtual ~Integrator();

	/**
	 * 
	 */
	 /// <summary>
	 /// Renders the specified scene it is passed a reference to the Scene
	 /// to use to compute an image of the scene or more generally,
	 /// a set of measurements of the scene lighting
	 /// </summary>
	 /// <param name="scene">The scene.</param>
	virtual void Render(const Scene &scene) = 0;
};

Spectrum UniformSampleAllLights(const Interaction &it, const Scene &scene,
	MemoryArena &arena, Sampler &sampler,
	const std::vector<int> &nLightSamples,
	bool handleMedia = false);

Spectrum UniformSampleOneLight(const Interaction &it, const Scene &scene,
	MemoryArena &arena, Sampler &sampler,
	bool handleMedia = false,
	const Distribution1D *lightDistrib = nullptr);

Spectrum EstimateDirect(const Interaction &it, const Point2f &uShading,
	const Light &light, const Point2f &uLight,
	const Scene &scene, Sampler &sampler,
	MemoryArena &arena, bool handleMedia = false,
	bool specular = false);

std::unique_ptr<Distribution1D> ComputeLightPowerDistribution(
	const Scene &scene);

// SamplerIntegrator Declarations

/// <summary>
/// 
/// </summary>
/// <seealso cref="Integrator" />
class SamplerIntegrator : public Integrator {
public:
	// SamplerIntegrator Public Methods

	/// <summary>
	/// Initializes a new instance of the <see cref="SamplerIntegrator"/> class.
	/// </summary>
	/// <param name="camera">The camera.</param>
	/// <param name="sampler">The sampler.</param>
	/// <param name="pixelBounds">The pixel bounds.</param>
	SamplerIntegrator(std::shared_ptr<const Camera> camera,
		std::shared_ptr<Sampler> sampler,
		const Bounds2i &pixelBounds)
		: camera(camera), sampler(sampler), pixelBounds(pixelBounds) {}

	/// <summary>
	/// Preprocesses the specified scene.
	/// </summary>
	/// <param name="scene">The scene.</param>
	/// <param name="sampler">The sampler.</param>
	virtual void Preprocess(const Scene & scene,
							Sampler &sampler) {}

	/// <summary>
	/// Renders the specified scene.
	/// </summary>
	/// <param name="scene">The scene.</param>
	void Render(const Scene & scene);


	/// <summary>
	/// Determine the amount of light arriving at the image plane along that ray.
	/// </summary>
	/// <param name="ray">The ray along which the incident radiance should be evaluated. </param>
	/// <param name="scene">The scene being rendered. The implementation will query the scene for information
	///						about the lights and geometry, and so on. </param>
	/// <param name="sampler">The sample geneartor used to solve the light transport equation via Monte Carlo
	///						  integration. </param>
	/// <param name="arena">The MemoryArena for efficient temporary memory allocation by the integrator.
	///						The integrator should assume that any memory it allocates with the arena will be 
	///						free shortly after the Li() method returns and thus should not use the arena to 
	///						allocate any memory that must persist for longer than is needed for the current ray. </param>
	/// <param name="depth">The number of ray bounces from the camera that have occurred up until the current call
	///						to Li(). </param>
	/// <returns></returns>
	virtual Spectrum Li(const RayDifferential & ray,
						const Scene &scene,
						Sampler &sampler,
						MemoryArena &arena,
						int depth = 0) const = 0;

	/// <summary>
	/// Speculars the refelect.
	/// </summary>
	/// <param name="ray">The ray.</param>
	/// <param name="isect">The isect.</param>
	/// <param name="scene">The scene.</param>
	/// <param name="sampler">The sampler.</param>
	/// <param name="arena">The arena.</param>
	/// <param name="depht">The depht.</param>
	/// <returns></returns>
	Spectrum SpecularReflect(const RayDifferential & ray,
								const SurfaceInteraction &isect,
								const Scene &scene,
								Sampler &sampler,
								Memory &arena,
								int depht) const;

	/// <summary>
	/// Speculars the transmit.
	/// </summary>
	/// <param name="ray">The ray.</param>
	/// <param name="isect">The isect.</param>
	/// <param name="scene">The scene.</param>
	/// <param name="sampler">The sampler.</param>
	/// <param name="arena">The arena.</param>
	/// <param name="depht">The depht.</param>
	/// <returns></returns>
	Spectrum SpecularTransmit(const RayDifferential & ray,
								const SurfaceInteraction &isect,
								const Scene &scene,
								Sampler &sampler,
								Memory &arena,
								int depht) const;

protected:
	// SamplerIntegrator Protected Data
	std::shared_ptr<const Camera> camera;

private:
	// SamplerIntegrator Private Data
	std::shared_ptr<Sampler> sampler;
	const Bounds2i pixelBounds;
};

}  // namespace pbrt

#endif // PBRT_CORE_INTEGRATOR_H