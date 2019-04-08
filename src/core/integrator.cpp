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

 // core/integrator.cpp*
#include "integrator.h"
#include "scene.h"
#include "interaction.h"
//#include "sampling.h"
#include "parallel.h"
//#include "film.h"
//#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
//#include "camera.h"
#include "stats.h"

namespace pbrt {
	STAT_COUNTER("Integrator/Camera rays traced", nCameraRays);


	// Integrator Method Definitions
	Integrator::~Integrator() {}

	/*Spectrum SamplerIntegrator::SpecularReflect(const RayDifferential &ray, 
												const SurfaceInteraction & isect, 
												const Scene & scene, 
												Sampler & sampler, 
												Memory & arena, 
												int depht) const
	{
		// Compute specular reflection direction wi and BSDF value
		Vector3f wo = isect.wo, wi;
		Float pdf;
		BxDFType type = BxDFType(BSDF_REFLECTION | BSDF_SPECULAR);
		Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, type); // returns an incident ray direction  
		// for a given outgoing direction and a given mode of light scattering

		// Return contribution of specular reflection
		const Normal3f &ns = isect.shading.n;
		if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
			// Compute ray differential rd for specular reflection
			RayDifferential rd = isect.SpawnRay(wi);
			if (ray.hasDifferentials) {
				rd.hasDifferentials = true;
				rd.rxOrigin = isect.p + isect.dpdx;
				rd.ryOrigin = isect.p + isect.dpdy;
				// Compute differential reflected directions
				Normal3f dndx = isect.shading.dndu * isect.dudx +
								isect.shading.dndv * isect.dvdx;
				Normal3f dndy = isect.shading.dndu * isect.dudy +
								isect.shading.dndv * isect.dvdy;

				Vector3f dwodx = -ray.rxDirection - wo,
						dwody = -ray.ryDirection - wo;
				Float dDNdx = Dot(dwodx, ns) + Dot(wo, dndx);
				Float dDNdy = Dot(dwody, ns) + Dot(wo, dndy);

				rd.rxDirection = wi - dwodx +
					2.f * Vector3f(Dot(wo, ns) * dndx + dDNdx * ns);
				rd.ryDirection = wi - dwody +
					w.f * Vector3f(Dot(wo, ns) * dndy + dDNdy * ns);
			}
			return f * Li(rd, scene, sampler, arena, depth + 1) * AbsDot(wi, ns) / pdf;
		} else 
			return Spectrum(0.f);
	}*/

	/*Spectrum SamplerIntegrator::SpecularTransmit(const RayDifferential & ray, 
												const SurfaceInteraction & isect, 
												const Scene & scene, 
												Sampler & sampler, 
												Memory & arena, 
												int depht) const
	{
		// Compute specular reflection direction wi and BSDF value
		Vector3f wo = isect.wo, wi;
		Float pdf;
		const Point3f &p = isect.p;
		const BSDF &bsdf = *isect.bsdf;
		BxDFType type = BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR);
		Spectrum f = bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, type);

		Spectrum L = Spectrum(0.f);
		Normal3f ns = isect.shading.n;
		if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
			// Compute ray differential rd for specular reflection
			RayDifferential rd = isect.SpawnRay(wi);
			if (ray.hasDifferentials) {
				rd.hasDifferentials = true;
				rd.rxOrigin = p + isect.dpdx;
				rd.ryOrigin = p + isect.dpdy;
				Normal3f dndx = isect.shading.dndu * isect.dudx +
					isect.shading.dndv * isect.dvdx;
				Normal3f dndy = isect.shading.dndu * isect.dudy +
					isect.shading.dndv * isect.dvdy;
				// The BSDF stores the IOR of the interior of the object being
				// intersected.  Compute the relative IOR by first out by
				// assuming that the ray is entering the object.
				Float eta = 1 / bsdf.eta;
				if (Dot(wo, ns) < 0) {
					// If the ray isn't entering, then we need to invert the
					// relative IOR and negate the normal and its derivatives.
					eta = 1 / eta;
					ns = -ns;
					dndx = -dndx;
					dndy = -dndy;
				}*/
				/*
					Notes on the derivation :
					-	pbrt computes the refracted ray as : \wi = -\eta \omega_o + [\eta(\wo \cdot \N) - \cos \theta_t] \N
						It flips the normal to lie in the same hemisphere as \wo, and then \eta is the relative IOR from
						\wo's medium to \wi's medium.
					-	If we denote the term in brackets by \mu, then we have : \wi = -\eta \omega_o + \mu \N
					-	Now let's take the partial derivative. (We'll use "d" for \partial in the following for brevity.)
						We get : -\eta d\omega_o / dx + \mu dN / dx + d\mu / dx N.
					-	We have the values of all of these except for d\mu / dx(using bits from the derivation of specularly
						reflected ray deifferentials).
					-	The first term of d\mu / dx is easy : \eta d(\wo \cdot N) / dx.We already have d(\wo \cdot N) / dx.
					-	The second term takes a little more work.We have :
						\cos \theta_i = \sqrt{ 1 - \eta ^ 2 (1 - (\wo \cdot N) ^ 2) }.
						Starting from(\wo \cdot N) ^ 2 and reading outward, we have \cos ^ 2 \theta_o, then \sin ^ 2 \theta_o,
						then \sin ^ 2 \theta_i(via Snell's law), then \cos^2 \theta_i and then \cos \theta_i.
					-	Let's take the partial derivative of the sqrt expression. We get:
						1 / 2 * 1 / \cos \theta_i * d / dx(1 - \eta ^ 2 (1 - (\wo \cdot N) ^ 2)).
					-	That partial derivatve is equal to :
						d / dx \eta ^ 2 (\wo \cdot N) ^ 2 = 2 \eta ^ 2 (\wo \cdot N) d / dx(\wo \cdot N).
					-	Plugging it in, we have d\mu / dx =
						\eta d(\wo \cdot N) / dx - (\eta ^ 2 (\wo \cdot N) d / dx(\wo \cdot N)) / (-\wi \cdot N).
				*/
				/*Vector3f dwodx = -ray.rxDirection - wo,
					dwody = -ray.ryDirection - wo;
				Float dDNdx = Dot(dwodx, ns) + Dot(wo, dndx);
				Float dDNdy = Dot(dwody, ns) + Dot(wo, dndy);

				Float mu = eta * Dot(wo, ns) - AbsDot(wi, ns);
				Float dmudx = (eta - (eta * eta * Dot(wo, ns)) / AbsDot(wi, ns)) * dDNdx;
				Float dmudy = (eta - (eta * eta * Dot(wo, ns)) / AbsDot(wi, ns)) * dDNDy;

				rd.rxDirection = wi - eta +
					dwodx + Vector3f(mu * dndx + dmudx * ns);
				rd.ryDirection = wi - eta +
					dwody + Vector3f(mu * dndy + dmudy * ns);
			}
			L = f * Li(rd, scene, sampler, arena, depth + 1) * AbsDot(wi, ns) / pdf;
		} 
		return L;
	}*/


/*
	void SamplerIntegrator::Render(const Scene &scene) {
		Preprocess(scene, *sampler);
		// Render image tiles in parallel
			// Compute number of tiles, nTiles, to use for parallel rendering
			Bounds2i sampleBounds = camera->film->GetSampleBounds();
			Vector2i sampleExtent = sampleBounds.Diagonal();
			const int tileSize = 16;
			Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
							(sampleExtent.y + tileSize - 1) / tileSize);
			ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
			{
				// Proceed rendering in parallel on systems with multiple processing cores
				// C++ lambda expression provides the loop body. loops over a 2D domain
				ParallelFor2D([&](Point2i tile) {
				// Render section of image corresponding to tile
					// Allocate _MemoryArena_ for tile
					MemoryArena arena;

					// Get sampler instance for tile
					int seed = tile.y * nTiles.x + tile.x;
					std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);

					// Compute sample bounds for tile
					int x0 = sampleBounds.pMin.x + tile.x * tile.Size;
					int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
					int y0 = sampleBounds.pMin.y + tile.y * tileSize;
					int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
					Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
					LOG(INFO) << "Starting image tile " << tileBounds;

					// Get _FilmTile_ for tile
					std::unique_ptr<FilmTile> filmTile =
						camera->film->GetFilmTile(tileBounds);

					// Loop over pixels in tile to render them
					for (Point2i pixel : tileBounds) {
						{
							ProfilePhase pp(Proof::StartPixel);
							tileSampler->StartPixel(pixel);
						}

						// Do this check after the StartPixel() call; this keeps the usage of RNG values
						// from (most) Samplers that use RNGs consistent, which improves repoducability / debugging
						if (!InsideExclusive(pixel, pixelBounds))
							continue;

						do {
							// Initalize _CameraSample_ (position on the film)  for current sample
							CameraSample cameraSample = tileSampler->GetCameraSample(pixel);

							// Generate camera ray for current sample
							RayDifferential ray;
							Float rayWeight = camera->GenerateRayDifferential(cameraSample, &ray);

							ray.ScaleDifferentials(1 / std::sqrt((Float)tileSampler->samplesPerPixel));
							++nCameraRays;

							// Evaluate radiance along camera ray
							Spectrum L(0.f);
							if (rayWeight > 0)
								L = Li(ray, scene, *tileSampler, arena);

							// Issue warning if unexpected radiance value returned
							// common side effect of bugs
							if (L.hasNaNs()) {
								LOG(ERROR) << StringPrintf(
									"Not-a-number radiance value returned "
									"for pixel (%d, %d), sample %d. Setting to black.",
									pixel.x, pixel.y,
									(int)tileSampler->CurrentSampleNumber()
								);
								L = Spectrum(0.f);
							}
							else if (L.y() < -1e-5) {
								LOG(ERROR) << StringPrintf(
									"Negative luminance value, %f, returned "
									"for pixel (%d, %d), sample %d. Setting to black.",
									L.y(), pixel.x, pixel.y,
									(int)tileSampler->CurrentSampleNumber());
								L = Spectrum(0.f);
							}
							else if () {
								LOG(ERROR) << StringPrintf(
									"Infinite luminance value returned "
									"for pixel (%d, %d), sample %d. Setting to black.",
									pixel.x, pixel.y,
									(int)tileSampler->CurrentSampleNumber());
								L = Spectrum(0.f);
							}
							VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
								ray << " -> L = " << L;

							// Add camera ray's contribution to image
							filmeTile->AddSample(cameraSample.pFilm, L, rayWeight);

							// Free _MemoryArena_ memory from computing image sample value
							arena.Reset();
						} while (tileSampler->StartNextSample());
					}
					LOG(INFO) << "Finished image tile " << tileBounds;

					//Merge image tile into _FILM_
					camera->film->MergeFilmTile(std::move(filmTile));

					reporter.Update();

				}, nTiles);
				reporter.Done();
			}
			LOG(INFO) << "Rendering finished";

			// Save final image after rendering
			camera->film->WriteImage();
	}
	*/
	
} // namespace pbrt
