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

#ifndef PBRT_CORE_PRIMITIVE_H
#define PBRT_CORE_PRIMITIVE_H

 // core/primitive.h*
#include "pbrt.h"
#include "shape.h"
#include "material.h"
#include "medium.h"
#include "transform.h"

namespace pbrt {

	// Primitive Declarations

	 /// <summary>
	 /// Bridge between the geometry processing and shading subsystems of pbrt.
	 /// </summary>
	class Primitive {
	public:
		// Primitive Interface
		virtual ~Primitive();
		
		/// <summary>
		/// Returns a box that encloses the primitive's geometry in world space.
		/// </summary>
		/// <returns>Box that encloses the primitive's geometry in world space.</returns>
		virtual Bounds3f WorldBound() const = 0;

		/// <summary>
		/// Intersection Test for a ray and a surface.
		/// Updates Ray::tMax with this value if an intersection is found.
		/// It is also responsible for initializing additional SurfaceInteraction member variables, including
		/// a pointer to the  Primtive the ray hit.
		/// </summary>
		/// <param name="r">Ray</param>
		/// <param name="">Surfaceinteraction</param>
		/// <returns>True if intersection is found.</returns>
		virtual bool Intersect(const Ray& r, SurfaceInteraction*) const = 0;

		/// <summary>
		/// Intersection Test for a ray.
		/// Updates Ray::tMax with this value if an intersection is found.
		/// </summary>
		/// <param name="r">Ray</param>
		/// <returns>True if intersection is found.</returns>
		virtual bool IntersectP(const Ray& r) const = 0;
		
		/// <summary>
		/// Gets the area light.
		/// If the primitive is not emissive, this method should return a nullptr.
		/// </summary>
		/// <returns> Pointer to the AreaLight that describes the primitive's emission destribution, if the primitive is itself a light source.
		/// </returns>
		virtual const AreaLight* GetAreaLight() const = 0;

		/// <summary>
		/// Gets the material.
		/// If nullptr is returned, ray intersections with the primtive should be ignored.
		/// The primitive only serves to delineate a volume of space for participating media.
		/// This method is also used to check if 2 rays have intersected the same object by comparing their Material pointers.
		/// </summary>
		/// <returns>Pointer to the material instance assigned to the primitive</returns>
		virtual const Material* GetMaterial() const = 0;

		/// <summary>
		/// Initializes representations of the light-scattering properties of the material at the intersection point on the surface.
		/// The BSDF object describes local light-scattering properties at the intersection point.
		/// If applicable, this method also initalizes a BSSRDF, which describes subsurface scattering inside the primitive - light 
		/// that enters the surface at points far from where it exits.
		/// </summary>
		/// <param name="isect">The isect.</param>
		/// <param name="arena">MemoryArena to allocate memory for the BSDF/BSSRDF.</param>
		/// <param name="mode">Transport enumerant that indicates whether the ray path that found this intersection point started from the camera or a light source.</param>
		/// <param name="allowMultipleLobes">if set to <c>true</c> allow multiple lobes and controls a detail of how some types of BRDF's are represented.</param>
		virtual void ComputeScatteringFunctions(SurfaceInteraction* isect,
			MemoryArena& arena,
			TransportMode mode,
			bool allowMultipleLobes) const = 0;
	};

	/// <summary>
	/// Represents a single shape (e.g. sphere) in the scene. One GeometricPrimitive is allocated for each shape
	/// in the scene description provided by the user.
	/// </summary>
	/// GeometricPrimitive Declarations
	/// <seealso cref="Primitive" />
	class GeometricPrimitive : public Primitive {
	public:
		// GeometricPrimitive Public Methods

		/// <summary>
		/// Returns a box that encloses the primitive's geometry in world space.
		/// Calls the Shape::WorldBound() method of its enclosed Shape to to the actual work
		/// </summary>
		/// <returns>
		/// Box that encloses the primitive's geometry in world space.
		/// </returns>
		virtual Bounds3f WorldBound() const;
		
		/// <summary>
		/// Intersection for the given ray with the given SurfaceInteraction.
		/// Calls the Shape::Intersect() method of its enclosed Shape to do the actual intersection test and initializes a
		/// SurfaceIntersaction to describe the intersection if any.
		/// </summary>
		/// <param name="r">Ray</param>
		/// <param name="isect">SurfaceInteraction.</param>
		/// <returns><c>true</c> if intersection is found.</returns>
		virtual bool Intersect(const Ray& r, SurfaceInteraction* isect) const;

		/// <summary>
		/// Intersection Test for a ray.
		/// Updates Ray::tMax with this value if an intersection is found.
		/// Calls the Shape::Intersect() method of its enclosed Shape to do the actual intersection testand initializes a
		/// SurfaceIntersaction to describe the intersection if any.
		/// </summary>
		/// <param name="r">Ray</param>
		/// <returns>
		/// True if intersection is found.
		/// </returns>
		virtual bool IntersectP(const Ray& r) const;

		/// <summary>
		/// Constructor
		/// Initializes a new instance of the <see cref="GeometricPrimitive"/> class.
		/// Initializes the member variables.
		/// </summary>
		/// <param name="shape">The shape.</param>
		/// <param name="material">The material.</param>
		/// <param name="areaLight">The area light.</param>
		/// <param name="mediumInterface">The medium interface.</param>
		GeometricPrimitive(const std::shared_ptr<Shape>& shape,
			const std::shared_ptr<Material>& material,
			const std::shared_ptr<AreaLight>& areaLight,
			const MediumInterface& mediumInterface);

		/// <summary>
		/// Gets the area light.
		/// If the primitive is not emissive, this method should return a nullptr.
		/// </summary>
		/// <returns>
		/// Returns the GeometricPrimitive::areaLight memeber
		/// </returns>
		const AreaLight* GetAreaLight() const;

		/// <summary>
		/// Gets the material.
		/// If nullptr is returned, ray intersections with the primtive should be ignored.
		/// The primitive only serves to delineate a volume of space for participating media.
		/// This method is also used to check if 2 rays have intersected the same object by comparing their Material pointers.
		/// </summary>
		/// <returns>
	/// Returns the GeometricPrimitive::material memeber
		/// </returns>
		const Material* GetMaterial() const;

		/// <summary>
		/// Initializes representations of the light-scattering properties of the material at the intersection point on the surface.
		/// The BSDF object describes local light-scattering properties at the intersection point.
		/// If applicable, this method also initalizes a BSSRDF, which describes subsurface scattering inside the primitive - light
		/// that enters the surface at points far from where it exits.
		/// It Forwars the requests on the Shape.
		/// </summary>
		/// <param name="isect">The isect.</param>
		/// <param name="arena">MemoryArena to allocate memory for the BSDF/BSSRDF.</param>
		/// <param name="mode">Transport enumerant that indicates whether the ray path that found this intersection point started from the camera or a light source.</param>
		/// <param name="allowMultipleLobes">if set to <c>true</c> allow multiple lobes and controls a detail of how some types of BRDF's are represented.</param>
		void ComputeScatteringFunctions(SurfaceInteraction* isect,
			MemoryArena& arena, TransportMode mode,
			bool allowMultipleLobes) const;

	private:
		// GeometricPrimitive Private Data
		std::shared_ptr<Shape> shape;			// Reference to a Shape
		std::shared_ptr<Material> material;		// Reference to its Shape's Material
		std::shared_ptr<AreaLight> areaLight;	// Pointer to an AreaLight object that describes its emission characteristics (nullptr if primitive does not emit light)
		MediumInterface mediumInterface;		// Encodes information about the participating media on the inside and outside of the primitive
	};

	// TransformedPrimitive Declarations
	/// <summary>
	/// Holds a single Primitivve and also includes an AnimatedTransform that is injected in between the underlying 
	/// primitive and its representation in the scene.
	/// THIS ENABLES: object instancing and primitives with animated transformations.
	/// </summary>
	/// <seealso cref="Primitive" />
	class TransformedPrimitive : public Primitive {
	public:
		// TransformedPrimitive Public Methods
		TransformedPrimitive(std::shared_ptr<Primitive>& primitive,
			const AnimatedTransform& PrimitiveToWorld);
		bool Intersect(const Ray& r, SurfaceInteraction* in) const;
		bool IntersectP(const Ray& r) const;
		const AreaLight* GetAreaLight() const { return nullptr; }
		const Material* GetMaterial() const { return nullptr; }
		void ComputeScatteringFunctions(SurfaceInteraction* isect,
			MemoryArena& arena, TransportMode mode,
			bool allowMultipleLobes) const {
			LOG(FATAL) <<
				"TransformedPrimitive::ComputeScatteringFunctions() shouldn't be "
				"called";
		}
		Bounds3f WorldBound() const {
			return PrimitiveToWorld.MotionBounds(primitive->WorldBound());
		}

	private:
		// TransformedPrimitive Private Data
		std::shared_ptr<Primitive> primitive;
		const AnimatedTransform PrimitiveToWorld;
	};

	// Aggregate Declarations
	class Aggregate : public Primitive {
	public:
		// Aggregate Public Methods
		const AreaLight* GetAreaLight() const;
		const Material* GetMaterial() const;
		void ComputeScatteringFunctions(SurfaceInteraction* isect,
			MemoryArena& arena, TransportMode mode,
			bool allowMultipleLobes) const;
	};

}  // namespace pbrt

#endif // PBRT_CORE_PRIMITIVE_H