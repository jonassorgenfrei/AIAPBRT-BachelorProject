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

#ifndef PBRT_SHAPES_TRIANGLE_H
#define PBRT_SHAPES_TRIANGLE_H

 // shapes/triangle.h*
#include "shape.h"
#include "stats.h"
#include <map>

namespace pbrt {

	STAT_MEMORY_COUNTER("Memory/Triangle meshes", triMeshBytes);

	// Triangle Declarations
	struct TriangleMesh {
		// TriangleMesh Public Methods
		/// <summary>
		/// Initializes a new instance of the <see cref="TriangleMesh"/> struct.
		/// </summary>
		/// <param name="ObjectToWorld">The object to world transformation for the mesh.</param>
		/// <param name="nTriangles">The total number of triangles in the mesh.</param>
		/// <param name="vertexIndices">A pointer to an array of verex indices. For the ith triangle, its three vertex positions are 
		///								P[vertexIndices[3*i]], P[vertexIndices[3*i+1]], and P[vertexIndices[3*i+2]]</param>
		/// <param name="nVertices">The total number of vertices in the mesh.</param>
		/// <param name="P">An array of nVertices vertex positions.</param>
		/// <param name="S">An optional array of tangent vectors, one per vertex in the mesh. These are used to compute shading tangents.</param>
		/// <param name="N">An optional array of normal vectors, one per vertex in the mesh. If present, these are interpolated across triangle faces to compute shading normals.</param>
		/// <param name="uv">An optional array of parametric (u,v) values, one for each vertex.</param>
		/// <param name="alphaMask">An optional alpha mask texture, which can be used to cut away parts of triangle surfaces.</param>
		/// <param name="shadowAlphaMask">The shadow alpha mask.</param>
		/// <param name="faceIndices">The face indices.</param>
		TriangleMesh(const Transform& ObjectToWorld,
						int nTriangles,
						const int* vertexIndices,
						int nVertices,
						const Point3f* P,
						const Vector3f* S,
						const Normal3f* N,
						const Point2f* uv,
						const std::shared_ptr<Texture<Float>>& alphaMask,
						const std::shared_ptr<Texture<Float>>& shadowAlphaMask,
						const int* faceIndices);

		// TriangleMesh Data
		const int nTriangles, nVertices;
		std::vector<int> vertexIndices;
		std::unique_ptr<Point3f[]> p;	
		std::unique_ptr<Normal3f[]> n;
		std::unique_ptr<Vector3f[]> s;
		std::unique_ptr<Point2f[]> uv;
		std::shared_ptr<Texture<Float>> alphaMask, shadowAlphaMask;
		std::vector<int> faceIndices;
	};

	/// <summary>
	/// Implements actually the Shape interface. It represents a single triangle
	/// </summary>
	/// <seealso cref="Shape" />
	class Triangle : public Shape {
		
	public:
		// Triangle Public Methods
		/// <summary>
		/// Initializes a new instance of the <see cref="Triangle"/> class.
		/// </summary>
		/// <param name="ObjectToWorld">The object to world.</param>
		/// <param name="WorldToObject">The world to object.</param>
		/// <param name="reverseOrientation">if set to <c>true</c> [reverse orientation].</param>
		/// <param name="mesh">The mesh.</param>
		/// <param name="triNumber">The tri number.</param>
		Triangle(const Transform* ObjectToWorld,
					const Transform* WorldToObject,
					bool reverseOrientation,
					const std::shared_ptr<TriangleMesh>& mesh,
					int triNumber)
			: Shape(ObjectToWorld, WorldToObject, reverseOrientation), mesh(mesh) {
			v = &mesh->vertexIndices[3 * triNumber];
			triMeshBytes += sizeof(*this);
			faceIndex = mesh->faceIndices.size() ? mesh->faceIndices[triNumber] : 0;
		}

		Bounds3f ObjectBound() const;

		Bounds3f WorldBound() const;

		/// <summary>
		/// Returns geometric information about a single ray-shape intersection corresponding to the first intersection
		/// if any in the (0, tMax) parametric range along the ray.
		/// </summary>
		/// <param name="ray">The ray.</param>
		/// <param name="tHit">Pointer to store the parametric distance along the ray if intersection is found. (for multiple intersection, the closest one is returned)</param>
		/// <param name="isect">Pointer to store information about the intersection. (Completly caputres the local geometric pproperties of a surface)</param>
		/// <param name="testAlphaTexture">if set to <c>true</c> the function shoudl perform an alpha cutting operation.</param>
		/// <returns></returns>
		bool Intersect(const Ray& ray,
						Float* tHit,
						SurfaceInteraction* isect,
						bool testAlphaTexture = true) const;

		/// <summary>
		/// Predicate function that determines whether or not an intersection occurs, without returning any details about the intersection itself.
		/// </summary>
		/// <param name="ray">The ray.</param>
		/// <param name="testAlphaTexture">if set to <c>true</c> [test alpha texture].</param>
		/// <returns></returns>
		bool IntersectP(const Ray& ray,
						bool testAlphaTexture = true) const;

		/// <summary>
		/// Returns the area of the surface. (e.g. to use as area lights)
		/// </summary>
		/// <returns></returns>
		Float Area() const;

		using Shape::Sample;  // Bring in the other Sample() overload.
		Interaction Sample(const Point2f& u, Float* pdf) const;

		// Returns the solid angle subtended by the triangle w.r.t. the given
		// reference point p.
		Float SolidAngle(const Point3f& p, int nSamples = 0) const;

	private:
		// Triangle Private Methods
		/// <summary>
		/// Utility Function, that returns the (u,v) coordinates for the three vertices of the triangle
		/// </summary>
		/// <param name="uv">The uv.</param>
		void GetUVs(Point2f uv[3]) const {
			if (mesh->uv) {
				uv[0] = mesh->uv[v[0]];
				uv[1] = mesh->uv[v[1]];
				uv[2] = mesh->uv[v[2]];
			}
			else {	// returning default values if explicit (u,v) coordinates were not specified with the mesh
				uv[0] = Point2f(0, 0);
				uv[1] = Point2f(1, 0);
				uv[2] = Point2f(1, 1);
			}
		}

		// Triangle Private Data
		std::shared_ptr<TriangleMesh> mesh;	// pointer to the parent TriangleMesh
		const int* v;	// Stores a pointer to the first vertex index.
		int faceIndex;
	};
	/// <summary>
	/// Creates the triangle mesh as well as a Triangle for each triangle in the mesh.
	/// </summary>
	/// <param name="ObjectToWorld">The object to world.</param>
	/// <param name="WorldToObject">The world to object.</param>
	/// <param name="reverseOrientation">if set to <c>true</c> [reverse orientation].</param>
	/// <param name="nTriangles">The n triangles.</param>
	/// <param name="vertexIndices">The vertex indices.</param>
	/// <param name="nVertices">The n vertices.</param>
	/// <param name="p">The p.</param>
	/// <param name="s">The s.</param>
	/// <param name="n">The n.</param>
	/// <param name="uv">The uv.</param>
	/// <param name="alphaMask">The alpha mask.</param>
	/// <param name="shadowAlphaMask">The shadow alpha mask.</param>
	/// <param name="faceIndices">The face indices.</param>
	/// <returns>Vector of triangle shapes</returns>
	std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(	const Transform* o2w, 
															const Transform* w2o,
															bool reverseOrientation,
															int nTriangles,
															const int* vertexIndices,
															int nVertices,
															const Point3f* p,
															const Vector3f* s, 
															const Normal3f* n, 
															const Point2f* uv,
															const std::shared_ptr<Texture<Float>>& alphaTexture,
															const std::shared_ptr<Texture<Float>>& shadowAlphaTexture,
															const int* faceIndices = nullptr);

	std::vector<std::shared_ptr<Shape>> CreateTriangleMeshShape(const Transform* o2w, 
																const Transform* w2o,
																bool reverseOrientation,
																const ParamSet& params,
																std::map<std::string, std::shared_ptr<Texture<Float>>>* floatTextures =	nullptr);

	bool WritePlyFile(	const std::string& filename,
						int nTriangles,
						const int* vertexIndices,
						int nVertices,
						const Point3f* P,
						const Vector3f* S,
						const Normal3f* N,
						const Point2f* UV,
						const int* faceIndices);

}  // namespace pbrt

#endif // PBRT_SHAPES_TRIANGLE_H