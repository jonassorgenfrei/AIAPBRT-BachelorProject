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

#ifndef PBRT_ACCELERATORS_BVH_H
#define PBRT_ACCELERATORS_BVH_H

 // accelerators/bvh.h*
#include "pbrt.h"
#include "primitive.h"
#include <atomic>

namespace pbrt {
	/// <summary>
	/// Structure that represents a node of the BVH
	/// </summary>
	struct BVHBuildNode;

	// BVHAccel Forward Declarations
	/// <summary>
	/// Structure to store information about Primitive in the BVH
	/// </summary>
	struct BVHPrimitiveInfo;
	struct MortonPrimitive;
	struct LinearBVHNode;

	// BVHAccel Declarations
	/// <summary>
	/// Bounding Volume Hierachies Acceleration.
	/// Total number of nodes => 2n-1, where n is the number of primitives
	/// N leaf nodes and n-1 interior nodes 
	/// If leaves store multiple primitives, fewer nodes are needed
	/// </summary>
	/// <seealso cref="Aggregate" />
	class BVHAccel : public Aggregate {
	public:
		// BVHAccel Public Types
		enum class SplitMethod {	SAH,		// Surface Area heuristic 
									HLBVH,		// Hierarchical Linear Bounding Volume Hierarchy (can be constructed more efficiently and easily parallized)
									Middle,		// (less computation to build the tree but create fairly low-quality trees)
									EqualCounts // (less computation to build the tree but create fairly low-quality trees)
								};

		// BVHAccel Public Methods
		/// <summary>
		/// Constructor
		/// Initializes a new instance of the <see cref="BVHAccel"/> class.
		/// </summary>
		/// <param name="p">The primitives to be stored.</param>
		/// <param name="maxPrimsInNode">The maximum number of primitives in any leaf node.</param>
		/// <param name="splitMethod">The split method enumerator, to determine the partitioning algorithm. Default: SAH.</param>
		BVHAccel(std::vector<std::shared_ptr<Primitive>> p,
				int maxPrimsInNode = 1,
				SplitMethod splitMethod = SplitMethod::SAH);

		/// <summary>
		/// Worlds the bound.
		/// </summary>
		/// <returns></returns>
		Bounds3f WorldBound() const;
		
		/// <summary>
		/// Finalizes an instance of the <see cref="BVHAccel"/> class.
		/// </summary>
		~BVHAccel();
		
		/// <summary>
		/// Intersects the specified ray.
		/// </summary>
		/// <param name="ray">The ray.</param>
		/// <param name="isect">The isect.</param>
		/// <returns></returns>
		bool Intersect(const Ray& ray, SurfaceInteraction* isect) const;

		/// <summary>
		/// Intersects the p.
		/// </summary>
		/// <param name="ray">The ray.</param>
		/// <returns></returns>
		bool IntersectP(const Ray& ray) const;

	private:
		// BVHAccel Private Methods
		
		/// <summary>
		/// Recursives the build.
		/// Builds the tree structure for the given range [start, end).
		/// If the range covers only a single primitive, than the recursion has bottomed out 
		/// and a leaf node is created.
		/// Otherwise the algorithm partitions the elements of the array in that range using
		/// one of the partitioning algorithms and reorders the array elements in the range accordingly, 
		/// so that the ranges from [start, mid) and [mid, end) represent the partitioned subets. 
		/// And create a new node for those children.
		/// </summary>
		/// <param name="arena">The memory area for memory management.</param>
		/// <param name="primitiveInfo">The primitive information.</param>
		/// <param name="start">The start of the subset Range.</param>
		/// <param name="end">The end of the subset Range.</param>
		/// <param name="totalNodes">Pointer to the number of total nodes of the BVH, that have been created.</param>
		/// <param name="orderedPrims">The ordered primitive array.</param>
		/// <returns>Pointer to the root of the tree</returns>
		BVHBuildNode* recursiveBuild(
			MemoryArena& arena, std::vector<BVHPrimitiveInfo>& primitiveInfo,
			int start, int end, int* totalNodes,
			std::vector<std::shared_ptr<Primitive>>& orderedPrims);
		
		/// <summary>
		/// Builds the tree using HLBVH Algorithm.
		/// </summary>
		/// <param name="arena">The memory area for memory management.</param>
		/// <param name="primitiveInfo">The primitive information.</param>
		/// <param name="totalNodes">Pointer to the number of total nodes of the BVH, that have been created.</param>
		/// <param name="orderedPrims">The ordered primitive array.</param>
		/// <returns>Pointer to the root of the tree</returns>
		BVHBuildNode* HLBVHBuild(
			MemoryArena& arena, const std::vector<BVHPrimitiveInfo>& primitiveInfo,
			int* totalNodes,
			std::vector<std::shared_ptr<Primitive>>& orderedPrims) const;


		BVHBuildNode* emitLBVH(
			BVHBuildNode*& buildNodes,
			const std::vector<BVHPrimitiveInfo>& primitiveInfo,
			MortonPrimitive* mortonPrims, int nPrimitives, int* totalNodes,
			std::vector<std::shared_ptr<Primitive>>& orderedPrims,
			std::atomic<int>* orderedPrimsOffset, int bitIndex) const;
		BVHBuildNode* buildUpperSAH(MemoryArena& arena,
			std::vector<BVHBuildNode*>& treeletRoots,
			int start, int end, int* totalNodes) const;
		
		int flattenBVHTree(BVHBuildNode* node, int* offset);

		// BVHAccel Private Data
		const int maxPrimsInNode;		// The maximum number of primitives in any leaf node
		const SplitMethod splitMethod;	//
		std::vector<std::shared_ptr<Primitive>> primitives;
		LinearBVHNode* nodes = nullptr;
	};

	std::shared_ptr<BVHAccel> CreateBVHAccelerator(
		std::vector<std::shared_ptr<Primitive>> prims, const ParamSet& ps);

}  // namespace pbrt

#endif // PBRT_ACCELERATORS_BVH_H