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

/**
 * Kd-Tree Accelerator
 * Binary space partioning tree (BSP).
 * BSPs can easily handle uneven distributions of geometry. 
 * (1) Splitting planes can be placed at arbitrary positions inside the overall bound 
 * 2) diff. parts of 3D space can be refined to different degrees)
 *  ==> diff Kd-tree / Octree
 *		Kd-tree:
 *			simply restrict the splitting planes to be perpendicular to one of the coordinate axes,
 *		Oc-tree:
 *			uses 3 axes-perpendicular lanes to simultaneously split the box into eight regions at each step
 *			typically by splitting down the center of the extent in each direction.
 */

#ifndef PBRT_ACCELERATORS_KDTREEACCEL_H
#define PBRT_ACCELERATORS_KDTREEACCEL_H

 // accelerators/kdtreeaccel.h*
#include "pbrt.h"
#include "primitive.h"

namespace pbrt {
	// KdTreeAccel Declarations

	struct KdAccelNode;
	struct BoundEdge;

	/// <summary>
	/// Kd-Tree for ray intersection acceleration.
	/// </summary>
	/// <seealso cref="Aggregate" />
	class KdTreeAccel : public Aggregate {
	public:
		// KdTreeAccel Public Methods
		/// <summary>
		/// Initializes a new instance of the <see cref="KdTreeAccel"/> class.
		/// </summary>
		/// <param name="p">The array of primitives.</param>
		/// <param name="isectCost">The isect cost.</param>
		/// <param name="traversalCost">The traversal cost.</param>
		/// <param name="emptyBonus">The empty bonus.</param>
		/// <param name="maxPrims">The maximum prims.</param>
		/// <param name="maxDepth">The maximum depth.</param>
		KdTreeAccel(std::vector<std::shared_ptr<Primitive>> p,
			int isectCost = 80, int traversalCost = 1,
			Float emptyBonus = 0.5, int maxPrims = 1, int maxDepth = -1);

		/// <summary>
		/// Worlds the bound.
		/// </summary>
		/// <returns></returns>
		Bounds3f WorldBound() const { return bounds; }

		/// <summary>
		/// Finalizes an instance of the <see cref="KdTreeAccel"/> class.
		/// </summary>
		~KdTreeAccel();
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
		// KdTreeAccel Private Methods
		/// <summary>
		/// Builds the tree with a recursive top-down algorithm.
		/// </summary>
		/// <param name="nodeNum">The offset into the array of KdAccelNodes to use for the node that it .</param>
		/// <param name="bounds">The bounding box that gives the region of space that the node covers.</param>
		/// <param name="primBounds">Array with all primitive bounds.</param>
		/// <param name="primNums">The indices of primitives that overlap the node..</param>
		/// <param name="nprims">The nprims.</param>
		/// <param name="depth">The depth.</param>
		/// <param name="edges">The edges.</param>
		/// <param name="prims0">The prims0.</param>
		/// <param name="prims1">The prims1.</param>
		/// <param name="badRefines">The bad refines.</param>
		void buildTree(int nodeNum, const Bounds3f& bounds,
			const std::vector<Bounds3f>& primBounds, int* primNums,
			int nprims, int depth,
			const std::unique_ptr<BoundEdge[]> edges[3], int* prims0,
			int* prims1, int badRefines = 0);

		// KdTreeAccel Private Data
		// Parameters to guide the decision that will be made as the tree is built
		const int isectCost,		// intersection cost for the SAH tree construction
					traversalCost,	// travel cost for the SAH tree construction 
					/* NOTE: in generall greater ratio --> visiting a kd-tree node is less expansive than visting a BVH node */
					maxPrims;

		/// <summary>
		/// A bonus in the SAH, for giving a slight preference to a a split where one of the children has no primitive overlapping it 
		/// Default: [0,1]
		/// </summary>
		const Float emptyBonus;
		std::vector<std::shared_ptr<Primitive>> primitives;

		/// <summary>
		/// Stores they the indices of more primitives overlapping in a segment.
		/// </summary>
		std::vector<int> primitiveIndices;
		/// <summary>
		/// The nodes
		/// </summary>
		KdAccelNode* nodes;
		
		/// <summary>
		/// The total number of nodes that have been allocated
		/// </summary>
		int nAllocedNodes;

		/// <summary>
		/// Records the next node that is available in the contiguous array of nodes
		/// </summary>
		int nextFreeNode;

		/// <summary>
		/// The bounding box of the Kd-Tree
		/// </summary>
		Bounds3f bounds;
	};

	struct KdToDo {
		const KdAccelNode* node;
		Float tMin, tMax;
	};

	std::shared_ptr<KdTreeAccel> CreateKdTreeAccelerator(
		std::vector<std::shared_ptr<Primitive>> prims, const ParamSet& ps);

}  // namespace pbrt

#endif // PBRT_ACCELERATORS_KDTREEACCEL_H