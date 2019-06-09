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


 // accelerators/kdtreeaccel.cpp*
#include "accelerators/kdtreeaccel.h"
#include "paramset.h"
#include "interaction.h"
#include "stats.h"
#include <algorithm>

namespace pbrt {

	// KdTreeAccel Local Declarations
	/// <summary>
	/// Representation for leafe and interior nodes of the Kd-Tree Structure
	/// leaf notes only store 8 Bytes of memory --> 8 nodes will fit into a 64-byte cache line (Speed Increase)
	/// </summary>
	struct KdAccelNode {
		// KdAccelNode Methods
		void InitLeaf(int* primNums, int np, std::vector<int>* primitiveIndices);

		/// <summary>
		/// Initializes the <see ref="KdAccelNode">  as an interior node
		/// </summary>
		/// <param name="axis">The split axis x=0, y=1, z=2.</param>
		/// <param name="ac">The amount of above children.</param>
		/// <param name="s">The position along the chosen split axis, where the node splits space.</param>
		void InitInterior(int axis, int ac, Float s) {
			split = s;
			flags = axis;	// store the split axis
			aboveChild |= (ac << 2);
		}

		/// <summary>
		/// Getter Function for the Split Pos
		/// </summary>
		/// <returns> Ths Split position.</returns>
		Float SplitPos() const { return split; }
	
		/// <summary>
		/// Getter Function for the amount of stored primitives.
		/// </summary>
		/// <returns>Number of primitives stored in that node</returns>
		int nPrimitives() const { return nPrims >> 2; }

		/// <summary>
		/// Getter Function for the Split axis.
		/// </summary>
		/// <returns>0 if split axis is x, 1 if split axis is y, 2 if split axis is z, if interior node it returns 3</returns>
		int SplitAxis() const { return flags & 3; }

		/// <summary>
		/// Determines whether this instance is leaf.
		/// </summary>
		/// <returns>
		///   <c>true</c> if this instance is leaf; otherwise, <c>false</c>.
		/// </returns>
		bool IsLeaf() const { return (flags & 3) == 3; }

		/// <summary>
		/// Getter Methods for the offset to the other above child.
		/// </summary>
		/// <returns>The offset to the other abve child.</returns>
		int AboveChild() const { return aboveChild >> 2; }

		union {	
			Float split;                 // Interior
			int onePrimitive;            // Leaf
			/// <summary>
			/// Stores the offset to the first index for the leaf; the indicies for the rest, directly follow
			/// </summary>
			int primitiveIndicesOffset;  // Leaf
		};

	private:
		union {	// leaf nodes and nPrimes share the same storage
			/// <summary>
			///  Two low-order bits ==> {0 - xsplit (interior), 1 - ysplit (interior), 2 - zsplit (interior); 3 - leaf node}
			/// </summary>
			int flags;       // Both 
			/// <summary>
			/// Upper 30 bits available to record the number of primitives overlapping it 
			/// </summary>
			int nPrims;      // Leaf
			/// <summary>
			/// Pointer to the other child of the leave. (the first one is stored immediately after the current node)
			/// </summary>
			int aboveChild;  // Interior
		};
	};

	/// <summary>
	/// Enum define if the projected Edge is Start or Endpoint on the particular axis
	/// </summary>
	enum class EdgeType { Start, End };

	/// <summary>
	/// Structure to represent points of edges of the bounding boxes projected onto one of the 3 axes
	/// </summary>
	struct BoundEdge {
		// BoundEdge Public Methods
		BoundEdge() {}
		BoundEdge(Float t, int primNum, bool starting) : t(t), primNum(primNum) {
			type = starting ? EdgeType::Start : EdgeType::End;
		}
		Float t;	// Value to define an ordering
		int primNum;
		EdgeType type;
	};

	// KdTreeAccel Method Definitions
	KdTreeAccel::KdTreeAccel(std::vector<std::shared_ptr<Primitive>> p,
		int isectCost, int traversalCost, Float emptyBonus,
		int maxPrims, int maxDepth)
		: isectCost(isectCost),
		traversalCost(traversalCost),
		maxPrims(maxPrims),
		emptyBonus(emptyBonus),
		primitives(std::move(p)) {

		// Build kd-tree for accelerator
		ProfilePhase _(Prof::AccelConstruction);
		nextFreeNode = nAllocedNodes = 0;	// ensure that allocation will be done immediately when the first node of the tree is initalized
		if (maxDepth <= 0)
			maxDepth = std::round(8 + 1.3f * Log2Int(int64_t(primitives.size())));	// determine the max tree depth if not manually set

		// Compute bounds for kd-tree construction
		std::vector<Bounds3f> primBounds;	// stores the bounding boxes of the primitives along the wy
		primBounds.reserve(primitives.size());
		for (const std::shared_ptr<Primitive>& prim : primitives) {
			Bounds3f b = prim->WorldBound();
			bounds = Union(bounds, b);	// copmutes the bounding box of the whole structure
			primBounds.push_back(b);
		}

		// Allocate working memory for kd-tree construction 
		std::unique_ptr<BoundEdge[]> edges[3];
		for (int i = 0; i < 3; ++i)
			edges[i].reset(new BoundEdge[2 * primitives.size()]);	// for the Bounding box projection on the axes, reuse this on each iteration
		// reserve memory of integers for the prims0 & prims1 for storing worst-case possible number of overlapping primitive numbers
		std::unique_ptr<int[]> prims0(new int[primitives.size()]);
		std::unique_ptr<int[]> prims1(new int[(maxDepth + 1) * primitives.size()]);	// more memory for prims1, cause prims 0 only stores prims for a single level at a time

		// Initialize _primNums_ for kd-tree construction
		std::unique_ptr<int[]> primNums(new int[primitives.size()]);
		for (size_t i = 0; i < primitives.size(); ++i) primNums[i] = i;	// because all nodes overlap the root node; --> inialize an array with all all primitive indices 

		// Start recursive construction of kd-tree
		buildTree(0, bounds, primBounds, primNums.get(), primitives.size(),
			maxDepth, edges, prims0.get(), prims1.get());	
		// calls the the build tree submethode; last three parameters: points to the data allocated in the allocate working memory for kd-tree construction fragment
	}

	/// <summary>
	/// Initializes the <see ref="KdAccelNode"> to a leaf node. 
	/// Setting the flag to 3.
	/// </summary>
	/// <param name="primNums">The prim nums.</param>
	/// <param name="np">The number of primitives in the leaf.</param>
	/// <param name="primitiveIndices">The primitive indices.</param>
	void KdAccelNode::InitLeaf(int* primNums, int np,
		std::vector<int> * primitiveIndices) {
		flags = 3;
		nPrims |= (np << 2);	// shift 2 bits to the, because of shared memory
		// Store primitive ids for leaf node
		if (np == 0)
			onePrimitive = 0;		// leaf node without overlapping primitives
		else if (np == 1)
			onePrimitive = primNums[0];	// leaf node with one overlapping primitives
		else {	// leaf node with more than one overlapping primitive
			primitiveIndicesOffset = primitiveIndices->size();	
			for (int i = 0; i < np; ++i) primitiveIndices->push_back(primNums[i]);	// allocate storage in the primitive indice array
		}
	}

	KdTreeAccel::~KdTreeAccel() { FreeAligned(nodes); }

	void KdTreeAccel::buildTree(int nodeNum, const Bounds3f& nodeBounds,
		const std::vector<Bounds3f> & allPrimBounds,
		int* primNums, int nPrimitives, int depth,
		const std::unique_ptr<BoundEdge[]> edges[3],
		int* prims0, int* prims1, int badRefines) {

		CHECK_EQ(nodeNum, nextFreeNode);
		// Get next free node from _nodes_ array
		if (nextFreeNode == nAllocedNodes) {	// 
			int nNewAllocNodes = std::max(2 * nAllocedNodes, 512);	// allocation size => twice as big as previous
			
			KdAccelNode* n = AllocAligned<KdAccelNode>(nNewAllocNodes);	// reallocate node memory with new size
			if (nAllocedNodes > 0) {		// only if initail block was previously allocated
				memcpy(n, nodes, nAllocedNodes * sizeof(KdAccelNode));	// copy old values
				FreeAligned(nodes);	// free old memory
			}
			nodes = n;
			nAllocedNodes = nNewAllocNodes;
		}
		++nextFreeNode;

		// Initialize leaf node if termination criteria met
		if (nPrimitives <= maxPrims || depth == 0) {	// check if the number of primitives in the region is sufficiently small or max depth has been reached
			nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);	// create leaf node
			return;
		}

		// Initialize interior node and continue recursion

		/* NOTE: USING Surface Area Heuristic */

		// Choose split axis position for interior node
		int bestAxis = -1, bestOffset = -1;		// record the axis & boundingbox edge index that have given the lowest cost so far
		Float bestCost = Infinity;				// the current best (lowest) cost
		Float oldCost = isectCost * Float(nPrimitives);
		Float totalSA = nodeBounds.SurfaceArea();	// nodes surface area
		Float invTotalSA = 1 / totalSA;			// reciprocal of the node's surface area (will be used when computing the probabilities of rays passing through the candidate children nodes)
		Vector3f d = nodeBounds.pMax - nodeBounds.pMin;

		// Choose which axis to split along
		int axis = nodeBounds.MaximumExtent();	// first tries to find a split along the largest spatial extent
		int retries = 0;
	retrySplit:

		// Initialize edges for _axis_
		// using the bounding boxes of the overlapping primitives
		for (int i = 0; i < nPrimitives; ++i) {
			int pn = primNums[i];
			const Bounds3f& bounds = allPrimBounds[pn];
			edges[axis][2 * i] = BoundEdge(bounds.pMin[axis], pn, true);
			edges[axis][2 * i + 1] = BoundEdge(bounds.pMax[axis], pn, false);
		}

		// Sort _edges_ for _axis_
		// from low to high along the axis so that it can sweep over the box edges from first to last
		std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
			[](const BoundEdge & e0, const BoundEdge & e1) -> bool {
				if (e0.t == e1.t)
					return (int)e0.type < (int)e1.type;		// try to break the tie by comparing the node's types
				else
					return e0.t < e1.t;
			});

		// Compute cost of all splits for _axis_ to find best
		int nBelow = 0, nAbove = nPrimitives;	// nBelow --> primitives that end up below the splitting plane; nAbove -->primitives that end up above the splitting plane
		for (int i = 0; i < 2 * nPrimitives; ++i) {
			if (edges[axis][i].type == EdgeType::End) --nAbove;		// update the primtive counts
			Float edgeT = edges[axis][i].t;
			if (edgeT > nodeBounds.pMin[axis] && edgeT < nodeBounds.pMax[axis]) {
				// Compute cost for split at _i_th edge

				// Compute child surface areas for split at _edgeT_
				int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
				Float belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
					(edgeT - nodeBounds.pMin[axis]) *
					(d[otherAxis0] + d[otherAxis1]));					// surface area of one child candidate bounds (computed by adding up the areas of the 6 faces)
				Float aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
					(nodeBounds.pMax[axis] - edgeT) *
					(d[otherAxis0] + d[otherAxis1]));					// surface area of other child candidate bounds
				// compute the the cost for this particular split 
				Float pBelow = belowSA * invTotalSA;
				Float pAbove = aboveSA * invTotalSA;
				Float eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0;		// bonus if one of the regions is completly empty
				Float cost =
					traversalCost +
					isectCost * (1 - eb) * (pBelow * nBelow + pAbove * nAbove);

				// Update best split if this is lowest cost so far
				if (cost < bestCost) {
					bestCost = cost;
					bestAxis = axis;
					bestOffset = i;
				}
			}
			if (edges[axis][i].type == EdgeType::Start) ++nBelow;	// update the primtive counts
		}
		CHECK(nBelow == nPrimitives && nAbove == 0);

		// Create leaf if no good splits were found
		if (bestAxis == -1 && retries < 2) {
			++retries;
			axis = (axis + 1) % 3;
			goto retrySplit;
		}
		if (bestCost > oldCost) ++badRefines;	// allowing a few slightly poor refinements
		if ((bestCost > 4 * oldCost && nPrimitives < 16)	|| // check if split cost is higher than not splitting the node && number of primitives is not too high
										bestAxis == -1		||
										badRefines == 3) {
			nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);	// make leaf node
			return;
		}

		// Classify primitives with respect to split being above or below or on both sides
		int n0 = 0, n1 = 0;
		for (int i = 0; i < bestOffset; ++i)		
			if (edges[bestAxis][i].type == EdgeType::Start)
				prims0[n0++] = edges[bestAxis][i].primNum;
		for (int i = bestOffset + 1; i < 2 * nPrimitives; ++i)		// skip bestOffset entry
			if (edges[bestAxis][i].type == EdgeType::End)
				prims1[n1++] = edges[bestAxis][i].primNum;

		// Recursively initialize children nodes
		Float tSplit = edges[bestAxis][bestOffset].t;
		Bounds3f bounds0 = nodeBounds, bounds1 = nodeBounds;
		bounds0.pMax[bestAxis] = bounds1.pMin[bestAxis] = tSplit;
		buildTree(nodeNum + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, edges,
			prims0, prims1 + nPrimitives, badRefines);
		int aboveChild = nextFreeNode;	// used for the above child
		nodes[nodeNum].InitInterior(bestAxis, aboveChild, tSplit);
		buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1, depth - 1, edges,
			prims0, prims1 + nPrimitives, badRefines);
	}

	bool KdTreeAccel::Intersect(const Ray & ray, SurfaceInteraction * isect) const {
		ProfilePhase p(Prof::AccelIntersect);
		// Compute initial parametric range of ray inside kd-tree extent
		Float tMin, tMax;	// overall parametric range  of the ray's overlap with the tree
		if (!bounds.IntersectP(ray, &tMin, &tMax)) {
			return false;		// return false if the ray misses the overall bounding box
		}

		// Prepare to traverse kd-tree for ray
		Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);	// precompute to avoid an division each iteration step
		PBRT_CONSTEXPR int maxTodo = 64;	// max number = max depth in the kd-tree
		KdToDo todo[maxTodo];	// ordered so that the last active entry in the array is the next node that should be considered
		int todoPos = 0;

		// Traverse kd-tree nodes in order (depth-first front-to-back traversal) for ray
		bool hit = false;
		const KdAccelNode * node = &nodes[0];
		while (node != nullptr) {
			// Bail out if we found a hit closer than the current node
			if (ray.tMax < tMin) break;	// only break where tMin is beyond the intersection, than it's certain that there is no closer intersection with some other primitive
			if (!node->IsLeaf()) {	
				// Process kd-tree interior node

				// Compute parametric distance along ray to split plane
				int axis = node->SplitAxis();
				Float tPlane = (node->SplitPos() - ray.o[axis]) * invDir[axis]; // intersect ray with node's splitting plane
				// this helps to determine the order of nodes to check

				// Get node children pointers for ray by determining the order in which the ray encounters the children nodes
				const KdAccelNode * firstChild,		// near node
								*secondChild;		// far node
				int belowFirst =
					(ray.o[axis] < node->SplitPos()) ||
					(ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0);	// handling the case that the origin lies on the splitting plane; then the direction is used to discriminate
				if (belowFirst) {
					firstChild = node + 1;
					secondChild = &nodes[node->AboveChild()];
				}
				else {
					firstChild = &nodes[node->AboveChild()];
					secondChild = node + 1;
				}

				// Advance to next child node, possibly enqueue other child
				if (tPlane > tMax || tPlane <= 0)	// only near node needs to be processed
					node = firstChild;				// becaus the ray doesnt overlap the far node due to its direction facing away from it or tsplit > tmax
				else if (tPlane < tMin)				// only far node needs to be processed 
					node = secondChild;				// since ray doesnt overlapp the near node
				else {
					// Enqueue _secondChild_ in todo list
					todo[todoPos].node = secondChild;	// add far node to the todo stack
					todo[todoPos].tMin = tPlane;
					todo[todoPos].tMax = tMax;
					++todoPos;
					node = firstChild;					// process near node next
					tMax = tPlane;
				}
			}
			else {
				// Check for intersections inside leaf node
				int nPrimitives = node->nPrimitives();
				if (nPrimitives == 1) {
					const std::shared_ptr<Primitive>& p =
						primitives[node->onePrimitive];
					// Check one primitive inside leaf node
					if (p->Intersect(ray, isect)) hit = true;	// passing intersection request on the primitive
				}
				else {
					for (int i = 0; i < nPrimitives; ++i) {
						int index =
							primitiveIndices[node->primitiveIndicesOffset + i];
						const std::shared_ptr<Primitive>& p = primitives[index];
						// Check one primitive inside leaf node
						if (p->Intersect(ray, isect)) hit = true;
					}
				}

				// Grab next node to process from todo list
				if (todoPos > 0) {
					--todoPos;
					node = todo[todoPos].node;
					// always hold the parametric range for the ray's overlap wiith the current node
					tMin = todo[todoPos].tMin;	
					tMax = todo[todoPos].tMax;
				}
				else
					break;	// end if todo array is empty
			}
		}
		return hit;
	}

	bool KdTreeAccel::IntersectP(const Ray & ray) const {
		ProfilePhase p(Prof::AccelIntersectP);
		// Compute initial parametric range of ray inside kd-tree extent
		Float tMin, tMax;
		if (!bounds.IntersectP(ray, &tMin, &tMax)) {
			return false;
		}

		// Prepare to traverse kd-tree for ray
		Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
		PBRT_CONSTEXPR int maxTodo = 64;
		KdToDo todo[maxTodo];
		int todoPos = 0;
		const KdAccelNode * node = &nodes[0];
		while (node != nullptr) {
			if (node->IsLeaf()) {
				// Check for shadow ray intersections inside leaf node
				int nPrimitives = node->nPrimitives();
				if (nPrimitives == 1) {
					const std::shared_ptr<Primitive>& p =
						primitives[node->onePrimitive];
					if (p->IntersectP(ray)) {
						return true;
					}
				}
				else {
					for (int i = 0; i < nPrimitives; ++i) {
						int primitiveIndex =
							primitiveIndices[node->primitiveIndicesOffset + i];
						const std::shared_ptr<Primitive>& prim =
							primitives[primitiveIndex];
						if (prim->IntersectP(ray)) {
							return true;
						}
					}
				}

				// Grab next node to process from todo list
				if (todoPos > 0) {
					--todoPos;
					node = todo[todoPos].node;
					tMin = todo[todoPos].tMin;
					tMax = todo[todoPos].tMax;
				}
				else
					break;
			}
			else {
				// Process kd-tree interior node

				// Compute parametric distance along ray to split plane
				int axis = node->SplitAxis();
				Float tPlane = (node->SplitPos() - ray.o[axis]) * invDir[axis];

				// Get node children pointers for ray
				const KdAccelNode * firstChild, *secondChild;
				int belowFirst =
					(ray.o[axis] < node->SplitPos()) ||
					(ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0);
				if (belowFirst) {
					firstChild = node + 1;
					secondChild = &nodes[node->AboveChild()];
				}
				else {
					firstChild = &nodes[node->AboveChild()];
					secondChild = node + 1;
				}

				// Advance to next child node, possibly enqueue other child
				if (tPlane > tMax || tPlane <= 0)
					node = firstChild;
				else if (tPlane < tMin)
					node = secondChild;
				else {
					// Enqueue _secondChild_ in todo list
					todo[todoPos].node = secondChild;
					todo[todoPos].tMin = tPlane;
					todo[todoPos].tMax = tMax;
					++todoPos;
					node = firstChild;
					tMax = tPlane;
				}
			}
		}
		return false;
	}

	std::shared_ptr<KdTreeAccel> CreateKdTreeAccelerator(
		std::vector<std::shared_ptr<Primitive>> prims, const ParamSet & ps) {
		int isectCost = ps.FindOneInt("intersectcost", 80);
		int travCost = ps.FindOneInt("traversalcost", 1);
		Float emptyBonus = ps.FindOneFloat("emptybonus", 0.5f);
		int maxPrims = ps.FindOneInt("maxprims", 1);
		int maxDepth = ps.FindOneInt("maxdepth", -1);
		return std::make_shared<KdTreeAccel>(std::move(prims), isectCost, travCost, emptyBonus,
			maxPrims, maxDepth);
	}

} // namespace pbrt