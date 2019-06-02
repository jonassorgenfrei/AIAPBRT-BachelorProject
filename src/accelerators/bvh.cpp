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


 // accelerators/bvh.cpp*
#include "accelerators/bvh.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "parallel.h"
#include <algorithm>

namespace pbrt {

	STAT_MEMORY_COUNTER("Memory/BVH tree", treeBytes);
	STAT_RATIO("BVH/Primitives per leaf node", totalPrimitives, totalLeafNodes);
	STAT_COUNTER("BVH/Interior nodes", interiorNodes);
	STAT_COUNTER("BVH/Leaf nodes", leafNodes);

	// BVHAccel Local Declarations
	struct BVHPrimitiveInfo {
		/// <summary>
		/// Initializes a new instance of the <see cref="BVHPrimitiveInfo"/> struct.
		/// </summary>
		BVHPrimitiveInfo() {}
		
		/// <summary>
		/// Initializes a new instance of the <see cref="BVHPrimitiveInfo"/> struct.
		/// </summary>
		/// <param name="primitiveNumber">The primitive number.</param>
		/// <param name="bounds">The bounds.</param>
		BVHPrimitiveInfo(size_t primitiveNumber, const Bounds3f& bounds)
			: primitiveNumber(primitiveNumber),
			bounds(bounds),
			centroid(.5f * bounds.pMin + .5f * bounds.pMax) {}

		size_t primitiveNumber;	// index in the primitives array
		Bounds3f bounds;		// complete bounding box
		Point3f centroid;		// centeroid of the bounding box
	};

	struct BVHBuildNode {
		// BVHBuildNode Public Methods
		/// <summary>
		/// Initializes a leaf node.
		/// </summary>
		/// <param name="first">The primitive first index</param>
		/// <param name="n">The number of primitives (from the (including) first index)</param>
		/// <param name="b">The bounding box of the leave</param>
		void InitLeaf(int first, int n, const Bounds3f& b) {
			firstPrimOffset = first;
			nPrimitives = n;
			bounds = b;
			children[0] = children[1] = nullptr;	// no children on a leaf node
			++leafNodes;
			++totalLeafNodes;
			totalPrimitives += n;
		}

		/// <summary>
		/// Initializes the interior node.
		/// </summary>
		/// <param name="axis">The coordinate axis along which primitives were partitioned for distribution to their children.</param>
		/// <param name="c0">Child Node 1</param>
		/// <param name="c1">Child Node 2</param>
		void InitInterior(int axis, BVHBuildNode* c0, BVHBuildNode* c1) {
			children[0] = c0;
			children[1] = c1;
			bounds = Union(c0->bounds, c1->bounds);	// computing the bounds of the interior node using the union function
			splitAxis = axis;
			nPrimitives = 0;
			++interiorNodes;
		}

		Bounds3f bounds;				// represents the bounds of all of the children beneath the node
		BVHBuildNode* children[2];		// pointers to its two children
		int splitAxis,					// coordinate axis along which primitives were partitioned for distribution to their children
			// Primitive Indices [firstPrimOffset, firstPrimOffset+nPrimtives-1]
			firstPrimOffset,			// first primitiv index
			nPrimitives;				// (not including) last primitiv index
	};

	struct MortonPrimitive {
		int primitiveIndex;		// index of the primitive in the primitive Info array
		uint32_t mortonCode;	// Morton code
	};

	/// <summary>
	/// LBVHTreelet representing each primitive cluster 
	/// </summary>
	struct LBVHTreelet {
		int startIndex,				// index in the mortonPrims array of the first primitive in the cluster 
			nPrimitives;			// number of following primitives
		BVHBuildNode* buildNodes;	// pointer to the root of the corresponding LBVH
	};

	struct LinearBVHNode {
		Bounds3f bounds;		// bounding box for the node 
		union {
			int primitivesOffset;   // leaf --> stores offset
			int secondChildOffset;  // interior --> stores offset to the second child
		};
		uint16_t nPrimitives;  // 0 -> interior node
		uint8_t axis;          // interior node: xyz (which of the coordinate axes the primitive were partitioned along when the hierarchy was built)
		uint8_t pad[1];        // ensure 32 byte total size (avoid nodes to straddle cache lines on modern CPU architectures)
	};

	// BVHAccel Utility Functions
	/// <summary>
	/// Shifts the ith bit to be in position 3i. In other words, shifts it 2i places to
	/// the left. All other bits are set to zero.
	/// </summary>
	/// <param name="x">A 32-bit value</param>
	/// <returns>Result of shifting the ith bit to be at the 3ith bit,leaving eros in other bits</returns>
	inline uint32_t LeftShift3(uint32_t x) {
		CHECK_LE(x, (1 << 10));
		if (x == (1 << 10)) --x;
#ifdef PBRT_HAVE_BINARY_CONSTANTS
		x = (x | (x << 16)) & 0b00000011000000000000000011111111;
		// x = ---- --98 ---- ---- ---- ---- 7654 3210
		x = (x | (x << 8)) & 0b00000011000000001111000000001111;
		// x = ---- --98 ---- ---- 7654 ---- ---- 3210
		x = (x | (x << 4)) & 0b00000011000011000011000011000011;
		// x = ---- --98 ---- 76-- --54 ---- 32-- --10
		x = (x | (x << 2)) & 0b00001001001001001001001001001001;
		// x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
#else
		x = (x | (x << 16)) & 0x30000ff;
		// x = ---- --98 ---- ---- ---- ---- 7654 3210
		x = (x | (x << 8)) & 0x300f00f;
		// x = ---- --98 ---- ---- 7654 ---- ---- 3210
		x = (x | (x << 4)) & 0x30c30c3;
		// x = ---- --98 ---- 76-- --54 ---- 32-- --10
		x = (x | (x << 2)) & 0x9249249;
		// x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
#endif // PBRT_HAVE_BINARY_CONSTANTS
		return x;
	}

	/// <summary>
	/// Converts a 3D coordinate value, where each component is a floating-point value between 0 and 2^10 
	/// to integers and then computes the Morton Code by expanding the three 10-bit quantized values
	/// so that their ith bits are at position 3i, then shifting the y bits over one more, the z bits
	/// over two more and ORing together the results.
	/// </summary>
	/// <param name="v">3D coordinate value</param>
	/// <returns></returns>
	inline uint32_t EncodeMorton3(const Vector3f& v) {
		CHECK_GE(v.x, 0);
		CHECK_GE(v.y, 0);
		CHECK_GE(v.z, 0);
		return (LeftShift3(v.z) << 2) | (LeftShift3(v.y) << 1) | LeftShift3(v.x);
	}

	/// <summary>
	/// Sort the Morton index value using a radix sort.
	/// Based on bucketing items based on some key.
	/// Radix sort can be used to sort integer values by sorting them one digit
	/// at a time, going from the rightmost digit to the leftmost.
	/// Noticable faster than using std::sort().
	/// </summary>
	/// <param name="v">Vector of MortonPrimtives</param>
	static void RadixSort(std::vector<MortonPrimitive>* v) {

		std::vector<MortonPrimitive> tempVector(v->size());
		PBRT_CONSTEXPR int bitsPerPass = 6;	// number of bits processed per pass
		PBRT_CONSTEXPR int nBits = 30;		// number of bits
		
		static_assert((nBits % bitsPerPass) == 0,
			"Radix sort bitsPerPass must evenly divide nBits");
		PBRT_CONSTEXPR int nPasses = nBits / bitsPerPass;

		for (int pass = 0; pass < nPasses; ++pass) {
			// Perform one pass of radix sort, sorting _bitsPerPass_ bits
			int lowBit = pass * bitsPerPass;	// starting

			// Set in and out vector pointers for radix sort pass
			// alternates each pass through the loop 
			std::vector<MortonPrimitive>& in = (pass & 1) ? tempVector : *v;
			std::vector<MortonPrimitive>& out = (pass & 1) ? *v : tempVector;

			// Count number of zero bits in array for current radix sort bit
			PBRT_CONSTEXPR int nBuckets = 1 << bitsPerPass;	// 2^n buckets
			int bucketCount[nBuckets] = { 0 };
			PBRT_CONSTEXPR int bitMask = (1 << bitsPerPass) - 1;
			for (const MortonPrimitive& mp : in) {	// count how many values will land in each bucket
				int bucket = (mp.mortonCode >> lowBit) & bitMask;	// to determine the index of the current value, shift the index so that the bit at index lowBit is at bit 0
																	// then mask off the low bitsPerPass bits
				CHECK_GE(bucket, 0);
				CHECK_LT(bucket, nBuckets);
				++bucketCount[bucket];
			}

			// Compute starting index in output array for each bucket
			int outIndex[nBuckets];
			outIndex[0] = 0;
			for (int i = 1; i < nBuckets; ++i)	// compute the offset in the output array where each bucket's values start
				outIndex[i] = outIndex[i - 1] + bucketCount[i - 1];	// sum of how many values land in the preceding buckets

			// Store sorted values in output array
			for (const MortonPrimitive& mp : in) {	// loop over the MortonPrimitives
				int bucket = (mp.mortonCode >> lowBit) & bitMask;	//  recompute the bucket that each one lands in
				out[outIndex[bucket]++] = mp;	// store MortonPrimitives in the output array
			}
		}
		// Copy final result from _tempVector_, if needed
		if (nPasses & 1)	//if an odd number of radix sort passes were performed
			std::swap(*v, tempVector);	// copy from the temp vector to the output vector that was originally passed to RadixSort()
	}

	// BVHAccel Method Definitions
	BVHAccel::BVHAccel(std::vector<std::shared_ptr<Primitive>> p,
		int maxPrimsInNode, SplitMethod splitMethod)
		: maxPrimsInNode(std::min(255, maxPrimsInNode)),
		splitMethod(splitMethod),
		primitives(std::move(p)) {

		ProfilePhase _(Prof::AccelConstruction);
		if (primitives.empty()) return;
		// Build BVH from _primitives_

		// Initialize _primitiveInfo_ array for primitives
		std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
		for (size_t i = 0; i < primitives.size(); ++i)
			primitiveInfo[i] = { i, primitives[i]->WorldBound() };	// compute bounding information about each primitive

		// Build BVH tree for primitives using _primitiveInfo_
		MemoryArena arena(1024 * 1024);	// memory area for memory management
		int totalNodes = 0;				// total number of nodes created
		std::vector<std::shared_ptr<Primitive>> orderedPrims;	// stores the primitives ordered
																// so that the primitives in
																// leaf nodes occupy contiguous 
																// ranges in the array
		orderedPrims.reserve(primitives.size());

		BVHBuildNode* root;	// Pointer to the Tree-Root Node

		if (splitMethod == SplitMethod::HLBVH)	// HLBVH-algorithm to build the tree
			root = HLBVHBuild(arena, primitiveInfo, &totalNodes, orderedPrims);
		else
			root = recursiveBuild(arena, primitiveInfo, 0, primitives.size(),
				&totalNodes, orderedPrims);

		primitives.swap(orderedPrims);	// swap original primitives with the ordered primitives

		primitiveInfo.resize(0);
		LOG(INFO) << StringPrintf("BVH created with %d nodes for %d "
			"primitives (%.2f MB), arena allocated %.2f MB",
			totalNodes, (int)primitives.size(),
			float(totalNodes * sizeof(LinearBVHNode)) /
			(1024.f * 1024.f),
			float(arena.TotalAllocated()) /
			(1024.f * 1024.f));

		// Compute representation of depth-first traversal of BVH tree
		treeBytes += totalNodes * sizeof(LinearBVHNode) + sizeof(*this) +
			primitives.size() * sizeof(primitives[0]);
		nodes = AllocAligned<LinearBVHNode>(totalNodes);
		int offset = 0;
		flattenBVHTree(root, &offset);	// store BVH in linear order (depth-first-order)
		CHECK_EQ(totalNodes, offset);
	}

	Bounds3f BVHAccel::WorldBound() const {
		return nodes ? nodes[0].bounds : Bounds3f();
	}

	/// <summary>
	/// Struct for Heuristic Buckets for SAH
	/// </summary>
	struct BucketInfo {
		int count = 0;		// count of primitives in the bucket
		Bounds3f bounds;	// bounding box of the bucket
	};

	BVHBuildNode* BVHAccel::recursiveBuild(
		MemoryArena& arena, std::vector<BVHPrimitiveInfo>& primitiveInfo, int start,
		int end, int* totalNodes,
		std::vector<std::shared_ptr<Primitive>>& orderedPrims) {

		CHECK_NE(start, end);
		BVHBuildNode* node = arena.Alloc<BVHBuildNode>();
		(*totalNodes)++;

		// Compute bounds of all primitives in BVH node
		Bounds3f bounds;
		for (int i = start; i < end; ++i)
			bounds = Union(bounds, primitiveInfo[i].bounds);
		
		int nPrimitives = end - start;

		if (nPrimitives == 1) {
			// Create leaf _BVHBuildNode_
			int firstPrimOffset = orderedPrims.size();
			for (int i = start; i < end; ++i) {		// the primitives overlapping the leaf are appended to the orderedPrims array
				int primNum = primitiveInfo[i].primitiveNumber;
				orderedPrims.push_back(primitives[primNum]);
			}
			node->InitLeaf(firstPrimOffset, nPrimitives, bounds);	// initalize a leaf node object
			return node;
		}
		else {	// interior nodes
			// Compute bound of primitive centroids, choose split dimension _dim_
			Bounds3f centroidBounds;
			for (int i = start; i < end; ++i)
				centroidBounds = Union(centroidBounds, primitiveInfo[i].centroid);
			int dim = centroidBounds.MaximumExtent();	// find axis with probably minimum overlapping

			// Partition primitives into two sets and build children
			int mid = (start + end) / 2;
			if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {	// check if all centeroid points are at the same position
				// Create leaf _BVHBuildNode_
				int firstPrimOffset = orderedPrims.size();
				for (int i = start; i < end; ++i) {
					int primNum = primitiveInfo[i].primitiveNumber;
					orderedPrims.push_back(primitives[primNum]);
				}
				node->InitLeaf(firstPrimOffset, nPrimitives, bounds);	// create a new leaf node
				return node;
			}
			else {
				// Partition primitives based on _splitMethod_
				switch (splitMethod) {
					case SplitMethod::Middle: {
						// Partition primitives through node's midpoint
						Float pmid =
							(centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;	// computes the midpoint of the primtives' centroids along the splitting axis
						
						// classify primitives into two sets, depending on whether their 
						// centroids are above or below the midpoint.
						BVHPrimitiveInfo* midPtr = std::partition(	// std function to order elements retunring true/false
														&primitiveInfo[start], 
														&primitiveInfo[end - 1] + 1,
														[dim, pmid](const BVHPrimitiveInfo& pi) {	// comparison function 
															return pi.centroid[dim] < pmid;
														});
						mid = midPtr - &primitiveInfo[0];
						// For lots of prims with large overlapping bounding boxes, this
						// may fail to partition; in that case don't break and fall
						// through
						// to EqualCounts.
						if (mid != start && mid != end)
							break;
					}
					case SplitMethod::EqualCounts: {
						// Partition primitives into equally-sized subsets
						mid = (start + end) / 2;	// n /2
						// split the primitives into two equal-sized subsets
						std::nth_element(	&primitiveInfo[start],			// start pointer
											&primitiveInfo[mid],			// middle pointer
											&primitiveInfo[end - 1] + 1,	// ending pointer
											[dim](const BVHPrimitiveInfo& a,	
												const BVHPrimitiveInfo& b) {	// comparing function 
												return a.centroid[dim] < b.centroid[dim];
											});
						break;
					}
					case SplitMethod::SAH:
					default: {	
						// Partition primitives using approximate SAH (Surface Area Heuristic)
						if (nPrimitives <= 2) {
							// Partition primitives into equally-sized subsets
							mid = (start + end) / 2;
							std::nth_element(&primitiveInfo[start],
											&primitiveInfo[mid],
											&primitiveInfo[end - 1] + 1,
											[dim](const BVHPrimitiveInfo& a,
												const BVHPrimitiveInfo& b) {
													return a.centroid[dim] < b.centroid[dim];
											});
						}
						else {
							// Allocate _BucketInfo_ for SAH partition buckets
							PBRT_CONSTEXPR int nBuckets = 12;
							BucketInfo buckets[nBuckets];	// containing the devided range along the axis 

							// Initialize _BucketInfo_ for SAH partition buckets
							for (int i = start; i < end; ++i) {	// for each primitive
								// determine the bucket that the centroid lies in
								int b = nBuckets *
									centroidBounds.Offset(primitiveInfo[i].centroid)[dim];

								if (b == nBuckets) 
									b = nBuckets - 1;
								CHECK_GE(b, 0);
								CHECK_LT(b, nBuckets);
								buckets[b].count++;	// increase count of primtives
								// update bucket bounds to include the primitive's bounds
								buckets[b].bounds =	Union(buckets[b].bounds, primitiveInfo[i].bounds);
							}

							// Compute costs for splitting after each bucket
							Float cost[nBuckets - 1];	// array of cost for splitting at each of the bucket boundaries
							for (int i = 0; i < nBuckets - 1; ++i) {	// loops over all buckets and and compute the cost
								// NOTE: doesn't split after the last buckt!
								Bounds3f b0, b1;
								int count0 = 0, count1 = 0;
								// TODO: current -> computation here has O(n^2) | linear time implementation based on forward scan 
								// over the buckets and backwards scan over the buckets that incrementally compute and store bounds and counts
								
								for (int j = 0; j <= i; ++j) { // sum buckets till current bucket
									b0 = Union(b0, buckets[j].bounds);
									count0 += buckets[j].count;
								}
								for (int j = i + 1; j < nBuckets; ++j) {	// sum buckets from current to last
									b1 = Union(b1, buckets[j].bounds);
									count1 += buckets[j].count;
								}
								// arbitrarily set estimated traversal cost and estimated
								// intersection cost to 1
								cost[i] = 1 +
									(count0 * b0.SurfaceArea() +
										count1 * b1.SurfaceArea()) /
										bounds.SurfaceArea();	
								// c(A,B) = ttrav + pA * sigma[1,NA](tisect(ai)) + pB * sigma[1,NB](tisect(bi))
								// probabilistic that a ray passing through A also passes through B and C are given by (SB+SC)/SA
							}

							// Find bucket to split at that minimizes SAH metric (linear scan)
							Float minCost = cost[0];
							int minCostSplitBucket = 0;	// find partition with minimum cost
							for (int i = 1; i < nBuckets - 1; ++i) {
								if (cost[i] < minCost) {
									minCost = cost[i];
									minCostSplitBucket = i;
								}
							}

							// Either create leaf or split primitives at selected SAH
							// bucket
							Float leafCost = nPrimitives;	// because the estimated intersection cost was arbitrarily set to 1
							
							// check if number of primitives is less than max number of primitives allowed in a node
							// or estimated bucket cost is lower than building a node with the exisiting primitives 
							if (nPrimitives > maxPrimsInNode || minCost < leafCost) {	
								// reordering nodes 
								BVHPrimitiveInfo* pmid = std::partition(	&primitiveInfo[start], 
																			&primitiveInfo[end - 1] + 1,
																			[=](const BVHPrimitiveInfo& pi) {	// comparison function
																				int b = nBuckets *
																					centroidBounds.Offset(pi.centroid)[dim];
																				if (b == nBuckets) 
																					b = nBuckets - 1;
																				CHECK_GE(b, 0);
																				CHECK_LT(b, nBuckets);
																				return b <= minCostSplitBucket;
																			});
								mid = pmid - &primitiveInfo[0];
							}
							else {
								// Create leaf _BVHBuildNode_
								int firstPrimOffset = orderedPrims.size();
								for (int i = start; i < end; ++i) {
									int primNum = primitiveInfo[i].primitiveNumber;
									orderedPrims.push_back(primitives[primNum]);
								}
								node->InitLeaf(firstPrimOffset, nPrimitives, bounds);	// init new leaf
								return node;
							}
						}
						break;
					}
				}
				node->InitInterior(dim,
					recursiveBuild(arena, primitiveInfo, start, mid,
						totalNodes, orderedPrims),
					recursiveBuild(arena, primitiveInfo, mid, end,
						totalNodes, orderedPrims));
			}
		}
		return node;
	}

	BVHBuildNode* BVHAccel::HLBVHBuild(	MemoryArena& arena, 
										const std::vector<BVHPrimitiveInfo>& primitiveInfo,
										int* totalNodes,
										std::vector<std::shared_ptr<Primitive>>& orderedPrims) const {

		// Compute bounding box of all primitive centroids
		Bounds3f bounds;	// bind the centroids of all primitives
		for (const BVHPrimitiveInfo& pi : primitiveInfo)
			bounds = Union(bounds, pi.centroid);

		// Compute Morton indices of primitives
		std::vector<MortonPrimitive> mortonPrims(primitiveInfo.size());
		ParallelFor([&](int i) {
			// Initialize _mortonPrims[i]_ for _i_th primitive
			PBRT_CONSTEXPR int mortonBits = 10;	// allows to fit into a single 32-bit variable
			PBRT_CONSTEXPR int mortonScale = 1 << mortonBits;
			mortonPrims[i].primitiveIndex = primitiveInfo[i].primitiveNumber;
			// quantize centroid positions with respect to the overall bounds
			Vector3f centroidOffset = bounds.Offset(primitiveInfo[i].centroid);	
			// scale floating-point centroid offsets inside the bounding box by 2^10
			mortonPrims[i].mortonCode = EncodeMorton3(centroidOffset * mortonScale);
			}, primitiveInfo.size(), 512);	// loop chunck size of 512 ==> worker threads are given groups of 512 primitives

		// Radix sort primitive Morton indices
		RadixSort(&mortonPrims);

		// Create LBVH treelets at bottom of BVH

		// Find intervals of primitives for each treelet
		std::vector<LBVHTreelet> treeletsToBuild;
		for (int start = 0, end = 1; end <= (int)mortonPrims.size(); ++end) {
#ifdef PBRT_HAVE_BINARY_CONSTANTS
			uint32_t mask = 0b00111111111111000000000000000000;
#else
			uint32_t mask = 0x3ffc0000;
#endif
			// find sets of primitives that have the same values for the 12 high bits of their 30-bit Morton codes
			if (end == (int)mortonPrims.size() ||
				((mortonPrims[start].mortonCode & mask) !=
				(mortonPrims[end].mortonCode & mask))) {	
				// Add entry to _treeletsToBuild_ for this treelet
				int nPrimitives = end - start;
				int maxBVHNodes = 2 * nPrimitives;
				BVHBuildNode* nodes = arena.Alloc<BVHBuildNode>(maxBVHNodes, false);
					// false - don't execute the constructors of the underlying objects because of significant overhead
				treeletsToBuild.push_back({ start, nPrimitives, nodes });
				start = end;
			}
		}

		// Create LBVHs for treelets in parallel
		std::atomic<int> atomicTotal(0),			// total number of nodes in all of the LBVHs
						orderedPrimsOffset(0);		// index of the next available entry in orderedPrims
		orderedPrims.resize(primitives.size());
		ParallelFor([&](int i) {
			// Generate _i_th LBVH treelet
			int nodesCreated = 0;
			const int firstBitIndex = 29 - 12;	// previously used high 12 bits to cluster the primitives
			LBVHTreelet& tr = treeletsToBuild[i];
			tr.buildNodes =
				emitLBVH(tr.buildNodes, primitiveInfo, &mortonPrims[tr.startIndex],
					tr.nPrimitives, &nodesCreated, orderedPrims,
					&orderedPrimsOffset, firstBitIndex);
			atomicTotal += nodesCreated;	// update atomic variable once per treelet when each treelet is done
			}, treeletsToBuild.size());
		*totalNodes = atomicTotal;

		// Create and return SAH BVH from LBVH treelets
		std::vector<BVHBuildNode*> finishedTreelets;
		finishedTreelets.reserve(treeletsToBuild.size());
		for (LBVHTreelet& treelet : treeletsToBuild)	// add all treelets to the finsihedTreelets vector
			finishedTreelets.push_back(treelet.buildNodes);
		return buildUpperSAH(arena, finishedTreelets, 0, finishedTreelets.size(),
			totalNodes); // create a BVH of all treelets
	}

	BVHBuildNode* BVHAccel::emitLBVH(
		BVHBuildNode*& buildNodes,
		const std::vector<BVHPrimitiveInfo>& primitiveInfo,
		MortonPrimitive* mortonPrims, int nPrimitives, int* totalNodes,
		std::vector<std::shared_ptr<Primitive>>& orderedPrims,
		std::atomic<int>* orderedPrimsOffset, int bitIndex) const {

		CHECK_GT(nPrimitives, 0);
		if (bitIndex == -1 || nPrimitives < maxPrimsInNode) {	// check if partitioned the primitive with the final low bit or if it's down a small number of primitives
			// Create and return leaf node of LBVH treelet
			(*totalNodes)++;
			BVHBuildNode* node = buildNodes++;
			Bounds3f bounds;
			int firstPrimOffset = orderedPrimsOffset->fetch_add(nPrimitives);	// atomically add the value of nPrimitives to orderedPrimsOffset and returns its olf value before the addition
			for (int i = 0; i < nPrimitives; ++i) {
				int primitiveIndex = mortonPrims[i].primitiveIndex;
				orderedPrims[firstPrimOffset + i] = primitives[primitiveIndex];
				bounds = Union(bounds, primitiveInfo[primitiveIndex].bounds);	// create bounding box around Primitives
			}
			node->InitLeaf(firstPrimOffset, nPrimitives, bounds);	// create new Leaf
			return node;
		}
		else {
			int mask = 1 << bitIndex;
			// Advance to next subtree level if there's no LBVH split for this bit
			if ((mortonPrims[0].mortonCode & mask) ==
				(mortonPrims[nPrimitives - 1].mortonCode & mask)) { // check if first and last primitive in the range both have the same bit value for this plane
				// if all of the primitives lie on the same side of the splitting plane
				// --> proceed to the next bit without unnecessarily creating a node.
				return emitLBVH(buildNodes, primitiveInfo, mortonPrims, nPrimitives,
					totalNodes, orderedPrims, orderedPrimsOffset,
					bitIndex - 1);
			}

			// Find LBVH split point for this dimension
			// using binary search
			int searchStart = 0, searchEnd = nPrimitives - 1;
			while (searchStart + 1 != searchEnd) {
				CHECK_NE(searchStart, searchEnd);
				int mid = (searchStart + searchEnd) / 2;
				// check if the bitIndexth bit goes from 0 to 1 in the current set of primitives
				if ((mortonPrims[searchStart].mortonCode & mask) ==
					(mortonPrims[mid].mortonCode & mask))
					searchStart = mid;
				else {
					CHECK_EQ(mortonPrims[mid].mortonCode & mask,
						mortonPrims[searchEnd].mortonCode & mask);
					searchEnd = mid;
				}
			}
			int splitOffset = searchEnd;	// splitting offset

			CHECK_LE(splitOffset, nPrimitives - 1);
			CHECK_NE(mortonPrims[splitOffset - 1].mortonCode & mask,
				mortonPrims[splitOffset].mortonCode & mask);

			// Create and return interior LBVH node
			(*totalNodes)++;
			// claim a node to use as an interior node
			BVHBuildNode* node = buildNodes++;
			// recursively build LBVHs for both partitioned sets of primitives
			BVHBuildNode* lbvh[2] = {
				emitLBVH(buildNodes, primitiveInfo, mortonPrims, splitOffset,
						 totalNodes, orderedPrims, orderedPrimsOffset,
						 bitIndex - 1),
				emitLBVH(buildNodes, primitiveInfo, &mortonPrims[splitOffset],
						 nPrimitives - splitOffset, totalNodes, orderedPrims,
						 orderedPrimsOffset, bitIndex - 1)
			};
			int axis = bitIndex % 3;
			node->InitInterior(axis, lbvh[0], lbvh[1]);
			return node;
		}
	}

	BVHBuildNode* BVHAccel::buildUpperSAH(MemoryArena& arena,
		std::vector<BVHBuildNode*>& treeletRoots,
		int start, int end,
		int* totalNodes) const {
		CHECK_LT(start, end);
		int nNodes = end - start;
		if (nNodes == 1) return treeletRoots[start];
		(*totalNodes)++;
		BVHBuildNode* node = arena.Alloc<BVHBuildNode>();

		// Compute bounds of all nodes under this HLBVH node
		Bounds3f bounds;
		for (int i = start; i < end; ++i)
			bounds = Union(bounds, treeletRoots[i]->bounds);

		// Compute bound of HLBVH node centroids, choose split dimension _dim_
		Bounds3f centroidBounds;
		for (int i = start; i < end; ++i) {
			Point3f centroid =
				(treeletRoots[i]->bounds.pMin + treeletRoots[i]->bounds.pMax) *
				0.5f;
			centroidBounds = Union(centroidBounds, centroid);
		}
		int dim = centroidBounds.MaximumExtent();
		// FIXME: if this hits, what do we need to do?
		// Make sure the SAH split below does something... ?
		CHECK_NE(centroidBounds.pMax[dim], centroidBounds.pMin[dim]);

		// Allocate _BucketInfo_ for SAH partition buckets
		PBRT_CONSTEXPR int nBuckets = 12;
		struct BucketInfo {
			int count = 0;
			Bounds3f bounds;
		};
		BucketInfo buckets[nBuckets];

		// Initialize _BucketInfo_ for HLBVH SAH partition buckets
		for (int i = start; i < end; ++i) {
			Float centroid = (treeletRoots[i]->bounds.pMin[dim] +
				treeletRoots[i]->bounds.pMax[dim]) *
				0.5f;
			int b =
				nBuckets * ((centroid - centroidBounds.pMin[dim]) /
				(centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
			if (b == nBuckets) b = nBuckets - 1;
			CHECK_GE(b, 0);
			CHECK_LT(b, nBuckets);
			buckets[b].count++;
			buckets[b].bounds = Union(buckets[b].bounds, treeletRoots[i]->bounds);
		}

		// Compute costs for splitting after each bucket
		Float cost[nBuckets - 1];
		for (int i = 0; i < nBuckets - 1; ++i) {
			Bounds3f b0, b1;
			int count0 = 0, count1 = 0;
			for (int j = 0; j <= i; ++j) {
				b0 = Union(b0, buckets[j].bounds);
				count0 += buckets[j].count;
			}
			for (int j = i + 1; j < nBuckets; ++j) {
				b1 = Union(b1, buckets[j].bounds);
				count1 += buckets[j].count;
			}
			cost[i] = .125f +
				(count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea()) /
				bounds.SurfaceArea();
		}

		// Find bucket to split at that minimizes SAH metric
		Float minCost = cost[0];
		int minCostSplitBucket = 0;
		for (int i = 1; i < nBuckets - 1; ++i) {
			if (cost[i] < minCost) {
				minCost = cost[i];
				minCostSplitBucket = i;
			}
		}

		// Split nodes and create interior HLBVH SAH node
		BVHBuildNode** pmid = std::partition(
			&treeletRoots[start], &treeletRoots[end - 1] + 1,
			[=](const BVHBuildNode* node) {
				Float centroid =
					(node->bounds.pMin[dim] + node->bounds.pMax[dim]) * 0.5f;
				int b = nBuckets *
					((centroid - centroidBounds.pMin[dim]) /
					(centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
				if (b == nBuckets) b = nBuckets - 1;
				CHECK_GE(b, 0);
				CHECK_LT(b, nBuckets);
				return b <= minCostSplitBucket;
			});
		int mid = pmid - &treeletRoots[0];
		CHECK_GT(mid, start);
		CHECK_LT(mid, end);
		node->InitInterior(
			dim, this->buildUpperSAH(arena, treeletRoots, start, mid, totalNodes),
			this->buildUpperSAH(arena, treeletRoots, mid, end, totalNodes));
		return node;
	}

	int BVHAccel::flattenBVHTree(BVHBuildNode* node, int* offset) {
		LinearBVHNode* linearNode = &nodes[*offset];
		linearNode->bounds = node->bounds;
		int myOffset = (*offset)++; //increase offset 
		if (node->nPrimitives > 0) {	// leaf node
			CHECK(!node->children[0] && !node->children[1]);
			CHECK_LT(node->nPrimitives, 65536);
			linearNode->primitivesOffset = node->firstPrimOffset;
			linearNode->nPrimitives = node->nPrimitives;
		}
		else { // interior node
			// Create interior flattened BVH node
			linearNode->axis = node->splitAxis;
			linearNode->nPrimitives = 0;
			flattenBVHTree(node->children[0], offset); //first child ends up immediately after the current node in the array 
			linearNode->secondChildOffset =
				flattenBVHTree(node->children[1], offset);
		}
		return myOffset;
	}

	BVHAccel::~BVHAccel() { FreeAligned(nodes); }

	bool BVHAccel::Intersect(const Ray& ray, SurfaceInteraction* isect) const {

		if (!nodes)			// return false if nodes are empty
			return false;
		
		ProfilePhase p(Prof::AccelIntersect);
		bool hit = false;
		// precompute some values related to the ray that will be used repeatedly
		Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
		int dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };

		int toVisitOffset = 0,		// offset to the next free element in the nodesToVisit Stack
			currentNodeIndex = 0;	// holds the offset into the nodes array of the node to be visited
		int nodesToVisit[64];		// nodes that still need to be visited	(like a stack)

		// Follow ray through BVH nodes to find primitive intersections	
		while (true) {
			const LinearBVHNode* node = &nodes[currentNodeIndex];	// current node
			// Check ray against BVH node
			if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {	// check if the ray intersects the node's bounding box (or starts inside of it)
				if (node->nPrimitives > 0) {	// leaf node
					// Intersect ray with primitives in leaf BVH node
					for (int i = 0; i < node->nPrimitives; ++i)
						if (primitives[node->primitivesOffset + i]->Intersect(ray, isect))
							hit = true;	// even if an intersection is found, the remaining nodes have to be visited in case one of them yields a closer intersection
					if (toVisitOffset == 0) 
						break;	// end if stack is empty
					currentNodeIndex = nodesToVisit[--toVisitOffset];	// pop node from stack
				}
				else {	// interior node
					// Put far BVH node on _nodesToVisit_ stack, advance to near
					// node; Using the sign of the ray's direction vector for the coordinate axis along which primitives were partitioned for the current node 
					if (dirIsNeg[node->axis]) {		// if sign is negativ: visit 2nd child first; since the primitives that went into the second child's subtree were on the upper side of the partition points
						nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
						currentNodeIndex = node->secondChildOffset;
					} else {	// if sign is positiv: visit 1st child first
						nodesToVisit[toVisitOffset++] = node->secondChildOffset;
						currentNodeIndex = currentNodeIndex + 1;		// depth-first layout: first child is immediately after the current node
					}
				}
			}
			else {	// no intersection found
				if (toVisitOffset == 0) break;	// end if stack is empty
				currentNodeIndex = nodesToVisit[--toVisitOffset];	// pop node from stack
			}
		}
		return hit;
	}

	bool BVHAccel::IntersectP(const Ray& ray) const {

		if (!nodes)			// return false if nodes are empty
			return false;

		ProfilePhase p(Prof::AccelIntersectP);

		// precompute some values related to the ray that will be used repeatedly
		Vector3f invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
		int dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };

		int nodesToVisit[64];		// nodes that still need to be visited	(like a stack)
		int toVisitOffset = 0,		// offset to the next free element in the nodesToVisit Stack
			currentNodeIndex = 0;	// holds the offset into the nodes array of the node to be visited

		// Follow ray through BVH nodes to find primitive intersections	
		while (true) {
			const LinearBVHNode* node = &nodes[currentNodeIndex];		// current node

			if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {	// check if the ray intersects the node's bounding box (or starts inside of it)
				// Process BVH node _node_ for traversal
				if (node->nPrimitives > 0) {	// leaf node
					// Intersect ray with primitives in leaf BVH node
					for (int i = 0; i < node->nPrimitives; ++i) {
						if (primitives[node->primitivesOffset + i]->IntersectP(ray)) {
							return true;	// traversal stops immediately when any intersection is founds
						}
					}
					if (toVisitOffset == 0) break;	// end if stack is empty
					currentNodeIndex = nodesToVisit[--toVisitOffset];	// pop node from stack
				}
				else {	// interior node
					// Put far BVH node on _nodesToVisit_ stack, advance to near
					// node; Using the sign of the ray's direction vector for the coordinate axis along which primitives were partitioned for the current node 
					if (dirIsNeg[node->axis]) {		// if sign is negativ: visit 2nd child first; since the primitives that went into the second child's subtree were on the upper side of the partition points
						/// second child first
						nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
						currentNodeIndex = node->secondChildOffset;
					}
					else {	// no intersection found
						nodesToVisit[toVisitOffset++] = node->secondChildOffset;
						currentNodeIndex = currentNodeIndex + 1;
					}
				}
			}
			else {
				if (toVisitOffset == 0) break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			}
		}
		return false;
	}

	std::shared_ptr<BVHAccel> CreateBVHAccelerator(
		std::vector<std::shared_ptr<Primitive>> prims, const ParamSet& ps) {

		std::string splitMethodName = ps.FindOneString("splitmethod", "sah");
		BVHAccel::SplitMethod splitMethod;
		if (splitMethodName == "sah")
			splitMethod = BVHAccel::SplitMethod::SAH;
		else if (splitMethodName == "hlbvh")
			splitMethod = BVHAccel::SplitMethod::HLBVH;
		else if (splitMethodName == "middle")
			splitMethod = BVHAccel::SplitMethod::Middle;
		else if (splitMethodName == "equal")
			splitMethod = BVHAccel::SplitMethod::EqualCounts;
		else {
			Warning("BVH split method \"%s\" unknown.  Using \"sah\".",
				splitMethodName.c_str());
			splitMethod = BVHAccel::SplitMethod::SAH;
		}

		int maxPrimsInNode = ps.FindOneInt("maxnodeprims", 4);
		return std::make_shared<BVHAccel>(std::move(prims), maxPrimsInNode, splitMethod);
	}

} // namespace pbrt