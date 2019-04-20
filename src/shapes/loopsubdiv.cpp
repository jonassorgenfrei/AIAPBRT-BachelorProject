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

 // shapes/loopsubdiv.cpp*
#include "shapes/loopsubdiv.h"
#include "shapes/triangle.h"
#include "paramset.h"
#include <set>
#include <map>

namespace pbrt {

	struct SDFace;
	struct SDVertex;

	// LoopSubdiv Macros
	// Determine the vertex and face indices before or adter a particular index
	// offsetsand copmute the result modulo three to handling cycling around
	#define NEXT(i) (((i) + 1) % 3)	
	#define PREV(i) (((i) + 2) % 3)	

	// LoopSubdiv Local Structures
	/// <summary>
	/// Subdivision Surface Vertex
	/// </summary>
	struct SDVertex {
		// SDVertex Constructor
		SDVertex(const Point3f& p = Point3f(0, 0, 0)) : p(p) {}

		// SDVertex Methods
		int valence();
		void oneRing(Point3f* p);
		//SDVertex Data
		Point3f p;						// Position
		SDFace* startFace = nullptr;	// pointer to an arbitrary face which is adjacent to it (gives a starting point for finding all of the adjacent faces to it)
		SDVertex* child = nullptr;		// pointer to store the corresponding SDVertex for the next level of subdivision (if any)
		bool regular = false,			// indicates whether it is a regular or an extraordinary vertex
			boundary = false;			// records if it lies on the boundary of the mesh. 
	};

	/// <summary>
	/// Subdivision Surface Face
	/// Maintaines most of the topological structure about the mesh
	/// </summary>
	struct SDFace {
		// SDFace Constructor
		SDFace() {
			for (int i = 0; i < 3; ++i) {
				v[i] = nullptr;
				f[i] = nullptr;
			}
			for (int i = 0; i < 4; ++i) children[i] = nullptr;
		}

		// SDFace Methods
		/// <summary>
		/// Finds the index of a given vertex pointer.
		/// It is a fatal error to pass a pointer to a vertex that isn't 
		/// part of the current face - this case would represent a bug elsewhere
		/// in the subdivision code.
		/// </summary>
		/// <param name="vert">Vertex Pointer</param>
		/// <returns></returns>
		int vnum(SDVertex* vert) const {
			for (int i = 0; i < 3; ++i)
				if (v[i] == vert) return i;
			LOG(FATAL) << "Basic logic error in SDFace::vnum()";
			return -1;
		}


		/// <summary>
		/// The next face is over the ith edge.
		/// </summary>
		/// <param name="vert">The vert.</param>
		/// <returns></returns>
		SDFace* nextFace(SDVertex* vert) { return f[vnum(vert)]; }
		
		/// <summary>
		/// The previous face is across the ede from PREV(i) to i
		/// </summary>
		/// <param name="vert">The vert.</param>
		/// <returns></returns>
		SDFace* prevFace(SDVertex* vert) { return f[PREV(vnum(vert))]; }

		/// <summary>
		/// Returns the next vertice
		/// </summary>
		/// <param name="vert">The vert.</param>
		/// <returns></returns>
		SDVertex* nextVert(SDVertex* vert) { return v[NEXT(vnum(vert))]; }
		/// <summary>
		/// Returns the previous vertice.
		/// </summary>
		/// <param name="vert">The vert.</param>
		/// <returns></returns>
		SDVertex* prevVert(SDVertex* vert) { return v[PREV(vnum(vert))]; }
		
		/// <summary>
		/// Returns the vertex opposite a given edge of a face.
		/// </summary>
		/// <param name="v0">The v0.</param>
		/// <param name="v1">The v1.</param>
		/// <returns></returns>
		SDVertex* otherVert(SDVertex* v0, SDVertex* v1) {
			for (int i = 0; i < 3; ++i)
				if (v[i] != v0 && v[i] != v1) return v[i];
			LOG(FATAL) << "Basic logic error in SDVertex::otherVert()";
			return nullptr;
		}
		SDVertex * v[3];	// pointers to their 3 vertices
		SDFace * f[3];		// pointers to their adjacent faces
		SDFace * children[4];	// pointers to child faces at the next level of subdivision
	};

	struct SDEdge {
		// SDEdge Constructor
		/// <summary>
		/// Initializes a new instance of the <see cref="SDEdge"/> struct.
		/// Takes pointers to the two vertices at each end of the edge.
		/// </summary>
		/// <param name="v0">The v0.</param>
		/// <param name="v1">The v1.</param>
		SDEdge(SDVertex* v0 = nullptr, SDVertex* v1 = nullptr) {
			v[0] = std::min(v0, v1);	// orders them, so v[0] is the first in memory
			v[1] = std::max(v0, v1);
			f[0] = f[1] = nullptr;
			f0edgeNum = -1;
		}

		// SDEdge Comparison Function
		/// <summary>
		/// Ordering operation for SDEdge objects so that they can be stored in
		/// other data structures that rely on ordering being well defined.
		/// </summary>
		/// <param name="e2">The e2.</param>
		/// <returns></returns>
		bool operator<(const SDEdge& e2) const {
			if (v[0] == e2.v[0]) return v[1] < e2.v[1];
			return v[0] < e2.v[0];
		}

		SDVertex* v[2];
		SDFace* f[2];
		int f0edgeNum;
	};

	// LoopSubdiv Local Declarations
	static Point3f weightOneRing(SDVertex* vert, Float beta);
	static Point3f weightBoundary(SDVertex* vert, Float beta);

	// LoopSubdiv Inline Functions
	/// <summary>
	/// Returns the valence of the vertex.
	/// </summary>
	/// <returns></returns>
	inline int SDVertex::valence() {
		SDFace* f = startFace;
		if (!boundary) {
			// Compute valence of interior (non-boundary) vertex
			int nf = 1; 
			while ((f = f->nextFace(this)) != startFace) ++nf;
			// copmute number of adjacent faces starting by following each face's
			// neighbor pointers around the vertex until it reaches the starting
			// face
			return nf;	// the valence is equal to the number of faces visited
		}
		else {
			// Compute valence of boundary vertex
			int nf = 1;
			while ((f = f->nextFace(this)) != nullptr) ++nf;
			// follows pointer to the next face around until it reaches the boundary
			f = startFace;
			while ((f = f->prevFace(this)) != nullptr) ++nf;
			// follows pointers to the previous face around until it reaches the boundary
			return nf + 1; // the valence is one more than the number of adjacent faces
		}
	}

	/// <summary>
	/// Computes a beta value based on the vertex's valence, that ensures smoothness.
	/// </summary>
	/// <param name="valence">The valence.</param>
	/// <returns></returns>
	inline Float beta(int valence) {
		if (valence == 3)
			return 3.f / 16.f;
		else
			return 3.f / (8.f * valence);
	}

	/// <summary>
	/// Computes appropriate vertex weights based on the valence of the vertex.
	/// </summary>
	/// <param name="valence">The valence.</param>
	/// <returns></returns>
	inline Float loopGamma(int valence) {
		return 1.f / (valence + 3.f / (8.f * beta(valence)));
	}

	// LoopSubdiv Function Definitions
	/// <summary>
	/// Applies the subdivision to a mesh represented by collection of vertices
	/// and vertex indices.
	/// The control mesh has to be manifold - no more than 2 faces share any given
	/// edge. It may be closed or open.
	/// Expects that the control mesh is consistently ordered. Each directed edge in the
	/// mesh can be present only once.
	/// </summary>
	/// <param name="ObjectToWorld">The object to world.</param>
	/// <param name="WorldToObject">The world to object.</param>
	/// <param name="reverseOrientation">if set to <c>true</c> [reverse orientation].</param>
	/// <param name="nLevels">The n levels.</param>
	/// <param name="nIndices">The n indices.</param>
	/// <param name="vertexIndices">The vertex indices.</param>
	/// <param name="nVertices">The n vertices.</param>
	/// <param name="p">The p.</param>
	/// <returns>a vector of Triangles that represents the final subdivided mesh</returns>
	static std::vector<std::shared_ptr<Shape>> LoopSubdivide(
		const Transform * ObjectToWorld,
		const Transform * WorldToObject,
		bool reverseOrientation,
		int nLevels,
		int nIndices,
		const int* vertexIndices,
		int nVertices,
		const Point3f * p) {
		
		std::vector<SDVertex*> vertices;
		std::vector<SDFace*> faces;

		// Allocate _LoopSubdiv_ vertices and faces
		std::unique_ptr<SDVertex[]> verts(new SDVertex[nVertices]);
		for (int i = 0; i < nVertices; ++i) {	// allocating one instance of SDVertex for each vertex in the mesh
			verts[i] = SDVertex(p[i]);
			vertices.push_back(&verts[i]);
		}
		int nFaces = nIndices / 3;
		std::unique_ptr<SDFace[]> fs(new SDFace[nFaces]);
		for (int i = 0; i < nFaces; ++i) faces.push_back(&fs[i]);
		//allocate one instance of SDFace for each face.

		// Set face to vertex pointers
		const int* vp = vertexIndices;
		for (int i = 0; i < nFaces; ++i, vp += 3) {	// loops over all faces
			SDFace* f = faces[i];
			for (int j = 0; j < 3; ++j) {
				SDVertex* v = vertices[vp[j]];	// pointer to their vertices
				f->v[j] = v;		//
				v->startFace = f; // pointer to one of the vertex's neighboring faces (reset every time)
			}
		}

		// Set neighbor pointers in _faces_
		std::set<SDEdge> edges;	// makes it possible to search for a particular edge in O(logn)
					// stores the edges that have only one adjacent face so far
		for (int i = 0; i < nFaces; ++i) {	//loops over the faces
			SDFace* f = faces[i];			
			for (int edgeNum = 0; edgeNum < 3; ++edgeNum) { // loop over their three edges
				// Update neighbor pointer for _edgeNum_
				// Updates neighbor pointers as it goes
				int v0 = edgeNum, v1 = NEXT(edgeNum);
				SDEdge e(f->v[v0], f->v[v1]);	// creates an SDEdge object
				if (edges.find(e) == edges.end()) {	// checks if the same edge has been seen previously
					// Handle new edge
					// adds the edge to the set of edges
					e.f[0] = f;				// store current face's pointer
					e.f0edgeNum = edgeNum;	// 
					edges.insert(e);
				}
				else {
					// Handle previously seen edge
					// initialize both faces' neighbor pointers across the edge
					e = *edges.find(e);
					e.f[0]->f[e.f0edgeNum] = f;	// set neighboring face for each face
					f->f[edgeNum] = e.f[0];	// set neighboring face for each face
					edges.erase(e); // remove edge from the edge set
				}
			}
		}

		// Finish vertex initialization
		for (int i = 0; i < nVertices; ++i) {
			SDVertex* v = vertices[i];
			SDFace* f = v->startFace;
			do {
				f = f->nextFace(v);
			} while (f && f != v->startFace); // determine if a vertex is a boundary
											// vertex following next facepointers around
											// the vertex
			v->boundary = (f == nullptr);
			if (!v->boundary && v->valence() == 6)
				v->regular = true; // if valence is 6 for an interior vertex 
			else if (v->boundary && v->valence() == 4)
				v->regular = true;	// if valence is 4 for a boundary vertex
			else
				v->regular = false;	// otherwise it's an extraordinary vertex
		}

		// Refine _LoopSubdiv_ into triangles ==> repeatedly applies the subdivision rules to the mesh, 
		// each time generating a new mesh to be used as the input to the next step
		std::vector<SDFace*> f = faces;		
		std::vector<SDVertex*> v = vertices;
		MemoryArena arena;		// allocate temporary storage through this process.
		for (int i = 0; i < nLevels; ++i) {
			// Update _f_ and _v_ for next level of subdivision
			std::vector<SDFace*> newFaces;
			std::vector<SDVertex*> newVertices;

			// Allocate (storage) next level of children in mesh tree
			for (SDVertex* vertex : v) {	// for vertices
				vertex->child = arena.Alloc<SDVertex>();
				vertex->child->regular = vertex->regular;
				vertex->child->boundary = vertex->boundary;
				newVertices.push_back(vertex->child);
			}
			for (SDFace* face : f) {	// for child faces
				for (int k = 0; k < 4; ++k) {
					face->children[k] = arena.Alloc<SDFace>();
					newFaces.push_back(face->children[k]);
				}
			}

			// Update vertex positions and create new edge vertices

			// Update vertex positions for even vertices
			for (SDVertex* vertex : v) {
				if (!vertex->boundary) {
					// Apply one-ring rule for even vertex
					if (vertex->regular)
						vertex->child->p = weightOneRing(vertex, 1.f / 16.f);
					else
						vertex->child->p =
						weightOneRing(vertex, beta(vertex->valence()));
				}
				else {
					// Apply boundary rule for even vertex
					vertex->child->p = weightBoundary(vertex, 1.f / 8.f);	// note: same weight: 1/8 is used for regular and extraordinary vertices
				}
			}

			// Compute new odd edge vertices
			std::map<SDEdge, SDVertex*> edgeVerts;	// associative array structure that performs efficient lookups
			for (SDFace* face : f) {
				for (int k = 0; k < 3; ++k) {	// loops over each edge of each face in the mesh, computing the new vertex that splits the edge
					// Compute odd vertex on _k_th edge
					SDEdge edge(face->v[k], face->v[NEXT(k)]);	// create an SDEdge-Object for the edge and checked to see if it is in the set of edges that has already been visited
					SDVertex* vert = edgeVerts[edge];		// compute the new vertex on this edge
					if (!vert) {
						// Create and initialize new odd vertex
						vert = arena.Alloc<SDVertex>();
						newVertices.push_back(vert);
						vert->regular = true;		// all (new) vertices will be regular
						vert->boundary = (face->f[k] == nullptr);	// initialize the boundary member
						vert->startFace = face->children[3];	// set new vertex's startFace pointer (child face number three (center) is guaranteed to be adjacent to the new vertex)

						// Apply edge rules to compute new vertex position
						if (vert->boundary) {	// boundary vertices's
							vert->p = 0.5f * edge.v[0]->p;		// just the average of the two adjacent vertices
							vert->p += 0.5f * edge.v[1]->p;
						}
						else {	// interior vertices's
							vert->p = 3.f / 8.f * edge.v[0]->p;	// 2 vertices at the end of the edge with weight 3/8
							vert->p += 3.f / 8.f * edge.v[1]->p;
							vert->p += 1.f / 8.f *				// 2 vertices opposite the edge with weight 1/8
								face->otherVert(edge.v[0], edge.v[1])->p;
							vert->p +=
								1.f / 8.f *
								face->f[k]->otherVert(edge.v[0], edge.v[1])->p;
						}
						edgeVerts[edge] = vert;
					}
				}
			}

			// Update new mesh topology

			// Update even vertex face pointers
			for (SDVertex* vertex : v) {	// loop through all the parent vertices in the mesh
				int vertNum = vertex->startFace->vnum(vertex);						// find its vertex indix in its startface
				vertex->child->startFace = vertex->startFace->children[vertNum];	// use index to find the child face adjacent to the new even vertex
			}

			// Update face neighbor pointers
			for (SDFace* face : f) {
				for (int j = 0; j < 3; ++j) {
					// Update children _f_ pointers for siblings
					face->children[3]->f[j] = face->children[NEXT(j)];		// the k+1st child face (k = 0,1,2) is across the kth edge of the interior face
																			// the interior face is across the k+1st edge of the kth face
					face->children[j]->f[NEXT(j)] = face->children[3];		// the interior face is always sotred in children[3]

					// Update children _f_ pointers for neighbor children
					SDFace* f2 = face->f[j];	// find the neighbor parent across that edge
					face->children[j]->f[j] =		//set the kth edges of the ith child
						f2 ? f2->children[f2->vnum(face->v[j])] : nullptr;	// check if f2 is NOT a boundary face
					f2 = face->f[PREV(j)];		// set the PREV(k)th of the ith child
					face->children[j]->f[PREV(j)] =
						f2 ? f2->children[f2->vnum(face->v[j])] : nullptr;
				}
			}

			// Update face vertex pointers
			for (SDFace* face : f) {
				for (int j = 0; j < 3; ++j) {
					// Update child vertex pointer to new even vertex
					// for the kth child face the kth vertex corresponds to the even vertex that is adjacent to the child face

					face->children[j]->v[j] = face->v[j]->child;	// find the even vertex for the non interior child faces[1 even; 2 odd]  (the interior faces have three odd vertices)

					// Update child vertex pointer to new odd vertex
					SDVertex* vert =
						edgeVerts[SDEdge(face->v[j], face->v[NEXT(j)])];		// to find the off vertex for each split edge of the parent face
					// three child faces have that vertex as an incident vertex
					face->children[j]->v[NEXT(j)] = vert;
					face->children[NEXT(j)]->v[j] = vert;
					face->children[3]->v[j] = vert;
				}
			}

			// Prepare for next level of subdivision
			// point to the next level of subdivision
			f = newFaces;	// move newly created vertices and faces into the v and f arrays
			v = newVertices;
		}

		// Push vertices to limit surface
		std::unique_ptr<Point3f[]> pLimit(new Point3f[v.size()]);	// initialize an array of limit surface positions
		for (size_t i = 0; i < v.size(); ++i) {
			if (v[i]->boundary)
				pLimit[i] = weightBoundary(v[i], 1.f / 5.f);		// limit rule for boundary vertex weights
			else
				pLimit[i] = weightOneRing(v[i], loopGamma(v[i]->valence()));	// limit rule for interior vertices
		}
		for (size_t i = 0; i < v.size(); ++i) v[i]->p = pLimit[i];

		// Compute vertex tangents on limit surface
		std::vector<Normal3f> Ns;
		Ns.reserve(v.size());
		std::vector<Point3f> pRing(16, Point3f());
		for (SDVertex* vertex : v) {
			Vector3f S(0, 0, 0), T(0, 0, 0);
			int valence = vertex->valence();
			if (valence > (int)pRing.size()) pRing.resize(valence);
			vertex->oneRing(&pRing[0]);
			// compute a pair of non-parallel tangent vectors
			if (!vertex->boundary) {
				// Compute tangents of interior face
				for (int j = 0; j < valence; ++j) {
					S += std::cos(2 * Pi * j / valence) * Vector3f(pRing[j]);
					T += std::sin(2 * Pi * j / valence) * Vector3f(pRing[j]);
				}
			}
			else {
				// Compute tangents of boundary face
				S = pRing[valence - 1] - pRing[0];	// across tangent is given by the vector between the two neighboring boundary vertices
				// computing transverse tangent, based on the vertex valence
				// the weights sum to zero for all values of i, the weighted sum is in fact a tangent vector
				if (valence == 2)
					T = Vector3f(pRing[0] + pRing[1] - 2 * vertex->p);
				else if (valence == 3)
					T = pRing[1] - vertex->p;
				else if (valence == 4)  // regular
					T = Vector3f(-1 * pRing[0] + 2 * pRing[1] + 2 * pRing[2] +
						-1 * pRing[3] + -2 * vertex->p);
				else {
					Float theta = Pi / float(valence - 1);
					T = Vector3f(std::sin(theta) * (pRing[0] + pRing[valence - 1]));
					for (int k = 1; k < valence - 1; ++k) {
						Float wt = (2 * std::cos(theta) - 2) * std::sin((k)* theta);
						T += Vector3f(wt * pRing[k]);
					}
					T = -T;
				}
			}
			Ns.push_back(Normal3f(Cross(S, T)));
		}

		// Create triangle mesh from subdivision mesh
		{
			// initializes a vector of Triangles corresponding to the triangulation of the limit surface
			// transformation of the subdivided mesh into an indexed triangle mesh
			size_t ntris = f.size();
			std::unique_ptr<int[]> verts(new int[3 * ntris]);
			int* vp = verts.get();
			size_t totVerts = v.size();
			std::map<SDVertex*, int> usedVerts;
			for (size_t i = 0; i < totVerts; ++i) usedVerts[v[i]] = i;
			for (size_t i = 0; i < ntris; ++i) {
				for (int j = 0; j < 3; ++j) {
					*vp = usedVerts[f[i]->v[j]];
					++vp;
				}
			}
			return CreateTriangleMesh(ObjectToWorld, WorldToObject,
				reverseOrientation, ntris, verts.get(),
				totVerts, pLimit.get(), nullptr, &Ns[0],
				nullptr, nullptr, nullptr);
		}
	}

	std::vector<std::shared_ptr<Shape>> CreateLoopSubdiv(const Transform * o2w,
		const Transform * w2o,
		bool reverseOrientation,
		const ParamSet & params) {
		int nLevels = params.FindOneInt("levels",
			params.FindOneInt("nlevels", 3));
		int nps, nIndices;
		const int* vertexIndices = params.FindInt("indices", &nIndices);
		const Point3f* P = params.FindPoint3f("P", &nps);
		if (!vertexIndices) {
			Error("Vertex indices \"indices\" not provided for LoopSubdiv shape.");
			return std::vector<std::shared_ptr<Shape>>();
		}
		if (!P) {
			Error("Vertex positions \"P\" not provided for LoopSubdiv shape.");
			return std::vector<std::shared_ptr<Shape>>();
		}

		// don't actually use this for now...
		std::string scheme = params.FindOneString("scheme", "loop");
		return LoopSubdivide(o2w, w2o, reverseOrientation, nLevels, nIndices,
			vertexIndices, nps, P);
	}

	/// <summary>
	/// Loops overt the one-ring of adjacent vertices and applies the given weight to compute a new vertex Position.
	/// </summary>
	/// <param name="vert">The vert.</param>
	/// <param name="beta">The beta.</param>
	/// <returns></returns>
	static Point3f weightOneRing(SDVertex* vert, Float beta) {
		// Put _vert_ one-ring in _pRing_
		int valence = vert->valence();
		Point3f* pRing = ALLOCA(Point3f, valence);
		vert->oneRing(pRing);
		Point3f p = (1 - valence * beta) * vert->p;
		for (int i = 0; i < valence; ++i) p += beta * pRing[i];
		return p;
	}

	/// <summary>
	/// Ones the ring.
	/// </summary>
	/// <param name="p">Pointer, which points to an area of memory large enough to hold the one-ring-around the vertex</param>
	void SDVertex::oneRing(Point3f* p) {
		if (!boundary) {
			// Get one-ring vertices for interior vertex
			SDFace* face = startFace;
			do {
				*p++ = face->nextVert(this)->p;
				face = face->nextFace(this);
			} while (face != startFace);
		}
		else {
			// Get one-ring vertices for boundary vertex
			SDFace* face = startFace, * f2;
			while ((f2 = face->nextFace(this)) != nullptr) face = f2;	// first looping around neighbors faces until boundary (to start with)
			*p++ = face->nextVert(this)->p;
			do {
				*p++ = face->prevVert(this)->p;
				face = face->prevFace(this);
			} while (face != nullptr);
		}
	}

	/// <summary>
	/// Applies the given widhts at a boundary vertex.
	/// </summary>
	/// <param name="vert">Vertex Pointer</param>
	/// <param name="beta">The beta value</param>
	/// <returns></returns>
	static Point3f weightBoundary(SDVertex* vert, Float beta) {
		// Put _vert_ one-ring in _pRing_
		int valence = vert->valence();
		Point3f* pRing = ALLOCA(Point3f, valence);	//efficiently allocate space to store their positions
		vert->oneRing(pRing);
		Point3f p = (1 - 2 * beta) * vert->p;
		p += beta * pRing[0];				// can be used since the oneRing function orders the boundary vertex's one-ring such 
		p += beta * pRing[valence - 1];		// that the 1. and last entries are boundary neighbors.
		return p;
	}

} // namespace pbrt