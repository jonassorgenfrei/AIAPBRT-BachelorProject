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

#ifndef PBRT_CORE_TRANSFORM_H
#define PBRT_CORE_TRANSFORM_H

 // core/transform.h*
#include "pbrt.h"
#include "stringprint.h"
#include "geometry.h"
#include "quaternion.h"

namespace pbrt {

	// Matrix4x4 Declarations
	// ----------------------
	// low-level representation of 4x4 matrices
	struct Matrix4x4 {
		// Matrix4x4 Public Methods
		// Creates identity matrix
		Matrix4x4() {
			m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.f;
			m[0][1] = m[0][2] = m[0][3] = 
			m[1][0] = m[1][2] = m[1][3] = 
			m[2][0] = m[2][1] = m[2][3] = 
			m[3][0] = m[3][1] = m[3][2] = 0.f;
		}
		Matrix4x4(Float mat[4][4]);
		Matrix4x4(Float t00, Float t01, Float t02, Float t03, 
				  Float t10, Float t11,	Float t12, Float t13, 
				  Float t20, Float t21, Float t22, Float t23,
				  Float t30, Float t31, Float t32, Float t33);
		
		// compares the given Matrix with the current Matrix; returns true if they are equal
		bool operator==(const Matrix4x4 &m2) const {
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					if (m[i][j] != m2.m[i][j]) return false;
			return true;
		}
		
		// compares the given Matrix with the current Matrix; returns true if they are not equal
		bool operator!=(const Matrix4x4 &m2) const {
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					if (m[i][j] != m2.m[i][j]) return true;
			return false;
		}

		// Transposes the current matrix and returns a new matrix
		friend Matrix4x4 Transpose(const Matrix4x4 &);

		void Print(FILE *f) const {
			fprintf(f, "[ ");
			for (int i = 0; i < 4; ++i) {
				fprintf(f, "  [ ");
				for (int j = 0; j < 4; ++j) {
					fprintf(f, "%f", m[i][j]);
					if (j != 3) fprintf(f, ", ");
				}
				fprintf(f, " ]\n");
			}
			fprintf(f, " ] ");
		}

		// Multiplies the two given matrices 
		static Matrix4x4 Mul(const Matrix4x4 &m1, const Matrix4x4 &m2) {
			Matrix4x4 r;
			// computetd by setting the (i,j)th element of the result to the inner product of the ith row 
			// of M1 with the jth element of M2
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					r.m[i][j] = m1.m[i][0] * m2.m[0][j] + m1.m[i][1] * m2.m[1][j] +
					m1.m[i][2] * m2.m[2][j] + m1.m[i][3] * m2.m[3][j];
			return r;
		}

		// inverts the current Matrix
		friend Matrix4x4 Inverse(const Matrix4x4 &);

		friend std::ostream &operator<<(std::ostream &os, const Matrix4x4 &m) {
			// clang-format off
			os << StringPrintf("[ [ %f, %f, %f, %f ] "
				"[ %f, %f, %f, %f ] "
				"[ %f, %f, %f, %f ] "
				"[ %f, %f, %f, %f ] ]",
				m.m[0][0], m.m[0][1], m.m[0][2], m.m[0][3],
				m.m[1][0], m.m[1][1], m.m[1][2], m.m[1][3],
				m.m[2][0], m.m[2][1], m.m[2][2], m.m[2][3],
				m.m[3][0], m.m[3][1], m.m[3][2], m.m[3][3]);
			// clang-format on
			return os;
		}

		// Matrix4x4 Public Data
		Float m[4][4];	// stored in row-major form
	};

	// Transform Declarations
	// ----------------------
	class Transform {
		public:
			// Transform Public Methods

			// Basic Operations
			// ----------------
			Transform() {}
			Transform(const Float mat[4][4]) {
				m = Matrix4x4(mat[0][0], mat[0][1], mat[0][2], mat[0][3], mat[1][0],
					mat[1][1], mat[1][2], mat[1][3], mat[2][0], mat[2][1],
					mat[2][2], mat[2][3], mat[3][0], mat[3][1], mat[3][2],
					mat[3][3]);
				mInv = Inverse(m);
			}


			Transform(const Matrix4x4 &m) : m(m), mInv(Inverse(m)) {}

			// taking matrix and (pre-computed) inverse to avoid the
			// expense and potential loss of precision from computing a general
			// 4x4 matrix inverse
			Transform(const Matrix4x4 &m, 
						const Matrix4x4 &mInv) : m(m), mInv(mInv) {}

			void Print(FILE *f) const;
			
			// inverse Transformation (swapping m and mInv)
			friend Transform Inverse(const Transform &t) {
				return Transform(t.mInv, t.m);
			}

			// transpose the two matrices in the the transform (to create a new transform)
			friend Transform Transpose(const Transform &t) {
				return Transform(Transpose(t.m), Transpose(t.mInv));
			}

			// compares the current matrix with the given matrix (and their inverse matrices), 
			// returns true if they are equal
			bool operator==(const Transform &t) const {
				return t.m == m && t.mInv == mInv;
			}

			// compares the current matrix with the given matrix (and their inverse matrices), 
			// returns true if they are not equal
			bool operator!=(const Transform &t) const {
				return t.m != m || t.mInv != mInv;
			}

			bool operator<(const Transform &t2) const {
				for (int i = 0; i < 4; ++i)
					for (int j = 0; j < 4; ++j) {
						if (m.m[i][j] < t2.m.m[i][j]) return true;
						if (m.m[i][j] > t2.m.m[i][j]) return false;
					}
				return false;
			}

			// checks if the given matrix is an Identity Matrix
			bool IsIdentity() const {
				return (m.m[0][0] == 1.f && m.m[0][1] == 0.f && m.m[0][2] == 0.f &&
					m.m[0][3] == 0.f && m.m[1][0] == 0.f && m.m[1][1] == 1.f &&
					m.m[1][2] == 0.f && m.m[1][3] == 0.f && m.m[2][0] == 0.f &&
					m.m[2][1] == 0.f && m.m[2][2] == 1.f && m.m[2][3] == 0.f &&
					m.m[3][0] == 0.f && m.m[3][1] == 0.f && m.m[3][2] == 0.f &&
					m.m[3][3] == 1.f);
			}

			// returns the current matrix
			const Matrix4x4 &GetMatrix() const { return m; }
			// returns the inverse of the current matrix
			const Matrix4x4 &GetInverseMatrix() const { return mInv; }

			/// <summary>
			/// Determines whether this instance has scale.
			/// Test if transformation has a scaling term in it
			/// by taking the three coordinate axes and see if any of their lengths
			/// are appreciably different from one
			/// </summary>
			/// <returns>
			///   <c>true</c> if this instance has scale; otherwise, <c>false</c>.
			/// </returns>
			bool HasScale() const {
				Float la2 = (*this)(Vector3f(1, 0, 0)).LengthSquared();
				Float lb2 = (*this)(Vector3f(0, 1, 0)).LengthSquared();
				Float lc2 = (*this)(Vector3f(0, 0, 1)).LengthSquared();

#define NOT_ONE(x) ((x) < .999f || (x) > 1.001f)
				return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
			};

			/// <summary>
			/// Point transformation routine
			/// Takes a point (x,y,z) and implicitly represents it as the homogeneous column vector
			/// [x y z 1]T. It then transforms the point by pre-multiplying this vector with the 
			/// transformation matrix.
			/// Finally it divides by w to convert back to a non homogeneous point representation.
			/// </summary>
			/// <param name="p">Point to transform</param>
			/// <returns></returns>
			template <typename T>
			inline Point3<T> operator()(const Point3<T> &p) const;

			/// <summary>
			/// Vector transformation routine
			/// Simplified since the implicit homogeneous w coordinate is zero.
			/// </summary>
			/// <param name="v">Vector to transform</param>
			/// <returns></returns>
			template <typename T>
			inline Vector3<T> operator()(const Vector3<T> &v) const;

			/// <summary>
			/// Normal transformation routine
			/// Multiplies the normal with the inverse transpose of the transformation matrix
			/// </summary>
			/// <param name="">The .</param>
			/// <returns></returns>
			template <typename T>
			inline Normal3<T> operator()(const Normal3<T> &) const;

			/// <summary>
			/// Ray transformation routine
			/// Transforming the origin and direction and copying the other data members
			/// </summary>
			/// <param name="r">The r.</param>
			/// <returns></returns>
			inline Ray operator()(const Ray &r) const;

			/// <summary>
			/// Ray transformation routine
			/// Transforming the origin and direction and copying the other data members
			/// </summary>
			/// <param name="r">The r.</param>
			/// <returns></returns>
			inline RayDifferential operator()(const RayDifferential &r) const;

			/// <summary>
			/// Bounding Box transformation routine
			/// Transform all eight of its corner vertices and then compute a new bounding box that encompasses
			/// those points.
			/// TODO: MAKE MORE EFFICIENT!
			/// </summary>
			/// <param name="b">The b.</param>
			/// <returns></returns>
			Bounds3f operator()(const Bounds3f &b) const;

			
			/// <summary>
			/// Composition of Transformations
			/// Multiplies the given Transformation to the current matrix.
			/// </summary>
			/// <param name="t2">Transformation</param>
			/// <returns></returns>
			Transform operator*(const Transform &t2) const;

			/// <summary>
			/// Checks if the handedness is changed by a transformation. 
			/// (happens only when the determinant of the transformation's upper-left 3x3 sub matrix is negative)
			/// </summary>
			/// <returns>
			///		<c>true</c> if the handedness is changed; otherwise, <c>false</c>.
			/// </returns>
			bool SwapsHandedness() const;

			SurfaceInteraction operator()(const SurfaceInteraction &si) const;
			template <typename T>
			inline Point3<T> operator()(const Point3<T> &pt,
				Vector3<T> *absError) const;
			template <typename T>
			inline Point3<T> operator()(const Point3<T> &p, const Vector3<T> &pError,
				Vector3<T> *pTransError) const;
			template <typename T>
			inline Vector3<T> operator()(const Vector3<T> &v,
				Vector3<T> *vTransError) const;
			template <typename T>
			inline Vector3<T> operator()(const Vector3<T> &v, const Vector3<T> &vError,
				Vector3<T> *vTransError) const;
			inline Ray operator()(const Ray &r, Vector3f *oError,
				Vector3f *dError) const;
			inline Ray operator()(const Ray &r, const Vector3f &oErrorIn,
				const Vector3f &dErrorIn, Vector3f *oErrorOut,
				Vector3f *dErrorOut) const;

			friend std::ostream &operator<<(std::ostream &os, const Transform &t) {
				os << "t=" << t.m << ", inv=" << t.mInv;
				return os;
			}
		private:
			// Transform Private Data
			Matrix4x4 m,
					mInv;	// inverse	
			friend class AnimatedTransform;
			friend struct Quaternion;
	};

	// Translate 
	// ---------
	Transform Translate(const Vector3f &delta);

	// Scaling
	// -------
	Transform Scale(Float x, Float y, Float z);

	// x,y,z-Rotations
	// ---------------
	Transform RotateX(Float theta);
	Transform RotateY(Float theta);
	Transform RotateZ(Float theta);

	// Rotation around arbitrary axis
	// ------------------------------
	Transform Rotate(Float theta, const Vector3f &axis);


	// Look-At Transformation
	// ----------------------
	Transform LookAt(const Point3f &pos, const Point3f &look, const Vector3f &up);
	

	// Transform Inline Functions	

	template <typename T>
	inline Point3<T> Transform::operator()(const Point3<T> &p) const {
		T x = p.x, y = p.y, z = p.z;
		T xp = m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z + m.m[0][3];
		T yp = m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z + m.m[1][3];
		T zp = m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z + m.m[2][3];
		T wp = m.m[3][0] * x + m.m[3][1] * y + m.m[3][2] * z + m.m[3][3];
		CHECK_NE(wp, 0);
		if (wp == 1)	// skips division if w is 1
			return Point3<T>(xp, yp, zp);
		else
			return Point3<T>(xp, yp, zp) / wp;
	}

	template <typename T>
	inline Vector3<T> Transform::operator()(const Vector3<T> &v) const {
		T x = v.x, y = v.y, z = v.z;
		return Vector3<T>(	m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
							m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
							m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
	}	

	template <typename T>
	inline Normal3<T> Transform::operator()(const Normal3<T> &n) const {
		T x = n.x, y = n.y, z = n.z;
		return Normal3<T>(mInv.m[0][0] * x + mInv.m[1][0] * y + mInv.m[2][0] * z,
							mInv.m[0][1] * x + mInv.m[1][1] * y + mInv.m[2][1] * z,
							mInv.m[0][2] * x + mInv.m[1][2] * y + mInv.m[2][2] * z);
	}	

	inline Ray Transform::operator()(const Ray &r) const {
		Vector3f oError;
		Point3f o = (*this)(r.o, &oError);
		Vector3f d = (*this)(r.d);

		// Offset ray origin to edge of error bounds and compute _tMax_
		Float lengthSquared = d.LengthSquared();
		Float tMax = r.tMax;
		if (lengthSquared > 0) {	// manage floating-point round-off errors
			Float dt = Dot(Abs(d), oError) / lengthSquared;
			o += d * dt;
			tMax -= dt;
		}

		return Ray(o, d, tMax, r.time, r.medium);
	}
	inline RayDifferential Transform::operator()(const RayDifferential &r) const {
		Ray tr = (*this)(Ray(r));
		RayDifferential ret(tr.o, tr.d, tr.tMax, tr.time, tr.medium);
		ret.hasDifferentials = r.hasDifferentials;
		ret.rxOrigin = (*this)(r.rxOrigin);
		ret.ryOrigin = (*this)(r.ryOrigin);
		ret.rxDirection = (*this)(r.rxDirection);
		ret.ryDirection = (*this)(r.ryDirection);
		return ret;
	}
	
	// AnimatedTransform Declarations
	class AnimatedTransform {
	public:
		// AnimatedTransform Public Methods
		AnimatedTransform(const Transform *startTransform, Float startTime,
			const Transform *endTransform, Float endTime);
		static void Decompose(const Matrix4x4 &m, Vector3f *T, Quaternion *R,
			Matrix4x4 *S);
		void Interpolate(Float time, Transform *t) const;
		Ray operator()(const Ray &r) const;
		RayDifferential operator()(const RayDifferential &r) const;
		Point3f operator()(Float time, const Point3f &p) const;
		Vector3f operator()(Float time, const Vector3f &v) const;
		bool HasScale() const {
			return startTransform->HasScale() || endTransform->HasScale();
		}
		Bounds3f MotionBounds(const Bounds3f &b) const;
		Bounds3f BoundPointMotion(const Point3f &p) const;

	private:
		// AnimatedTransform Private Data
		const Transform *startTransform, *endTransform;
		const Float startTime, endTime;
		const bool actuallyAnimated;
		Vector3f T[2];
		Quaternion R[2];
		Matrix4x4 S[2];
		bool hasRotation;
		struct DerivativeTerm {
			DerivativeTerm() {}
			DerivativeTerm(Float c, Float x, Float y, Float z)
				: kc(c), kx(x), ky(y), kz(z) {}
			Float kc, kx, ky, kz;
			Float Eval(const Point3f &p) const {
				return kc + kx * p.x + ky * p.y + kz * p.z;
			}
		};
		DerivativeTerm c1[3], c2[3], c3[3], c4[3], c5[3];
	};


}  // namespace pbrt

#endif // PBRT_CORE_TRANSFORM_H