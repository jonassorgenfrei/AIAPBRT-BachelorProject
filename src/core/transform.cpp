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

 // core/transform.cpp*
#include "transform.h"
//#include "interaction.h"

namespace pbrt {
	// Matrix4x4 Method Definitions
	bool SolveLinearSystem2x2(const Float A[2][2], 
								const Float B[2], 
								Float *x0,
								Float *x1) {
		Float det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
		if (std::abs(det) < 1e-10f) return false;
		*x0 = (A[1][1] * B[0] - A[0][1] * B[1]) / det;
		*x1 = (A[0][0] * B[1] - A[1][0] * B[0]) / det;
		if (std::isnan(*x0) || std::isnan(*x1)) return false;
		return true;
	}

	
	Matrix4x4::Matrix4x4(Float mat[4][4]) { memcpy(m, mat, 16 * sizeof(Float)); }

	Matrix4x4::Matrix4x4(Float t00, Float t01, Float t02, Float t03, Float t10,
		Float t11, Float t12, Float t13, Float t20, Float t21,
		Float t22, Float t23, Float t30, Float t31, Float t32,
		Float t33) {
		m[0][0] = t00;
		m[0][1] = t01;
		m[0][2] = t02;
		m[0][3] = t03;
		m[1][0] = t10;
		m[1][1] = t11;
		m[1][2] = t12;
		m[1][3] = t13;
		m[2][0] = t20;
		m[2][1] = t21;
		m[2][2] = t22;
		m[2][3] = t23;
		m[3][0] = t30;
		m[3][1] = t31;
		m[3][2] = t32;
		m[3][3] = t33;
	}

	Matrix4x4 Transpose(const Matrix4x4 &m) {
		return Matrix4x4(m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0], m.m[0][1],
			m.m[1][1], m.m[2][1], m.m[3][1], m.m[0][2], m.m[1][2],
			m.m[2][2], m.m[3][2], m.m[0][3], m.m[1][3], m.m[2][3],
			m.m[3][3]);
	}

	Matrix4x4 Inverse(const Matrix4x4 &m) {
		// numerically stable Gauss-Jordan elimination routine
		int indxc[4], indxr[4];
		int ipiv[4] = { 0, 0, 0, 0 };
		Float minv[4][4];
		memcpy(minv, m.m, 4 * 4 * sizeof(Float));
		for (int i = 0; i < 4; i++) {
			int irow = 0, icol = 0;
			Float big = 0.f;
			// Choose pivot
			for (int j = 0; j < 4; j++) {
				if (ipiv[j] != 1) {
					for (int k = 0; k < 4; k++) {
						if (ipiv[k] == 0) {
							if (std::abs(minv[j][k]) >= big) {
								big = Float(std::abs(minv[j][k]));
								irow = j;
								icol = k;
							}
						}
						else if (ipiv[k] > 1)
							 Error("Singular matrix in MatrixInvert");
					}
				}
			}
			++ipiv[icol];
			// Swap rows _irow_ and _icol_ for pivot
			if (irow != icol) {
				for (int k = 0; k < 4; ++k) std::swap(minv[irow][k], minv[icol][k]);
			}
			indxr[i] = irow;
			indxc[i] = icol;
			if (minv[icol][icol] == 0.f) Error("Singular matrix in MatrixInvert");

			// Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
			Float pivinv = 1. / minv[icol][icol];
			minv[icol][icol] = 1.;
			for (int j = 0; j < 4; j++) minv[icol][j] *= pivinv;

			// Subtract this row from others to zero out their columns
			for (int j = 0; j < 4; j++) {
				if (j != icol) {
					Float save = minv[j][icol];
					minv[j][icol] = 0;
					for (int k = 0; k < 4; k++) minv[j][k] -= minv[icol][k] * save;
				}
			}
		}
		// Swap columns to reflect permutation
		for (int j = 3; j >= 0; j--) {
			if (indxr[j] != indxc[j]) {
				for (int k = 0; k < 4; k++)
					std::swap(minv[k][indxr[j]], minv[k][indxc[j]]);
			}
		}
		return Matrix4x4(minv);
	}

	// Transform Method Definitions
	void Transform::Print(FILE *f) const { m.Print(f); }

	// Translation
	// -----------
	Transform Translate(const Vector3f &delta) {
		Matrix4x4 m(1, 0, 0, delta.x, 
					0, 1, 0, delta.y, 
					0, 0, 1, delta.z, 
					0, 0, 0,		1);
		Matrix4x4 minv( 1, 0, 0, -delta.x, 
						0, 1, 0, -delta.y, 
						0, 0, 1, -delta.z, 
						0, 0, 0, 1);
		return Transform(m, minv);
	}

	// Scaling
	// -------
	Transform Scale(Float x, Float y, Float z) {
		Matrix4x4 m(x, 0, 0, 0, 
					0, y, 0, 0,
					0, 0, z, 0,
					0, 0, 0, 1);
		Matrix4x4 minv(	1 / x,	0,		0,		0,
						0,		1 / y,	0,		0,
						0,		0,		1 / z,	0,
						0,		0,		0,		1);
		return Transform(m, minv);
	}

	// x,y,z-Rotations
	// ---------------
	Transform RotateX(Float theta) {
		Float sinTheta = std::sin(Radians(theta));
		Float cosTheta = std::cos(Radians(theta));
		Matrix4x4 m(1,			0,			0,	0, 
					0,	cosTheta,	-sinTheta,	0,
					0, sinTheta,	cosTheta,	0,
					0,			0,	 0,			1);
		// inverse equals the transpose for orthogonal matrices!
		return Transform(m, Transpose(m));
	}

	Transform RotateY(Float theta) {
		Float sinTheta = std::sin(Radians(theta));
		Float cosTheta = std::cos(Radians(theta));
		Matrix4x4 m(cosTheta,	0,	sinTheta,	0,
							0,	1,			0,	0,
					-sinTheta,	0,	cosTheta,	0,
							0,	0,			0,	1);
		// inverse equals the transpose for orthogonal matrices!
		return Transform(m, Transpose(m));
	}

	Transform RotateZ(Float theta) {
		Float sinTheta = std::sin(Radians(theta));
		Float cosTheta = std::cos(Radians(theta));
		Matrix4x4 m(cosTheta,	-sinTheta,	0,	0,
					sinTheta,	cosTheta,	0,	0,
							0,			0,	1,	0,
							0,			0,	0,	1);
		// inverse equals the transpose for orthogonal matrices!
		return Transform(m, Transpose(m));
	}

	// Rotation around arbitrary axis
	// ------------------------------
	Transform Rotate(Float theta, const Vector3f &axis) {
		// normalize(axis)
		// Vector vc through the end point of v and parallel to as
		// alpha -> angle between axis and v
		// vc = a ||v|| cos(alpha) = a(v dot a)
		// pair of basis vectors v1 and v2
		// v1 = v - vc
		// v2 = (v1 x a)
		// v' = vc + v1 * cos(theta) + v2 * sin(theta)

		Vector3f a = Normalize(axis);
		Float sinTheta = std::sin(Radians(theta));
		Float cosTheta = std::cos(Radians(theta));
		Matrix4x4 m;
		// Compute rotation of first basis vector
		m.m[0][0] = a.x * a.x + (1 - a.x * a.x) * cosTheta;
		m.m[0][1] = a.x * a.y * (1 - cosTheta) - a.z * sinTheta;
		m.m[0][2] = a.x * a.z * (1 - cosTheta) + a.y * sinTheta;
		m.m[0][3] = 0;

		// Compute rotations of second and third basis vectors
		m.m[1][0] = a.x * a.y * (1 - cosTheta) + a.z * sinTheta;
		m.m[1][1] = a.y * a.y + (1 - a.y * a.y) * cosTheta;
		m.m[1][2] = a.y * a.z * (1 - cosTheta) - a.x * sinTheta;
		m.m[1][3] = 0;

		m.m[2][0] = a.x * a.z * (1 - cosTheta) - a.y * sinTheta;
		m.m[2][1] = a.y * a.z * (1 - cosTheta) + a.x * sinTheta;
		m.m[2][2] = a.z * a.z + (1 - a.z * a.z) * cosTheta;
		m.m[2][3] = 0;

		// as with the other rotation matrices, the inverse is equal to the transpose
		return Transform(m, Transpose(m));
	}

	// Look-At Transformation
	// ----------------------
	Transform LookAt(const Point3f &pos, const Point3f &look, const Vector3f &up) {
		Matrix4x4 cameraToWorld;
		// Initialize fourth column of viewing matrix
		cameraToWorld.m[0][3] = pos.x;
		cameraToWorld.m[1][3] = pos.y;
		cameraToWorld.m[2][3] = pos.z;
		cameraToWorld.m[3][3] = 1;

		// Initialize first three columns of viewing matrix
		Vector3f dir = Normalize(look - pos);		// compute normalized direction (3rd column)
		// left-handed coordinate system: camera space is defined with the viewing direction down the +z
		// axis
		if (Cross(Normalize(up), dir).Length() == 0) {
			Error(
				"\"up\" vector (%f, %f, %f) and viewing direction (%f, %f, %f) "
				"passed to LookAt are pointing in the same direction.  Using "
				"the identity transformation.",
				up.x, up.y, up.z, dir.x, dir.y, dir.z);
			return Transform();
		}

		Vector3f right = Normalize(Cross(Normalize(up), dir)); // world space direction => +x axis in camera space
		Vector3f newUp = Cross(dir, right); // up-vector 
		cameraToWorld.m[0][0] = right.x;
		cameraToWorld.m[1][0] = right.y;
		cameraToWorld.m[2][0] = right.z;
		cameraToWorld.m[3][0] = 0.;
		cameraToWorld.m[0][1] = newUp.x;
		cameraToWorld.m[1][1] = newUp.y;
		cameraToWorld.m[2][1] = newUp.z;
		cameraToWorld.m[3][1] = 0.;
		cameraToWorld.m[0][2] = dir.x;
		cameraToWorld.m[1][2] = dir.y;
		cameraToWorld.m[2][2] = dir.z;
		cameraToWorld.m[3][2] = 0.;

		return Transform(Inverse(cameraToWorld), cameraToWorld);
	}

} // namespace pbrt
