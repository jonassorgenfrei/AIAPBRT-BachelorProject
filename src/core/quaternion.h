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

#ifndef PBRT_CORE_QUATERNION_H
#define PBRT_CORE_QUATERNION_H

 // core/quaternion.h*
#include "pbrt.h"
#include "stringprint.h"
#include "geometry.h"

namespace pbrt {
	// Quaternion Declarations
	// four-tuple
	struct Quaternion {
		// Quaternion Public Methods
		// Default constructor
		Quaternion() : v(0, 0, 0), w(1) {}

		// component-wise addition of quaternions
		Quaternion &operator+=(const Quaternion &q) {
			v += q.v;
			w += q.w;
			return *this;
		}

		// component-wise addition of quaternions
		friend Quaternion operator+(const Quaternion &q1, const Quaternion &q2) {
			Quaternion ret = q1;
			return ret += q2;
		}

		// component-wise subtraction of quaternions
		Quaternion &operator-=(const Quaternion &q) {
			v -= q.v;
			w -= q.w;
			return *this;
		}

		// component-wise subtraction of quaternions
		friend Quaternion operator-(const Quaternion &q1, const Quaternion &q2) {
			Quaternion ret = q1;
			return ret -= q2;
		}

		// component-wise inverse of quaternions
		Quaternion operator-() const {
			Quaternion ret;
			ret.v = -v;
			ret.w = -w;
			return ret;
		}

		// component-wise multiplication of quaternions
		Quaternion &operator*=(Float f) {
			v *= f;
			w *= f;
			return *this;
		}

		// component-wise multiplication of quaternions
		Quaternion operator*(Float f) const {
			Quaternion ret = *this;
			ret.v *= f;
			ret.w *= f;
			return ret;
		}

		// component-wise division of quaternions
		Quaternion &operator/=(Float f) {
			v /= f;
			w /= f;
			return *this;
		}

		// component-wise division of quaternions
		Quaternion operator/(Float f) const {
			Quaternion ret = *this;
			ret.v /= f;
			ret.w /= f;
			return ret;
		}

		
		/// <summary>
		/// Converts the current quaternion to a transform.
		/// </summary>
		/// <returns>Transformation</returns>
		Transform ToTransform() const;

		/// <summary>
		/// Initializes a new instance of the <see cref="Quaternion"/> struct from a <see cref="Transformation"/>.
		/// </summary>
		/// <param name="t">Transformation</param>
		Quaternion(const Transform &t);

		friend std::ostream &operator<<(std::ostream &os, const Quaternion &q) {
			os << StringPrintf("[ %f, %f, %f, %f ]", q.v.x, q.v.y, q.v.z,
				q.w);
			return os;
		}

		// Quaternion Public Data
		Vector3f v;	// imaginary 3-vector
		Float w;	// real part
	};

	/// <summary>
	/// Slerps the specified t.
	/// </summary>
	/// <param name="t">The t.</param>
	/// <param name="q1">The q1.</param>
	/// <param name="q2">The q2.</param>
	/// <returns></returns>
	Quaternion Slerp(Float t, const Quaternion &q1, const Quaternion &q2);

	// Quaternion Inline Functions
	inline Quaternion operator*(Float f, const Quaternion &q) { return q * f; }

	/// <summary>
	/// Returns the dot product of two quaternions.
	/// </summary>
	/// <param name="q1">Quaternion 1</param>
	/// <param name="q2">Quaternion 2</param>
	/// <returns>Dot product of two given quaternions.</returns>
	inline Float Dot(const Quaternion &q1, const Quaternion &q2) {
		return Dot(q1.v, q2.v) + q1.w * q2.w;
	}

	/// <summary>
	/// Normalizes the specified Quaternion q
	/// by dividing by its length.
	/// </summary>
	/// <param name="q">Quaternion</param>
	/// <returns>Normalized Quaternion</returns>
	inline Quaternion Normalize(const Quaternion &q) {
		return q / std::sqrt(Dot(q, q));
	}

}  // namespace pbrt

#endif // PBRT_CORE_QUATERNION_H