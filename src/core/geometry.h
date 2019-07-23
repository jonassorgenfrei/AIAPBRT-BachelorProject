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

#ifndef PBRT_CORE_GEOMETRY_H
#define PBRT_CORE_GEOMETRY_H

 // core/geometry.h*
#include "pbrt.h"
#include "stringprint.h"
#include <iterator>

namespace pbrt {

	template <typename T>
	inline bool isNaN(const T x) {
		return std::isnan(x);
	}
	template <>
	inline bool isNaN(const int x) {
		return false;
	}

	// Vector Declarations
	// -------------------
	template <typename T> class Vector2 {
		public:
			// Vector2 Public Methods
			Vector2() { x = y = 0; }
			Vector2(T xx, T yy) : x(xx), y(yy) { DCHECK(!HasNaNs()); }
			bool HasNaNs() const { return isNaN(x) || isNaN(y); }
			explicit Vector2(const Point2<T>& p);
			explicit Vector2(const Point3<T>& p);
#ifndef NDEBUG
			// The default versions of these are fine for release builds; for debug
			// we define them so that we can add the Assert checks.
			Vector2(const Vector2<T>& v) {
				DCHECK(!v.HasNaNs());
				x = v.x;
				y = v.y;
			}
			Vector2<T>& operator=(const Vector2<T>& v) {
				DCHECK(!v.HasNaNs());
				x = v.x;
				y = v.y;
				return *this;
			}
#endif // !NDEBUG

			Vector2<T> operator+(const Vector2<T>& v) const {
				DCHECK(!v.HasNaNs());
				return Vector2(x + v.x, y + v.y);
			}

			Vector2<T>& operator+=(const Vector2<T>& v) {
				DCHECK(!v.HasNaNs());
				x += v.x;
				y += v.y;
				return *this;
			}
			Vector2<T> operator-(const Vector2<T>& v) const {
				DCHECK(!v.HasNaNs());
				return Vector2(x - v.x, y - v.y);
			}

			Vector2<T>& operator-=(const Vector2<T>& v) {
				DCHECK(!v.HasNaNs());
				x -= v.x;
				y -= v.y;
				return *this;
			}

			bool operator==(const Vector2<T>& v) const { return x == v.x && y == v.y; }
			bool operator!=(const Vector2<T>& v) const { return x != v.x || y != v.y; }
			template <typename U>
			Vector2<T> operator*(U f) const {
				return Vector2<T>(f * x, f * y);
			}

			template <typename U>
			Vector2<T>& operator*=(U f) {
				DCHECK(!isNaN(f));
				x *= f;
				y *= f;
				return *this;
			}
			template <typename U>
			Vector2<T> operator/(U f) const {
				CHECK_NE(f, 0);
				Float inv = (Float)1 / f;
				return Vector2<T>(x * inv, y * inv);
			}

			template <typename U>
			Vector2<T>& operator/=(U f) {
				CHECK_NE(f, 0);
				Float inv = (Float)1 / f;
				x *= inv;
				y *= inv;
				return *this;
			}
			Vector2<T> operator-() const { return Vector2<T>(-x, -y); }
			T operator[](int i) const {
				DCHECK(i >= 0 && i <= 1);
				if (i == 0) return x;
				return y;
			}

			T& operator[](int i) {
				DCHECK(i >= 0 && i <= 1);
				if (i == 0) return x;
				return y;
			}
			Float LengthSquared() const { return x * x + y * y; }
			Float Length() const { return std::sqrt(LengthSquared()); }

			// Vector2 Public Data
			T x, y;
	};

	template <typename T>
	inline std::ostream& operator<<(std::ostream& os, const Vector2<T>& v) {
		os << "[ " << v.x << ", " << v.y << " ]";
		return os;
	}

	template <>
	inline std::ostream& operator<<(std::ostream& os, const Vector2<Float>& v) {
		os << StringPrintf("[ %f, %f ]", v.x, v.y);
		return os;
	}

	template <typename T> class Vector3 {
	public:
		// Vector3 Public Methods
		T operator[](int i) const {
			DCHECK(i >= 0 && i <= 2);
			if (i == 0) return x;
			if (i == 1) return y;
			return z;
		}

		T& operator[](int i) {
			DCHECK(i >= 0 && i <= 2);
			if (i == 0) return x;
			if (i == 1) return y;
			return z;
		}

		// Default Constructor 
		Vector3() { x = y = z = 0; }
		Vector3(T x, T y, T z) : x(x), y(y), z(z) { DCHECK(!HasNaNs()); }
		bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z);	}
		explicit Vector3(const Point3<T> &p);
#ifndef NDEBUG
		// The default versions of these are fine for release builds; for debug
		// we define them so that we can add the Assert checks.
		Vector3(const Vector3<T> &v) {
			DCHECK(!v.HasNaNs());
			x = v.x;
			y = v.y;
			z = v.z;
		}

		Vector3<T> &operator=(const Vector3<T> &v) {
			DCHECK(!v.HasNaNs());
			x = v.x;
			y = v.y;
			z = v.z;
			return *this;
		}
#endif // !NDEBUG

		Vector3<T> operator+(const Vector3<T> &v) const {
			DCHECK(!HasNaNs());
			return Vector3(x + v.x, y + v.y, z + v.z);
		}

		Vector3<T> operator+=(const Vector3<T> &v) {
			DCHECK(!v.HasNaNs());
			x += v.x;
			y += v.y;
			z += v.z;
			return *this;
		}

		Vector3<T> operator-(const Vector3<T> &v) const {
			DCHECK(!v.HasNaNs());
			return Vector3(x - v.x, y - v.y, z - v.z);
		}

		Vector3<T>& operator-=(const Vector3<T> &v) {
			x -= v.x;
			y -= v.y;
			z -= v.z;
			return *this;
		}

		bool operator==(const Vector3<T> &v) const {
			return x == v.x && y == v.y && z == v.z;
		}

		bool operator!=(const Vector3<T> &v) const {
			return x != v.x || y != v.y || z != v.z;
		}

		// Component-wise scalar multiplication
		template <typename U>
		Vector3<T> operator*(U s) const { 
			return Vector3<T>(s*x, s*y, s*z);
		}

		template <typename U>
		Vector3<T> &operator*=(U s) {
			DCHECK(!isNaN(s));
			x *= s;
			y *= s;
			z *= s;
			return *this;
		}

		// Component-wise scalar division
		template <typename U>
		Vector3<T> operator/(U f) const {
			CHECK_NE(f, 0);
			Float inv = (Float)1 / f;
			return Vector3<T>(x*inv, y*inv, z*inv);
		}

		Vector3<T> &operator/=(T f) {
			CHECK_NE(f, 0);
			Float inv = (Float)1 / f;
			x *= inv;
			y *= inv;
			z *= inv;
			return *this;
		}

		// unary negation operator (returns a bew vector pointing in the opposite direction of the original one)
		Vector3<T> operator-() const { return Vector3<T>(-x, -y, -z); }

		Float LengthSquared() const { return x * x + y * y + z*z; }
		Float Length() const { return std::sqrt(LengthSquared()); }
		explicit Vector3(const Normal3<T> &n);

		// Vector3 Public Data
		T x, y, z;
	};

	template <typename T>
	inline std::ostream& operator<<(std::ostream& os, const Vector3<T> &v) {
		os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
		return os;
	}

	template <>
	inline std::ostream &operator<<(std::ostream &os, const Vector3<Float> &v) {
		os << StringPrintf("[ %f, %f, %f ]", v.x, v.y, v.z);
		return os;
	}

	typedef Vector2<Float> Vector2f;
	typedef Vector2<int> Vector2i;
	typedef Vector3<Float> Vector3f;
	typedef Vector3<int> Vector3i;


	// Point Declarations
	// ------------------
	template <typename T>
	class Point2 {
		public:
			// Point 2 Public Methods
			
			/// <summary>
			/// Initializes a new instance of the <see cref="Point2"/> class from a <see cref="Point3"/>
			/// by dropping the z coordinate.
			/// Explicit qualifier so that this conversion can't happen without an explicit cast, lest
			/// it happen unintentionally
			/// </summary>
			/// <param name="p">Point3</param>
			explicit Point2(const Point3<T> &p) : x(p.x), y(p.y) { DCHECK(!HasNaNs()); }
			// Default constructor
			Point2() { x = y = 0; }
			Point2(T xx, T yy) : x(xx), y(yy) { DCHECK(!HasNaNs()); }

			template <typename U>
			explicit Point2(const Point2<U> &p) {
				x = (T)p.x;
				y = (T)p.y;
				DCHECK(!HasNaNs());
			}

			template <typename U>
			explicit Point2(const Vector2<U> &p) {
				x = (T)p.x;
				y = (T)p.y;
				DCHECK(!HasNaNs());
			}

			template <typename U>
			explicit operator Vector2<U>() const {
				return Vector2<U>(x, y);
			}

#ifndef NDEBUG
			Point2(const Point2<T> &p) {
				DCHECK(!p.HasNaNs());
				x = p.x;
				y = p.y;
			}

			Point2<T> &operator=(const Point2<T> &p) {
				DCHECK(!p.HasNaNs());
				x = p.x;
				y = p.y;
				return *this;
			}
#endif // !NDEBUG
			Point2<T> operator+(const Vector2<T> &v) const {
				DCHECK(!v.HasNaNs());
				return Point2<T>(x + v.x, y + v.y);
			}

			Point2<T> &operator+=(const Vector2<T> &v) {
				DCHECK(!v.HasNaNs());
				x += v.x;
				y += v.y;
				return *this;
			}
			Vector2<T> operator-(const Point2<T> &p) const {
				DCHECK(!p.HasNaNs());
				return Vector2<T>(x - p.x, y - p.y);
			}

			Point2<T> operator-(const Vector2<T> &v) const {
				DCHECK(!v.HasNaNs());
				return Point2<T>(x - v.x, y - v.y);
			}
			Point2<T> operator-() const { return Point2<T>(-x, -y); }
			Point2<T> &operator-=(const Vector2<T> &v) {
				DCHECK(!v.HasNaNs());
				x -= v.x;
				y -= v.y;
				return *this;
			}
			Point2<T> &operator+=(const Point2<T> &p) {
				DCHECK(!p.HasNaNs());
				x += p.x;
				y += p.y;
				return *this;
			}
			Point2<T> operator+(const Point2<T> &p) const {
				DCHECK(!p.HasNaNs());
				return Point2<T>(x + p.x, y + p.y);
			}
			template <typename U>
			Point2<T> operator*(U f) const {
				return Point2<T>(f * x, f * y);
			}
			template <typename U>
			Point2<T> &operator*=(U f) {
				x *= f;
				y *= f;
				return *this;
			}
			template <typename U>
			Point2<T> operator/(U f) const {
				CHECK_NE(f, 0);
				Float inv = (Float)1 / f;
				return Point2<T>(inv * x, inv * y);
			}
			template <typename U>
			Point2<T> &operator/=(U f) {
				CHECK_NE(f, 0);
				Float inv = (Float)1 / f;
				x *= inv;
				y *= inv;
				return *this;
			}
			T operator[](int i) const {
				DCHECK(i >= 0 && i <= 1);
				if (i == 0) return x;
				return y;
			}

			T &operator[](int i) {
				DCHECK(i >= 0 && i <= 1);
				if (i == 0) return x;
				return y;
			}
			bool operator==(const Point2<T> &p) const { return x == p.x && y == p.y; }
			bool operator!=(const Point2<T> &p) const { return x != p.x || y != p.y; }
			bool HasNaNs() const { return isNaN(x) || isNaN(y); }

			// Point2 Public Data
			T x, y;
	};

	template <typename T>
	inline std::ostream &operator<<(std::ostream &os, const Point2<T> &v) {
		os << "[ " << v.x << ", " << v.y << " ]";
		return os;
	}

	template <>
	inline std::ostream &operator<<(std::ostream &os, const Point2<Float> &v) {
		os << StringPrintf("[ %f, %f ]", v.x, v.y);
		return os;
	}

	template <typename T>
	class Point3 {
	public:
		// Point3 Public Methods
		// Default constructor
		Point3() { x = y = z = 0; }
		Point3(T x, T y, T z) : x(x), y(y), z(z) { DCHECK(!HasNaNs()); }
		
		// Converts a point with one element type to a point of another one
		template <typename U>
		explicit Point3(const Point3<U> &p)
			: x((T)p.x), y((T)p.y), z((T)p.z) {
			DCHECK(!HasNaNs());
		}
		// Converts a point to a vector with a different underlying element type
		template <typename U>
		explicit operator Vector3<U>() const {
			return Vector3<U>(x, y, z);
		}
#ifndef NDEBUG
		Point3(const Point3<T> &p) {
			DCHECK(!p.HasNaNs());
			x = p.x;
			y = p.y;
			z = p.z;
		}

		Point3<T> &operator=(const Point3<T> &p) {
			DCHECK(!p.HasNaNs());
			x = p.x;
			y = p.y;
			z = p.z;
			return *this;
		}
#endif // !NDEBUG
		// Add a vector to a point, offsetting it in the given direction to obtain a new point.
		Point3<T> operator+(const Vector3<T> &v) const {
			DCHECK(!v.HasNaNs());
			return Point3<T>(x + v.x, y + v.y, z + v.z);
		}
		Point3<T> &operator+=(const Vector3<T> &v) {
			DCHECK(!v.HasNaNs());
			x += v.x;
			y += v.y;
			z += v.z;
			return *this;
		}
		// Subtracts one point from another obtaining the vector between them.
		Vector3<T> operator-(const Point3<T> &p) const {
			DCHECK(!p.HasNaNs());
			return Vector3<T>(x - p.x, y - p.y, z - p.z);
		}
		// Subtracting a vector from a point gives a new point
		Point3<T> operator-(const Vector3<T> &v) const {
			DCHECK(!v.HasNaNs());
			return Point3<T>(x - v.x, y - v.y, z - v.z);
		}
		Point3<T> &operator-=(const Vector3<T> &v) {
			DCHECK(!v.HasNaNs());
			x -= v.x;
			y -= v.y;
			z -= v.z;
			return *this;
		}
		Point3<T> &operator+=(const Point3<T> &p) {
			DCHECK(!p.HasNaNs());
			x += p.x;
			y += p.y;
			z += p.z;
			return *this;
		}
		Point3<T> operator+(const Point3<T> &p) const {
			DCHECK(!p.HasNaNs());
			return Point3<T>(x + p.x, y + p.y, z + p.z);
		}
		template <typename U>
		Point3<T> operator*(U f) const {
			return Point3<T>(f * x, f * y, f * z);
		}
		template <typename U>
		Point3<T> &operator*=(U f) {
			x *= f;
			y *= f;
			z *= f;
			return *this;
		}
		template <typename U>
		Point3<T> operator/(U f) const {
			CHECK_NE(f, 0);
			Float inv = (Float)1 / f;
			return Point3<T>(inv * x, inv * y, inv * z);
		}
		template <typename U>
		Point3<T> &operator/=(U f) {
			CHECK_NE(f, 0);
			Float inv = (Float)1 / f;
			x *= inv;
			y *= inv;
			z *= inv;
			return *this;
		}
		T operator[](int i) const {
			DCHECK(i >= 0 && i <= 2);
			if (i == 0) return x;
			if (i == 1) return y;
			return z;
		}

		T &operator[](int i) {
			DCHECK(i >= 0 && i <= 2);
			if (i == 0) return x;
			if (i == 1) return y;
			return z;
		}
		bool operator==(const Point3<T> &p) const {
			return x == p.x && y == p.y && z == p.z;
		}
		bool operator!=(const Point3<T> &p) const {
			return x != p.x || y != p.y || z != p.z;
		}
		bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }
		Point3<T> operator-() const { return Point3<T>(-x, -y, -z); }

		// Point3 Public Data
		T x, y, z;
	};

	template <typename T>
	inline std::ostream &operator<<(std::ostream &os, const Point3<T> &v) {
		os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
		return os;
	}

	template <>
	inline std::ostream &operator<<(std::ostream &os, const Point3<Float> &v) {
		os << StringPrintf("[ %f, %f, %f ]", v.x, v.y, v.z);
		return os;
	}

	typedef Point2<Float> Point2f;
	typedef Point2<int> Point2i;
	typedef Point3<Float> Point3f;
	typedef Point3<int> Point3i;


	// Normal Declarations 
	// -------------------
	// Very similar to a Vector3
	// can't take the cross product of two normals, can't add normal to a point
	// NOT necessarily normalized
	template <typename T> class Normal3 {
	public:
		// Normal3 Public Methods
		// Default Constructor
		Normal3() { x = y = z = 0; }
		Normal3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) { DCHECK(!HasNaNs()); }
		Normal3<T> operator-() const { return Normal3(-x, -y, -z); }
		Normal3<T> operator+(const Normal3<T> &n) const {
			DCHECK(!n.HasNaNs());
			return Normal3<T>(x + n.x, y + n.y, z + n.z);
		}

		Normal3<T> &operator+=(const Normal3<T> &n) {
			DCHECK(!n.HasNaNs());
			x += n.x;
			y += n.y;
			z += n.z;
			return *this;
		}
		Normal3<T> operator-(const Normal3<T> &n) const {
			DCHECK(!n.HasNaNs());
			return Normal3<T>(x - n.x, y - n.y, z - n.z);
		}

		Normal3<T> &operator-=(const Normal3<T> &n) {
			DCHECK(!n.HasNaNs());
			x -= n.x;
			y -= n.y;
			z -= n.z;
			return *this;
		}
		bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }
		template <typename U>
		Normal3<T> operator*(U f) const {
			return Normal3<T>(f * x, f * y, f * z);
		}

		template <typename U>
		Normal3<T> &operator*=(U f) {
			x *= f;
			y *= f;
			z *= f;
			return *this;
		}
		template <typename U>
		Normal3<T> operator/(U f) const {
			CHECK_NE(f, 0);
			Float inv = (Float)1 / f;
			return Normal3<T>(x * inv, y * inv, z * inv);
		}

		template <typename U>
		Normal3<T> &operator/=(U f) {
			CHECK_NE(f, 0);
			Float inv = (Float)1 / f;
			x *= inv;
			y *= inv;
			z *= inv;
			return *this;
		}

		Float LengthSquared() const { return x * x + y * y + z * z; }
		Float Length() const { return std::sqrt(LengthSquared()); }

#ifndef NDEBUG
		Normal3<T>(const Normal3<T> &n) {
			DCHECK(!n.HasNaNs());
			x = n.x;
			y = n.y;
			z = n.z;
		}

		Normal3<T> &operator=(const Normal3<T> &n) {
			DCHECK(!n.HasNaNs());
			x = n.x;
			y = n.y;
			z = n.z;
			return *this;
		}
#endif  // !NDEBUG
		// Initializes a Normal3 from a Vector3
		explicit Normal3<T>(const Vector3<T> &v) : x(v.x), y(v.y), z(v.z) {
			DCHECK(!v.HasNaNs());
		}
		bool operator==(const Normal3<T> &n) const {
			return x == n.x && y == n.y && z == n.z;
		}
		bool operator!=(const Normal3<T> &n) const {
			return x != n.x || y != n.y || z != n.z;
		}

		T operator[](int i) const {
			DCHECK(i >= 0 && i <= 2);
			if (i == 0) return x;
			if (i == 1) return y;
			return z;
		}

		T &operator[](int i) {
			DCHECK(i >= 0 && i <= 2);
			if (i == 0) return x;
			if (i == 1) return y;
			return z;
		}

		
		// Normal3 Public Data 
		T x, y, z;
	};

	template <typename T>
	inline std::ostream &operator<<(std::ostream &os, const Normal3<T> &v) {
		os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
		return os;
	}

	template <>
	inline std::ostream &operator<<(std::ostream &os, const Normal3<Float> &v) {
		os << StringPrintf("[ %f, %f, %f ]", v.x, v.y, v.z);
		return os;
	}


	typedef Normal3<Float> Normal3f;


	// Bounds Declarations
	// -------------------
	// Represent the extent of these sorts of regions.
	// Both parameterized by a Type T (used to represent the coordinates of its extents)
	template <typename T>
	class Bounds2 {
		public:
		// Bounds2 Public Methods
		// Default Constructor 
		Bounds2() {
			T minNum = std::numeric_limits<T>::lowest();
			T maxNum = std::numeric_limits<T>::max();
			pMin = Point2<T>(maxNum, maxNum);
			pMax = Point2<T>(minNum, minNum);
		}
		// Initialize a Bounds2 to enclose a single point
		explicit Bounds2(const Point2<T> &p) : pMin(p), pMax(p) {}
		Bounds2(const Point2<T> &p1, const Point2<T> &p2) {
			pMin = Point2<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y));
			pMax = Point2<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y));
		}
		template <typename U>
		explicit operator Bounds2<U>() const {
			return Bounds2<U>((Point2<U>)pMin, (Point2<U>)pMax);
		}

		// returns the vector along the box diagonal from the minimum point to the maximum point
		Vector2<T> Diagonal() const { return pMax - pMin; }
		// returns the area of the box
		T Area() const {
			Vector2<T> d = pMax - pMin;
			return (d.x * d.y);
		}

		// returns the index of which of the three axes is longest
		int MaximumExtent() const {
			Vector2<T> diag = Diagonal();
			if (diag.x > diag.y)
				return 0;
			else
				return 1;
		}

		// Array indexing to select between the two points at the corners of the box
		inline const Point2<T> &operator[](int i) const {
			DCHECK(i == 0 || i == 1);
			return (i == 0) ? pMin : pMax;
		}
		inline Point2<T> &operator[](int i) {
			DCHECK(i == 0 || i == 1);
			return (i == 0) ? pMin : pMax;
		}

		bool operator==(const Bounds2<T> &b) const {
			return b.pMin == pMin && b.pMax == pMax;
		}
		bool operator!=(const Bounds2<T> &b) const {
			return b.pMin != pMin || b.pMax != pMax;
		}

		// linearly interpolates between the corners of the box by the given amount in each dimension
		Point2<T> Lerp(const Point2f &t) const {
			return Point2<T>(pbrt::Lerp(t.x, pMin.x, pMax.x),
				pbrt::Lerp(t.y, pMin.y, pMax.y));
		}

		// returns the continuous positions of a relative to the corners box, where a point at the minimum corner has 
		// offset (0, 0), a point at the maximum corner has offset (1, 1) an so forth
		Vector2<T> Offset(const Point2<T> &p) const {
			Vector2<T> o = p - pMin;
			if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
			if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
			return o;
		}

		// returns the the center and radius of a circle that bounds the bounding box
		// This may give a far looser fit than a circle that bounded the original contents of the Bounds2 directly.
		void BoundingSphere(Point2<T> *c, Float *rad) const {
			*c = (pMin + pMax) / 2;
			*rad = Inside(*c, *this) ? Distance(*c, pMax) : 0;
		}
		
		friend std::ostream &operator<<(std::ostream &os, const Bounds2<T> &b) {
			os << "[ " << b.pMin << " - " << b.pMax << " ]";
			return os;
		}

		// Bounds2 Public Data
		Point2<T> pMin, pMax;
	};

	template <typename T>
	class Bounds3 {
		public:
		// Bounds3 Public Methods
		// Default Constructor
		Bounds3() {
			T minNum = std::numeric_limits<T>::lowest();
			T maxNum = std::numeric_limits<T>::max();
			pMin = Point3<T>(maxNum, maxNum, maxNum);
			pMax = Point3<T>(minNum, minNum, minNum);
		
		}
		// Initialize a Bounds3 to enclose a single point
		explicit Bounds3(const Point3<T> &p) : pMin(p), pMax(p) {}
		Bounds3(const Point3<T> &p1, const Point3<T> &p2) 
			: pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y), 
				std::min(p1.z, p2.z)),
			 pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y), 
				std::max(p1.z, p2.z)) { }

		// Array indexing to select between the two points at the corners of the box
		const Point3<T> &operator[](int i) const;
		Point3<T> &operator[](int i);

		bool operator==(const Bounds3<T> &b) const {
			return b.pMin == pMin && b.pMax == pMax;
		}
		bool operator!=(const Bounds3<T> &b) const {
			return b.pMin != pMin || b.pMax != pMax;
		}

		// Returns the coordinantes of one of the eight corners of the bounding box
		Point3<T> Corner(int corner) const {
			DCHECK(corner >= 0 && corner < 8);
			return Point3<T>((*this)[(corner & 1)].x,
							(*this)[(corner & 2) ? 1 : 0].y,
							(*this)[(corner & 4) ? 1 : 0].z);
		}


		// returns the vector along the box diagonal from the minimum point to the maximum point
		Vector3<T> Diagonal() const { return pMax - pMin; }
		// returns the surface area of the box
		T SurfaceArea() const {
			Vector3<T> d = Diagonal();
			return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
		}
		// returns the volume inside of the box
		T Volume() const {
			Vector3<T> d = Diagonal();
			return d.x * d.y * d.z;
		}

		// returns the index of which of the three axes is longest
		int MaximumExtent() const {
			Vector3<T> d = Diagonal();
			if (d.x > d.y && d.x > d.z)
				return 0;
			else if (d.y > d.z)
				return 1;
			else
				return 2;
		}

		// linearly interpolates between the corners of the box by the given amount in each dimension
		Point3<T> Lerp(const Point3f &t) const {
			return Point3<T>(pbrt::Lerp(t.x, pMin.x, pMax.x),
				pbrt::Lerp(t.y, pMin.y, pMax.y),
				pbrt::Lerp(t.z, pMin.z, pMax.z));
		}

		// returns the continuous positions of a relative to the corners box, where a point at the minimum corner has 
		// offset (0, 0, 0), a point at the maximum corner has offset (1, 1, 1) an so forth
		Vector3<T> Offset(const Point3<T> &p) const {
			Vector3<T> o = p - pMin;
			if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
			if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
			if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;
			return o;
		}

		// returns the the center and radius of a sphere that bounds the bounding box
		// This may give a far looser fit than a sphere that bounded the original contents of the Bounds3 directly.
		void BoundingSphere(Point3<T> *center, Float *radius) const {
			*center = (pMin + pMax) / 2;
			*radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
		}

		template <typename U>
		explicit operator Bounds3<U>() const {
			return Bounds3<U>((Point3<U>)pMin, (Point3<U>)pMax);
		}

		/// <summary>
		/// Checks for a ray-box intersection and returns the two parametric t values of the intersection if any.
		/// </summary>
		/// <param name="ray">The ray.</param>
		/// <param name="hitt0">The hitt0.</param>
		/// <param name="hitt1">The hitt1.</param>
		/// <returns>returns true if the intersections's parametric range is returned in the optional arguments hitt 0 and
		//  hitt1.</returns>
		bool IntersectP(const Ray &ray,
						Float *hitt0 = nullptr,
						Float *hitt1 = nullptr) const;

		/// <summary>
		/// Specialized method that takes the reciprocal of the ray's direction as an additional parameter, so that 
		/// the three reciprocals don't need to be computed each time IntersectP() is called
		/// </summary>
		/// <param name="ray">The ray.</param>
		/// <param name="invDir">Precomputed value, indicates whether each direction component is negative (avoids comparisons of the computed tNear and tFar values)</param>
		/// <param name="dirIsNeg">Precomputed value, indicates whether each direction component is negative (avoids comparisons of the computed tNear and tFar values)</param>
		/// <returns>Returns true if the ray segment is entirely inside the bounding box, even if the intersections are not within the ray's (0, tMax) range</returns>
		inline bool IntersectP(const Ray &ray,
								const Vector3f &invDir,
								const int dirIsNeg[3]) const;

		friend std::ostream &operator<<(std::ostream &os, const Bounds3<T> &b) {
			os << "[ " << b.pMin << " - " << b.pMax << " ]";
			return os;
		}

		// Bounds2 Public Data
		Point3<T> pMin, pMax;
	};

	typedef Bounds2<Float> Bounds2f;
	typedef Bounds2<int> Bounds2i;
	typedef Bounds3<Float> Bounds3f;
	typedef Bounds3<int> Bounds3i;
	

	// Bounds2iIterator Declarations
	// ---------------------------
	class Bounds2iIterator : public std::forward_iterator_tag {
		public:
			Bounds2iIterator(const Bounds2i &b, const Point2i &pt)
				: p(pt), bounds(&b) {}
			Bounds2iIterator operator++() {
				advance();
				return *this;
			}
			Bounds2iIterator operator++(int) {
				Bounds2iIterator old = *this;
				advance();
				return old;
			}
			bool operator==(const Bounds2iIterator &bi) const {
				return p == bi.p && bounds == bi.bounds;
			}
			bool operator!=(const Bounds2iIterator &bi) const {
				return p != bi.p || bounds != bi.bounds;
			}

			Point2i operator*() const { return p; }

		private:
			void advance() {
				++p.x;
				if (p.x == bounds->pMax.x) {
					p.x = bounds->pMin.x;
					++p.y;
				}
			}
			Point2i p;
			const Bounds2i *bounds;
	};


	// Ray Declarations
	// ----------------
	class Ray {	// semi-infinite line specified by its origin(Point3f) and direction(Vector3f).

		public: 
			// Ray public Methods
			// Default Constructor
			Ray() : tMax(Infinity), time(0.f), medium(nullptr) {}
			Ray(const Point3f &o, const Vector3f &d, Float tMax = Infinity,
				Float time = 0.f, const Medium *medium = nullptr)
				: o(o), d(d), tMax(tMax), time(time), medium(medium) {}
			// parametric form of a ray: r(t) = o + td	0<=t<INFINITY
			Point3f operator()(Float t) const { return o + d * t; }
			bool HasNaNs() const { return (o.HasNaNs() || d.HasNaNs() || isNaN(tMax)); }
			friend std::ostream &operator<<(std::ostream &os, const Ray &r) {
				os << "[o=" << r.o << ", d=" << r.d << ", tMax=" << r.tMax
					<< ", time=" << r.time << "]";
				return os;
			}

			// Ray public Data
			Point3f o;				// origin
			Vector3f d;				// direction
			mutable Float tMax;		// limits the ray to a segment along its infinite extent
			Float time;				// time value of the ray, in scenes with animated objects, the rendering system constructs a representation of the scene at the appropriate time for each ray
			const Medium *medium;	// record of the medium containing its origin
	};

	// Subset of Ray
	// contains additional information about two auxiliary rays
	// represent: camera rays offset by one sample in the x and y direction from the main ray on the film plane
	class RayDifferential : public Ray {
	public:
		// RayDifferential Public Methods
		// constructors mirror the Ray's constructors
		RayDifferential() { hasDifferentials = false; }
		RayDifferential(const Point3f &o, const Vector3f &d, Float tMax = Infinity,
			Float time = 0.f, const Medium *medium = nullptr)
			: Ray(o, d, tMax, time, medium) {
			hasDifferentials = false;
		}
		// constructor to create RayDifferentials from Rays
		// sets hasDifferentials to false initially because the neighboring rays, if any, are not known
		RayDifferential(const Ray &ray) : Ray(ray) { hasDifferentials = false; }
		bool HasNaNs() const {
			return Ray::HasNaNs() ||
				(hasDifferentials &&
				(rxOrigin.HasNaNs() || ryOrigin.HasNaNs() ||
					rxDirection.HasNaNs() || ryDirection.HasNaNs()));
		}
		// update differential rays for an estimated sample spacing of s
		void ScaleDifferentials(Float s) {
			rxOrigin = o + (rxOrigin - o) * s;
			ryOrigin = o + (ryOrigin - o) * s;
			rxDirection = d + (rxDirection - d) * s;
			ryDirection = d + (ryDirection - d) * s;
		}
		friend std::ostream &operator<<(std::ostream &os, const RayDifferential &r) {
			os << "[ " << (Ray &)r << " has differentials: " <<
				(r.hasDifferentials ? "true" : "false") << ", xo = " << r.rxOrigin <<
				", xd = " << r.rxDirection << ", yo = " << r.ryOrigin << ", yd = " <<
				r.ryDirection;
			return os;
		}

		// RayDifferential Public Data
		bool hasDifferentials;
		Point3f rxOrigin, ryOrigin;
		Vector3f rxDirection, ryDirection;
	};



	// Geometry Inline Functions
	// -------------------------
	template <typename T>
	inline Vector3<T>::Vector3(const Point3<T> &p)
		: x(p.x), y(p.y), z(p.z) {
		DCHECK(!HasNaNs());
	}

	template <typename T, typename U>
	inline Vector3<T> operator*(U s, const Vector3<T> &v) {
		return v * s;
	}

	// returns a vector with the absolute value operation applied to its components
	template <typename T>
	Vector3<T> Abs(const Vector3<T> &v) {
		return Vector3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
	}

	// returns the Dot product of two vectors
	template <typename T>
	inline T Dot(const Vector3<T> &v1, const Vector3<T> &v2) {
		DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
		return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	}

	// returns the absolute value of the dot product
	template <typename T>
	inline T AbsDot(const Vector3<T> &v1, const Vector3<T> &v2) {
		DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
		return std::abs(Dot(v1, v2));
	}

	// returns the cross-product of two vectors
	template <typename T>
	inline Vector3<T> Cross(const Vector3<T> &v1, const Vector3<T> &v2) {
		DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
		double v1x = v1.x, v1y = v1.y, v1z = v1.z;
		double v2x = v2.x, v2y = v2.y, v2z = v2.z;
		return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
			(v1x * v2y) - (v1y * v2x));
	}

	template <typename T>
	inline Vector3<T> Cross(const Vector3<T> &v1, const Normal3<T> &v2) {
		DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
		double v1x = v1.x, v1y = v1.y, v1z = v1.z;
		double v2x = v2.x, v2y = v2.y, v2z = v2.z;
		return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
			(v1x * v2y) - (v1y * v2x));
	}

	template <typename T>
	inline Vector3<T> Cross(const Normal3<T> &v1, const Vector3<T> &v2) {
		DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
		double v1x = v1.x, v1y = v1.y, v1z = v1.z;
		double v2x = v2.x, v2y = v2.y, v2z = v2.z;
		return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
			(v1x * v2y) - (v1y * v2x));
	}

	// returns the normalized (unit-length) vector (by dividing each component by the length of the vector)
	template <typename T>
	inline Vector3<T> Normalize(const Vector3<T> &v) { return v / v.Length(); }

	// returns the smallest coordinate value
	template <typename T>
	T MinComponent(const Vector3<T> &v) {
		return std::min(v.x, std::min(v.y, v.z));
	}

	// returns the largest coordinate value
	template <typename T>
	T MaxComponent(const Vector3<T> &v) {
		return std::max(v.x, std::max(v.y, v.z));
	}

	// returns the index of the component with the largest values
	template <typename T>
	int MaxDimension(const Vector3<T> &v) {
		return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
	}

	// returns component-wise minimum
	template <typename T>
	Vector3<T> Min(const Vector3<T> &p1, const Vector3<T> &p2) {
		return Vector3<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z));
	}

	// returns component-wise maximum
	template <typename T>
	Vector3<T> Max(const Vector3<T> &p1, const Vector3<T> &p2) {
		return Vector3<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z));
	}

	// returns the permutation of the coordinante values according to the index values provided
	template <typename T>
	Vector3<T> Permute(const Vector3<T> &v, int x, int y, int z) {
		return Vector3<T>(v[x], v[y], v[z]);
	}
	
	/// <summary>
	/// Constructs a local coordinate system given only a single 3D Vector
	/// </summary>
	/// <param name="v1">The v1 - (normalized) Vector used to construct local coordinate system</param>
	/// <param name="v2">The v2 - constructed perpendicular vector</param>
	/// <param name="v3">The v3 - constructed perpendicular vector</param>
	template <typename T>
	inline void CoordinateSystem(const Vector3<T> &v1,
										Vector3<T> *v2,
										Vector3<T> *v3) {
		if (std::abs(v1.x) > std::abs(v1.y))
			*v2 = Vector3<T>(-v1.z, 0, v1.x) /
			std::sqrt(v1.x * v1.x + v1.z * v1.z);
		else
			*v2 = Vector3<T>(0, v1.z, -v1.y) /
			std::sqrt(v1.y * v1.y + v1.z * v1.z);
		*v3 = Cross(v1, *v2);
	}



	template <typename T>
	Vector2<T>::Vector2(const Point2<T> &p)
		: x(p.x), y(p.y) {
		DCHECK(!HasNaNs());
	}

	template <typename T>
	Vector2<T>::Vector2(const Point3<T> &p)
		: x(p.x), y(p.y) {
		DCHECK(!HasNaNs());
	}

	template <typename T, typename U>
	inline Vector2<T> operator*(U f, const Vector2<T> &v) {
		return v * f;
	}
	
	// returns the Dot product of two vectors
	template <typename T>
	inline Float Dot(const Vector2<T> &v1, const Vector2<T> &v2) {
		DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
		return v1.x * v2.x + v1.y * v2.y;
	}

	// returns the absolute value of the dot product
	template <typename T>
	inline Float AbsDot(const Vector2<T> &v1, const Vector2<T> &v2) {
		DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
		return std::abs(Dot(v1, v2));
	}

	// returns the normalized (unit-length) vector (by dividing each component by the length of the vector)
	template <typename T>
	inline Vector2<T> Normalize(const Vector2<T> &v) {
		return v / v.Length();
	}

	// returns a vector with the absolute value operation applied to its components
	template <typename T>
	Vector2<T> Abs(const Vector2<T> &v) {
		return Vector2<T>(std::abs(v.x), std::abs(v.y));
	}



	// Calculates distance between two points by subtracting them to compute
	// the vector between them and then finding the length of that vector
	template <typename T>
	inline Float Distance(const Point3<T> &p1, const Point3<T> &p2) {
		return (p1 - p2).Length();
	}

	template <typename T>
	inline Float DistanceSquared(const Point3<T> &p1, const Point3<T> &p2) {
		return (p1 - p2).LengthSquared();
	}

	template <typename T, typename U>
	inline Point3<T> operator*(U f, const Point3<T> &p) {
		DCHECK(!p.HasNaNs());
		return p * f;
	}

	/// <summary>
	/// Linearly interpolate between two points between them at other values of t 
	/// Returns p0 at t == 0, p1 at t == 1
	/// Extrapolates for t<0 or t>1
	/// </summary>
	/// <param name="t">The t.</param>
	/// <param name="p0">The p0.</param>
	/// <param name="p1">The p1.</param>
	/// <returns></returns>
	template <typename T>
	Point3<T> Lerp(Float t, const Point3<T> &p0, const Point3<T> &p1) {
		return (1 - t) * p0 + t * p1;
	}

	// Return points representing the component-wise minimums and maximums of the 2 given points
	template <typename T>
	Point3<T> Min(const Point3<T> &p1, const Point3<T> &p2) {
		return Point3<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
			std::min(p1.z, p2.z));
	}

	template <typename T>
	Point3<T> Max(const Point3<T> &p1, const Point3<T> &p2) {
		return Point3<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
			std::max(p1.z, p2.z));
	}

	template <typename T>
	Point3<T> Floor(const Point3<T> &p) {
		return Point3<T>(std::floor(p.x), std::floor(p.y), std::floor(p.z));
	}

	template <typename T>
	Point3<T> Ceil(const Point3<T> &p) {
		return Point3<T>(std::ceil(p.x), std::ceil(p.y), std::ceil(p.z));
	}

	template <typename T>
	Point3<T> Abs(const Point3<T> &p) {
		return Point3<T>(std::abs(p.x), std::abs(p.y), std::abs(p.z));
	}

	template <typename T>
	inline Float Distance(const Point2<T> &p1, const Point2<T> &p2) {
		return (p1 - p2).Length();
	}

	template <typename T>
	inline Float DistanceSquared(const Point2<T> &p1, const Point2<T> &p2) {
		return (p1 - p2).LengthSquared();
	}

	template <typename T, typename U>
	inline Point2<T> operator*(U f, const Point2<T> &p) {
		DCHECK(!p.HasNaNs());
		return p * f;
	}

	template <typename T>
	Point2<T> Floor(const Point2<T> &p) {
		return Point2<T>(std::floor(p.x), std::floor(p.y));
	}

	template <typename T>
	Point2<T> Ceil(const Point2<T> &p) {
		return Point2<T>(std::ceil(p.x), std::ceil(p.y));
	}

	/// <summary>
	/// Linearly interpolate between two points between them at other values of t 
	/// Returns p0 at t == 0, p1 at t == 1
	/// Extrapolates for t<0 or t>1
	/// </summary>
	/// <param name="t">The t.</param>
	/// <param name="p0">The p0.</param>
	/// <param name="p1">The p1.</param>
	template <typename T>
	Point2<T> Lerp(Float t, const Point2<T> &v0, const Point2<T> &v1) {
		return (1 - t) * v0 + t * v1;
	}

	// Return points representing the component-wise minimums and maximums of the 2 given points
	template <typename T>
	Point2<T> Min(const Point2<T> &pa, const Point2<T> &pb) {
		return Point2<T>(std::min(pa.x, pb.x), std::min(pa.y, pb.y));
	}

	template <typename T>
	Point2<T> Max(const Point2<T> &pa, const Point2<T> &pb) {
		return Point2<T>(std::max(pa.x, pb.x), std::max(pa.y, pb.y));
	}

	// Permutes the coordinate values according to the provided permutation.
	template <typename T>
	Point3<T> Permute(const Point3<T> &p, int x, int y, int z) {
		return Point3<T>(p[x], p[y], p[z]);
	}


	template <typename T, typename U>
	inline Normal3<T> operator*(U f, const Normal3<T> &n) {
		return Normal3<T>(f * n.x, f * n.y, f * n.z);
	}

	template <typename T>
	inline Normal3<T> Normalize(const Normal3<T> &n) {
		return n / n.Length();
	}

	// Constructs a Vector3 from a Normal3
	template <typename T>
	inline Vector3<T>::Vector3(const Normal3<T> &n)
		: x(n.x), y(n.y), z(n.z) {
		DCHECK(!n.HasNaNs());
	}

	// Dot product calculation
	template <typename T>
	inline T Dot(const Normal3<T> &n1, const Vector3<T> &v2) {
		DCHECK(!n1.HasNaNs() && !v2.HasNaNs());
		return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
	}

	template <typename T>
	inline T Dot(const Vector3<T> &v1, const Normal3<T> &n2) {
		DCHECK(!v1.HasNaNs() && !n2.HasNaNs());
		return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
	}

	template <typename T>
	inline T Dot(const Normal3<T> &n1, const Normal3<T> &n2) {
		DCHECK(!n1.HasNaNs() && !n2.HasNaNs());
		return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
	}

	// Abs of Dot product calculation
	template <typename T>
	inline T AbsDot(const Normal3<T> &n1, const Vector3<T> &v2) {
		DCHECK(!n1.HasNaNs() && !v2.HasNaNs());
		return std::abs(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
	}

	template <typename T>
	inline T AbsDot(const Vector3<T> &v1, const Normal3<T> &n2) {
		DCHECK(!v1.HasNaNs() && !n2.HasNaNs());
		return std::abs(v1.x * n2.x + v1.y * n2.y + v1.z * n2.z);
	}

	template <typename T>
	inline T AbsDot(const Normal3<T> &n1, const Normal3<T> &n2) {
		DCHECK(!n1.HasNaNs() && !n2.HasNaNs());
		return std::abs(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
	}

	// Flip a surface normal so that it lies in the same hemisphere as a 
	// given vector 
	// First parameter is the one that should be returned (possibly flipped)
	// and second is the one to test against
	template <typename T>
	inline Normal3<T> Faceforward(const Normal3<T> &n, const Vector3<T> &v) {
		return (Dot(n, v) < 0.f) ? -n : n;
	}

	template <typename T>
	inline Normal3<T> Faceforward(const Normal3<T> &n, const Normal3<T> &n2) {
		return (Dot(n, n2) < 0.f) ? -n : n;
	}

	template <typename T>
	inline Vector3<T> Faceforward(const Vector3<T> &v, const Vector3<T> &v2) {
		return (Dot(v, v2) < 0.f) ? -v : v;
	}

	template <typename T>
	inline Vector3<T> Faceforward(const Vector3<T> &v, const Normal3<T> &n2) {
		return (Dot(v, n2) < 0.f) ? -v : v;
	}

	template <typename T>
	Normal3<T> Abs(const Normal3<T> &v) {
		return Normal3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
	}

	template <typename T>
	inline const Point3<T> &Bounds3<T>::operator[](int i) const {
		DCHECK(i == 0 || i == 1);
		return (i == 0) ? pMin : pMax;
	}

	template <typename T>
	inline Point3<T> &Bounds3<T>::operator[](int i) {
		DCHECK(i == 0 || i == 1);
		return (i == 0) ? pMin : pMax;
	}

	// returns a new bounding box that encompasses that point as well as the original box
	template <typename T>
	Bounds3<T> Union(const Bounds3<T> &b, const Point3<T> &p) {
		Bounds3<T> ret;
		ret.pMin = Min(b.pMin, p);
		ret.pMax = Max(b.pMax, p);
		return ret;
	}

	// returns a new bounding box that encompasses two bounding boxes
	template <typename T>
	Bounds3<T> Union(const Bounds3<T> &b1, const Bounds3<T> &b2) {
		Bounds3<T> ret;
		ret.pMin = Min(b1.pMin, b2.pMin);
		ret.pMax = Max(b1.pMax, b2.pMax);
		return ret;
	}

	// returns the intersection bounding box of two bounding boxes
	template <typename T>
	Bounds3<T> Intersect(const Bounds3<T> &b1, const Bounds3<T> &b2) {
		// Important: assign to pMin/pMax directly and don't run the Bounds2()
		// constructor, since it takes min/max of the points passed to it.  In
		// turn, that breaks returning an invalid bound for the case where we
		// intersect non-overlapping bounds (as we'd like to happen).
		Bounds3<T> ret;
		ret.pMin = Max(b1.pMin, b2.pMin);
		ret.pMax = Min(b1.pMax, b2.pMax);
		return ret;
	}

	// returns true if two bounding boxes overlap
	template <typename T>
	bool Overlaps(const Bounds3<T> &b1, const Bounds3<T> &b2) {
		bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
		bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
		bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
		return (x && y && z);
	}

	// returns true if a point is inside (or on the boundary of) the bounding box (based on three 1D containment tests)
	template <typename T>
	bool Inside(const Point3<T> &p, const Bounds3<T> &b) {
		return (p.x >= b.pMin.x && p.x <= b.pMax.x && p.y >= b.pMin.y &&
			p.y <= b.pMax.y && p.z >= b.pMin.z && p.z <= b.pMax.z);
	}

	// returns true if a point is inside the bounding box (based on three 1D containment tests)
	// if its on the boundary or outside it returns false
	template <typename T>
	bool InsideExclusive(const Point3<T> &p, const Bounds3<T> &b) {
		return (p.x >= b.pMin.x && p.x < b.pMax.x && p.y >= b.pMin.y &&
			p.y < b.pMax.y && p.z >= b.pMin.z && p.z < b.pMax.z);
	}

	// pads the bounding box by a constant factor in all dimensions
	template <typename T, typename U>
	inline Bounds3<T> Expand(const Bounds3<T> &b, U delta) {
		return Bounds3<T>(b.pMin - Vector3<T>(delta, delta, delta),
			b.pMax + Vector3<T>(delta, delta, delta));
	}

	// Minimum squared distance from point to box; returns zero if point is
	// inside.
	template <typename T, typename U>
	inline Float DistanceSquared(const Point3<T> &p, const Bounds3<U> &b) {
		Float dx = std::max({ Float(0), b.pMin.x - p.x, p.x - b.pMax.x });
		Float dy = std::max({ Float(0), b.pMin.y - p.y, p.y - b.pMax.y });
		Float dz = std::max({ Float(0), b.pMin.z - p.z, p.z - b.pMax.z });
		return dx * dx + dy * dy + dz * dz;
	}

	template <typename T, typename U>
	inline Float Distance(const Point3<T> &p, const Bounds3<U> &b) {
		return std::sqrt(DistanceSquared(p, b));
	}

	inline Bounds2iIterator begin(const Bounds2i &b) {
		return Bounds2iIterator(b, b.pMin);
	}

	inline Bounds2iIterator end(const Bounds2i &b) {
		// Normally, the ending point is at the minimum x value and one past
		// the last valid y value.
		Point2i pEnd(b.pMin.x, b.pMax.y);
		// However, if the bounds are degenerate, override the end point to
		// equal the start point so that any attempt to iterate over the bounds
		// exits out immediately.
		if (b.pMin.x >= b.pMax.x || b.pMin.y >= b.pMax.y)
			pEnd = b.pMin;
		return Bounds2iIterator(b, pEnd);
	}

	// returns a new bounding box that encompasses two bounding boxes
	template <typename T>
	Bounds2<T> Union(const Bounds2<T> &b, const Point2<T> &p) {
		Bounds2<T> ret;
		ret.pMin = Min(b.pMin, p);
		ret.pMax = Max(b.pMax, p);
		return ret;
	}

	template <typename T>
	Bounds2<T> Union(const Bounds2<T> &b, const Bounds2<T> &b2) {
		Bounds2<T> ret;
		ret.pMin = Min(b.pMin, b2.pMin);
		ret.pMax = Max(b.pMax, b2.pMax);
		return ret;
	}

	// returns the intersection bounding box of two bounding boxes
	template <typename T>
	Bounds2<T> Intersect(const Bounds2<T> &b1, const Bounds2<T> &b2) {
		// Important: assign to pMin/pMax directly and don't run the Bounds2()
		// constructor, since it takes min/max of the points passed to it.  In
		// turn, that breaks returning an invalid bound for the case where we
		// intersect non-overlapping bounds (as we'd like to happen).
		Bounds2<T> ret;
		ret.pMin = Max(b1.pMin, b2.pMin);
		ret.pMax = Min(b1.pMax, b2.pMax);
		return ret;
	}
	
	// returns true if two bounding boxes overlap
	template <typename T>
	bool Overlaps(const Bounds2<T> &ba, const Bounds2<T> &bb) {
		bool x = (ba.pMax.x >= bb.pMin.x) && (ba.pMin.x <= bb.pMax.x);
		bool y = (ba.pMax.y >= bb.pMin.y) && (ba.pMin.y <= bb.pMax.y);
		return (x && y);
	}

	// returns true if a point is inside (or on the boundary of) the bounding box (based on two 1D containment tests)
	template <typename T>
	bool Inside(const Point2<T> &pt, const Bounds2<T> &b) {
		return (pt.x >= b.pMin.x && pt.x <= b.pMax.x && pt.y >= b.pMin.y &&
			pt.y <= b.pMax.y);
	}

	// returns true if a point is inside the bounding box (based on two 1D containment tests)
	// if its on the boundary or outside it returns false
	template <typename T>
	bool InsideExclusive(const Point2<T> &pt, const Bounds2<T> &b) {
		return (pt.x >= b.pMin.x && pt.x < b.pMax.x && pt.y >= b.pMin.y &&
			pt.y < b.pMax.y);
	}

	// pads the bounding box by a constant factor in all dimensions
	template <typename T, typename U>
	Bounds2<T> Expand(const Bounds2<T> &b, U delta) {
		return Bounds2<T>(b.pMin - Vector2<T>(delta, delta),
			b.pMax + Vector2<T>(delta, delta));
	}

	template <typename T>
	inline bool Bounds3<T>::IntersectP(const Ray &ray,
										Float *hitt0,
										Float *hitt1) const {
		Float t0 = 0, t1 = ray.tMax; // (0, Ray::tMax) range others are ignored
		for (int i = 0; i < 3; ++i) { // three slabs of the bounding box
			// Update interval for _i_th bounding box slab
			Float invRayDir = 1 / ray.d[i];		// computing the reciprocal of the corresponding component of the ray direction
			Float tNear = (pMin[i] - ray.o[i]) * invRayDir;	// compute t values of slab interscetion
			Float tFar = (pMax[i] - ray.o[i]) * invRayDir;	// compute t values of slab interscetion

			// Update parametric interval from slab intersection $t$ values
			if (tNear > tFar) std::swap(tNear, tFar); // reorder so that tNear holds the closer intersection and tFar the farther one

			// Update _tFar_ to ensure robust ray--bounds intersection
			// in case where the ray origin is in the plane of one of the bounding box slabs and the ray lies in the plane of the slab it is
			// possible that tNear or tFar will be computed by an expression of the form 0/0 (NaN)
			tFar *= 1 + 2 * gamma(3);
			t0 = tNear > t0 ? tNear : t0;	// avoid NAN on t0 and t1
			t1 = tFar < t1 ? tFar : t1;
			if (t0 > t1) return false;
		}
		if (hitt0) *hitt0 = t0;
		if (hitt1) *hitt1 = t1;
		return true;
	}

	template <typename T>
	inline bool Bounds3<T>::IntersectP(const Ray &ray, const Vector3f &invDir,
										const int dirIsNeg[3]) const {
		// provides a 15% performance improvement in overall rendering time compared to using the other variant
		const Bounds3f &bounds = *this;
		// Check for ray intersection against $x$ and $y$ slabs
		// if the ray dir vector is negative, the near parametric intersection will be found with the slab
		// with the larger of the two bounding values
		// the far intersection will be found with the slab with the smaller of them
		// ==> directly compute the near and far parametric values in each direction
		Float tMin = (bounds[dirIsNeg[0]].x - ray.o.x) * invDir.x;
		Float tMax = (bounds[1 - dirIsNeg[0]].x - ray.o.x) * invDir.x;
		Float tyMin = (bounds[dirIsNeg[1]].y - ray.o.y) * invDir.y;
		Float tyMax = (bounds[1 - dirIsNeg[1]].y - ray.o.y) * invDir.y;

		// Update _tMax_ and _tyMax_ to ensure robust bounds intersection
		tMax *= 1 + 2 * gamma(3);
		tyMax *= 1 + 2 * gamma(3);
		if (tMin > tyMax || tyMin > tMax) return false;
		if (tyMin > tMin) tMin = tyMin;
		if (tyMax < tMax) tMax = tyMax;

		// Check for ray intersection against $z$ slab
		Float tzMin = (bounds[dirIsNeg[2]].z - ray.o.z) * invDir.z;
		Float tzMax = (bounds[1 - dirIsNeg[2]].z - ray.o.z) * invDir.z;

		// Update _tzMax_ to ensure robust bounds intersection
		tzMax *= 1 + 2 * gamma(3);
		if (tMin > tzMax || tzMin > tMax) return false;
		if (tzMin > tMin) tMin = tzMin;
		if (tzMax < tMax) tMax = tzMax;
		return (tMin < ray.tMax) && (tMax > 0);
	}

	/// <summary>
	/// Offsets the ray origin along the normal by minimizing the distance.
	/// </summary>
	/// <param name="p">The p.</param>
	/// <param name="pError">The p error.</param>
	/// <param name="n">The n.</param>
	/// <param name="w">The w.</param>
	/// <returns></returns>
	inline Point3f OffsetRayOrigin(const Point3f& p, const Vector3f& pError,
									const Normal3f &n, const Vector3f &w) {
		Float d = Dot(Abs(n), pError);
#ifdef PBRT_FLOAT_AS_DOUBLE
		// We have tons of precision; for now bump up the offset a bunch just
		// to be extra sure that we start on the right side of the surface
		// (In case of any bugs in the epsilons code...)
		d *= 1024.;
#endif
		Vector3f offset = d * Vector3f(n);
		if (Dot(w, n) < 0) offset = -offset;
		Point3f po = p + offset;
		// Round offset point _po_ away from _p_
		// advancing each coordinate of the computed point one 
		// floating-point value away from p ensures that it is outside of
		// the error box
		for (int i = 0; i < 3; ++i) {
			if (offset[i] > 0)
				po[i] = NextFloatUp(po[i]);
			else if (offset[i] < 0)
				po[i] = NextFloatDown(po[i]);
		}
		return po;
	}

	/// <summary>
	/// Converts spherical coordinates (theata, phi) to a direction vector.
	/// </summary>
	/// <param name="sinTheta">The sin of the theta angle.</param>
	/// <param name="cosTheta">The cos of the theta angle.</param>
	/// <param name="phi">The phi angle.</param>
	/// <returns>The corresponding direction vector.</returns>
	inline Vector3f SphericalDirection(Float sinTheta, Float cosTheta, Float phi) {
		return Vector3f(sinTheta * std::cos(phi), 
						sinTheta * std::sin(phi),
						cosTheta);
	}

	/// <summary>
	/// Converts basis vectors to the appropriate direction vector with resprect to the coordinate drame defined by them.
	/// </summary>
	/// <param name="sinTheta">The sin theta.</param>
	/// <param name="cosTheta">The cos theta.</param>
	/// <param name="phi">The phi.</param>
	/// <param name="x">The x basis vector.</param>
	/// <param name="y">The y basis vector.</param>
	/// <param name="z">The z basis vector.</param>
	/// <returns>The appropriate direction vector.</returns>
	inline Vector3f SphericalDirection(Float sinTheta, Float cosTheta, Float phi,
										const Vector3f &x, const Vector3f &y,
										const Vector3f &z) {
		return	sinTheta * std::cos(phi) * x + 
				sinTheta * std::sin(phi) * y +
				cosTheta * z;
	}

	/// <summary>
	/// Converts a dierection to the spherical angle Theta.
	/// Not the function assumes that the vector v has been normalized before being passed in.
	/// </summary>
	/// <param name="v">The direction vector.</param>
	/// <returns></returns>
	inline Float SphericalTheta(const Vector3f& v) {
		return std::acos(Clamp(v.z, -1, 1));	// clamp purley avoids errors from std::acos if |v.z| is slightly greater than 1 due to floating point round-off error
	}

	/// <summary>
	/// Converts a dierection to the spherical angle Phi.
	/// </summary>
	/// <param name="v">The direction vector.</param>
	/// <returns></returns>
	inline Float SphericalPhi(const Vector3f& v) {
		Float p = std::atan2(v.y, v.x);
		return (p < 0) ? (p + 2 * Pi) : p;
	}

} // namespace pbrt
#endif // PBRT_CORE_GEOMETRY_H