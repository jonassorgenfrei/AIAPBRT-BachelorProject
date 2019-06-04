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

#ifndef PBRT_CORE_SPECTRUM_H
#define PBRT_CORE_SPECTRUM_H

 // core/spectrum.h*
#include "pbrt.h"
#include "stringprint.h"

namespace pbrt {

	// Spectrum Utility Declarations
	static const int sampledLambdaStart = 400;	// wavelength range minimum
	static const int sampledLambdaEnd = 700;	// wavelength range maximum
	static const int nSpectralSamples = 60;		// number of samples in the Spectrum


	extern bool SpectrumSamplesSorted(const Float* lambda, const Float* vals,
		int n);
	extern void SortSpectrumSamples(Float* lambda, Float* vals, int n);
	
	extern Float AverageSpectrumSamples(const Float* lambda, const Float* vals,
		int n, Float lambdaStart, Float lambdaEnd);
	inline void XYZToRGB(const Float xyz[3], Float rgb[3]) {
		rgb[0] = 3.240479f * xyz[0] - 1.537150f * xyz[1] - 0.498535f * xyz[2];
		rgb[1] = -0.969256f * xyz[0] + 1.875991f * xyz[1] + 0.041556f * xyz[2];
		rgb[2] = 0.055648f * xyz[0] - 0.204043f * xyz[1] + 1.057311f * xyz[2];
	}

	inline void RGBToXYZ(const Float rgb[3], Float xyz[3]) {
		xyz[0] = 0.412453f * rgb[0] + 0.357580f * rgb[1] + 0.180423f * rgb[2];
		xyz[1] = 0.212671f * rgb[0] + 0.715160f * rgb[1] + 0.072169f * rgb[2];
		xyz[2] = 0.019334f * rgb[0] + 0.119193f * rgb[1] + 0.950227f * rgb[2];
	}

	enum class SpectrumType { Reflectance, Illuminant };
	extern Float InterpolateSpectrumSamples(const Float* lambda, const Float* vals,
		int n, Float l);
	extern void Blackbody(const Float* lambda, int n, Float T, Float* Le);
	extern void BlackbodyNormalized(const Float* lambda, int n, Float T,
		Float* vals);

	// Spectral Data Declarations
	static const int nCIESamples = 471;
	extern const Float CIE_X[nCIESamples];
	extern const Float CIE_Y[nCIESamples];
	extern const Float CIE_Z[nCIESamples];
	extern const Float CIE_lambda[nCIESamples];
	static const Float CIE_Y_integral = 106.856895;
	static const int nRGB2SpectSamples = 32;
	extern const Float RGB2SpectLambda[nRGB2SpectSamples];
	extern const Float RGBRefl2SpectWhite[nRGB2SpectSamples];
	extern const Float RGBRefl2SpectCyan[nRGB2SpectSamples];
	extern const Float RGBRefl2SpectMagenta[nRGB2SpectSamples];
	extern const Float RGBRefl2SpectYellow[nRGB2SpectSamples];
	extern const Float RGBRefl2SpectRed[nRGB2SpectSamples];
	extern const Float RGBRefl2SpectGreen[nRGB2SpectSamples];
	extern const Float RGBRefl2SpectBlue[nRGB2SpectSamples];
	extern const Float RGBIllum2SpectWhite[nRGB2SpectSamples];
	extern const Float RGBIllum2SpectCyan[nRGB2SpectSamples];
	extern const Float RGBIllum2SpectMagenta[nRGB2SpectSamples];
	extern const Float RGBIllum2SpectYellow[nRGB2SpectSamples];
	extern const Float RGBIllum2SpectRed[nRGB2SpectSamples];
	extern const Float RGBIllum2SpectGreen[nRGB2SpectSamples];
	extern const Float RGBIllum2SpectBlue[nRGB2SpectSamples];

	/// Spectrum Declarations

	/// <summary>
	/// Representation of a spectrum as a
	/// particular number of samples given as the nSpectrumSample
	/// template parameter.
	/// Note: implicit assumption that the spectral representation is a set of 
	/// coefficients that linearly scale a fixed set of basis functions.
	/// </summary>
	template <int nSpectrumSamples>
	class CoefficientSpectrum {
	public:
		// CoefficientSpectrum Public Methods

		/// <summary>
		/// Initializes a new instance of the <see cref="CoefficientSpectrum"/> class.
		/// </summary>
		/// <param name="v">Constant value across all wavelengths.</param>
		CoefficientSpectrum(Float v = 0.f) {
			for (int i = 0; i < nSpectrumSamples; ++i) c[i] = v;
			DCHECK(!HasNaNs());
		}

#ifdef DEBUG
		CoefficientSpectrum(const CoefficientSpectrum & s) {
			DCHECK(!s.HasNaNs());
			for (int i = 0; i < nSpectrumSamples; ++i) c[i] = s.c[i];
		}

		CoefficientSpectrum& operator=(const CoefficientSpectrum & s) {
			DCHECK(!s.HasNaNs());
			for (int i = 0; i < nSpectrumSamples; ++i) c[i] = s.c[i];
			return *this;
		}
#endif  // DEBUG

		/// <summary>
		/// Prints the Spectrum on the specified filestream.
		/// </summary>
		/// <param name="f">The filestream</param>
		void Print(FILE* f) const {
			fprintf(f, "[ ");
			for (int i = 0; i < nSpectrumSamples; ++i) {
				fprintf(f, "%f", c[i]);
				if (i != nSpectrumSamples - 1) fprintf(f, ", ");
			}
			fprintf(f, "]");
		}
		
		/// <summary>
		/// Operation to add pairs of spectral distribution.
		/// Component-wise.
		/// </summary>
		/// <param name="s2">The Spectrum to add.</param>
		/// <returns>This</returns>
		CoefficientSpectrum& operator+=(const CoefficientSpectrum& s2) {
			DCHECK(!s2.HasNaNs());
			for (int i = 0; i < nSpectrumSamples; ++i) c[i] += s2.c[i];
			return *this;
		}

		/// <summary>
		/// Operation to add pairs of spectral distribution.
		/// Component-wise.
		/// </summary>
		/// <param name="s2">The Spectrum to add.</param>
		/// <returns>The Result of the add operation.</returns>
		CoefficientSpectrum operator+(const CoefficientSpectrum& s2) const {
			DCHECK(!s2.HasNaNs());
			CoefficientSpectrum ret = *this;
			for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] += s2.c[i];
			return ret;
		}

		/// <summary>
		/// Operator to subtract a pair of spectral distribution.
		/// Component-wise.
		/// </summary>
		/// <param name="s2">The Spectrum to subtract.</param>
		/// <returns>The Result of the subtract operation.</returns>
		CoefficientSpectrum operator-(const CoefficientSpectrum& s2) const {
			DCHECK(!s2.HasNaNs());
			CoefficientSpectrum ret = *this;
			for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] -= s2.c[i];
			return ret;
		}

		/// <summary>
		/// Operator to divide a pair of spectral distribution.
		/// Component-wise.
		/// </summary>
		/// <param name="s2">The Spectrum which is the divisior.</param>
		/// <returns>The Result of the division operation.</returns>
		CoefficientSpectrum operator/(const CoefficientSpectrum & s2) const {
			DCHECK(!s2.HasNaNs());
			CoefficientSpectrum ret = *this;
			for (int i = 0; i < nSpectrumSamples; ++i) {
				CHECK_NE(s2.c[i], 0);
				ret.c[i] /= s2.c[i];
			}
			return ret;
		}

		/// <summary>
		/// Operator to multiplicate a pair of spectral distribution.
		/// Component-wise.
		/// </summary>
		/// <param name="sp">The Spectrum which is the multiplicator.</param>
		/// <returns>The Result of the multiplication operation.</returns>
		CoefficientSpectrum operator*(const CoefficientSpectrum & sp) const {
			DCHECK(!sp.HasNaNs());
			CoefficientSpectrum ret = *this;
			for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] *= sp.c[i];
			return ret;
		}

		/// <summary>
		/// Operator to multiplicate a pair of spectral distribution.
		/// Component-wise.
		/// </summary>
		/// <param name="sp">The Spectrum which is the multiplicator.</param>
		/// <returns>This</returns>
		CoefficientSpectrum& operator*=(const CoefficientSpectrum & sp) {
			DCHECK(!sp.HasNaNs());
			for (int i = 0; i < nSpectrumSamples; ++i) c[i] *= sp.c[i];
			return *this;
		}

		/// <summary>
		/// Operator to multiplicate a spectral distribution with a scalar.
		/// Component-wise.
		/// </summary>
		/// <param name="s2">The scalar to multiply onto.</param>
		/// <returns>The Result of the multiplication operation.</returns>
		CoefficientSpectrum operator*(Float a) const {
			CoefficientSpectrum ret = *this;
			for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] *= a;
			DCHECK(!ret.HasNaNs());
			return ret;
		}

		/// <summary>
		/// Operator to multiplicate a spectral distribution with a scalar.
		/// Component-wise.
		/// </summary>
		/// <param name="s2">The scalar to multiply.</param>
		/// <returns>This.</returns>
		CoefficientSpectrum& operator*=(Float a) {
			for (int i = 0; i < nSpectrumSamples; ++i) c[i] *= a;
			DCHECK(!HasNaNs());
			return *this;
		}

		/// <summary>
		/// Operator to multiplicate a scalar with a spectral distribution.
		/// Component-wise.
		/// </summary>
		/// <param name="a">The scalar</param>
		/// <param name="s">The spectrum</param>
		/// <returns>
		/// The Result of the multiplication operation.
		/// </returns>
		friend inline CoefficientSpectrum operator*(Float a,
			const CoefficientSpectrum & s) {
			DCHECK(!std::isnan(a) && !s.HasNaNs());
			return s * a;
		}


		/// <summary>
		/// Operator to divide the spectral distribution by a scalar.
		/// Component-wise.
		/// </summary>
		/// <param name="a">The scalar (divisior)</param>
		/// <returns>
		/// The result of the division operation.
		/// </returns>
		CoefficientSpectrum operator/(Float a) const {
			CHECK_NE(a, 0);
			DCHECK(!std::isnan(a));
			CoefficientSpectrum ret = *this;
			for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] /= a;
			DCHECK(!ret.HasNaNs());
			return ret;
		}

		/// <summary>
		/// Operator to divide the spectral distribution by a scalar.
		/// Component-wise.
		/// </summary>
		/// <param name="a">The scalar (divisior)</param>
		/// <returns>This.</returns>
		CoefficientSpectrum& operator/=(Float a) {
			CHECK_NE(a, 0);
			DCHECK(!std::isnan(a));
			for (int i = 0; i < nSpectrumSamples; ++i) c[i] /= a;
			return *this;
		}

		/// <summary>
		/// Checks component wise if the specified spectrum is equal to this 
		/// spectrum.
		/// </summary>
		/// <param name="sp">Specific spectrum.</param>
		/// <returns><c>true</c> if the spectra are equal, otherwise <c>false</c></returns>
		bool operator==(const CoefficientSpectrum& sp) const {
			for (int i = 0; i < nSpectrumSamples; ++i)
				if (c[i] != sp.c[i]) return false;
			return true;
		}

		/// <summary>
		/// Checks component wise if the specified spectrum is unequal to this 
		/// spectrum.
		/// </summary>
		/// <param name="sp">Specific spectrum.</param>
		/// <returns><c>true</c> if the spectra are unequal, otherwise <c>false</c></returns>
		bool operator!=(const CoefficientSpectrum & sp) const {
			return !(*this == sp);
		}

		/// <summary>
		/// Determines whether this instance is black.
		/// It Means, the the surface has values zero everywhere.
		/// </summary>
		/// <returns>
		///   <c>true</c> if this instance is black; otherwise, <c>false</c>.
		/// </returns>
		bool IsBlack() const {
			for (int i = 0; i < nSpectrumSamples; ++i)
				if (c[i] != 0.) return false;
			return true;
		}

		/// <summary>
		/// Takes the Square root of spectrum component-wise.
		/// </summary>
		/// <param name="s">The s.</param>
		/// <returns></returns>
		friend CoefficientSpectrum Sqrt(const CoefficientSpectrum& s) {
			CoefficientSpectrum ret;
			for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = std::sqrt(s.c[i]);
			DCHECK(!ret.HasNaNs());
			return ret;
		}

		/// <summary>
		/// Raises the SPD to a given power.
		/// </summary>
		/// <param name="s">The spectrum.</param>
		/// <param name="e">The given power.</param>
		/// <returns>Result of the raise operation.</returns>
		template <int n>
		friend inline CoefficientSpectrum<n> Pow(const CoefficientSpectrum<n>& s,
			Float e);

		/// <summary>
		/// Inverts the spectrum.
		/// Component-wise.
		/// </summary>
		/// <returns>The inverted spectrum.</returns>
		CoefficientSpectrum operator-() const {
			CoefficientSpectrum ret;
			for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = -c[i];
			return ret;
		}

		/// <summary>
		/// Computes e raised to the given power component-wise.
		/// </summary>
		/// <param name="s">The spectrum.</param>
		/// <returns>Result of the e raise operation.</returns>
		friend CoefficientSpectrum Exp(const CoefficientSpectrum& s) {
			CoefficientSpectrum ret;
			for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = std::exp(s.c[i]);
			DCHECK(!ret.HasNaNs());
			return ret;
		}

		/// <summary>
		/// Operators the specified os.
		/// </summary>
		/// <param name="os">The os.</param>
		/// <param name="s">The spectrum.</param>
		/// <returns>The stream of the spectrum written on the osstream.</returns>
		friend std::ostream& operator<<(std::ostream& os,
			const CoefficientSpectrum & s) {
			return os << s.ToString();
		}

		/// <summary>
		/// Converts to the Spectrumto a string.
		/// </summary>
		/// <returns>Spectrum String.</returns>
		std::string ToString() const {
			std::string str = "[ ";
			for (int i = 0; i < nSpectrumSamples; ++i) {
				str += StringPrintf("%f", c[i]);
				if (i + 1 < nSpectrumSamples) str += ", ";
			}
			str += " ]";
			return str;
		}
		
		/// <summary>
		/// Clamps the Spectrum to the given parameters.
		/// </summary>
		/// <param name="low">The low clamp border. Default 0.</param>
		/// <param name="high">The high clamp border.Default Infinity.</param>
		/// <returns>The clamped spectrum.</returns>
		CoefficientSpectrum Clamp(Float low = 0, Float high = Infinity) const {
			CoefficientSpectrum ret;
			for (int i = 0; i < nSpectrumSamples; ++i)
				ret.c[i] = pbrt::Clamp(c[i], low, high);
			DCHECK(!ret.HasNaNs());
			return ret;
		}
		
		/// <summary>
		/// Returns the maximum component value in the spectrum.
		/// </summary>
		/// <returns>Maximum component value in the spectrum.</returns>
		Float MaxComponentValue() const {
			Float m = c[0];
			for (int i = 1; i < nSpectrumSamples; ++i)
				m = std::max(m, c[i]);
			return m;
		}

		/// <summary>
		/// Debuging Routine that determines whether the spectrum has NaNs.
		/// </summary>
		/// <returns>
		///   <c>true</c> if it has NaNs otherwise, <c>false</c>.
		/// </returns>
		bool HasNaNs() const {
			for (int i = 0; i < nSpectrumSamples; ++i)
				if (std::isnan(c[i])) return true;
			return false;
		}

		/// <summary>
		/// Writes the specified file pointer.
		/// </summary>
		/// <param name="f">The file pointer.</param>
		/// <returns>
		///  <c>true</c> if the write was successfull, otherwise <c>false</c>
		///</returns>
		bool Write(FILE* f) const {
			for (int i = 0; i < nSpectrumSamples; ++i)
				if (fprintf(f, "%f ", c[i]) < 0) return false;
			return true;
		}

		/// <summary>
		/// Reads the specified file pointer.
		/// </summary>
		/// <param name="f">The file pointer.</param>
		/// <returns>
		///  <c>true</c> if the write was successfull, otherwise <c>false</c>
		///</returns>
		bool Read(FILE* f) {
			for (int i = 0; i < nSpectrumSamples; ++i) {
				double v;
				if (fscanf(f, "%lf ", &v) != 1) return false;
				c[i] = v;
			}
			return true;
		}

		/// <summary>
		/// Operator to access individual samples values.
		/// </summary>
		/// <param name="i">The index.</param>
		/// <returns>The pointer to the samples at the specified index.</returns>
		Float& operator[](int i) {
			DCHECK(i >= 0 && i < nSpectrumSamples);
			return c[i];
		}
		
		/// <summary>
		/// Operator to access individual samples values.
		/// </summary>
		/// <param name="i">The index.</param>
		/// <returns>The samples at the specified index.</returns>
		Float operator[](int i) const {
			DCHECK(i >= 0 && i < nSpectrumSamples);
			return c[i];
		}

		// CoefficientSpectrum Public Data
		/// <summary>
		/// Gives the number of samples used to represent the SPD
		/// </summary>
		static const int nSamples = nSpectrumSamples;

	protected:
		// CoefficientSpectrum Protected Data
		Float c[nSpectrumSamples];
	};

	/// <summary>
	/// 
	/// It uses the <see cref="Coefficientspectrum"\> infrastructure to represent 
	/// an SPD with uniformly spaced samples between a starting and an ending wavelength.
	/// </summary>
	/// <seealso cref="CoefficientSpectrum{nSpectralSamples}" />
	class SampledSpectrum : public CoefficientSpectrum<nSpectralSamples> {
	public:
		// SampledSpectrum Public Methods

		/// <summary>
		/// Initializes a new instance of the <see cref="SampledSpectrum"/> class.
		/// </summary>
		/// <param name="v">The Spectrum.</param>
		SampledSpectrum(Float v = 0.f) : CoefficientSpectrum(v) {}
		
		/// <summary>
		/// Initializes a new instance of the <see cref="SampledSpectrum"/> class.
		/// </summary>
		/// <param name="v">The Spectrum with the templated amount of Samples.</param>
		SampledSpectrum(const CoefficientSpectrum<nSpectralSamples>& v)
			: CoefficientSpectrum<nSpectralSamples>(v) {}

		/// <summary>
		/// Takes arrays of SPD sample values v at given wavelengths lambda and 
		/// uses them to define a piecewise linear function to represent the SPD.
		/// </summary>
		/// <param name="lambda">The wavelength values.</param>
		/// <param name="v">The value coressponding to the wavelength lambda at the same index..</param>
		/// <param name="n">The number of sample values.</param>
		/// <returns></returns>
		static SampledSpectrum FromSampled(const Float* lambda, const Float* v,
			int n) {
			// Sort samples if unordered, use sorted for returned spectrum
			if (!SpectrumSamplesSorted(lambda, v, n)) {
				// if samples are not sorted allocate new storage
				std::vector<Float> slambda(&lambda[0], &lambda[n]);
				std::vector<Float> sv(&v[0], &v[n]);
				SortSpectrumSamples(&slambda[0], &sv[0], n);
				return FromSampled(&slambda[0], &sv[0], n);
			}
			SampledSpectrum r;
			for (int i = 0; i < nSpectralSamples; ++i) {
				// Compute average value of given SPD over $i$th sample's range
				// 1D instance of sampling and reconstruction
				Float lambda0 = Lerp(Float(i) / Float(nSpectralSamples),
					sampledLambdaStart, sampledLambdaEnd);
				Float lambda1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
					sampledLambdaStart, sampledLambdaEnd);
				// compute the average of the piecewise linear function over the range of
				// wavelengths that each SPD sample is responsible for
				r.c[i] = AverageSpectrumSamples(lambda, v, n, lambda0, lambda1);
			}
			return r;
		}

		static void Init() {
			// Compute XYZ matching functions for _SampledSpectrum_
			for (int i = 0; i < nSpectralSamples; ++i) {
				Float wl0 = Lerp(Float(i) / Float(nSpectralSamples),
					sampledLambdaStart, sampledLambdaEnd);
				Float wl1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
					sampledLambdaStart, sampledLambdaEnd);
				X.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples, wl0,
					wl1);
				Y.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples, wl0,
					wl1);
				Z.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples, wl0,
					wl1);
			}

			// Compute RGB to spectrum functions for _SampledSpectrum_
			for (int i = 0; i < nSpectralSamples; ++i) {
				Float wl0 = Lerp(Float(i) / Float(nSpectralSamples),
					sampledLambdaStart, sampledLambdaEnd);
				Float wl1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
					sampledLambdaStart, sampledLambdaEnd);
				rgbRefl2SpectWhite.c[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectWhite,
						nRGB2SpectSamples, wl0, wl1);
				rgbRefl2SpectCyan.c[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectCyan,
						nRGB2SpectSamples, wl0, wl1);
				rgbRefl2SpectMagenta.c[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectMagenta,
						nRGB2SpectSamples, wl0, wl1);
				rgbRefl2SpectYellow.c[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectYellow,
						nRGB2SpectSamples, wl0, wl1);
				rgbRefl2SpectRed.c[i] = AverageSpectrumSamples(
					RGB2SpectLambda, RGBRefl2SpectRed, nRGB2SpectSamples, wl0, wl1);
				rgbRefl2SpectGreen.c[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectGreen,
						nRGB2SpectSamples, wl0, wl1);
				rgbRefl2SpectBlue.c[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectBlue,
						nRGB2SpectSamples, wl0, wl1);

				rgbIllum2SpectWhite.c[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectWhite,
						nRGB2SpectSamples, wl0, wl1);
				rgbIllum2SpectCyan.c[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectCyan,
						nRGB2SpectSamples, wl0, wl1);
				rgbIllum2SpectMagenta.c[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectMagenta,
						nRGB2SpectSamples, wl0, wl1);
				rgbIllum2SpectYellow.c[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectYellow,
						nRGB2SpectSamples, wl0, wl1);
				rgbIllum2SpectRed.c[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectRed,
						nRGB2SpectSamples, wl0, wl1);
				rgbIllum2SpectGreen.c[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectGreen,
						nRGB2SpectSamples, wl0, wl1);
				rgbIllum2SpectBlue.c[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectBlue,
						nRGB2SpectSamples, wl0, wl1);
			}
		}
		
		void ToXYZ(Float xyz[3]) const {
			xyz[0] = xyz[1] = xyz[2] = 0.f;
			for (int i = 0; i < nSpectralSamples; ++i) {
				xyz[0] += X.c[i] * c[i];
				xyz[1] += Y.c[i] * c[i];
				xyz[2] += Z.c[i] * c[i];
			}
			Float scale = Float(sampledLambdaEnd - sampledLambdaStart) /
				Float(CIE_Y_integral * nSpectralSamples);
			xyz[0] *= scale;
			xyz[1] *= scale;
			xyz[2] *= scale;
		}
		
		Float y() const {
			Float yy = 0.f;
			for (int i = 0; i < nSpectralSamples; ++i) yy += Y.c[i] * c[i];
			return yy * Float(sampledLambdaEnd - sampledLambdaStart) /
				Float(CIE_Y_integral * nSpectralSamples);
		}
		
		void ToRGB(Float rgb[3]) const {
			Float xyz[3];
			ToXYZ(xyz);
			XYZToRGB(xyz, rgb);
		}
		
		RGBSpectrum ToRGBSpectrum() const;
		
		static SampledSpectrum FromRGB(
			const Float rgb[3], SpectrumType type = SpectrumType::Illuminant);
		static SampledSpectrum FromXYZ(
			const Float xyz[3], SpectrumType type = SpectrumType::Reflectance) {
			Float rgb[3];
			XYZToRGB(xyz, rgb);
			return FromRGB(rgb, type);
		}
		SampledSpectrum(const RGBSpectrum & r,
			SpectrumType type = SpectrumType::Reflectance);

	private:
		// SampledSpectrum Private Data
		static SampledSpectrum X, Y, Z;
		static SampledSpectrum rgbRefl2SpectWhite, rgbRefl2SpectCyan;
		static SampledSpectrum rgbRefl2SpectMagenta, rgbRefl2SpectYellow;
		static SampledSpectrum rgbRefl2SpectRed, rgbRefl2SpectGreen;
		static SampledSpectrum rgbRefl2SpectBlue;
		static SampledSpectrum rgbIllum2SpectWhite, rgbIllum2SpectCyan;
		static SampledSpectrum rgbIllum2SpectMagenta, rgbIllum2SpectYellow;
		static SampledSpectrum rgbIllum2SpectRed, rgbIllum2SpectGreen;
		static SampledSpectrum rgbIllum2SpectBlue;
	};

	/// <summary>
	/// More efficient but less accurate Representation of SPD
	/// </summary>
	/// <seealso cref="CoefficientSpectrum{3}" />
	class RGBSpectrum : public CoefficientSpectrum<3> {
		using CoefficientSpectrum<3>::c;

	public:
		// RGBSpectrum Public Methods
		RGBSpectrum(Float v = 0.f) : CoefficientSpectrum<3>(v) {}
		RGBSpectrum(const CoefficientSpectrum<3>& v) : CoefficientSpectrum<3>(v) {}
		RGBSpectrum(const RGBSpectrum& s,
			SpectrumType type = SpectrumType::Reflectance) {
			*this = s;
		}
		static RGBSpectrum FromRGB(const Float rgb[3],
			SpectrumType type = SpectrumType::Reflectance) {
			RGBSpectrum s;
			s.c[0] = rgb[0];
			s.c[1] = rgb[1];
			s.c[2] = rgb[2];
			DCHECK(!s.HasNaNs());
			return s;
		}
		void ToRGB(Float* rgb) const {
			rgb[0] = c[0];
			rgb[1] = c[1];
			rgb[2] = c[2];
		}
		const RGBSpectrum& ToRGBSpectrum() const { return *this; }
		void ToXYZ(Float xyz[3]) const { RGBToXYZ(c, xyz); }
		static RGBSpectrum FromXYZ(const Float xyz[3],
			SpectrumType type = SpectrumType::Reflectance) {
			RGBSpectrum r;
			XYZToRGB(xyz, r.c);
			return r;
		}
		Float y() const {
			const Float YWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
			return YWeight[0] * c[0] + YWeight[1] * c[1] + YWeight[2] * c[2];
		}
		static RGBSpectrum FromSampled(const Float* lambda, const Float* v, int n) {
			// Sort samples if unordered, use sorted for returned spectrum
			if (!SpectrumSamplesSorted(lambda, v, n)) {
				std::vector<Float> slambda(&lambda[0], &lambda[n]);
				std::vector<Float> sv(&v[0], &v[n]);
				SortSpectrumSamples(&slambda[0], &sv[0], n);
				return FromSampled(&slambda[0], &sv[0], n);
			}
			Float xyz[3] = { 0, 0, 0 };
			for (int i = 0; i < nCIESamples; ++i) {
				Float val = InterpolateSpectrumSamples(lambda, v, n, CIE_lambda[i]);
				xyz[0] += val * CIE_X[i];
				xyz[1] += val * CIE_Y[i];
				xyz[2] += val * CIE_Z[i];
			}
			Float scale = Float(CIE_lambda[nCIESamples - 1] - CIE_lambda[0]) /
				Float(CIE_Y_integral * nCIESamples);
			xyz[0] *= scale;
			xyz[1] *= scale;
			xyz[2] *= scale;
			return FromXYZ(xyz);
		}
	};

	// Spectrum Inline Functions
	template <int nSpectrumSamples>
	inline CoefficientSpectrum<nSpectrumSamples> Pow(
		const CoefficientSpectrum<nSpectrumSamples> & s, Float e) {
		CoefficientSpectrum<nSpectrumSamples> ret;
		for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = std::pow(s.c[i], e);
		DCHECK(!ret.HasNaNs());
		return ret;
	}

	/// <summary>
	/// Linearly interpolate between two SPDs with a parameter t.
	/// </summary>
	/// <param name="t">The t value to interpolate between.</param>
	/// <param name="s1">The spectrum 1.</param>
	/// <param name="s2">The spectrum 2.</param>
	/// <returns>The interpolated Spectrum.</returns>
	inline RGBSpectrum Lerp(Float t, const RGBSpectrum& s1, const RGBSpectrum& s2) {
		return (1 - t) * s1 + t * s2;
	}

	/// <summary>
	/// Linearly interpolate between two SPDs with a parameter t.
	/// </summary>
	/// <param name="t">The t.</param>
	/// <param name="s1">The spectrum 1.</param>
	/// <param name="s2">The spectrum 2.</param>
	/// <returns>The interpolated Spectrum.</returns>
	inline SampledSpectrum Lerp(Float t, const SampledSpectrum& s1,
		const SampledSpectrum & s2) {
		return (1 - t) * s1 + t * s2;
	}

	void ResampleLinearSpectrum(const Float * lambdaIn, const Float * vIn, int nIn,
		Float lambdaMin, Float lambdaMax, int nOut,
		Float * vOut);

}  // namespace pbrt

#endif // PBRT_CORE_SPECTRUM_H