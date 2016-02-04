/////////////////////////////////////////////////////////////////////
// = NMatrix
//
// A linear algebra library for scientific computation in Ruby.
// NMatrix is part of SciRuby.
//
// NMatrix was originally inspired by and derived from NArray, by
// Masahiro Tanaka: http://narray.rubyforge.org
//
// == Copyright Information
//
// SciRuby is Copyright (c) 2010 - 2014, Ruby Science Foundation
// NMatrix is Copyright (c) 2012 - 2014, John Woods and the Ruby Science Foundation
//
// Please see LICENSE.txt for additional copyright notices.
//
// == Contributing
//
// By contributing source code to SciRuby, you agree to be bound by
// our Contributor Agreement:
//
// * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
//
// == nrm2.h
//
// CBLAS nrm2 function
//

/*
 *             Automatically Tuned Linear Algebra Software v3.8.4
 *                    (C) Copyright 1999 R. Clint Whaley
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef NRM2_H
# define NRM2_H

#include "math/long_dtype.h"


namespace nm { namespace math {

/*
 * Level 1 BLAS routine which returns the 2-norm of an n-vector x.
 #
 * Based on input types, these are the valid return types:
 *    int -> int
 *    float -> float or double
 *    double -> double
 *    complex64 -> float or double
 *    complex128 -> double
 */
template <typename ReturnDType, typename DType>
ReturnDType nrm2(const int N, const DType* X, const int incX) {
  const DType ONE = 1, ZERO = 0;
  typename LongDType<DType>::type scale = 0, ssq = 1, absxi, temp;


  if ((N < 1) || (incX < 1))    return ZERO;
  else if (N == 1)              return std::abs(X[0]);

  for (int i = 0; i < N; ++i) {
    absxi = std::abs(X[i*incX]);
    if (scale < absxi) {
      temp  = scale / absxi;
      scale = absxi;
      ssq   = ONE + ssq * (temp * temp);
    } else {
      temp = absxi / scale;
      ssq += temp * temp;
    }
  }

  return scale * std::sqrt( ssq );
}


template <typename FloatDType>
static inline void nrm2_complex_helper(const FloatDType& xr, const FloatDType& xi, double& scale, double& ssq) {
  double absx = std::abs(xr);
  if (scale < absx) {
    double temp  = scale / absx;
    scale = absx;
    ssq   = 1.0 + ssq * (temp * temp);
  } else {
    double temp = absx / scale;
    ssq += temp * temp;
  }

  absx = std::abs(xi);
  if (scale < absx) {
    double temp  = scale / absx;
    scale = absx;
    ssq   = 1.0 + ssq * (temp * temp);
  } else {
    double temp = absx / scale;
    ssq += temp * temp;
  }
}

template <>
float nrm2(const int N, const Complex64* X, const int incX) {
  double scale = 0, ssq = 1;

  if ((N < 1) || (incX < 1))    return 0.0;

  for (int i = 0; i < N; ++i) {
    nrm2_complex_helper<float>(X[i*incX].r, X[i*incX].i, scale, ssq);
  }

  return scale * std::sqrt( ssq );
}

template <>
double nrm2(const int N, const Complex128* X, const int incX) {
  double scale = 0, ssq = 1;

  if ((N < 1) || (incX < 1))    return 0.0;

  for (int i = 0; i < N; ++i) {
    nrm2_complex_helper<double>(X[i*incX].r, X[i*incX].i, scale, ssq);
  }

  return scale * std::sqrt( ssq );
}

template <typename ReturnDType, typename DType>
inline void cblas_nrm2(const int N, const void* X, const int incX, void* result) {
  *reinterpret_cast<ReturnDType*>( result ) = nrm2<ReturnDType, DType>( N, reinterpret_cast<const DType*>(X), incX );
}



}} // end of namespace nm::math

#endif // NRM2_H
