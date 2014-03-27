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
// == rot.h
//
// BLAS rot function in native C++.
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

#ifndef ROT_H
# define ROT_H

namespace nm { namespace math {

// FIXME: This is not working properly for rational numbers. Do we need some kind of symbolic
// type to handle square roots?


// TODO: Test this to see if it works properly on complex. ATLAS has a separate algorithm for complex, which looks like
// TODO: it may actually be the same one.
//
// This function is called ATL_rot in ATLAS 3.8.4.
template <typename DType>
inline void rot_helper(const int N, DType* X, const int incX, DType* Y, const int incY, const DType c, const DType s) {
  if (c != 1 || s != 0) {
    if (incX == 1 && incY == 1) {
      for (int i = 0; i != N; ++i) {
        DType tmp = X[i] * c + Y[i] * s;
        Y[i]      = Y[i] * c - X[i] * s;
        X[i]      = tmp;
      }
    } else {
      for (int i = N; i > 0; --i, Y += incY, X += incX) {
        DType tmp = *X * c + *Y * s;
        *Y  = *Y * c - *X * s;
        *X  = tmp;
      }
    }
  }
}


/* Applies a plane rotation. From ATLAS 3.8.4. */
template <typename DType, typename CSDType>
inline void rot(const int N, DType* X, const int incX, DType* Y, const int incY, const CSDType c, const CSDType s) {
  int incx = incX, incy = incY;
  DType *x = X, *y = Y;

  if (N > 0) {
    if (incX < 0) {
      if (incY < 0) { incx = -incx; incy = -incy; }
      else x += -incX * (N-1);
    } else if (incY < 0) {
      incy = -incy;
      incx = -incx;
      x += (N-1) * incX;
    }
    rot_helper<DType>(N, x, incx, y, incy, c, s);
  }
}

template <>
inline void rot(const int N, float* X, const int incX, float* Y, const int incY, const float c, const float s) {
  cblas_srot(N, X, incX, Y, incY, (float)c, (float)s);
}

template <>
inline void rot(const int N, double* X, const int incX, double* Y, const int incY, const double c, const double s) {
  cblas_drot(N, X, incX, Y, incY, c, s);
}

template <>
inline void rot(const int N, Complex64* X, const int incX, Complex64* Y, const int incY, const float c, const float s) {
  cblas_csrot(N, X, incX, Y, incY, c, s);
}

template <>
inline void rot(const int N, Complex128* X, const int incX, Complex128* Y, const int incY, const double c, const double s) {
  cblas_zdrot(N, X, incX, Y, incY, c, s);
}


template <typename DType, typename CSDType>
inline void cblas_rot(const int N, void* X, const int incX, void* Y, const int incY, const void* c, const void* s) {
  rot<DType,CSDType>(N, reinterpret_cast<DType*>(X), incX, reinterpret_cast<DType*>(Y), incY,
                       *reinterpret_cast<const CSDType*>(c), *reinterpret_cast<const CSDType*>(s));
}


} } //nm::math

#endif // ROT_H