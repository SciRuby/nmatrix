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
// == swap.h
//
// BLAS level 2 swap function in native C++.
//

#ifndef SWAP_H
#define SWAP_H

namespace nm { namespace math {
/*
template <typename DType>
inline void swap(int n, DType *dx, int incx, DType *dy, int incy) {

  if (n <= 0) return;

  // For negative increments, start at the end of the array.
  int ix = incx < 0 ? (-n+1)*incx : 0,
      iy = incy < 0 ? (-n+1)*incy : 0;

  if (incx < 0) ix = (-n + 1) * incx;
  if (incy < 0) iy = (-n + 1) * incy;

  for (size_t i = 0; i < n; ++i, ix += incx, iy += incy) {
    DType dtemp = dx[ix];
    dx[ix]      = dy[iy];
    dy[iy]      = dtemp;
  }
  return;
} /* dswap */

// This is the old BLAS version of this function. ATLAS has an optimized version, but
// it's going to be tough to translate.
template <typename DType>
static void swap(const int N, DType* X, const int incX, DType* Y, const int incY) {
  if (N > 0) {
    int ix = 0, iy = 0;
    for (int i = 0; i < N; ++i) {
      DType temp = X[i];
      X[i]       = Y[i];
      Y[i]       = temp;

      ix += incX;
      iy += incY;
    }
  }
}

}} // end nm::math

#endif