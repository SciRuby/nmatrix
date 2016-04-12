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
// SciRuby is Copyright (c) 2010 - present, Ruby Science Foundation
// NMatrix is Copyright (c) 2012 - present, John Woods and the Ruby Science Foundation
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
// == imax.h
//
// BLAS level 1 function imax.
//

#ifndef IMAX_H
#define IMAX_H

#include "math/magnitude.h"

namespace nm { namespace math {


template<typename DType>
inline int imax(const int n, const DType *x, const int incx) {

  if (n < 1 || incx <= 0) {
    return -1;
  }
  if (n == 1) {
    return 0;
  }

  typename MagnitudeDType<DType>::type dmax;
  int imax = 0;

  if (incx == 1) { // if incrementing by 1

    dmax = magnitude(x[0]);

    for (int i = 1; i < n; ++i) {
      if (magnitude(x[i]) > dmax) {
        imax = i;
        dmax = magnitude(x[i]);
      }
    }

  } else { // if incrementing by more than 1

    dmax = magnitude(x[0]);

    for (int i = 1, ix = incx; i < n; ++i, ix += incx) {
      if (magnitude(x[ix]) > dmax) {
        imax = i;
        dmax = magnitude(x[ix]);
      }
    }
  }
  return imax;
}

template<typename DType>
inline int cblas_imax(const int n, const void* x, const int incx) {
  return imax<DType>(n, reinterpret_cast<const DType*>(x), incx);
}

}} // end of namespace nm::math

#endif /* IMAX_H */
