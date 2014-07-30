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
// == idamax.h
//
// LAPACK idamax function in native C.
//

#ifndef IDAMAX_H
#define IDAMAX_H

namespace nm { namespace math {

/*  Purpose */
/*  ======= */

/*     IDAMAX finds the index of element having max. absolute value. */

/*  Further Details */
/*  =============== */

/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */

/*  ===================================================================== */

template <typename DType>
inline int idamax(size_t n, DType *dx, int incx) {

  /* Function Body */
  if (n < 1 || incx <= 0) return -1;
  if (n == 1)             return 0;

  DType dmax;
  size_t imax = 0;

  if (incx == 1) { // if incrementing by 1

    dmax = abs(dx[0]);

    for (size_t i = 1; i < n; ++i) {
      if (std::abs(dx[i]) > dmax) {
        imax = i;
        dmax = std::abs(dx[i]);
      }
    }

  } else { // if incrementing by more than 1

    dmax = std::abs(dx[0]);

    for (size_t i = 1, ix = incx; i < n; ++i, ix += incx) {
      if (std::abs(dx[ix]) > dmax) {
        imax = i;
        dmax = std::abs(dx[ix]);
      }
    }
  }
  return imax;
} /* idamax_ */

}} // end of namespace nm::math

#endif

