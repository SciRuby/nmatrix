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
// == scal.h
//
// LAPACK scal function in native C.
//

#ifndef SCAL_H
#define SCAL_H

namespace nm { namespace math {

/*  Purpose */
/*  ======= */

/*     DSCAL scales a vector by a constant. */
/*     uses unrolled loops for increment equal to one. */

/*  Further Details */
/*  =============== */

/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */

/*  ===================================================================== */

template <typename DType>
inline void scal(const int n, const DType da, DType* dx,	const int incx) {

  // This used to have unrolled loops, like dswap. They were in the way.

  if (n <= 0 || incx <= 0) return;

  for (int i = 0; incx < 0 ? i > n*incx : i < n*incx; i += incx) {
    dx[i] = da * dx[i];
  }
} /* scal */


/*
 * Function signature conversion for LAPACK's scal function.
 */
template <typename DType>
inline void clapack_scal(const int n, const void* da, void* dx, const int incx) {
  // FIXME: See if we can call the clapack version instead of our C++ version.
  scal<DType>(n, *reinterpret_cast<const DType*>(da), reinterpret_cast<DType*>(dx), incx);
}

}} // end of nm::math

#endif
