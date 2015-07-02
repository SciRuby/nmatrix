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
// == getrf.h
//

#ifndef GETRF_LAPACK_H
#define GETRF_LAPACK_H

#include "math/getrf.h"

namespace nm { namespace math { namespace lapack {

template <typename DType>
inline int getrf(const enum CBLAS_ORDER order, const int m, const int n, DType* a, const int lda, int* ipiv) {
  //unless overridden below, use the internal implementation of getrf
  //from ext/nmatrix/math/getrf.h
  //THIS IS PROBLEMATIC because the internal implementation uses the ATLAS
  //interpretation of getrf which behaves differently for row-major matrices
  //(different in terms of whether upper or lower matrix has the unit diagonals).
  return nm::math::getrf<DType>(order, m, n, a, lda, ipiv);
}

template <>
inline int getrf(const enum CBLAS_ORDER order, const int m, const int n, float* a, const int lda, int* ipiv) {
  return LAPACKE_sgetrf(order, m, n, a, lda, ipiv);
}

template <>
inline int getrf(const enum CBLAS_ORDER order, const int m, const int n, double* a, const int lda, int* ipiv) {
  return LAPACKE_dgetrf(order, m, n, a, lda, ipiv);
}

template <>
inline int getrf(const enum CBLAS_ORDER order, const int m, const int n, Complex64* a, const int lda, int* ipiv) {
  return LAPACKE_cgetrf(order, m, n, a, lda, ipiv);
}

template <>
inline int getrf(const enum CBLAS_ORDER order, const int m, const int n, Complex128* a, const int lda, int* ipiv) {
  return LAPACKE_zgetrf(order, m, n, a, lda, ipiv);
}

/*
* Function signature conversion for calling LAPACK's getrf functions as directly as possible.
*
* For documentation: http://www.netlib.org/lapack/double/dgetrf.f
*
* This function should normally go in math.cpp, but we need it to be available to nmatrix.cpp.
*/
template <typename DType>
inline int clapack_getrf(const enum CBLAS_ORDER order, const int m, const int n, void* a, const int lda, int* ipiv) {
  return getrf<DType>(order, m, n, static_cast<DType*>(a), lda, ipiv);
}

} } } // end nm::math::lapack

#endif // GETRF_LAPACK_H