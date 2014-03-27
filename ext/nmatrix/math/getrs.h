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
// == getrs.h
//
// getrs function in native C++.
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

#ifndef GETRS_H
#define GETRS_H

extern "C" {
#if defined HAVE_CBLAS_H
  #include <cblas.h>
#elif defined HAVE_ATLAS_CBLAS_H
  #include <atlas/cblas.h>
#endif
}

namespace nm { namespace math {


/*
 * Solves a system of linear equations A*X = B with a general NxN matrix A using the LU factorization computed by GETRF.
 *
 * From ATLAS 3.8.0.
 */
template <typename DType>
int getrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans, const int N, const int NRHS, const DType* A,
           const int lda, const int* ipiv, DType* B, const int ldb)
{
  // enum CBLAS_DIAG Lunit, Uunit; // These aren't used. Not sure why they're declared in ATLAS' src.

  if (!N || !NRHS) return 0;

  const DType ONE = 1;

  if (Order == CblasColMajor) {
    if (Trans == CblasNoTrans) {
      nm::math::laswp<DType>(NRHS, B, ldb, 0, N, ipiv, 1);
      nm::math::trsm<DType>(Order, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, N, NRHS, ONE, A, lda, B, ldb);
      nm::math::trsm<DType>(Order, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, NRHS, ONE, A, lda, B, ldb);
    } else {
      nm::math::trsm<DType>(Order, CblasLeft, CblasUpper, Trans, CblasNonUnit, N, NRHS, ONE, A, lda, B, ldb);
      nm::math::trsm<DType>(Order, CblasLeft, CblasLower, Trans, CblasUnit, N, NRHS, ONE, A, lda, B, ldb);
      nm::math::laswp<DType>(NRHS, B, ldb, 0, N, ipiv, -1);
    }
  } else {
    if (Trans == CblasNoTrans) {
      nm::math::trsm<DType>(Order, CblasRight, CblasLower, CblasTrans, CblasNonUnit, NRHS, N, ONE, A, lda, B, ldb);
      nm::math::trsm<DType>(Order, CblasRight, CblasUpper, CblasTrans, CblasUnit, NRHS, N, ONE, A, lda, B, ldb);
      nm::math::laswp<DType>(NRHS, B, ldb, 0, N, ipiv, -1);
    } else {
      nm::math::laswp<DType>(NRHS, B, ldb, 0, N, ipiv, 1);
      nm::math::trsm<DType>(Order, CblasRight, CblasUpper, CblasNoTrans, CblasUnit, NRHS, N, ONE, A, lda, B, ldb);
      nm::math::trsm<DType>(Order, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, NRHS, N, ONE, A, lda, B, ldb);
    }
  }
  return 0;
}


/*
* Function signature conversion for calling LAPACK's getrs functions as directly as possible.
*
* For documentation: http://www.netlib.org/lapack/double/dgetrs.f
*
* This function should normally go in math.cpp, but we need it to be available to nmatrix.cpp.
*/
template <typename DType>
inline int clapack_getrs(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE trans, const int n, const int nrhs,
                         const void* a, const int lda, const int* ipiv, void* b, const int ldb) {
  return getrs<DType>(order, trans, n, nrhs, reinterpret_cast<const DType*>(a), lda, ipiv, reinterpret_cast<DType*>(b), ldb);
}


} } // end nm::math

#endif // GETRS_H
