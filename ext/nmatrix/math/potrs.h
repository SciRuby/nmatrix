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

#ifndef POTRS_H
#define POTRS_H

extern "C" {
#if defined HAVE_CBLAS_H
  #include <cblas.h>
#elif defined HAVE_ATLAS_CBLAS_H
  #include <atlas/cblas.h>
#endif
}

namespace nm { namespace math {

/*
 * Solves a system of linear equations A*X = B with a symmetric positive definite matrix A using the Cholesky factorization computed by POTRF.
 *
 * From ATLAS 3.8.0.
 */
template <typename DType, bool is_complex>
int potrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const int N, const int NRHS, const DType* A,
           const int lda, DType* B, const int ldb)
{
  // enum CBLAS_DIAG Lunit, Uunit; // These aren't used. Not sure why they're declared in ATLAS' src.

  CBLAS_TRANSPOSE MyTrans = is_complex ? CblasConjTrans : CblasTrans;

  if (!N || !NRHS) return 0;

  const DType ONE = 1;

  if (Order == CblasColMajor) {
    if (Uplo == CblasUpper) {
      nm::math::trsm<DType>(Order, CblasLeft, CblasUpper, MyTrans, CblasNonUnit, N, NRHS, ONE, A, lda, B, ldb);
      nm::math::trsm<DType>(Order, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, NRHS, ONE, A, lda, B, ldb);
    } else {
      nm::math::trsm<DType>(Order, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, N, NRHS, ONE, A, lda, B, ldb);
      nm::math::trsm<DType>(Order, CblasLeft, CblasLower, MyTrans, CblasNonUnit, N, NRHS, ONE, A, lda, B, ldb);
    }
  } else {
    // There's some kind of scaling operation that normally happens here in ATLAS. Not sure what it does, so we'll only
    // worry if something breaks. It probably has to do with their non-templated code and doesn't apply to us.

    if (Uplo == CblasUpper) {
      nm::math::trsm<DType>(Order, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, NRHS, N, ONE, A, lda, B, ldb);
      nm::math::trsm<DType>(Order, CblasRight, CblasUpper, MyTrans, CblasNonUnit, NRHS, N, ONE, A, lda, B, ldb);
    } else {
      nm::math::trsm<DType>(Order, CblasRight, CblasLower, MyTrans, CblasNonUnit, NRHS, N, ONE, A, lda, B, ldb);
      nm::math::trsm<DType>(Order, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, NRHS, N, ONE, A, lda, B, ldb);
    }
  }
  return 0;
}


/*
* Function signature conversion for calling LAPACK's potrs functions as directly as possible.
*
* For documentation: http://www.netlib.org/lapack/double/dpotrs.f
*
* This function should normally go in math.cpp, but we need it to be available to nmatrix.cpp.
*/
template <typename DType, bool is_complex>
inline int clapack_potrs(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, const int nrhs,
                         const void* a, const int lda, void* b, const int ldb) {
  return potrs<DType,is_complex>(order, uplo, n, nrhs, reinterpret_cast<const DType*>(a), lda, reinterpret_cast<DType*>(b), ldb);
}


} } // end nm::math

#endif // POTRS_H
