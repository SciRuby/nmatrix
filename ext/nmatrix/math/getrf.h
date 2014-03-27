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
// getrf function in native C++.
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

#ifndef GETRF_H
#define GETRF_H

namespace nm { namespace math {

/* Numeric inverse -- usually just 1 / f, but a little more complicated for complex. */
template <typename DType>
inline DType numeric_inverse(const DType& n) {
  return n.inverse();
}
template <> inline float numeric_inverse(const float& n) { return 1 / n; }
template <> inline double numeric_inverse(const double& n) { return 1 / n; }



/*
 * Templated version of row-order and column-order getrf, derived from ATL_getrfR.c (from ATLAS 3.8.0).
 *
 * 1. Row-major factorization of form
 *   A = L * U * P
 * where P is a column-permutation matrix, L is lower triangular (lower
 * trapazoidal if M > N), and U is upper triangular with unit diagonals (upper
 * trapazoidal if M < N).  This is the recursive Level 3 BLAS version.
 *
 * 2. Column-major factorization of form
 *   A = P * L * U
 * where P is a row-permutation matrix, L is lower triangular with unit diagonal
 * elements (lower trapazoidal if M > N), and U is upper triangular (upper
 * trapazoidal if M < N).  This is the recursive Level 3 BLAS version.
 *
 * Template argument determines whether 1 or 2 is utilized.
 */
template <bool RowMajor, typename DType>
inline int getrf_nothrow(const int M, const int N, DType* A, const int lda, int* ipiv) {
  const int MN = std::min(M, N);
  int ierr = 0;

  // Symbols used by ATLAS in the several versions of this function:
  // Row   Col      Us
  // Nup   Nleft    N_ul
  // Ndown Nright   N_dr
  // We're going to use N_ul, N_dr

  DType neg_one = -1, one = 1;

  if (MN > 1) {
    int N_ul = MN >> 1;

    // FIXME: Figure out how ATLAS #defines NB
#ifdef NB
    if (N_ul > NB) N_ul = ATL_MulByNB(ATL_DivByNB(N_ul));
#endif

    int N_dr = M - N_ul;

    int i = RowMajor ? getrf_nothrow<true,DType>(N_ul, N, A, lda, ipiv) : getrf_nothrow<false,DType>(M, N_ul, A, lda, ipiv);

    if (i) if (!ierr) ierr = i;

    DType *Ar, *Ac, *An;
    if (RowMajor) {
      Ar = &(A[N_ul * lda]),
      Ac = &(A[N_ul]);
      An = &(Ar[N_ul]);

      nm::math::laswp<DType>(N_dr, Ar, lda, 0, N_ul, ipiv, 1);

      nm::math::trsm<DType>(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasUnit, N_dr, N_ul, one, A, lda, Ar, lda);
      nm::math::gemm<DType>(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_dr, N-N_ul, N_ul, &neg_one, Ar, lda, Ac, lda, &one, An, lda);

      i = getrf_nothrow<true,DType>(N_dr, N-N_ul, An, lda, ipiv+N_ul);
    } else {
      Ar = NULL;
      Ac = &(A[N_ul * lda]);
      An = &(Ac[N_ul]);

      nm::math::laswp<DType>(N_dr, Ac, lda, 0, N_ul, ipiv, 1);

      nm::math::trsm<DType>(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, N_ul, N_dr, one, A, lda, Ac, lda);
      nm::math::gemm<DType>(CblasColMajor, CblasNoTrans, CblasNoTrans, M-N_ul, N_dr, N_ul, &neg_one, An, lda, Ac, lda, &one, An, lda);

      i = getrf_nothrow<false,DType>(M-N_ul, N_dr, An, lda, ipiv+N_ul);
    }

    if (i) if (!ierr) ierr = N_ul + i;

    for (i = N_ul; i != MN; i++) {
      ipiv[i] += N_ul;
    }

    nm::math::laswp<DType>(N_ul, A, lda, N_ul, MN, ipiv, 1);  /* apply pivots */

  } else if (MN == 1) { // there's another case for the colmajor version, but i don't know that it's that critical. Calls ATLAS LU2, who knows what that does.

    int i = *ipiv = nm::math::idamax<DType>(N, A, 1); // cblas_iamax(N, A, 1);

    DType tmp = A[i];
    if (tmp != 0) {

      nm::math::scal<DType>((RowMajor ? N : M), nm::math::numeric_inverse(tmp), A, 1);
      A[i] = *A;
      *A   = tmp;

    } else ierr = 1;

  }
  return(ierr);
}


/*
 * From ATLAS 3.8.0:
 *
 * Computes one of two LU factorizations based on the setting of the Order
 * parameter, as follows:
 * ----------------------------------------------------------------------------
 *                       Order == CblasColMajor
 * Column-major factorization of form
 *   A = P * L * U
 * where P is a row-permutation matrix, L is lower triangular with unit
 * diagonal elements (lower trapazoidal if M > N), and U is upper triangular
 * (upper trapazoidal if M < N).
 *
 * ----------------------------------------------------------------------------
 *                       Order == CblasRowMajor
 * Row-major factorization of form
 *   A = P * L * U
 * where P is a column-permutation matrix, L is lower triangular (lower
 * trapazoidal if M > N), and U is upper triangular with unit diagonals (upper
 * trapazoidal if M < N).
 *
 * ============================================================================
 * Let IERR be the return value of the function:
 *    If IERR == 0, successful exit.
 *    If (IERR < 0) the -IERR argument had an illegal value
 *    If (IERR > 0 && Order == CblasColMajor)
 *       U(i-1,i-1) is exactly zero.  The factorization has been completed,
 *       but the factor U is exactly singular, and division by zero will
 *       occur if it is used to solve a system of equations.
 *    If (IERR > 0 && Order == CblasRowMajor)
 *       L(i-1,i-1) is exactly zero.  The factorization has been completed,
 *       but the factor L is exactly singular, and division by zero will
 *       occur if it is used to solve a system of equations.
 */
template <typename DType>
inline int getrf(const enum CBLAS_ORDER Order, const int M, const int N, DType* A, int lda, int* ipiv) {
  if (Order == CblasRowMajor) {
    if (lda < std::max(1,N)) {
      rb_raise(rb_eArgError, "GETRF: lda must be >= MAX(N,1): lda=%d N=%d", lda, N);
      return -6;
    }

    return getrf_nothrow<true,DType>(M, N, A, lda, ipiv);
  } else {
    if (lda < std::max(1,M)) {
      rb_raise(rb_eArgError, "GETRF: lda must be >= MAX(M,1): lda=%d M=%d", lda, M);
      return -6;
    }

    return getrf_nothrow<false,DType>(M, N, A, lda, ipiv);
    //rb_raise(rb_eNotImpError, "column major getrf not implemented");
  }
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
  return getrf<DType>(order, m, n, reinterpret_cast<DType*>(a), lda, ipiv);
}


} } // end nm::math

#endif