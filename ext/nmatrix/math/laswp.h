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
// == laswp.h
//
// laswp function in native C++.
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

#ifndef LASWP_H
#define LASWP_H

namespace nm { namespace math {


/*
 * ATLAS function which performs row interchanges on a general rectangular matrix. Modeled after the LAPACK LASWP function.
 *
 * This version is templated for use by template <> getrf().
 */
template <typename DType>
inline void laswp(const int N, DType* A, const int lda, const int K1, const int K2, const int *piv, const int inci) {
  //const int n = K2 - K1; // not sure why this is declared. commented it out because it's unused.

  int nb = N >> 5;

  const int mr = N - (nb<<5);
  const int incA = lda << 5;

  if (K2 < K1) return;

  int i1, i2;
  if (inci < 0) {
    piv -= (K2-1) * inci;
    i1 = K2 - 1;
    i2 = K1;
  } else {
    piv += K1 * inci;
    i1 = K1;
    i2 = K2-1;
  }

  if (nb) {

    do {
      const int* ipiv = piv;
      int i           = i1;
      int KeepOn;

      do {
        int ip = *ipiv; ipiv += inci;

        if (ip != i) {
          DType *a0 = &(A[i]),
                *a1 = &(A[ip]);

          for (register int h = 32; h; h--) {
            DType r   = *a0;
            *a0       = *a1;
            *a1       = r;

            a0 += lda;
            a1 += lda;
          }

        }
        if (inci > 0) KeepOn = (++i <= i2);
        else          KeepOn = (--i >= i2);

      } while (KeepOn);
      A += incA;
    } while (--nb);
  }

  if (mr) {
    const int* ipiv = piv;
    int i           = i1;
    int KeepOn;

    do {
      int ip = *ipiv; ipiv += inci;
      if (ip != i) {
        DType *a0 = &(A[i]),
              *a1 = &(A[ip]);

        for (register int h = mr; h; h--) {
          DType r   = *a0;
          *a0       = *a1;
          *a1       = r;

          a0 += lda;
          a1 += lda;
        }
      }

      if (inci > 0) KeepOn = (++i <= i2);
      else          KeepOn = (--i >= i2);

    } while (KeepOn);
  }
}


/*
* Function signature conversion for calling LAPACK's laswp functions as directly as possible.
*
* For documentation: http://www.netlib.org/lapack/double/dlaswp.f
*
* This function should normally go in math.cpp, but we need it to be available to nmatrix.cpp.
*/
template <typename DType>
inline void clapack_laswp(const int n, void* a, const int lda, const int k1, const int k2, const int* ipiv, const int incx) {
  laswp<DType>(n, reinterpret_cast<DType*>(a), lda, k1, k2, ipiv, incx);
}

} }  // namespace nm::math
#endif // LASWP_H