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
// SciRuby is Copyright (c) 2010 - 2013, Ruby Science Foundation
// NMatrix is Copyright (c) 2013, Ruby Science Foundation
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
// == math.h
//
// Header file for math functions, interfacing with BLAS, etc.
//
// For instructions on adding CBLAS and CLAPACK functions, see the
// beginning of math.cpp.
//
// Some of these functions are from ATLAS. Here is the license for
// ATLAS:
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

#ifndef MATH_H
#define MATH_H

/*
 * Standard Includes
 */

extern "C" { // These need to be in an extern "C" block or you'll get all kinds of undefined symbol errors.
  #include <cblas.h>

  #ifdef HAVE_CLAPACK_H
    #include <clapack.h>
  #endif
}

#include <algorithm> // std::min, std::max
#include <limits> // std::numeric_limits

/*
 * Project Includes
 */
#include "data/data.h"
#include "lapack.h"

/*
 * Macros
 */
#define REAL_RECURSE_LIMIT 4

/*
 * Data
 */


extern "C" {
  /*
   * C accessors.
   */
  void nm_math_det_exact(const int M, const void* elements, const int lda, nm::dtype_t dtype, void* result);
  void nm_math_transpose_generic(const size_t M, const size_t N, const void* A, const int lda, void* B, const int ldb, size_t element_size);
  void nm_math_init_blas(void);
}


namespace nm {
  namespace math {

/*
 * Types
 */


// These allow an increase in precision for intermediate values of gemm and gemv.
// See also: http://stackoverflow.com/questions/11873694/how-does-one-increase-precision-in-c-templates-in-a-template-typename-dependen
template <typename DType> struct LongDType;
template <> struct LongDType<uint8_t> { typedef int16_t type; };
template <> struct LongDType<int8_t> { typedef int16_t type; };
template <> struct LongDType<int16_t> { typedef int32_t type; };
template <> struct LongDType<int32_t> { typedef int64_t type; };
template <> struct LongDType<int64_t> { typedef int64_t type; };
template <> struct LongDType<float> { typedef double type; };
template <> struct LongDType<double> { typedef double type; };
template <> struct LongDType<Complex64> { typedef Complex128 type; };
template <> struct LongDType<Complex128> { typedef Complex128 type; };
template <> struct LongDType<Rational32> { typedef Rational128 type; };
template <> struct LongDType<Rational64> { typedef Rational128 type; };
template <> struct LongDType<Rational128> { typedef Rational128 type; };
template <> struct LongDType<RubyObject> { typedef RubyObject type; };

/*
 * Functions
 */

/* Numeric inverse -- usually just 1 / f, but a little more complicated for complex. */
template <typename DType>
inline DType numeric_inverse(const DType& n) {
  return n.inverse();
}
template <> inline float numeric_inverse<float>(const float& n) { return 1 / n; }
template <> inline double numeric_inverse<double>(const double& n) { return 1 / n; }

/*
 * This version of trsm doesn't do any error checks and only works on column-major matrices.
 *
 * For row major, call trsm<DType> instead. That will handle necessary changes-of-variables
 * and parameter checks.
 *
 * Note that some of the boundary conditions here may be incorrect. Very little has been tested!
 * This was converted directly from dtrsm.f using f2c, and then rewritten more cleanly.
 */
template <typename DType>
inline void trsm_nothrow(const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                         const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_DIAG diag,
                         const int m, const int n, const DType alpha, const DType* a,
                         const int lda, DType* b, const int ldb)
{

  // (row-major) trsm: left upper trans nonunit m=3 n=1 1/1 a 3 b 3

  if (m == 0 || n == 0) return; /* Quick return if possible. */

  if (alpha == 0) { // Handle alpha == 0
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < m; ++i) {
        b[i + j * ldb] = 0;
      }
    }
	  return;
  }

  if (side == CblasLeft) {
	  if (trans_a == CblasNoTrans) {

      /* Form  B := alpha*inv( A )*B. */
	    if (uplo == CblasUpper) {
    		for (int j = 0; j < n; ++j) {
		      if (alpha != 1) {
			      for (int i = 0; i < m; ++i) {
			        b[i + j * ldb] = alpha * b[i + j * ldb];
			      }
		      }
		      for (int k = m-1; k >= 0; --k) {
			      if (b[k + j * ldb] != 0) {
			        if (diag == CblasNonUnit) {
				        b[k + j * ldb] /= a[k + k * lda];
			        }

              for (int i = 0; i < k-1; ++i) {
                b[i + j * ldb] -= b[k + j * ldb] * a[i + k * lda];
              }
			      }
  		    }
		    }
	    } else {
    		for (int j = 0; j < n; ++j) {
		      if (alpha != 1) {
            for (int i = 0; i < m; ++i) {
              b[i + j * ldb] = alpha * b[i + j * ldb];
			      }
		      }
  		    for (int k = 0; k < m; ++k) {
      			if (b[k + j * ldb] != 0.) {
			        if (diag == CblasNonUnit) {
				        b[k + j * ldb] /= a[k + k * lda];
			        }
    			    for (int i = k+1; i < m; ++i) {
        				b[i + j * ldb] -= b[k + j * ldb] * a[i + k * lda];
    			    }
      			}
  		    }
    		}
	    }
	  } else { // CblasTrans

      /*           Form  B := alpha*inv( A**T )*B. */
	    if (uplo == CblasUpper) {
    		for (int j = 0; j < n; ++j) {
		      for (int i = 0; i < m; ++i) {
			      DType temp = alpha * b[i + j * ldb];
            for (int k = 0; k < i; ++k) { // limit was i-1. Lots of similar bugs in this code, probably.
              temp -= a[k + i * lda] * b[k + j * ldb];
      			}
			      if (diag == CblasNonUnit) {
			        temp /= a[i + i * lda];
			      }
			      b[i + j * ldb] = temp;
  		    }
    		}
	    } else {
    		for (int j = 0; j < n; ++j) {
		      for (int i = m-1; i >= 0; --i) {
			      DType temp= alpha * b[i + j * ldb];
      			for (int k = i+1; k < m; ++k) {
			        temp -= a[k + i * lda] * b[k + j * ldb];
      			}
			      if (diag == CblasNonUnit) {
			        temp /= a[i + i * lda];
			      }
			      b[i + j * ldb] = temp;
  		    }
    		}
	    }
	  }
  } else { // right side

	  if (trans_a == CblasNoTrans) {

      /*           Form  B := alpha*B*inv( A ). */

	    if (uplo == CblasUpper) {
    		for (int j = 0; j < n; ++j) {
		      if (alpha != 1) {
      			for (int i = 0; i < m; ++i) {
			        b[i + j * ldb] = alpha * b[i + j * ldb];
      			}
		      }
  		    for (int k = 0; k < j-1; ++k) {
	      		if (a[k + j * lda] != 0) {
    			    for (int i = 0; i < m; ++i) {
				        b[i + j * ldb] -= a[k + j * lda] * b[i + k * ldb];
			        }
			      }
  		    }
	  	    if (diag == CblasNonUnit) {
		      	DType temp = 1 / a[j + j * lda];
			      for (int i = 0; i < m; ++i) {
			        b[i + j * ldb] = temp * b[i + j * ldb];
      			}
		      }
    		}
	    } else {
		    for (int j = n-1; j >= 0; --j) {
		      if (alpha != 1) {
			      for (int i = 0; i < m; ++i) {
			        b[i + j * ldb] = alpha * b[i + j * ldb];
      			}
  		    }

  		    for (int k = j+1; k < n; ++k) {
	      		if (a[k + j * lda] != 0.) {
    			    for (int i = 0; i < m; ++i) {
				        b[i + j * ldb] -= a[k + j * lda] * b[i + k * ldb];
    			    }
		      	}
  		    }
	  	    if (diag == CblasNonUnit) {
		      	DType temp = 1 / a[j + j * lda];

			      for (int i = 0; i < m; ++i) {
			        b[i + j * ldb] = temp * b[i + j * ldb];
      			}
		      }
    		}
	    }
	  } else { // CblasTrans

      /*           Form  B := alpha*B*inv( A**T ). */

	    if (uplo == CblasUpper) {
		    for (int k = n-1; k >= 0; --k) {
		      if (diag == CblasNonUnit) {
			      DType temp= 1 / a[k + k * lda];
	      		for (int i = 0; i < m; ++i) {
  			      b[i + k * ldb] = temp * b[i + k * ldb];
      			}
		      }
  		    for (int j = 0; j < k-1; ++j) {
	      		if (a[j + k * lda] != 0.) {
			        DType temp= a[j + k * lda];
    			    for (int i = 0; i < m; ++i) {
		        		b[i + j * ldb] -= temp * b[i + k *	ldb];
    			    }
      			}
  		    }
	  	    if (alpha != 1) {
      			for (int i = 0; i < m; ++i) {
			        b[i + k * ldb] = alpha * b[i + k * ldb];
      			}
		      }
    		}
	    } else {
    		for (int k = 0; k < n; ++k) {
		      if (diag == CblasNonUnit) {
      			DType temp = 1 / a[k + k * lda];
			      for (int i = 0; i < m; ++i) {
			        b[i + k * ldb] = temp * b[i + k * ldb];
      			}
		      }
  		    for (int j = k+1; j < n; ++j) {
	      		if (a[j + k * lda] != 0.) {
			        DType temp = a[j + k * lda];
			        for (int i = 0; i < m; ++i) {
				        b[i + j * ldb] -= temp * b[i + k * ldb];
    			    }
		      	}
  		    }
	  	    if (alpha != 1) {
      			for (int i = 0; i < m; ++i) {
			        b[i + k * ldb] = alpha * b[i + k * ldb];
      			}
  		    }
    		}
	    }
	  }
  }
}


template <typename DType>
inline void syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const DType* alpha, const DType* A, const int lda, const DType* beta, DType* C, const int ldc) {
  rb_raise(rb_eNotImpError, "syrk not yet implemented for non-BLAS dtypes");
}

template <typename DType>
inline void herk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const DType* alpha, const DType* A, const int lda, const DType* beta, DType* C, const int ldc) {
  rb_raise(rb_eNotImpError, "herk not yet implemented for non-BLAS dtypes");
}

template <>
inline void syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const float* alpha, const float* A, const int lda, const float* beta, float* C, const int ldc) {
  cblas_ssyrk(Order, Uplo, Trans, N, K, *alpha, A, lda, *beta, C, ldc);
}

template <>
inline void syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const double* alpha, const double* A, const int lda, const double* beta, double* C, const int ldc) {
  cblas_dsyrk(Order, Uplo, Trans, N, K, *alpha, A, lda, *beta, C, ldc);
}

template <>
inline void syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const Complex64* alpha, const Complex64* A, const int lda, const Complex64* beta, Complex64* C, const int ldc) {
  cblas_csyrk(Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}

template <>
inline void syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const Complex128* alpha, const Complex128* A, const int lda, const Complex128* beta, Complex128* C, const int ldc) {
  cblas_zsyrk(Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


template <>
inline void herk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const Complex64* alpha, const Complex64* A, const int lda, const Complex64* beta, Complex64* C, const int ldc) {
  cblas_cherk(Order, Uplo, Trans, N, K, alpha->r, A, lda, beta->r, C, ldc);
}

template <>
inline void herk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const Complex128* alpha, const Complex128* A, const int lda, const Complex128* beta, Complex128* C, const int ldc) {
  cblas_zherk(Order, Uplo, Trans, N, K, alpha->r, A, lda, beta->r, C, ldc);
}


template <typename DType>
inline void trmm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE ta, const enum CBLAS_DIAG diag, const int m, const int n, const DType* alpha,
                 const DType* A, const int lda, DType* B, const int ldb) {
  rb_raise(rb_eNotImpError, "trmm not yet implemented for non-BLAS dtypes");
}

template <>
inline void trmm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE ta, const enum CBLAS_DIAG diag, const int m, const int n, const float* alpha,
                 const float* A, const int lda, float* B, const int ldb) {
  cblas_strmm(order, side, uplo, ta, diag, m, n, *alpha, A, lda, B, ldb);
}

template <>
inline void trmm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE ta, const enum CBLAS_DIAG diag, const int m, const int n, const double* alpha,
                 const double* A, const int lda, double* B, const int ldb) {
  cblas_dtrmm(order, side, uplo, ta, diag, m, n, *alpha, A, lda, B, ldb);
}

template <>
inline void trmm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE ta, const enum CBLAS_DIAG diag, const int m, const int n, const Complex64* alpha,
                 const Complex64* A, const int lda, Complex64* B, const int ldb) {
  cblas_ctrmm(order, side, uplo, ta, diag, m, n, alpha, A, lda, B, ldb);
}

template <>
inline void trmm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE ta, const enum CBLAS_DIAG diag, const int m, const int n, const Complex128* alpha,
                 const Complex128* A, const int lda, Complex128* B, const int ldb) {
  cblas_ztrmm(order, side, uplo, ta, diag, m, n, alpha, A, lda, B, ldb);
}


/*
 * BLAS' DTRSM function, generalized.
 */
template <typename DType, typename = typename std::enable_if<!std::is_integral<DType>::value>::type>
inline void trsm(const enum CBLAS_ORDER order,
                 const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_DIAG diag,
                 const int m, const int n, const DType alpha, const DType* a,
                 const int lda, DType* b, const int ldb)
{
  /*using std::cerr;
  using std::endl;*/

  int                     num_rows_a = n;
  if (side == CblasLeft)  num_rows_a = m;

  if (lda < std::max(1,num_rows_a)) {
    fprintf(stderr, "TRSM: num_rows_a = %d; got lda=%d\n", num_rows_a, lda);
    rb_raise(rb_eArgError, "TRSM: Expected lda >= max(1, num_rows_a)");
  }

  // Test the input parameters.
  if (order == CblasRowMajor) {
    if (ldb < std::max(1,n)) {
      fprintf(stderr, "TRSM: M=%d; got ldb=%d\n", m, ldb);
      rb_raise(rb_eArgError, "TRSM: Expected ldb >= max(1,N)");
    }

    // For row major, need to switch side and uplo
    enum CBLAS_SIDE side_ = side == CblasLeft  ? CblasRight : CblasLeft;
    enum CBLAS_UPLO uplo_ = uplo == CblasUpper ? CblasLower : CblasUpper;

/*
    cerr << "(row-major) trsm: " << (side_ == CblasLeft ? "left " : "right ")
         << (uplo_ == CblasUpper ? "upper " : "lower ")
         << (trans_a == CblasTrans ? "trans " : "notrans ")
         << (diag == CblasNonUnit ? "nonunit " : "unit ")
         << n << " " << m << " " << alpha << " a " << lda << " b " << ldb << endl;
*/
    trsm_nothrow<DType>(side_, uplo_, trans_a, diag, n, m, alpha, a, lda, b, ldb);

  } else { // CblasColMajor

    if (ldb < std::max(1,m)) {
      fprintf(stderr, "TRSM: M=%d; got ldb=%d\n", m, ldb);
      rb_raise(rb_eArgError, "TRSM: Expected ldb >= max(1,M)");
    }
/*
    cerr << "(col-major) trsm: " << (side == CblasLeft ? "left " : "right ")
         << (uplo == CblasUpper ? "upper " : "lower ")
         << (trans_a == CblasTrans ? "trans " : "notrans ")
         << (diag == CblasNonUnit ? "nonunit " : "unit ")
         << m << " " << n << " " << alpha << " a " << lda << " b " << ldb << endl;
*/
    trsm_nothrow<DType>(side, uplo, trans_a, diag, m, n, alpha, a, lda, b, ldb);

  }

}


template <>
inline void trsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_DIAG diag,
                 const int m, const int n, const float alpha, const float* a,
                 const int lda, float* b, const int ldb)
{
  cblas_strsm(order, side, uplo, trans_a, diag, m, n, alpha, a, lda, b, ldb);
}

template <>
inline void trsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_DIAG diag,
                 const int m, const int n, const double alpha, const double* a,
                 const int lda, double* b, const int ldb)
{
/*  using std::cerr;
  using std::endl;
  cerr << "(row-major) dtrsm: " << (side == CblasLeft ? "left " : "right ")
       << (uplo == CblasUpper ? "upper " : "lower ")
       << (trans_a == CblasTrans ? "trans " : "notrans ")
       << (diag == CblasNonUnit ? "nonunit " : "unit ")
       << m << " " << n << " " << alpha << " a " << lda << " b " << ldb << endl;
*/
  cblas_dtrsm(order, side, uplo, trans_a, diag, m, n, alpha, a, lda, b, ldb);
}


template <>
inline void trsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_DIAG diag,
                 const int m, const int n, const Complex64 alpha, const Complex64* a,
                 const int lda, Complex64* b, const int ldb)
{
  cblas_ctrsm(order, side, uplo, trans_a, diag, m, n, (const void*)(&alpha), (const void*)(a), lda, (void*)(b), ldb);
}

template <>
inline void trsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_DIAG diag,
                 const int m, const int n, const Complex128 alpha, const Complex128* a,
                 const int lda, Complex128* b, const int ldb)
{
  cblas_ztrsm(order, side, uplo, trans_a, diag, m, n, (const void*)(&alpha), (const void*)(a), lda, (void*)(b), ldb);
}


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
 * GEneral Matrix Multiplication: based on dgemm.f from Netlib.
 *
 * This is an extremely inefficient algorithm. Recommend using ATLAS' version instead.
 *
 * Template parameters: LT -- long version of type T. Type T is the matrix dtype.
 *
 * This version throws no errors. Use gemm<DType> instead for error checking.
 */
template <typename DType>
inline void gemm_nothrow(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
                 const DType* alpha, const DType* A, const int lda, const DType* B, const int ldb, const DType* beta, DType* C, const int ldc)
{

  typename LongDType<DType>::type temp;

  // Quick return if possible
  if (!M or !N or ((*alpha == 0 or !K) and *beta == 1)) return;

  // For alpha = 0
  if (*alpha == 0) {
    if (*beta == 0) {
      for (int j = 0; j < N; ++j)
        for (int i = 0; i < M; ++i) {
          C[i+j*ldc] = 0;
        }
    } else {
      for (int j = 0; j < N; ++j)
        for (int i = 0; i < M; ++i) {
          C[i+j*ldc] *= *beta;
        }
    }
    return;
  }

  // Start the operations
  if (TransB == CblasNoTrans) {
    if (TransA == CblasNoTrans) {
      // C = alpha*A*B+beta*C
      for (int j = 0; j < N; ++j) {
        if (*beta == 0) {
          for (int i = 0; i < M; ++i) {
            C[i+j*ldc] = 0;
          }
        } else if (*beta != 1) {
          for (int i = 0; i < M; ++i) {
            C[i+j*ldc] *= *beta;
          }
        }

        for (int l = 0; l < K; ++l) {
          if (B[l+j*ldb] != 0) {
            temp = *alpha * B[l+j*ldb];
            for (int i = 0; i < M; ++i) {
              C[i+j*ldc] += A[i+l*lda] * temp;
            }
          }
        }
      }

    } else {

      // C = alpha*A**DType*B + beta*C
      for (int j = 0; j < N; ++j) {
        for (int i = 0; i < M; ++i) {
          temp = 0;
          for (int l = 0; l < K; ++l) {
            temp += A[l+i*lda] * B[l+j*ldb];
          }

          if (*beta == 0) {
            C[i+j*ldc] = *alpha*temp;
          } else {
            C[i+j*ldc] = *alpha*temp + *beta*C[i+j*ldc];
          }
        }
      }

    }

  } else if (TransA == CblasNoTrans) {

    // C = alpha*A*B**T + beta*C
    for (int j = 0; j < N; ++j) {
      if (*beta == 0) {
        for (int i = 0; i < M; ++i) {
          C[i+j*ldc] = 0;
        }
      } else if (*beta != 1) {
        for (int i = 0; i < M; ++i) {
          C[i+j*ldc] *= *beta;
        }
      }

      for (int l = 0; l < K; ++l) {
        if (B[j+l*ldb] != 0) {
          temp = *alpha * B[j+l*ldb];
          for (int i = 0; i < M; ++i) {
            C[i+j*ldc] += A[i+l*lda] * temp;
          }
        }
      }

    }

  } else {

    // C = alpha*A**DType*B**T + beta*C
    for (int j = 0; j < N; ++j) {
      for (int i = 0; i < M; ++i) {
        temp = 0;
        for (int l = 0; l < K; ++l) {
          temp += A[l+i*lda] * B[j+l*ldb];
        }

        if (*beta == 0) {
          C[i+j*ldc] = *alpha*temp;
        } else {
          C[i+j*ldc] = *alpha*temp + *beta*C[i+j*ldc];
        }
      }
    }

  }

  return;
}



template <typename DType>
inline void gemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
                 const DType* alpha, const DType* A, const int lda, const DType* B, const int ldb, const DType* beta, DType* C, const int ldc)
{
  if (Order == CblasRowMajor) {
    if (TransA == CblasNoTrans) {
      if (lda < std::max(K,1)) {
        rb_raise(rb_eArgError, "lda must be >= MAX(K,1): lda=%d K=%d", lda, K);
      }
    } else {
      if (lda < std::max(M,1)) { // && TransA == CblasTrans
        rb_raise(rb_eArgError, "lda must be >= MAX(M,1): lda=%d M=%d", lda, M);
      }
    }

    if (TransB == CblasNoTrans) {
      if (ldb < std::max(N,1)) {
        rb_raise(rb_eArgError, "ldb must be >= MAX(N,1): ldb=%d N=%d", ldb, N);
      }
    } else {
      if (ldb < std::max(K,1)) {
        rb_raise(rb_eArgError, "ldb must be >= MAX(K,1): ldb=%d K=%d", ldb, K);
      }
    }

    if (ldc < std::max(N,1)) {
      rb_raise(rb_eArgError, "ldc must be >= MAX(N,1): ldc=%d N=%d", ldc, N);
    }
  } else { // CblasColMajor
    if (TransA == CblasNoTrans) {
      if (lda < std::max(M,1)) {
        rb_raise(rb_eArgError, "lda must be >= MAX(M,1): lda=%d M=%d", lda, M);
      }
    } else {
      if (lda < std::max(K,1)) { // && TransA == CblasTrans
        rb_raise(rb_eArgError, "lda must be >= MAX(K,1): lda=%d K=%d", lda, K);
      }
    }

    if (TransB == CblasNoTrans) {
      if (ldb < std::max(K,1)) {
        rb_raise(rb_eArgError, "ldb must be >= MAX(K,1): ldb=%d N=%d", ldb, K);
      }
    } else {
      if (ldb < std::max(N,1)) { // NOTE: This error message is actually wrong in the ATLAS source currently. Or are we wrong?
        rb_raise(rb_eArgError, "ldb must be >= MAX(N,1): ldb=%d N=%d", ldb, N);
      }
    }

    if (ldc < std::max(M,1)) {
      rb_raise(rb_eArgError, "ldc must be >= MAX(M,1): ldc=%d N=%d", ldc, M);
    }
  }

  /*
   * Call SYRK when that's what the user is actually asking for; just handle beta=0, because beta=X requires
   * we copy C and then subtract to preserve asymmetry.
   */

  if (A == B && M == N && TransA != TransB && lda == ldb && beta == 0) {
    rb_raise(rb_eNotImpError, "syrk and syreflect not implemented");
    /*syrk<DType>(CblasUpper, (Order == CblasColMajor) ? TransA : TransB, N, K, alpha, A, lda, beta, C, ldc);
    syreflect(CblasUpper, N, C, ldc);
    */
  }

  if (Order == CblasRowMajor)    gemm_nothrow<DType>(TransB, TransA, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
  else                           gemm_nothrow<DType>(TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);

}


template <>
inline void gemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const float* alpha, const float* A, const int lda, const float* B, const int ldb, const float* beta, float* C, const int ldc) {
  cblas_sgemm(Order, TransA, TransB, M, N, K, *alpha, A, lda, B, ldb, *beta, C, ldc);
}

template <>
inline void gemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const double* alpha, const double* A, const int lda, const double* B, const int ldb, const double* beta, double* C, const int ldc) {
  cblas_dgemm(Order, TransA, TransB, M, N, K, *alpha, A, lda, B, ldb, *beta, C, ldc);
}

template <>
inline void gemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const Complex64* alpha, const Complex64* A, const int lda, const Complex64* B, const int ldb, const Complex64* beta, Complex64* C, const int ldc) {
  cblas_cgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

template <>
inline void gemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const Complex128* alpha, const Complex128* A, const int lda, const Complex128* B, const int ldb, const Complex128* beta, Complex128* C, const int ldc) {
  cblas_zgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


/*
 * GEneral Matrix-Vector multiplication: based on dgemv.f from Netlib.
 *
 * This is an extremely inefficient algorithm. Recommend using ATLAS' version instead.
 *
 * Template parameters: LT -- long version of type T. Type T is the matrix dtype.
 */
template <typename DType>
inline bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const DType* alpha, const DType* A, const int lda,
          const DType* X, const int incX, const DType* beta, DType* Y, const int incY) {
  int lenX, lenY, i, j;
  int kx, ky, iy, jx, jy, ix;

  typename LongDType<DType>::type temp;

  // Test the input parameters
  if (Trans < 111 || Trans > 113) {
    rb_raise(rb_eArgError, "GEMV: TransA must be CblasNoTrans, CblasTrans, or CblasConjTrans");
    return false;
  } else if (lda < std::max(1, N)) {
    fprintf(stderr, "GEMV: N = %d; got lda=%d", N, lda);
    rb_raise(rb_eArgError, "GEMV: Expected lda >= max(1, N)");
    return false;
  } else if (incX == 0) {
    rb_raise(rb_eArgError, "GEMV: Expected incX != 0\n");
    return false;
  } else if (incY == 0) {
    rb_raise(rb_eArgError, "GEMV: Expected incY != 0\n");
    return false;
  }

  // Quick return if possible
  if (!M or !N or (*alpha == 0 and *beta == 1)) return true;

  if (Trans == CblasNoTrans) {
    lenX = N;
    lenY = M;
  } else {
    lenX = M;
    lenY = N;
  }

  if (incX > 0) kx = 0;
  else          kx = (lenX - 1) * -incX;

  if (incY > 0) ky = 0;
  else          ky =  (lenY - 1) * -incY;

  // Start the operations. In this version, the elements of A are accessed sequentially with one pass through A.
  if (*beta != 1) {
    if (incY == 1) {
      if (*beta == 0) {
        for (i = 0; i < lenY; ++i) {
          Y[i] = 0;
        }
      } else {
        for (i = 0; i < lenY; ++i) {
          Y[i] *= *beta;
        }
      }
    } else {
      iy = ky;
      if (*beta == 0) {
        for (i = 0; i < lenY; ++i) {
          Y[iy] = 0;
          iy += incY;
        }
      } else {
        for (i = 0; i < lenY; ++i) {
          Y[iy] *= *beta;
          iy += incY;
        }
      }
    }
  }

  if (*alpha == 0) return false;

  if (Trans == CblasNoTrans) {

    // Form  y := alpha*A*x + y.
    jx = kx;
    if (incY == 1) {
      for (j = 0; j < N; ++j) {
        if (X[jx] != 0) {
          temp = *alpha * X[jx];
          for (i = 0; i < M; ++i) {
            Y[i] += A[j+i*lda] * temp;
          }
        }
        jx += incX;
      }
    } else {
      for (j = 0; j < N; ++j) {
        if (X[jx] != 0) {
          temp = *alpha * X[jx];
          iy = ky;
          for (i = 0; i < M; ++i) {
            Y[iy] += A[j+i*lda] * temp;
            iy += incY;
          }
        }
        jx += incX;
      }
    }

  } else { // TODO: Check that indices are correct! They're switched for C.

    // Form  y := alpha*A**DType*x + y.
    jy = ky;

    if (incX == 1) {
      for (j = 0; j < N; ++j) {
        temp = 0;
        for (i = 0; i < M; ++i) {
          temp += A[j+i*lda]*X[j];
        }
        Y[jy] += *alpha * temp;
        jy += incY;
      }
    } else {
      for (j = 0; j < N; ++j) {
        temp = 0;
        ix = kx;
        for (i = 0; i < M; ++i) {
          temp += A[j+i*lda] * X[ix];
          ix += incX;
        }

        Y[jy] += *alpha * temp;
        jy += incY;
      }
    }
  }

  return true;
}  // end of GEMV

template <>
inline bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const float* alpha, const float* A, const int lda,
          const float* X, const int incX, const float* beta, float* Y, const int incY) {
  cblas_sgemv(CblasRowMajor, Trans, M, N, *alpha, A, lda, X, incX, *beta, Y, incY);
  return true;
}

template <>
inline bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const double* alpha, const double* A, const int lda,
          const double* X, const int incX, const double* beta, double* Y, const int incY) {
  cblas_dgemv(CblasRowMajor, Trans, M, N, *alpha, A, lda, X, incX, *beta, Y, incY);
  return true;
}

template <>
inline bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const Complex64* alpha, const Complex64* A, const int lda,
          const Complex64* X, const int incX, const Complex64* beta, Complex64* Y, const int incY) {
  cblas_cgemv(CblasRowMajor, Trans, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  return true;
}

template <>
inline bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const Complex128* alpha, const Complex128* A, const int lda,
          const Complex128* X, const int incX, const Complex128* beta, Complex128* Y, const int incY) {
  cblas_zgemv(CblasRowMajor, Trans, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  return true;
}


// Yale: numeric matrix multiply c=a*b
template <typename DType, typename IType>
inline void numbmm(const unsigned int n, const unsigned int m, const unsigned int l, const IType* ia, const IType* ja, const DType* a, const bool diaga,
            const IType* ib, const IType* jb, const DType* b, const bool diagb, IType* ic, IType* jc, DType* c, const bool diagc) {
  const unsigned int max_lmn = std::max(std::max(m, n), l);
  IType next[max_lmn];
  DType sums[max_lmn];

  DType v;

  IType head, length, temp, ndnz = 0;
  IType minmn = std::min(m,n);
  IType minlm = std::min(l,m);

  for (IType idx = 0; idx < max_lmn; ++idx) { // initialize scratch arrays
    next[idx] = std::numeric_limits<IType>::max();
    sums[idx] = 0;
  }

  for (IType i = 0; i < n; ++i) { // walk down the rows
    head = std::numeric_limits<IType>::max()-1; // head gets assigned as whichever column of B's row j we last visited
    length = 0;

    for (IType jj = ia[i]; jj <= ia[i+1]; ++jj) { // walk through entries in each row
      IType j;

      if (jj == ia[i+1]) { // if we're in the last entry for this row:
        if (!diaga || i >= minmn) continue;
        j   = i;      // if it's a new Yale matrix, and last entry, get the diagonal position (j) and entry (ajj)
        v   = a[i];
      } else {
        j   = ja[jj]; // if it's not the last entry for this row, get the column (j) and entry (ajj)
        v   = a[jj];
      }

      for (IType kk = ib[j]; kk <= ib[j+1]; ++kk) {

        IType k;

        if (kk == ib[j+1]) { // Get the column id for that entry
          if (!diagb || j >= minlm) continue;
          k  = j;
          sums[k] += v*b[k];
        } else {
          k  = jb[kk];
          sums[k] += v*b[kk];
        }

        if (next[k] == std::numeric_limits<IType>::max()) {
          next[k] = head;
          head    = k;
          ++length;
        }
      } // end of kk loop
    } // end of jj loop

    for (IType jj = 0; jj < length; ++jj) {
      if (sums[head] != 0) {
        if (diagc && head == i) {
          c[head] = sums[head];
        } else {
          jc[n+1+ndnz] = head;
          c[n+1+ndnz]  = sums[head];
          ++ndnz;
        }
      }

      temp = head;
      head = next[head];

      next[temp] = std::numeric_limits<IType>::max();
      sums[temp] = 0;
    }

    ic[i+1] = n+1+ndnz;
  }
} /* numbmm_ */


/*
template <typename DType, typename IType>
inline void new_yale_matrix_multiply(const unsigned int m, const IType* ija, const DType* a, const IType* ijb, const DType* b, YALE_STORAGE* c_storage) {
  unsigned int n = c_storage->shape[0],
               l = c_storage->shape[1];

  // Create a working vector of dimension max(m,l,n) and initial value IType::max():
  std::vector<IType> mask(std::max(std::max(m,l),n), std::numeric_limits<IType>::max());

  for (IType i = 0; i < n; ++i) { // A.rows.each_index do |i|

    IType j, k;
    size_t ndnz;

    for (IType jj = ija[i]; jj <= ija[i+1]; ++jj) { // walk through column pointers for row i of A
      j = (jj == ija[i+1]) ? i : ija[jj];   // Get the current column index (handle diagonals last)

      if (j >= m) {
        if (j == ija[jj]) rb_raise(rb_eIndexError, "ija array for left-hand matrix contains an out-of-bounds column index %u at position %u", jj, j);
        else              break;
      }

      for (IType kk = ijb[j]; kk <= ijb[j+1]; ++kk) { // walk through column pointers for row j of B
        if (j >= m) continue; // first of all, does B *have* a row j?
        k = (kk == ijb[j+1]) ? j : ijb[kk];   // Get the current column index (handle diagonals last)

        if (k >= l) {
          if (k == ijb[kk]) rb_raise(rb_eIndexError, "ija array for right-hand matrix contains an out-of-bounds column index %u at position %u", kk, k);
          else              break;
        }

        if (mask[k] == )
      }

    }
  }
}
*/

// Yale: Symbolic matrix multiply c=a*b
template <typename IType>
inline size_t symbmm(const unsigned int n, const unsigned int m, const unsigned int l, const IType* ia, const IType* ja, const bool diaga,
            const IType* ib, const IType* jb, const bool diagb, IType* ic, const bool diagc) {
  unsigned int max_lmn = std::max(std::max(m,n), l);
  IType mask[max_lmn];  // INDEX in the SMMP paper.
  IType j, k; /* Local variables */
  size_t ndnz = n;

  for (IType idx = 0; idx < max_lmn; ++idx)
    mask[idx] = std::numeric_limits<IType>::max();

  if (ic) { // Only write to ic if it's supplied; otherwise, we're just counting.
    if (diagc)  ic[0] = n+1;
    else        ic[0] = 0;
  }

  IType minmn = std::min(m,n);
  IType minlm = std::min(l,m);

  for (IType i = 0; i < n; ++i) { // MAIN LOOP: through rows

    for (IType jj = ia[i]; jj <= ia[i+1]; ++jj) { // merge row lists, walking through columns in each row

      // j <- column index given by JA[jj], or handle diagonal.
      if (jj == ia[i+1]) { // Don't really do it the last time -- just handle diagonals in a new yale matrix.
        if (!diaga || i >= minmn) continue;
        j = i;
      } else j = ja[jj];

      for (IType kk = ib[j]; kk <= ib[j+1]; ++kk) { // Now walk through columns K of row J in matrix B.
        if (kk == ib[j+1]) {
          if (!diagb || j >= minlm) continue;
          k = j;
        } else k = jb[kk];

        if (mask[k] != i) {
          mask[k] = i;
          ++ndnz;
        }
      }
    }

    if (diagc && mask[i] == std::numeric_limits<IType>::max()) --ndnz;

    if (ic) ic[i+1] = ndnz;
  }

  return ndnz;
} /* symbmm_ */


// In-place quicksort (from Wikipedia) -- called by smmp_sort_columns, below. All functions are inclusive of left, right.
namespace smmp_sort {
  const size_t THRESHOLD = 4;  // switch to insertion sort for 4 elements or fewer

  template <typename DType, typename IType>
  void print_array(DType* vals, IType* array, IType left, IType right) {
    for (IType i = left; i <= right; ++i) {
      std::cerr << array[i] << ":" << vals[i] << "  ";
    }
    std::cerr << std::endl;
  }

  template <typename DType, typename IType>
  IType partition(DType* vals, IType* array, IType left, IType right, IType pivot) {
    IType pivotJ = array[pivot];
    DType pivotV = vals[pivot];

    // Swap pivot and right
    array[pivot] = array[right];
    vals[pivot]  = vals[right];
    array[right] = pivotJ;
    vals[right]  = pivotV;

    IType store = left;
    for (IType idx = left; idx < right; ++idx) {
      if (array[idx] <= pivotJ) {
        // Swap i and store
        std::swap(array[idx], array[store]);
        std::swap(vals[idx],  vals[store]);
        ++store;
      }
    }

    std::swap(array[store], array[right]);
    std::swap(vals[store],  vals[right]);

    return store;
  }

  // Recommended to use the median of left, right, and mid for the pivot.
  template <typename IType>
  IType median(IType a, IType b, IType c) {
    if (a < b) {
      if (b < c) return b; // a b c
      if (a < c) return c; // a c b
                 return a; // c a b

    } else { // a > b
      if (a < c) return a; // b a c
      if (b < c) return c; // b c a
                 return b; // c b a
    }
  }


  // Insertion sort is more efficient than quicksort for small N
  template <typename DType, typename IType>
  void insertion_sort(DType* vals, IType* array, IType left, IType right) {
    for (IType idx = left; idx <= right; ++idx) {
      IType col_to_insert = array[idx];
      DType val_to_insert = vals[idx];

      IType hole_pos = idx;
      for (; hole_pos > left && col_to_insert < array[hole_pos-1]; --hole_pos) {
        array[hole_pos] = array[hole_pos - 1];  // shift the larger column index up
        vals[hole_pos]  = vals[hole_pos - 1];   // value goes along with it
      }

      array[hole_pos] = col_to_insert;
      vals[hole_pos]  = val_to_insert;
    }
  }


  template <typename DType, typename IType>
  void quicksort(DType* vals, IType* array, IType left, IType right) {

    if (left < right) {
      if (right - left < THRESHOLD) {
        insertion_sort(vals, array, left, right);
      } else {
        // choose any pivot such that left < pivot < right
        IType pivot = median(left, right, (IType)(((unsigned long)left + (unsigned long)right) / 2));
        pivot = partition(vals, array, left, right, pivot);

        // recursively sort elements smaller than the pivot
        quicksort<DType,IType>(vals, array, left, pivot-1);

        // recursively sort elements at least as big as the pivot
        quicksort<DType,IType>(vals, array, pivot+1, right);
      }
    }
  }


}; // end of namespace smmp_sort


/*
 * For use following symbmm and numbmm. Sorts the matrix entries in each row according to the column index.
 * This utilizes quicksort, which is an in-place unstable sort (since there are no duplicate entries, we don't care
 * about stability).
 *
 * TODO: It might be worthwhile to do a test for free memory, and if available, use an unstable sort that isn't in-place.
 *
 * TODO: It's actually probably possible to write an even faster sort, since symbmm/numbmm are not producing a random
 * ordering. If someone is doing a lot of Yale matrix multiplication, it might benefit them to consider even insertion
 * sort.
 */
template <typename DType, typename IType>
inline void smmp_sort_columns(const size_t n, const IType* ia, IType* ja, DType* a) {
  for (size_t i = 0; i < n; ++i) {
    if (ia[i+1] - ia[i] < 2) continue; // no need to sort rows containing only one or two elements.
    else if (ia[i+1] - ia[i] <= smmp_sort::THRESHOLD) {
      smmp_sort::insertion_sort<DType, IType>(a, ja, ia[i], ia[i+1]-1); // faster for small rows
    } else {
      smmp_sort::quicksort<DType, IType>(a, ja, ia[i], ia[i+1]-1);      // faster for large rows (and may call insertion_sort as well)
    }
  }
}



/*
 * Transposes a generic Yale matrix (old or new). Specify new by setting diaga = true.
 *
 * Based on transp from SMMP (same as symbmm and numbmm).
 *
 * This is not named in the same way as most yale_storage functions because it does not act on a YALE_STORAGE
 * object.
 */
template <typename DType, typename IType>
void transpose_yale(const size_t n, const size_t m, const void* ia_, const void* ja_, const void* a_,
                             const bool diaga, void* ib_, void* jb_, void* b_, const bool move)
{
  const IType *ia = reinterpret_cast<const IType*>(ia_),
              *ja = reinterpret_cast<const IType*>(ja_);
  const DType *a  = reinterpret_cast<const DType*>(a_);

  IType *ib = reinterpret_cast<IType*>(ib_),
        *jb = reinterpret_cast<IType*>(jb_);
  DType *b  = reinterpret_cast<DType*>(b_);



  size_t index;

  // Clear B
  for (size_t i = 0; i < m+1; ++i) ib[i] = 0;

  if (move)
    for (size_t i = 0; i < m+1; ++i) b[i] = 0;

  if (diaga) ib[0] = m + 1;
  else       ib[0] = 0;

  /* count indices for each column */

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = ia[i]; j < ia[i+1]; ++j) {
      ++(ib[ja[j]+1]);
    }
  }

  for (size_t i = 0; i < m; ++i) {
    ib[i+1] = ib[i] + ib[i+1];
  }

  /* now make jb */

  for (size_t i = 0; i < n; ++i) {

    for (size_t j = ia[i]; j < ia[i+1]; ++j) {
      index = ja[j];
      jb[ib[index]] = i;

      if (move)
        b[ib[index]] = a[j];

      ++(ib[index]);
    }
  }

  /* now fixup ib */

  for (size_t i = m; i >= 1; --i) {
    ib[i] = ib[i-1];
  }


  if (diaga) {
    if (move) {
      size_t j = std::min(n,m);

      for (size_t i = 0; i < j; ++i) {
        b[i] = a[i];
      }
    }
    ib[0] = m + 1;

  } else {
    ib[0] = 0;
  }
}


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

  // Symbols used by ATLAS:
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

    int i = *ipiv = nm::math::lapack::idamax<DType>(N, A, 1); // cblas_iamax(N, A, 1);

    DType tmp = A[i];
    if (tmp != 0) {

      nm::math::lapack::scal<DType>((RowMajor ? N : M), nm::math::numeric_inverse(tmp), A, 1);
      A[i] = *A;
      *A   = tmp;

    } else ierr = 1;

  }
  return(ierr);
}

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
inline int potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, DType* A, const int lda) {
#ifdef HAVE_CLAPACK_H
  rb_raise(rb_eNotImpError, "not yet implemented for non-BLAS dtypes");
#else
  rb_raise(rb_eNotImpError, "only LAPACK version implemented thus far");
#endif
  return 0;
}

#ifdef HAVE_CLAPACK_H
template <>
inline int potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, float* A, const int lda) {
  return clapack_spotrf(order, uplo, N, A, lda);
}

template <>
inline int potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, double* A, const int lda) {
  return clapack_dpotrf(order, uplo, N, A, lda);
}

template <>
inline int potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, Complex64* A, const int lda) {
  return clapack_cpotrf(order, uplo, N, reinterpret_cast<void*>(A), lda);
}

template <>
inline int potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, Complex128* A, const int lda) {
  return clapack_zpotrf(order, uplo, N, reinterpret_cast<void*>(A), lda);
}
#endif


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


// Copies an upper row-major array from U, zeroing U; U is unit, so diagonal is not copied.
//
// From ATLAS 3.8.0.
template <typename DType>
static inline void trcpzeroU(const int M, const int N, DType* U, const int ldu, DType* C, const int ldc) {

  for (int i = 0; i != M; ++i) {
    for (int j = i+1; j < N; ++j) {
      C[j] = U[j];
      U[j] = 0;
    }

    C += ldc;
    U += ldu;
  }
}


/*
 * Un-comment the following lines when we figure out how to calculate NB for each of the ATLAS-derived
 * functions. This is probably really complicated.
 *
 * Also needed: ATL_MulByNB, ATL_DivByNB (both defined in the build process for ATLAS), and ATL_mmMU.
 *
 */

/*

template <bool RowMajor, bool Upper, typename DType>
static int trtri_4(const enum CBLAS_DIAG Diag, DType* A, const int lda) {

  if (RowMajor) {
    DType *pA0 = A, *pA1 = A+lda, *pA2 = A+2*lda, *pA3 = A+3*lda;
    DType tmp;
    if (Upper) {
      DType A01 = pA0[1], A02 = pA0[2], A03 = pA0[3],
                          A12 = pA1[2], A13 = pA1[3],
                                        A23 = pA2[3];

      if (Diag == CblasNonUnit) {
        pA0->inverse();
        (pA1+1)->inverse();
        (pA2+2)->inverse();
        (pA3+3)->inverse();

        pA0[1] = -A01 * pA1[1] * pA0[0];
        pA1[2] = -A12 * pA2[2] * pA1[1];
        pA2[3] = -A23 * pA3[3] * pA2[2];

        pA0[2] = -(A01 * pA1[2] + A02 * pA2[2]) * pA0[0];
        pA1[3] = -(A12 * pA2[3] + A13 * pA3[3]) * pA1[1];

        pA0[3] = -(A01 * pA1[3] + A02 * pA2[3] + A03 * pA3[3]) * pA0[0];

      } else {

        pA0[1] = -A01;
        pA1[2] = -A12;
        pA2[3] = -A23;

        pA0[2] = -(A01 * pA1[2] + A02);
        pA1[3] = -(A12 * pA2[3] + A13);

        pA0[3] = -(A01 * pA1[3] + A02 * pA2[3] + A03);
      }

    } else { // Lower
      DType A10 = pA1[0],
            A20 = pA2[0], A21 = pA2[1],
            A30 = PA3[0], A31 = pA3[1], A32 = pA3[2];
      DType *B10 = pA1,
            *B20 = pA2,
            *B30 = pA3,
            *B21 = pA2+1,
            *B31 = pA3+1,
            *B32 = pA3+2;


      if (Diag == CblasNonUnit) {
        pA0->inverse();
        (pA1+1)->inverse();
        (pA2+2)->inverse();
        (pA3+3)->inverse();

        *B10 = -A10 * pA0[0] * pA1[1];
        *B21 = -A21 * pA1[1] * pA2[2];
        *B32 = -A32 * pA2[2] * pA3[3];
        *B20 = -(A20 * pA0[0] + A21 * (*B10)) * pA2[2];
        *B31 = -(A31 * pA1[1] + A32 * (*B21)) * pA3[3];
        *B30 = -(A30 * pA0[0] + A31 * (*B10) + A32 * (*B20)) * pA3;
      } else {
        *B10 = -A10;
        *B21 = -A21;
        *B32 = -A32;
        *B20 = -(A20 + A21 * (*B10));
        *B31 = -(A31 + A32 * (*B21));
        *B30 = -(A30 + A31 * (*B10) + A32 * (*B20));
      }
    }

  } else {
    rb_raise(rb_eNotImpError, "only row-major implemented at this time");
  }

  return 0;

}


template <bool RowMajor, bool Upper, typename DType>
static int trtri_3(const enum CBLAS_DIAG Diag, DType* A, const int lda) {

  if (RowMajor) {

    DType tmp;

    if (Upper) {
      DType A01 = pA0[1], A02 = pA0[2], A03 = pA0[3],
                          A12 = pA1[2], A13 = pA1[3];

      DType *B01 = pA0 + 1,
            *B02 = pA0 + 2,
            *B12 = pA1 + 2;

      if (Diag == CblasNonUnit) {
        pA0->inverse();
        (pA1+1)->inverse();
        (pA2+2)->inverse();

        *B01 = -A01 * pA1[1] * pA0[0];
        *B12 = -A12 * pA2[2] * pA1[1];
        *B02 = -(A01 * (*B12) + A02 * pA2[2]) * pA0[0];
      } else {
        *B01 = -A01;
        *B12 = -A12;
        *B02 = -(A01 * (*B12) + A02);
      }

    } else { // Lower
      DType *pA0=A, *pA1=A+lda, *pA2=A+2*lda;
      DType A10=pA1[0],
            A20=pA2[0], A21=pA2[1];

      DType *B10 = pA1,
            *B20 = pA2;
            *B21 = pA2+1;

      if (Diag == CblasNonUnit) {
        pA0->inverse();
        (pA1+1)->inverse();
        (pA2+2)->inverse();
        *B10 = -A10 * pA0[0] * pA1[1];
        *B21 = -A21 * pA1[1] * pA2[2];
        *B20 = -(A20 * pA0[0] + A21 * (*B10)) * pA2[2];
      } else {
        *B10 = -A10;
        *B21 = -A21;
        *B20 = -(A20 + A21 * (*B10));
      }
    }


  } else {
    rb_raise(rb_eNotImpError, "only row-major implemented at this time");
  }

  return 0;

}

template <bool RowMajor, bool Upper, bool Real, typename DType>
static void trtri(const enum CBLAS_DIAG Diag, const int N, DType* A, const int lda) {
  DType *Age, *Atr;
  DType tmp;
  int Nleft, Nright;

  int ierr = 0;

  static const DType ONE = 1;
  static const DType MONE -1;
  static const DType NONE = -1;

  if (RowMajor) {

    // FIXME: Use REAL_RECURSE_LIMIT here for float32 and float64 (instead of 1)
    if ((Real && N > REAL_RECURSE_LIMIT) || (N > 1)) {
      Nleft = N >> 1;
#ifdef NB
      if (Nleft > NB) NLeft = ATL_MulByNB(ATL_DivByNB(Nleft));
#endif

      Nright = N - Nleft;

      if (Upper) {
        Age = A + Nleft;
        Atr = A + (Nleft * (lda+1));

        nm::math::trsm<DType>(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, Diag,
                              Nleft, Nright, ONE, Atr, lda, Age, lda);

        nm::math::trsm<DType>(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, Diag,
                              Nleft, Nright, MONE, A, lda, Age, lda);

      } else { // Lower
        Age = A + ((Nleft*lda));
        Atr = A + (Nleft * (lda+1));

        nm::math::trsm<DType>(CblasRowMajor, CblasRight, CblasLower, CblasNoTrans, Diag,
                              Nright, Nleft, ONE, A, lda, Age, lda);
        nm::math::trsm<DType>(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, Diag,
                              Nright, Nleft, MONE, Atr, lda, Age, lda);
      }

      ierr = trtri<RowMajor,Upper,Real,DType>(Diag, Nleft, A, lda);
      if (ierr) return ierr;

      ierr = trtri<RowMajor,Upper,Real,DType>(Diag, Nright, Atr, lda);
      if (ierr) return ierr + Nleft;

    } else {
      if (Real) {
        if (N == 4) {
          return trtri_4<RowMajor,Upper,Real,DType>(Diag, A, lda);
        } else if (N == 3) {
          return trtri_3<RowMajor,Upper,Real,DType>(Diag, A, lda);
        } else if (N == 2) {
          if (Diag == CblasNonUnit) {
            A->inverse();
            (A+(lda+1))->inverse();

            if (Upper) {
              *(A+1)     *=   *A;         // TRI_MUL
              *(A+1)     *=   *(A+lda+1); // TRI_MUL
            } else {
              *(A+lda)   *=   *A;         // TRI_MUL
              *(A+lda)   *=   *(A+lda+1); // TRI_MUL
            }
          }

          if (Upper) *(A+1)   = -*(A+1);      // TRI_NEG
          else       *(A+lda) = -*(A+lda);    // TRI_NEG
        } else if (Diag == CblasNonUnit) A->inverse();
      } else { // not real
        if (Diag == CblasNonUnit) A->inverse();
      }
    }

  } else {
    rb_raise(rb_eNotImpError, "only row-major implemented at this time");
  }

  return ierr;
}


template <bool RowMajor, bool Real, typename DType>
int getri(const int N, DType* A, const int lda, const int* ipiv, DType* wrk, const int lwrk) {

  if (!RowMajor) rb_raise(rb_eNotImpError, "only row-major implemented at this time");

  int jb, nb, I, ndown, iret;

  const DType ONE = 1, NONE = -1;

  int iret = trtri<RowMajor,false,Real,DType>(CblasNonUnit, N, A, lda);
  if (!iret && N > 1) {
    jb = lwrk / N;
    if (jb >= NB) nb = ATL_MulByNB(ATL_DivByNB(jb));
    else if (jb >= ATL_mmMU) nb = (jb/ATL_mmMU)*ATL_mmMU;
    else nb = jb;
    if (!nb) return -6; // need at least 1 row of workspace

    // only first iteration will have partial block, unroll it

    jb = N - (N/nb) * nb;
    if (!jb) jb = nb;
    I = N - jb;
    A += lda * I;
    trcpzeroU<DType>(jb, jb, A+I, lda, wrk, jb);
    nm::math::trsm<DType>(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
                          jb, N, ONE, wrk, jb, A, lda);

    if (I) {
      do {
        I -= nb;
        A -= nb * lda;
        ndown = N-I;
        trcpzeroU<DType>(nb, ndown, A+I, lda, wrk, ndown);
        nm::math::gemm<DType>(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
                              nb, N, ONE, wrk, ndown, A, lda);
      } while (I);
    }

    // Apply row interchanges

    for (I = N - 2; I >= 0; --I) {
      jb = ipiv[I];
      if (jb != I) nm::math::swap<DType>(N, A+I*lda, 1, A+jb*lda, 1);
    }
  }

  return iret;
}
*/


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


/* Givens plane rotation. From ATLAS 3.8.4. */
// FIXME: Need a specialized algorithm for Rationals. BLAS' algorithm simply will not work for most values due to the
// FIXME: sqrt.
template <typename DType>
inline void rotg(DType* a, DType* b, DType* c, DType* s) {
  DType aa    = std::abs(*a), ab = std::abs(*b);
  DType roe   = aa > ab ? *a : *b;
  DType scal  = aa + ab;

  if (scal == 0) {
    *c =  1;
    *s = *a = *b = 0;
  } else {
    DType t0  = aa / scal, t1 = ab / scal;
    DType r   = scal * std::sqrt(t0 * t0 + t1 * t1);
    if (roe < 0) r = -r;
    *c = *a / r;
    *s = *b / r;
    DType z   = (*c != 0) ? (1 / *c) : DType(1);
    *a = r;
    *b = z;
  }
}

template <>
inline void rotg(float* a, float* b, float* c, float* s) {
  cblas_srotg(a, b, c, s);
}

template <>
inline void rotg(double* a, double* b, double* c, double* s) {
  cblas_drotg(a, b, c, s);
}

template <>
inline void rotg(Complex64* a, Complex64* b, Complex64* c, Complex64* s) {
  cblas_crotg(reinterpret_cast<void*>(a), reinterpret_cast<void*>(b), reinterpret_cast<void*>(c), reinterpret_cast<void*>(s));
}

template <>
inline void rotg(Complex128* a, Complex128* b, Complex128* c, Complex128* s) {
  cblas_zrotg(reinterpret_cast<void*>(a), reinterpret_cast<void*>(b), reinterpret_cast<void*>(c), reinterpret_cast<void*>(s));
}

template <typename DType>
inline void cblas_rotg(void* a, void* b, void* c, void* s) {
  rotg<DType>(reinterpret_cast<DType*>(a), reinterpret_cast<DType*>(b), reinterpret_cast<DType*>(c), reinterpret_cast<DType*>(s));
}


/* Applies a plane rotation. From ATLAS 3.8.4. */
template <typename DType, typename CSDType>
inline void rot(const int N, DType* X, const int incX, DType* Y, const int incY, const CSDType c, const CSDType s) {
  DType *x = X, *y = Y;
  int incx = incX, incy = incY;

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
  rot<DType,CSDType>(N, reinterpret_cast<DType*>(X), incX, reinterpret_cast<DType*>(Y), incY, *reinterpret_cast<const CSDType*>(c), *reinterpret_cast<const CSDType*>(s));
}


template <bool is_complex, typename DType>
inline void lauum(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, DType* A, const int lda) {

  int Nleft, Nright;
  const DType ONE = 1;
  DType *G, *U0 = A, *U1;

  if (N > 1) {
    Nleft = N >> 1;
    #ifdef NB
      if (Nleft > NB) Nleft = ATL_MulByNB(ATL_DivByNB(Nleft));
    #endif

    Nright = N - Nleft;

    // FIXME: There's a simpler way to write this next block, but I'm way too tired to work it out right now.
    if (uplo == CblasUpper) {
      if (order == CblasRowMajor) {
        G = A + Nleft;
        U1 = G + Nleft * lda;
      } else {
        G = A + Nleft * lda;
        U1 = G + Nleft;
      }
    } else {
      if (order == CblasRowMajor) {
        G = A + Nleft * lda;
        U1 = G + Nleft;
      } else {
        G = A + Nleft;
        U1 = G + Nleft * lda;
      }
    }

    lauum<is_complex, DType>(order, uplo, Nleft, U0, lda);

    if (is_complex) {

      nm::math::herk<DType>(order, uplo,
                            uplo == CblasLower ? CblasConjTrans : CblasNoTrans,
                            Nleft, Nright, &ONE, G, lda, &ONE, U0, lda);

      nm::math::trmm<DType>(order, CblasLeft, uplo, CblasConjTrans, CblasNonUnit, Nright, Nleft, &ONE, U1, lda, G, lda);
    } else {
      nm::math::syrk<DType>(order, uplo,
                            uplo == CblasLower ? CblasTrans : CblasNoTrans,
                            Nleft, Nright, &ONE, G, lda, &ONE, U0, lda);

      nm::math::trmm<DType>(order, CblasLeft, uplo, CblasTrans, CblasNonUnit, Nright, Nleft, &ONE, U1, lda, G, lda);
    }
    lauum<is_complex, DType>(order, uplo, Nright, U1, lda);

  } else {
    *A = *A * *A;
  }
}


#ifdef HAVE_CLAPACK_H
template <bool is_complex>
inline void lauum(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, float* A, const int lda) {
  clapack_slauum(order, uplo, N, A, lda);
}

template <bool is_complex>
inline void lauum(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, double* A, const int lda) {
  clapack_dlauum(order, uplo, N, A, lda);
}

template <bool is_complex>
inline void lauum(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, Complex64* A, const int lda) {
  clapack_clauum(order, uplo, N, A, lda);
}

template <bool is_complex>
inline void lauum(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, Complex128* A, const int lda) {
  clapack_zlauum(order, uplo, N, A, lda);
}
#endif


/*
* Function signature conversion for calling LAPACK's lauum functions as directly as possible.
*
* For documentation: http://www.netlib.org/lapack/double/dlauum.f
*
* This function should normally go in math.cpp, but we need it to be available to nmatrix.cpp.
*/
template <bool is_complex, typename DType>
inline int clapack_lauum(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, void* a, const int lda) {
  if (n < 0)              rb_raise(rb_eArgError, "n cannot be less than zero, is set to %d", n);
  if (lda < n || lda < 1) rb_raise(rb_eArgError, "lda must be >= max(n,1); lda=%d, n=%d\n", lda, n);

  if (uplo == CblasUpper) lauum<is_complex, DType>(order, uplo, n, reinterpret_cast<DType*>(a), lda);
  else                    lauum<is_complex, DType>(order, uplo, n, reinterpret_cast<DType*>(a), lda);

  return 0;
}




/*
 * Macro for declaring LAPACK specializations of the getrf function.
 *
 * type is the DType; call is the specific function to call; cast_as is what the DType* should be
 * cast to in order to pass it to LAPACK.
 */
#define LAPACK_GETRF(type, call, cast_as)                                     \
template <>                                                                   \
inline int getrf(const enum CBLAS_ORDER Order, const int M, const int N, type * A, const int lda, int* ipiv) { \
  int info = call(Order, M, N, reinterpret_cast<cast_as *>(A), lda, ipiv);    \
  if (!info) return info;                                                     \
  else {                                                                      \
    rb_raise(rb_eArgError, "getrf: problem with argument %d\n", info);        \
    return info;                                                              \
  }                                                                           \
}

/* Specialize for ATLAS types */
/*LAPACK_GETRF(float,      clapack_sgetrf, float)
LAPACK_GETRF(double,     clapack_dgetrf, double)
LAPACK_GETRF(Complex64,  clapack_cgetrf, void)
LAPACK_GETRF(Complex128, clapack_zgetrf, void)
*/


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


/*
* Function signature conversion for calling LAPACK's potrf functions as directly as possible.
*
* For documentation: http://www.netlib.org/lapack/double/dpotrf.f
*
* This function should normally go in math.cpp, but we need it to be available to nmatrix.cpp.
*/
template <typename DType>
inline int clapack_potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, void* a, const int lda) {
  return potrf<DType>(order, uplo, n, reinterpret_cast<DType*>(a), lda);
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

template <typename DType>
inline int getri(const enum CBLAS_ORDER order, const int n, DType* a, const int lda, const int* ipiv) {
  rb_raise(rb_eNotImpError, "getri not yet implemented for non-BLAS dtypes");
  return 0;
}

#ifdef HAVE_CLAPACK_H
template <>
inline int getri(const enum CBLAS_ORDER order, const int n, float* a, const int lda, const int* ipiv) {
  return clapack_sgetri(order, n, a, lda, ipiv);
}

template <>
inline int getri(const enum CBLAS_ORDER order, const int n, double* a, const int lda, const int* ipiv) {
  return clapack_dgetri(order, n, a, lda, ipiv);
}

template <>
inline int getri(const enum CBLAS_ORDER order, const int n, Complex64* a, const int lda, const int* ipiv) {
  return clapack_cgetri(order, n, reinterpret_cast<void*>(a), lda, ipiv);
}

template <>
inline int getri(const enum CBLAS_ORDER order, const int n, Complex128* a, const int lda, const int* ipiv) {
  return clapack_zgetri(order, n, reinterpret_cast<void*>(a), lda, ipiv);
}
#endif


template <typename DType>
inline int potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, DType* a, const int lda) {
  rb_raise(rb_eNotImpError, "potri not yet implemented for non-BLAS dtypes");
  return 0;
}


#ifdef HAVE_CLAPACK_H
template <>
inline int potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, float* a, const int lda) {
  return clapack_spotri(order, uplo, n, a, lda);
}

template <>
inline int potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, double* a, const int lda) {
  return clapack_dpotri(order, uplo, n, a, lda);
}

template <>
inline int potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, Complex64* a, const int lda) {
  return clapack_cpotri(order, uplo, n, reinterpret_cast<void*>(a), lda);
}

template <>
inline int potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, Complex128* a, const int lda) {
  return clapack_zpotri(order, uplo, n, reinterpret_cast<void*>(a), lda);
}
#endif

/*
 * Function signature conversion for calling LAPACK's getri functions as directly as possible.
 *
 * For documentation: http://www.netlib.org/lapack/double/dgetri.f
 *
 * This function should normally go in math.cpp, but we need it to be available to nmatrix.cpp.
 */
template <typename DType>
inline int clapack_getri(const enum CBLAS_ORDER order, const int n, void* a, const int lda, const int* ipiv) {
  return getri<DType>(order, n, reinterpret_cast<DType*>(a), lda, ipiv);
}


/*
 * Function signature conversion for calling LAPACK's potri functions as directly as possible.
 *
 * For documentation: http://www.netlib.org/lapack/double/dpotri.f
 *
 * This function should normally go in math.cpp, but we need it to be available to nmatrix.cpp.
 */
template <typename DType>
inline int clapack_potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, void* a, const int lda) {
  return potri<DType>(order, uplo, n, reinterpret_cast<DType*>(a), lda);
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



}} // end namespace nm::math


#endif // MATH_H
