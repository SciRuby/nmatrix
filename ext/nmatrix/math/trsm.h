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
// == trsm.h
//
// trsm function in native C++.
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

#ifndef TRSM_H
#define TRSM_H


extern "C" {
#if defined HAVE_CBLAS_H
  #include <cblas.h>
#elif defined HAVE_ATLAS_CBLAS_H
  #include <atlas/cblas.h>
#endif
}

namespace nm { namespace math {


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


} }  // namespace nm::math
#endif // TRSM_H
