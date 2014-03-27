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
// == gemm.h
//
// Header file for interface with ATLAS's CBLAS gemm functions and
// native templated version of LAPACK's gemm function.
//

#ifndef GEMM_H
# define GEMM_H

extern "C" { // These need to be in an extern "C" block or you'll get all kinds of undefined symbol errors.
#if defined HAVE_CBLAS_H
  #include <cblas.h>
#elif defined HAVE_ATLAS_CBLAS_H
  #include <atlas/cblas.h>
#endif
}


namespace nm { namespace math {
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

}} // end of namespace nm::math

#endif // GEMM_H
