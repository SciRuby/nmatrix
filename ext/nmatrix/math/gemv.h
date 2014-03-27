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
// == gemv.h
//
// Header file for interface with ATLAS's CBLAS gemv functions and
// native templated version of LAPACK's gemv function.
//

#ifndef GEMV_H
# define GEMV_H

extern "C" { // These need to be in an extern "C" block or you'll get all kinds of undefined symbol errors.
#if defined HAVE_CBLAS_H
  #include <cblas.h>
#elif defined HAVE_ATLAS_CBLAS_H
  #include <atlas/cblas.h>
#endif
}


namespace nm { namespace math {

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

}} // end of namespace nm::math

#endif // GEMM_H
