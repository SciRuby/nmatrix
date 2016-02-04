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
// == clapack_templates.h
//
// Collection of functions used to call ATLAS CLAPACK functions
// directly.
//

#ifndef CLAPACK_TEMPLATES_H
#define CLAPACK_TEMPLATES_H

//needed to get access to internal implementations
#include "math/getrf.h"
#include "math/getrs.h"

namespace nm { namespace math { namespace atlas {
//The first group of functions are those for which we have internal implementations.
//The internal implementations are defined in the ext/nmatrix/math directory
//and are the non-specialized
//forms of the template functions nm::math::whatever().
//They are are called below for non-BLAS
//types in the non-specialized form of the template nm::math::atlas::whatever().
//The specialized forms call the appropriate clapack functions.

//We also define the clapack_whatever() template
//functions below, which just cast
//their arguments to the appropriate types.


//getrf
template <typename DType>
inline int getrf(const enum CBLAS_ORDER order, const int m, const int n, DType* a, const int lda, int* ipiv) {
  return nm::math::getrf<DType>(order, m, n, a, lda, ipiv);
}

//Apparently CLAPACK isn't available on OS X, so we only define these
//specializations if available,
#if defined (HAVE_CLAPACK_H) || defined (HAVE_ATLAS_CLAPACK_H)
template <>
inline int getrf(const enum CBLAS_ORDER order, const int m, const int n, float* a, const int lda, int* ipiv) {
  return clapack_sgetrf(order, m, n, a, lda, ipiv);
}

template <>
inline int getrf(const enum CBLAS_ORDER order, const int m, const int n, double* a, const int lda, int* ipiv) {
  return clapack_dgetrf(order, m, n, a, lda, ipiv);
}

template <>
inline int getrf(const enum CBLAS_ORDER order, const int m, const int n, Complex64* a, const int lda, int* ipiv) {
  return clapack_cgetrf(order, m, n, a, lda, ipiv);
}

template <>
inline int getrf(const enum CBLAS_ORDER order, const int m, const int n, Complex128* a, const int lda, int* ipiv) {
  return clapack_zgetrf(order, m, n, a, lda, ipiv);
}
#endif

template <typename DType>
inline int clapack_getrf(const enum CBLAS_ORDER order, const int m, const int n, void* a, const int lda, int* ipiv) {
  return getrf<DType>(order, m, n, static_cast<DType*>(a), lda, ipiv);
}

//getrs
/*
 * Solves a system of linear equations A*X = B with a general NxN matrix A using the LU factorization computed by GETRF.
 *
 * From ATLAS 3.8.0.
 */
template <typename DType>
inline int getrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans, const int N, const int NRHS, const DType* A,
           const int lda, const int* ipiv, DType* B, const int ldb)
{
  return nm::math::getrs<DType>(Order, Trans, N, NRHS, A, lda, ipiv, B, ldb);
}

#if defined (HAVE_CLAPACK_H) || defined (HAVE_ATLAS_CLAPACK_H)
template <>
inline int getrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans, const int N, const int NRHS, const float* A,
           const int lda, const int* ipiv, float* B, const int ldb)
{
  return clapack_sgetrs(Order, Trans, N, NRHS, A, lda, ipiv, B, ldb);
}

template <>
inline int getrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans, const int N, const int NRHS, const double* A,
           const int lda, const int* ipiv, double* B, const int ldb)
{
  return clapack_dgetrs(Order, Trans, N, NRHS, A, lda, ipiv, B, ldb);
}

template <>
inline int getrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans, const int N, const int NRHS, const Complex64* A,
           const int lda, const int* ipiv, Complex64* B, const int ldb)
{
  return clapack_cgetrs(Order, Trans, N, NRHS, A, lda, ipiv, static_cast<void*>(B), ldb);
}

template <>
inline int getrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans, const int N, const int NRHS, const Complex128* A,
           const int lda, const int* ipiv, Complex128* B, const int ldb)
{
  return clapack_zgetrs(Order, Trans, N, NRHS, A, lda, ipiv, static_cast<void*>(B), ldb);
}
#endif

template <typename DType>
inline int clapack_getrs(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE trans, const int n, const int nrhs,
                         const void* a, const int lda, const int* ipiv, void* b, const int ldb) {
  return getrs<DType>(order, trans, n, nrhs, static_cast<const DType*>(a), lda, ipiv, static_cast<DType*>(b), ldb);
}


//Functions without internal implementations below:

//getri
template <typename DType>
inline int getri(const enum CBLAS_ORDER order, const int n, DType* a, const int lda, const int* ipiv) {
  rb_raise(rb_eNotImpError, "getri not yet implemented for non-BLAS dtypes");
  return 0;
}

#if defined (HAVE_CLAPACK_H) || defined (HAVE_ATLAS_CLAPACK_H)
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
  return clapack_cgetri(order, n, a, lda, ipiv);
}

template <>
inline int getri(const enum CBLAS_ORDER order, const int n, Complex128* a, const int lda, const int* ipiv) {
  return clapack_zgetri(order, n, a, lda, ipiv);
}
#endif

template <typename DType>
inline int clapack_getri(const enum CBLAS_ORDER order, const int n, void* a, const int lda, const int* ipiv) {
  return getri<DType>(order, n, static_cast<DType*>(a), lda, ipiv);
}

//potrf
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
#if defined HAVE_CLAPACK_H || defined HAVE_ATLAS_CLAPACK_H
  rb_raise(rb_eNotImpError, "not yet implemented for non-BLAS dtypes");
#else
  rb_raise(rb_eNotImpError, "only CLAPACK version implemented thus far");
#endif
  return 0;
}

#if defined HAVE_CLAPACK_H || defined HAVE_ATLAS_CLAPACK_H
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
  return clapack_cpotrf(order, uplo, N, A, lda);
}

template <>
inline int potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, Complex128* A, const int lda) {
  return clapack_zpotrf(order, uplo, N, A, lda);
}
#endif

template <typename DType>
inline int clapack_potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, void* a, const int lda) {
  return potrf<DType>(order, uplo, n, static_cast<DType*>(a), lda);
}

//potri
template <typename DType>
inline int potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, DType* a, const int lda) {
  rb_raise(rb_eNotImpError, "potri not yet implemented for non-BLAS dtypes");
  return 0;
}


#if defined HAVE_CLAPACK_H || defined HAVE_ATLAS_CLAPACK_H
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
  return clapack_cpotri(order, uplo, n, a, lda);
}

template <>
inline int potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, Complex128* a, const int lda) {
  return clapack_zpotri(order, uplo, n, a, lda);
}
#endif

template <typename DType>
inline int clapack_potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, void* a, const int lda) {
  return potri<DType>(order, uplo, n, static_cast<DType*>(a), lda);
}

//potrs
/*
 * Solves a system of linear equations A*X = B with a symmetric positive definite matrix A using the Cholesky factorization computed by POTRF.
 */
template <typename DType>
inline int potrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const int N, const int NRHS, const DType* A,
           const int lda, DType* B, const int ldb)
{
#if defined HAVE_CLAPACK_H || defined HAVE_ATLAS_CLAPACK_H
  rb_raise(rb_eNotImpError, "not yet implemented for non-BLAS dtypes");
#else
  rb_raise(rb_eNotImpError, "only CLAPACK version implemented thus far");
#endif
}

#if defined (HAVE_CLAPACK_H) || defined (HAVE_ATLAS_CLAPACK_H)
template <>
inline int potrs<float> (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const int N, const int NRHS, const float* A,
           const int lda, float* B, const int ldb)
{
  return clapack_spotrs(Order, Uplo, N, NRHS, A, lda, B, ldb);
}

template <>
inline int potrs<double>(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const int N, const int NRHS, const double* A,
           const int lda, double* B, const int ldb)
{
  return clapack_dpotrs(Order, Uplo, N, NRHS, A, lda, B, ldb);
}

template <>
inline int potrs<Complex64>(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const int N, const int NRHS, const Complex64* A,
           const int lda, Complex64* B, const int ldb)
{
  return clapack_cpotrs(Order, Uplo, N, NRHS, A, lda, static_cast<void *>(B), ldb);
}

template <>
inline int potrs<Complex128>(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const int N, const int NRHS, const Complex128* A,
           const int lda, Complex128* B, const int ldb)
{
  return clapack_zpotrs(Order, Uplo, N, NRHS, A, lda, static_cast<void *>(B), ldb);
}
#endif

template <typename DType>
inline int clapack_potrs(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, const int nrhs,
                         const void* a, const int lda, void* b, const int ldb) {
  return potrs<DType>(order, uplo, n, nrhs, static_cast<const DType*>(a), lda, static_cast<DType*>(b), ldb);
}

}}}

#endif
