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
// == lapacke_templates.h
//
// Templated functions for calling LAPACKE functions directly.
//

#ifndef LAPACKE_TEMPLATES_H
#define LAPACKE_TEMPLATES_H

namespace nm { namespace math { namespace lapacke {

//getrf
template <typename DType>
inline int getrf(const enum CBLAS_ORDER order, const int m, const int n, DType* a, const int lda, int* ipiv) {
  //We don't want to call the internal implementation since the the CLAPACK interface is slightly different than the LAPACKE.
  rb_raise(rb_eNotImpError, "lapacke_getrf not implemented for non_BLAS dtypes. Try clapack_getrf instead.");
  return 0;
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

template <typename DType>
inline int lapacke_getrf(const enum CBLAS_ORDER order, const int m, const int n, void* a, const int lda, int* ipiv) {
  return getrf<DType>(order, m, n, static_cast<DType*>(a), lda, ipiv);
}

//getri
template <typename DType>
inline int getri(const enum CBLAS_ORDER order, const int n, DType* a, const int lda, const int* ipiv) {
  rb_raise(rb_eNotImpError, "getri not yet implemented for non-BLAS dtypes");
  return 0;
}

template <>
inline int getri(const enum CBLAS_ORDER order, const int n, float* a, const int lda, const int* ipiv) {
  return LAPACKE_sgetri(order, n, a, lda, ipiv);
}

template <>
inline int getri(const enum CBLAS_ORDER order, const int n, double* a, const int lda, const int* ipiv) {
  return LAPACKE_dgetri(order, n, a, lda, ipiv);
}

template <>
inline int getri(const enum CBLAS_ORDER order, const int n, Complex64* a, const int lda, const int* ipiv) {
  return LAPACKE_cgetri(order, n, a, lda, ipiv);
}

template <>
inline int getri(const enum CBLAS_ORDER order, const int n, Complex128* a, const int lda, const int* ipiv) {
  return LAPACKE_zgetri(order, n, a, lda, ipiv);
}

template <typename DType>
inline int lapacke_getri(const enum CBLAS_ORDER order, const int n, void* a, const int lda, const int* ipiv) {
  return getri<DType>(order, n, static_cast<DType*>(a), lda, ipiv);
}

//getrs
template <typename DType>
inline int getrs(const enum CBLAS_ORDER Order, char Trans, const int N, const int NRHS, const DType* A,
           const int lda, const int* ipiv, DType* B, const int ldb)
{
  rb_raise(rb_eNotImpError, "lapacke_getrs not implemented for non_BLAS dtypes. Try clapack_getrs instead.");
  return 0;
}

template <>
inline int getrs(const enum CBLAS_ORDER Order, char Trans, const int N, const int NRHS, const float* A,
           const int lda, const int* ipiv, float* B, const int ldb)
{
  return LAPACKE_sgetrs(Order, Trans, N, NRHS, A, lda, ipiv, B, ldb);
}

template <>
inline int getrs(const enum CBLAS_ORDER Order, char Trans, const int N, const int NRHS, const double* A,
           const int lda, const int* ipiv, double* B, const int ldb)
{
  return LAPACKE_dgetrs(Order, Trans, N, NRHS, A, lda, ipiv, B, ldb);
}

template <>
inline int getrs(const enum CBLAS_ORDER Order, char Trans, const int N, const int NRHS, const Complex64* A,
           const int lda, const int* ipiv, Complex64* B, const int ldb)
{
  return LAPACKE_cgetrs(Order, Trans, N, NRHS, A, lda, ipiv, B, ldb);
}

template <>
inline int getrs(const enum CBLAS_ORDER Order, char Trans, const int N, const int NRHS, const Complex128* A,
           const int lda, const int* ipiv, Complex128* B, const int ldb)
{
  return LAPACKE_zgetrs(Order, Trans, N, NRHS, A, lda, ipiv, B, ldb);
}

template <typename DType>
inline int lapacke_getrs(const enum CBLAS_ORDER order, char trans, const int n, const int nrhs,
                         const void* a, const int lda, const int* ipiv, void* b, const int ldb) {
  return getrs<DType>(order, trans, n, nrhs, static_cast<const DType*>(a), lda, ipiv, static_cast<DType*>(b), ldb);
}

//potrf
template <typename DType>
inline int potrf(const enum CBLAS_ORDER order, char uplo, const int N, DType* A, const int lda) {
  rb_raise(rb_eNotImpError, "not implemented for non-BLAS dtypes");
  return 0;
}

template <>
inline int potrf(const enum CBLAS_ORDER order, char uplo, const int N, float* A, const int lda) {
  return LAPACKE_spotrf(order, uplo, N, A, lda);
}

template <>
inline int potrf(const enum CBLAS_ORDER order, char uplo, const int N, double* A, const int lda) {
  return LAPACKE_dpotrf(order, uplo, N, A, lda);
}

template <>
inline int potrf(const enum CBLAS_ORDER order, char uplo, const int N, Complex64* A, const int lda) {
  return LAPACKE_cpotrf(order, uplo, N, A, lda);
}

template <>
inline int potrf(const enum CBLAS_ORDER order, char uplo, const int N, Complex128* A, const int lda) {
  return LAPACKE_zpotrf(order, uplo, N, A, lda);
}

template <typename DType>
inline int lapacke_potrf(const enum CBLAS_ORDER order, char uplo, const int n, void* a, const int lda) {
  return potrf<DType>(order, uplo, n, static_cast<DType*>(a), lda);
}

//potrs
template <typename DType>
inline int potrs(const enum CBLAS_ORDER Order, char Uplo, const int N, const int NRHS, const DType* A,
           const int lda, DType* B, const int ldb)
{
  rb_raise(rb_eNotImpError, "not implemented for non-BLAS dtypes");
  return 0;
}

template <>
inline int potrs<float> (const enum CBLAS_ORDER Order, char Uplo, const int N, const int NRHS, const float* A,
           const int lda, float* B, const int ldb)
{
  return LAPACKE_spotrs(Order, Uplo, N, NRHS, A, lda, B, ldb);
}

template <>
inline int potrs<double>(const enum CBLAS_ORDER Order, char Uplo, const int N, const int NRHS, const double* A,
           const int lda, double* B, const int ldb)
{
  return LAPACKE_dpotrs(Order, Uplo, N, NRHS, A, lda, B, ldb);
}

template <>
inline int potrs<Complex64>(const enum CBLAS_ORDER Order, char Uplo, const int N, const int NRHS, const Complex64* A,
           const int lda, Complex64* B, const int ldb)
{
  return LAPACKE_cpotrs(Order, Uplo, N, NRHS, A, lda, B, ldb);
}

template <>
inline int potrs<Complex128>(const enum CBLAS_ORDER Order, char Uplo, const int N, const int NRHS, const Complex128* A,
           const int lda, Complex128* B, const int ldb)
{
  return LAPACKE_zpotrs(Order, Uplo, N, NRHS, A, lda, B, ldb);
}

template <typename DType>
inline int lapacke_potrs(const enum CBLAS_ORDER order, char uplo, const int n, const int nrhs,
                         const void* a, const int lda, void* b, const int ldb) {
  return potrs<DType>(order, uplo, n, nrhs, static_cast<const DType*>(a), lda, static_cast<DType*>(b), ldb);
}

//potri
template <typename DType>
inline int potri(const enum CBLAS_ORDER order, char uplo, const int n, DType* a, const int lda) {
  rb_raise(rb_eNotImpError, "potri not yet implemented for non-BLAS dtypes");
  return 0;
}

template <>
inline int potri(const enum CBLAS_ORDER order, char uplo, const int n, float* a, const int lda) {
  return LAPACKE_spotri(order, uplo, n, a, lda);
}

template <>
inline int potri(const enum CBLAS_ORDER order, char uplo, const int n, double* a, const int lda) {
  return LAPACKE_dpotri(order, uplo, n, a, lda);
}

template <>
inline int potri(const enum CBLAS_ORDER order, char uplo, const int n, Complex64* a, const int lda) {
  return LAPACKE_cpotri(order, uplo, n, a, lda);
}

template <>
inline int potri(const enum CBLAS_ORDER order, char uplo, const int n, Complex128* a, const int lda) {
  return LAPACKE_zpotri(order, uplo, n, a, lda);
}

template <typename DType>
inline int lapacke_potri(const enum CBLAS_ORDER order, char uplo, const int n, void* a, const int lda) {
  return potri<DType>(order, uplo, n, static_cast<DType*>(a), lda);
}

//gesvd
template <typename DType, typename CType>
inline int gesvd(int matrix_layout, char jobu, char jobvt, int m, int n, DType* a, int lda, CType* s, DType* u, int ldu, DType* vt, int ldvt, CType* superb) {
  rb_raise(rb_eNotImpError, "gesvd not yet implemented for non-BLAS dtypes");
  return 0;
}

template <>
inline int gesvd<float, float>(int matrix_layout, char jobu, char jobvt, int m, int n, float* a, int lda, float* s, float* u, int ldu, float* vt, int ldvt, float* superb) {
  return LAPACKE_sgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}

template <>
inline int gesvd<double, double>(int matrix_layout, char jobu, char jobvt, int m, int n, double* a, int lda, double* s, double* u, int ldu, double* vt, int ldvt, double* superb) {
  return LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}

template <>
inline int gesvd<nm::Complex64, float>(int matrix_layout, char jobu, char jobvt, int m, int n, nm::Complex64* a, int lda, float* s, nm::Complex64* u, int ldu, nm::Complex64* vt, int ldvt, float* superb) {
  return LAPACKE_cgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}

template <>
inline int gesvd<nm::Complex128, double>(int matrix_layout, char jobu, char jobvt, int m, int n, nm::Complex128* a, int lda, double* s, nm::Complex128* u, int ldu, nm::Complex128* vt, int ldvt, double* superb) {
  return LAPACKE_zgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}

template <typename DType, typename CType>
inline int lapacke_gesvd(int matrix_layout, char jobu, char jobvt, int m, int n, void* a, int lda, void* s, void* u, int ldu, void* vt, int ldvt, void* superb) {
  return gesvd<DType,CType>(matrix_layout, jobu, jobvt, m, n, static_cast<DType*>(a), lda, static_cast<CType*>(s), static_cast<DType*>(u), ldu, static_cast<DType*>(vt), ldvt, static_cast<CType*>(superb));
}

//gesdd
template <typename DType, typename CType>
inline int gesdd(int matrix_layout, char jobz, int m, int n, DType* a, int lda, CType* s, DType* u, int ldu, DType* vt, int ldvt) {
  rb_raise(rb_eNotImpError, "gesdd not yet implemented for non-BLAS dtypes");
  return 0;
}

template <>
inline int gesdd<float, float>(int matrix_layout, char jobz, int m, int n, float* a, int lda, float* s, float* u, int ldu, float* vt, int ldvt) {
  return LAPACKE_sgesdd(matrix_layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
}

template <>
inline int gesdd<double, double>(int matrix_layout, char jobz, int m, int n, double* a, int lda, double* s, double* u, int ldu, double* vt, int ldvt) {
  return LAPACKE_dgesdd(matrix_layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
}

template <>
inline int gesdd<nm::Complex64, float>(int matrix_layout, char jobz, int m, int n, nm::Complex64* a, int lda, float* s, nm::Complex64* u, int ldu, nm::Complex64* vt, int ldvt) {
  return LAPACKE_cgesdd(matrix_layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
}

template <>
inline int gesdd<nm::Complex128, double>(int matrix_layout, char jobz, int m, int n, nm::Complex128* a, int lda, double* s, nm::Complex128* u, int ldu, nm::Complex128* vt, int ldvt) {
  return LAPACKE_zgesdd(matrix_layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
}

template <typename DType, typename CType>
inline int lapacke_gesdd(int matrix_layout, char jobz, int m, int n, void* a, int lda, void* s, void* u, int ldu, void* vt, int ldvt) {
  return gesdd<DType,CType>(matrix_layout, jobz, m, n, static_cast<DType*>(a), lda, static_cast<CType*>(s), static_cast<DType*>(u), ldu, static_cast<DType*>(vt), ldvt);
}

//geev
//This one is a little tricky. The signature is different for the complex
//versions than for the real ones. This is because real matrices can have
//complex eigenvalues. For the complex types, the eigenvalues are just
//returned in argument that's a complex array, but for real types the real
//parts of the eigenvalues are returned
//in one (array) argument, and the complex parts in a separate argument.
//The solution is that the template takes an vi argument, but it is just
//ignored in the specializations for complex types.

template <typename DType>
inline int geev(int matrix_layout, char jobvl, char jobvr, int n, DType* a, int lda, DType* w, DType* wi, DType* vl, int ldvl, DType* vr, int ldvr) {
  rb_raise(rb_eNotImpError, "not yet implemented for non-BLAS dtypes");
  return -1;
}

template <>
inline int geev(int matrix_layout, char jobvl, char jobvr, int n, float* a, int lda, float* w, float* wi, float* vl, int ldvl, float* vr, int ldvr) {
  return LAPACKE_sgeev(matrix_layout, jobvl, jobvr, n, a, lda, w, wi, vl, ldvl, vr, ldvr);
}

template <>
inline int geev(int matrix_layout, char jobvl, char jobvr, int n, double* a, int lda, double* w, double* wi, double* vl, int ldvl, double* vr, int ldvr) {
  return LAPACKE_dgeev(matrix_layout, jobvl, jobvr, n, a, lda, w, wi, vl, ldvl, vr, ldvr);
}

template <>
inline int geev(int matrix_layout, char jobvl, char jobvr, int n, Complex64* a, int lda, Complex64* w, Complex64* wi, Complex64* vl, int ldvl, Complex64* vr, int ldvr) {
  return LAPACKE_cgeev(matrix_layout, jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr);
}

template <>
inline int geev(int matrix_layout, char jobvl, char jobvr, int n, Complex128* a, int lda, Complex128* w, Complex128* wi, Complex128* vl, int ldvl, Complex128* vr, int ldvr) {
  return LAPACKE_zgeev(matrix_layout, jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr);
}

template <typename DType>
inline int lapacke_geev(int matrix_layout, char jobvl, char jobvr, int n, void* a, int lda, void* w, void* wi, void* vl, int ldvl, void* vr, int ldvr) {
  return geev<DType>(matrix_layout, jobvl, jobvr, n, static_cast<DType*>(a), lda, static_cast<DType*>(w), static_cast<DType*>(wi), static_cast<DType*>(vl), ldvl, static_cast<DType*>(vr), ldvr);
}

}}}

#endif
