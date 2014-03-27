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
// == geev.h
//
// Header file for interface with LAPACK's xGEEV functions.
//

#ifndef GEEV_H
# define GEEV_H

extern "C" {
  void sgeev_(char* jobvl, char* jobvr, int* n, float* a,          int* lda, float* wr,  float* wi,  float* vl,          int* ldvl, float* vr,          int* ldvr, float* work,          int* lwork,                int* info);
  void dgeev_(char* jobvl, char* jobvr, int* n, double* a,         int* lda, double* wr, double* wi, double* vl,         int* ldvl, double* vr,         int* ldvr, double* work,         int* lwork,                int* info);
  void cgeev_(char* jobvl, char* jobvr, int* n, nm::Complex64* a,  int* lda, nm::Complex64* w,       nm::Complex64* vl,  int* ldvl, nm::Complex64* vr,  int* ldvr, nm::Complex64* work,  int* lwork, float* rwork,  int* info);
  void zgeev_(char* jobvl, char* jobvr, int* n, nm::Complex128* a, int* lda, nm::Complex128* w,      nm::Complex128* vl, int* ldvl, nm::Complex128* vr, int* ldvr, nm::Complex128* work, int* lwork, double* rwork, int* info);
}

namespace nm { namespace math {

template <typename DType, typename CType>                         // wr
inline int geev(char jobvl, char jobvr, int n, DType* a, int lda, DType* w, DType* wi, DType* vl, int ldvl, DType* vr, int ldvr, DType* work, int lwork, CType* rwork) {
  rb_raise(rb_eNotImpError, "not yet implemented for non-BLAS dtypes");
  return -1;
}

template <>
inline int geev(char jobvl, char jobvr, int n, float* a, int lda, float* w, float* wi, float* vl, int ldvl, float* vr, int ldvr, float* work, int lwork, float* rwork) {
  int info;
  sgeev_(&jobvl, &jobvr, &n, a, &lda, w, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
  return info;
}

template <>
inline int geev(char jobvl, char jobvr, int n, double* a, int lda, double* w, double* wi, double* vl, int ldvl, double* vr, int ldvr, double* work, int lwork, double* rwork) {
  int info;
  dgeev_(&jobvl, &jobvr, &n, a, &lda, w, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
  return info;
}

template <>
inline int geev(char jobvl, char jobvr, int n, Complex64* a, int lda, Complex64* w, Complex64* wi, Complex64* vl, int ldvl, Complex64* vr, int ldvr, Complex64* work, int lwork, float* rwork) {
  int info;
  cgeev_(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
  return info;
}

template <>
inline int geev(char jobvl, char jobvr, int n, Complex128* a, int lda, Complex128* w, Complex128* wi, Complex128* vl, int ldvl, Complex128* vr, int ldvr, Complex128* work, int lwork, double* rwork) {
  int info;
  zgeev_(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
  return info;
}

template <typename DType, typename CType>
inline int lapack_geev(char jobvl, char jobvr, int n, void* a, int lda, void* w, void* wi, void* vl, int ldvl, void* vr, int ldvr, void* work, int lwork, void* rwork) {
  return geev<DType,CType>(jobvl, jobvr, n, reinterpret_cast<DType*>(a), lda, reinterpret_cast<DType*>(w), reinterpret_cast<DType*>(wi), reinterpret_cast<DType*>(vl), ldvl, reinterpret_cast<DType*>(vr), ldvr, reinterpret_cast<DType*>(work), lwork, reinterpret_cast<CType*>(rwork));
}

}} // end nm::math

#endif // GEEV_H
