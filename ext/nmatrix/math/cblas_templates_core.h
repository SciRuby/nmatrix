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
// == cblas_templates_core.h
//
// This header files is not used by the main nmatrix gem but has to be stored
// in this directory so that it can be shared between nmatrix-atlas and
// nmatrix-lapack.
//

//This is not a normal header file so we don't use an include guard.
//See ext/nmatrix_atlas/math_atlas/cblas_templates_atlas.h for how
//to use.

//Below are the BLAS functions for which we have internal implementations.
//The internal implementations are defined in the ext/nmatrix/math directory
//and are the non-specialized
//forms of the template functions nm::math::whatever().
//They are are called below for non-BLAS
//types in the non-specialized form of the template nm::math::something_else::whatever().
//The specialized forms call the appropriate cblas functions.

//For all functions besides herk, we also define the cblas_whatever() template
//functions below, which just cast
//their arguments to the appropriate types.

//rotg
template <typename DType>
inline void rotg(DType* a, DType* b, DType* c, DType* s) {
  nm::math::rotg(a, b, c, s);
}

template <>
inline void rotg(float* a, float* b, float* c, float* s) {
  cblas_srotg(a, b, c, s);
}

template <>
inline void rotg(double* a, double* b, double* c, double* s) {
  cblas_drotg(a, b, c, s);
}

//Complex versions of rot and rotg are available in the ATLAS (and Intel)
//version of CBLAS, but not part
//of the reference implementation or OpenBLAS, so we omit them here
//and fall back to the generic internal implementation.
//Another options would be to directly call the fortran functions, e.g. ZROTG,
//which for some reason are a part of the standard.
//We can still define complex specializations of these functions in an ATLAS-specific
//header.

template <typename DType>
inline void cblas_rotg(void* a, void* b, void* c, void* s) {
  rotg<DType>(static_cast<DType*>(a), static_cast<DType*>(b), static_cast<DType*>(c), static_cast<DType*>(s));
}

//rot
template <typename DType, typename CSDType>
inline void rot(const int N, DType* X, const int incX, DType* Y, const int incY, const CSDType c, const CSDType s) {
  nm::math::rot<DType,CSDType>(N, X, incX, Y, incY, c, s);
}

template <>
inline void rot(const int N, float* X, const int incX, float* Y, const int incY, const float c, const float s) {
  cblas_srot(N, X, incX, Y, incY, (float)c, (float)s);
}

template <>
inline void rot(const int N, double* X, const int incX, double* Y, const int incY, const double c, const double s) {
  cblas_drot(N, X, incX, Y, incY, c, s);
}

template <typename DType, typename CSDType>
inline void cblas_rot(const int N, void* X, const int incX, void* Y, const int incY, const void* c, const void* s) {
  rot<DType,CSDType>(N, static_cast<DType*>(X), incX, static_cast<DType*>(Y), incY,
                       *static_cast<const CSDType*>(c), *static_cast<const CSDType*>(s));
}

/*
 * Level 1 BLAS routine which sums the absolute values of a vector's contents. If the vector consists of complex values,
 * the routine sums the absolute values of the real and imaginary components as well.
 *
 * So, based on input types, these are the valid return types:
 *    int -> int
 *    float -> float or double
 *    double -> double
 *    complex64 -> float or double
 *    complex128 -> double
 */
template <typename ReturnDType, typename DType>
inline ReturnDType asum(const int N, const DType* X, const int incX) {
  return nm::math::asum<ReturnDType,DType>(N,X,incX);
}


template <>
inline float asum(const int N, const float* X, const int incX) {
  return cblas_sasum(N, X, incX);
}

template <>
inline double asum(const int N, const double* X, const int incX) {
  return cblas_dasum(N, X, incX);
}

template <>
inline float asum(const int N, const Complex64* X, const int incX) {
  return cblas_scasum(N, X, incX);
}

template <>
inline double asum(const int N, const Complex128* X, const int incX) {
  return cblas_dzasum(N, X, incX);
}


template <typename ReturnDType, typename DType>
inline void cblas_asum(const int N, const void* X, const int incX, void* sum) {
  *static_cast<ReturnDType*>( sum ) = asum<ReturnDType, DType>( N, static_cast<const DType*>(X), incX );
}

/*
 * Level 1 BLAS routine which returns the 2-norm of an n-vector x.
 #
 * Based on input types, these are the valid return types:
 *    int -> int
 *    float -> float or double
 *    double -> double
 *    complex64 -> float or double
 *    complex128 -> double
 */
template <typename ReturnDType, typename DType>
inline ReturnDType nrm2(const int N, const DType* X, const int incX) {
  return nm::math::nrm2<ReturnDType,DType>(N, X, incX);
}


template <>
inline float nrm2(const int N, const float* X, const int incX) {
  return cblas_snrm2(N, X, incX);
}

template <>
inline double nrm2(const int N, const double* X, const int incX) {
  return cblas_dnrm2(N, X, incX);
}

template <>
inline float nrm2(const int N, const Complex64* X, const int incX) {
  return cblas_scnrm2(N, X, incX);
}

template <>
inline double nrm2(const int N, const Complex128* X, const int incX) {
  return cblas_dznrm2(N, X, incX);
}

template <typename ReturnDType, typename DType>
inline void cblas_nrm2(const int N, const void* X, const int incX, void* result) {
  *static_cast<ReturnDType*>( result ) = nrm2<ReturnDType, DType>( N, static_cast<const DType*>(X), incX );
}

//imax
template<typename DType>
inline int imax(const int n, const DType *x, const int incx) {
  return nm::math::imax(n, x, incx);
}

template<>
inline int imax(const int n, const float* x, const int incx) {
  return cblas_isamax(n, x, incx);
}

template<>
inline int imax(const int n, const double* x, const int incx) {
  return cblas_idamax(n, x, incx);
}

template<>
inline int imax(const int n, const Complex64* x, const int incx) {
  return cblas_icamax(n, x, incx);
}

template <>
inline int imax(const int n, const Complex128* x, const int incx) {
  return cblas_izamax(n, x, incx);
}

template<typename DType>
inline int cblas_imax(const int n, const void* x, const int incx) {
  return imax<DType>(n, static_cast<const DType*>(x), incx);
}

//scal
template <typename DType>
inline void scal(const int n, const DType scalar, DType* x, const int incx) {
  nm::math::scal(n, scalar, x, incx);
}

template <>
inline void scal(const int n, const float scalar, float* x, const int incx) {
  cblas_sscal(n, scalar, x, incx);
}

template <>
inline void scal(const int n, const double scalar, double* x, const int incx) {
  cblas_dscal(n, scalar, x, incx);
}

template <>
inline void scal(const int n, const Complex64 scalar, Complex64* x, const int incx) {
  cblas_cscal(n, &scalar, x, incx);
}

template <>
inline void scal(const int n, const Complex128 scalar, Complex128* x, const int incx) {
  cblas_zscal(n, &scalar, x, incx);
}

template <typename DType>
inline void cblas_scal(const int n, const void* scalar, void* x, const int incx) {
  scal<DType>(n, *static_cast<const DType*>(scalar), static_cast<DType*>(x), incx);
}

//gemv
template <typename DType>
inline bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const DType* alpha, const DType* A, const int lda,
          const DType* X, const int incX, const DType* beta, DType* Y, const int incY) {
  return nm::math::gemv(Trans, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}

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

template <typename DType>
inline static bool cblas_gemv(const enum CBLAS_TRANSPOSE trans,
                              const int m, const int n,
                              const void* alpha,
                              const void* a, const int lda,
                              const void* x, const int incx,
                              const void* beta,
                              void* y, const int incy)
{
  return gemv<DType>(trans,
                     m, n, static_cast<const DType*>(alpha),
                     static_cast<const DType*>(a), lda,
                     static_cast<const DType*>(x), incx, static_cast<const DType*>(beta),
                     static_cast<DType*>(y), incy);
}

//gemm
template <typename DType>
inline void gemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
                 const DType* alpha, const DType* A, const int lda, const DType* B, const int ldb, const DType* beta, DType* C, const int ldc)
{
  nm::math::gemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
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

template <typename DType>
inline static void cblas_gemm(const enum CBLAS_ORDER order,
                              const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_TRANSPOSE trans_b,
                              int m, int n, int k,
                              void* alpha,
                              void* a, int lda,
                              void* b, int ldb,
                              void* beta,
                              void* c, int ldc)
{
  gemm<DType>(order, trans_a, trans_b, m, n, k, static_cast<DType*>(alpha),
              static_cast<DType*>(a), lda,
              static_cast<DType*>(b), ldb, static_cast<DType*>(beta),
              static_cast<DType*>(c), ldc);
}

//trsm
template <typename DType, typename = typename std::enable_if<!std::is_integral<DType>::value>::type>
inline void trsm(const enum CBLAS_ORDER order,
                 const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_DIAG diag,
                 const int m, const int n, const DType alpha, const DType* a,
                 const int lda, DType* b, const int ldb)
{
  nm::math::trsm(order, side, uplo, trans_a, diag, m, n, alpha, a, lda, b, ldb);
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
  cblas_dtrsm(order, side, uplo, trans_a, diag, m, n, alpha, a, lda, b, ldb);
}


template <>
inline void trsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_DIAG diag,
                 const int m, const int n, const Complex64 alpha, const Complex64* a,
                 const int lda, Complex64* b, const int ldb)
{
  cblas_ctrsm(order, side, uplo, trans_a, diag, m, n, &alpha, a, lda, b, ldb);
}

template <>
inline void trsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_DIAG diag,
                 const int m, const int n, const Complex128 alpha, const Complex128* a,
                 const int lda, Complex128* b, const int ldb)
{
  cblas_ztrsm(order, side, uplo, trans_a, diag, m, n, &alpha, a, lda, b, ldb);
}

template <typename DType>
inline static void cblas_trsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                               const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_DIAG diag,
                               const int m, const int n, const void* alpha, const void* a,
                               const int lda, void* b, const int ldb)
{
  trsm<DType>(order, side, uplo, trans_a, diag, m, n, *static_cast<const DType*>(alpha),
              static_cast<const DType*>(a), lda, static_cast<DType*>(b), ldb);
}

//Below are BLAS functions that we don't have an internal implementation for.
//In this case the non-specialized form just raises an error.

//syrk
template <typename DType>
inline void syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const DType* alpha, const DType* A, const int lda, const DType* beta, DType* C, const int ldc) {
  rb_raise(rb_eNotImpError, "syrk not yet implemented for non-BLAS dtypes");
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

template <typename DType>
inline static void cblas_syrk(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const enum CBLAS_TRANSPOSE trans,
                              const int n, const int k, const void* alpha,
                              const void* A, const int lda, const void* beta, void* C, const int ldc)
{
  syrk<DType>(order, uplo, trans, n, k, static_cast<const DType*>(alpha),
              static_cast<const DType*>(A), lda, static_cast<const DType*>(beta), static_cast<DType*>(C), ldc);
}

//herk
template <typename DType>
inline void herk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const DType* alpha, const DType* A, const int lda, const DType* beta, DType* C, const int ldc) {
  rb_raise(rb_eNotImpError, "herk not yet implemented for non-BLAS dtypes");
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

//trmm
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

template <typename DType>
inline static void cblas_trmm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                              const enum CBLAS_TRANSPOSE ta, const enum CBLAS_DIAG diag, const int m, const int n, const void* alpha,
                              const void* A, const int lda, void* B, const int ldb)
{
  trmm<DType>(order, side, uplo, ta, diag, m, n, static_cast<const DType*>(alpha),
              static_cast<const DType*>(A), lda, static_cast<DType*>(B), ldb);
}
