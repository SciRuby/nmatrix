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
// == math_lapacke.cpp
//
// Ruby-exposed CBLAS and LAPACK functions that call BLAS
// and LAPACKE functions.
//

#include "data/data.h"

#include "lapacke_nmatrix.h"

#include "math_lapacke/cblas_local.h"

#include "math/util.h"

#include "math_lapacke/cblas_templates_lapacke.h"

#include "math_lapacke/lapacke_templates.h"


/*
 * Forward Declarations
 */

extern "C" {
  /* BLAS Level 1. */
  static VALUE nm_lapacke_cblas_scal(VALUE self, VALUE n, VALUE scale, VALUE vector, VALUE incx);
  static VALUE nm_lapacke_cblas_nrm2(VALUE self, VALUE n, VALUE x, VALUE incx);
  static VALUE nm_lapacke_cblas_asum(VALUE self, VALUE n, VALUE x, VALUE incx);
  static VALUE nm_lapacke_cblas_rot(VALUE self, VALUE n, VALUE x, VALUE incx, VALUE y, VALUE incy, VALUE c, VALUE s);
  static VALUE nm_lapacke_cblas_rotg(VALUE self, VALUE ab);
  static VALUE nm_lapacke_cblas_imax(VALUE self, VALUE n, VALUE x, VALUE incx);

  /* BLAS Level 2. */
  static VALUE nm_lapacke_cblas_gemv(VALUE self, VALUE trans_a, VALUE m, VALUE n, VALUE vAlpha, VALUE a, VALUE lda,
                             VALUE x, VALUE incx, VALUE vBeta, VALUE y, VALUE incy);

  /* BLAS Level 3. */
  static VALUE nm_lapacke_cblas_gemm(VALUE self, VALUE order, VALUE trans_a, VALUE trans_b, VALUE m, VALUE n, VALUE k, VALUE vAlpha,
                             VALUE a, VALUE lda, VALUE b, VALUE ldb, VALUE vBeta, VALUE c, VALUE ldc);
  static VALUE nm_lapacke_cblas_trsm(VALUE self, VALUE order, VALUE side, VALUE uplo, VALUE trans_a, VALUE diag, VALUE m, VALUE n,
                             VALUE vAlpha, VALUE a, VALUE lda, VALUE b, VALUE ldb);
  static VALUE nm_lapacke_cblas_trmm(VALUE self, VALUE order, VALUE side, VALUE uplo, VALUE trans_a, VALUE diag, VALUE m, VALUE n,
                             VALUE alpha, VALUE a, VALUE lda, VALUE b, VALUE ldb);
  static VALUE nm_lapacke_cblas_herk(VALUE self, VALUE order, VALUE uplo, VALUE trans, VALUE n, VALUE k, VALUE alpha, VALUE a,
                             VALUE lda, VALUE beta, VALUE c, VALUE ldc);
  static VALUE nm_lapacke_cblas_syrk(VALUE self, VALUE order, VALUE uplo, VALUE trans, VALUE n, VALUE k, VALUE alpha, VALUE a,
                             VALUE lda, VALUE beta, VALUE c, VALUE ldc);

  /* LAPACK. */
  static VALUE nm_lapacke_lapacke_getrf(VALUE self, VALUE order, VALUE m, VALUE n, VALUE a, VALUE lda);
  static VALUE nm_lapacke_lapacke_getrs(VALUE self, VALUE order, VALUE trans, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE ipiv, VALUE b, VALUE ldb);
  static VALUE nm_lapacke_lapacke_getri(VALUE self, VALUE order, VALUE n, VALUE a, VALUE lda, VALUE ipiv);
  static VALUE nm_lapacke_lapacke_potrf(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda);
  static VALUE nm_lapacke_lapacke_potrs(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE b, VALUE ldb);
  static VALUE nm_lapacke_lapacke_potri(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda);

  static VALUE nm_lapacke_lapacke_geqrf(VALUE self, VALUE order, VALUE m, VALUE n, VALUE a, VALUE lda, VALUE tau);
  static VALUE nm_lapacke_lapacke_ormqr(VALUE self, VALUE order, VALUE side, VALUE trans, VALUE m, VALUE n, VALUE k, VALUE a, VALUE lda, VALUE tau, VALUE c, VALUE ldc);
  static VALUE nm_lapacke_lapacke_unmqr(VALUE self, VALUE order, VALUE side, VALUE trans, VALUE m, VALUE n, VALUE k, VALUE a, VALUE lda, VALUE tau, VALUE c, VALUE ldc);


  static VALUE nm_lapacke_lapacke_gesvd(VALUE self, VALUE order, VALUE jobu, VALUE jobvt, VALUE m, VALUE n, VALUE a, VALUE lda, VALUE s, VALUE u, VALUE ldu, VALUE vt, VALUE ldvt, VALUE superb);
  static VALUE nm_lapacke_lapacke_gesdd(VALUE self, VALUE order, VALUE jobz, VALUE m, VALUE n, VALUE a, VALUE lda, VALUE s, VALUE u, VALUE ldu, VALUE vt, VALUE ldvt);
  static VALUE nm_lapacke_lapacke_geev(VALUE self, VALUE order, VALUE jobvl, VALUE jobvr, VALUE n, VALUE a, VALUE lda, VALUE w, VALUE wi, VALUE vl, VALUE ldvl, VALUE vr, VALUE ldvr);
}

extern "C" {

///////////////////
// Ruby Bindings //
///////////////////

void nm_math_init_lapack() {

  VALUE cNMatrix_LAPACKE = rb_define_module_under(cNMatrix, "LAPACKE");

  VALUE cNMatrix_LAPACKE_LAPACK = rb_define_module_under(cNMatrix_LAPACKE, "LAPACK");
  VALUE cNMatrix_LAPACKE_BLAS = rb_define_module_under(cNMatrix_LAPACKE, "BLAS");

  //BLAS Level 1
  rb_define_singleton_method(cNMatrix_LAPACKE_BLAS, "cblas_scal", (METHOD)nm_lapacke_cblas_scal, 4);
  rb_define_singleton_method(cNMatrix_LAPACKE_BLAS, "cblas_nrm2", (METHOD)nm_lapacke_cblas_nrm2, 3);
  rb_define_singleton_method(cNMatrix_LAPACKE_BLAS, "cblas_asum", (METHOD)nm_lapacke_cblas_asum, 3);
  rb_define_singleton_method(cNMatrix_LAPACKE_BLAS, "cblas_rot",  (METHOD)nm_lapacke_cblas_rot,  7);
  rb_define_singleton_method(cNMatrix_LAPACKE_BLAS, "cblas_rotg", (METHOD)nm_lapacke_cblas_rotg, 1);
  rb_define_singleton_method(cNMatrix_LAPACKE_BLAS, "cblas_imax", (METHOD)nm_lapacke_cblas_imax, 3);

  //BLAS Level 2
  rb_define_singleton_method(cNMatrix_LAPACKE_BLAS, "cblas_gemv", (METHOD)nm_lapacke_cblas_gemv, 11);

  //BLAS Level 3
  rb_define_singleton_method(cNMatrix_LAPACKE_BLAS, "cblas_gemm", (METHOD)nm_lapacke_cblas_gemm, 14);
  rb_define_singleton_method(cNMatrix_LAPACKE_BLAS, "cblas_trsm", (METHOD)nm_lapacke_cblas_trsm, 12);
  rb_define_singleton_method(cNMatrix_LAPACKE_BLAS, "cblas_trmm", (METHOD)nm_lapacke_cblas_trmm, 12);
  rb_define_singleton_method(cNMatrix_LAPACKE_BLAS, "cblas_syrk", (METHOD)nm_lapacke_cblas_syrk, 11);
  rb_define_singleton_method(cNMatrix_LAPACKE_BLAS, "cblas_herk", (METHOD)nm_lapacke_cblas_herk, 11);

  /* LAPACK Functions */
  rb_define_singleton_method(cNMatrix_LAPACKE_LAPACK, "lapacke_getrf", (METHOD)nm_lapacke_lapacke_getrf, 5);
  rb_define_singleton_method(cNMatrix_LAPACKE_LAPACK, "lapacke_getrs", (METHOD)nm_lapacke_lapacke_getrs, 9);
  rb_define_singleton_method(cNMatrix_LAPACKE_LAPACK, "lapacke_getri", (METHOD)nm_lapacke_lapacke_getri, 5);
  rb_define_singleton_method(cNMatrix_LAPACKE_LAPACK, "lapacke_potrf", (METHOD)nm_lapacke_lapacke_potrf, 5);
  rb_define_singleton_method(cNMatrix_LAPACKE_LAPACK, "lapacke_potrs", (METHOD)nm_lapacke_lapacke_potrs, 8);
  rb_define_singleton_method(cNMatrix_LAPACKE_LAPACK, "lapacke_potri", (METHOD)nm_lapacke_lapacke_potri, 5);

  rb_define_singleton_method(cNMatrix_LAPACKE_LAPACK, "lapacke_geqrf", (METHOD)nm_lapacke_lapacke_geqrf, 6);
  rb_define_singleton_method(cNMatrix_LAPACKE_LAPACK, "lapacke_ormqr", (METHOD)nm_lapacke_lapacke_ormqr, 11);
  rb_define_singleton_method(cNMatrix_LAPACKE_LAPACK, "lapacke_unmqr", (METHOD)nm_lapacke_lapacke_unmqr, 11);

  rb_define_singleton_method(cNMatrix_LAPACKE_LAPACK, "lapacke_gesvd", (METHOD)nm_lapacke_lapacke_gesvd, 13);
  rb_define_singleton_method(cNMatrix_LAPACKE_LAPACK, "lapacke_gesdd", (METHOD)nm_lapacke_lapacke_gesdd, 11);
  rb_define_singleton_method(cNMatrix_LAPACKE_LAPACK, "lapacke_geev", (METHOD)nm_lapacke_lapacke_geev, 12);
}

/*
 * call-seq:
 *     NMatrix::BLAS.cblas_scal(n, alpha, vector, inc) -> NMatrix
 *
 * BLAS level 1 function +scal+. Works with all dtypes.
 *
 * Scale +vector+ in-place by +alpha+ and also return it. The operation is as
 * follows:
 *  x <- alpha * x
 *
 * - +n+ -> Number of elements of +vector+.
 * - +alpha+ -> Scalar value used in the operation.
 * - +vector+ -> NMatrix of shape [n,1] or [1,n]. Modified in-place.
 * - +inc+ -> Increment used in the scaling function. Should generally be 1.
 */
static VALUE nm_lapacke_cblas_scal(VALUE self, VALUE n, VALUE alpha, VALUE vector, VALUE incx) {
  nm::dtype_t dtype = NM_DTYPE(vector);

  void* scalar = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, scalar);

  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::lapacke::cblas_scal, void, const int n,
      const void* scalar, void* x, const int incx);

  ttable[dtype](FIX2INT(n), scalar, NM_STORAGE_DENSE(vector)->elements,
      FIX2INT(incx));

  return vector;
}

/*
 * Call any of the cblas_xrotg functions as directly as possible.
 *
 * xROTG computes the elements of a Givens plane rotation matrix such that:
 *
 *  |  c s |   | a |   | r |
 *  | -s c | * | b | = | 0 |
 *
 * where r = +- sqrt( a**2 + b**2 ) and c**2 + s**2 = 1.
 *
 * The Givens plane rotation can be used to introduce zero elements into a matrix selectively.
 *
 * This function differs from most of the other raw BLAS accessors. Instead of
 * providing a, b, c, s as arguments, you should only provide a and b (the
 * inputs), and you should provide them as the first two elements of any dense
 * NMatrix type.
 *
 * The outputs [c,s] will be returned in a Ruby Array at the end; the input
 * NMatrix will also be modified in-place.
 *
 * This function, like the other cblas_ functions, does minimal type-checking.
 */
static VALUE nm_lapacke_cblas_rotg(VALUE self, VALUE ab) {
  static void (*ttable[nm::NUM_DTYPES])(void* a, void* b, void* c, void* s) = {
      NULL, NULL, NULL, NULL, NULL, // can't represent c and s as integers, so no point in having integer operations.
      nm::math::lapacke::cblas_rotg<float>,
      nm::math::lapacke::cblas_rotg<double>,
      nm::math::lapacke::cblas_rotg<nm::Complex64>,
      nm::math::lapacke::cblas_rotg<nm::Complex128>,
      NULL //nm::math::lapacke::cblas_rotg<nm::RubyObject>
  };

  nm::dtype_t dtype = NM_DTYPE(ab);

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this operation undefined for integer vectors");
    return Qnil;

  } else {
    NM_CONSERVATIVE(nm_register_value(&self));
    NM_CONSERVATIVE(nm_register_value(&ab));
    void *pC = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]),
         *pS = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);

    // extract A and B from the NVector (first two elements)
    void* pA = NM_STORAGE_DENSE(ab)->elements;
    void* pB = (char*)(NM_STORAGE_DENSE(ab)->elements) + DTYPE_SIZES[dtype];
    // c and s are output

    ttable[dtype](pA, pB, pC, pS);

    VALUE result = rb_ary_new2(2);

    if (dtype == nm::RUBYOBJ) {
      rb_ary_store(result, 0, *reinterpret_cast<VALUE*>(pC));
      rb_ary_store(result, 1, *reinterpret_cast<VALUE*>(pS));
    } else {
      rb_ary_store(result, 0, nm::rubyobj_from_cval(pC, dtype).rval);
      rb_ary_store(result, 1, nm::rubyobj_from_cval(pS, dtype).rval);
    }
    NM_CONSERVATIVE(nm_unregister_value(&ab));
    NM_CONSERVATIVE(nm_unregister_value(&self));
    return result;
  }
}


/*
 * Call any of the cblas_xrot functions as directly as possible.
 *
 * xROT is a BLAS level 1 routine (taking two vectors) which applies a plane rotation.
 *
 * It's tough to find documentation on xROT. Here are what we think the arguments are for:
 *  * n     :: number of elements to consider in x and y
 *  * x     :: a vector (expects an NVector)
 *  * incx  :: stride of x
 *  * y     :: a vector (expects an NVector)
 *  * incy  :: stride of y
 *  * c     :: cosine of the angle of rotation
 *  * s     :: sine of the angle of rotation
 *
 * Note that c and s will be the same dtype as x and y, except when x and y are complex. If x and y are complex, c and s
 * will be float for Complex64 or double for Complex128.
 *
 * You probably don't want to call this function. Instead, why don't you try rot, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_lapacke_cblas_rot(VALUE self, VALUE n, VALUE x, VALUE incx, VALUE y, VALUE incy, VALUE c, VALUE s) {
  static void (*ttable[nm::NUM_DTYPES])(const int N, void*, const int, void*, const int, const void*, const void*) = {
      NULL, NULL, NULL, NULL, NULL, // can't represent c and s as integers, so no point in having integer operations.
      nm::math::lapacke::cblas_rot<float,float>,
      nm::math::lapacke::cblas_rot<double,double>,
      nm::math::lapacke::cblas_rot<nm::Complex64,float>,
      nm::math::lapacke::cblas_rot<nm::Complex128,double>,
      nm::math::lapacke::cblas_rot<nm::RubyObject,nm::RubyObject>
  };

  nm::dtype_t dtype = NM_DTYPE(x);


  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this operation undefined for integer vectors");
    return Qfalse;
  } else {
    void *pC, *pS;

    // We need to ensure the cosine and sine arguments are the correct dtype -- which may differ from the actual dtype.
    if (dtype == nm::COMPLEX64) {
      pC = NM_ALLOCA_N(float,1);
      pS = NM_ALLOCA_N(float,1);
      rubyval_to_cval(c, nm::FLOAT32, pC);
      rubyval_to_cval(s, nm::FLOAT32, pS);
    } else if (dtype == nm::COMPLEX128) {
      pC = NM_ALLOCA_N(double,1);
      pS = NM_ALLOCA_N(double,1);
      rubyval_to_cval(c, nm::FLOAT64, pC);
      rubyval_to_cval(s, nm::FLOAT64, pS);
    } else {
      pC = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
      pS = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
      rubyval_to_cval(c, dtype, pC);
      rubyval_to_cval(s, dtype, pS);
    }


    ttable[dtype](FIX2INT(n), NM_STORAGE_DENSE(x)->elements, FIX2INT(incx), NM_STORAGE_DENSE(y)->elements, FIX2INT(incy), pC, pS);

    return Qtrue;
  }
}


/*
 * Call any of the cblas_xnrm2 functions as directly as possible.
 *
 * xNRM2 is a BLAS level 1 routine which calculates the 2-norm of an n-vector x.
 *
 * Arguments:
 *  * n     :: length of x, must be at least 0
 *  * x     :: pointer to first entry of input vector
 *  * incx  :: stride of x, must be POSITIVE (ATLAS says non-zero, but 3.8.4 code only allows positive)
 *
 * You probably don't want to call this function. Instead, why don't you try nrm2, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_lapacke_cblas_nrm2(VALUE self, VALUE n, VALUE x, VALUE incx) {

  static void (*ttable[nm::NUM_DTYPES])(const int N, const void* X, const int incX, void* sum) = {
      NULL, NULL, NULL, NULL, NULL, // no help for integers
      nm::math::lapacke::cblas_nrm2<float32_t>,
      nm::math::lapacke::cblas_nrm2<float64_t>,
      nm::math::lapacke::cblas_nrm2<nm::Complex64>,
      nm::math::lapacke::cblas_nrm2<nm::Complex128>,
      nm::math::lapacke::cblas_nrm2<nm::RubyObject>
  };

  nm::dtype_t dtype  = NM_DTYPE(x);

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this operation undefined for integer vectors");
    return Qnil;

  } else {
    // Determine the return dtype and allocate it
    nm::dtype_t rdtype = dtype;
    if      (dtype == nm::COMPLEX64)  rdtype = nm::FLOAT32;
    else if (dtype == nm::COMPLEX128) rdtype = nm::FLOAT64;

    void *Result = NM_ALLOCA_N(char, DTYPE_SIZES[rdtype]);

    ttable[dtype](FIX2INT(n), NM_STORAGE_DENSE(x)->elements, FIX2INT(incx), Result);

    return nm::rubyobj_from_cval(Result, rdtype).rval;
  }
}



/*
 * Call any of the cblas_xasum functions as directly as possible.
 *
 * xASUM is a BLAS level 1 routine which calculates the sum of absolute values of the entries
 * of a vector x.
 *
 * Arguments:
 *  * n     :: length of x, must be at least 0
 *  * x     :: pointer to first entry of input vector
 *  * incx  :: stride of x, must be POSITIVE (ATLAS says non-zero, but 3.8.4 code only allows positive)
 *
 * You probably don't want to call this function. Instead, why don't you try asum, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_lapacke_cblas_asum(VALUE self, VALUE n, VALUE x, VALUE incx) {

  static void (*ttable[nm::NUM_DTYPES])(const int N, const void* X, const int incX, void* sum) = {
      nm::math::lapacke::cblas_asum<uint8_t>,
      nm::math::lapacke::cblas_asum<int8_t>,
      nm::math::lapacke::cblas_asum<int16_t>,
      nm::math::lapacke::cblas_asum<int32_t>,
      nm::math::lapacke::cblas_asum<int64_t>,
      nm::math::lapacke::cblas_asum<float32_t>,
      nm::math::lapacke::cblas_asum<float64_t>,
      nm::math::lapacke::cblas_asum<nm::Complex64>,
      nm::math::lapacke::cblas_asum<nm::Complex128>,
      nm::math::lapacke::cblas_asum<nm::RubyObject>
  };

  nm::dtype_t dtype  = NM_DTYPE(x);

  // Determine the return dtype and allocate it
  nm::dtype_t rdtype = dtype;
  if      (dtype == nm::COMPLEX64)  rdtype = nm::FLOAT32;
  else if (dtype == nm::COMPLEX128) rdtype = nm::FLOAT64;

  void *Result = NM_ALLOCA_N(char, DTYPE_SIZES[rdtype]);

  ttable[dtype](FIX2INT(n), NM_STORAGE_DENSE(x)->elements, FIX2INT(incx), Result);

  return nm::rubyobj_from_cval(Result, rdtype).rval;
}

/*
 * call-seq:
 *    NMatrix::BLAS.cblas_imax(n, vector, inc) -> Fixnum
 *
 * BLAS level 1 routine.
 *
 * Return the index of the largest element of +vector+.
 *
 * - +n+ -> Vector's size. Generally, you can use NMatrix#rows or NMatrix#cols.
 * - +vector+ -> A NMatrix of shape [n,1] or [1,n] with any dtype.
 * - +inc+ -> It's the increment used when searching. Use 1 except if you know
 *   what you're doing.
 */
static VALUE nm_lapacke_cblas_imax(VALUE self, VALUE n, VALUE x, VALUE incx) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::lapacke::cblas_imax, int, const int n, const void* x, const int incx);

  nm::dtype_t dtype = NM_DTYPE(x);

  int index = ttable[dtype](FIX2INT(n), NM_STORAGE_DENSE(x)->elements, FIX2INT(incx));

  // Convert to Ruby's Int value.
  return INT2FIX(index);
}

/* Call any of the cblas_xgemv functions as directly as possible.
 *
 * The cblas_xgemv functions (dgemv, sgemv, cgemv, and zgemv) define the following operation:
 *
 *    y = alpha*op(A)*x + beta*y
 *
 * where op(A) is one of <tt>op(A) = A</tt>, <tt>op(A) = A**T</tt>, or the complex conjugate of A.
 *
 * Note that this will only work for dense matrices that are of types :float32, :float64, :complex64, and :complex128.
 * Other types are not implemented in BLAS, and while they exist in NMatrix, this method is intended only to
 * expose the ultra-optimized ATLAS versions.
 *
 * == Arguments
 * See: http://www.netlib.org/blas/dgemm.f
 *
 * You probably don't want to call this function. Instead, why don't you try cblas_gemv, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_lapacke_cblas_gemv(VALUE self,
                           VALUE trans_a,
                           VALUE m, VALUE n,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE x, VALUE incx,
                           VALUE beta,
                           VALUE y, VALUE incy)
{
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::lapacke::cblas_gemv, bool, const enum CBLAS_TRANSPOSE, const int, const int, const void*, const void*, const int, const void*, const int, const void*, void*, const int)

  nm::dtype_t dtype = NM_DTYPE(a);

  void *pAlpha = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]),
       *pBeta  = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, pAlpha);
  rubyval_to_cval(beta, dtype, pBeta);

  return ttable[dtype](blas_transpose_sym(trans_a), FIX2INT(m), FIX2INT(n), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(x)->elements, FIX2INT(incx), pBeta, NM_STORAGE_DENSE(y)->elements, FIX2INT(incy)) ? Qtrue : Qfalse;
}

/* Call any of the cblas_xgemm functions as directly as possible.
 *
 * The cblas_xgemm functions (dgemm, sgemm, cgemm, and zgemm) define the following operation:
 *
 *    C = alpha*op(A)*op(B) + beta*C
 *
 * where op(X) is one of <tt>op(X) = X</tt>, <tt>op(X) = X**T</tt>, or the complex conjugate of X.
 *
 * Note that this will only work for dense matrices that are of types :float32, :float64, :complex64, and :complex128.
 * Other types are not implemented in BLAS, and while they exist in NMatrix, this method is intended only to
 * expose the ultra-optimized ATLAS versions.
 *
 * == Arguments
 * See: http://www.netlib.org/blas/dgemm.f
 *
 * You probably don't want to call this function. Instead, why don't you try gemm, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_lapacke_cblas_gemm(VALUE self,
                           VALUE order,
                           VALUE trans_a, VALUE trans_b,
                           VALUE m, VALUE n, VALUE k,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb,
                           VALUE beta,
                           VALUE c, VALUE ldc)
{
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::lapacke::cblas_gemm, void, const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_TRANSPOSE trans_b, int m, int n, int k, void* alpha, void* a, int lda, void* b, int ldb, void* beta, void* c, int ldc);

  nm::dtype_t dtype = NM_DTYPE(a);

  void *pAlpha = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]),
       *pBeta  = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, pAlpha);
  rubyval_to_cval(beta, dtype, pBeta);

  ttable[dtype](blas_order_sym(order), blas_transpose_sym(trans_a), blas_transpose_sym(trans_b), FIX2INT(m), FIX2INT(n), FIX2INT(k), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb), pBeta, NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc));

  return c;
}


static VALUE nm_lapacke_cblas_trsm(VALUE self,
                           VALUE order,
                           VALUE side, VALUE uplo,
                           VALUE trans_a, VALUE diag,
                           VALUE m, VALUE n,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb)
{
  static void (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const enum CBLAS_SIDE, const enum CBLAS_UPLO,
                                        const enum CBLAS_TRANSPOSE, const enum CBLAS_DIAG,
                                        const int m, const int n, const void* alpha, const void* a,
                                        const int lda, void* b, const int ldb) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::lapacke::cblas_trsm<float>,
      nm::math::lapacke::cblas_trsm<double>,
      cblas_ctrsm, cblas_ztrsm, // call directly, same function signature!
      nm::math::lapacke::cblas_trsm<nm::RubyObject>
  };

  nm::dtype_t dtype = NM_DTYPE(a);

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    void *pAlpha = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
    rubyval_to_cval(alpha, dtype, pAlpha);

    ttable[dtype](blas_order_sym(order), blas_side_sym(side), blas_uplo_sym(uplo), blas_transpose_sym(trans_a), blas_diag_sym(diag), FIX2INT(m), FIX2INT(n), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb));
  }

  return Qtrue;
}

static VALUE nm_lapacke_cblas_trmm(VALUE self,
                           VALUE order,
                           VALUE side, VALUE uplo,
                           VALUE trans_a, VALUE diag,
                           VALUE m, VALUE n,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb)
{
  static void (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER,
                                        const enum CBLAS_SIDE, const enum CBLAS_UPLO,
                                        const enum CBLAS_TRANSPOSE, const enum CBLAS_DIAG,
                                        const int m, const int n, const void* alpha, const void* a,
                                        const int lda, void* b, const int ldb) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::lapacke::cblas_trmm<float>,
      nm::math::lapacke::cblas_trmm<double>,
      cblas_ctrmm, cblas_ztrmm, // call directly, same function signature!
      NULL
  };

  nm::dtype_t dtype = NM_DTYPE(a);

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this matrix operation not yet defined for non-BLAS dtypes");
  } else {
    void *pAlpha = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
    rubyval_to_cval(alpha, dtype, pAlpha);

    ttable[dtype](blas_order_sym(order), blas_side_sym(side), blas_uplo_sym(uplo), blas_transpose_sym(trans_a), blas_diag_sym(diag), FIX2INT(m), FIX2INT(n), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb));
  }

  return b;
}

static VALUE nm_lapacke_cblas_syrk(VALUE self,
                           VALUE order,
                           VALUE uplo,
                           VALUE trans,
                           VALUE n, VALUE k,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE beta,
                           VALUE c, VALUE ldc)
{
  static void (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const enum CBLAS_UPLO, const enum CBLAS_TRANSPOSE,
                                        const int n, const int k, const void* alpha, const void* a,
                                        const int lda, const void* beta, void* c, const int ldc) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::lapacke::cblas_syrk<float>,
      nm::math::lapacke::cblas_syrk<double>,
      cblas_csyrk, cblas_zsyrk, // call directly, same function signature!
      NULL
  };

  nm::dtype_t dtype = NM_DTYPE(a);

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    void *pAlpha = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]),
         *pBeta = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
    rubyval_to_cval(alpha, dtype, pAlpha);
    rubyval_to_cval(beta, dtype, pBeta);

    ttable[dtype](blas_order_sym(order), blas_uplo_sym(uplo), blas_transpose_sym(trans), FIX2INT(n), FIX2INT(k), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), pBeta, NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc));
  }

  return Qtrue;
}

static VALUE nm_lapacke_cblas_herk(VALUE self,
                           VALUE order,
                           VALUE uplo,
                           VALUE trans,
                           VALUE n, VALUE k,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE beta,
                           VALUE c, VALUE ldc)
{

  nm::dtype_t dtype = NM_DTYPE(a);

  if (dtype == nm::COMPLEX64) {
    cblas_cherk(blas_order_sym(order), blas_uplo_sym(uplo), blas_transpose_sym(trans), FIX2INT(n), FIX2INT(k), NUM2DBL(alpha), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NUM2DBL(beta), NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc));
  } else if (dtype == nm::COMPLEX128) {
    cblas_zherk(blas_order_sym(order), blas_uplo_sym(uplo), blas_transpose_sym(trans), FIX2INT(n), FIX2INT(k), NUM2DBL(alpha), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NUM2DBL(beta), NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc));
  } else
    rb_raise(rb_eNotImpError, "this matrix operation undefined for non-complex dtypes");
  return Qtrue;
}

/* Call any of the lapacke_xgetri functions as directly as possible.
 *
 * This version (the LAPACKE version) differs from the CLAPACK version in terms of the
 * input it expects (which is the output of getrf). See getrf for details.
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 *
 * Returns an array giving the pivot indices (normally these are argument #5).
 */
static VALUE nm_lapacke_lapacke_getri(VALUE self, VALUE order, VALUE n, VALUE a, VALUE lda, VALUE ipiv) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const int n, void* a, const int lda, const int* ipiv) = {
      NULL, NULL, NULL, NULL, NULL,
      nm::math::lapacke::lapacke_getri<float>,
      nm::math::lapacke::lapacke_getri<double>,
      nm::math::lapacke::lapacke_getri<nm::Complex64>,
      nm::math::lapacke::lapacke_getri<nm::Complex128>,
      NULL
  };

  // Allocate the C version of the pivot index array
  int* ipiv_;
  if (!RB_TYPE_P(ipiv, T_ARRAY)) {
    rb_raise(rb_eArgError, "ipiv must be of type Array");
  } else {
    ipiv_ = NM_ALLOCA_N(int, RARRAY_LEN(ipiv));
    for (int index = 0; index < RARRAY_LEN(ipiv); ++index) {
      ipiv_[index] = FIX2INT( RARRAY_AREF(ipiv, index) );
    }
  }

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
  } else {
    ttable[NM_DTYPE(a)](blas_order_sym(order), FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), ipiv_);
  }

  return a;
}

/* Call any of the lapacke_xgetrf functions as directly as possible.
 *
 * The lapacke_getrf functions (dgetrf, sgetrf, cgetrf, and zgetrf) compute an LU factorization of a general M-by-N
 * matrix A using partial pivoting with row interchanges.
 *
 * The factorization has the form:
 *    A = P * L * U
 * where P is a permutation matrix, L is lower triangular with unit diagonal elements (lower trapezoidal if m > n),
 * and U is upper triangular (upper trapezoidal if m < n).
 *
 * This version of getrf (the LAPACKE one) differs from the CLAPACK version. The CLAPACK has
 * different behavior for row-major matrices (the upper matrix has unit diagonals instead of
 * the lower and it uses column permutations instead of rows).
 *
 * This is the right-looking level 3 BLAS version of the algorithm.
 *
 * == Arguments
 * See: http://www.netlib.org/lapack/double/dgetrf.f
 * (You don't need argument 5; this is the value returned by this function.)
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 *
 * Returns an array giving the pivot indices (normally these are argument #5).
 */
static VALUE nm_lapacke_lapacke_getrf(VALUE self, VALUE order, VALUE m, VALUE n, VALUE a, VALUE lda) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const int m, const int n, void* a, const int lda, int* ipiv) = {
      NULL, NULL, NULL, NULL, NULL,
      nm::math::lapacke::lapacke_getrf<float>,
      nm::math::lapacke::lapacke_getrf<double>,
      nm::math::lapacke::lapacke_getrf<nm::Complex64>,
      nm::math::lapacke::lapacke_getrf<nm::Complex128>,
      NULL
  };

  int M = FIX2INT(m),
      N = FIX2INT(n);

  // Allocate the pivot index array, which is of size MIN(M, N).
  size_t ipiv_size = std::min(M,N);
  int* ipiv = NM_ALLOCA_N(int, ipiv_size);

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    ttable[NM_DTYPE(a)](blas_order_sym(order), M, N, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), ipiv);
  }

  // Result will be stored in a. We return ipiv as an array.
  VALUE ipiv_array = rb_ary_new2(ipiv_size);
  for (size_t i = 0; i < ipiv_size; ++i) {
    rb_ary_store(ipiv_array, i, INT2FIX(ipiv[i]));
  }

  return ipiv_array;
}

/*
 * Call any of the lapacke_xgetrs functions as directly as possible.
 */
static VALUE nm_lapacke_lapacke_getrs(VALUE self, VALUE order, VALUE trans, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE ipiv, VALUE b, VALUE ldb) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER Order, char Trans, const int N,
                                       const int NRHS, const void* A, const int lda, const int* ipiv, void* B,
                                       const int ldb) = {
      NULL, NULL, NULL, NULL, NULL,
      nm::math::lapacke::lapacke_getrs<float>,
      nm::math::lapacke::lapacke_getrs<double>,
      nm::math::lapacke::lapacke_getrs<nm::Complex64>,
      nm::math::lapacke::lapacke_getrs<nm::Complex128>,
      NULL
  };

  // Allocate the C version of the pivot index array
  int* ipiv_;
  if (!RB_TYPE_P(ipiv, T_ARRAY)) {
    rb_raise(rb_eArgError, "ipiv must be of type Array");
  } else {
    ipiv_ = NM_ALLOCA_N(int, RARRAY_LEN(ipiv));
    for (int index = 0; index < RARRAY_LEN(ipiv); ++index) {
      ipiv_[index] = FIX2INT( RARRAY_AREF(ipiv, index) );
    }
  }

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    ttable[NM_DTYPE(a)](blas_order_sym(order), lapacke_transpose_sym(trans), FIX2INT(n), FIX2INT(nrhs), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda),
                        ipiv_, NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb));
  }

  // b is both returned and modified directly in the argument list.
  return b;
}

/* Call any of the LAPACKE_xpotrf functions as directly as possible.
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_lapacke_lapacke_potrf(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda) {

  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, char, const int n, void* a, const int lda) = {
      NULL, NULL, NULL, NULL, NULL,
      nm::math::lapacke::lapacke_potrf<float>,
      nm::math::lapacke::lapacke_potrf<double>,
      nm::math::lapacke::lapacke_potrf<nm::Complex64>,
      nm::math::lapacke::lapacke_potrf<nm::Complex128>,
      NULL
  };

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
  } else {
    ttable[NM_DTYPE(a)](blas_order_sym(order), lapacke_uplo_sym(uplo), FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda));
  }

  return a;
}

/*
 * Call any of the LAPACKE_xpotrs functions as directly as possible.
 */
static VALUE nm_lapacke_lapacke_potrs(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE b, VALUE ldb) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER Order, char Uplo, const int N,
                                       const int NRHS, const void* A, const int lda, void* B, const int ldb) = {
      NULL, NULL, NULL, NULL, NULL,
      nm::math::lapacke::lapacke_potrs<float>,
      nm::math::lapacke::lapacke_potrs<double>,
      nm::math::lapacke::lapacke_potrs<nm::Complex64>,
      nm::math::lapacke::lapacke_potrs<nm::Complex128>,
      NULL
  };


  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {

    ttable[NM_DTYPE(a)](blas_order_sym(order), lapacke_uplo_sym(uplo), FIX2INT(n), FIX2INT(nrhs), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda),
                        NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb));
  }

  // b is both returned and modified directly in the argument list.
  return b;
}

/* Call any of the lapacke_xpotri functions as directly as possible.
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_lapacke_lapacke_potri(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda) {

  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, char, const int n, void* a, const int lda) = {
      NULL, NULL, NULL, NULL, NULL,
      nm::math::lapacke::lapacke_potri<float>,
      nm::math::lapacke::lapacke_potri<double>,
      nm::math::lapacke::lapacke_potri<nm::Complex64>,
      nm::math::lapacke::lapacke_potri<nm::Complex128>,
      NULL
  };

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
  } else {
    ttable[NM_DTYPE(a)](blas_order_sym(order), lapacke_uplo_sym(uplo), FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda));
  }

  return a;
}

//badly need docs for gesvd, gesdd because of the real/complex mixing

/*
 * xGESVD computes the singular value decomposition (SVD) of a real
 * M-by-N matrix A, optionally computing the left and/or right singular
 * vectors. The SVD is written
 *
 *      A = U * SIGMA * transpose(V)
 *
 * where SIGMA is an M-by-N matrix which is zero except for its
 * min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
 * V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
 * are the singular values of A; they are real and non-negative, and
 * are returned in descending order.  The first min(m,n) columns of
 * U and V are the left and right singular vectors of A.
 *
 * Note that the routine returns V**T, not V.
 */
static VALUE nm_lapacke_lapacke_gesvd(VALUE self, VALUE order, VALUE jobu, VALUE jobvt, VALUE m, VALUE n, VALUE a, VALUE lda, VALUE s, VALUE u, VALUE ldu, VALUE vt, VALUE ldvt, VALUE superb) {
  static int (*gesvd_table[nm::NUM_DTYPES])(int, char, char, int, int, void* a, int, void* s, void* u, int, void* vt, int, void* superb) = {
    NULL, NULL, NULL, NULL, NULL, // no integer ops
    nm::math::lapacke::lapacke_gesvd<float,float>,
    nm::math::lapacke::lapacke_gesvd<double,double>,
    nm::math::lapacke::lapacke_gesvd<nm::Complex64,float>,
    nm::math::lapacke::lapacke_gesvd<nm::Complex128,double>,
    NULL // no Ruby objects
  };

  nm::dtype_t dtype = NM_DTYPE(a);


  if (!gesvd_table[dtype]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
    return Qfalse;
  } else {
    int M = FIX2INT(m),
        N = FIX2INT(n);

    char JOBU = lapack_svd_job_sym(jobu),
         JOBVT = lapack_svd_job_sym(jobvt);

    int info = gesvd_table[dtype](blas_order_sym(order),JOBU, JOBVT, M, N, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda),
      NM_STORAGE_DENSE(s)->elements, NM_STORAGE_DENSE(u)->elements, FIX2INT(ldu), NM_STORAGE_DENSE(vt)->elements, FIX2INT(ldvt),
      NM_STORAGE_DENSE(superb)->elements);
    return INT2FIX(info);
  }
}

static VALUE nm_lapacke_lapacke_gesdd(VALUE self, VALUE order, VALUE jobz, VALUE m, VALUE n, VALUE a, VALUE lda, VALUE s, VALUE u, VALUE ldu, VALUE vt, VALUE ldvt) {
  static int (*gesdd_table[nm::NUM_DTYPES])(int, char, int, int, void* a, int, void* s, void* u, int, void* vt, int) = {
    NULL, NULL, NULL, NULL, NULL, // no integer ops
    nm::math::lapacke::lapacke_gesdd<float,float>,
    nm::math::lapacke::lapacke_gesdd<double,double>,
    nm::math::lapacke::lapacke_gesdd<nm::Complex64,float>,
    nm::math::lapacke::lapacke_gesdd<nm::Complex128,double>,
    NULL // no Ruby objects
  };

  nm::dtype_t dtype = NM_DTYPE(a);


  if (!gesdd_table[dtype]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
    return Qfalse;
  } else {
    int M = FIX2INT(m),
        N = FIX2INT(n);

    char JOBZ = lapack_svd_job_sym(jobz);

    int info = gesdd_table[dtype](blas_order_sym(order),JOBZ, M, N, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda),
      NM_STORAGE_DENSE(s)->elements, NM_STORAGE_DENSE(u)->elements, FIX2INT(ldu), NM_STORAGE_DENSE(vt)->elements, FIX2INT(ldvt));
    return INT2FIX(info);
  }
}

/*
 * GEEV computes for an N-by-N real nonsymmetric matrix A, the
 * eigenvalues and, optionally, the left and/or right eigenvectors.
 *
 * The right eigenvector v(j) of A satisfies
 *                    A * v(j) = lambda(j) * v(j)
 * where lambda(j) is its eigenvalue.
 *
 * The left eigenvector u(j) of A satisfies
 *                 u(j)**H * A = lambda(j) * u(j)**H
 * where u(j)**H denotes the conjugate transpose of u(j).
 *
 * The computed eigenvectors are normalized to have Euclidean norm
 * equal to 1 and largest component real.
 */
//note on wi
static VALUE nm_lapacke_lapacke_geev(VALUE self, VALUE order, VALUE jobvl, VALUE jobvr, VALUE n, VALUE a, VALUE lda, VALUE w, VALUE wi, VALUE vl, VALUE ldvl, VALUE vr, VALUE ldvr) {
  static int (*geev_table[nm::NUM_DTYPES])(int, char, char, int, void* a, int, void* w, void* wi, void* vl, int, void* vr, int) = {
    NULL, NULL, NULL, NULL, NULL, // no integer ops
    nm::math::lapacke::lapacke_geev<float>,
    nm::math::lapacke::lapacke_geev<double>,
    nm::math::lapacke::lapacke_geev<nm::Complex64>,
    nm::math::lapacke::lapacke_geev<nm::Complex128>,
    NULL // no Ruby objects
  };

  nm::dtype_t dtype = NM_DTYPE(a);


  if (!geev_table[dtype]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
    return Qfalse;
  } else {
    int N = FIX2INT(n);

    char JOBVL = lapack_evd_job_sym(jobvl),
         JOBVR = lapack_evd_job_sym(jobvr);

    void* A  = NM_STORAGE_DENSE(a)->elements;
    void* W = NM_STORAGE_DENSE(w)->elements;
    void* WI = wi == Qnil ? NULL : NM_STORAGE_DENSE(wi)->elements; //For complex, wi should be nil
    void* VL = JOBVL == 'V' ? NM_STORAGE_DENSE(vl)->elements : NULL;
    void* VR = JOBVR == 'V' ? NM_STORAGE_DENSE(vr)->elements : NULL;

    // Perform the actual calculation.
    int info = geev_table[dtype](blas_order_sym(order), JOBVL, JOBVR, N, A, FIX2INT(lda), W, WI, VL, FIX2INT(ldvl), VR, FIX2INT(ldvr));

    return INT2FIX(info);
  }
}

/* 
 * GEQRF calculates the QR factorization for an MxN real or complex matrix.
 *  
 * The QR factorization is A = QR, where Q is orthogonal and R is Upper Triangular
 * +A+ is overwritten with the elements of R and Q with Q being represented by the 
 * elements below A's diagonal and an array of scalar factors in the output NMatrix. 
 *
 * The matrix Q is represented as a product of elementary reflectors
 *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
 *
 * Each H(i) has the form
 *
 *     H(i) = I - tau * v * v'
 *
 * http://www.netlib.org/lapack/explore-html/d3/d69/dgeqrf_8f.html
 */

static VALUE nm_lapacke_lapacke_geqrf(VALUE self, VALUE order, VALUE m, VALUE n, VALUE a, VALUE lda, VALUE tau) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER order, const int m, const int n, void* a, const int lda, void* tau) = {
      NULL, NULL, NULL, NULL, NULL,
      nm::math::lapacke::lapacke_geqrf<float>,
      nm::math::lapacke::lapacke_geqrf<double>,
      nm::math::lapacke::lapacke_geqrf<nm::Complex64>,
      nm::math::lapacke::lapacke_geqrf<nm::Complex128>,
      NULL
  };

  int M = FIX2INT(m),
      N = FIX2INT(n);
 
  nm::dtype_t dtype = NM_DTYPE(a);

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this matrix operation is undefined for integer matrices");
    return Qfalse;
  } else {
    int info = ttable[dtype](blas_order_sym(order), M, N, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(tau)->elements);
    return INT2FIX(info);
  }
}

/* ORMQR calculates the orthogonal matrix Q from TAU and A after calling GEQRF on a real matrix
 *  
 *
 * The matrix Q is represented as a product of elementary reflectors
 *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
 *
 * Each H(i) has the form
 *
 *     H(i) = I - tau * v * v'
 *  
 *  v is contained in the matrix passed to GEQRF     
 *
 *  www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
 */

static VALUE nm_lapacke_lapacke_ormqr(VALUE self, VALUE order, VALUE side, VALUE trans, VALUE m, VALUE n, VALUE k, VALUE a, VALUE lda, VALUE tau, VALUE c, VALUE ldc) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER order, char side, char trans, const int m, const int n, const int k, void* a, const int lda, void* tau, void* c, const int ldc) = {
      NULL, NULL, NULL, NULL, NULL,
      nm::math::lapacke::lapacke_ormqr<float>,
      nm::math::lapacke::lapacke_ormqr<double>,
      NULL,NULL,NULL // no complex or Ruby objects
  };

  int M = FIX2INT(m),
      N = FIX2INT(n),
      K = FIX2INT(k); 

  char SIDE  = lapacke_side_sym(side),
       TRANS = lapacke_transpose_sym(trans);

  nm::dtype_t dtype = NM_DTYPE(a);


  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this matrix operation is undefined for integer matrices");
    return Qfalse;
  } else {
    int info = ttable[dtype](blas_order_sym(order), SIDE, TRANS, M, N, K, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(tau)->elements, NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc));
    return INT2FIX(info);
  }
}

/* UNMQR calculates the orthogonal matrix Q from TAU and A after calling GEQRF on a complex matrix.
 *  
 *
 * The matrix Q is represented as a product of elementary reflectors
 *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
 *
 * Each H(i) has the form
 *
 *     H(i) = I - tau * v * v'
 *  
 *  v is contained in the matrix passed to GEQRF     
 *
 *  http://www.netlib.org/lapack/explore-html/d5/d65/zunmqr_8f.html
 */

static VALUE nm_lapacke_lapacke_unmqr(VALUE self, VALUE order, VALUE side, VALUE trans, VALUE m, VALUE n, VALUE k, VALUE a, VALUE lda, VALUE tau, VALUE c, VALUE ldc) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER order, char side, char trans, const int m, const int n, const int k, void* a, const int lda, void* tau, void* c, const int ldc) = {
      NULL, NULL, NULL, NULL, NULL,NULL,NULL, // no non-complex ops
      nm::math::lapacke::lapacke_unmqr<nm::Complex64>,
      nm::math::lapacke::lapacke_unmqr<nm::Complex128>,
      NULL // no Ruby objects
  };

  int M = FIX2INT(m),
      N = FIX2INT(n),
      K = FIX2INT(k); 

  char SIDE  = lapacke_side_sym(side),
       TRANS = lapacke_transpose_sym(trans);

  nm::dtype_t dtype = NM_DTYPE(a);

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this matrix operation is valid only for complex datatypes");
    return Qfalse;
  } else {
    int info = ttable[dtype](blas_order_sym(order), SIDE, TRANS, M, N, K, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(tau)->elements, NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc));
    return INT2FIX(info);
  }
}

}
