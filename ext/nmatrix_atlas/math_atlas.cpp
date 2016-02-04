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
// == math_atlas.cpp
//
// Ruby-exposed CBLAS and LAPACK functions that call ATLAS
// functions.
//

/*
 * Project Includes
 */

#include "data/data.h"

#include "math_atlas/inc.h"

#include "math/util.h"

//BLAS
#include "math_atlas/cblas_templates_atlas.h"

//LAPACK
#include "math/laswp.h"
#include "math_atlas/clapack_templates.h"

#include "math_atlas/gesvd.h"
#include "math_atlas/gesdd.h"
#include "math_atlas/geev.h"


/*
 * Forward Declarations
 */

extern "C" {
  /* BLAS Level 1. */
  static VALUE nm_atlas_cblas_scal(VALUE self, VALUE n, VALUE scale, VALUE vector, VALUE incx);
  static VALUE nm_atlas_cblas_nrm2(VALUE self, VALUE n, VALUE x, VALUE incx);
  static VALUE nm_atlas_cblas_asum(VALUE self, VALUE n, VALUE x, VALUE incx);
  static VALUE nm_atlas_cblas_rot(VALUE self, VALUE n, VALUE x, VALUE incx, VALUE y, VALUE incy, VALUE c, VALUE s);
  static VALUE nm_atlas_cblas_rotg(VALUE self, VALUE ab);
  static VALUE nm_atlas_cblas_imax(VALUE self, VALUE n, VALUE x, VALUE incx);

  /* BLAS Level 2. */
  static VALUE nm_atlas_cblas_gemv(VALUE self, VALUE trans_a, VALUE m, VALUE n, VALUE vAlpha, VALUE a, VALUE lda,
                             VALUE x, VALUE incx, VALUE vBeta, VALUE y, VALUE incy);

  /* BLAS Level 3. */
  static VALUE nm_atlas_cblas_gemm(VALUE self, VALUE order, VALUE trans_a, VALUE trans_b, VALUE m, VALUE n, VALUE k, VALUE vAlpha,
                             VALUE a, VALUE lda, VALUE b, VALUE ldb, VALUE vBeta, VALUE c, VALUE ldc);
  static VALUE nm_atlas_cblas_trsm(VALUE self, VALUE order, VALUE side, VALUE uplo, VALUE trans_a, VALUE diag, VALUE m, VALUE n,
                             VALUE vAlpha, VALUE a, VALUE lda, VALUE b, VALUE ldb);
  static VALUE nm_atlas_cblas_trmm(VALUE self, VALUE order, VALUE side, VALUE uplo, VALUE trans_a, VALUE diag, VALUE m, VALUE n,
                             VALUE alpha, VALUE a, VALUE lda, VALUE b, VALUE ldb);
  static VALUE nm_atlas_cblas_herk(VALUE self, VALUE order, VALUE uplo, VALUE trans, VALUE n, VALUE k, VALUE alpha, VALUE a,
                             VALUE lda, VALUE beta, VALUE c, VALUE ldc);
  static VALUE nm_atlas_cblas_syrk(VALUE self, VALUE order, VALUE uplo, VALUE trans, VALUE n, VALUE k, VALUE alpha, VALUE a,
                             VALUE lda, VALUE beta, VALUE c, VALUE ldc);

  /* LAPACK. */
  static VALUE nm_atlas_has_clapack(VALUE self);
  static VALUE nm_atlas_clapack_getrf(VALUE self, VALUE order, VALUE m, VALUE n, VALUE a, VALUE lda);
  static VALUE nm_atlas_clapack_potrf(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda);
  static VALUE nm_atlas_clapack_getrs(VALUE self, VALUE order, VALUE trans, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE ipiv, VALUE b, VALUE ldb);
  static VALUE nm_atlas_clapack_potrs(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE b, VALUE ldb);
  static VALUE nm_atlas_clapack_getri(VALUE self, VALUE order, VALUE n, VALUE a, VALUE lda, VALUE ipiv);
  static VALUE nm_atlas_clapack_potri(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda);
  static VALUE nm_atlas_clapack_laswp(VALUE self, VALUE n, VALUE a, VALUE lda, VALUE k1, VALUE k2, VALUE ipiv, VALUE incx);

  static VALUE nm_atlas_lapack_gesvd(VALUE self, VALUE jobu, VALUE jobvt, VALUE m, VALUE n, VALUE a, VALUE lda, VALUE s, VALUE u, VALUE ldu, VALUE vt, VALUE ldvt, VALUE lworkspace_size);
  static VALUE nm_atlas_lapack_gesdd(VALUE self, VALUE jobz, VALUE m, VALUE n, VALUE a, VALUE lda, VALUE s, VALUE u, VALUE ldu, VALUE vt, VALUE ldvt, VALUE lworkspace_size);
  static VALUE nm_atlas_lapack_geev(VALUE self, VALUE compute_left, VALUE compute_right, VALUE n, VALUE a, VALUE lda, VALUE w, VALUE wi, VALUE vl, VALUE ldvl, VALUE vr, VALUE ldvr, VALUE lwork);
}

////////////////////
// Math Functions //
////////////////////

namespace nm { 
  namespace math {
  namespace atlas {

    /*
     * Function signature conversion for calling CBLAS' gesvd functions as directly as possible.
     */
    template <typename DType, typename CType>
    inline static int lapack_gesvd(char jobu, char jobvt, int m, int n, void* a, int lda, void* s, void* u, int ldu, void* vt, int ldvt, void* work, int lwork, void* rwork) {
      return gesvd<DType,CType>(jobu, jobvt, m, n, reinterpret_cast<DType*>(a), lda, reinterpret_cast<CType*>(s), reinterpret_cast<DType*>(u), ldu, reinterpret_cast<DType*>(vt), ldvt, reinterpret_cast<DType*>(work), lwork, reinterpret_cast<CType*>(rwork));
    }

    /*
     * Function signature conversion for calling CBLAS' gesdd functions as directly as possible.
     */
    template <typename DType, typename CType>
    inline static int lapack_gesdd(char jobz, int m, int n, void* a, int lda, void* s, void* u, int ldu, void* vt, int ldvt, void* work, int lwork, int* iwork, void* rwork) {
      return gesdd<DType,CType>(jobz, m, n, reinterpret_cast<DType*>(a), lda, reinterpret_cast<CType*>(s), reinterpret_cast<DType*>(u), ldu, reinterpret_cast<DType*>(vt), ldvt, reinterpret_cast<DType*>(work), lwork, iwork, reinterpret_cast<CType*>(rwork));
    }


  }
  }
}

extern "C" {

///////////////////
// Ruby Bindings //
///////////////////

void nm_math_init_atlas() {
  VALUE cNMatrix_ATLAS = rb_define_module_under(cNMatrix, "ATLAS");

  rb_define_singleton_method(cNMatrix, "has_clapack?", (METHOD)nm_atlas_has_clapack, 0);

  VALUE cNMatrix_ATLAS_LAPACK = rb_define_module_under(cNMatrix_ATLAS, "LAPACK");

  /* ATLAS-CLAPACK Functions */
  rb_define_singleton_method(cNMatrix_ATLAS_LAPACK, "clapack_getrf", (METHOD)nm_atlas_clapack_getrf, 5);
  rb_define_singleton_method(cNMatrix_ATLAS_LAPACK, "clapack_potrf", (METHOD)nm_atlas_clapack_potrf, 5);
  rb_define_singleton_method(cNMatrix_ATLAS_LAPACK, "clapack_getrs", (METHOD)nm_atlas_clapack_getrs, 9);
  rb_define_singleton_method(cNMatrix_ATLAS_LAPACK, "clapack_potrs", (METHOD)nm_atlas_clapack_potrs, 8);
  rb_define_singleton_method(cNMatrix_ATLAS_LAPACK, "clapack_getri", (METHOD)nm_atlas_clapack_getri, 5);
  rb_define_singleton_method(cNMatrix_ATLAS_LAPACK, "clapack_potri", (METHOD)nm_atlas_clapack_potri, 5);
  rb_define_singleton_method(cNMatrix_ATLAS_LAPACK, "clapack_laswp", (METHOD)nm_atlas_clapack_laswp, 7);

  /* Non-ATLAS regular LAPACK Functions called via Fortran interface */
  rb_define_singleton_method(cNMatrix_ATLAS_LAPACK, "lapack_gesvd", (METHOD)nm_atlas_lapack_gesvd, 12);
  rb_define_singleton_method(cNMatrix_ATLAS_LAPACK, "lapack_gesdd", (METHOD)nm_atlas_lapack_gesdd, 11);
  rb_define_singleton_method(cNMatrix_ATLAS_LAPACK, "lapack_geev",  (METHOD)nm_atlas_lapack_geev,  12);

  VALUE cNMatrix_ATLAS_BLAS = rb_define_module_under(cNMatrix_ATLAS, "BLAS");

  //BLAS Level 1
  rb_define_singleton_method(cNMatrix_ATLAS_BLAS, "cblas_scal", (METHOD)nm_atlas_cblas_scal, 4);
  rb_define_singleton_method(cNMatrix_ATLAS_BLAS, "cblas_nrm2", (METHOD)nm_atlas_cblas_nrm2, 3);
  rb_define_singleton_method(cNMatrix_ATLAS_BLAS, "cblas_asum", (METHOD)nm_atlas_cblas_asum, 3);
  rb_define_singleton_method(cNMatrix_ATLAS_BLAS, "cblas_rot",  (METHOD)nm_atlas_cblas_rot,  7);
  rb_define_singleton_method(cNMatrix_ATLAS_BLAS, "cblas_rotg", (METHOD)nm_atlas_cblas_rotg, 1);
  rb_define_singleton_method(cNMatrix_ATLAS_BLAS, "cblas_imax", (METHOD)nm_atlas_cblas_imax, 3);

  //BLAS Level 2
  rb_define_singleton_method(cNMatrix_ATLAS_BLAS, "cblas_gemv", (METHOD)nm_atlas_cblas_gemv, 11);

  //BLAS Level 3
  rb_define_singleton_method(cNMatrix_ATLAS_BLAS, "cblas_gemm", (METHOD)nm_atlas_cblas_gemm, 14);
  rb_define_singleton_method(cNMatrix_ATLAS_BLAS, "cblas_trsm", (METHOD)nm_atlas_cblas_trsm, 12);
  rb_define_singleton_method(cNMatrix_ATLAS_BLAS, "cblas_trmm", (METHOD)nm_atlas_cblas_trmm, 12);
  rb_define_singleton_method(cNMatrix_ATLAS_BLAS, "cblas_syrk", (METHOD)nm_atlas_cblas_syrk, 11);
  rb_define_singleton_method(cNMatrix_ATLAS_BLAS, "cblas_herk", (METHOD)nm_atlas_cblas_herk, 11);

}

/*
 * Simple way to check from within Ruby code if clapack functions are available, without
 * having to wait around for an exception to be thrown.
 */
static VALUE nm_atlas_has_clapack(VALUE self) {
#if defined (HAVE_CLAPACK_H) || defined (HAVE_ATLAS_CLAPACK_H)
  return Qtrue;
#else
  return Qfalse;
#endif
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
static VALUE nm_atlas_cblas_scal(VALUE self, VALUE n, VALUE alpha, VALUE vector, VALUE incx) {
  nm::dtype_t dtype = NM_DTYPE(vector);

  void* scalar = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, scalar);

  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::atlas::cblas_scal, void, const int n,
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
static VALUE nm_atlas_cblas_rotg(VALUE self, VALUE ab) {
  static void (*ttable[nm::NUM_DTYPES])(void* a, void* b, void* c, void* s) = {
      NULL, NULL, NULL, NULL, NULL, // can't represent c and s as integers, so no point in having integer operations.
      nm::math::atlas::cblas_rotg<float>,
      nm::math::atlas::cblas_rotg<double>,
      nm::math::atlas::cblas_rotg<nm::Complex64>,
      nm::math::atlas::cblas_rotg<nm::Complex128>,
      NULL //nm::math::atlas::cblas_rotg<nm::RubyObject>
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
static VALUE nm_atlas_cblas_rot(VALUE self, VALUE n, VALUE x, VALUE incx, VALUE y, VALUE incy, VALUE c, VALUE s) {
  static void (*ttable[nm::NUM_DTYPES])(const int N, void*, const int, void*, const int, const void*, const void*) = {
      NULL, NULL, NULL, NULL, NULL, // can't represent c and s as integers, so no point in having integer operations.
      nm::math::atlas::cblas_rot<float,float>,
      nm::math::atlas::cblas_rot<double,double>,
      nm::math::atlas::cblas_rot<nm::Complex64,float>,
      nm::math::atlas::cblas_rot<nm::Complex128,double>,
      nm::math::atlas::cblas_rot<nm::RubyObject,nm::RubyObject>
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
static VALUE nm_atlas_cblas_nrm2(VALUE self, VALUE n, VALUE x, VALUE incx) {

  static void (*ttable[nm::NUM_DTYPES])(const int N, const void* X, const int incX, void* sum) = {
      NULL, NULL, NULL, NULL, NULL, // no help for integers
      nm::math::atlas::cblas_nrm2<float32_t,float32_t>,
      nm::math::atlas::cblas_nrm2<float64_t,float64_t>,
      nm::math::atlas::cblas_nrm2<float32_t,nm::Complex64>,
      nm::math::atlas::cblas_nrm2<float64_t,nm::Complex128>,
      nm::math::atlas::cblas_nrm2<nm::RubyObject,nm::RubyObject>
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
static VALUE nm_atlas_cblas_asum(VALUE self, VALUE n, VALUE x, VALUE incx) {

  static void (*ttable[nm::NUM_DTYPES])(const int N, const void* X, const int incX, void* sum) = {
      nm::math::atlas::cblas_asum<uint8_t,uint8_t>,
      nm::math::atlas::cblas_asum<int8_t,int8_t>,
      nm::math::atlas::cblas_asum<int16_t,int16_t>,
      nm::math::atlas::cblas_asum<int32_t,int32_t>,
      nm::math::atlas::cblas_asum<int64_t,int64_t>,
      nm::math::atlas::cblas_asum<float32_t,float32_t>,
      nm::math::atlas::cblas_asum<float64_t,float64_t>,
      nm::math::atlas::cblas_asum<float32_t,nm::Complex64>,
      nm::math::atlas::cblas_asum<float64_t,nm::Complex128>,
      nm::math::atlas::cblas_asum<nm::RubyObject,nm::RubyObject>
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
static VALUE nm_atlas_cblas_imax(VALUE self, VALUE n, VALUE x, VALUE incx) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::atlas::cblas_imax, int, const int n, const void* x, const int incx);

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
static VALUE nm_atlas_cblas_gemv(VALUE self,
                           VALUE trans_a,
                           VALUE m, VALUE n,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE x, VALUE incx,
                           VALUE beta,
                           VALUE y, VALUE incy)
{
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::atlas::cblas_gemv, bool, const enum CBLAS_TRANSPOSE, const int, const int, const void*, const void*, const int, const void*, const int, const void*, void*, const int)

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
static VALUE nm_atlas_cblas_gemm(VALUE self,
                           VALUE order,
                           VALUE trans_a, VALUE trans_b,
                           VALUE m, VALUE n, VALUE k,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb,
                           VALUE beta,
                           VALUE c, VALUE ldc)
{
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::atlas::cblas_gemm, void, const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_TRANSPOSE trans_b, int m, int n, int k, void* alpha, void* a, int lda, void* b, int ldb, void* beta, void* c, int ldc);

  nm::dtype_t dtype = NM_DTYPE(a);

  void *pAlpha = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]),
       *pBeta  = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, pAlpha);
  rubyval_to_cval(beta, dtype, pBeta);

  ttable[dtype](blas_order_sym(order), blas_transpose_sym(trans_a), blas_transpose_sym(trans_b), FIX2INT(m), FIX2INT(n), FIX2INT(k), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb), pBeta, NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc));

  return c;
}


static VALUE nm_atlas_cblas_trsm(VALUE self,
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
      nm::math::atlas::cblas_trsm<float>,
      nm::math::atlas::cblas_trsm<double>,
      cblas_ctrsm, cblas_ztrsm, // call directly, same function signature!
      nm::math::atlas::cblas_trsm<nm::RubyObject>
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

static VALUE nm_atlas_cblas_trmm(VALUE self,
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
      nm::math::atlas::cblas_trmm<float>,
      nm::math::atlas::cblas_trmm<double>,
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

static VALUE nm_atlas_cblas_syrk(VALUE self,
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
      nm::math::atlas::cblas_syrk<float>,
      nm::math::atlas::cblas_syrk<double>,
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

static VALUE nm_atlas_cblas_herk(VALUE self,
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

/*
 * Function signature conversion for calling CBLAS' gesvd functions as directly as possible.
 *
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
static VALUE nm_atlas_lapack_gesvd(VALUE self, VALUE jobu, VALUE jobvt, VALUE m, VALUE n, VALUE a, VALUE lda, VALUE s, VALUE u, VALUE ldu, VALUE vt, VALUE ldvt, VALUE lwork) {
  static int (*gesvd_table[nm::NUM_DTYPES])(char, char, int, int, void* a, int, void* s, void* u, int, void* vt, int, void* work, int, void* rwork) = {
    NULL, NULL, NULL, NULL, NULL, // no integer ops
    nm::math::atlas::lapack_gesvd<float,float>,
    nm::math::atlas::lapack_gesvd<double,double>,
    nm::math::atlas::lapack_gesvd<nm::Complex64,float>,
    nm::math::atlas::lapack_gesvd<nm::Complex128,double>,
    NULL // no Ruby objects
  };

  nm::dtype_t dtype = NM_DTYPE(a);


  if (!gesvd_table[dtype]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
    return Qfalse;
  } else {
    int M = FIX2INT(m),
        N = FIX2INT(n);

    int min_mn  = NM_MIN(M,N);
    int max_mn  = NM_MAX(M,N);

    char JOBU = lapack_svd_job_sym(jobu),
         JOBVT = lapack_svd_job_sym(jobvt);

    // only need rwork for complex matrices
    int rwork_size  = (dtype == nm::COMPLEX64 || dtype == nm::COMPLEX128) ? 5 * min_mn : 0;
    void* rwork     = rwork_size > 0 ? NM_ALLOCA_N(char, DTYPE_SIZES[dtype] * rwork_size) : NULL;
    int work_size   = FIX2INT(lwork);

    // ignore user argument for lwork if it's too small.
    work_size       = NM_MAX((dtype == nm::COMPLEX64 || dtype == nm::COMPLEX128 ? 2 * min_mn + max_mn : NM_MAX(3*min_mn + max_mn, 5*min_mn)), work_size);
    void* work      = NM_ALLOCA_N(char, DTYPE_SIZES[dtype] * work_size);

    int info = gesvd_table[dtype](JOBU, JOBVT, M, N, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda),
      NM_STORAGE_DENSE(s)->elements, NM_STORAGE_DENSE(u)->elements, FIX2INT(ldu), NM_STORAGE_DENSE(vt)->elements, FIX2INT(ldvt),
      work, work_size, rwork);
    return INT2FIX(info);
  }
}

/*
 * Function signature conversion for calling CBLAS' gesdd functions as directly as possible.
 *
 * xGESDD uses a divide-and-conquer strategy to compute the singular value decomposition (SVD) of a real
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
static VALUE nm_atlas_lapack_gesdd(VALUE self, VALUE jobz, VALUE m, VALUE n, VALUE a, VALUE lda, VALUE s, VALUE u, VALUE ldu, VALUE vt, VALUE ldvt, VALUE lwork) {
  static int (*gesdd_table[nm::NUM_DTYPES])(char, int, int, void* a, int, void* s, void* u, int, void* vt, int, void* work, int, int* iwork, void* rwork) = {
    NULL, NULL, NULL, NULL, NULL, // no integer ops
    nm::math::atlas::lapack_gesdd<float,float>,
    nm::math::atlas::lapack_gesdd<double,double>,
    nm::math::atlas::lapack_gesdd<nm::Complex64,float>,
    nm::math::atlas::lapack_gesdd<nm::Complex128,double>,
    NULL // no Ruby objects
  };

  nm::dtype_t dtype = NM_DTYPE(a);

  if (!gesdd_table[dtype]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
    return Qfalse;
  } else {
    int M = FIX2INT(m),
        N = FIX2INT(n);

    int min_mn  = NM_MIN(M,N);
    int max_mn  = NM_MAX(M,N);

    char JOBZ = lapack_svd_job_sym(jobz);

    // only need rwork for complex matrices
    void* rwork = NULL;

    int work_size = FIX2INT(lwork); // Make sure we allocate enough work, regardless of the user request.
    if (dtype == nm::COMPLEX64 || dtype == nm::COMPLEX128) {
      int rwork_size = min_mn * (JOBZ == 'N' ? 5 : NM_MAX(5*min_mn + 7, 2*max_mn + 2*min_mn + 1));
      rwork = NM_ALLOCA_N(char, DTYPE_SIZES[dtype] * rwork_size);

      if (JOBZ == 'N')      work_size = NM_MAX(work_size, 3*min_mn + NM_MAX(max_mn, 6*min_mn));
      else if (JOBZ == 'O') work_size = NM_MAX(work_size, 3*min_mn*min_mn + NM_MAX(max_mn, 5*min_mn*min_mn + 4*min_mn));
      else                  work_size = NM_MAX(work_size, 3*min_mn*min_mn + NM_MAX(max_mn, 4*min_mn*min_mn + 4*min_mn));
    } else {
      if (JOBZ == 'N')      work_size = NM_MAX(work_size, 2*min_mn + max_mn);
      else if (JOBZ == 'O') work_size = NM_MAX(work_size, 2*min_mn*min_mn + max_mn + 2*min_mn);
      else                  work_size = NM_MAX(work_size, min_mn*min_mn + max_mn + 2*min_mn);
    }
    void* work  = NM_ALLOCA_N(char, DTYPE_SIZES[dtype] * work_size);
    int* iwork  = NM_ALLOCA_N(int, 8*min_mn);

    int info = gesdd_table[dtype](JOBZ, M, N, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda),
      NM_STORAGE_DENSE(s)->elements, NM_STORAGE_DENSE(u)->elements, FIX2INT(ldu), NM_STORAGE_DENSE(vt)->elements, FIX2INT(ldvt),
      work, work_size, iwork, rwork);
    return INT2FIX(info);
  }
}

/*
 * Function signature conversion for calling CBLAS' geev functions as directly as possible.
 *
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
static VALUE nm_atlas_lapack_geev(VALUE self, VALUE compute_left, VALUE compute_right, VALUE n, VALUE a, VALUE lda, VALUE w, VALUE wi, VALUE vl, VALUE ldvl, VALUE vr, VALUE ldvr, VALUE lwork) {
  static int (*geev_table[nm::NUM_DTYPES])(char, char, int, void* a, int, void* w, void* wi, void* vl, int, void* vr, int, void* work, int, void* rwork) = {
    NULL, NULL, NULL, NULL, NULL, // no integer ops
    nm::math::atlas::lapack_geev<float,float>,
    nm::math::atlas::lapack_geev<double,double>,
    nm::math::atlas::lapack_geev<nm::Complex64,float>,
    nm::math::atlas::lapack_geev<nm::Complex128,double>,
    NULL // no Ruby objects
  };

  nm::dtype_t dtype = NM_DTYPE(a);


  if (!geev_table[dtype]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
    return Qfalse;
  } else {
    int N = FIX2INT(n);

    char JOBVL = lapack_evd_job_sym(compute_left),
         JOBVR = lapack_evd_job_sym(compute_right);

    void* A  = NM_STORAGE_DENSE(a)->elements;
    void* WR = NM_STORAGE_DENSE(w)->elements;
    void* WI = wi == Qnil ? NULL : NM_STORAGE_DENSE(wi)->elements;
    void* VL = JOBVL == 'V' ? NM_STORAGE_DENSE(vl)->elements : NULL;
    void* VR = JOBVR == 'V' ? NM_STORAGE_DENSE(vr)->elements : NULL;

    // only need rwork for complex matrices (wi == Qnil for complex)
    int rwork_size  = dtype == nm::COMPLEX64 || dtype == nm::COMPLEX128 ? N * DTYPE_SIZES[dtype] : 0; // 2*N*floattype for complex only, otherwise 0
    void* rwork     = rwork_size > 0 ? NM_ALLOCA_N(char, rwork_size) : NULL;
    int work_size   = FIX2INT(lwork);
    void* work;

    int info;

    // if work size is 0 or -1, query.
    if (work_size <= 0) {
      work_size = -1;
      work = NM_ALLOC_N(char, DTYPE_SIZES[dtype]); //2*N * DTYPE_SIZES[dtype]);
      info = geev_table[dtype](JOBVL, JOBVR, N, A, FIX2INT(lda), WR, WI, VL, FIX2INT(ldvl), VR, FIX2INT(ldvr), work, work_size, rwork);
      work_size = (int)(dtype == nm::COMPLEX64 || dtype == nm::FLOAT32 ? reinterpret_cast<float*>(work)[0] : reinterpret_cast<double*>(work)[0]);
      // line above is basically: work_size = (int)(work[0]); // now have new work_size
      NM_FREE(work);
      if (info == 0)
        rb_warn("geev: calculated optimal lwork of %d; to eliminate this message, use a positive value for lwork (at least 2*shape[i])", work_size);
      else return INT2FIX(info); // error of some kind on query!
    }

    // if work size is < 2*N, just set it to 2*N
    if (work_size < 2*N) work_size = 2*N;
    if (work_size < 3*N && (dtype == nm::FLOAT32 || dtype == nm::FLOAT64)) {
      work_size = JOBVL == 'V' || JOBVR == 'V' ? 4*N : 3*N;
    }

    // Allocate work array for actual run
    work = NM_ALLOCA_N(char, work_size * DTYPE_SIZES[dtype]);

    // Perform the actual calculation.
    info = geev_table[dtype](JOBVL, JOBVR, N, A, FIX2INT(lda), WR, WI, VL, FIX2INT(ldvl), VR, FIX2INT(ldvr), work, work_size, rwork);

    return INT2FIX(info);
  }
}

/* Call any of the clapack_xgetrf functions as directly as possible.
 *
 * The clapack_getrf functions (dgetrf, sgetrf, cgetrf, and zgetrf) compute an LU factorization of a general M-by-N
 * matrix A using partial pivoting with row interchanges.
 *
 * The factorization has the form:
 *    A = P * L * U
 * where P is a permutation matrix, L is lower triangular with unit diagonal elements (lower trapezoidal if m > n),
 * and U is upper triangular (upper trapezoidal if m < n).
 *
 * This is the right-looking level 3 BLAS version of the algorithm.
 *
 * == Arguments
 * See: http://www.netlib.org/lapack/double/dgetrf.f
 * (You don't need argument 5; this is the value returned by this function.)
 *
 * You probably don't want to call this function. Instead, why don't you try clapack_getrf, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 *
 * Returns an array giving the pivot indices (normally these are argument #5).
 */
static VALUE nm_atlas_clapack_getrf(VALUE self, VALUE order, VALUE m, VALUE n, VALUE a, VALUE lda) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const int m, const int n, void* a, const int lda, int* ipiv) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::atlas::clapack_getrf<float>,
      nm::math::atlas::clapack_getrf<double>,
#if defined (HAVE_CLAPACK_H) || defined (HAVE_ATLAS_CLAPACK_H)
      clapack_cgetrf, clapack_zgetrf, // call directly, same function signature!
#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::atlas::clapack_getrf<nm::Complex64>,
      nm::math::atlas::clapack_getrf<nm::Complex128>,
#endif
      nm::math::atlas::clapack_getrf<nm::RubyObject>
  };

  int M = FIX2INT(m),
      N = FIX2INT(n);

  // Allocate the pivot index array, which is of size MIN(M, N).
  size_t ipiv_size = std::min(M,N);
  int* ipiv = NM_ALLOCA_N(int, ipiv_size);

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    // Call either our version of getrf or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), M, N, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), ipiv);
  }

  // Result will be stored in a. We return ipiv as an array.
  VALUE ipiv_array = rb_ary_new2(ipiv_size);
  for (size_t i = 0; i < ipiv_size; ++i) {
    rb_ary_store(ipiv_array, i, INT2FIX(ipiv[i]));
  }

  return ipiv_array;
}


/* Call any of the clapack_xpotrf functions as directly as possible.
 *
 * You probably don't want to call this function. Instead, why don't you try clapack_potrf, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_atlas_clapack_potrf(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda) {
#if !defined(HAVE_CLAPACK_H) && !defined(HAVE_ATLAS_CLAPACK_H)
  rb_raise(rb_eNotImpError, "potrf currently requires CLAPACK");
#endif

  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const enum CBLAS_UPLO, const int n, void* a, const int lda) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::atlas::clapack_potrf<float>,
      nm::math::atlas::clapack_potrf<double>,
#if defined (HAVE_CLAPACK_H) || defined (HAVE_ATLAS_CLAPACK_H)
      clapack_cpotrf, clapack_zpotrf, // call directly, same function signature!
#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::atlas::clapack_potrf<nm::Complex64>,
      nm::math::atlas::clapack_potrf<nm::Complex128>,
#endif
      NULL
  };

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
    // FIXME: Once BLAS dtypes are implemented, replace error above with the error below.
    //rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    // Call either our version of potrf or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), blas_uplo_sym(uplo), FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda));
  }

  return a;
}


/*
 * Call any of the clapack_xgetrs functions as directly as possible.
 */
static VALUE nm_atlas_clapack_getrs(VALUE self, VALUE order, VALUE trans, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE ipiv, VALUE b, VALUE ldb) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans, const int N,
                                       const int NRHS, const void* A, const int lda, const int* ipiv, void* B,
                                       const int ldb) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::atlas::clapack_getrs<float>,
      nm::math::atlas::clapack_getrs<double>,
#if defined (HAVE_CLAPACK_H) || defined (HAVE_ATLAS_CLAPACK_H)
      clapack_cgetrs, clapack_zgetrs, // call directly, same function signature!
#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::atlas::clapack_getrs<nm::Complex64>,
      nm::math::atlas::clapack_getrs<nm::Complex128>,
#endif
      nm::math::atlas::clapack_getrs<nm::RubyObject>
  };

  // Allocate the C version of the pivot index array
  int* ipiv_;
  if (TYPE(ipiv) != T_ARRAY) {
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

    // Call either our version of getrs or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), blas_transpose_sym(trans), FIX2INT(n), FIX2INT(nrhs), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda),
                        ipiv_, NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb));
  }

  // b is both returned and modified directly in the argument list.
  return b;
}


/*
 * Call any of the clapack_xpotrs functions as directly as possible.
 */
static VALUE nm_atlas_clapack_potrs(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE b, VALUE ldb) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const int N,
                                       const int NRHS, const void* A, const int lda, void* B, const int ldb) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::atlas::clapack_potrs<float>,
      nm::math::atlas::clapack_potrs<double>,
#if defined (HAVE_CLAPACK_H) || defined (HAVE_ATLAS_CLAPACK_H)
      clapack_cpotrs, clapack_zpotrs, // call directly, same function signature!
#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::atlas::clapack_potrs<nm::Complex64>,
      nm::math::atlas::clapack_potrs<nm::Complex128>,
#endif
      nm::math::atlas::clapack_potrs<nm::RubyObject>
  };


  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {

    // Call either our version of potrs or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), blas_uplo_sym(uplo), FIX2INT(n), FIX2INT(nrhs), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda),
                        NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb));
  }

  // b is both returned and modified directly in the argument list.
  return b;
}

/* Call any of the clapack_xgetri functions as directly as possible.
 *
 * You probably don't want to call this function. Instead, why don't you try clapack_getri, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 *
 * Returns an array giving the pivot indices (normally these are argument #5).
 */
static VALUE nm_atlas_clapack_getri(VALUE self, VALUE order, VALUE n, VALUE a, VALUE lda, VALUE ipiv) {
#if !defined (HAVE_CLAPACK_H) && !defined (HAVE_ATLAS_CLAPACK_H)
  rb_raise(rb_eNotImpError, "getri currently requires CLAPACK");
#endif

  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const int n, void* a, const int lda, const int* ipiv) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::atlas::clapack_getri<float>,
      nm::math::atlas::clapack_getri<double>,
#if defined (HAVE_CLAPACK_H) || defined (HAVE_ATLAS_CLAPACK_H)
      clapack_cgetri, clapack_zgetri, // call directly, same function signature!
#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::atlas::clapack_getri<nm::Complex64>,
      nm::math::atlas::clapack_getri<nm::Complex128>,
#endif
      NULL
  };

  // Allocate the C version of the pivot index array
  int* ipiv_;
  if (TYPE(ipiv) != T_ARRAY) {
    rb_raise(rb_eArgError, "ipiv must be of type Array");
  } else {
    ipiv_ = NM_ALLOCA_N(int, RARRAY_LEN(ipiv));
    for (int index = 0; index < RARRAY_LEN(ipiv); ++index) {
      ipiv_[index] = FIX2INT( RARRAY_AREF(ipiv, index) );
    }
  }

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
    // FIXME: Once non-BLAS dtypes are implemented, replace error above with the error below.
    //rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    // Call either our version of getri or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), ipiv_);
  }

  return a;
}


/* Call any of the clapack_xpotri functions as directly as possible.
 *
 * You probably don't want to call this function. Instead, why don't you try clapack_potri, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_atlas_clapack_potri(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda) {
#if !defined (HAVE_CLAPACK_H) && !defined (HAVE_ATLAS_CLAPACK_H)
  rb_raise(rb_eNotImpError, "getri currently requires CLAPACK");
#endif

  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const enum CBLAS_UPLO, const int n, void* a, const int lda) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::atlas::clapack_potri<float>,
      nm::math::atlas::clapack_potri<double>,
#if defined (HAVE_CLAPACK_H) || defined (HAVE_ATLAS_CLAPACK_H)
      clapack_cpotri, clapack_zpotri, // call directly, same function signature!
#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::atlas::clapack_potri<nm::Complex64>,
      nm::math::atlas::clapack_potri<nm::Complex128>,
#endif
      NULL
  };

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
    // FIXME: Once BLAS dtypes are implemented, replace error above with the error below.
    //rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    // Call either our version of getri or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), blas_uplo_sym(uplo), FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda));
  }

  return a;
}


/*
 * Call any of the clapack_xlaswp functions as directly as possible.
 *
 * Note that LAPACK's xlaswp functions accept a column-order matrix, but NMatrix uses row-order. Thus, n should be the
 * number of rows and lda should be the number of columns, no matter what it says in the documentation for dlaswp.f.
 */
static VALUE nm_atlas_clapack_laswp(VALUE self, VALUE n, VALUE a, VALUE lda, VALUE k1, VALUE k2, VALUE ipiv, VALUE incx) {
  //We have actually never used the ATLAS version of laswp. For the time being
  //I will leave it like that and just always call the internal implementation.
  //I don't know if there is a good reason for this or not.
  //Maybe because our internal version swaps columns instead of rows.
  static void (*ttable[nm::NUM_DTYPES])(const int n, void* a, const int lda, const int k1, const int k2, const int* ipiv, const int incx) = {
      nm::math::clapack_laswp<uint8_t>,
      nm::math::clapack_laswp<int8_t>,
      nm::math::clapack_laswp<int16_t>,
      nm::math::clapack_laswp<int32_t>,
      nm::math::clapack_laswp<int64_t>,
      nm::math::clapack_laswp<float>,
      nm::math::clapack_laswp<double>,
      nm::math::clapack_laswp<nm::Complex64>,
      nm::math::clapack_laswp<nm::Complex128>,
      nm::math::clapack_laswp<nm::RubyObject>
  };

  // Allocate the C version of the pivot index array
  int* ipiv_;
  if (TYPE(ipiv) != T_ARRAY) {
    rb_raise(rb_eArgError, "ipiv must be of type Array");
  } else {
    ipiv_ = NM_ALLOCA_N(int, RARRAY_LEN(ipiv));
    for (int index = 0; index < RARRAY_LEN(ipiv); ++index) {
      ipiv_[index] = FIX2INT( RARRAY_AREF(ipiv, index) );
    }
  }

  // Call either our version of laswp or the LAPACK version.
  ttable[NM_DTYPE(a)](FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), FIX2INT(k1), FIX2INT(k2), ipiv_, FIX2INT(incx));

  // a is both returned and modified directly in the argument list.
  return a;
}


}
