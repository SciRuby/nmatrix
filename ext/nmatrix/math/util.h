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
// == util.h
//
// Collect a few utility functions which convert ruby symbols into arguments
// that CBLAS or LAPACK can understand: either enum's for CBLAS or char's
// for LAPACK.
//

#ifndef UTIL_H
#define UTIL_H

/* Interprets cblas argument which could be any of false/:no_transpose, :transpose, or :complex_conjugate,
 * into an enum recognized by cblas.
 *
 * Called by nm_cblas_gemm -- basically inline.
 *
 */
static inline enum CBLAS_TRANSPOSE blas_transpose_sym(VALUE op) {
  if (op == Qfalse || rb_to_id(op) == nm_rb_no_transpose) return CblasNoTrans;
  else if (rb_to_id(op) == nm_rb_transpose) return CblasTrans;
  else if (rb_to_id(op) == nm_rb_complex_conjugate) return CblasConjTrans;
  else rb_raise(rb_eArgError, "Expected false, :transpose, or :complex_conjugate");
  return CblasNoTrans;
}

/* Interprets transpose argument which could be any of false/:no_transpose, :transpose, or :complex_conjugate,
 * into an character recognized by LAPACKE. LAPACKE uses a different system than CBLAS for this.
 *
 */
static inline char lapacke_transpose_sym(VALUE op) {
  if (op == Qfalse || rb_to_id(op) == nm_rb_no_transpose) return 'N';
  else if (rb_to_id(op) == nm_rb_transpose) return 'T';
  else if (rb_to_id(op) == nm_rb_complex_conjugate) return 'C';
  else rb_raise(rb_eArgError, "Expected false, :transpose, or :complex_conjugate");
  return 'N';
}

/*
 * Interprets cblas argument which could be :left or :right
 *
 * Called by nm_cblas_trsm -- basically inline
 */
static inline enum CBLAS_SIDE blas_side_sym(VALUE op) {
  ID op_id = rb_to_id(op);
  if (op_id == nm_rb_left)  return CblasLeft;
  if (op_id == nm_rb_right) return CblasRight;
  rb_raise(rb_eArgError, "Expected :left or :right for side argument");
  return CblasLeft;
}

/*
 * Interprets cblas argument which could be :upper or :lower
 *
 * Called by nm_cblas_trsm -- basically inline
 */
static inline enum CBLAS_UPLO blas_uplo_sym(VALUE op) {
  ID op_id = rb_to_id(op);
  if (op_id == nm_rb_upper) return CblasUpper;
  if (op_id == nm_rb_lower) return CblasLower;
  rb_raise(rb_eArgError, "Expected :upper or :lower for uplo argument");
  return CblasUpper;
}

/*
 * Interprets argument which could be :upper or :lower for LAPACKE
 *
 * Called by nm_cblas_trsm -- basically inline
 */
static inline char lapacke_uplo_sym(VALUE op) {
  ID op_id = rb_to_id(op);
  if (op_id == nm_rb_upper) return 'U';
  if (op_id == nm_rb_lower) return 'L';
  rb_raise(rb_eArgError, "Expected :upper or :lower for uplo argument");
  return 'U';
}

/*
 * Interprets cblas argument which could be :unit (true) or :nonunit (false or anything other than true/:unit)
 *
 * Called by nm_cblas_trsm -- basically inline
 */
static inline enum CBLAS_DIAG blas_diag_sym(VALUE op) {
  if (rb_to_id(op) == nm_rb_unit || op == Qtrue) return CblasUnit;
  return CblasNonUnit;
}

/*
 * Interprets cblas argument which could be :row or :col
 * 
 * This function, unlike the other ones, works for LAPACKE as well as for CBLAS/CLAPACK.
 * Although LAPACKE calls this an int instead of a enum, the magic values are the same
 * (101 for row-major, 102 for column-major).
 */
static inline enum CBLAS_ORDER blas_order_sym(VALUE op) {
  if (rb_to_id(op) == rb_intern("row") || rb_to_id(op) == rb_intern("row_major")) return CblasRowMajor;
  else if (rb_to_id(op) == rb_intern("col") || rb_to_id(op) == rb_intern("col_major") ||
           rb_to_id(op) == rb_intern("column") || rb_to_id(op) == rb_intern("column_major")) return CblasColMajor;
  rb_raise(rb_eArgError, "Expected :row or :col for order argument");
  return CblasRowMajor;
}

/*
 * Interprets lapack jobu and jobvt arguments, for which LAPACK needs character values A, S, O, or N.
 *
 * Called by lapack_gesvd -- basically inline. svd stands for singular value decomposition.
 */
static inline char lapack_svd_job_sym(VALUE op) {
  if (rb_to_id(op) == rb_intern("all") || rb_to_id(op) == rb_intern("a")) return 'A';
  else if (rb_to_id(op) == rb_intern("return") || rb_to_id(op) == rb_intern("s")) return 'S';
  else if (rb_to_id(op) == rb_intern("overwrite") || rb_to_id(op) == rb_intern("o")) return 'O';
  else if (rb_to_id(op) == rb_intern("none") || rb_to_id(op) == rb_intern("n")) return 'N';
  else rb_raise(rb_eArgError, "Expected :all, :return, :overwrite, :none (or :a, :s, :o, :n, respectively)");
  return 'a';
}

/*
 * Interprets lapack jobvl and jobvr arguments, for which LAPACK needs character values N or V.
 *
 * Called by lapack_geev -- basically inline. evd stands for eigenvalue decomposition.
 */
static inline char lapack_evd_job_sym(VALUE op) {
  if (op == Qfalse || op == Qnil || rb_to_id(op) == rb_intern("n")) return 'N';
  else return 'V';
}

#endif
