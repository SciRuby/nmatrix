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
// == ruby_symbols.cpp
//
// Ruby symbols used throught the NMatrix project.

/*
 * Standard Includes
 */

#include <ruby.h>

/*
 * Global Variables
 */

ID  nm_rb_dtype,
    nm_rb_stype,

    nm_rb_capacity,
    nm_rb_default,

    nm_rb_real,
    nm_rb_imag,

    nm_rb_numer,
    nm_rb_denom,

    nm_rb_complex_conjugate,
    nm_rb_transpose,
    nm_rb_no_transpose,
    nm_rb_left,
    nm_rb_right,
    nm_rb_upper,
    nm_rb_lower,
    nm_rb_unit,
    nm_rb_nonunit,

    nm_rb_dense,
    nm_rb_list,
    nm_rb_yale,

    nm_rb_row,
    nm_rb_column,
    nm_rb_add,
    nm_rb_sub,
    nm_rb_mul,
    nm_rb_div,
    nm_rb_both,
    nm_rb_none,

    nm_rb_negate,

    nm_rb_percent,
    nm_rb_gt,
    nm_rb_lt,
    nm_rb_eql,
    nm_rb_neql,
    nm_rb_gte,
    nm_rb_lte,

    nm_rb_hash;

VALUE cNMatrix,
      cNMatrix_IO,
      cNMatrix_IO_Matlab,
      cNMatrix_YaleFunctions,

      cNMatrix_GC_holder,

      nm_eDataTypeError,
      nm_eConvergenceError,
      nm_eStorageTypeError,
      nm_eShapeError,
      nm_eNotInvertibleError;

/*
 * Functions
 */

void nm_init_ruby_constants(void) {
  nm_rb_dtype             = rb_intern("dtype");
  nm_rb_stype             = rb_intern("stype");

  nm_rb_capacity          = rb_intern("capacity");
  nm_rb_default           = rb_intern("default");

  nm_rb_real              = rb_intern("real");
  nm_rb_imag              = rb_intern("imag");

  nm_rb_numer              = rb_intern("numerator");
  nm_rb_denom              = rb_intern("denominator");

  nm_rb_complex_conjugate  = rb_intern("complex_conjugate");
  nm_rb_transpose          = rb_intern("transpose");
  nm_rb_no_transpose      = rb_intern("no_transpose");

  nm_rb_dense             = rb_intern("dense");
  nm_rb_list              = rb_intern("list");
  nm_rb_yale              = rb_intern("yale");

  nm_rb_add                = rb_intern("+");
  nm_rb_sub                = rb_intern("-");
  nm_rb_mul                = rb_intern("*");
  nm_rb_div                = rb_intern("/");

  nm_rb_negate            = rb_intern("-@");

  nm_rb_percent            = rb_intern("%");
  nm_rb_gt                = rb_intern(">");
  nm_rb_lt                = rb_intern("<");
  nm_rb_eql                = rb_intern("==");
  nm_rb_neql              = rb_intern("!=");
  nm_rb_gte                = rb_intern(">=");
  nm_rb_lte                = rb_intern("<=");

  nm_rb_left              = rb_intern("left");
  nm_rb_right             = rb_intern("right");
  nm_rb_upper             = rb_intern("upper");
  nm_rb_lower             = rb_intern("lower");
  nm_rb_unit              = rb_intern("unit");
  nm_rb_nonunit           = rb_intern("nonunit");
  nm_rb_hash              = rb_intern("hash");

  nm_rb_column            = rb_intern("column");
  nm_rb_row               = rb_intern("row");

  nm_rb_both              = rb_intern("both");
  nm_rb_none              = rb_intern("none");
}
