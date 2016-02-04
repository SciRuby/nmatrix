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
// == nmatrix_fftw.cpp
//
// Main file for nmatrix_fftw extension
//

#include <ruby.h>
#include <complex.h>
#include <fftw3.h>
#include "storage/common.h"
#include "nmatrix.h"
#include <iostream>

#define TYPE_COMPLEX_COMPLEX 0
#define TYPE_REAL_COMPLEX    1
#define TYPE_COMPLEX_REAL    2
#define TYPE_REAL_REAL       3

// @private Used internally by the C API.
static VALUE cNMatrix_FFTW_Plan_Data;

// @private Used internally by the C API.
//
// ADT for encapsulating various data structures required for sucessfully planning
//   and executing a fourier transform with FFTW. Uses void* pointers because 
//   input/output can be either double or fftw_complex depending on the type of
//   FFT being planned.
struct fftw_data {
  void* input; 
  void* output;
  fftw_plan plan;
};

// @private Used internally by the C API.
// Method used by Ruby GC for freeing memory allocated by FFTW.
static void nm_fftw_cleanup(fftw_data* d)
{
  xfree(d->input);
  xfree(d->output);
  fftw_destroy_plan(d->plan);
  xfree(d);
}

// @private Used internally by the C API.
// Used for converting a Ruby Array containing the shape to a C++ array of ints.
static int* nm_fftw_interpret_shape(VALUE rb_shape, const int dimensions)
{
  Check_Type(rb_shape, T_ARRAY);

  int *shape = new int[dimensions];
  const VALUE *arr = RARRAY_CONST_PTR(rb_shape);

  for (int i = 0; i < dimensions; ++i) {
    shape[i] = FIX2INT(arr[i]);
  }

  return shape;
}

// @private Used internally by the C API.
// Convert values passed in Ruby Array containing kinds of real-real transforms 
//   to a C array of ints. 
static void
nm_fftw_interpret_real_real_kind(VALUE real_real_kind, int *r2r_kinds)
{
  int size = RARRAY_LEN(real_real_kind);
  const VALUE *a = RARRAY_CONST_PTR(real_real_kind);
  for (int i = 0; i < size; ++i) { 
    r2r_kinds[i] = FIX2INT(a[i]); 
  }
}

// @private Used internally by the C API.
// Actually calls the FFTW planner routines based on the input/output and the
//   type of routine selected. Also allocates memory for input and output pointers.
static void nm_fftw_actually_create_plan(fftw_data* data, 
  size_t size, const int dimensions, const int* shape, int sign, unsigned flags, 
  VALUE rb_type, VALUE real_real_kind)
{
  switch (FIX2INT(rb_type))
  {
    case TYPE_COMPLEX_COMPLEX:
      data->input  = ALLOC_N(fftw_complex, size);
      data->output = ALLOC_N(fftw_complex, size);
      data->plan   = fftw_plan_dft(dimensions, shape, (fftw_complex*)data->input, 
        (fftw_complex*)data->output, sign, flags);
      break;
    case TYPE_REAL_COMPLEX:
      data->input  = ALLOC_N(double      , size);
      data->output = ALLOC_N(fftw_complex, size);
      data->plan   = fftw_plan_dft_r2c(dimensions, shape, (double*)data->input, 
        (fftw_complex*)data->output, flags);
      break;
    case TYPE_COMPLEX_REAL:
      data->input  = ALLOC_N(fftw_complex,  size);
      data->output = ALLOC_N(double      ,  size);
      data->plan   = fftw_plan_dft_c2r(dimensions, shape, (fftw_complex*)data->input, 
        (double*)data->output, flags);
      break;
    case TYPE_REAL_REAL:
      int* r2r_kinds = ALLOC_N(int, FIX2INT(real_real_kind));
      nm_fftw_interpret_real_real_kind(real_real_kind, r2r_kinds);
      data->input  = ALLOC_N(double, size);
      data->output = ALLOC_N(double, size);
      data->plan   = fftw_plan_r2r(dimensions, shape, (double*)data->input, 
        (double*)data->output, (fftw_r2r_kind*)r2r_kinds, flags);
      xfree(r2r_kinds);
      break;
  }
}

/** \brief Create a plan for performing the fourier transform based on input,
 * output pointers and the underlying hardware.
 *
 * @param[in] self          Object on which the function is called
 * @param[in] rb_shape      Shape of the plan.
 * @param[in] rb_size       Size of the plan.
 * @param[in] rb_dim        Dimension of the FFT to be performed.
 * @param[in] rb_flags      Number denoting the planner flags.
 * @param[in] rb_direction  Direction of FFT (can be -1 or +1). Specifies the
 *   sign of the exponent.
 * @param[in] rb_type       Number specifying the type of FFT being planned (one
 *    of :complex_complex, :complex_real, :real_complex and :real_real)
 * @param[in] rb_real_real_kind    Ruby Array specifying the kind of DFT to perform over
 *   each axis in case of a real input/real output FFT.
 *
 * \returns An object of type NMatrix::FFTW::Plan::Data that encapsulates the
 * plan and relevant input/output arrays.
 */
static VALUE nm_fftw_create_plan(VALUE self, VALUE rb_shape, VALUE rb_size,
  VALUE rb_dim, VALUE rb_flags, VALUE rb_direction, VALUE rb_type, VALUE rb_real_real_kind)
{ 
  const int dimensions = FIX2INT(rb_dim);
  const int* shape     = nm_fftw_interpret_shape(rb_shape, dimensions);
  size_t size          = FIX2INT(rb_size);
  int sign             = FIX2INT(rb_direction);
  unsigned flags       = FIX2INT(rb_flags);
  fftw_data *data      = ALLOC(fftw_data);

  nm_fftw_actually_create_plan(data, size, dimensions, shape, 
    sign, flags, rb_type, rb_real_real_kind);
  
  return Data_Wrap_Struct(cNMatrix_FFTW_Plan_Data, NULL, nm_fftw_cleanup, data);
}

// @private Used internally by the C API.
template <typename InputType>
static void nm_fftw_actually_set(VALUE nmatrix, VALUE plan_data)
{
  fftw_data* data;
  Data_Get_Struct(plan_data, fftw_data, data);
  memcpy((InputType*)data->input, (InputType*)NM_DENSE_ELEMENTS(nmatrix), 
    sizeof(InputType)*NM_DENSE_COUNT(nmatrix));
}

/** \brief Here is a brief description of what this function does.
 *
 * @param[in,out] self       Object on which the function is called.
 * @param[in]     plan_data  An internal data structure of type 
 *   NMatrix::FFTW::Plan::Data that is created by Data_Wrap_Struct in 
 *   nm_fftw_create_plan and which encapsulates the FFTW plan in a Ruby object.
 * @param[in]     nmatrix    An NMatrix object (pre-allocated) which contains the
 *   input elements for the fourier transform.
 * @param[in]     type       A number representing the type of fourier transform 
 *   being performed. (:complex_complex, :real_complex, :complex_real or :real_real).
 *
 * \returns self
 */
static VALUE nm_fftw_set_input(VALUE self, VALUE nmatrix, VALUE plan_data, 
  VALUE type)
{
  switch(FIX2INT(type))
  {
    case TYPE_COMPLEX_COMPLEX:
    case TYPE_COMPLEX_REAL:
      nm_fftw_actually_set<fftw_complex>(nmatrix, plan_data);
      break;
    case TYPE_REAL_COMPLEX:
    case TYPE_REAL_REAL:
      nm_fftw_actually_set<double>(nmatrix, plan_data);
      break;
    default:
      rb_raise(rb_eArgError, "Invalid type of DFT.");
  }

  return self;
}

// @private Used internally by the C API.
// Call fftw_execute and copy the resulting data into the nmatrix object.
template <typename OutputType>
static void nm_fftw_actually_execute(VALUE nmatrix, VALUE plan_data)
{
  fftw_data *data;
  Data_Get_Struct(plan_data, fftw_data, data);
  fftw_execute(data->plan);
  memcpy((OutputType*)NM_DENSE_ELEMENTS(nmatrix), (OutputType*)data->output, 
    sizeof(OutputType)*NM_DENSE_COUNT(nmatrix));
}

/** \brief Executes the fourier transform by calling the fftw_execute function 
 * and copies the output to the output nmatrix object, which can be accessed from
 * Ruby.
 *
 * @param[in] self       Object on which the function is called.
 * @param[in] plan_data  An internal data structure of type 
 *   NMatrix::FFTW::Plan::Data that is created by Data_Wrap_Struct in 
 *   nm_fftw_create_plan and which encapsulates the FFTW plan in a Ruby object.
 * @param[in] nmatrix    An NMatrix object (pre-allocated) into which the computed
 *   data will be copied.
 * @param[in] type       A number representing the type of fourier transform being
 *   performed. (:complex_complex, :real_complex, :complex_real or :real_real).
 *
 * \returns TrueClass if computation completed without errors.
 */
static VALUE nm_fftw_execute(VALUE self, VALUE nmatrix, VALUE plan_data, VALUE type)
{
  switch(FIX2INT(type))
  {
    case TYPE_COMPLEX_COMPLEX:
    case TYPE_REAL_COMPLEX:
      nm_fftw_actually_execute<fftw_complex>(nmatrix, plan_data);
      break;
    case TYPE_COMPLEX_REAL:
    case TYPE_REAL_REAL:
      nm_fftw_actually_execute<double>(nmatrix, plan_data);
      break;
    default:
      rb_raise(rb_eTypeError, "Invalid type of DFT.");
  }

  return Qtrue;
}

extern "C" {
  void Init_nmatrix_fftw() 
  {
    VALUE cNMatrix                = rb_define_class("NMatrix", rb_cObject);
    VALUE cNMatrix_FFTW           = rb_define_module_under(cNMatrix, "FFTW");
    VALUE cNMatrix_FFTW_Plan      = rb_define_class_under(cNMatrix_FFTW, "Plan", 
      rb_cObject);
    VALUE cNMatrix_FFTW_Plan_Data = rb_define_class_under(
      cNMatrix_FFTW_Plan, "Data", rb_cObject);

    rb_define_private_method(cNMatrix_FFTW_Plan, "c_create_plan", 
      (METHOD)nm_fftw_create_plan, 7);
    rb_define_private_method(cNMatrix_FFTW_Plan, "c_set_input",
      (METHOD)nm_fftw_set_input, 3);
    rb_define_private_method(cNMatrix_FFTW_Plan, "c_execute",
      (METHOD)nm_fftw_execute, 3);
  }
}
