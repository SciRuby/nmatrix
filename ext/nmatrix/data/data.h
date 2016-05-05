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
// == data.h
//
// Header file for dealing with data types.

#ifndef DATA_H
#define DATA_H

/*
 * Standard Includes
 */

#include <ruby.h>
#include <string>

/*
 * Project Includes
 */

#include "nmatrix.h"

#include "types.h"

#include "complex.h"
#include "ruby_object.h"

namespace nm {


  /*
   * Constants
   */

  const int NUM_DTYPES = 10;
  const int NUM_ITYPES = 4;
  const int NUM_EWOPS = 12;
  const int NUM_UNARYOPS = 24;
  const int NUM_NONCOM_EWOPS = 3;

  enum ewop_t {
    EW_ADD,
    EW_SUB,
    EW_MUL,
    EW_DIV,
    EW_POW,
    EW_MOD,
    EW_EQEQ,
    EW_NEQ,
    EW_LT,
    EW_GT,
    EW_LEQ,
    EW_GEQ,
  };

  enum noncom_ewop_t {
    NONCOM_EW_ATAN2,
    NONCOM_EW_LDEXP,
    NONCOM_EW_HYPOT
  };

  enum unaryop_t {
    UNARY_SIN,
    UNARY_COS,
    UNARY_TAN,
    UNARY_ASIN,
    UNARY_ACOS,
    UNARY_ATAN,
    UNARY_SINH,
    UNARY_COSH,
    UNARY_TANH,
    UNARY_ASINH,
    UNARY_ACOSH,
    UNARY_ATANH,
    UNARY_EXP,
    UNARY_LOG2,
    UNARY_LOG10,
    UNARY_SQRT,
    UNARY_ERF,
    UNARY_ERFC,
    UNARY_CBRT,
    UNARY_GAMMA,
    UNARY_NEGATE,
    UNARY_FLOOR,
    UNARY_CEIL,
    UNARY_ROUND
  };

  // element-wise and scalar operators
  extern const char* const  EWOP_OPS[nm::NUM_EWOPS];
  extern const std::string  EWOP_NAMES[nm::NUM_EWOPS];
  extern const std::string  UNARYOPS[nm::NUM_UNARYOPS];
  extern const std::string  NONCOM_EWOP_NAMES[nm::NUM_NONCOM_EWOPS];


  template <typename Type>
  Complex<Type>::Complex(const RubyObject& other) {
    *this = other;
  }

  template <typename Type>
  Complex<Type>& Complex<Type>::operator=(const RubyObject& other) {
    if (RB_TYPE_P(other.rval, T_COMPLEX)) {
      this->r = NUM2DBL(rb_funcall(other.rval, rb_intern("real"), 0));
      this->i = NUM2DBL(rb_funcall(other.rval, rb_intern("imag"), 0));
    }
    else if (RB_TYPE_P(other.rval, T_FLOAT) ||
             RB_TYPE_P(other.rval, T_FIXNUM) ||
             RB_TYPE_P(other.rval, T_BIGNUM)) {
      this->r = NUM2DBL(other.rval);
      this->i = 0.0;
    }
    else {
      rb_raise(rb_eTypeError, "not sure how to convert this type of VALUE to a complex");
    }
    return *this;
  }

  template<typename Type>
  Complex<Type>::operator RubyObject () const {
    return RubyObject(*this);
  }

  nm::RubyObject  rubyobj_from_cval(void* val, nm::dtype_t dtype);
} // end of namespace nm

/*
 * Macros
 */

#define STYPE_MARK_TABLE(name)                  \
  static void (*(name)[nm::NUM_STYPES])(STORAGE*) = {  \
    nm_dense_storage_mark,                      \
    nm_list_storage_mark,                        \
    nm_yale_storage_mark                        \
  };

#define STYPE_REGISTER_TABLE(name)              \
  static void (*(name)[nm::NUM_STYPES])(const STORAGE*) = { \
    nm_dense_storage_register,                  \
    nm_list_storage_register,                   \
    nm_yale_storage_register                    \
  };

#define STYPE_UNREGISTER_TABLE(name)              \
  static void (*(name)[nm::NUM_STYPES])(const STORAGE*) = { \
    nm_dense_storage_unregister,                \
    nm_list_storage_unregister,                 \
    nm_yale_storage_unregister                  \
  };

#define CAST_TABLE(name)                                                   \
  static STORAGE* (*(name)[nm::NUM_STYPES][nm::NUM_STYPES])(const STORAGE*, nm::dtype_t, void*) = {      \
    { nm_dense_storage_cast_copy,  nm_dense_storage_from_list,  nm_dense_storage_from_yale },  \
    { nm_list_storage_from_dense,  nm_list_storage_cast_copy,   nm_list_storage_from_yale  },  \
    { nm_yale_storage_from_dense,  nm_yale_storage_from_list,   nm_yale_storage_cast_copy  }   \
  };

/*
 * Defines a static array that hold function pointers to dtype templated
 * versions of the specified function.
 */
#define DTYPE_TEMPLATE_TABLE(fun, ret, ...) NAMED_DTYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)

#define NAMED_DTYPE_TEMPLATE_TABLE(name, fun, ret, ...) \
  static ret (*(name)[nm::NUM_DTYPES])(__VA_ARGS__) =  { \
    fun<uint8_t>,                                        \
    fun<int8_t>,                                        \
    fun<int16_t>,                                        \
    fun<int32_t>,                                        \
    fun<int64_t>,                                        \
    fun<float32_t>,                                      \
    fun<float64_t>,                                      \
    fun<nm::Complex64>,                                  \
    fun<nm::Complex128>,                                \
    fun<nm::RubyObject>                                 \
  };

#define DTYPE_OBJECT_STATIC_TABLE(obj, fun, ret, ...)     \
  static ret (*(ttable)[nm::NUM_DTYPES])(__VA_ARGS__) =  { \
    obj<uint8_t>::fun,                                  \
    obj<int8_t>::fun,                                    \
    obj<int16_t>::fun,                                  \
    obj<int32_t>::fun,                                  \
    obj<int64_t>::fun,                                  \
    obj<float32_t>::fun,                                \
    obj<float64_t>::fun,                                \
    obj<nm::Complex64>::fun,                            \
    obj<nm::Complex128>::fun,                            \
    obj<nm::RubyObject>::fun                            \
  };

#define NAMED_DTYPE_TEMPLATE_TABLE_NO_ROBJ(name, fun, ret, ...) \
  static ret (*(name)[nm::NUM_DTYPES])(__VA_ARGS__) =  {      \
    fun<uint8_t>,                                        \
    fun<int8_t>,                                        \
    fun<int16_t>,                                        \
    fun<int32_t>,                                        \
    fun<int64_t>,                                        \
    fun<float32_t>,                                      \
    fun<float64_t>,                                      \
    fun<nm::Complex64>,                                  \
    fun<nm::Complex128>                                \
  };


/*
 * Same as DTYPE_TEMPLATE_TABLE but for functions that have two template
 * parameters.
 *
 * The left-hand DType is used as the first index, and the right-hand side is
 * the second index.  Not all left- and right-hand side combinations are valid,
 * and an invalid combination will result in a NULL pointer.
 */
#define LR_DTYPE_TEMPLATE_TABLE(fun, ret, ...) NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)

#define NAMED_LR_DTYPE_TEMPLATE_TABLE(name, fun, ret, ...)                                                                                                                \
  static ret (*(name)[nm::NUM_DTYPES][nm::NUM_DTYPES])(__VA_ARGS__) = {  \
    {fun<uint8_t, uint8_t>, fun<uint8_t, int8_t>, fun<uint8_t, int16_t>, fun<uint8_t, int32_t>, fun<uint8_t, int64_t>, fun<uint8_t, float32_t>, fun<uint8_t, float64_t>, fun<uint8_t, nm::Complex64>, fun<uint8_t, nm::Complex128>, fun<uint8_t, nm::RubyObject>}, \
    {fun<int8_t, uint8_t>, fun<int8_t, int8_t>, fun<int8_t, int16_t>, fun<int8_t, int32_t>, fun<int8_t, int64_t>, fun<int8_t, float32_t>, fun<int8_t, float64_t>, fun<int8_t, nm::Complex64>, fun<int8_t, nm::Complex128>, fun<int8_t, nm::RubyObject>},               \
    {fun<int16_t, uint8_t>, fun<int16_t, int8_t>, fun<int16_t, int16_t>, fun<int16_t, int32_t>, fun<int16_t, int64_t>, fun<int16_t, float32_t>, fun<int16_t, float64_t>, fun<int16_t, nm::Complex64>, fun<int16_t, nm::Complex128>, fun<int16_t, nm::RubyObject>},  \
    {fun<int32_t, uint8_t>, fun<int32_t, int8_t>, fun<int32_t, int16_t>, fun<int32_t, int32_t>, fun<int32_t, int64_t>, fun<int32_t, float32_t>, fun<int32_t, float64_t>, fun<int32_t, nm::Complex64>, fun<int32_t, nm::Complex128>, fun<int32_t, nm::RubyObject>},  \
    {fun<int64_t, uint8_t>, fun<int64_t, int8_t>, fun<int64_t, int16_t>, fun<int64_t, int32_t>, fun<int64_t, int64_t>, fun<int64_t, float32_t>, fun<int64_t, float64_t>, fun<int64_t, nm::Complex64>, fun<int64_t, nm::Complex128>, fun<int64_t, nm::RubyObject>},  \
    {fun<float32_t, uint8_t>, fun<float32_t, int8_t>, fun<float32_t, int16_t>, fun<float32_t, int32_t>, fun<float32_t, int64_t>, fun<float32_t, float32_t>, fun<float32_t, float64_t>, fun<float32_t, nm::Complex64>, fun<float32_t, nm::Complex128>, fun<float32_t, nm::RubyObject>},  \
    {fun<float64_t, uint8_t>, fun<float64_t, int8_t>, fun<float64_t, int16_t>, fun<float64_t, int32_t>, fun<float64_t, int64_t>, fun<float64_t, float32_t>, fun<float64_t, float64_t>, fun<float64_t, nm::Complex64>, fun<float64_t, nm::Complex128>, fun<float64_t, nm::RubyObject>},  \
    {fun<nm::Complex64, uint8_t>, fun<nm::Complex64, int8_t>, fun<nm::Complex64, int16_t>, fun<nm::Complex64, int32_t>, fun<nm::Complex64, int64_t>, fun<nm::Complex64, float32_t>, fun<nm::Complex64, float64_t>, fun<nm::Complex64, nm::Complex64>, fun<nm::Complex64, nm::Complex128>, fun<nm::Complex64, nm::RubyObject>},               \
    {fun<nm::Complex128, uint8_t>, fun<nm::Complex128, int8_t>, fun<nm::Complex128, int16_t>, fun<nm::Complex128, int32_t>, fun<nm::Complex128, int64_t>, fun<nm::Complex128, float32_t>, fun<nm::Complex128, float64_t>, fun<nm::Complex128, nm::Complex64>, fun<nm::Complex128, nm::Complex128>, fun<nm::Complex128, nm::RubyObject>},  \
    {fun<nm::RubyObject, uint8_t>, fun<nm::RubyObject, int8_t>, fun<nm::RubyObject, int16_t>, fun<nm::RubyObject, int32_t>, fun<nm::RubyObject, int64_t>, fun<nm::RubyObject, float32_t>, fun<nm::RubyObject, float64_t>, fun<nm::RubyObject, nm::Complex64>, fun<nm::RubyObject, nm::Complex128>, fun<nm::RubyObject, nm::RubyObject>}   \
  };

/*
 * Defines a static array that holds function pointers to operation, and left-
 * and right-side dtype templated version sof the specified function.
 */
#define OP_LR_DTYPE_TEMPLATE_TABLE(fun, ret, ...) NAMED_OP_LR_DTYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)

#define NAMED_OP_LR_DTYPE_TEMPLATE_TABLE(name, fun, ret, ...)                                                                                                               \
  static ret (*(name)[nm::NUM_EWOPS][nm::NUM_DTYPES][nm::NUM_DTYPES])(__VA_ARGS__) = {                                                                                                  \
    {                                                                                                                                                                        \
      {fun<nm::EW_ADD, uint8_t, uint8_t>, fun<nm::EW_ADD, uint8_t, int8_t>, fun<nm::EW_ADD, uint8_t, int16_t>, fun<nm::EW_ADD, uint8_t, int32_t>, fun<nm::EW_ADD, uint8_t, int64_t>,            \
        fun<nm::EW_ADD, uint8_t, float32_t>, fun<nm::EW_ADD, uint8_t, float64_t>, fun<nm::EW_ADD, uint8_t, nm::Complex64>, fun<nm::EW_ADD, uint8_t, nm::Complex128>,                        \
        fun<nm::EW_ADD, int8_t, float32_t>, fun<nm::EW_ADD, int8_t, float64_t>, fun<nm::EW_ADD, int8_t, nm::Complex64>, fun<nm::EW_ADD, int8_t, nm::Complex128>,                            \
         NULL},                                              \
                                                                                                                                                                            \
      {fun<nm::EW_ADD, int16_t, uint8_t>, fun<nm::EW_ADD, int16_t, int8_t>, fun<nm::EW_ADD, int16_t, int16_t>, fun<nm::EW_ADD, int16_t, int32_t>, fun<nm::EW_ADD, int16_t, int64_t>,            \
        fun<nm::EW_ADD, int16_t, float32_t>, fun<nm::EW_ADD, int16_t, float64_t>, fun<nm::EW_ADD, int16_t, nm::Complex64>, fun<nm::EW_ADD, int16_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_ADD, int32_t, uint8_t>, fun<nm::EW_ADD, int32_t, int8_t>, fun<nm::EW_ADD, int32_t, int16_t>, fun<nm::EW_ADD, int32_t, int32_t>, fun<nm::EW_ADD, int32_t, int64_t>,            \
        fun<nm::EW_ADD, int32_t, float32_t>, fun<nm::EW_ADD, int32_t, float64_t>, fun<nm::EW_ADD, int32_t, nm::Complex64>, fun<nm::EW_ADD, int32_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_ADD, int64_t, uint8_t>, fun<nm::EW_ADD, int64_t, int8_t>, fun<nm::EW_ADD, int64_t, int16_t>, fun<nm::EW_ADD, int64_t, int32_t>, fun<nm::EW_ADD, int64_t, int64_t>,            \
        fun<nm::EW_ADD, int64_t, float32_t>, fun<nm::EW_ADD, int64_t, float64_t>, fun<nm::EW_ADD, int64_t, nm::Complex64>, fun<nm::EW_ADD, int64_t, nm::Complex128>,                        \
         NULL},                                           \
                                                                                                                                                                            \
      {fun<nm::EW_ADD, float32_t, uint8_t>, fun<nm::EW_ADD, float32_t, int8_t>, fun<nm::EW_ADD, float32_t, int16_t>, fun<nm::EW_ADD, float32_t, int32_t>, fun<nm::EW_ADD, float32_t, int64_t>,  \
        fun<nm::EW_ADD, float32_t, float32_t>, fun<nm::EW_ADD, float32_t, float64_t>, fun<nm::EW_ADD, float32_t, nm::Complex64>, fun<nm::EW_ADD, float32_t, nm::Complex128>,                \
         NULL},                                      \
                                                                                                                                                                            \
      {fun<nm::EW_ADD, float64_t, uint8_t>, fun<nm::EW_ADD, float64_t, int8_t>, fun<nm::EW_ADD, float64_t, int16_t>, fun<nm::EW_ADD, float64_t, int32_t>, fun<nm::EW_ADD, float64_t, int64_t>,  \
        fun<nm::EW_ADD, float64_t, float32_t>, fun<nm::EW_ADD, float64_t, float64_t>, fun<nm::EW_ADD, float64_t, nm::Complex64>, fun<nm::EW_ADD, float64_t, nm::Complex128>,                \
         NULL},                                      \
                                                                                                                                                                            \
      {fun<nm::EW_ADD, nm::Complex64, uint8_t>, fun<nm::EW_ADD, nm::Complex64, int8_t>, fun<nm::EW_ADD, nm::Complex64, int16_t>, fun<nm::EW_ADD, nm::Complex64, int32_t>,                    \
        fun<nm::EW_ADD, nm::Complex64, int64_t>, fun<nm::EW_ADD, nm::Complex64, float32_t>, fun<nm::EW_ADD, nm::Complex64, float64_t>, fun<nm::EW_ADD, nm::Complex64, nm::Complex64>,        \
        fun<nm::EW_ADD, nm::Complex64, nm::Complex128>,                                  \
         NULL},                                                                                                                  \
                                                                                                                                                                            \
      {fun<nm::EW_ADD, nm::Complex128, uint8_t>, fun<nm::EW_ADD, nm::Complex128, int8_t>, fun<nm::EW_ADD, nm::Complex128, int16_t>, fun<nm::EW_ADD, nm::Complex128, int32_t>,                \
        fun<nm::EW_ADD, nm::Complex128, int64_t>, fun<nm::EW_ADD, nm::Complex128, float32_t>, fun<nm::EW_ADD, nm::Complex128, float64_t>, fun<nm::EW_ADD, nm::Complex128, nm::Complex64>,    \
        fun<nm::EW_ADD, nm::Complex128, nm::Complex128>,                                \
         NULL},                                                                                                                \
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_ADD, nm::RubyObject, nm::RubyObject>}                                                  \
    },                                                                                                                                                                      \
                                                                                                                                                                            \
    {                                                                                                                                                                        \
      {fun<nm::EW_SUB, uint8_t, uint8_t>, fun<nm::EW_SUB, uint8_t, int8_t>, fun<nm::EW_SUB, uint8_t, int16_t>, fun<nm::EW_SUB, uint8_t, int32_t>, fun<nm::EW_SUB, uint8_t, int64_t>,            \
        fun<nm::EW_SUB, uint8_t, float32_t>, fun<nm::EW_SUB, uint8_t, float64_t>, fun<nm::EW_SUB, uint8_t, nm::Complex64>, fun<nm::EW_SUB, uint8_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_SUB, int8_t, uint8_t>, fun<nm::EW_SUB, int8_t, int8_t>, fun<nm::EW_SUB, int8_t, int16_t>, fun<nm::EW_SUB, int8_t, int32_t>, fun<nm::EW_SUB, int8_t, int64_t>,                  \
        fun<nm::EW_SUB, int8_t, float32_t>, fun<nm::EW_SUB, int8_t, float64_t>, fun<nm::EW_SUB, int8_t, nm::Complex64>, fun<nm::EW_SUB, int8_t, nm::Complex128>,                            \
         NULL},                                              \
                                                                                                                                                                            \
      {fun<nm::EW_SUB, int16_t, uint8_t>, fun<nm::EW_SUB, int16_t, int8_t>, fun<nm::EW_SUB, int16_t, int16_t>, fun<nm::EW_SUB, int16_t, int32_t>, fun<nm::EW_SUB, int16_t, int64_t>,            \
        fun<nm::EW_SUB, int16_t, float32_t>, fun<nm::EW_SUB, int16_t, float64_t>, fun<nm::EW_SUB, int16_t, nm::Complex64>, fun<nm::EW_SUB, int16_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_SUB, int32_t, uint8_t>, fun<nm::EW_SUB, int32_t, int8_t>, fun<nm::EW_SUB, int32_t, int16_t>, fun<nm::EW_SUB, int32_t, int32_t>, fun<nm::EW_SUB, int32_t, int64_t>,            \
        fun<nm::EW_SUB, int32_t, float32_t>, fun<nm::EW_SUB, int32_t, float64_t>, fun<nm::EW_SUB, int32_t, nm::Complex64>, fun<nm::EW_SUB, int32_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_SUB, int64_t, uint8_t>, fun<nm::EW_SUB, int64_t, int8_t>, fun<nm::EW_SUB, int64_t, int16_t>, fun<nm::EW_SUB, int64_t, int32_t>, fun<nm::EW_SUB, int64_t, int64_t>,            \
        fun<nm::EW_SUB, int64_t, float32_t>, fun<nm::EW_SUB, int64_t, float64_t>, fun<nm::EW_SUB, int64_t, nm::Complex64>, fun<nm::EW_SUB, int64_t, nm::Complex128>,                        \
         NULL},                                           \
                                                                                                                                                                            \
      {fun<nm::EW_SUB, float32_t, uint8_t>, fun<nm::EW_SUB, float32_t, int8_t>, fun<nm::EW_SUB, float32_t, int16_t>, fun<nm::EW_SUB, float32_t, int32_t>, fun<nm::EW_SUB, float32_t, int64_t>,  \
        fun<nm::EW_SUB, float32_t, float32_t>, fun<nm::EW_SUB, float32_t, float64_t>, fun<nm::EW_SUB, float32_t, nm::Complex64>, fun<nm::EW_SUB, float32_t, nm::Complex128>,                \
         NULL},                                      \
                                                                                                                                                                            \
      {fun<nm::EW_SUB, float64_t, uint8_t>, fun<nm::EW_SUB, float64_t, int8_t>, fun<nm::EW_SUB, float64_t, int16_t>, fun<nm::EW_SUB, float64_t, int32_t>, fun<nm::EW_SUB, float64_t, int64_t>,  \
        fun<nm::EW_SUB, float64_t, float32_t>, fun<nm::EW_SUB, float64_t, float64_t>, fun<nm::EW_SUB, float64_t, nm::Complex64>, fun<nm::EW_SUB, float64_t, nm::Complex128>,                \
         NULL},                                      \
                                                                                                                                                                            \
      {fun<nm::EW_SUB, nm::Complex64, uint8_t>, fun<nm::EW_SUB, nm::Complex64, int8_t>, fun<nm::EW_SUB, nm::Complex64, int16_t>, fun<nm::EW_SUB, nm::Complex64, int32_t>,                    \
        fun<nm::EW_SUB, nm::Complex64, int64_t>, fun<nm::EW_SUB, nm::Complex64, float32_t>, fun<nm::EW_SUB, nm::Complex64, float64_t>, fun<nm::EW_SUB, nm::Complex64, nm::Complex64>,        \
        fun<nm::EW_SUB, nm::Complex64, nm::Complex128>,                                  \
         NULL},                                                                                                                  \
                                                                                                                                                                            \
      {fun<nm::EW_SUB, nm::Complex128, uint8_t>, fun<nm::EW_SUB, nm::Complex128, int8_t>, fun<nm::EW_SUB, nm::Complex128, int16_t>, fun<nm::EW_SUB, nm::Complex128, int32_t>,                \
        fun<nm::EW_SUB, nm::Complex128, int64_t>, fun<nm::EW_SUB, nm::Complex128, float32_t>, fun<nm::EW_SUB, nm::Complex128, float64_t>, fun<nm::EW_SUB, nm::Complex128, nm::Complex64>,    \
        fun<nm::EW_SUB, nm::Complex128, nm::Complex128>,                                \
                                                                                                                                                                            \
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_SUB, nm::RubyObject, nm::RubyObject>}                                                  \
    },                                                                                                                                                                      \
                                                                                                                                                                            \
    {                                                                                                                                                                        \
      {fun<nm::EW_MUL, uint8_t, uint8_t>, fun<nm::EW_MUL, uint8_t, int8_t>, fun<nm::EW_MUL, uint8_t, int16_t>, fun<nm::EW_MUL, uint8_t, int32_t>, fun<nm::EW_MUL, uint8_t, int64_t>,            \
        fun<nm::EW_MUL, uint8_t, float32_t>, fun<nm::EW_MUL, uint8_t, float64_t>, fun<nm::EW_MUL, uint8_t, nm::Complex64>, fun<nm::EW_MUL, uint8_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_MUL, int8_t, uint8_t>, fun<nm::EW_MUL, int8_t, int8_t>, fun<nm::EW_MUL, int8_t, int16_t>, fun<nm::EW_MUL, int8_t, int32_t>, fun<nm::EW_MUL, int8_t, int64_t>,                  \
        fun<nm::EW_MUL, int8_t, float32_t>, fun<nm::EW_MUL, int8_t, float64_t>, fun<nm::EW_MUL, int8_t, nm::Complex64>, fun<nm::EW_MUL, int8_t, nm::Complex128>,                            \
         NULL},                                              \
                                                                                                                                                                            \
      {fun<nm::EW_MUL, int16_t, uint8_t>, fun<nm::EW_MUL, int16_t, int8_t>, fun<nm::EW_MUL, int16_t, int16_t>, fun<nm::EW_MUL, int16_t, int32_t>, fun<nm::EW_MUL, int16_t, int64_t>,            \
        fun<nm::EW_MUL, int16_t, float32_t>, fun<nm::EW_MUL, int16_t, float64_t>, fun<nm::EW_MUL, int16_t, nm::Complex64>, fun<nm::EW_MUL, int16_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_MUL, int32_t, uint8_t>, fun<nm::EW_MUL, int32_t, int8_t>, fun<nm::EW_MUL, int32_t, int16_t>, fun<nm::EW_MUL, int32_t, int32_t>, fun<nm::EW_MUL, int32_t, int64_t>,            \
        fun<nm::EW_MUL, int32_t, float32_t>, fun<nm::EW_MUL, int32_t, float64_t>, fun<nm::EW_MUL, int32_t, nm::Complex64>, fun<nm::EW_MUL, int32_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_MUL, int64_t, uint8_t>, fun<nm::EW_MUL, int64_t, int8_t>, fun<nm::EW_MUL, int64_t, int16_t>, fun<nm::EW_MUL, int64_t, int32_t>, fun<nm::EW_MUL, int64_t, int64_t>,            \
        fun<nm::EW_MUL, int64_t, float32_t>, fun<nm::EW_MUL, int64_t, float64_t>, fun<nm::EW_MUL, int64_t, nm::Complex64>, fun<nm::EW_MUL, int64_t, nm::Complex128>,                        \
         NULL},                                           \
                                                                                                                                                                            \
      {fun<nm::EW_MUL, float32_t, uint8_t>, fun<nm::EW_MUL, float32_t, int8_t>, fun<nm::EW_MUL, float32_t, int16_t>, fun<nm::EW_MUL, float32_t, int32_t>, fun<nm::EW_MUL, float32_t, int64_t>,  \
        fun<nm::EW_MUL, float32_t, float32_t>, fun<nm::EW_MUL, float32_t, float64_t>, fun<nm::EW_MUL, float32_t, nm::Complex64>, fun<nm::EW_MUL, float32_t, nm::Complex128>,                \
         NULL},                                      \
                                                                                                                                                                            \
      {fun<nm::EW_MUL, float64_t, uint8_t>, fun<nm::EW_MUL, float64_t, int8_t>, fun<nm::EW_MUL, float64_t, int16_t>, fun<nm::EW_MUL, float64_t, int32_t>, fun<nm::EW_MUL, float64_t, int64_t>,  \
        fun<nm::EW_MUL, float64_t, float32_t>, fun<nm::EW_MUL, float64_t, float64_t>, fun<nm::EW_MUL, float64_t, nm::Complex64>, fun<nm::EW_MUL, float64_t, nm::Complex128>,                \
         NULL},                                      \
                                                                                                                                                                            \
      {fun<nm::EW_MUL, nm::Complex64, uint8_t>, fun<nm::EW_MUL, nm::Complex64, int8_t>, fun<nm::EW_MUL, nm::Complex64, int16_t>, fun<nm::EW_MUL, nm::Complex64, int32_t>,                    \
        fun<nm::EW_MUL, nm::Complex64, int64_t>, fun<nm::EW_MUL, nm::Complex64, float32_t>, fun<nm::EW_MUL, nm::Complex64, float64_t>, fun<nm::EW_MUL, nm::Complex64, nm::Complex64>,        \
        fun<nm::EW_MUL, nm::Complex64, nm::Complex128>,                                  \
         NULL},                                                                                                                  \
                                                                                                                                                                            \
      {fun<nm::EW_MUL, nm::Complex128, uint8_t>, fun<nm::EW_MUL, nm::Complex128, int8_t>, fun<nm::EW_MUL, nm::Complex128, int16_t>, fun<nm::EW_MUL, nm::Complex128, int32_t>,                \
        fun<nm::EW_MUL, nm::Complex128, int64_t>, fun<nm::EW_MUL, nm::Complex128, float32_t>, fun<nm::EW_MUL, nm::Complex128, float64_t>, fun<nm::EW_MUL, nm::Complex128, nm::Complex64>,    \
        fun<nm::EW_MUL, nm::Complex128, nm::Complex128>,                                \
                                                                                                                                                                            \
    {                                                                                                                                                                        \
      {fun<nm::EW_DIV, uint8_t, uint8_t>, fun<nm::EW_DIV, uint8_t, int8_t>, fun<nm::EW_DIV, uint8_t, int16_t>, fun<nm::EW_DIV, uint8_t, int32_t>, fun<nm::EW_DIV, uint8_t, int64_t>,            \
        fun<nm::EW_DIV, uint8_t, float32_t>, fun<nm::EW_DIV, uint8_t, float64_t>, fun<nm::EW_DIV, uint8_t, nm::Complex64>, fun<nm::EW_DIV, uint8_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_DIV, int8_t, uint8_t>, fun<nm::EW_DIV, int8_t, int8_t>, fun<nm::EW_DIV, int8_t, int16_t>, fun<nm::EW_DIV, int8_t, int32_t>, fun<nm::EW_DIV, int8_t, int64_t>,                  \
        fun<nm::EW_DIV, int8_t, float32_t>, fun<nm::EW_DIV, int8_t, float64_t>, fun<nm::EW_DIV, int8_t, nm::Complex64>, fun<nm::EW_DIV, int8_t, nm::Complex128>,                            \
         NULL},                                              \
                                                                                                                                                                            \
      {fun<nm::EW_DIV, int16_t, uint8_t>, fun<nm::EW_DIV, int16_t, int8_t>, fun<nm::EW_DIV, int16_t, int16_t>, fun<nm::EW_DIV, int16_t, int32_t>, fun<nm::EW_DIV, int16_t, int64_t>,            \
        fun<nm::EW_DIV, int16_t, float32_t>, fun<nm::EW_DIV, int16_t, float64_t>, fun<nm::EW_DIV, int16_t, nm::Complex64>, fun<nm::EW_DIV, int16_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_DIV, int32_t, uint8_t>, fun<nm::EW_DIV, int32_t, int8_t>, fun<nm::EW_DIV, int32_t, int16_t>, fun<nm::EW_DIV, int32_t, int32_t>, fun<nm::EW_DIV, int32_t, int64_t>,            \
        fun<nm::EW_DIV, int32_t, float32_t>, fun<nm::EW_DIV, int32_t, float64_t>, fun<nm::EW_DIV, int32_t, nm::Complex64>, fun<nm::EW_DIV, int32_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_DIV, int64_t, uint8_t>, fun<nm::EW_DIV, int64_t, int8_t>, fun<nm::EW_DIV, int64_t, int16_t>, fun<nm::EW_DIV, int64_t, int32_t>, fun<nm::EW_DIV, int64_t, int64_t>,            \
        fun<nm::EW_DIV, int64_t, float32_t>, fun<nm::EW_DIV, int64_t, float64_t>, fun<nm::EW_DIV, int64_t, nm::Complex64>, fun<nm::EW_DIV, int64_t, nm::Complex128>,                        \
         NULL},                                           \
                                                                                                                                                                            \
      {fun<nm::EW_DIV, float32_t, uint8_t>, fun<nm::EW_DIV, float32_t, int8_t>, fun<nm::EW_DIV, float32_t, int16_t>, fun<nm::EW_DIV, float32_t, int32_t>, fun<nm::EW_DIV, float32_t, int64_t>,  \
        fun<nm::EW_DIV, float32_t, float32_t>, fun<nm::EW_DIV, float32_t, float64_t>, fun<nm::EW_DIV, float32_t, nm::Complex64>, fun<nm::EW_DIV, float32_t, nm::Complex128>,                \
         NULL},                                      \
                                                                                                                                                                            \
      {fun<nm::EW_DIV, float64_t, uint8_t>, fun<nm::EW_DIV, float64_t, int8_t>, fun<nm::EW_DIV, float64_t, int16_t>, fun<nm::EW_DIV, float64_t, int32_t>, fun<nm::EW_DIV, float64_t, int64_t>,  \
        fun<nm::EW_DIV, float64_t, float32_t>, fun<nm::EW_DIV, float64_t, float64_t>, fun<nm::EW_DIV, float64_t, nm::Complex64>, fun<nm::EW_DIV, float64_t, nm::Complex128>,                \
         NULL},                                      \
                                                                                                                                                                            \
      {fun<nm::EW_DIV, nm::Complex64, uint8_t>, fun<nm::EW_DIV, nm::Complex64, int8_t>, fun<nm::EW_DIV, nm::Complex64, int16_t>, fun<nm::EW_DIV, nm::Complex64, int32_t>,                    \
        fun<nm::EW_DIV, nm::Complex64, int64_t>, fun<nm::EW_DIV, nm::Complex64, float32_t>, fun<nm::EW_DIV, nm::Complex64, float64_t>, fun<nm::EW_DIV, nm::Complex64, nm::Complex64>,        \
        fun<nm::EW_DIV, nm::Complex64, nm::Complex128>,                                  \
         NULL},                                                                                                                  \
                                                                                                                                                                            \
      {fun<nm::EW_DIV, nm::Complex128, uint8_t>, fun<nm::EW_DIV, nm::Complex128, int8_t>, fun<nm::EW_DIV, nm::Complex128, int16_t>, fun<nm::EW_DIV, nm::Complex128, int32_t>,                \
        fun<nm::EW_DIV, nm::Complex128, int64_t>, fun<nm::EW_DIV, nm::Complex128, float32_t>, fun<nm::EW_DIV, nm::Complex128, float64_t>, fun<nm::EW_DIV, nm::Complex128, nm::Complex64>,    \
        fun<nm::EW_DIV, nm::Complex128, nm::Complex128>,                                \
         NULL},                                                                                                                \
\
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_DIV, nm::RubyObject, nm::RubyObject>}                                                  \
    },                                                                                                                                                                      \
      \
    { \
      {fun<nm::EW_POW, uint8_t, uint8_t>, fun<nm::EW_POW, uint8_t, int8_t>, fun<nm::EW_POW, uint8_t, int16_t>, fun<nm::EW_POW, uint8_t, int32_t>, fun<nm::EW_POW, uint8_t, int64_t>,            \
        fun<nm::EW_POW, uint8_t, float32_t>, fun<nm::EW_POW, uint8_t, float64_t>, fun<nm::EW_POW, uint8_t, nm::Complex64>, fun<nm::EW_POW, uint8_t, nm::Complex128>,                        \
 NULL},                                            \
\
      {fun<nm::EW_POW, int8_t, uint8_t>, fun<nm::EW_POW, int8_t, int8_t>, fun<nm::EW_POW, int8_t, int16_t>, fun<nm::EW_POW, int8_t, int32_t>, fun<nm::EW_POW, int8_t, int64_t>,                  \
        fun<nm::EW_POW, int8_t, float32_t>, fun<nm::EW_POW, int8_t, float64_t>, fun<nm::EW_POW, int8_t, nm::Complex64>, fun<nm::EW_POW, int8_t, nm::Complex128>,                            \
 NULL},                                              \
\
      {fun<nm::EW_POW, int16_t, uint8_t>, fun<nm::EW_POW, int16_t, int8_t>, fun<nm::EW_POW, int16_t, int16_t>, fun<nm::EW_POW, int16_t, int32_t>, fun<nm::EW_POW, int16_t, int64_t>,            \
        fun<nm::EW_POW, int16_t, float32_t>, fun<nm::EW_POW, int16_t, float64_t>, fun<nm::EW_POW, int16_t, nm::Complex64>, fun<nm::EW_POW, int16_t, nm::Complex128>,                        \
 NULL},                                            \
\
      {fun<nm::EW_POW, int32_t, uint8_t>, fun<nm::EW_POW, int32_t, int8_t>, fun<nm::EW_POW, int32_t, int16_t>, fun<nm::EW_POW, int32_t, int32_t>, fun<nm::EW_POW, int32_t, int64_t>,            \
        fun<nm::EW_POW, int32_t, float32_t>, fun<nm::EW_POW, int32_t, float64_t>, fun<nm::EW_POW, int32_t, nm::Complex64>, fun<nm::EW_POW, int32_t, nm::Complex128>,                        \
 NULL},                                            \
\
      {fun<nm::EW_POW, int64_t, uint8_t>, fun<nm::EW_POW, int64_t, int8_t>, fun<nm::EW_POW, int64_t, int16_t>, fun<nm::EW_POW, int64_t, int32_t>, fun<nm::EW_POW, int64_t, int64_t>,            \
        fun<nm::EW_POW, int64_t, float32_t>, fun<nm::EW_POW, int64_t, float64_t>, fun<nm::EW_POW, int64_t, nm::Complex64>, fun<nm::EW_POW, int64_t, nm::Complex128>,                        \
 NULL},                                           \
\
      {fun<nm::EW_POW, float32_t, uint8_t>, fun<nm::EW_POW, float32_t, int8_t>, fun<nm::EW_POW, float32_t, int16_t>, fun<nm::EW_POW, float32_t, int32_t>, fun<nm::EW_POW, float32_t, int64_t>,  \
        fun<nm::EW_POW, float32_t, float32_t>, fun<nm::EW_POW, float32_t, float64_t>, fun<nm::EW_POW, float32_t, nm::Complex64>, fun<nm::EW_POW, float32_t, nm::Complex128>,                \
 NULL},                                      \
\
      {fun<nm::EW_POW, float64_t, uint8_t>, fun<nm::EW_POW, float64_t, int8_t>, fun<nm::EW_POW, float64_t, int16_t>, fun<nm::EW_POW, float64_t, int32_t>, fun<nm::EW_POW, float64_t, int64_t>,  \
        fun<nm::EW_POW, float64_t, float32_t>, fun<nm::EW_POW, float64_t, float64_t>, fun<nm::EW_POW, float64_t, nm::Complex64>, fun<nm::EW_POW, float64_t, nm::Complex128>,                \
 NULL},                                      \
\
      {fun<nm::EW_POW, nm::Complex64, uint8_t>, fun<nm::EW_POW, nm::Complex64, int8_t>, fun<nm::EW_POW, nm::Complex64, int16_t>, fun<nm::EW_POW, nm::Complex64, int32_t>,                    \
        fun<nm::EW_POW, nm::Complex64, int64_t>, fun<nm::EW_POW, nm::Complex64, float32_t>, fun<nm::EW_POW, nm::Complex64, float64_t>, fun<nm::EW_POW, nm::Complex64, nm::Complex64>,        \
        fun<nm::EW_POW, nm::Complex64, nm::Complex128>,                                  \
 NULL},                                                                                                                  \
\
      {fun<nm::EW_POW, nm::Complex128, uint8_t>, fun<nm::EW_POW, nm::Complex128, int8_t>, fun<nm::EW_POW, nm::Complex128, int16_t>, fun<nm::EW_POW, nm::Complex128, int32_t>,                \
        fun<nm::EW_POW, nm::Complex128, int64_t>, fun<nm::EW_POW, nm::Complex128, float32_t>, fun<nm::EW_POW, nm::Complex128, float64_t>, fun<nm::EW_POW, nm::Complex128, nm::Complex64>,    \
        fun<nm::EW_POW, nm::Complex128, nm::Complex128>,                                \
 NULL},                                                                                                                \
\
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_POW, nm::RubyObject, nm::RubyObject>}                                                  \
    },  \
\
    {                                                                                                                                                                        \
      {fun<nm::EW_MOD, uint8_t, uint8_t>, fun<nm::EW_MOD, uint8_t, int8_t>, fun<nm::EW_MOD, uint8_t, int16_t>, fun<nm::EW_MOD, uint8_t, int32_t>, fun<nm::EW_MOD, uint8_t, int64_t>,            \
        fun<nm::EW_MOD, uint8_t, float32_t>, fun<nm::EW_MOD, uint8_t, float64_t>, fun<nm::EW_MOD, uint8_t, nm::Complex64>, fun<nm::EW_MOD, uint8_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_MOD, int8_t, uint8_t>, fun<nm::EW_MOD, int8_t, int8_t>, fun<nm::EW_MOD, int8_t, int16_t>, fun<nm::EW_MOD, int8_t, int32_t>, fun<nm::EW_MOD, int8_t, int64_t>,                  \
        fun<nm::EW_MOD, int8_t, float32_t>, fun<nm::EW_MOD, int8_t, float64_t>, fun<nm::EW_MOD, int8_t, nm::Complex64>, fun<nm::EW_MOD, int8_t, nm::Complex128>,                            \
         NULL},                                              \
                                                                                                                                                                            \
      {fun<nm::EW_MOD, int16_t, uint8_t>, fun<nm::EW_MOD, int16_t, int8_t>, fun<nm::EW_MOD, int16_t, int16_t>, fun<nm::EW_MOD, int16_t, int32_t>, fun<nm::EW_MOD, int16_t, int64_t>,            \
        fun<nm::EW_MOD, int16_t, float32_t>, fun<nm::EW_MOD, int16_t, float64_t>, fun<nm::EW_MOD, int16_t, nm::Complex64>, fun<nm::EW_MOD, int16_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_MOD, int32_t, uint8_t>, fun<nm::EW_MOD, int32_t, int8_t>, fun<nm::EW_MOD, int32_t, int16_t>, fun<nm::EW_MOD, int32_t, int32_t>, fun<nm::EW_MOD, int32_t, int64_t>,            \
        fun<nm::EW_MOD, int32_t, float32_t>, fun<nm::EW_MOD, int32_t, float64_t>, fun<nm::EW_MOD, int32_t, nm::Complex64>, fun<nm::EW_MOD, int32_t, nm::Complex128>,                        \
         NULL},                                            \
                                                                                                                                                                            \
      {fun<nm::EW_MOD, int64_t, uint8_t>, fun<nm::EW_MOD, int64_t, int8_t>, fun<nm::EW_MOD, int64_t, int16_t>, fun<nm::EW_MOD, int64_t, int32_t>, fun<nm::EW_MOD, int64_t, int64_t>,            \
        fun<nm::EW_MOD, int64_t, float32_t>, fun<nm::EW_MOD, int64_t, float64_t>, fun<nm::EW_MOD, int64_t, nm::Complex64>, fun<nm::EW_MOD, int64_t, nm::Complex128>,                        \
         NULL},                                           \
                                                                                                                                                                            \
      {fun<nm::EW_MOD, float32_t, uint8_t>, fun<nm::EW_MOD, float32_t, int8_t>, fun<nm::EW_MOD, float32_t, int16_t>, fun<nm::EW_MOD, float32_t, int32_t>, fun<nm::EW_MOD, float32_t, int64_t>,  \
        fun<nm::EW_MOD, float32_t, float32_t>, fun<nm::EW_MOD, float32_t, float64_t>, fun<nm::EW_MOD, float32_t, nm::Complex64>, fun<nm::EW_MOD, float32_t, nm::Complex128>,                \
         NULL},                                      \
                                                                                                                                                                            \
      {fun<nm::EW_MOD, float64_t, uint8_t>, fun<nm::EW_MOD, float64_t, int8_t>, fun<nm::EW_MOD, float64_t, int16_t>, fun<nm::EW_MOD, float64_t, int32_t>, fun<nm::EW_MOD, float64_t, int64_t>,  \
        fun<nm::EW_MOD, float64_t, float32_t>, fun<nm::EW_MOD, float64_t, float64_t>, fun<nm::EW_MOD, float64_t, nm::Complex64>, fun<nm::EW_MOD, float64_t, nm::Complex128>,                \
         NULL},                                      \
                                                                                                                                                                            \
      {fun<nm::EW_MOD, nm::Complex64, uint8_t>, fun<nm::EW_MOD, nm::Complex64, int8_t>, fun<nm::EW_MOD, nm::Complex64, int16_t>, fun<nm::EW_MOD, nm::Complex64, int32_t>,                    \
        fun<nm::EW_MOD, nm::Complex64, int64_t>, fun<nm::EW_MOD, nm::Complex64, float32_t>, fun<nm::EW_MOD, nm::Complex64, float64_t>, fun<nm::EW_MOD, nm::Complex64, nm::Complex64>,        \
        fun<nm::EW_MOD, nm::Complex64, nm::Complex128>,                                  \
         NULL},                                                                                                                  \
                                                                                                                                                                            \
      {fun<nm::EW_MOD, nm::Complex128, uint8_t>, fun<nm::EW_MOD, nm::Complex128, int8_t>, fun<nm::EW_MOD, nm::Complex128, int16_t>, fun<nm::EW_MOD, nm::Complex128, int32_t>,                \
        fun<nm::EW_MOD, nm::Complex128, int64_t>, fun<nm::EW_MOD, nm::Complex128, float32_t>, fun<nm::EW_MOD, nm::Complex128, float64_t>, fun<nm::EW_MOD, nm::Complex128, nm::Complex64>,    \
        fun<nm::EW_MOD, nm::Complex128, nm::Complex128>,                                \
         NULL},                                                                                                                \
\
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_MOD, nm::RubyObject, nm::RubyObject>}                                                  \
    },                                                                                                                                                                      \
                                                                                                                                                                            \
    {                                                                                                                                                                       \
      {fun<nm::EW_EQEQ, uint8_t, uint8_t>, fun<nm::EW_EQEQ, uint8_t, int8_t>, fun<nm::EW_EQEQ, uint8_t, int16_t>, fun<nm::EW_EQEQ, uint8_t, int32_t>, \
        fun<nm::EW_EQEQ, uint8_t, int64_t>, fun<nm::EW_EQEQ, uint8_t, float32_t>, fun<nm::EW_EQEQ, uint8_t, float64_t>, fun<nm::EW_EQEQ, uint8_t, nm::Complex64>, \
        fun<nm::EW_EQEQ, uint8_t, nm::Complex128>, \
 NULL}, \
      {fun<nm::EW_EQEQ, int8_t, uint8_t>, fun<nm::EW_EQEQ, int8_t, int8_t>, fun<nm::EW_EQEQ, int8_t, int16_t>, fun<nm::EW_EQEQ, int8_t, int32_t>, fun<nm::EW_EQEQ, int8_t, int64_t>, fun<nm::EW_EQEQ, int8_t, float32_t>, fun<nm::EW_EQEQ, int8_t, float64_t>, fun<nm::EW_EQEQ, int8_t, nm::Complex64>, fun<nm::EW_EQEQ, int8_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_EQEQ, int16_t, uint8_t>, fun<nm::EW_EQEQ, int16_t, int8_t>, fun<nm::EW_EQEQ, int16_t, int16_t>, fun<nm::EW_EQEQ, int16_t, int32_t>, fun<nm::EW_EQEQ, int16_t, int64_t>, fun<nm::EW_EQEQ, int16_t, float32_t>, fun<nm::EW_EQEQ, int16_t, float64_t>, fun<nm::EW_EQEQ, int16_t, nm::Complex64>, fun<nm::EW_EQEQ, int16_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_EQEQ, int32_t, uint8_t>, fun<nm::EW_EQEQ, int32_t, int8_t>, fun<nm::EW_EQEQ, int32_t, int16_t>, fun<nm::EW_EQEQ, int32_t, int32_t>, fun<nm::EW_EQEQ, int32_t, int64_t>, fun<nm::EW_EQEQ, int32_t, float32_t>, fun<nm::EW_EQEQ, int32_t, float64_t>, fun<nm::EW_EQEQ, int32_t, nm::Complex64>, fun<nm::EW_EQEQ, int32_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_EQEQ, int64_t, uint8_t>, fun<nm::EW_EQEQ, int64_t, int8_t>, fun<nm::EW_EQEQ, int64_t, int16_t>, fun<nm::EW_EQEQ, int64_t, int32_t>, fun<nm::EW_EQEQ, int64_t, int64_t>, fun<nm::EW_EQEQ, int64_t, float32_t>, fun<nm::EW_EQEQ, int64_t, float64_t>, fun<nm::EW_EQEQ, int64_t, nm::Complex64>, fun<nm::EW_EQEQ, int64_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_EQEQ, float32_t, uint8_t>, fun<nm::EW_EQEQ, float32_t, int8_t>, fun<nm::EW_EQEQ, float32_t, int16_t>, fun<nm::EW_EQEQ, float32_t, int32_t>, fun<nm::EW_EQEQ, float32_t, int64_t>, fun<nm::EW_EQEQ, float32_t, float32_t>, fun<nm::EW_EQEQ, float32_t, float64_t>, fun<nm::EW_EQEQ, float32_t, nm::Complex64>, fun<nm::EW_EQEQ, float32_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_EQEQ, float64_t, uint8_t>, fun<nm::EW_EQEQ, float64_t, int8_t>, fun<nm::EW_EQEQ, float64_t, int16_t>, fun<nm::EW_EQEQ, float64_t, int32_t>, fun<nm::EW_EQEQ, float64_t, int64_t>, fun<nm::EW_EQEQ, float64_t, float32_t>, fun<nm::EW_EQEQ, float64_t, float64_t>, fun<nm::EW_EQEQ, float64_t, nm::Complex64>, fun<nm::EW_EQEQ, float64_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_EQEQ, nm::Complex64, uint8_t>, fun<nm::EW_EQEQ, nm::Complex64, int8_t>, fun<nm::EW_EQEQ, nm::Complex64, int16_t>, fun<nm::EW_EQEQ, nm::Complex64, int32_t>, fun<nm::EW_EQEQ, nm::Complex64, int64_t>, fun<nm::EW_EQEQ, nm::Complex64, float32_t>, fun<nm::EW_EQEQ, nm::Complex64, float64_t>, fun<nm::EW_EQEQ, nm::Complex64, nm::Complex64>, fun<nm::EW_EQEQ, nm::Complex64, nm::Complex128>, NULL}, \
      {fun<nm::EW_EQEQ, nm::Complex128, uint8_t>, fun<nm::EW_EQEQ, nm::Complex128, int8_t>, fun<nm::EW_EQEQ, nm::Complex128, int16_t>, fun<nm::EW_EQEQ, nm::Complex128, int32_t>, fun<nm::EW_EQEQ, nm::Complex128, int64_t>, fun<nm::EW_EQEQ, nm::Complex128, float32_t>, fun<nm::EW_EQEQ, nm::Complex128, float64_t>, fun<nm::EW_EQEQ, nm::Complex128, nm::Complex64>, fun<nm::EW_EQEQ, nm::Complex128, nm::Complex128>, NULL}, \
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_EQEQ, nm::RubyObject, nm::RubyObject>}  \
    }, \
    {{fun<nm::EW_NEQ, uint8_t, uint8_t>, fun<nm::EW_NEQ, uint8_t, int8_t>, fun<nm::EW_NEQ, uint8_t, int16_t>, fun<nm::EW_NEQ, uint8_t, int32_t>, fun<nm::EW_NEQ, uint8_t, int64_t>, fun<nm::EW_NEQ, uint8_t, float32_t>, fun<nm::EW_NEQ, uint8_t, float64_t>, fun<nm::EW_NEQ, uint8_t, nm::Complex64>, fun<nm::EW_NEQ, uint8_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_NEQ, int8_t, uint8_t>, fun<nm::EW_NEQ, int8_t, int8_t>, fun<nm::EW_NEQ, int8_t, int16_t>, fun<nm::EW_NEQ, int8_t, int32_t>, fun<nm::EW_NEQ, int8_t, int64_t>, fun<nm::EW_NEQ, int8_t, float32_t>, fun<nm::EW_NEQ, int8_t, float64_t>, fun<nm::EW_NEQ, int8_t, nm::Complex64>, fun<nm::EW_NEQ, int8_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_NEQ, int16_t, uint8_t>, fun<nm::EW_NEQ, int16_t, int8_t>, fun<nm::EW_NEQ, int16_t, int16_t>, fun<nm::EW_NEQ, int16_t, int32_t>, fun<nm::EW_NEQ, int16_t, int64_t>, fun<nm::EW_NEQ, int16_t, float32_t>, fun<nm::EW_NEQ, int16_t, float64_t>, fun<nm::EW_NEQ, int16_t, nm::Complex64>, fun<nm::EW_NEQ, int16_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_NEQ, int32_t, uint8_t>, fun<nm::EW_NEQ, int32_t, int8_t>, fun<nm::EW_NEQ, int32_t, int16_t>, fun<nm::EW_NEQ, int32_t, int32_t>, fun<nm::EW_NEQ, int32_t, int64_t>, fun<nm::EW_NEQ, int32_t, float32_t>, fun<nm::EW_NEQ, int32_t, float64_t>, fun<nm::EW_NEQ, int32_t, nm::Complex64>, fun<nm::EW_NEQ, int32_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_NEQ, int64_t, uint8_t>, fun<nm::EW_NEQ, int64_t, int8_t>, fun<nm::EW_NEQ, int64_t, int16_t>, fun<nm::EW_NEQ, int64_t, int32_t>, fun<nm::EW_NEQ, int64_t, int64_t>, fun<nm::EW_NEQ, int64_t, float32_t>, fun<nm::EW_NEQ, int64_t, float64_t>, fun<nm::EW_NEQ, int64_t, nm::Complex64>, fun<nm::EW_NEQ, int64_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_NEQ, float32_t, uint8_t>, fun<nm::EW_NEQ, float32_t, int8_t>, fun<nm::EW_NEQ, float32_t, int16_t>, fun<nm::EW_NEQ, float32_t, int32_t>, fun<nm::EW_NEQ, float32_t, int64_t>, fun<nm::EW_NEQ, float32_t, float32_t>, fun<nm::EW_NEQ, float32_t, float64_t>, fun<nm::EW_NEQ, float32_t, nm::Complex64>, fun<nm::EW_NEQ, float32_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_NEQ, float64_t, uint8_t>, fun<nm::EW_NEQ, float64_t, int8_t>, fun<nm::EW_NEQ, float64_t, int16_t>, fun<nm::EW_NEQ, float64_t, int32_t>, fun<nm::EW_NEQ, float64_t, int64_t>, fun<nm::EW_NEQ, float64_t, float32_t>, fun<nm::EW_NEQ, float64_t, float64_t>, fun<nm::EW_NEQ, float64_t, nm::Complex64>, fun<nm::EW_NEQ, float64_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_NEQ, nm::Complex64, uint8_t>, fun<nm::EW_NEQ, nm::Complex64, int8_t>, fun<nm::EW_NEQ, nm::Complex64, int16_t>, fun<nm::EW_NEQ, nm::Complex64, int32_t>, fun<nm::EW_NEQ, nm::Complex64, int64_t>, fun<nm::EW_NEQ, nm::Complex64, float32_t>, fun<nm::EW_NEQ, nm::Complex64, float64_t>, fun<nm::EW_NEQ, nm::Complex64, nm::Complex64>, fun<nm::EW_NEQ, nm::Complex64, nm::Complex128>, NULL}, \
      {fun<nm::EW_NEQ, nm::Complex128, uint8_t>, fun<nm::EW_NEQ, nm::Complex128, int8_t>, fun<nm::EW_NEQ, nm::Complex128, int16_t>, fun<nm::EW_NEQ, nm::Complex128, int32_t>, fun<nm::EW_NEQ, nm::Complex128, int64_t>, fun<nm::EW_NEQ, nm::Complex128, float32_t>, fun<nm::EW_NEQ, nm::Complex128, float64_t>, fun<nm::EW_NEQ, nm::Complex128, nm::Complex64>, fun<nm::EW_NEQ, nm::Complex128, nm::Complex128>, NULL}, \
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_NEQ, nm::RubyObject, nm::RubyObject>}}, \
    {{fun<nm::EW_LT, uint8_t, uint8_t>, fun<nm::EW_LT, uint8_t, int8_t>, fun<nm::EW_LT, uint8_t, int16_t>, fun<nm::EW_LT, uint8_t, int32_t>, fun<nm::EW_LT, uint8_t, int64_t>, fun<nm::EW_LT, uint8_t, float32_t>, fun<nm::EW_LT, uint8_t, float64_t>, fun<nm::EW_LT, uint8_t, nm::Complex64>, fun<nm::EW_LT, uint8_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LT, int8_t, uint8_t>, fun<nm::EW_LT, int8_t, int8_t>, fun<nm::EW_LT, int8_t, int16_t>, fun<nm::EW_LT, int8_t, int32_t>, fun<nm::EW_LT, int8_t, int64_t>, fun<nm::EW_LT, int8_t, float32_t>, fun<nm::EW_LT, int8_t, float64_t>, fun<nm::EW_LT, int8_t, nm::Complex64>, fun<nm::EW_LT, int8_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LT, int16_t, uint8_t>, fun<nm::EW_LT, int16_t, int8_t>, fun<nm::EW_LT, int16_t, int16_t>, fun<nm::EW_LT, int16_t, int32_t>, fun<nm::EW_LT, int16_t, int64_t>, fun<nm::EW_LT, int16_t, float32_t>, fun<nm::EW_LT, int16_t, float64_t>, fun<nm::EW_LT, int16_t, nm::Complex64>, fun<nm::EW_LT, int16_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LT, int32_t, uint8_t>, fun<nm::EW_LT, int32_t, int8_t>, fun<nm::EW_LT, int32_t, int16_t>, fun<nm::EW_LT, int32_t, int32_t>, fun<nm::EW_LT, int32_t, int64_t>, fun<nm::EW_LT, int32_t, float32_t>, fun<nm::EW_LT, int32_t, float64_t>, fun<nm::EW_LT, int32_t, nm::Complex64>, fun<nm::EW_LT, int32_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LT, int64_t, uint8_t>, fun<nm::EW_LT, int64_t, int8_t>, fun<nm::EW_LT, int64_t, int16_t>, fun<nm::EW_LT, int64_t, int32_t>, fun<nm::EW_LT, int64_t, int64_t>, fun<nm::EW_LT, int64_t, float32_t>, fun<nm::EW_LT, int64_t, float64_t>, fun<nm::EW_LT, int64_t, nm::Complex64>, fun<nm::EW_LT, int64_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LT, float32_t, uint8_t>, fun<nm::EW_LT, float32_t, int8_t>, fun<nm::EW_LT, float32_t, int16_t>, fun<nm::EW_LT, float32_t, int32_t>, fun<nm::EW_LT, float32_t, int64_t>, fun<nm::EW_LT, float32_t, float32_t>, fun<nm::EW_LT, float32_t, float64_t>, fun<nm::EW_LT, float32_t, nm::Complex64>, fun<nm::EW_LT, float32_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LT, float64_t, uint8_t>, fun<nm::EW_LT, float64_t, int8_t>, fun<nm::EW_LT, float64_t, int16_t>, fun<nm::EW_LT, float64_t, int32_t>, fun<nm::EW_LT, float64_t, int64_t>, fun<nm::EW_LT, float64_t, float32_t>, fun<nm::EW_LT, float64_t, float64_t>, fun<nm::EW_LT, float64_t, nm::Complex64>, fun<nm::EW_LT, float64_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LT, nm::Complex64, uint8_t>, fun<nm::EW_LT, nm::Complex64, int8_t>, fun<nm::EW_LT, nm::Complex64, int16_t>, fun<nm::EW_LT, nm::Complex64, int32_t>, fun<nm::EW_LT, nm::Complex64, int64_t>, fun<nm::EW_LT, nm::Complex64, float32_t>, fun<nm::EW_LT, nm::Complex64, float64_t>, fun<nm::EW_LT, nm::Complex64, nm::Complex64>, fun<nm::EW_LT, nm::Complex64, nm::Complex128>, NULL}, \
      {fun<nm::EW_LT, nm::Complex128, uint8_t>, fun<nm::EW_LT, nm::Complex128, int8_t>, fun<nm::EW_LT, nm::Complex128, int16_t>, fun<nm::EW_LT, nm::Complex128, int32_t>, fun<nm::EW_LT, nm::Complex128, int64_t>, fun<nm::EW_LT, nm::Complex128, float32_t>, fun<nm::EW_LT, nm::Complex128, float64_t>, fun<nm::EW_LT, nm::Complex128, nm::Complex64>, fun<nm::EW_LT, nm::Complex128, nm::Complex128>, NULL}, \
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_LT, nm::RubyObject, nm::RubyObject>}}, \
    {{fun<nm::EW_GT, uint8_t, uint8_t>, fun<nm::EW_GT, uint8_t, int8_t>, fun<nm::EW_GT, uint8_t, int16_t>, fun<nm::EW_GT, uint8_t, int32_t>, fun<nm::EW_GT, uint8_t, int64_t>, fun<nm::EW_GT, uint8_t, float32_t>, fun<nm::EW_GT, uint8_t, float64_t>, fun<nm::EW_GT, uint8_t, nm::Complex64>, fun<nm::EW_GT, uint8_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GT, int8_t, uint8_t>, fun<nm::EW_GT, int8_t, int8_t>, fun<nm::EW_GT, int8_t, int16_t>, fun<nm::EW_GT, int8_t, int32_t>, fun<nm::EW_GT, int8_t, int64_t>, fun<nm::EW_GT, int8_t, float32_t>, fun<nm::EW_GT, int8_t, float64_t>, fun<nm::EW_GT, int8_t, nm::Complex64>, fun<nm::EW_GT, int8_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GT, int16_t, uint8_t>, fun<nm::EW_GT, int16_t, int8_t>, fun<nm::EW_GT, int16_t, int16_t>, fun<nm::EW_GT, int16_t, int32_t>, fun<nm::EW_GT, int16_t, int64_t>, fun<nm::EW_GT, int16_t, float32_t>, fun<nm::EW_GT, int16_t, float64_t>, fun<nm::EW_GT, int16_t, nm::Complex64>, fun<nm::EW_GT, int16_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GT, int32_t, uint8_t>, fun<nm::EW_GT, int32_t, int8_t>, fun<nm::EW_GT, int32_t, int16_t>, fun<nm::EW_GT, int32_t, int32_t>, fun<nm::EW_GT, int32_t, int64_t>, fun<nm::EW_GT, int32_t, float32_t>, fun<nm::EW_GT, int32_t, float64_t>, fun<nm::EW_GT, int32_t, nm::Complex64>, fun<nm::EW_GT, int32_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GT, int64_t, uint8_t>, fun<nm::EW_GT, int64_t, int8_t>, fun<nm::EW_GT, int64_t, int16_t>, fun<nm::EW_GT, int64_t, int32_t>, fun<nm::EW_GT, int64_t, int64_t>, fun<nm::EW_GT, int64_t, float32_t>, fun<nm::EW_GT, int64_t, float64_t>, fun<nm::EW_GT, int64_t, nm::Complex64>, fun<nm::EW_GT, int64_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GT, float32_t, uint8_t>, fun<nm::EW_GT, float32_t, int8_t>, fun<nm::EW_GT, float32_t, int16_t>, fun<nm::EW_GT, float32_t, int32_t>, fun<nm::EW_GT, float32_t, int64_t>, fun<nm::EW_GT, float32_t, float32_t>, fun<nm::EW_GT, float32_t, float64_t>, fun<nm::EW_GT, float32_t, nm::Complex64>, fun<nm::EW_GT, float32_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GT, float64_t, uint8_t>, fun<nm::EW_GT, float64_t, int8_t>, fun<nm::EW_GT, float64_t, int16_t>, fun<nm::EW_GT, float64_t, int32_t>, fun<nm::EW_GT, float64_t, int64_t>, fun<nm::EW_GT, float64_t, float32_t>, fun<nm::EW_GT, float64_t, float64_t>, fun<nm::EW_GT, float64_t, nm::Complex64>, fun<nm::EW_GT, float64_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GT, nm::Complex64, uint8_t>, fun<nm::EW_GT, nm::Complex64, int8_t>, fun<nm::EW_GT, nm::Complex64, int16_t>, fun<nm::EW_GT, nm::Complex64, int32_t>, fun<nm::EW_GT, nm::Complex64, int64_t>, fun<nm::EW_GT, nm::Complex64, float32_t>, fun<nm::EW_GT, nm::Complex64, float64_t>, fun<nm::EW_GT, nm::Complex64, nm::Complex64>, fun<nm::EW_GT, nm::Complex64, nm::Complex128>, NULL}, \
      {fun<nm::EW_GT, nm::Complex128, uint8_t>, fun<nm::EW_GT, nm::Complex128, int8_t>, fun<nm::EW_GT, nm::Complex128, int16_t>, fun<nm::EW_GT, nm::Complex128, int32_t>, fun<nm::EW_GT, nm::Complex128, int64_t>, fun<nm::EW_GT, nm::Complex128, float32_t>, fun<nm::EW_GT, nm::Complex128, float64_t>, fun<nm::EW_GT, nm::Complex128, nm::Complex64>, fun<nm::EW_GT, nm::Complex128, nm::Complex128>, NULL}, \
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_GT, nm::RubyObject, nm::RubyObject>}}, \
    {{fun<nm::EW_LEQ, uint8_t, uint8_t>, fun<nm::EW_LEQ, uint8_t, int8_t>, fun<nm::EW_LEQ, uint8_t, int16_t>, fun<nm::EW_LEQ, uint8_t, int32_t>, fun<nm::EW_LEQ, uint8_t, int64_t>, fun<nm::EW_LEQ, uint8_t, float32_t>, fun<nm::EW_LEQ, uint8_t, float64_t>, fun<nm::EW_LEQ, uint8_t, nm::Complex64>, fun<nm::EW_LEQ, uint8_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LEQ, int8_t, uint8_t>, fun<nm::EW_LEQ, int8_t, int8_t>, fun<nm::EW_LEQ, int8_t, int16_t>, fun<nm::EW_LEQ, int8_t, int32_t>, fun<nm::EW_LEQ, int8_t, int64_t>, fun<nm::EW_LEQ, int8_t, float32_t>, fun<nm::EW_LEQ, int8_t, float64_t>, fun<nm::EW_LEQ, int8_t, nm::Complex64>, fun<nm::EW_LEQ, int8_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LEQ, int16_t, uint8_t>, fun<nm::EW_LEQ, int16_t, int8_t>, fun<nm::EW_LEQ, int16_t, int16_t>, fun<nm::EW_LEQ, int16_t, int32_t>, fun<nm::EW_LEQ, int16_t, int64_t>, fun<nm::EW_LEQ, int16_t, float32_t>, fun<nm::EW_LEQ, int16_t, float64_t>, fun<nm::EW_LEQ, int16_t, nm::Complex64>, fun<nm::EW_LEQ, int16_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LEQ, int32_t, uint8_t>, fun<nm::EW_LEQ, int32_t, int8_t>, fun<nm::EW_LEQ, int32_t, int16_t>, fun<nm::EW_LEQ, int32_t, int32_t>, fun<nm::EW_LEQ, int32_t, int64_t>, fun<nm::EW_LEQ, int32_t, float32_t>, fun<nm::EW_LEQ, int32_t, float64_t>, fun<nm::EW_LEQ, int32_t, nm::Complex64>, fun<nm::EW_LEQ, int32_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LEQ, int64_t, uint8_t>, fun<nm::EW_LEQ, int64_t, int8_t>, fun<nm::EW_LEQ, int64_t, int16_t>, fun<nm::EW_LEQ, int64_t, int32_t>, fun<nm::EW_LEQ, int64_t, int64_t>, fun<nm::EW_LEQ, int64_t, float32_t>, fun<nm::EW_LEQ, int64_t, float64_t>, fun<nm::EW_LEQ, int64_t, nm::Complex64>, fun<nm::EW_LEQ, int64_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LEQ, float32_t, uint8_t>, fun<nm::EW_LEQ, float32_t, int8_t>, fun<nm::EW_LEQ, float32_t, int16_t>, fun<nm::EW_LEQ, float32_t, int32_t>, fun<nm::EW_LEQ, float32_t, int64_t>, fun<nm::EW_LEQ, float32_t, float32_t>, fun<nm::EW_LEQ, float32_t, float64_t>, fun<nm::EW_LEQ, float32_t, nm::Complex64>, fun<nm::EW_LEQ, float32_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LEQ, float64_t, uint8_t>, fun<nm::EW_LEQ, float64_t, int8_t>, fun<nm::EW_LEQ, float64_t, int16_t>, fun<nm::EW_LEQ, float64_t, int32_t>, fun<nm::EW_LEQ, float64_t, int64_t>, fun<nm::EW_LEQ, float64_t, float32_t>, fun<nm::EW_LEQ, float64_t, float64_t>, fun<nm::EW_LEQ, float64_t, nm::Complex64>, fun<nm::EW_LEQ, float64_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_LEQ, nm::Complex64, uint8_t>, fun<nm::EW_LEQ, nm::Complex64, int8_t>, fun<nm::EW_LEQ, nm::Complex64, int16_t>, fun<nm::EW_LEQ, nm::Complex64, int32_t>, fun<nm::EW_LEQ, nm::Complex64, int64_t>, fun<nm::EW_LEQ, nm::Complex64, float32_t>, fun<nm::EW_LEQ, nm::Complex64, float64_t>, fun<nm::EW_LEQ, nm::Complex64, nm::Complex64>, fun<nm::EW_LEQ, nm::Complex64, nm::Complex128>, NULL}, \
      {fun<nm::EW_LEQ, nm::Complex128, uint8_t>, fun<nm::EW_LEQ, nm::Complex128, int8_t>, fun<nm::EW_LEQ, nm::Complex128, int16_t>, fun<nm::EW_LEQ, nm::Complex128, int32_t>, fun<nm::EW_LEQ, nm::Complex128, int64_t>, fun<nm::EW_LEQ, nm::Complex128, float32_t>, fun<nm::EW_LEQ, nm::Complex128, float64_t>, fun<nm::EW_LEQ, nm::Complex128, nm::Complex64>, fun<nm::EW_LEQ, nm::Complex128, nm::Complex128>, NULL}, \
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_LEQ, nm::RubyObject, nm::RubyObject>}}, \
    {{fun<nm::EW_GEQ, uint8_t, uint8_t>, fun<nm::EW_GEQ, uint8_t, int8_t>, fun<nm::EW_GEQ, uint8_t, int16_t>, fun<nm::EW_GEQ, uint8_t, int32_t>, fun<nm::EW_GEQ, uint8_t, int64_t>, fun<nm::EW_GEQ, uint8_t, float32_t>, fun<nm::EW_GEQ, uint8_t, float64_t>, fun<nm::EW_GEQ, uint8_t, nm::Complex64>, fun<nm::EW_GEQ, uint8_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GEQ, int8_t, uint8_t>, fun<nm::EW_GEQ, int8_t, int8_t>, fun<nm::EW_GEQ, int8_t, int16_t>, fun<nm::EW_GEQ, int8_t, int32_t>, fun<nm::EW_GEQ, int8_t, int64_t>, fun<nm::EW_GEQ, int8_t, float32_t>, fun<nm::EW_GEQ, int8_t, float64_t>, fun<nm::EW_GEQ, int8_t, nm::Complex64>, fun<nm::EW_GEQ, int8_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GEQ, int16_t, uint8_t>, fun<nm::EW_GEQ, int16_t, int8_t>, fun<nm::EW_GEQ, int16_t, int16_t>, fun<nm::EW_GEQ, int16_t, int32_t>, fun<nm::EW_GEQ, int16_t, int64_t>, fun<nm::EW_GEQ, int16_t, float32_t>, fun<nm::EW_GEQ, int16_t, float64_t>, fun<nm::EW_GEQ, int16_t, nm::Complex64>, fun<nm::EW_GEQ, int16_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GEQ, int32_t, uint8_t>, fun<nm::EW_GEQ, int32_t, int8_t>, fun<nm::EW_GEQ, int32_t, int16_t>, fun<nm::EW_GEQ, int32_t, int32_t>, fun<nm::EW_GEQ, int32_t, int64_t>, fun<nm::EW_GEQ, int32_t, float32_t>, fun<nm::EW_GEQ, int32_t, float64_t>, fun<nm::EW_GEQ, int32_t, nm::Complex64>, fun<nm::EW_GEQ, int32_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GEQ, int64_t, uint8_t>, fun<nm::EW_GEQ, int64_t, int8_t>, fun<nm::EW_GEQ, int64_t, int16_t>, fun<nm::EW_GEQ, int64_t, int32_t>, fun<nm::EW_GEQ, int64_t, int64_t>, fun<nm::EW_GEQ, int64_t, float32_t>, fun<nm::EW_GEQ, int64_t, float64_t>, fun<nm::EW_GEQ, int64_t, nm::Complex64>, fun<nm::EW_GEQ, int64_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GEQ, float32_t, uint8_t>, fun<nm::EW_GEQ, float32_t, int8_t>, fun<nm::EW_GEQ, float32_t, int16_t>, fun<nm::EW_GEQ, float32_t, int32_t>, fun<nm::EW_GEQ, float32_t, int64_t>, fun<nm::EW_GEQ, float32_t, float32_t>, fun<nm::EW_GEQ, float32_t, float64_t>, fun<nm::EW_GEQ, float32_t, nm::Complex64>, fun<nm::EW_GEQ, float32_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GEQ, float64_t, uint8_t>, fun<nm::EW_GEQ, float64_t, int8_t>, fun<nm::EW_GEQ, float64_t, int16_t>, fun<nm::EW_GEQ, float64_t, int32_t>, fun<nm::EW_GEQ, float64_t, int64_t>, fun<nm::EW_GEQ, float64_t, float32_t>, fun<nm::EW_GEQ, float64_t, float64_t>, fun<nm::EW_GEQ, float64_t, nm::Complex64>, fun<nm::EW_GEQ, float64_t, nm::Complex128>, NULL}, \
      {fun<nm::EW_GEQ, nm::Complex64, uint8_t>, fun<nm::EW_GEQ, nm::Complex64, int8_t>, fun<nm::EW_GEQ, nm::Complex64, int16_t>, fun<nm::EW_GEQ, nm::Complex64, int32_t>, fun<nm::EW_GEQ, nm::Complex64, int64_t>, fun<nm::EW_GEQ, nm::Complex64, float32_t>, fun<nm::EW_GEQ, nm::Complex64, float64_t>, fun<nm::EW_GEQ, nm::Complex64, nm::Complex64>, fun<nm::EW_GEQ, nm::Complex64, nm::Complex128>, NULL}, \
      {fun<nm::EW_GEQ, nm::Complex128, uint8_t>, fun<nm::EW_GEQ, nm::Complex128, int8_t>, fun<nm::EW_GEQ, nm::Complex128, int16_t>, fun<nm::EW_GEQ, nm::Complex128, int32_t>, fun<nm::EW_GEQ, nm::Complex128, int64_t>, fun<nm::EW_GEQ, nm::Complex128, float32_t>, fun<nm::EW_GEQ, nm::Complex128, float64_t>, fun<nm::EW_GEQ, nm::Complex128, nm::Complex64>, fun<nm::EW_GEQ, nm::Complex128, nm::Complex128>, NULL}, \
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_GEQ, nm::RubyObject, nm::RubyObject>} \
    } \
  };

/*
 * Defines a static array that holds function pointers to an elementwise op,
 * itype, dtype templated versions of the specified function.
 */
#define OP_ITYPE_DTYPE_TEMPLATE_TABLE(fun, ret, ...) NAMED_OP_ITYPE_DTYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)

#define NAMED_OP_ITYPE_DTYPE_TEMPLATE_TABLE(name,  fun,  ret,  ...) \
  static ret (*(name)[nm::NUM_EWOPS][nm::NUM_ITYPES][nm::NUM_DTYPES])(__VA_ARGS__) = \
    {{{fun<nm::EW_ADD, uint8_t, uint8_t>,fun<nm::EW_ADD, uint8_t, int8_t>,fun<nm::EW_ADD, uint8_t, int16_t>,fun<nm::EW_ADD, uint8_t, int32_t>,fun<nm::EW_ADD, uint8_t, int64_t>,fun<nm::EW_ADD, uint8_t, float32_t>,fun<nm::EW_ADD, uint8_t, float64_t>,fun<nm::EW_ADD, uint8_t, nm::Complex64>,fun<nm::EW_ADD, uint8_t, nm::Complex128>,fun<nm::EW_ADD, uint8_t, nm::RubyObject>},\
{fun<nm::EW_ADD, uint16_t, uint8_t>,fun<nm::EW_ADD, uint16_t, int8_t>,fun<nm::EW_ADD, uint16_t, int16_t>,fun<nm::EW_ADD, uint16_t, int32_t>,fun<nm::EW_ADD, uint16_t, int64_t>,fun<nm::EW_ADD, uint16_t, float32_t>,fun<nm::EW_ADD, uint16_t, float64_t>,fun<nm::EW_ADD, uint16_t, nm::Complex64>,fun<nm::EW_ADD, uint16_t, nm::Complex128>,fun<nm::EW_ADD, uint16_t, nm::RubyObject>},\
{fun<nm::EW_ADD, uint32_t, uint8_t>,fun<nm::EW_ADD, uint32_t, int8_t>,fun<nm::EW_ADD, uint32_t, int16_t>,fun<nm::EW_ADD, uint32_t, int32_t>,fun<nm::EW_ADD, uint32_t, int64_t>,fun<nm::EW_ADD, uint32_t, float32_t>,fun<nm::EW_ADD, uint32_t, float64_t>,fun<nm::EW_ADD, uint32_t, nm::Complex64>,fun<nm::EW_ADD, uint32_t, nm::Complex128>,fun<nm::EW_ADD, uint32_t, nm::RubyObject>},\
{fun<nm::EW_ADD, uint64_t, uint8_t>,fun<nm::EW_ADD, uint64_t, int8_t>,fun<nm::EW_ADD, uint64_t, int16_t>,fun<nm::EW_ADD, uint64_t, int32_t>,fun<nm::EW_ADD, uint64_t, int64_t>,fun<nm::EW_ADD, uint64_t, float32_t>,fun<nm::EW_ADD, uint64_t, float64_t>,fun<nm::EW_ADD, uint64_t, nm::Complex64>,fun<nm::EW_ADD, uint64_t, nm::Complex128>,fun<nm::EW_ADD, uint64_t, nm::RubyObject>}},\
{{fun<nm::EW_SUB, uint8_t, uint8_t>,fun<nm::EW_SUB, uint8_t, int8_t>,fun<nm::EW_SUB, uint8_t, int16_t>,fun<nm::EW_SUB, uint8_t, int32_t>,fun<nm::EW_SUB, uint8_t, int64_t>,fun<nm::EW_SUB, uint8_t, float32_t>,fun<nm::EW_SUB, uint8_t, float64_t>,fun<nm::EW_SUB, uint8_t, nm::Complex64>,fun<nm::EW_SUB, uint8_t, nm::Complex128>,fun<nm::EW_SUB, uint8_t, nm::RubyObject>},\
{fun<nm::EW_SUB, uint16_t, uint8_t>,fun<nm::EW_SUB, uint16_t, int8_t>,fun<nm::EW_SUB, uint16_t, int16_t>,fun<nm::EW_SUB, uint16_t, int32_t>,fun<nm::EW_SUB, uint16_t, int64_t>,fun<nm::EW_SUB, uint16_t, float32_t>,fun<nm::EW_SUB, uint16_t, float64_t>,fun<nm::EW_SUB, uint16_t, nm::Complex64>,fun<nm::EW_SUB, uint16_t, nm::Complex128>,fun<nm::EW_SUB, uint16_t, nm::RubyObject>},\
{fun<nm::EW_SUB, uint32_t, uint8_t>,fun<nm::EW_SUB, uint32_t, int8_t>,fun<nm::EW_SUB, uint32_t, int16_t>,fun<nm::EW_SUB, uint32_t, int32_t>,fun<nm::EW_SUB, uint32_t, int64_t>,fun<nm::EW_SUB, uint32_t, float32_t>,fun<nm::EW_SUB, uint32_t, float64_t>,fun<nm::EW_SUB, uint32_t, nm::Complex64>,fun<nm::EW_SUB, uint32_t, nm::Complex128>,fun<nm::EW_SUB, uint32_t, nm::RubyObject>},\
{fun<nm::EW_SUB, uint64_t, uint8_t>,fun<nm::EW_SUB, uint64_t, int8_t>,fun<nm::EW_SUB, uint64_t, int16_t>,fun<nm::EW_SUB, uint64_t, int32_t>,fun<nm::EW_SUB, uint64_t, int64_t>,fun<nm::EW_SUB, uint64_t, float32_t>,fun<nm::EW_SUB, uint64_t, float64_t>,fun<nm::EW_SUB, uint64_t, nm::Complex64>,fun<nm::EW_SUB, uint64_t, nm::Complex128>,fun<nm::EW_SUB, uint64_t, nm::RubyObject>}},\
{{fun<nm::EW_MUL, uint8_t, uint8_t>,fun<nm::EW_MUL, uint8_t, int8_t>,fun<nm::EW_MUL, uint8_t, int16_t>,fun<nm::EW_MUL, uint8_t, int32_t>,fun<nm::EW_MUL, uint8_t, int64_t>,fun<nm::EW_MUL, uint8_t, float32_t>,fun<nm::EW_MUL, uint8_t, float64_t>,fun<nm::EW_MUL, uint8_t, nm::Complex64>,fun<nm::EW_MUL, uint8_t, nm::Complex128>,fun<nm::EW_MUL, uint8_t, nm::RubyObject>},\
{fun<nm::EW_MUL, uint16_t, uint8_t>,fun<nm::EW_MUL, uint16_t, int8_t>,fun<nm::EW_MUL, uint16_t, int16_t>,fun<nm::EW_MUL, uint16_t, int32_t>,fun<nm::EW_MUL, uint16_t, int64_t>,fun<nm::EW_MUL, uint16_t, float32_t>,fun<nm::EW_MUL, uint16_t, float64_t>,fun<nm::EW_MUL, uint16_t, nm::Complex64>,fun<nm::EW_MUL, uint16_t, nm::Complex128>,fun<nm::EW_MUL, uint16_t, nm::RubyObject>},\
{fun<nm::EW_MUL, uint32_t, uint8_t>,fun<nm::EW_MUL, uint32_t, int8_t>,fun<nm::EW_MUL, uint32_t, int16_t>,fun<nm::EW_MUL, uint32_t, int32_t>,fun<nm::EW_MUL, uint32_t, int64_t>,fun<nm::EW_MUL, uint32_t, float32_t>,fun<nm::EW_MUL, uint32_t, float64_t>,fun<nm::EW_MUL, uint32_t, nm::Complex64>,fun<nm::EW_MUL, uint32_t, nm::Complex128>,fun<nm::EW_MUL, uint32_t, nm::RubyObject>},\
{fun<nm::EW_MUL, uint64_t, uint8_t>,fun<nm::EW_MUL, uint64_t, int8_t>,fun<nm::EW_MUL, uint64_t, int16_t>,fun<nm::EW_MUL, uint64_t, int32_t>,fun<nm::EW_MUL, uint64_t, int64_t>,fun<nm::EW_MUL, uint64_t, float32_t>,fun<nm::EW_MUL, uint64_t, float64_t>,fun<nm::EW_MUL, uint64_t, nm::Complex64>,fun<nm::EW_MUL, uint64_t, nm::Complex128>,fun<nm::EW_MUL, uint64_t, nm::RubyObject>}},\
{{fun<nm::EW_DIV, uint8_t, uint8_t>,fun<nm::EW_DIV, uint8_t, int8_t>,fun<nm::EW_DIV, uint8_t, int16_t>,fun<nm::EW_DIV, uint8_t, int32_t>,fun<nm::EW_DIV, uint8_t, int64_t>,fun<nm::EW_DIV, uint8_t, float32_t>,fun<nm::EW_DIV, uint8_t, float64_t>,fun<nm::EW_DIV, uint8_t, nm::Complex64>,fun<nm::EW_DIV, uint8_t, nm::Complex128>,fun<nm::EW_DIV, uint8_t, nm::RubyObject>},\
{fun<nm::EW_DIV, uint16_t, uint8_t>,fun<nm::EW_DIV, uint16_t, int8_t>,fun<nm::EW_DIV, uint16_t, int16_t>,fun<nm::EW_DIV, uint16_t, int32_t>,fun<nm::EW_DIV, uint16_t, int64_t>,fun<nm::EW_DIV, uint16_t, float32_t>,fun<nm::EW_DIV, uint16_t, float64_t>,fun<nm::EW_DIV, uint16_t, nm::Complex64>,fun<nm::EW_DIV, uint16_t, nm::Complex128>,fun<nm::EW_DIV, uint16_t, nm::RubyObject>},\
{fun<nm::EW_DIV, uint32_t, uint8_t>,fun<nm::EW_DIV, uint32_t, int8_t>,fun<nm::EW_DIV, uint32_t, int16_t>,fun<nm::EW_DIV, uint32_t, int32_t>,fun<nm::EW_DIV, uint32_t, int64_t>,fun<nm::EW_DIV, uint32_t, float32_t>,fun<nm::EW_DIV, uint32_t, float64_t>,fun<nm::EW_DIV, uint32_t, nm::Complex64>,fun<nm::EW_DIV, uint32_t, nm::Complex128>,fun<nm::EW_DIV, uint32_t, nm::RubyObject>},\
{fun<nm::EW_DIV, uint64_t, uint8_t>,fun<nm::EW_DIV, uint64_t, int8_t>,fun<nm::EW_DIV, uint64_t, int16_t>,fun<nm::EW_DIV, uint64_t, int32_t>,fun<nm::EW_DIV, uint64_t, int64_t>,fun<nm::EW_DIV, uint64_t, float32_t>,fun<nm::EW_DIV, uint64_t, float64_t>,fun<nm::EW_DIV, uint64_t, nm::Complex64>,fun<nm::EW_DIV, uint64_t, nm::Complex128>,fun<nm::EW_DIV, uint64_t, nm::RubyObject>}},\
{{fun<nm::EW_MOD, uint8_t, uint8_t>,fun<nm::EW_MOD, uint8_t, int8_t>,fun<nm::EW_MOD, uint8_t, int16_t>,fun<nm::EW_MOD, uint8_t, int32_t>,fun<nm::EW_MOD, uint8_t, int64_t>,fun<nm::EW_MOD, uint8_t, float32_t>,fun<nm::EW_MOD, uint8_t, float64_t>,fun<nm::EW_MOD, uint8_t, nm::Complex64>,fun<nm::EW_MOD, uint8_t, nm::Complex128>,fun<nm::EW_MOD, uint8_t, nm::RubyObject>},\
{fun<nm::EW_MOD, uint16_t, uint8_t>,fun<nm::EW_MOD, uint16_t, int8_t>,fun<nm::EW_MOD, uint16_t, int16_t>,fun<nm::EW_MOD, uint16_t, int32_t>,fun<nm::EW_MOD, uint16_t, int64_t>,fun<nm::EW_MOD, uint16_t, float32_t>,fun<nm::EW_MOD, uint16_t, float64_t>,fun<nm::EW_MOD, uint16_t, nm::Complex64>,fun<nm::EW_MOD, uint16_t, nm::Complex128>,fun<nm::EW_MOD, uint16_t, nm::RubyObject>},\
{fun<nm::EW_MOD, uint32_t, uint8_t>,fun<nm::EW_MOD, uint32_t, int8_t>,fun<nm::EW_MOD, uint32_t, int16_t>,fun<nm::EW_MOD, uint32_t, int32_t>,fun<nm::EW_MOD, uint32_t, int64_t>,fun<nm::EW_MOD, uint32_t, float32_t>,fun<nm::EW_MOD, uint32_t, float64_t>,fun<nm::EW_MOD, uint32_t, nm::Complex64>,fun<nm::EW_MOD, uint32_t, nm::Complex128>,fun<nm::EW_MOD, uint32_t, nm::RubyObject>},\
{fun<nm::EW_MOD, uint64_t, uint8_t>,fun<nm::EW_MOD, uint64_t, int8_t>,fun<nm::EW_MOD, uint64_t, int16_t>,fun<nm::EW_MOD, uint64_t, int32_t>,fun<nm::EW_MOD, uint64_t, int64_t>,fun<nm::EW_MOD, uint64_t, float32_t>,fun<nm::EW_MOD, uint64_t, float64_t>,fun<nm::EW_MOD, uint64_t, nm::Complex64>,fun<nm::EW_MOD, uint64_t, nm::Complex128>,fun<nm::EW_MOD, uint64_t, nm::RubyObject>}},\
{{fun<nm::EW_EQEQ, uint8_t, uint8_t>,fun<nm::EW_EQEQ, uint8_t, int8_t>,fun<nm::EW_EQEQ, uint8_t, int16_t>,fun<nm::EW_EQEQ, uint8_t, int32_t>,fun<nm::EW_EQEQ, uint8_t, int64_t>,fun<nm::EW_EQEQ, uint8_t, float32_t>,fun<nm::EW_EQEQ, uint8_t, float64_t>,fun<nm::EW_EQEQ, uint8_t, nm::Complex64>,fun<nm::EW_EQEQ, uint8_t, nm::Complex128>,fun<nm::EW_EQEQ, uint8_t, nm::RubyObject>},\
{fun<nm::EW_EQEQ, uint16_t, uint8_t>,fun<nm::EW_EQEQ, uint16_t, int8_t>,fun<nm::EW_EQEQ, uint16_t, int16_t>,fun<nm::EW_EQEQ, uint16_t, int32_t>,fun<nm::EW_EQEQ, uint16_t, int64_t>,fun<nm::EW_EQEQ, uint16_t, float32_t>,fun<nm::EW_EQEQ, uint16_t, float64_t>,fun<nm::EW_EQEQ, uint16_t, nm::Complex64>,fun<nm::EW_EQEQ, uint16_t, nm::Complex128>,fun<nm::EW_EQEQ, uint16_t, nm::RubyObject>},\
{fun<nm::EW_EQEQ, uint32_t, uint8_t>,fun<nm::EW_EQEQ, uint32_t, int8_t>,fun<nm::EW_EQEQ, uint32_t, int16_t>,fun<nm::EW_EQEQ, uint32_t, int32_t>,fun<nm::EW_EQEQ, uint32_t, int64_t>,fun<nm::EW_EQEQ, uint32_t, float32_t>,fun<nm::EW_EQEQ, uint32_t, float64_t>,fun<nm::EW_EQEQ, uint32_t, nm::Complex64>,fun<nm::EW_EQEQ, uint32_t, nm::Complex128>,fun<nm::EW_EQEQ, uint32_t, nm::RubyObject>},\
{fun<nm::EW_EQEQ, uint64_t, uint8_t>,fun<nm::EW_EQEQ, uint64_t, int8_t>,fun<nm::EW_EQEQ, uint64_t, int16_t>,fun<nm::EW_EQEQ, uint64_t, int32_t>,fun<nm::EW_EQEQ, uint64_t, int64_t>,fun<nm::EW_EQEQ, uint64_t, float32_t>,fun<nm::EW_EQEQ, uint64_t, float64_t>,fun<nm::EW_EQEQ, uint64_t, nm::Complex64>,fun<nm::EW_EQEQ, uint64_t, nm::Complex128>,fun<nm::EW_EQEQ, uint64_t, nm::RubyObject>}},\
{{fun<nm::EW_NEQ, uint8_t, uint8_t>,fun<nm::EW_NEQ, uint8_t, int8_t>,fun<nm::EW_NEQ, uint8_t, int16_t>,fun<nm::EW_NEQ, uint8_t, int32_t>,fun<nm::EW_NEQ, uint8_t, int64_t>,fun<nm::EW_NEQ, uint8_t, float32_t>,fun<nm::EW_NEQ, uint8_t, float64_t>,fun<nm::EW_NEQ, uint8_t, nm::Complex64>,fun<nm::EW_NEQ, uint8_t, nm::Complex128>,fun<nm::EW_NEQ, uint8_t, nm::RubyObject>},\
{fun<nm::EW_NEQ, uint16_t, uint8_t>,fun<nm::EW_NEQ, uint16_t, int8_t>,fun<nm::EW_NEQ, uint16_t, int16_t>,fun<nm::EW_NEQ, uint16_t, int32_t>,fun<nm::EW_NEQ, uint16_t, int64_t>,fun<nm::EW_NEQ, uint16_t, float32_t>,fun<nm::EW_NEQ, uint16_t, float64_t>,fun<nm::EW_NEQ, uint16_t, nm::Complex64>,fun<nm::EW_NEQ, uint16_t, nm::Complex128>,fun<nm::EW_NEQ, uint16_t, nm::RubyObject>},\
{fun<nm::EW_NEQ, uint32_t, uint8_t>,fun<nm::EW_NEQ, uint32_t, int8_t>,fun<nm::EW_NEQ, uint32_t, int16_t>,fun<nm::EW_NEQ, uint32_t, int32_t>,fun<nm::EW_NEQ, uint32_t, int64_t>,fun<nm::EW_NEQ, uint32_t, float32_t>,fun<nm::EW_NEQ, uint32_t, float64_t>,fun<nm::EW_NEQ, uint32_t, nm::Complex64>,fun<nm::EW_NEQ, uint32_t, nm::Complex128>,fun<nm::EW_NEQ, uint32_t, nm::RubyObject>},\
{fun<nm::EW_NEQ, uint64_t, uint8_t>,fun<nm::EW_NEQ, uint64_t, int8_t>,fun<nm::EW_NEQ, uint64_t, int16_t>,fun<nm::EW_NEQ, uint64_t, int32_t>,fun<nm::EW_NEQ, uint64_t, int64_t>,fun<nm::EW_NEQ, uint64_t, float32_t>,fun<nm::EW_NEQ, uint64_t, float64_t>,fun<nm::EW_NEQ, uint64_t, nm::Complex64>,fun<nm::EW_NEQ, uint64_t, nm::Complex128>,fun<nm::EW_NEQ, uint64_t, nm::RubyObject>}},\
{{fun<nm::EW_LT, uint8_t, uint8_t>,fun<nm::EW_LT, uint8_t, int8_t>,fun<nm::EW_LT, uint8_t, int16_t>,fun<nm::EW_LT, uint8_t, int32_t>,fun<nm::EW_LT, uint8_t, int64_t>,fun<nm::EW_LT, uint8_t, float32_t>,fun<nm::EW_LT, uint8_t, float64_t>,fun<nm::EW_LT, uint8_t, nm::Complex64>,fun<nm::EW_LT, uint8_t, nm::Complex128>,fun<nm::EW_LT, uint8_t, nm::RubyObject>},\
{fun<nm::EW_LT, uint16_t, uint8_t>,fun<nm::EW_LT, uint16_t, int8_t>,fun<nm::EW_LT, uint16_t, int16_t>,fun<nm::EW_LT, uint16_t, int32_t>,fun<nm::EW_LT, uint16_t, int64_t>,fun<nm::EW_LT, uint16_t, float32_t>,fun<nm::EW_LT, uint16_t, float64_t>,fun<nm::EW_LT, uint16_t, nm::Complex64>,fun<nm::EW_LT, uint16_t, nm::Complex128>,fun<nm::EW_LT, uint16_t, nm::RubyObject>},\
{fun<nm::EW_LT, uint32_t, uint8_t>,fun<nm::EW_LT, uint32_t, int8_t>,fun<nm::EW_LT, uint32_t, int16_t>,fun<nm::EW_LT, uint32_t, int32_t>,fun<nm::EW_LT, uint32_t, int64_t>,fun<nm::EW_LT, uint32_t, float32_t>,fun<nm::EW_LT, uint32_t, float64_t>,fun<nm::EW_LT, uint32_t, nm::Complex64>,fun<nm::EW_LT, uint32_t, nm::Complex128>,fun<nm::EW_LT, uint32_t, nm::RubyObject>},\
{fun<nm::EW_LT, uint64_t, uint8_t>,fun<nm::EW_LT, uint64_t, int8_t>,fun<nm::EW_LT, uint64_t, int16_t>,fun<nm::EW_LT, uint64_t, int32_t>,fun<nm::EW_LT, uint64_t, int64_t>,fun<nm::EW_LT, uint64_t, float32_t>,fun<nm::EW_LT, uint64_t, float64_t>,fun<nm::EW_LT, uint64_t, nm::Complex64>,fun<nm::EW_LT, uint64_t, nm::Complex128>,fun<nm::EW_LT, uint64_t, nm::RubyObject>}},\
{{fun<nm::EW_GT, uint8_t, uint8_t>,fun<nm::EW_GT, uint8_t, int8_t>,fun<nm::EW_GT, uint8_t, int16_t>,fun<nm::EW_GT, uint8_t, int32_t>,fun<nm::EW_GT, uint8_t, int64_t>,fun<nm::EW_GT, uint8_t, float32_t>,fun<nm::EW_GT, uint8_t, float64_t>,fun<nm::EW_GT, uint8_t, nm::Complex64>,fun<nm::EW_GT, uint8_t, nm::Complex128>,fun<nm::EW_GT, uint8_t, nm::RubyObject>},\
{fun<nm::EW_GT, uint16_t, uint8_t>,fun<nm::EW_GT, uint16_t, int8_t>,fun<nm::EW_GT, uint16_t, int16_t>,fun<nm::EW_GT, uint16_t, int32_t>,fun<nm::EW_GT, uint16_t, int64_t>,fun<nm::EW_GT, uint16_t, float32_t>,fun<nm::EW_GT, uint16_t, float64_t>,fun<nm::EW_GT, uint16_t, nm::Complex64>,fun<nm::EW_GT, uint16_t, nm::Complex128>,fun<nm::EW_GT, uint16_t, nm::RubyObject>},\
{fun<nm::EW_GT, uint32_t, uint8_t>,fun<nm::EW_GT, uint32_t, int8_t>,fun<nm::EW_GT, uint32_t, int16_t>,fun<nm::EW_GT, uint32_t, int32_t>,fun<nm::EW_GT, uint32_t, int64_t>,fun<nm::EW_GT, uint32_t, float32_t>,fun<nm::EW_GT, uint32_t, float64_t>,fun<nm::EW_GT, uint32_t, nm::Complex64>,fun<nm::EW_GT, uint32_t, nm::Complex128>,fun<nm::EW_GT, uint32_t, nm::RubyObject>},\
{fun<nm::EW_GT, uint64_t, uint8_t>,fun<nm::EW_GT, uint64_t, int8_t>,fun<nm::EW_GT, uint64_t, int16_t>,fun<nm::EW_GT, uint64_t, int32_t>,fun<nm::EW_GT, uint64_t, int64_t>,fun<nm::EW_GT, uint64_t, float32_t>,fun<nm::EW_GT, uint64_t, float64_t>,fun<nm::EW_GT, uint64_t, nm::Complex64>,fun<nm::EW_GT, uint64_t, nm::Complex128>,fun<nm::EW_GT, uint64_t, nm::RubyObject>}},\
{{fun<nm::EW_LEQ, uint8_t, uint8_t>,fun<nm::EW_LEQ, uint8_t, int8_t>,fun<nm::EW_LEQ, uint8_t, int16_t>,fun<nm::EW_LEQ, uint8_t, int32_t>,fun<nm::EW_LEQ, uint8_t, int64_t>,fun<nm::EW_LEQ, uint8_t, float32_t>,fun<nm::EW_LEQ, uint8_t, float64_t>,fun<nm::EW_LEQ, uint8_t, nm::Complex64>,fun<nm::EW_LEQ, uint8_t, nm::Complex128>,fun<nm::EW_LEQ, uint8_t, nm::RubyObject>},\
{fun<nm::EW_LEQ, uint16_t, uint8_t>,fun<nm::EW_LEQ, uint16_t, int8_t>,fun<nm::EW_LEQ, uint16_t, int16_t>,fun<nm::EW_LEQ, uint16_t, int32_t>,fun<nm::EW_LEQ, uint16_t, int64_t>,fun<nm::EW_LEQ, uint16_t, float32_t>,fun<nm::EW_LEQ, uint16_t, float64_t>,fun<nm::EW_LEQ, uint16_t, nm::Complex64>,fun<nm::EW_LEQ, uint16_t, nm::Complex128>,fun<nm::EW_LEQ, uint16_t, nm::RubyObject>},\
{fun<nm::EW_LEQ, uint32_t, uint8_t>,fun<nm::EW_LEQ, uint32_t, int8_t>,fun<nm::EW_LEQ, uint32_t, int16_t>,fun<nm::EW_LEQ, uint32_t, int32_t>,fun<nm::EW_LEQ, uint32_t, int64_t>,fun<nm::EW_LEQ, uint32_t, float32_t>,fun<nm::EW_LEQ, uint32_t, float64_t>,fun<nm::EW_LEQ, uint32_t, nm::Complex64>,fun<nm::EW_LEQ, uint32_t, nm::Complex128>,fun<nm::EW_LEQ, uint32_t, nm::RubyObject>},\
{fun<nm::EW_LEQ, uint64_t, uint8_t>,fun<nm::EW_LEQ, uint64_t, int8_t>,fun<nm::EW_LEQ, uint64_t, int16_t>,fun<nm::EW_LEQ, uint64_t, int32_t>,fun<nm::EW_LEQ, uint64_t, int64_t>,fun<nm::EW_LEQ, uint64_t, float32_t>,fun<nm::EW_LEQ, uint64_t, float64_t>,fun<nm::EW_LEQ, uint64_t, nm::Complex64>,fun<nm::EW_LEQ, uint64_t, nm::Complex128>,fun<nm::EW_LEQ, uint64_t, nm::RubyObject>}},\
{{fun<nm::EW_GEQ, uint8_t, uint8_t>,fun<nm::EW_GEQ, uint8_t, int8_t>,fun<nm::EW_GEQ, uint8_t, int16_t>,fun<nm::EW_GEQ, uint8_t, int32_t>,fun<nm::EW_GEQ, uint8_t, int64_t>,fun<nm::EW_GEQ, uint8_t, float32_t>,fun<nm::EW_GEQ, uint8_t, float64_t>,fun<nm::EW_GEQ, uint8_t, nm::Complex64>,fun<nm::EW_GEQ, uint8_t, nm::Complex128>,fun<nm::EW_GEQ, uint8_t, nm::RubyObject>},\
{fun<nm::EW_GEQ, uint16_t, uint8_t>,fun<nm::EW_GEQ, uint16_t, int8_t>,fun<nm::EW_GEQ, uint16_t, int16_t>,fun<nm::EW_GEQ, uint16_t, int32_t>,fun<nm::EW_GEQ, uint16_t, int64_t>,fun<nm::EW_GEQ, uint16_t, float32_t>,fun<nm::EW_GEQ, uint16_t, float64_t>,fun<nm::EW_GEQ, uint16_t, nm::Complex64>,fun<nm::EW_GEQ, uint16_t, nm::Complex128>,fun<nm::EW_GEQ, uint16_t, nm::RubyObject>},\
{fun<nm::EW_GEQ, uint32_t, uint8_t>,fun<nm::EW_GEQ, uint32_t, int8_t>,fun<nm::EW_GEQ, uint32_t, int16_t>,fun<nm::EW_GEQ, uint32_t, int32_t>,fun<nm::EW_GEQ, uint32_t, int64_t>,fun<nm::EW_GEQ, uint32_t, float32_t>,fun<nm::EW_GEQ, uint32_t, float64_t>,fun<nm::EW_GEQ, uint32_t, nm::Complex64>,fun<nm::EW_GEQ, uint32_t, nm::Complex128>,fun<nm::EW_GEQ, uint32_t, nm::RubyObject>},\
{fun<nm::EW_GEQ, uint64_t, uint8_t>,fun<nm::EW_GEQ, uint64_t, int8_t>,fun<nm::EW_GEQ, uint64_t, int16_t>,fun<nm::EW_GEQ, uint64_t, int32_t>,fun<nm::EW_GEQ, uint64_t, int64_t>,fun<nm::EW_GEQ, uint64_t, float32_t>,fun<nm::EW_GEQ, uint64_t, float64_t>,fun<nm::EW_GEQ, uint64_t, nm::Complex64>,fun<nm::EW_GEQ, uint64_t, nm::Complex128>,fun<nm::EW_GEQ, uint64_t, nm::RubyObject>}}};


extern "C" {


/*
 * Data
 */

// regular data types
extern const char* const  DTYPE_NAMES[nm::NUM_DTYPES];
extern const size_t       DTYPE_SIZES[nm::NUM_DTYPES];

extern const nm::dtype_t Upcast[nm::NUM_DTYPES][nm::NUM_DTYPES];


/*
 * Functions
 */


void*            rubyobj_to_cval(VALUE val, nm::dtype_t dtype);
void            rubyval_to_cval(VALUE val, nm::dtype_t dtype, void* loc);

void nm_init_data();

} // end of extern "C" block

#endif // DATA_H
