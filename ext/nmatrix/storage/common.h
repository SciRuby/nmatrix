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
// == common.h
//
// Header file for code common to all storage types.

#ifndef STORAGE_COMMON_H
#define STORAGE_COMMON_H

/*
 * Standard Includes
 */

#include <cmath> // pow().

/*
 * Project Includes
 */
#include "data/data.h"
#include "nmatrix.h"

/*
 * Macros
 */

extern "C" {

/*
 * Types
 */

// For binary operations involving matrices that need to be casted.
struct STORAGE_PAIR {
  STORAGE* left;
  STORAGE* right;
};

struct SLICE {
  size_t*	coords; // Coordinate of first element
  size_t*	lengths; // Lengths of slice
  bool  	single; // true if all lengths equal to 1 (represents single matrix element)
};

/*
 * Data
 */

/*
 * Functions
 */

  size_t nm_storage_count_max_elements(const STORAGE* storage);
  VALUE nm_enumerator_length(VALUE nmatrix);

} // end of extern "C" block

namespace nm {

  /*
   * Templated helper function for element-wise operations, used by dense, yale, and list.
   */
  template <ewop_t op, typename LDType, typename RDType>
  inline VALUE ew_op_switch(LDType left, RDType right) {
    switch (op) {
      case EW_ADD:
        return RubyObject(left + right).rval;

      case EW_SUB:
        return RubyObject(left - right).rval;

      case EW_MUL:
        return RubyObject(left * right).rval;

      case EW_DIV:
        return RubyObject(left / right).rval;

      case EW_POW:
        return RubyObject(pow(left, right)).rval;

      case EW_MOD:
        rb_raise(rb_eNotImpError, "Element-wise modulo is currently not supported.");
        break;

      default:
        rb_raise(rb_eStandardError, "This should not happen.");
    }
    return Qnil;
  }

  #define EWOP_INT_INT_DIV(ltype, rtype)       template <>       \
  inline VALUE ew_op_switch<EW_DIV>( ltype left, rtype right) { \
    if (right == 0) rb_raise(rb_eZeroDivError, "cannot divide type by 0, would throw SIGFPE");  \
    if ((left > 0 && right > 0) || (left < 0 && right < 0)) \
      return left / right;  \
    else \
      return ( ltype )(std::floor((double)(left) / (double)(right)));  \
  }

  #define EWOP_UINT_UINT_DIV(ltype, rtype)       template <>       \
  inline VALUE ew_op_switch<EW_DIV>( ltype left, rtype right) { \
    if (right == 0) rb_raise(rb_eZeroDivError, "cannot divide type by 0, would throw SIGFPE");  \
    return left / right;  \
  }

  #define EWOP_INT_UINT_DIV(ltype, rtype)       template <>       \
  inline VALUE ew_op_switch<EW_DIV>( ltype left, rtype right) { \
    if (right == 0) rb_raise(rb_eZeroDivError, "cannot divide type by 0, would throw SIGFPE");  \
    if (left > 0 )  return left / right;  \
    else            return ( ltype )(std::floor((double)(left) / (double)(right)));  \
  }

  #define EWOP_UINT_INT_DIV(ltype, rtype)       template <>       \
  inline VALUE ew_op_switch<EW_DIV>( ltype left, rtype right) { \
    if (right == 0) rb_raise(rb_eZeroDivError, "cannot divide type by 0, would throw SIGFPE");  \
    if (right > 0)  return left / right;  \
    else            return ( ltype )(std::floor((double)(left) / (double)(right)));  \
  }

  #define EWOP_FLOAT_INT_DIV(ltype, rtype)       template <>       \
  inline VALUE ew_op_switch<EW_DIV>( ltype left, rtype right) { \
    return left / (ltype)(right);  \
  }

  // Ensure that divisions are done in the Ruby way, and that (int)x/0 always raises a Ruby error instead
  // of throwing a SIGFPE.
  EWOP_INT_INT_DIV(int64_t, int64_t)
  EWOP_INT_INT_DIV(int32_t, int32_t)
  EWOP_INT_INT_DIV(int32_t, int64_t)
  EWOP_INT_INT_DIV(int16_t, int16_t)
  EWOP_INT_INT_DIV(int16_t, int32_t)
  EWOP_INT_INT_DIV(int16_t, int64_t)
  EWOP_INT_INT_DIV(int8_t, int8_t)
  EWOP_INT_UINT_DIV(int8_t, u_int8_t)
  EWOP_INT_INT_DIV(int8_t, int16_t)
  EWOP_INT_INT_DIV(int8_t, int32_t)
  EWOP_INT_INT_DIV(int8_t, int64_t)
  EWOP_UINT_UINT_DIV(uint8_t, uint8_t)
  EWOP_UINT_INT_DIV(uint8_t, int8_t)
  EWOP_UINT_INT_DIV(uint8_t, int16_t)
  EWOP_UINT_INT_DIV(uint8_t, int32_t)
  EWOP_UINT_INT_DIV(uint8_t, int64_t)
  EWOP_FLOAT_INT_DIV(float, int8_t)
  EWOP_FLOAT_INT_DIV(float, u_int8_t)
  EWOP_FLOAT_INT_DIV(float, int16_t)
  EWOP_FLOAT_INT_DIV(float, int32_t)
  EWOP_FLOAT_INT_DIV(float, int64_t)
  EWOP_FLOAT_INT_DIV(double, int8_t)
  EWOP_FLOAT_INT_DIV(double, u_int8_t)
  EWOP_FLOAT_INT_DIV(double, int16_t)
  EWOP_FLOAT_INT_DIV(double, int32_t)
  EWOP_FLOAT_INT_DIV(double, int64_t)

}

#endif // STORAGE_COMMON_H
