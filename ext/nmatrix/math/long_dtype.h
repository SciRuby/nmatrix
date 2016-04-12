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
// SciRuby is Copyright (c) 2010 - present, Ruby Science Foundation
// NMatrix is Copyright (c) 2012 - present, John Woods and the Ruby Science Foundation
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
// == long_dtype.h
//
// Declarations necessary for the native versions of GEMM and GEMV,
// as well as for IMAX.
//

#ifndef LONG_DTYPE_H
#define LONG_DTYPE_H

namespace nm { namespace math {
  // These allow an increase in precision for intermediate values of gemm and gemv.
  // See also: http://stackoverflow.com/questions/11873694/how-does-one-increase-precision-in-c-templates-in-a-template-typename-dependen
  template <typename DType> struct LongDType;
  template <> struct LongDType<uint8_t> { typedef int16_t type; };
  template <> struct LongDType<int8_t> { typedef int16_t type; };
  template <> struct LongDType<int16_t> { typedef int32_t type; };
  template <> struct LongDType<int32_t> { typedef int64_t type; };
  template <> struct LongDType<int64_t> { typedef int64_t type; };
  template <> struct LongDType<float> { typedef double type; };
  template <> struct LongDType<double> { typedef double type; };
  template <> struct LongDType<Complex64> { typedef Complex128 type; };
  template <> struct LongDType<Complex128> { typedef Complex128 type; };
  template <> struct LongDType<RubyObject> { typedef RubyObject type; };

  template <typename DType> struct MagnitudeDType;
  template <> struct MagnitudeDType<uint8_t> { typedef uint8_t type; };
  template <> struct MagnitudeDType<int8_t> { typedef int8_t type; };
  template <> struct MagnitudeDType<int16_t> { typedef int16_t type; };
  template <> struct MagnitudeDType<int32_t> { typedef int32_t type; };
  template <> struct MagnitudeDType<int64_t> { typedef int64_t type; };
  template <> struct MagnitudeDType<float> { typedef float type; };
  template <> struct MagnitudeDType<double> { typedef double type; };
  template <> struct MagnitudeDType<Complex64> { typedef float type; };
  template <> struct MagnitudeDType<Complex128> { typedef double type; };
  template <> struct MagnitudeDType<RubyObject> { typedef RubyObject type; };
  
}} // end of namespace nm::math

#endif
