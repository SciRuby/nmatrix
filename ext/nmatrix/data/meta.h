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
// == meta.h
//
// Header file for dealing with template metaprogramming.

#ifndef META_H
# define META_H

namespace nm {
  /*
   * Template Metaprogramming
   */
  template <typename T> struct ctype_to_dtype_enum {
    static const nm::dtype_t value_type = nm::BYTE;
  };
  template <> struct ctype_to_dtype_enum<uint8_t> { static const nm::dtype_t value_type = nm::BYTE; };
  template <> struct ctype_to_dtype_enum<int8_t>  { static const nm::dtype_t value_type = nm::INT8; };
  template <> struct ctype_to_dtype_enum<int16_t> { static const nm::dtype_t value_type = nm::INT16; };
  template <> struct ctype_to_dtype_enum<int32_t> { static const nm::dtype_t value_type = nm::INT32; };
  template <> struct ctype_to_dtype_enum<int64_t> { static const nm::dtype_t value_type = nm::INT64; };
  template <> struct ctype_to_dtype_enum<float>   { static const nm::dtype_t value_type = nm::FLOAT32; };
  template <> struct ctype_to_dtype_enum<double>  { static const nm::dtype_t value_type = nm::FLOAT64; };
  template <> struct ctype_to_dtype_enum<Complex64>   { static const nm::dtype_t value_type = nm::COMPLEX64; };
  template <> struct ctype_to_dtype_enum<Complex128>  { static const nm::dtype_t value_type = nm::COMPLEX128; };
  template <> struct ctype_to_dtype_enum<Rational32>   { static const nm::dtype_t value_type = nm::RATIONAL32; };
  template <> struct ctype_to_dtype_enum<Rational64>  { static const nm::dtype_t value_type = nm::RATIONAL64; };
  template <> struct ctype_to_dtype_enum<Rational128>  { static const nm::dtype_t value_type = nm::RATIONAL128; };
  template <> struct ctype_to_dtype_enum<RubyObject>  { static const nm::dtype_t value_type = nm::RUBYOBJ; };


  template <nm::dtype_t Enum> struct dtype_enum_T;
  template <> struct dtype_enum_T<nm::BYTE> { typedef uint8_t type; };
  template <> struct dtype_enum_T<nm::INT8>  { typedef int8_t type;  };
  template <> struct dtype_enum_T<nm::INT16> { typedef int16_t type; };
  template <> struct dtype_enum_T<nm::INT32> { typedef int32_t type; };
  template <> struct dtype_enum_T<nm::INT64> { typedef int64_t type; };
  template <> struct dtype_enum_T<nm::FLOAT32> { typedef float type; };
  template <> struct dtype_enum_T<nm::FLOAT64> { typedef double type; };
  template <> struct dtype_enum_T<nm::COMPLEX64> { typedef nm::Complex64 type; };
  template <> struct dtype_enum_T<nm::COMPLEX128> { typedef nm::Complex128 type; };
  template <> struct dtype_enum_T<nm::RATIONAL32> { typedef nm::Rational32 type; };
  template <> struct dtype_enum_T<nm::RATIONAL64> { typedef nm::Rational64 type; };
  template <> struct dtype_enum_T<nm::RATIONAL128> { typedef nm::Rational128 type; };
  template <> struct dtype_enum_T<nm::RUBYOBJ> { typedef nm::RubyObject type; };

} // end namespace nm

#endif