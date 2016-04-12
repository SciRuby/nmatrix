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
// == math/magnitude.h
//
// Takes the absolute value (meaning magnitude) of each DType.
// Needed for a variety of BLAS/LAPACK functions.
//

#ifndef MAGNITUDE_H
#define MAGNITUDE_H

#include "math/long_dtype.h"

namespace nm { namespace math {

/* Magnitude -- may be complicated for unsigned types, and need to call the correct STL abs for floats/doubles */ 
template <typename DType, typename MDType = typename MagnitudeDType<DType>::type>
inline MDType magnitude(const DType& v) {
  return v.abs();
}
template <> inline float magnitude(const float& v) { return std::abs(v); }
template <> inline double magnitude(const double& v) { return std::abs(v); }
template <> inline uint8_t magnitude(const uint8_t& v) { return v; }
template <> inline int8_t magnitude(const int8_t& v) { return std::abs(v); }
template <> inline int16_t magnitude(const int16_t& v) { return std::abs(v); }
template <> inline int32_t magnitude(const int32_t& v) { return std::abs(v); }
template <> inline int64_t magnitude(const int64_t& v) { return std::abs(v); }
template <> inline float magnitude(const nm::Complex64& v) { return std::sqrt(v.r * v.r + v.i * v.i); }
template <> inline double magnitude(const nm::Complex128& v) { return std::sqrt(v.r * v.r + v.i * v.i); } 
    
}}

#endif // MAGNITUDE_H
