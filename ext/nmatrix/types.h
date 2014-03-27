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
// == types.h
//
// Definition of simple types used throughout NMatrix.

#ifndef NMATRIX_TYPES_H
#define NMATRIX_TYPES_H

/*
 * Standard Includes
 */

#include <cstdint>

/*
 * Macros
 */

#define EPSILON 1E-10
#define FP_IS_ZERO(n) (-EPSILON < n && n < EPSILON)
#define FP_EQUAL(a, b) FP_IS_ZERO((a - b))

/*
 * Types
 */

typedef float		float32_t;
typedef double	float64_t;

typedef size_t  IType;

#endif
