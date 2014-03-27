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
// Header file for utility functions and data.

#ifndef UTIL_H
#define UTIL_H

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

#include "types.h"

/*
 * Macros
 */

/*
 * Types
 */

/*
 * Data
 */

/*
 * Functions
 */
namespace nm {
  template <typename Type>
  inline Type gcf(Type x, Type y) {
    Type t;

    if (x < 0) x = -x;
    if (y < 0) y = -y;

    if (x == 0) return y;
    if (y == 0) return x;

    while (x > 0) {
      t = x;
      x = y % x;
      y = t;
    }

    return y;
  }
} // end of namespace nm


#endif // UTIL_H
