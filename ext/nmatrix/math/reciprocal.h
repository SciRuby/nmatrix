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
// == math/reciprocal.h
//
// Helper function for computing the reciprocal of a template type.
// Called by getrf and possibly other things.
//

#ifndef RECIPROCAL_H
#define RECIPROCAL_H


namespace nm { namespace math {

/* Numeric inverse -- basically just 1 / n, which in Ruby is 1.quo(n). */
template <typename DType>
inline DType reciprocal(const DType& n) {
  return n.reciprocal();
}
template <> inline float reciprocal(const float& n) { return 1 / n; }
template <> inline double reciprocal(const double& n) { return 1 / n; }

}}


#endif // RECIPROCAL_H

