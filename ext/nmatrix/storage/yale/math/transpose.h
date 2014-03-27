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
// == transpose.h
//
// Functions for Yale math: transposing
//

#ifndef YALE_MATH_TRANSPOSE_H
# define YALE_MATH_TRANSPOSE_H

namespace nm { namespace yale_storage {

/*
 * Transposes a generic Yale matrix (old or new). Specify new by setting RDiag = true.
 *
 * Based on transp from SMMP (same as symbmm and numbmm).
 *
 * This is not named in the same way as most yale_storage functions because it does not act on a YALE_STORAGE
 * object.
 */

template <typename AD, typename BD, bool DiagA, bool Move>
void transpose_yale(const size_t n, const size_t m,
                    const size_t* ia, const size_t* ja, const AD* a, const AD& a_default,
                    size_t* ib, size_t* jb, BD* b, const BD& b_default) {

  size_t index;

  // Clear B
  for (size_t i = 0; i < m+1; ++i) ib[i] = 0;

  if (Move)
    for (size_t i = 0; i < m+1; ++i) b[i] = b_default;

  if (DiagA) ib[0] = m + 1;
  else       ib[0] = 0;

  /* count indices for each column */

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = ia[i]; j < ia[i+1]; ++j) {
      ++(ib[ja[j]+1]);
    }
  }

  for (size_t i = 0; i < m; ++i) {
    ib[i+1] = ib[i] + ib[i+1];
  }

  /* now make jb */

  for (size_t i = 0; i < n; ++i) {

    for (size_t j = ia[i]; j < ia[i+1]; ++j) {
      index = ja[j];
      jb[ib[index]] = i;

      if (Move && a[j] != a_default)
        b[ib[index]] = a[j];

      ++(ib[index]);
    }
  }

  /* now fixup ib */

  for (size_t i = m; i >= 1; --i) {
    ib[i] = ib[i-1];
  }


  if (DiagA) {
    if (Move) {
      size_t j = std::min(n,m);

      for (size_t i = 0; i < j; ++i) {
        b[i] = a[i];
      }
    }
    ib[0] = m + 1;

  } else {
    ib[0] = 0;
  }
}

} } // end of namespace nm::yale_storage

#endif