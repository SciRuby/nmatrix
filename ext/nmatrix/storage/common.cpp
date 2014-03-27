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
// == common.cpp
//
// Code for the STORAGE struct that is common to all storage types.

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

#include "common.h"

/*
 * Macros
 */

/*
 * Global Variables
 */

/*
 * Forward Declarations
 */

/*
 * Functions
 */

extern "C" {
  /*
   * Calculate the number of elements in the dense storage structure, based on
   * shape and dim.
   */
  size_t nm_storage_count_max_elements(const STORAGE* storage) {
    unsigned int i;
    size_t count = 1;

    for (i = storage->dim; i-- > 0;) {
      count *= storage->shape[i];
    }

    return count;
  }

  // Helper function used only for the RETURN_SIZED_ENUMERATOR macro. Returns the length of
  // the matrix's storage.
  VALUE nm_enumerator_length(VALUE nmatrix) {
    long len = nm_storage_count_max_elements(NM_STORAGE_DENSE(nmatrix));
    return LONG2NUM(len);
  }

} // end of extern "C" block
