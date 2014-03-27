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
// == dense.h
//
// Dense n-dimensional matrix storage.

#ifndef DENSE_H
#define DENSE_H

/*
 * Standard Includes
 */

#include <stdlib.h>

/*
 * Project Includes
 */

#include "types.h"
//#include "util/math.h"

#include "data/data.h"

#include "../common.h"

#include "nmatrix.h"

/*
 * Macros
 */

/*
 * Types
 */

/*
 * Data
 */

extern "C" {

/*
 * Functions
 */

///////////////
// Lifecycle //
///////////////

DENSE_STORAGE*	nm_dense_storage_create(nm::dtype_t dtype, size_t* shape, size_t dim, void* elements, size_t elements_length);
void						nm_dense_storage_delete(STORAGE* s);
void						nm_dense_storage_delete_ref(STORAGE* s);
void						nm_dense_storage_mark(STORAGE*);
void            nm_dense_storage_register(const STORAGE* s);
void            nm_dense_storage_unregister(const STORAGE* s);


///////////////
// Accessors //
///////////////


VALUE nm_dense_map_pair(VALUE self, VALUE right);
VALUE nm_dense_map(VALUE self);
VALUE nm_dense_each(VALUE nmatrix);
VALUE nm_dense_each_with_indices(VALUE nmatrix);
void*	nm_dense_storage_get(const STORAGE* s, SLICE* slice);
void*	nm_dense_storage_ref(const STORAGE* s, SLICE* slice);
void  nm_dense_storage_set(VALUE left, SLICE* slice, VALUE right);

///////////
// Tests //
///////////

bool nm_dense_storage_eqeq(const STORAGE* left, const STORAGE* right);
bool nm_dense_storage_is_symmetric(const DENSE_STORAGE* mat, int lda);
bool nm_dense_storage_is_hermitian(const DENSE_STORAGE* mat, int lda);

//////////
// Math //
//////////

STORAGE* nm_dense_storage_matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector);

/////////////
// Utility //
/////////////

size_t nm_dense_storage_pos(const DENSE_STORAGE* s, const size_t* coords);
void nm_dense_storage_coords(const DENSE_STORAGE* s, const size_t slice_pos, size_t* coords_out);

/////////////////////////
// Copying and Casting //
/////////////////////////

DENSE_STORAGE*  nm_dense_storage_copy(const DENSE_STORAGE* rhs);
STORAGE*        nm_dense_storage_copy_transposed(const STORAGE* rhs_base);
STORAGE*        nm_dense_storage_cast_copy(const STORAGE* rhs, nm::dtype_t new_dtype, void*);

} // end of extern "C" block

namespace nm {
  std::pair<NMATRIX*,bool> interpret_arg_as_dense_nmatrix(VALUE right, nm::dtype_t dtype);
} // end of namespace nm

#endif // DENSE_H
