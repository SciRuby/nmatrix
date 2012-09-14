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
// SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
// NMatrix is Copyright (c) 2012, Ruby Science Foundation
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
// == yale.h
//
// "new yale" storage format for 2D matrices (like yale, but with
// the diagonal pulled out for O(1) access).
//
// Specifications:
// * dtype and index dtype must necessarily differ
//      * index dtype is defined by whatever unsigned type can store
//        max(rows,cols)
//      * that means vector ija stores only index dtype, but a stores
//        dtype
// * vectors must be able to grow as necessary
//      * maximum size is rows*cols+1

#ifndef YALE_H
#define YALE_H

/*
 * Standard Includes
 */

#include <limits> // for std::numeric_limits<T>::max()

/*
 * Project Includes
 */

#include "types.h"

#include "data/data.h"

#include "common.h"

#include "nmatrix.h"

extern "C" {

  /*
   * Macros
   */

  #define NM_YALE_MINIMUM(sptr)               (((YALE_STORAGE*)(sptr))->shape[0] + 1) // arbitrarily defined

  #ifndef NM_CHECK_ALLOC
   #define NM_CHECK_ALLOC(x) if (!x) rb_raise(rb_eNoMemError, "insufficient memory");
  #endif

  /*
   * Types
   */


  /*
   * Data
   */


  /*
   * Functions
   */

  ///////////////
  // Lifecycle //
  ///////////////

  YALE_STORAGE* nm_yale_storage_create(dtype_t dtype, size_t* shape, size_t dim, size_t init_capacity);
  YALE_STORAGE* nm_yale_storage_create_from_old_yale(dtype_t dtype, size_t* shape, void* ia, void* ja, void* a, dtype_t from_dtype);
  YALE_STORAGE*	nm_yale_storage_create_merged(const YALE_STORAGE* merge_template, const YALE_STORAGE* other);
  void          nm_yale_storage_delete(STORAGE* s);
  void					nm_yale_storage_init(YALE_STORAGE* s);
  void					nm_yale_storage_mark(void*);

  ///////////////
  // Accessors //
  ///////////////

  void* nm_yale_storage_get(STORAGE* s, SLICE* slice);
  void*	nm_yale_storage_ref(STORAGE* s, SLICE* slice);
  char  nm_yale_storage_set(STORAGE* storage, SLICE* slice, void* v);

  inline size_t  nm_yale_storage_get_size(const YALE_STORAGE* storage);

  ///////////
  // Tests //
  ///////////

  bool nm_yale_storage_eqeq(const STORAGE* left, const STORAGE* right);

  //////////
  // Math //
  //////////
	
	STORAGE* nm_yale_storage_ew_op(nm::ewop_t op, const STORAGE* left, const STORAGE* right);
  STORAGE* nm_yale_storage_matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector);

  /////////////
  // Utility //
  /////////////

  /*
   * Calculates the itype a YALE_STORAGE object would need without actually needing
   * to see the YALE_STORAGE object. Does this just by looking at the shape.
   *
   * Useful for creating Yale Storage by other means than NMatrix.new(:yale, ...),
   * e.g., from a MATLAB v5 .mat file.
   */
  inline itype_t nm_yale_storage_itype_by_shape(const size_t* shape) {
    uint64_t yale_max_size = shape[0] * (shape[1]+1);

    if (yale_max_size < static_cast<uint64_t>(std::numeric_limits<uint8_t>::max()) - 2) {
      return UINT8;

    } else if (yale_max_size < static_cast<uint64_t>(std::numeric_limits<uint16_t>::max()) - 2) {
      return UINT16;

    } else if (yale_max_size < std::numeric_limits<uint32_t>::max() - 2) {
      return UINT32;

    } else {
      return UINT64;
    }
  }

  /*
   * Determine the index dtype (which will be used for the ija vector). This is
   * determined by matrix shape, not IJA/A vector capacity. Note that it's MAX-2
   * because UINTX_MAX and UINTX_MAX-1 are both reserved for sparse matrix
   * multiplication.
   */
  inline itype_t nm_yale_storage_itype(const YALE_STORAGE* s) {
    return nm_yale_storage_itype_by_shape(s->shape);
  }


  /////////////////////////
  // Copying and Casting //
  /////////////////////////

  STORAGE*      nm_yale_storage_cast_copy(const STORAGE* rhs, dtype_t new_dtype);
  STORAGE*      nm_yale_storage_copy_transposed(const STORAGE* rhs_base);



  void nm_init_yale_functions(void);


} // end of extern "C" block

namespace nm { namespace yale_storage {

  /*
   * Constants
   */
  const float GROWTH_CONSTANT = 1.5;


  /*
   * Templated Functions
   */

  template <typename IType>
  int binary_search(YALE_STORAGE* s, IType left, IType right, IType key);

  /*
   * Clear out the D portion of the A vector (clearing the diagonal and setting
   * the zero value).
   *
   * Note: This sets a literal 0 value. If your dtype is RUBYOBJ (a Ruby object),
   * it'll actually be INT2FIX(0) instead of a string of NULLs.
   */
  template <typename DType>
  inline void clear_diagonal_and_zero(YALE_STORAGE* s) {
    DType* a = reinterpret_cast<DType*>(s->a);

    // Clear out the diagonal + one extra entry
    for (size_t i = 0; i < s->shape[0]+1; ++i) // insert Ruby zeros
      a[i] = 0;
  }

  template <typename DType, typename IType>
  void init(YALE_STORAGE* s);

}} // end of namespace nm::yale_storage

#endif // YALE_H
