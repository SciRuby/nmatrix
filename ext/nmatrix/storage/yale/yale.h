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
#include <stdexcept>

/*
 * Project Includes
 */

#include "../../types.h"
#include "../../data/data.h"
#include "../common.h"
#include "../../nmatrix.h"

extern "C" {

  /*
   * Macros
   */

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

  YALE_STORAGE* nm_yale_storage_create(nm::dtype_t dtype, size_t* shape, size_t dim, size_t init_capacity);
  YALE_STORAGE* nm_yale_storage_create_from_old_yale(nm::dtype_t dtype, size_t* shape, char* ia, char* ja, char* a, nm::dtype_t from_dtype);
  YALE_STORAGE*	nm_yale_storage_create_merged(const YALE_STORAGE* merge_template, const YALE_STORAGE* other);
  void          nm_yale_storage_delete(STORAGE* s);
  void          nm_yale_storage_delete_ref(STORAGE* s);
  void					nm_yale_storage_init(YALE_STORAGE* s, void* default_val);
  void					nm_yale_storage_mark(STORAGE*);
  void          nm_yale_storage_register(const STORAGE* s);
  void          nm_yale_storage_unregister(const STORAGE* s);
  void		nm_yale_storage_register_a(void* a, size_t size);
  void		nm_yale_storage_unregister_a(void* a, size_t size); 
    
  ///////////////
  // Accessors //
  ///////////////

  VALUE nm_yale_each_with_indices(VALUE nmatrix);
  VALUE nm_yale_each_stored_with_indices(VALUE nmatrix);
  VALUE nm_yale_stored_diagonal_each_with_indices(VALUE nmatrix);
  VALUE nm_yale_stored_nondiagonal_each_with_indices(VALUE nmatrix);
  VALUE nm_yale_each_ordered_stored_with_indices(VALUE nmatrix);
  void* nm_yale_storage_get(const STORAGE* s, SLICE* slice);
  void*	nm_yale_storage_ref(const STORAGE* s, SLICE* slice);
  void  nm_yale_storage_set(VALUE left, SLICE* slice, VALUE right);

  //char  nm_yale_storage_vector_insert(YALE_STORAGE* s, size_t pos, size_t* js, void* vals, size_t n, bool struct_only, nm::dtype_t dtype, nm::itype_t itype);
  //void  nm_yale_storage_increment_ia_after(YALE_STORAGE* s, size_t ija_size, size_t i, size_t n);

  size_t  nm_yale_storage_get_size(const YALE_STORAGE* storage);
  VALUE   nm_yale_default_value(VALUE self);
  VALUE   nm_yale_map_stored(VALUE self);
  VALUE   nm_yale_map_merged_stored(VALUE left, VALUE right, VALUE init);

  ///////////
  // Tests //
  ///////////

  bool nm_yale_storage_eqeq(const STORAGE* left, const STORAGE* right);

  //////////
  // Math //
  //////////

  STORAGE* nm_yale_storage_matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector);

  /////////////
  // Utility //
  /////////////



  /////////////////////////
  // Copying and Casting //
  /////////////////////////

  STORAGE*      nm_yale_storage_cast_copy(const STORAGE* rhs, nm::dtype_t new_dtype, void*);
  STORAGE*      nm_yale_storage_copy_transposed(const STORAGE* rhs_base);



  void nm_init_yale_functions(void);

  VALUE nm_vector_set(int argc, VALUE* argv, VALUE self);


} // end of extern "C" block

namespace nm {

namespace yale_storage {

  /*
   * Typedefs
   */

  typedef size_t IType;


  /*
   * Templated Functions
   */

  int binary_search(YALE_STORAGE* s, IType left, IType right, IType key);

  /*
   * Clear out the D portion of the A vector (clearing the diagonal and setting
   * the zero value).
   *
   * Note: This sets a literal 0 value. If your dtype is RUBYOBJ (a Ruby object),
   * it'll actually be INT2FIX(0) instead of a string of NULLs. You can actually
   * set a default for Ruby objects other than zero -- you generally want it to
   * be Qfalse, Qnil, or INT2FIX(0). The last is the default.
   */
  template <typename DType>
  inline void clear_diagonal_and_zero(YALE_STORAGE* s, void* init_val) {
    DType* a = reinterpret_cast<DType*>(s->a);

    // Clear out the diagonal + one extra entry
    if (init_val) {
      for (size_t i = 0; i <= s->shape[0]; ++i) // insert Ruby zeros, falses, or whatever else.
        a[i] = *reinterpret_cast<DType*>(init_val);
    } else {
      for (size_t i = 0; i <= s->shape[0]; ++i) // insert zeros.
        a[i] = 0;
    }
  }

  template <typename DType>
  void init(YALE_STORAGE* s, void* init_val);

  size_t  get_size(const YALE_STORAGE* storage);

  IType binary_search_left_boundary(const YALE_STORAGE* s, IType left, IType right, IType bound);


}} // end of namespace nm::yale_storage

#endif // YALE_H
