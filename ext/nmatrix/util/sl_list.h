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
// == sl_list.h
//
// Singly-linked list implementation used for List Storage.

#ifndef SL_LIST_H
#define SL_LIST_H


/*
 * Standard Includes
 */

#include <type_traits>
#include <cstdlib>

/*
 * Project Includes
 */

#include "types.h"

#include "data/data.h"

#include "nmatrix.h"

namespace nm { namespace list {

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
 
////////////////
// Lifecycle //
///////////////

LIST*	create(void);
void	del(LIST* list, size_t recursions);
void	mark(LIST* list, size_t recursions);

///////////////
// Accessors //
///////////////

NODE* insert(LIST* list, bool replace, size_t key, void* val);
NODE* insert_copy(LIST *list, bool replace, size_t key, void *val, size_t size);
NODE* insert_first_node(LIST* list, size_t key, void* val, size_t val_size);
NODE* insert_first_list(LIST* list, size_t key, LIST* l);
NODE* insert_after(NODE* node, size_t key, void* val);
NODE* replace_insert_after(NODE* node, size_t key, void* val, bool copy, size_t copy_size);
void* remove(LIST* list, size_t key);
void* remove_by_node(LIST* list, NODE* prev, NODE* rm);
bool remove_recursive(LIST* list, const size_t* coords, const size_t* offset, const size_t* lengths, size_t r, const size_t& dim);
bool node_is_within_slice(NODE* n, size_t coord, size_t len);

template <typename Type>
inline NODE* insert_helper(LIST* list, NODE* node, size_t key, Type val) {
	Type* val_mem = NM_ALLOC(Type);
	*val_mem = val;
	
	if (node == NULL) {
		return insert(list, false, key, val_mem);
		
	} else {
		return insert_after(node, key, val_mem);
	}
}

template <typename Type>
inline NODE* insert_helper(LIST* list, NODE* node, size_t key, Type* ptr) {
	if (node == NULL) {
		return insert(list, false, key, ptr);
		
	} else {
		return insert_after(node, key, ptr);
	}
}

///////////
// Tests //
///////////


/////////////
// Utility //
/////////////

NODE* find(LIST* list, size_t key);
NODE* find_preceding_from_node(NODE* prev, size_t key);
NODE* find_preceding_from_list(LIST* l, size_t key);
NODE* find_nearest(LIST* list, size_t key);
NODE* find_nearest_from(NODE* prev, size_t key);

/////////////////////////
// Copying and Casting //
/////////////////////////

template <typename LDType, typename RDType>
void cast_copy_contents(LIST* lhs, const LIST* rhs, size_t recursions);

}} // end of namespace nm::list

extern "C" {
  void nm_list_cast_copy_contents(LIST* lhs, const LIST* rhs, nm::dtype_t lhs_dtype, nm::dtype_t rhs_dtype, size_t recursions);
  VALUE nm_list_copy_to_hash(const LIST* l, const nm::dtype_t dtype, size_t recursions, VALUE default_value);
} // end of extern "C" block

#endif // SL_LIST_H
