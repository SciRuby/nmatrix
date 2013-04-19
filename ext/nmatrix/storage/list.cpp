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
// == list.c
//
// List-of-lists n-dimensional matrix storage. Uses singly-linked
// lists.

/*
 * Standard Includes
 */

#include <ruby.h>
#include <algorithm> // std::min
#include <iostream>

/*
 * Project Includes
 */

#include "types.h"

#include "data/data.h"

#include "common.h"
#include "list.h"

#include "util/math.h"
#include "util/sl_list.h"

/*
 * Macros
 */

/*
 * Global Variables
 */

namespace nm { namespace list_storage {

/*
 * Forward Declarations
 */

template <typename LDType, typename RDType>
static LIST_STORAGE* cast_copy(const LIST_STORAGE* rhs, dtype_t new_dtype);

template <typename LDType, typename RDType>
static bool eqeq(const LIST_STORAGE* left, const LIST_STORAGE* right);

template <ewop_t op, typename LDType, typename RDType>
static void* ew_op(LIST* dest, const LIST* left, const void* l_default, const LIST* right, const void* r_default, const size_t* shape, size_t dim);

template <ewop_t op, typename LDType, typename RDType>
static void ew_op_prime(LIST* dest, LDType d_default, const LIST* left, LDType l_default, const LIST* right, RDType r_default, const size_t* shape, size_t last_level, size_t level);

template <ewop_t op, typename LDType, typename RDType>
static void ew_comp_prime(LIST* dest, uint8_t d_default, const LIST* left, LDType l_default, const LIST* right, RDType r_default, const size_t* shape, size_t last_level, size_t level);

} // end of namespace list_storage

extern "C" {

/*
 * Functions
 */


////////////////
// Lifecycle //
///////////////

/*
 * Creates a list-of-lists(-of-lists-of-lists-etc) storage framework for a
 * matrix.
 *
 * Note: The pointers you pass in for shape and init_val become property of our
 * new storage. You don't need to free them, and you shouldn't re-use them.
 */
LIST_STORAGE* nm_list_storage_create(dtype_t dtype, size_t* shape, size_t dim, void* init_val) {
  LIST_STORAGE* s;

  s = ALLOC( LIST_STORAGE );

  s->dim   = dim;
  s->shape = shape;
  s->dtype = dtype;

  s->offset = ALLOC_N(size_t, s->dim);
  memset(s->offset, 0, s->dim * sizeof(size_t));

  s->rows  = list::create();
  s->default_val = init_val;
  s->count = 1;
  s->src = s;

  return s;
}

/*
 * Documentation goes here.
 */
void nm_list_storage_delete(STORAGE* s) {
  if (s) {
    LIST_STORAGE* storage = (LIST_STORAGE*)s;
    if (storage->count-- == 1) {
      list::del( storage->rows, storage->dim - 1 );

      free(storage->shape);
      free(storage->offset);
      free(storage->default_val);
      free(s);
    }
  }
}

/*
 * Documentation goes here.
 */
void nm_list_storage_delete_ref(STORAGE* s) {
  if (s) {
    LIST_STORAGE* storage = (LIST_STORAGE*)s;

    nm_list_storage_delete( reinterpret_cast<STORAGE*>(storage->src ) );
    free(storage->shape);
    free(storage->offset);
    free(s);
  }
}

/*
 * Documentation goes here.
 */
void nm_list_storage_mark(void* storage_base) {
  LIST_STORAGE* storage = (LIST_STORAGE*)storage_base;

  if (storage && storage->dtype == RUBYOBJ) {
    rb_gc_mark(*((VALUE*)(storage->default_val)));
    list::mark(storage->rows, storage->dim - 1);
  }
}

///////////////
// Accessors //
///////////////

/*
 * Documentation goes here.
 */
NODE* list_storage_get_single_node(LIST_STORAGE* s, SLICE* slice)
{
  size_t r;
  LIST*  l = s->rows;
  NODE*  n;

  for (r = 0; r < s->dim; r++) {
    n = list::find(l, s->offset[r] + slice->coords[r]);
    if (n)  l = reinterpret_cast<LIST*>(n->val);
    else return NULL;
  }

  return n;
}
 /*
   * Each stored iterator, brings along the indices
   */
  VALUE nm_list_each_stored_with_indices(VALUE nmatrix) {

    // If we don't have a block, return an enumerator.
    RETURN_ENUMERATOR(nmatrix, 0, 0);
    LIST_STORAGE* s = NM_STORAGE_LIST(nmatrix);
    // Create indices and initialize to zero
    size_t* coords = ALLOCA_N(size_t, s->dim);
    memset(coords, 0, sizeof(size_t) * s->dim);

    // Set up the LIST and NODE
    LIST* l = s->rows;
    NODE* curr = l->first;

    // Levels...
    while (curr) { // LOOPS through the rows
      NODE* subnode = reinterpret_cast<LIST*>(curr->val)->first;
      size_t row = curr->key;
      // Iterate along each row, yielding the val, row, col for each non-default entry
      while (subnode) {
        VALUE ary = rb_ary_new();
        size_t col = subnode->key;
        
        // Conditional type handling
        if (NM_DTYPE(nmatrix) == nm::RUBYOBJ) {
          rb_ary_push(ary, *reinterpret_cast<VALUE*>(subnode->val));
        } else {
          rb_ary_push(ary, rubyobj_from_cval((char*)(subnode->val ) , NM_DTYPE(nmatrix)).rval);
        }

        // Push the coordinate values into the ary
        rb_ary_push(ary, INT2FIX(row));
        rb_ary_push(ary, INT2FIX(col));

        // Yield the ary
        rb_yield(ary);

        // Update the col position
        subnode = subnode->next;
      }
      // Update the row node 
      curr = curr->next;
    }
    return nmatrix;
  }


static LIST* slice_copy(const LIST_STORAGE *src, LIST *src_rows, size_t *coords, size_t *lengths, size_t n)
{
  NODE *src_node;
  LIST *dst_rows = NULL;
  void *val = NULL;
  int key;
  
  dst_rows = list::create();
  src_node = src_rows->first;

  while (src_node) {
    key = src_node->key - (src->offset[n] + coords[n]);
    
    if (key >= 0 && (size_t)key < lengths[n]) {
      if (src->dim - n > 1) {
        val = slice_copy(src,  
          reinterpret_cast<LIST*>(src_node->val), 
          coords,
          lengths,
          n + 1);  

        if (val) 
          list::insert_with_copy(dst_rows, key, val, sizeof(LIST));          
        
      }
      else {
        list::insert_with_copy(dst_rows, key, src_node->val, DTYPE_SIZES[src->dtype]);
      }
    }

    src_node = src_node->next;
  }

  return dst_rows;
}

/*
 * Documentation goes here.
 */
void* nm_list_storage_get(STORAGE* storage, SLICE* slice) {
  LIST_STORAGE* s = (LIST_STORAGE*)storage;
  LIST_STORAGE* ns = NULL;
  NODE* n;

  if (slice->single) {
    n = list_storage_get_single_node(s, slice); 
    return (n ? n->val : s->default_val);
  } 
  else {
    void *init_val = ALLOC_N(char, DTYPE_SIZES[s->dtype]);
    memcpy(init_val, s->default_val, DTYPE_SIZES[s->dtype]);

    size_t *shape = ALLOC_N(size_t, s->dim);
    memcpy(shape, slice->lengths, sizeof(size_t) * s->dim);

    ns = nm_list_storage_create(s->dtype, shape, s->dim, init_val);
    
    ns->rows = slice_copy(s, s->rows, slice->coords, slice->lengths, 0);
    return ns;
  }
}

/*
 * Get the contents of some set of coordinates. Note: Does not make a copy!
 * Don't free!
 */
void* nm_list_storage_ref(STORAGE* storage, SLICE* slice) {
  LIST_STORAGE* s = (LIST_STORAGE*)storage;
  LIST_STORAGE* ns = NULL;
  NODE* n;

  //TODO: It needs a refactoring.
  if (slice->single) {
    n = list_storage_get_single_node(s, slice); 
    return (n ? n->val : s->default_val);
  } 
  else {
    ns = ALLOC( LIST_STORAGE );
    
    ns->dim = s->dim;
    ns->dtype = s->dtype;
    ns->offset     = ALLOC_N(size_t, ns->dim);
    ns->shape      = ALLOC_N(size_t, ns->dim);

    for (size_t i = 0; i < ns->dim; ++i) {
      ns->offset[i] = slice->coords[i] + s->offset[i];
      ns->shape[i]  = slice->lengths[i];
    }

    ns->rows = s->rows;
    ns->default_val = s->default_val;
    
    s->src->count++;
    ns->src = s->src;
    
    return ns;
  }
}

/*
 * Documentation goes here.
 *
 * TODO: Allow this function to accept an entire row and not just one value -- for slicing
 */
void* nm_list_storage_insert(STORAGE* storage, SLICE* slice, void* val) {
  LIST_STORAGE* s = (LIST_STORAGE*)storage;
  // Pretend dims = 2
  // Then coords is going to be size 2
  // So we need to find out if some key already exists
  size_t r;
  NODE*  n;
  LIST*  l = s->rows;

  // drill down into the structure
  for (r = s->dim; r > 1; --r) {
    n = list::insert(l, false, s->offset[s->dim - r] + slice->coords[s->dim - r], list::create());
    l = reinterpret_cast<LIST*>(n->val);
  }

  n = list::insert(l, true, s->offset[s->dim - r] + slice->coords[s->dim - r], val);
  return n->val;
}

/*
 * Remove an item from list storage.
 */
void* nm_list_storage_remove(STORAGE* storage, SLICE* slice) {
  LIST_STORAGE* s = (LIST_STORAGE*)storage;
  void* rm = NULL;

  // This returns a boolean, which will indicate whether s->rows is empty.
  // We can safely ignore it, since we never want to delete s->rows until
  // it's time to destroy the LIST_STORAGE object.
  list::remove_recursive(s->rows, slice->coords, s->offset, 0, s->dim, rm);

  return rm;
}

///////////
// Tests //
///////////

/*
 * Comparison of contents for list storage.
 */
bool nm_list_storage_eqeq(const STORAGE* left, const STORAGE* right) {
	NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::list_storage::eqeq, bool, const LIST_STORAGE* left, const LIST_STORAGE* right);

	return ttable[left->dtype][right->dtype]((const LIST_STORAGE*)left, (const LIST_STORAGE*)right);
}

//////////
// Math //
//////////

/*
 * Element-wise operations for list storage.
 *
 * If a scalar is given, a temporary matrix is created with that scalar as a default value.
 */
STORAGE* nm_list_storage_ew_op(nm::ewop_t op, const STORAGE* left, const STORAGE* right, VALUE scalar) {
  // rb_raise(rb_eNotImpError, "elementwise operations for list storage currently broken");

  bool cleanup = false;
  LIST_STORAGE *r, *new_l;
  const LIST_STORAGE* l = reinterpret_cast<const LIST_STORAGE*>(left);

  if (!right) { // need to build a right-hand matrix temporarily, with default value of 'scalar'

    dtype_t scalar_dtype  = nm_dtype_guess(scalar);
    void* scalar_init     = rubyobj_to_cval(scalar, scalar_dtype);

    size_t* shape         = ALLOC_N(size_t, l->dim);
    memcpy(shape, left->shape, sizeof(size_t) * l->dim);

    r = nm_list_storage_create(scalar_dtype, shape, l->dim, scalar_init);

    cleanup = true;

  } else {

    r = reinterpret_cast<LIST_STORAGE*>(const_cast<STORAGE*>(right));

  }

  // We may need to upcast our arguments to the same type.
  dtype_t new_dtype = Upcast[left->dtype][r->dtype];

	// Make sure we allocate a byte-storing matrix for comparison operations; otherwise, use the argument dtype (new_dtype)
	dtype_t result_dtype = static_cast<uint8_t>(op) < NUM_NONCOMP_EWOPS ? new_dtype : BYTE;

	OP_LR_DTYPE_TEMPLATE_TABLE(nm::list_storage::ew_op, void*, LIST* dest, const LIST* left, const void* l_default, const LIST* right, const void* r_default, const size_t* shape, size_t dim);
	
	// Allocate a new shape array for the resulting matrix.
	size_t* new_shape = ALLOC_N(size_t, l->dim);
	memcpy(new_shape, left->shape, sizeof(size_t) * l->dim);
	
	// Create the result matrix.
	LIST_STORAGE* result = nm_list_storage_create(result_dtype, new_shape, left->dim, NULL);
	
	/*
	 * Call the templated elementwise multiplication function and set the default
	 * value for the resulting matrix.
	 */
	if (new_dtype != left->dtype) {
		// Upcast the left-hand side if necessary.
		new_l = reinterpret_cast<LIST_STORAGE*>(nm_list_storage_cast_copy(l, new_dtype));
		
		result->default_val =
			ttable[op][new_l->dtype][r->dtype](result->rows, new_l->rows, new_l->default_val, r->rows, r->default_val, result->shape, result->dim);
		
		// Delete the temporary left-hand side matrix.
		nm_list_storage_delete(reinterpret_cast<STORAGE*>(new_l));

	} else {
		result->default_val =
			ttable[op][left->dtype][r->dtype](result->rows, l->rows, l->default_val, r->rows, r->default_val, result->shape, result->dim);
	}

  // If we created a temporary scalar matrix (for matrix-scalar operations), we now need to delete it.
	if (cleanup) {
	  nm_list_storage_delete(reinterpret_cast<STORAGE*>(r));
	}
	
	return result;
}


/*
 * List storage matrix multiplication.
 */
STORAGE* nm_list_storage_matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector) {
  free(resulting_shape);
  rb_raise(rb_eNotImpError, "multiplication not implemented for list-of-list matrices");
  return NULL;
  //DTYPE_TEMPLATE_TABLE(dense_storage::matrix_multiply, NMATRIX*, STORAGE_PAIR, size_t*, bool);

  //return ttable[reinterpret_cast<DENSE_STORAGE*>(casted_storage.left)->dtype](casted_storage, resulting_shape, vector);
}


/*
 * List storage to Hash conversion. Uses Hashes with default values, so you can continue to pretend
 * it's a sparse matrix.
 */
VALUE nm_list_storage_to_hash(const LIST_STORAGE* s, const dtype_t dtype) {

  // Get the default value for the list storage.
  VALUE default_value = rubyobj_from_cval(s->default_val, dtype).rval;

  // Recursively copy each dimension of the matrix into a nested hash.
  return nm_list_copy_to_hash(s->rows, dtype, s->dim - 1, default_value);
}

/////////////
// Utility //
/////////////

/*
 * Recursively count the non-zero elements in a list storage object.
 */
size_t nm_list_storage_count_elements_r(const LIST* l, size_t recursions) {
  size_t count = 0;
  NODE* curr = l->first;
  
  if (recursions) {
    while (curr) {
      count += nm_list_storage_count_elements_r(reinterpret_cast<const LIST*>(curr->val), recursions - 1);
      curr   = curr->next;
    }
    
  } else {
    while (curr) {
      ++count;
      curr = curr->next;
    }
  }
  
  return count;
}

/*
 * Count non-diagonal non-zero elements.
 */
size_t nm_list_storage_count_nd_elements(const LIST_STORAGE* s) {
  NODE *i_curr, *j_curr;
  size_t count = 0;
  
  if (s->dim != 2) {
  	rb_raise(rb_eNotImpError, "non-diagonal element counting only defined for dim = 2");
  }

  for (i_curr = s->rows->first; i_curr; i_curr = i_curr->next) {
    int i = i_curr->key - s->offset[0];
    if (i < 0 || i >= (int)s->shape[0]) continue;

    for (j_curr = ((LIST*)(i_curr->val))->first; j_curr; j_curr = j_curr->next) {
      int j = j_curr->key - s->offset[1];
      if (j < 0 || j >= (int)s->shape[1]) continue;

      if (i != j)  	++count;
    }
  }
  
  return count;
}

/////////////////////////
// Copying and Casting //
/////////////////////////
//
/*
 * List storage copy constructor C access.
 */

LIST_STORAGE* nm_list_storage_copy(const LIST_STORAGE* rhs)
{
  size_t *shape = ALLOC_N(size_t, rhs->dim);
  memcpy(shape, rhs->shape, sizeof(size_t) * rhs->dim);
  
  void *init_val = ALLOC_N(char, DTYPE_SIZES[rhs->dtype]);
  memcpy(init_val, rhs->default_val, DTYPE_SIZES[rhs->dtype]);

  LIST_STORAGE* lhs = nm_list_storage_create(rhs->dtype, shape, rhs->dim, init_val);
  
  lhs->rows = slice_copy(rhs, rhs->rows, lhs->offset, lhs->shape, 0);

  return lhs;
}

/*
 * List storage copy constructor C access with casting.
 */
STORAGE* nm_list_storage_cast_copy(const STORAGE* rhs, dtype_t new_dtype) {
  NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::list_storage::cast_copy, LIST_STORAGE*, const LIST_STORAGE* rhs, dtype_t new_dtype);

  return (STORAGE*)ttable[new_dtype][rhs->dtype]((LIST_STORAGE*)rhs, new_dtype);
}


/*
 * List storage copy constructor for transposing.
 */
STORAGE* nm_list_storage_copy_transposed(const STORAGE* rhs_base) {
  rb_raise(rb_eNotImpError, "list storage transpose not yet implemented");
  return NULL;
}


} // end of extern "C" block


/////////////////////////
// Templated Functions //
/////////////////////////

namespace list_storage {

/*
 * List storage copy constructor for changing dtypes.
 */
template <typename LDType, typename RDType>
static LIST_STORAGE* cast_copy(const LIST_STORAGE* rhs, dtype_t new_dtype) {

  // allocate and copy shape
  size_t* shape = ALLOC_N(size_t, rhs->dim);
  memcpy(shape, rhs->shape, rhs->dim * sizeof(size_t));

  // copy default value
  LDType* default_val = ALLOC_N(LDType, 1);
  *default_val = *reinterpret_cast<RDType*>(rhs->default_val);

  LIST_STORAGE* lhs = nm_list_storage_create(new_dtype, shape, rhs->dim, default_val);
  lhs->rows         = list::create();

  // TODO: Needs optimization. When matrix is reference it is copped twice.
  if (rhs->src == rhs) 
    list::cast_copy_contents<LDType, RDType>(lhs->rows, rhs->rows, rhs->dim - 1);
  else {
    LIST_STORAGE *tmp = nm_list_storage_copy(rhs);
    list::cast_copy_contents<LDType, RDType>(lhs->rows, tmp->rows, rhs->dim - 1);
    nm_list_storage_delete(tmp);
  }

  return lhs;
}


/*
 * Do these two dense matrices of the same dtype have exactly the same
 * contents?
 */
template <typename LDType, typename RDType>
bool eqeq(const LIST_STORAGE* left, const LIST_STORAGE* right) {
  bool result;

  // in certain cases, we need to keep track of the number of elements checked.
  size_t num_checked  = 0,
	max_elements = nm_storage_count_max_elements(left);
  LIST_STORAGE *tmp1 = NULL, *tmp2 = NULL;

  if (!left->rows->first) {
    // Easy: both lists empty -- just compare default values
    if (!right->rows->first) {
    	return *reinterpret_cast<LDType*>(left->default_val) == *reinterpret_cast<RDType*>(right->default_val);
    	
    } else if (!list::eqeq_value<RDType,LDType>(right->rows, reinterpret_cast<LDType*>(left->default_val), left->dim-1, num_checked)) {
    	// Left empty, right not empty. Do all values in right == left->default_val?
    	return false;
    	
    } else if (num_checked < max_elements) {
    	// If the matrix isn't full, we also need to compare default values.
    	return *reinterpret_cast<LDType*>(left->default_val) == *reinterpret_cast<RDType*>(right->default_val);
    }

  } else if (!right->rows->first) {
    // fprintf(stderr, "!right->rows true\n");
    // Right empty, left not empty. Do all values in left == right->default_val?
    if (!list::eqeq_value<LDType,RDType>(left->rows, reinterpret_cast<RDType*>(right->default_val), left->dim-1, num_checked)) {
    	return false;
    	
    } else if (num_checked < max_elements) {
   		// If the matrix isn't full, we also need to compare default values.
    	return *reinterpret_cast<LDType*>(left->default_val) == *reinterpret_cast<RDType*>(right->default_val);
    }

  } else {
    // fprintf(stderr, "both matrices have entries\n");
    // Hardest case. Compare lists node by node. Let's make it simpler by requiring that both have the same default value
    
    // left is reference
    if (left->src != left && right->src == right) {
      tmp1 = nm_list_storage_copy(left);
      result = list::eqeq<LDType,RDType>(tmp1->rows, right->rows, reinterpret_cast<LDType*>(left->default_val), reinterpret_cast<RDType*>(right->default_val), left->dim-1, num_checked);
      nm_list_storage_delete(tmp1);
    } 
    // right is reference
    if (left->src == left && right->src != right) {
      tmp2 = nm_list_storage_copy(right);
      result = list::eqeq<LDType,RDType>(left->rows, tmp2->rows, reinterpret_cast<LDType*>(left->default_val), reinterpret_cast<RDType*>(right->default_val), left->dim-1, num_checked);
      nm_list_storage_delete(tmp2);
    } 
    // both are references
    if (left->src != left && right->src != right) {
      tmp1 = nm_list_storage_copy(left);
      tmp2 = nm_list_storage_copy(right);
      result = list::eqeq<LDType,RDType>(tmp1->rows, tmp2->rows, reinterpret_cast<LDType*>(left->default_val), reinterpret_cast<RDType*>(right->default_val), left->dim-1, num_checked);
      nm_list_storage_delete(tmp1);
      nm_list_storage_delete(tmp2);
    }
    // both are normal matricies
    if (left->src == left && right->src == right) 
      result = list::eqeq<LDType,RDType>(left->rows, right->rows, reinterpret_cast<LDType*>(left->default_val), reinterpret_cast<RDType*>(right->default_val), left->dim-1, num_checked);

    if (!result){
      return result;
    } else if (num_checked < max_elements) {
      return *reinterpret_cast<LDType*>(left->default_val) == *reinterpret_cast<RDType*>(right->default_val);
    }
  }

  return true;
}

/*
 * List storage element-wise operations (including comparisons).
 */
template <ewop_t op, typename LDType, typename RDType>
static void* ew_op(LIST* dest, const LIST* left, const void* l_default, const LIST* right, const void* r_default, const size_t* shape, size_t dim) {

	if (static_cast<uint8_t>(op) < NUM_NONCOMP_EWOPS) {

    /*
     * Allocate space for, and calculate, the default value for the destination
     * matrix.
     */
    LDType* d_default_mem = ALLOC(LDType);
    *d_default_mem = ew_op_switch<op, LDType, RDType>(*reinterpret_cast<const LDType*>(l_default), *reinterpret_cast<const RDType*>(r_default));

    // Now that setup is done call the actual elementwise operation function.
    ew_op_prime<op, LDType, RDType>(dest, *reinterpret_cast<const LDType*>(d_default_mem),
                                    left, *reinterpret_cast<const LDType*>(l_default),
                                    right, *reinterpret_cast<const RDType*>(r_default),
                                    shape, dim - 1, 0);

    // Return a pointer to the destination matrix's default value.
    return d_default_mem;

	} else { // Handle comparison operations in a similar manner.
    /*
     * Allocate a byte for default, and set default value to 0.
     */
    uint8_t* d_default_mem = ALLOC(uint8_t);
    *d_default_mem = 0;
    switch (op) {
      case EW_EQEQ:
        *d_default_mem = *reinterpret_cast<const LDType*>(l_default) == *reinterpret_cast<const RDType*>(r_default);
        break;

      case EW_NEQ:
        *d_default_mem = *reinterpret_cast<const LDType*>(l_default) != *reinterpret_cast<const RDType*>(r_default);
        break;

      case EW_LT:
        *d_default_mem = *reinterpret_cast<const LDType*>(l_default) < *reinterpret_cast<const RDType*>(r_default);
        break;

      case EW_GT:
        *d_default_mem = *reinterpret_cast<const LDType*>(l_default) > *reinterpret_cast<const RDType*>(r_default);
        break;

      case EW_LEQ:
        *d_default_mem = *reinterpret_cast<const LDType*>(l_default) <= *reinterpret_cast<const RDType*>(r_default);
        break;

      case EW_GEQ:
        *d_default_mem = *reinterpret_cast<const LDType*>(l_default) >= *reinterpret_cast<const RDType*>(r_default);
        break;

      default:
        rb_raise(rb_eStandardError, "this should not happen");
    }

    // Now that setup is done call the actual elementwise comparison function.
    ew_comp_prime<op, LDType, RDType>(dest, *reinterpret_cast<const uint8_t*>(d_default_mem),
                                      left, *reinterpret_cast<const LDType*>(l_default),
                                      right, *reinterpret_cast<const RDType*>(r_default),
                                      shape, dim - 1, 0);

    // Return a pointer to the destination matrix's default value.
    return d_default_mem;
	}
}


/*
 * List storage element-wise comparisons, recursive helper.
 */
template <ewop_t op, typename LDType, typename RDType>
static void ew_comp_prime(LIST* dest, uint8_t d_default, const LIST* left, LDType l_default, const LIST* right, RDType r_default, const size_t* shape, size_t last_level, size_t level) {

	static LIST EMPTY_LIST = {NULL};

	size_t index;

	uint8_t tmp_result;

	LIST* new_level = NULL;

	NODE* l_node		= left->first,
			* r_node		= right->first,
			* dest_node	= NULL;

	for (index = 0; index < shape[level]; ++index) {
		if (l_node == NULL and r_node == NULL) {
			/*
			 * Both source lists are now empty.  Because the default value of the
			 * destination is already set appropriately we can now return.
			 */

			return;

		} else {
			// At least one list still has entries.

			if (l_node == NULL and r_node->key == index) {
				/*
				 * One source list is empty, but the index has caught up to the key of
				 * the other list.
				 */

				if (level == last_level) {
					switch (op) {
						case EW_EQEQ:
							tmp_result = (uint8_t)(l_default == *reinterpret_cast<RDType*>(r_node->val));
							break;

						case EW_NEQ:
							tmp_result = (uint8_t)(l_default != *reinterpret_cast<RDType*>(r_node->val));
							break;

						case EW_LT:
							tmp_result = (uint8_t)(l_default < *reinterpret_cast<RDType*>(r_node->val));
							break;

						case EW_GT:
							tmp_result = (uint8_t)(l_default > *reinterpret_cast<RDType*>(r_node->val));
							break;

						case EW_LEQ:
							tmp_result = (uint8_t)(l_default <= *reinterpret_cast<RDType*>(r_node->val));
							break;

						case EW_GEQ:
							tmp_result = (uint8_t)(l_default >= *reinterpret_cast<RDType*>(r_node->val));
							break;

            default:
              rb_raise(rb_eStandardError, "This should not happen.");
					}

					if (tmp_result != d_default) {
						dest_node = nm::list::insert_helper(dest, dest_node, index, tmp_result);
					}

				} else {
					new_level = nm::list::create();
					dest_node = nm::list::insert_helper(dest, dest_node, index, new_level);

					ew_comp_prime<op, LDType, RDType>(new_level, d_default, &EMPTY_LIST, l_default,
                                            reinterpret_cast<LIST*>(r_node->val), r_default,
                                            shape, last_level, level + 1);
				}

				r_node = r_node->next;

			} else if (r_node == NULL and l_node->key == index) {
				/*
				 * One source list is empty, but the index has caught up to the key of
				 * the other list.
				 */

				if (level == last_level) {
					switch (op) {
						case EW_EQEQ:
							tmp_result = (uint8_t)(*reinterpret_cast<LDType*>(l_node->val) == r_default);
							break;

						case EW_NEQ:
							tmp_result = (uint8_t)(*reinterpret_cast<LDType*>(l_node->val) != r_default);
							break;

						case EW_LT:
							tmp_result = (uint8_t)(*reinterpret_cast<LDType*>(l_node->val) < r_default);
							break;

						case EW_GT:
							tmp_result = (uint8_t)(*reinterpret_cast<LDType*>(l_node->val) > r_default);
							break;

						case EW_LEQ:
							tmp_result = (uint8_t)(*reinterpret_cast<LDType*>(l_node->val) <= r_default);
							break;

						case EW_GEQ:
							tmp_result = (uint8_t)(*reinterpret_cast<LDType*>(l_node->val) >= r_default);
							break;

            default:
              rb_raise(rb_eStandardError, "this should not happen");
					}

					if (tmp_result != d_default) {
						dest_node = nm::list::insert_helper(dest, dest_node, index, tmp_result);
					}

				} else {
					new_level = nm::list::create();
					dest_node = nm::list::insert_helper(dest, dest_node, index, new_level);

					ew_comp_prime<op, LDType, RDType>(new_level, d_default,
                                            reinterpret_cast<LIST*>(l_node->val), l_default,
                                            &EMPTY_LIST, r_default,
                                            shape, last_level, level + 1);
				}

				l_node = l_node->next;

			} else if (l_node != NULL and r_node != NULL and index == std::min(l_node->key, r_node->key)) {
				/*
				 * Neither list is empty and our index has caught up to one of the
				 * source lists.
				 */

				if (l_node->key == r_node->key) {

					if (level == last_level) {
						switch (op) {
							case EW_EQEQ:
								tmp_result = (uint8_t)(*reinterpret_cast<LDType*>(l_node->val) == *reinterpret_cast<RDType*>(r_node->val));
								break;

							case EW_NEQ:
								tmp_result = (uint8_t)(*reinterpret_cast<LDType*>(l_node->val) != *reinterpret_cast<RDType*>(r_node->val));
								break;

							case EW_LT:
								tmp_result = (uint8_t)(*reinterpret_cast<LDType*>(l_node->val) < *reinterpret_cast<RDType*>(r_node->val));
								break;

							case EW_GT:
								tmp_result = (uint8_t)(*reinterpret_cast<LDType*>(l_node->val) > *reinterpret_cast<RDType*>(r_node->val));
								break;

							case EW_LEQ:
								tmp_result = (uint8_t)(*reinterpret_cast<LDType*>(l_node->val) <= *reinterpret_cast<RDType*>(r_node->val));
								break;

							case EW_GEQ:
								tmp_result = (uint8_t)(*reinterpret_cast<LDType*>(l_node->val) >= *reinterpret_cast<RDType*>(r_node->val));
								break;

              default:
                rb_raise(rb_eStandardError, "this should not happen");
						}

						if (tmp_result != d_default) {
							dest_node = nm::list::insert_helper(dest, dest_node, index, tmp_result);
						}

					} else {
						new_level = nm::list::create();
						dest_node = nm::list::insert_helper(dest, dest_node, index, new_level);

						ew_comp_prime<op, LDType, RDType>(new_level, d_default,
							reinterpret_cast<LIST*>(l_node->val), l_default,
							reinterpret_cast<LIST*>(r_node->val), r_default,
							shape, last_level, level + 1);
					}

					l_node = l_node->next;
					r_node = r_node->next;

				} else if (l_node->key < r_node->key) {
					// Advance the left node knowing that the default value is OK.

					l_node = l_node->next;

				} else /* if (l_node->key > r_node->key) */ {
					// Advance the right node knowing that the default value is OK.

					r_node = r_node->next;
				}

			} else {
				/*
				 * Our index needs to catch up but the default value is OK.  This
				 * conditional is here only for documentation and should be optimized
				 * out.
				 */
			}
		}
	}
}



/*
 * List storage element-wise operations, recursive helper.
 */
template <ewop_t op, typename LDType, typename RDType>
static void ew_op_prime(LIST* dest, LDType d_default, const LIST* left, LDType l_default, const LIST* right, RDType r_default, const size_t* shape, size_t last_level, size_t level) {
	
	static LIST EMPTY_LIST = {NULL};
	
	size_t index;
	
	LDType tmp_result;
	
	LIST* new_level = NULL;
	
	NODE* l_node		= left->first,
			* r_node		= right->first,
			* dest_node	= NULL;
	
	for (index = 0; index < shape[level]; ++index) {
		if (l_node == NULL and r_node == NULL) {
			/*
			 * Both source lists are now empty.  Because the default value of the
			 * destination is already set appropriately we can now return.
			 */
			
			return;
			
		} else {
			// At least one list still has entries.
			
			if (op == EW_MUL) {
				// Special cases for multiplication.
				
				if (l_node == NULL and (l_default == 0 and d_default == 0)) {
					/* 
					 * The left hand list has run out of elements.  We don't need to add new
					 * values to the destination if l_default and d_default are both 0.
					 */
				
					return;
			
				} else if (r_node == NULL and (r_default == 0 and d_default == 0)) {
					/*
					 * The right hand list has run out of elements.  We don't need to add new
					 * values to the destination if r_default and d_default are both 0.
					 */
				
					return;
				}
				
			} else if (op == EW_DIV) {
				// Special cases for division.
				
				if (l_node == NULL and (l_default == 0 and d_default == 0)) {
					/* 
					 * The left hand list has run out of elements.  We don't need to add new
					 * values to the destination if l_default and d_default are both 0.
					 */
				
					return;
			
				} else if (r_node == NULL and (r_default == 0 and d_default == 0)) {
					/*
					 * The right hand list has run out of elements.  If the r_default
					 * value is 0 any further division will result in a SIGFPE.
					 */
				
					rb_raise(rb_eZeroDivError, "Cannot divide type by 0, would throw SIGFPE.");
				}
				
				// TODO: Add optimizations for addition and subtraction.
				
			}
			
			// We need to continue processing the lists.
			
			if (l_node == NULL and r_node->key == index) {
				/*
				 * One source list is empty, but the index has caught up to the key of
				 * the other list.
				 */
				
				if (level == last_level) {
				  tmp_result = ew_op_switch<op, LDType, RDType>(l_default, *reinterpret_cast<RDType*>(r_node->val));
				  std::cerr << "1. tmp_result = " << tmp_result << std::endl;
					
					if (tmp_result != d_default) {
						dest_node = nm::list::insert_helper(dest, dest_node, index, tmp_result);
					}
					
				} else {
					new_level = nm::list::create();
					dest_node = nm::list::insert_helper(dest, dest_node, index, new_level);
				
					ew_op_prime<op, LDType, RDType>(new_level, d_default,	&EMPTY_LIST, l_default,
						                              reinterpret_cast<LIST*>(r_node->val), r_default,
						                              shape, last_level, level + 1);
				}
				
				r_node = r_node->next;
				
			} else if (r_node == NULL and l_node->key == index) {
				/*
				 * One source list is empty, but the index has caught up to the key of
				 * the other list.
				 */
				
				if (level == last_level) {
				  tmp_result = ew_op_switch<op, LDType, RDType>(*reinterpret_cast<LDType*>(l_node->val), r_default);
				  std::cerr << "2. tmp_result = " << tmp_result << std::endl;

					if (tmp_result != d_default) {
						dest_node = nm::list::insert_helper(dest, dest_node, index, tmp_result);
					}
					
				} else {
					new_level = nm::list::create();
					dest_node = nm::list::insert_helper(dest, dest_node, index, new_level);
				
					ew_op_prime<op, LDType, RDType>(new_level, d_default,	reinterpret_cast<LIST*>(l_node->val), l_default,
						                              &EMPTY_LIST, r_default,	shape, last_level, level + 1);
				}
				
				l_node = l_node->next;
				
			} else if (l_node != NULL and r_node != NULL and index == std::min(l_node->key, r_node->key)) {
				/*
				 * Neither list is empty and our index has caught up to one of the
				 * source lists.
				 */
				
				if (l_node->key == r_node->key) {
					
					if (level == last_level) {
					  tmp_result = ew_op_switch<op, LDType, RDType>(*reinterpret_cast<LDType*>(l_node->val),*reinterpret_cast<RDType*>(r_node->val));
					  std::cerr << "3. tmp_result = " << tmp_result << std::endl;
						
						if (tmp_result != d_default) {
							dest_node = nm::list::insert_helper(dest, dest_node, index, tmp_result);
						}
						
					} else {
						new_level = nm::list::create();
						dest_node = nm::list::insert_helper(dest, dest_node, index, new_level);
					
						ew_op_prime<op, LDType, RDType>(new_level, d_default,
                                            reinterpret_cast<LIST*>(l_node->val), l_default,
                                            reinterpret_cast<LIST*>(r_node->val), r_default,
                                            shape, last_level, level + 1);
					}
				
					l_node = l_node->next;
					r_node = r_node->next;
			
				} else if (l_node->key < r_node->key) {
					// Advance the left node knowing that the default value is OK.
			
					l_node = l_node->next;
					 
				} else /* if (l_node->key > r_node->key) */ {
					// Advance the right node knowing that the default value is OK.
			
					r_node = r_node->next;
				}
				
			} else {
				/*
				 * Our index needs to catch up but the default value is OK.  This
				 * conditional is here only for documentation and should be optimized
				 * out.
				 */
			}
		}
	}
}

}} // end of namespace nm::list_storage

