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
// == storage.cpp
//
// Code that is used by or involves more then one storage type.

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

#include "data/data.h"

#include "storage.h"

#include "common.h"

/*
 * Macros
 */

/*
 * Global Variables
 */

extern "C" {

const char* const STYPE_NAMES[nm::NUM_STYPES] = {
	"dense",
	"list",
	"yale"
};

} // end extern "C" block

/*
 * Forward Declarations
 */

namespace nm {


/*
 * Functions
 */

/////////////////////////
// Templated Functions //
/////////////////////////

namespace dense_storage {

template <typename LDType, typename RDType>
static void cast_copy_list_contents(LDType* lhs, const LIST* rhs, RDType* default_val,
  size_t& pos, const size_t* shape, size_t dim, size_t max_elements, size_t recursions);

template <typename LDType, typename RDType>
static void cast_copy_list_default(LDType* lhs, RDType* default_val, size_t& pos,
  const size_t* shape, size_t dim, size_t max_elements, size_t recursions);

/*
 * Convert (by creating a copy) from list storage to dense storage.
 */
template <typename LDType, typename RDType>
DENSE_STORAGE* create_from_list_storage(const LIST_STORAGE* rhs, dtype_t l_dtype) {
  nm_list_storage_register(rhs);
  // allocate and copy shape
  size_t* shape = NM_ALLOC_N(size_t, rhs->dim);
  memcpy(shape, rhs->shape, rhs->dim * sizeof(size_t));

  DENSE_STORAGE* lhs = nm_dense_storage_create(l_dtype, shape, rhs->dim, NULL, 0);

  // Position in lhs->elements.
  size_t pos = 0;
  size_t max_elements = nm_storage_count_max_elements(rhs);

//static void dense_storage_cast_copy_list_contents_template(LDType* lhs, const LIST* rhs, RDType* default_val, size_t& pos, const size_t* shape, size_t dim, size_t max_elements, size_t recursions)
  // recursively copy the contents
  if (rhs->src == rhs)
    cast_copy_list_contents<LDType,RDType>(reinterpret_cast<LDType*>(lhs->elements),
                                         rhs->rows,
                                         reinterpret_cast<RDType*>(rhs->default_val),
                                         pos, shape, lhs->dim, max_elements, rhs->dim-1);
  else {
    LIST_STORAGE *tmp = nm_list_storage_copy(rhs);
    cast_copy_list_contents<LDType,RDType>(reinterpret_cast<LDType*>(lhs->elements),
                                         tmp->rows,
                                         reinterpret_cast<RDType*>(tmp->default_val),
                                         pos, shape, lhs->dim, max_elements, tmp->dim-1);
    nm_list_storage_delete(tmp);

  }
  nm_list_storage_unregister(rhs);

  return lhs;
}




/*
 * Create/allocate dense storage, copying into it the contents of a Yale matrix.
 */
template <typename LDType, typename RDType>
DENSE_STORAGE* create_from_yale_storage(const YALE_STORAGE* rhs, dtype_t l_dtype) {

  nm_yale_storage_register(rhs);
  // Position in rhs->elements.
  IType*  rhs_ija = reinterpret_cast<YALE_STORAGE*>(rhs->src)->ija;
  RDType* rhs_a   = reinterpret_cast<RDType*>(reinterpret_cast<YALE_STORAGE*>(rhs->src)->a);

  // Allocate and set shape.
  size_t* shape = NM_ALLOC_N(size_t, rhs->dim);
  shape[0] = rhs->shape[0];
  shape[1] = rhs->shape[1];

  DENSE_STORAGE* lhs = nm_dense_storage_create(l_dtype, shape, rhs->dim, NULL, 0);
  LDType* lhs_elements = reinterpret_cast<LDType*>(lhs->elements);

  // Position in dense to write to.
  size_t pos = 0;

  LDType LCAST_ZERO = rhs_a[rhs->src->shape[0]];

  // Walk through rows. For each entry we set in dense, increment pos.
  for (size_t i = 0; i < shape[0]; ++i) {
    IType ri = i + rhs->offset[0];

    if (rhs_ija[ri] == rhs_ija[ri+1]) { // Check boundaries of row: is row empty? (Yes.)

			// Write zeros in each column.
			for (size_t j = 0; j < shape[1]; ++j) { // Move to next dense position.

        // Fill in zeros and copy the diagonal entry for this empty row.
        if (ri == j + rhs->offset[1]) lhs_elements[pos] = static_cast<LDType>(rhs_a[ri]);
				else                          lhs_elements[pos] = LCAST_ZERO;

				++pos;
      }

    } else {  // Row contains entries: write those in each column, interspersed with zeros.

      // Get the first ija position of the row (as sliced)
      IType ija = nm::yale_storage::binary_search_left_boundary(rhs, rhs_ija[ri], rhs_ija[ri+1]-1, rhs->offset[1]);

      // What column is it?
      IType next_stored_rj = rhs_ija[ija];

			for (size_t j = 0; j < shape[1]; ++j) {
			  IType rj = j + rhs->offset[1];

        if (rj == ri) { // at a diagonal in RHS
          lhs_elements[pos] = static_cast<LDType>(rhs_a[ri]);

        } else if (rj == next_stored_rj) { // column ID was found in RHS
          lhs_elements[pos] = static_cast<LDType>(rhs_a[ija]); // Copy from rhs.

          // Get next.
          ++ija;

          // Increment to next column ID (or go off the end).
          if (ija < rhs_ija[ri+1]) next_stored_rj = rhs_ija[ija];
          else               	     next_stored_rj = rhs->src->shape[1];

        } else { // rj < next_stored_rj

          // Insert zero.
          lhs_elements[pos] = LCAST_ZERO;
        }

        // Move to next dense position.
        ++pos;
      }
    }
  }
  nm_yale_storage_unregister(rhs);

  return lhs;
}


/*
 * Copy list contents into dense recursively.
 */
template <typename LDType, typename RDType>
static void cast_copy_list_contents(LDType* lhs, const LIST* rhs, RDType* default_val, size_t& pos, const size_t* shape, size_t dim, size_t max_elements, size_t recursions) {

  NODE *curr = rhs->first;
  int last_key = -1;

  nm_list_storage_register_list(rhs, recursions);

  for (size_t i = 0; i < shape[dim - 1 - recursions]; ++i, ++pos) {

    if (!curr || (curr->key > (size_t)(last_key+1))) {

      if (recursions == 0)  lhs[pos] = static_cast<LDType>(*default_val);
      else               		cast_copy_list_default<LDType,RDType>(lhs, default_val, pos, shape, dim, max_elements, recursions-1);

      ++last_key;

    } else {

      if (recursions == 0)  lhs[pos] = static_cast<LDType>(*reinterpret_cast<RDType*>(curr->val));
      else                	cast_copy_list_contents<LDType,RDType>(lhs, (const LIST*)(curr->val),
                                                                                         default_val, pos, shape, dim, max_elements, recursions-1);

      last_key = curr->key;
      curr     = curr->next;
    }
  }

  nm_list_storage_unregister_list(rhs, recursions);

  --pos;
}

/*
 * Copy a set of default values into dense.
 */
template <typename LDType,typename RDType>
static void cast_copy_list_default(LDType* lhs, RDType* default_val, size_t& pos, const size_t* shape, size_t dim, size_t max_elements, size_t recursions) {
  for (size_t i = 0; i < shape[dim - 1 - recursions]; ++i, ++pos) {

    if (recursions == 0)    lhs[pos] = static_cast<LDType>(*default_val);
    else                  	cast_copy_list_default<LDType,RDType>(lhs, default_val, pos, shape, dim, max_elements, recursions-1);

  }

  --pos;
}


} // end of namespace dense_storage

namespace list_storage {


template <typename LDType, typename RDType>
static bool cast_copy_contents_dense(LIST* lhs, const RDType* rhs, RDType* zero, size_t& pos, size_t* coords, const size_t* shape, size_t dim, size_t recursions);

/*
 * Creation of list storage from dense storage.
 */
template <typename LDType, typename RDType>
LIST_STORAGE* create_from_dense_storage(const DENSE_STORAGE* rhs, dtype_t l_dtype, void* init) {
  nm_dense_storage_register(rhs);

  LDType* l_default_val = NM_ALLOC_N(LDType, 1);
  RDType* r_default_val = NM_ALLOCA_N(RDType, 1); // clean up when finished with this function

  // allocate and copy shape and coords
  size_t *shape  = NM_ALLOC_N(size_t, rhs->dim),
         *coords = NM_ALLOC_N(size_t, rhs->dim);

  memcpy(shape, rhs->shape, rhs->dim * sizeof(size_t));
  memset(coords, 0, rhs->dim * sizeof(size_t));

  // set list default_val to 0
  if (init) *l_default_val = *reinterpret_cast<LDType*>(init);
  else {
    if (l_dtype == RUBYOBJ)  	*l_default_val = INT2FIX(0);
    else    	                *l_default_val = 0;
  }

  // need test default value for comparing to elements in dense matrix
  if (rhs->dtype == l_dtype || rhs->dtype != RUBYOBJ) *r_default_val = static_cast<RDType>(*l_default_val);
  else                                                *r_default_val = rubyobj_from_cval(l_default_val, l_dtype);


  LIST_STORAGE* lhs = nm_list_storage_create(l_dtype, shape, rhs->dim, l_default_val);

  nm_list_storage_register(lhs);

  size_t pos = 0;

  if (rhs->src == rhs)
    list_storage::cast_copy_contents_dense<LDType,RDType>(lhs->rows,
                                                          reinterpret_cast<const RDType*>(rhs->elements),
                                                        r_default_val,
                                                        pos, coords, rhs->shape, rhs->dim, rhs->dim - 1);
  else {
    DENSE_STORAGE* tmp = nm_dense_storage_copy(rhs);
    list_storage::cast_copy_contents_dense<LDType,RDType>(lhs->rows,
                                                          reinterpret_cast<const RDType*>(tmp->elements),
                                                        r_default_val,
                                                        pos, coords, rhs->shape, rhs->dim, rhs->dim - 1);

    nm_dense_storage_delete(tmp);
  }

  nm_list_storage_unregister(lhs);
  nm_dense_storage_unregister(rhs);

  return lhs;
}



/*
 * Creation of list storage from yale storage.
 */
template <typename LDType, typename RDType>
LIST_STORAGE* create_from_yale_storage(const YALE_STORAGE* rhs, dtype_t l_dtype) {
  // allocate and copy shape
  nm_yale_storage_register(rhs);

  size_t *shape = NM_ALLOC_N(size_t, rhs->dim);
  shape[0] = rhs->shape[0]; shape[1] = rhs->shape[1];

  RDType* rhs_a    = reinterpret_cast<RDType*>(reinterpret_cast<YALE_STORAGE*>(rhs->src)->a);
  RDType R_ZERO    = rhs_a[ rhs->src->shape[0] ];

  // copy default value from the zero location in the Yale matrix
  LDType* default_val = NM_ALLOC_N(LDType, 1);
  *default_val        = static_cast<LDType>(R_ZERO);

  LIST_STORAGE* lhs = nm_list_storage_create(l_dtype, shape, rhs->dim, default_val);

  if (rhs->dim != 2)    rb_raise(nm_eStorageTypeError, "Can only convert matrices of dim 2 from yale.");

  IType* rhs_ija  = reinterpret_cast<YALE_STORAGE*>(rhs->src)->ija;

  NODE *last_row_added = NULL;
  // Walk through rows and columns as if RHS were a dense matrix
  for (IType i = 0; i < shape[0]; ++i) {
    IType ri = i + rhs->offset[0];

    NODE *last_added = NULL;

    // Get boundaries of beginning and end of row
    IType ija      = rhs_ija[ri],
          ija_next = rhs_ija[ri+1];

    // Are we going to need to add a diagonal for this row?
    bool add_diag = false;
    if (rhs_a[ri] != R_ZERO) add_diag = true; // non-zero and located within the bounds of the slice

    if (ija < ija_next || add_diag) {
      ija = nm::yale_storage::binary_search_left_boundary(rhs, ija, ija_next-1, rhs->offset[1]);

      LIST* curr_row = list::create();

      LDType* insert_val;

      while (ija < ija_next) {
        // Find first column in slice
        IType rj = rhs_ija[ija];
        IType j  = rj - rhs->offset[1];

        // Is there a nonzero diagonal item between the previously added item and the current one?
        if (rj > ri && add_diag) {
          // Allocate and copy insertion value
          insert_val  = NM_ALLOC_N(LDType, 1);
          *insert_val = static_cast<LDType>(rhs_a[ri]);

          // Insert the item in the list at the appropriate location.
          // What is the appropriate key? Well, it's definitely right(i)==right(j), but the
          // rj index has already been advanced past ri. So we should treat ri as the column and
          // subtract offset[1].
          if (last_added) 	last_added = list::insert_after(last_added, ri - rhs->offset[1], insert_val);
          else            	last_added = list::insert(curr_row, false,  ri - rhs->offset[1], insert_val);

					// don't add again!
          add_diag = false;
        }

        // now allocate and add the current item
        insert_val  = NM_ALLOC_N(LDType, 1);
        *insert_val = static_cast<LDType>(rhs_a[ija]);

        if (last_added)    	last_added = list::insert_after(last_added, j, insert_val);
        else              	last_added = list::insert(curr_row, false, j, insert_val);

        ++ija; // move to next entry in Yale matrix
      }

      if (add_diag) {

      	// still haven't added the diagonal.
        insert_val         = NM_ALLOC_N(LDType, 1);
        *insert_val        = static_cast<LDType>(rhs_a[ri]);

        // insert the item in the list at the appropriate location
        if (last_added)    	last_added = list::insert_after(last_added, ri - rhs->offset[1], insert_val);
        else              	last_added = list::insert(curr_row, false, ri - rhs->offset[1], insert_val);

        // no need to set add_diag to false because it'll be reset automatically in next iteration.
      }

      // Now add the list at the appropriate location
      if (last_row_added)   last_row_added = list::insert_after(last_row_added, i, curr_row);
      else                  last_row_added = list::insert(lhs->rows, false, i, curr_row);
    }

		// end of walk through rows
  }

  nm_yale_storage_unregister(rhs);

  return lhs;
}


/* Copy dense into lists recursively
 *
 * FIXME: This works, but could probably be cleaner (do we really need to pass coords around?)
 */
template <typename LDType, typename RDType>
static bool cast_copy_contents_dense(LIST* lhs, const RDType* rhs, RDType* zero, size_t& pos, size_t* coords, const size_t* shape, size_t dim, size_t recursions) {

  nm_list_storage_register_list(lhs, recursions);

  NODE *prev = NULL;
  LIST *sub_list;
  bool added = false, added_list = false;
  //void* insert_value;

  for (coords[dim-1-recursions] = 0; coords[dim-1-recursions] < shape[dim-1-recursions]; ++coords[dim-1-recursions], ++pos) {

    if (recursions == 0) {
    	// create nodes

      if (rhs[pos] != *zero) {
      	// is not zero

        // Create a copy of our value that we will insert in the list
        LDType* insert_value = NM_ALLOC_N(LDType, 1);
        *insert_value        = static_cast<LDType>(rhs[pos]);

        if (!lhs->first)    prev = list::insert(lhs, false, coords[dim-1-recursions], insert_value);
        else               	prev = list::insert_after(prev, coords[dim-1-recursions], insert_value);

        added = true;
      }
      // no need to do anything if the element is zero

    } else { // create lists
      // create a list as if there's something in the row in question, and then delete it if nothing turns out to be there
      sub_list = list::create();

      added_list = list_storage::cast_copy_contents_dense<LDType,RDType>(sub_list, rhs, zero, pos, coords, shape, dim, recursions-1);

      if (!added_list)      	list::del(sub_list, recursions-1);
      else if (!lhs->first)  	prev = list::insert(lhs, false, coords[dim-1-recursions], sub_list);
      else                  	prev = list::insert_after(prev, coords[dim-1-recursions], sub_list);

      // added = (added || added_list);
    }
  }

  nm_list_storage_unregister_list(lhs, recursions);

  coords[dim-1-recursions] = 0;
  --pos;

  return added;
}

} // end of namespace list_storage


namespace yale_storage { // FIXME: Move to yale.cpp
  /*
   * Creation of yale storage from dense storage.
   */
  template <typename LDType, typename RDType>
  YALE_STORAGE* create_from_dense_storage(const DENSE_STORAGE* rhs, dtype_t l_dtype, void* init) {

    if (rhs->dim != 2) rb_raise(nm_eStorageTypeError, "can only convert matrices of dim 2 to yale");

    nm_dense_storage_register(rhs);

    IType pos = 0;
    IType ndnz = 0;

    // We need a zero value. This should nearly always be zero, but sometimes you might want false or nil.
    LDType    L_INIT(0);
    if (init) {
      if (l_dtype == RUBYOBJ) L_INIT = *reinterpret_cast<VALUE*>(init);
      else                    L_INIT = *reinterpret_cast<LDType*>(init);
    }
    RDType R_INIT = static_cast<RDType>(L_INIT);

    RDType* rhs_elements = reinterpret_cast<RDType*>(rhs->elements);

    // First, count the non-diagonal nonzeros
    for (size_t i = rhs->shape[0]; i-- > 0;) {
      for (size_t j = rhs->shape[1]; j-- > 0;) {
        pos = rhs->stride[0]*(i + rhs->offset[0]) + rhs->stride[1]*(j + rhs->offset[1]);
        if (i != j && rhs_elements[pos] != R_INIT)	++ndnz;

        // move forward 1 position in dense matrix elements array
      }
    }

    // Copy shape for yale construction
    size_t* shape = NM_ALLOC_N(size_t, 2);
    shape[0] = rhs->shape[0];
    shape[1] = rhs->shape[1];

    size_t request_capacity = shape[0] + ndnz + 1;

    // Create with minimum possible capacity -- just enough to hold all of the entries
    YALE_STORAGE* lhs = nm_yale_storage_create(l_dtype, shape, 2, request_capacity);

    if (lhs->capacity < request_capacity)
      rb_raise(nm_eStorageTypeError, "conversion failed; capacity of %ld requested, max allowable is %ld", (unsigned long)request_capacity, (unsigned long)(lhs->capacity));

    LDType* lhs_a     = reinterpret_cast<LDType*>(lhs->a);
    IType* lhs_ija    = lhs->ija;

    // Set the zero position in the yale matrix
    lhs_a[shape[0]]   = L_INIT;

    // Start just after the zero position.
    IType ija = shape[0]+1;
    pos       = 0;

    // Copy contents
    for (IType i = 0; i < rhs->shape[0]; ++i) {
      // indicate the beginning of a row in the IJA array
      lhs_ija[i] = ija;

      for (IType j = 0; j < rhs->shape[1];  ++j) {
        pos = rhs->stride[0] * (i + rhs->offset[0]) + rhs->stride[1] * (j + rhs->offset[1]); // calc position with offsets

        if (i == j) { // copy to diagonal
          lhs_a[i]     = static_cast<LDType>(rhs_elements[pos]);
        } else if (rhs_elements[pos] != R_INIT) { // copy nonzero to LU
          lhs_ija[ija] = j; // write column index
          lhs_a[ija]   = static_cast<LDType>(rhs_elements[pos]);

          ++ija;
        }
      }
    }

    lhs_ija[shape[0]] = ija; // indicate the end of the last row
    lhs->ndnz = ndnz;

    nm_dense_storage_unregister(rhs);

    return lhs;
  }

  /*
   * Creation of yale storage from list storage.
   */
  template <typename LDType, typename RDType>
  YALE_STORAGE* create_from_list_storage(const LIST_STORAGE* rhs, nm::dtype_t l_dtype) {
    if (rhs->dim != 2) rb_raise(nm_eStorageTypeError, "can only convert matrices of dim 2 to yale");

    if (rhs->dtype == RUBYOBJ) {
      VALUE init_val = *reinterpret_cast<VALUE*>(rhs->default_val);
      if (rb_funcall(init_val, rb_intern("!="), 1, Qnil) == Qtrue && rb_funcall(init_val, rb_intern("!="), 1, Qfalse) == Qtrue && rb_funcall(init_val, rb_intern("!="), 1, INT2FIX(0)) == Qtrue)
        rb_raise(nm_eStorageTypeError, "list matrix of Ruby objects must have default value equal to 0, nil, or false to convert to yale");
    } else if (strncmp(reinterpret_cast<const char*>(rhs->default_val), "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", DTYPE_SIZES[rhs->dtype]))
      rb_raise(nm_eStorageTypeError, "list matrix of non-Ruby objects must have default value of 0 to convert to yale");

    nm_list_storage_register(rhs);

    size_t ndnz = nm_list_storage_count_nd_elements(rhs);
    // Copy shape for yale construction
    size_t* shape = NM_ALLOC_N(size_t, 2);
    shape[0] = rhs->shape[0];
    shape[1] = rhs->shape[1];

    size_t request_capacity = shape[0] + ndnz + 1;
    YALE_STORAGE* lhs = nm_yale_storage_create(l_dtype, shape, 2, request_capacity);

    if (lhs->capacity < request_capacity)
      rb_raise(nm_eStorageTypeError, "conversion failed; capacity of %ld requested, max allowable is %ld", (unsigned long)request_capacity, (unsigned long)(lhs->capacity));

    // Initialize the A and IJA arrays
    init<LDType>(lhs, rhs->default_val);

    IType*  lhs_ija = lhs->ija;
    LDType* lhs_a   = reinterpret_cast<LDType*>(lhs->a);

    IType ija = lhs->shape[0]+1;

    // Copy contents 
    for (NODE* i_curr = rhs->rows->first; i_curr; i_curr = i_curr->next) {

      // Shrink reference
      int i = i_curr->key - rhs->offset[0];
      if (i < 0 || i >= (int)rhs->shape[0]) continue;

      for (NODE* j_curr = ((LIST*)(i_curr->val))->first; j_curr; j_curr = j_curr->next) {
        
        // Shrink reference
        int j = j_curr->key - rhs->offset[1];
        if (j < 0 || j >= (int)rhs->shape[1]) continue;

        LDType cast_jcurr_val = *reinterpret_cast<RDType*>(j_curr->val);
        if (i_curr->key - rhs->offset[0] == j_curr->key - rhs->offset[1])
          lhs_a[i_curr->key - rhs->offset[0]] = cast_jcurr_val; // set diagonal
        else {
          lhs_ija[ija] = j_curr->key - rhs->offset[1];    // set column value

          lhs_a[ija]   = cast_jcurr_val;                      // set cell value

          ++ija;
          // indicate the beginning of a row in the IJA array
          for (size_t i = i_curr->key - rhs->offset[0] + 1; i < rhs->shape[0] + rhs->offset[0]; ++i) {
            lhs_ija[i] = ija;
          }

        }
      }

    }
    
    lhs_ija[rhs->shape[0]] = ija; // indicate the end of the last row
    lhs->ndnz = ndnz;

    nm_list_storage_unregister(rhs);

    return lhs;
  }

} // end of namespace yale_storage
} // end of namespace nm

extern "C" {

  /*
   * The following functions represent stype casts -- conversions from one
   * stype to another. Each of these is the C accessor for a templated C++
   * function.
   */


  STORAGE* nm_yale_storage_from_dense(const STORAGE* right, nm::dtype_t l_dtype, void* init) {
    NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::create_from_dense_storage, YALE_STORAGE*, const DENSE_STORAGE* rhs, nm::dtype_t l_dtype, void*);

    if (!ttable[l_dtype][right->dtype]) {
      rb_raise(nm_eDataTypeError, "casting between these dtypes is undefined");
      return NULL;
    }

    return (STORAGE*)ttable[l_dtype][right->dtype]((const DENSE_STORAGE*)right, l_dtype, init);
  }

  STORAGE* nm_yale_storage_from_list(const STORAGE* right, nm::dtype_t l_dtype, void* dummy) {
    NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::create_from_list_storage, YALE_STORAGE*, const LIST_STORAGE* rhs, nm::dtype_t l_dtype);

    if (!ttable[l_dtype][right->dtype]) {
      rb_raise(nm_eDataTypeError, "casting between these dtypes is undefined");
      return NULL;
    }

    return (STORAGE*)ttable[l_dtype][right->dtype]((const LIST_STORAGE*)right, l_dtype);
  }

  STORAGE* nm_dense_storage_from_list(const STORAGE* right, nm::dtype_t l_dtype, void* dummy) {
    NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::dense_storage::create_from_list_storage, DENSE_STORAGE*, const LIST_STORAGE* rhs, nm::dtype_t l_dtype);

    if (!ttable[l_dtype][right->dtype]) {
      rb_raise(nm_eDataTypeError, "casting between these dtypes is undefined");
      return NULL;
    }

    return (STORAGE*)ttable[l_dtype][right->dtype]((const LIST_STORAGE*)right, l_dtype);
  }

  STORAGE* nm_dense_storage_from_yale(const STORAGE* right, nm::dtype_t l_dtype, void* dummy) {
    NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::dense_storage::create_from_yale_storage, DENSE_STORAGE*, const YALE_STORAGE* rhs, nm::dtype_t l_dtype);

    const YALE_STORAGE* casted_right = reinterpret_cast<const YALE_STORAGE*>(right);

    if (!ttable[l_dtype][right->dtype]) {
      rb_raise(nm_eDataTypeError, "casting between these dtypes is undefined");
      return NULL;
    }

    return reinterpret_cast<STORAGE*>(ttable[l_dtype][right->dtype](casted_right, l_dtype));
  }

  STORAGE* nm_list_storage_from_dense(const STORAGE* right, nm::dtype_t l_dtype, void* init) {
    NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::list_storage::create_from_dense_storage, LIST_STORAGE*, const DENSE_STORAGE*, nm::dtype_t, void*);

    if (!ttable[l_dtype][right->dtype]) {
      rb_raise(nm_eDataTypeError, "casting between these dtypes is undefined");
      return NULL;
    }

    return (STORAGE*)ttable[l_dtype][right->dtype]((DENSE_STORAGE*)right, l_dtype, init);
  }

  STORAGE* nm_list_storage_from_yale(const STORAGE* right, nm::dtype_t l_dtype, void* dummy) {
    NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::list_storage::create_from_yale_storage, LIST_STORAGE*, const YALE_STORAGE* rhs, nm::dtype_t l_dtype);

    const YALE_STORAGE* casted_right = reinterpret_cast<const YALE_STORAGE*>(right);

    if (!ttable[l_dtype][right->dtype]) {
      rb_raise(nm_eDataTypeError, "casting between these dtypes is undefined");
      return NULL;
    }

    return (STORAGE*)ttable[l_dtype][right->dtype](casted_right, l_dtype);
  }

} // end of extern "C"

