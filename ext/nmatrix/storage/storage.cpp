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

void (* const STYPE_MARK[nm::NUM_STYPES])(void*) = {
	nm_dense_storage_mark,
	nm_list_storage_mark,
	nm_yale_storage_mark
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

  // allocate and copy shape
  size_t* shape = ALLOC_N(size_t, rhs->dim);
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

  return lhs;
}




/*
 * Create/allocate dense storage, copying into it the contents of a Yale matrix.
 */
template <typename LDType, typename RDType, typename RIType>
DENSE_STORAGE* create_from_yale_storage(const YALE_STORAGE* rhs, dtype_t l_dtype) {

  // Position in rhs->elements.
  RIType* rhs_ija = reinterpret_cast<RIType*>(rhs->ija);
  RDType* rhs_a   = reinterpret_cast<RDType*>(rhs->a);

  // Allocate and set shape.
  size_t* shape = ALLOC_N(size_t, rhs->dim);
  memcpy(shape, rhs->shape, rhs->dim * sizeof(size_t));

  DENSE_STORAGE* lhs = nm_dense_storage_create(l_dtype, shape, rhs->dim, NULL, 0);
  LDType* lhs_elements = reinterpret_cast<LDType*>(lhs->elements);

  // Position in dense to write to.
  size_t pos = 0;

  LDType LCAST_ZERO = rhs_a[rhs->shape[0]];

  // Walk through rows. For each entry we set in dense, increment pos.
  for (RIType i = 0; i < rhs->shape[0]; ++i) {

    // Position in yale array
    RIType ija = rhs_ija[i];

    if (ija == rhs_ija[i+1]) { // Check boundaries of row: is row empty?

			// Write zeros in each column.
			for (RIType j = 0; j < rhs->shape[1]; ++j) { // Move to next dense position.

        // Fill in zeros (except for diagonal)
        if (i == j) lhs_elements[pos] = rhs_a[i];
				else        lhs_elements[pos] = LCAST_ZERO;

				++pos;
      }

    } else {
      // Row contains entries: write those in each column, interspersed with zeros.
      RIType jj = rhs_ija[ija];

			for (size_t j = 0; j < rhs->shape[1]; ++j) {
        if (i == j) {
          lhs_elements[pos] = rhs_a[i];

        } else if (j == jj) {
          lhs_elements[pos] = rhs_a[ija]; // Copy from rhs.

          // Get next.
          ++ija;

          // Increment to next column ID (or go off the end).
          if (ija < rhs_ija[i+1]) jj = rhs_ija[ija];
          else               	    jj = rhs->shape[1];

        } else { // j < jj

          // Insert zero.
          lhs_elements[pos] = LCAST_ZERO;
        }

        // Move to next dense position.
        ++pos;
      }
    }
  }

  return lhs;
}


/*
 * Copy list contents into dense recursively.
 */
template <typename LDType, typename RDType>
static void cast_copy_list_contents(LDType* lhs, const LIST* rhs, RDType* default_val, size_t& pos, const size_t* shape, size_t dim, size_t max_elements, size_t recursions) {

  NODE *curr = rhs->first;
  int last_key = -1;

	for (size_t i = 0; i < shape[dim - 1 - recursions]; ++i, ++pos) {

    if (!curr || (curr->key > (size_t)(last_key+1))) {

      if (recursions == 0)  lhs[pos] = *default_val;
      else               		cast_copy_list_default<LDType,RDType>(lhs, default_val, pos, shape, dim, max_elements, recursions-1);

      ++last_key;

    } else {

      if (recursions == 0)  lhs[pos] = *reinterpret_cast<RDType*>(curr->val);
      else                	cast_copy_list_contents<LDType,RDType>(lhs, (const LIST*)(curr->val),
                                                                                         default_val, pos, shape, dim, max_elements, recursions-1);

      last_key = curr->key;
      curr     = curr->next;
    }
  }

  --pos;
}

/*
 * Copy a set of default values into dense.
 */
template <typename LDType,typename RDType>
static void cast_copy_list_default(LDType* lhs, RDType* default_val, size_t& pos, const size_t* shape, size_t dim, size_t max_elements, size_t recursions) {
	for (size_t i = 0; i < shape[dim - 1 - recursions]; ++i, ++pos) {

    if (recursions == 0)    lhs[pos] = *default_val;
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
LIST_STORAGE* create_from_dense_storage(const DENSE_STORAGE* rhs, dtype_t l_dtype) {

  LDType* l_default_val = ALLOC_N(LDType, 1);
  RDType* r_default_val = ALLOCA_N(RDType, 1); // clean up when finished with this function

  // allocate and copy shape and coords
  size_t *shape  = ALLOC_N(size_t, rhs->dim),
         *coords = ALLOC_N(size_t, rhs->dim);

  memcpy(shape, rhs->shape, rhs->dim * sizeof(size_t));
  memset(coords, 0, rhs->dim * sizeof(size_t));

  // set list default_val to 0
  if (l_dtype == RUBYOBJ)  	*l_default_val = INT2FIX(0);
  else    	                *l_default_val = 0;

  // need test default value for comparing to elements in dense matrix
  if (rhs->dtype == l_dtype)  	  *r_default_val = *l_default_val;
  else if (rhs->dtype == RUBYOBJ) *r_default_val = INT2FIX(0);
  else  	                        *r_default_val = 0;

  LIST_STORAGE* lhs = nm_list_storage_create(l_dtype, shape, rhs->dim, l_default_val);

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

  return lhs;
}



/*
 * Creation of list storage from yale storage.
 */
template <typename LDType, typename RDType, typename RIType>
LIST_STORAGE* create_from_yale_storage(const YALE_STORAGE* rhs, dtype_t l_dtype) {
  // allocate and copy shape
  size_t *shape = ALLOC_N(size_t, rhs->dim);
  shape[0] = rhs->shape[0]; shape[1] = rhs->shape[1];

  RDType* rhs_a    = reinterpret_cast<RDType*>(rhs->a);
  RDType R_ZERO    = rhs_a[ rhs->shape[0] ];

  // copy default value from the zero location in the Yale matrix
  LDType* default_val = ALLOC_N(LDType, 1);
  *default_val        = R_ZERO;

  LIST_STORAGE* lhs = nm_list_storage_create(l_dtype, shape, rhs->dim, default_val);

  if (rhs->dim != 2)    rb_raise(nm_eStorageTypeError, "Can only convert matrices of dim 2 from yale.");

  RIType* rhs_ija  = reinterpret_cast<RIType*>(rhs->ija);

  NODE *last_added = NULL, *last_row_added = NULL;

  // Walk through rows and columns as if RHS were a dense matrix
  for (RIType i = 0; i < rhs->shape[0]; ++i) {

    // Get boundaries of beginning and end of row
    RIType ija      = rhs_ija[i],
           ija_next = rhs_ija[i+1];

    // Are we going to need to add a diagonal for this row?
    bool add_diag = false;
    if (rhs_a[i] != R_ZERO) add_diag = true;

    if (ija < ija_next || add_diag) {

      LIST* curr_row = list::create();

      LDType* insert_val;

      while (ija < ija_next) {
        RIType jj = rhs_ija[ija]; // what column number is this?

        // Is there a nonzero diagonal item between the previously added item and the current one?
        if (jj > i && add_diag) {
          // Allocate and copy insertion value
          insert_val = ALLOC_N(LDType, 1);
          *insert_val        = rhs_a[i];

          // insert the item in the list at the appropriate location
          if (last_added) 	last_added = list::insert_after(last_added, i, insert_val);
          else            	last_added = list::insert(curr_row, false, i, insert_val);

					// don't add again!
          add_diag = false;
        }

        // now allocate and add the current item
        insert_val  = ALLOC_N(LDType, 1);
        *insert_val = rhs_a[ija];

        if (last_added)    	last_added = list::insert_after(last_added, jj, insert_val);
        else              	last_added = list::insert(curr_row, false, jj, insert_val);

        ++ija; // move to next entry in Yale matrix
      }

      if (add_diag) {
      	// still haven't added the diagonal.
        insert_val = ALLOC_N(LDType, 1);
        *insert_val        = rhs_a[i];

        // insert the item in the list at the appropriate location
        if (last_added)    	last_added = list::insert_after(last_added, i, insert_val);
        else              	last_added = list::insert(curr_row, false, i, insert_val);
      }

      // Now add the list at the appropriate location
      if (last_row_added)  	last_row_added = list::insert_after(last_row_added, i, curr_row);
      else                 	last_row_added = list::insert(lhs->rows, false, i, curr_row);
    }

		// end of walk through rows
  }

  return lhs;
}


/* Copy dense into lists recursively
 *
 * FIXME: This works, but could probably be cleaner (do we really need to pass coords around?)
 */
template <typename LDType, typename RDType>
static bool cast_copy_contents_dense(LIST* lhs, const RDType* rhs, RDType* zero, size_t& pos, size_t* coords, const size_t* shape, size_t dim, size_t recursions) {
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
        LDType* insert_value = ALLOC_N(LDType, 1);
        *insert_value        = (LDType)(rhs[pos]);

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

  coords[dim-1-recursions] = 0;
  --pos;

  return added;
}

} // end of namespace list_storage


namespace yale_storage { // FIXME: Move to yale.cpp
  /*
   * Creation of yale storage from dense storage.
   */
  template <typename LDType, typename RDType, typename LIType>
  YALE_STORAGE* create_from_dense_storage(const DENSE_STORAGE* rhs, dtype_t l_dtype) {
    if (rhs->dim != 2) rb_raise(nm_eStorageTypeError, "can only convert matrices of dim 2 to yale");

    LIType pos = 0;
    LIType ndnz = 0;

    RDType R_ZERO; // need zero for easier comparisons
    if (rhs->dtype == RUBYOBJ)  R_ZERO = INT2FIX(0);
    else                        R_ZERO = 0;


    RDType* rhs_elements = reinterpret_cast<RDType*>(rhs->elements);

    // First, count the non-diagonal nonzeros
    for (size_t i = rhs->shape[0]; i-- > 0;) {
      for (size_t j = rhs->shape[1]; j-- > 0;) {
        pos = rhs->stride[0]*(i + rhs->offset[0]) + rhs->stride[1]*(j + rhs->offset[1]);
        if (i != j && rhs_elements[pos] != R_ZERO)	++ndnz;

        // move forward 1 position in dense matrix elements array
      }
    }

    // Copy shape for yale construction
    size_t* shape = ALLOC_N(size_t, 2);
    shape[0] = rhs->shape[0];
    shape[1] = rhs->shape[1];

    size_t request_capacity = shape[0] + ndnz + 1;

    // Create with minimum possible capacity -- just enough to hold all of the entries
    YALE_STORAGE* lhs = nm_yale_storage_create(l_dtype, shape, 2, request_capacity);

    if (lhs->capacity < request_capacity)
      rb_raise(nm_eStorageTypeError, "conversion failed; capacity of %ld requested, max allowable is %ld", request_capacity, lhs->capacity);

    LDType* lhs_a     = reinterpret_cast<LDType*>(lhs->a);
    LIType* lhs_ija   = reinterpret_cast<LIType*>(lhs->ija);

    // Set the zero position in the yale matrix
    lhs_a[shape[0]] = R_ZERO;

    // Start just after the zero position.
    LIType ija = shape[0]+1;
    LIType i;
    pos        = 0;

    // Copy contents
    for (i = 0; i < rhs->shape[0]; ++i) {
      // indicate the beginning of a row in the IJA array
      lhs_ija[i]= ija;

      for (LIType j = 0; j < rhs->shape[1];  ++j) {
        pos = rhs->stride[0]*(i + rhs->offset[0]) + rhs->stride[1]*(j + rhs->offset[1]); // calc position with offsets

        if (i == j) { // copy to diagonal
          lhs_a[i]  = rhs_elements[pos];
        } else if (rhs_elements[pos] != R_ZERO) { // copy nonzero to LU
          lhs_ija[ija] = j; // write column index
          
          lhs_a[ija] = rhs_elements[pos];

          ++ija;
        }
      }
    }

    lhs_ija[shape[0]] = ija; // indicate the end of the last row
    lhs->ndnz = ndnz;

    return lhs;
  }

  /*
   * Creation of yale storage from list storage.
   */
  template <typename LDType, typename RDType, typename LIType>
  YALE_STORAGE* create_from_list_storage(const LIST_STORAGE* rhs, dtype_t l_dtype) {
    if (rhs->dim != 2) rb_raise(nm_eStorageTypeError, "can only convert matrices of dim 2 to yale");

    if ((rhs->dtype == RUBYOBJ and (*reinterpret_cast<RubyObject*>(rhs->default_val)) == RubyObject(INT2FIX(0)))
        || strncmp(reinterpret_cast<const char*>(rhs->default_val), "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", DTYPE_SIZES[rhs->dtype]))
      rb_raise(nm_eStorageTypeError, "list matrix must have default value of 0 to convert to yale");


    size_t ndnz = nm_list_storage_count_nd_elements(rhs);
    // Copy shape for yale construction
    size_t* shape = ALLOC_N(size_t, 2);
    shape[0] = rhs->shape[0];
    shape[1] = rhs->shape[1];

    size_t request_capacity = shape[0] + ndnz + 1;
    YALE_STORAGE* lhs = nm_yale_storage_create(l_dtype, shape, 2, request_capacity);

    if (lhs->capacity < request_capacity)
      rb_raise(nm_eStorageTypeError, "conversion failed; capacity of %ld requested, max allowable is %ld", request_capacity, lhs->capacity);

    // Initialize the A and IJA arrays
    init<LDType,LIType>(lhs);

    LIType* lhs_ija = reinterpret_cast<LIType*>(lhs->ija);
    LDType* lhs_a   = reinterpret_cast<LDType*>(lhs->a);

    LIType ija = lhs->shape[0]+1;

    // Copy contents 
    for (NODE* i_curr = rhs->rows->first; i_curr; i_curr = i_curr->next) {

      // Shrink refernce
      int i = i_curr->key - rhs->offset[0];
      if (i < 0 || i >= (int)rhs->shape[0]) continue;

      for (NODE* j_curr = ((LIST*)(i_curr->val))->first; j_curr; j_curr = j_curr->next) {
        
        // Shrink refernce
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



  STORAGE* nm_yale_storage_from_dense(const STORAGE* right, dtype_t l_dtype) {
    NAMED_LRI_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::create_from_dense_storage, YALE_STORAGE*, const DENSE_STORAGE* rhs, dtype_t l_dtype);

    itype_t itype = nm_yale_storage_itype((const YALE_STORAGE*)right);

    return (STORAGE*)ttable[l_dtype][right->dtype][itype]((const DENSE_STORAGE*)right, l_dtype);
  }

  STORAGE* nm_yale_storage_from_list(const STORAGE* right, dtype_t l_dtype) {
    NAMED_LRI_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::create_from_list_storage, YALE_STORAGE*, const LIST_STORAGE* rhs, dtype_t l_dtype);

    itype_t itype = nm_yale_storage_itype((const YALE_STORAGE*)right);

    return (STORAGE*)ttable[l_dtype][right->dtype][itype]((const LIST_STORAGE*)right, l_dtype);
  }

  STORAGE* nm_dense_storage_from_list(const STORAGE* right, dtype_t l_dtype) {
    NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::dense_storage::create_from_list_storage, DENSE_STORAGE*, const LIST_STORAGE* rhs, dtype_t l_dtype);

    return (STORAGE*)ttable[l_dtype][right->dtype]((const LIST_STORAGE*)right, l_dtype);
  }

  STORAGE* nm_dense_storage_from_yale(const STORAGE* right, dtype_t l_dtype) {
    NAMED_LRI_DTYPE_TEMPLATE_TABLE(ttable, nm::dense_storage::create_from_yale_storage, DENSE_STORAGE*, const YALE_STORAGE* rhs, dtype_t l_dtype);

    const YALE_STORAGE* casted_right = reinterpret_cast<const YALE_STORAGE*>(right);
    return reinterpret_cast<STORAGE*>(ttable[l_dtype][right->dtype][casted_right->itype](casted_right, l_dtype));
  }

  STORAGE* nm_list_storage_from_dense(const STORAGE* right, dtype_t l_dtype) {
    NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::list_storage::create_from_dense_storage, LIST_STORAGE*, const DENSE_STORAGE*, dtype_t);

    return (STORAGE*)ttable[l_dtype][right->dtype]((DENSE_STORAGE*)right, l_dtype);
  }

  STORAGE* nm_list_storage_from_yale(const STORAGE* right, dtype_t l_dtype) {
    NAMED_LRI_DTYPE_TEMPLATE_TABLE(ttable, nm::list_storage::create_from_yale_storage, LIST_STORAGE*, const YALE_STORAGE* rhs, dtype_t l_dtype);

    const YALE_STORAGE* casted_right = reinterpret_cast<const YALE_STORAGE*>(right);

    return (STORAGE*)ttable[l_dtype][right->dtype][casted_right->itype](casted_right, l_dtype);
  }

} // end of extern "C"

