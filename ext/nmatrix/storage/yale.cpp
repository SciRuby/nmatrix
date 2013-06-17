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
// SciRuby is Copyright (c) 2010 - 2013, Ruby Science Foundation
// NMatrix is Copyright (c) 2013, Ruby Science Foundation
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
// == yale.c
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

/*
 * Standard Includes
 */

#include <ruby.h>
#include <algorithm>  // std::min
#include <cstdio>     // std::fprintf
#include <iostream>

/*
 * Project Includes
 */

// #include "types.h"
#include "util/math.h"

#include "data/data.h"

#include "common.h"
#include "yale.h"

#include "nmatrix.h"
#include "ruby_constants.h"

/*
 * Macros
 */

#ifndef NM_MAX
#define NM_MAX(a,b) (((a)>(b))?(a):(b))
#define NM_MIN(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef NM_MAX_ITYPE
#define NM_MAX_ITYPE(a,b) ((static_cast<int8_t>(a) > static_cast<int8_t>(b)) ? static_cast<nm::itype_t>(a) : static_cast<nm::itype_t>(b))
#define NM_MIN_ITYPE(a,b) ((static_cast<int8_t>(a) < static_cast<int8_t>(b)) ? static_cast<nm::itype_t>(a) : static_cast<nm::itype_t>(b))
#endif

/*
 * Forward Declarations
 */

extern "C" {
  static YALE_STORAGE*  nm_copy_alloc_struct(const YALE_STORAGE* rhs, const nm::dtype_t new_dtype, const size_t new_capacity, const size_t new_size);
  static YALE_STORAGE*	alloc(nm::dtype_t dtype, size_t* shape, size_t dim, nm::itype_t min_itype);

  /* Ruby-accessible functions */
  static VALUE nm_size(VALUE self);
  static VALUE nm_a(int argc, VALUE* argv, VALUE self);
  static VALUE nm_d(int argc, VALUE* argv, VALUE self);
  static VALUE nm_lu(VALUE self);
  static VALUE nm_ia(VALUE self);
  static VALUE nm_ja(VALUE self);
  static VALUE nm_ija(int argc, VALUE* argv, VALUE self);

  static VALUE nm_nd_row(int argc, VALUE* argv, VALUE self);
  static VALUE nm_vector_insert(int argc, VALUE* argv, VALUE self);


} // end extern "C" block

namespace nm { namespace yale_storage {

template <typename DType, typename IType>
static bool						ndrow_is_empty(const YALE_STORAGE* s, IType ija, const IType ija_next);

template <typename LDType, typename RDType, typename IType>
static bool						ndrow_eqeq_ndrow(const YALE_STORAGE* l, const YALE_STORAGE* r, IType l_ija, const IType l_ija_next, IType r_ija, const IType r_ija_next);

template <typename LDType, typename RDType, typename IType>
static bool           eqeq(const YALE_STORAGE* left, const YALE_STORAGE* right);

template <typename IType>
static YALE_STORAGE*	copy_alloc_struct(const YALE_STORAGE* rhs, const dtype_t new_dtype, const size_t new_capacity, const size_t new_size);

template <typename IType>
static void						increment_ia_after(YALE_STORAGE* s, IType ija_size, IType i, IType n);

template <typename IType>
static void           c_increment_ia_after(YALE_STORAGE* s, size_t ija_size, size_t i, size_t n) {
  increment_ia_after<IType>(s, ija_size, i, n);
}

template <typename IType>
static IType				  insert_search(YALE_STORAGE* s, IType left, IType right, IType key, bool* found);

template <typename DType, typename IType>
static char           vector_insert(YALE_STORAGE* s, size_t pos, size_t* j, void* val_, size_t n, bool struct_only);

template <typename DType, typename IType>
static char           vector_insert_resize(YALE_STORAGE* s, size_t current_size, size_t pos, size_t* j, size_t n, bool struct_only);

template <typename nm::ewop_t op, typename IType, typename DType>
YALE_STORAGE* ew_op(const YALE_STORAGE* left, const YALE_STORAGE* right, dtype_t dtype);

/*
 * Functions
 */

/*
 * Copy a vector from one IType or DType to another.
 */
template <typename LType, typename RType>
static inline void copy_recast_vector(const void* in_, void* out_, size_t length) {
  const RType* in = reinterpret_cast<const RType*>(in_);
  LType* out      = reinterpret_cast<LType*>(out_);
  for (size_t i = 0; i < length; ++i) {
    out[i] = in[i];
  }
  out;
}


static inline void copy_recast_itype_vector(const void* in, nm::itype_t in_itype, void* out, nm::itype_t out_itype, size_t length) {
  NAMED_LR_ITYPE_TEMPLATE_TABLE(ttable, copy_recast_vector, void, const void* in_, void* out_, size_t length);

  ttable[out_itype][in_itype](in, out, length);
}


/*
 * Create Yale storage from IA, JA, and A vectors given in Old Yale format (probably from a file, since NMatrix only uses
 * new Yale for its storage).
 *
 * This function is needed for Matlab .MAT v5 IO.
 */
template <typename LDType, typename RDType, typename IType>
YALE_STORAGE* create_from_old_yale(dtype_t dtype, size_t* shape, void* r_ia, void* r_ja, void* r_a) {
  IType*  ir = reinterpret_cast<IType*>(r_ia);
  IType*  jr = reinterpret_cast<IType*>(r_ja);
  RDType* ar = reinterpret_cast<RDType*>(r_a);

  // Read through ia and ja and figure out the ndnz (non-diagonal non-zeros) count.
  size_t ndnz = 0, i, p, p_next;

  for (i = 0; i < shape[0]; ++i) { // Walk down rows
    for (p = ir[i], p_next = ir[i+1]; p < p_next; ++p) { // Now walk through columns

      if (i != jr[p]) ++ndnz; // entry is non-diagonal and probably nonzero

    }
  }

  // Having walked through the matrix, we now go about allocating the space for it.
  YALE_STORAGE* s = alloc(dtype, shape, 2, UINT8);

  s->capacity = shape[0] + ndnz + 1;
  s->ndnz     = ndnz;

  // Setup IJA and A arrays
  s->ija = ALLOC_N( IType, s->capacity );
  s->a   = ALLOC_N( LDType, s->capacity );
  IType* ijl    = reinterpret_cast<IType*>(s->ija);
  LDType* al    = reinterpret_cast<LDType*>(s->a);

  // set the diagonal to zero -- this prevents uninitialized values from popping up.
  for (size_t index = 0; index < shape[0]; ++index) {
    al[index] = 0;
  }

  // Figure out where to start writing JA in IJA:
  size_t pp = s->shape[0]+1;

  // Find beginning of first row
  p = ir[0];

  // Now fill the arrays
  for (i = 0; i < s->shape[0]; ++i) {

    // Set the beginning of the row (of output)
    ijl[i] = pp;

    // Now walk through columns, starting at end of row (of input)
    for (size_t p_next = ir[i+1]; p < p_next; ++p, ++pp) {

      if (i == jr[p]) { // diagonal

        al[i] = ar[p];
        --pp;

      } else {          // nondiagonal

        ijl[pp] = jr[p];
        al[pp]  = ar[p];

      }
    }
  }

  ijl[i] = pp; // Set the end of the last row

  // Set the zero position for our output matrix
  al[i] = 0;

  return s;
}


/*
 * Take two Yale storages and merge them into a new Yale storage.
 *
 * Uses the left as a template for the creation of a new one.
 */
template <typename DType, typename IType>
YALE_STORAGE* create_merged(const YALE_STORAGE* left, const YALE_STORAGE* right) {
  char ins_type;

  size_t size = get_size<IType>(left);

  // s represents the resulting storage
  YALE_STORAGE* s = copy_alloc_struct<IType>(left, left->dtype, NM_MAX(left->capacity, right->capacity), size);

  IType* sija = reinterpret_cast<IType*>(s->ija);
  IType* rija = reinterpret_cast<IType*>(right->ija);

  // set the element between D and LU (the boundary in A), which should be 0.
  reinterpret_cast<DType*>(s->a)[s->shape[0]] = reinterpret_cast<DType*>(left->a)[left->shape[0]];

  if (right && right != left) {
  	// some operations are unary and don't need this; others are x+x and don't need this

  	for (IType i = 0; i < s->shape[0]; ++i) {

      IType ija       = sija[i];
      IType ija_next  = sija[i+1];

      for (IType r_ija = rija[i]; r_ija < rija[i+1]; ++r_ija) {

        size_t ja = sija[ija]; // insert expects a size_t

        if (ija == ija_next) {
        	// destination row is empty
          ins_type = vector_insert<DType,IType>(s, ija, &ja, NULL, 1, true);
          increment_ia_after<IType>(s, s->shape[0], i, 1);
          ++(s->ndnz);
          ++ija;

          if (ins_type == 'i') ++ija_next;

        } else {
          bool found;

        	// merge positions into destination row
          IType pos = insert_search<IType>(s, ija, ija_next-1, sija[ija], &found);

          if (!found) {
            vector_insert<DType,IType>(s, pos, &ja, NULL, 1, true);
            increment_ia_after<IType>(s, s->shape[0], i, 1);
            ++(s->ndnz);

            if (ins_type == 'i') ++ija_next;
          }

          // can now set a left boundary for the next search
          ija = pos + 1;
        }
      }
    }
  }

  return s;
}


/*
 * Empty the matrix by initializing the IJA vector and setting the diagonal to 0.
 *
 * Called when most YALE_STORAGE objects are created.
 */
template <typename DType, typename IType>
void init(YALE_STORAGE* s) {
  IType IA_INIT = s->shape[0] + 1;

  IType* ija = reinterpret_cast<IType*>(s->ija);
  // clear out IJA vector
  for (IType i = 0; i < IA_INIT; ++i) {
    ija[i] = IA_INIT; // set initial values for IJA
  }

  clear_diagonal_and_zero<DType>(s);
}

size_t max_size(YALE_STORAGE* s) {
  size_t result = s->shape[0]*s->shape[1] + 1;
  if (s->shape[0] > s->shape[1])
    result += s->shape[0] - s->shape[1];

  return result;
}
///////////////
// Accessors //
///////////////

/*
 * Returns a slice of YALE_STORAGE object by copy
 */
template <typename DType,typename IType>
void* get(YALE_STORAGE* storage, SLICE* slice) {
  
  size_t *offset = slice->coords;
  // Copy shape for yale construction
  size_t* shape = ALLOC_N(size_t, 2);
  shape[0] = slice->lengths[0];
  shape[1] = slice->lengths[1];

  IType* src_ija = reinterpret_cast<IType*>(storage->ija);
  DType* src_a   = reinterpret_cast<DType*>(storage->a);

  // Calc ndnz for the destination
  size_t ndnz  = 0;
  size_t i,j; // indexes of destination matrix
  size_t k,l; // indexes of source matrix
  for (i = 0; i < shape[0]; i++) {
    k = i + offset[0];
    for (j = 0; j < shape[1]; j++) {
      l = j + offset[1];

      if (j == i)  continue;

      if (k == l) { // for diagonal element of source
        if (src_a[k] != 0) ++ndnz;
      } else { // for non-diagonal element
        for (size_t c = src_ija[k]; c < src_ija[k+1]; c++) {
          if (src_ija[c] == l) {
            ++ndnz;
            break;
          }
        }
      }

    }
  }

  size_t request_capacity = shape[0] + ndnz + 1;
  //fprintf(stderr, "yale get copy: shape0=%d, shape1=%d, ndnz=%d, request_capacity=%d\n", shape[0], shape[1], ndnz, request_capacity);
  YALE_STORAGE* ns = nm_yale_storage_create(storage->dtype, shape, 2, request_capacity, storage->itype);

  if (ns->capacity < request_capacity)
    rb_raise(nm_eStorageTypeError, "conversion failed; capacity of %ld requested, max allowable is %ld", request_capacity, ns->capacity);

   // Initialize the A and IJA arrays
  init<DType,IType>(ns);
  IType* dst_ija = reinterpret_cast<IType*>(ns->ija);
  DType* dst_a   = reinterpret_cast<DType*>(ns->a);
 
  size_t ija = shape[0] + 1;
  DType val = src_a[storage->shape[0]]; // use 0 as the default for copy
  for (i = 0; i < shape[0]; ++i) {
    k = i + offset[0];
    for (j = 0; j < shape[1]; ++j) {
      bool found = false;
      l = j + offset[1];
    
      // Get value from source matrix
      if (k == l) { // source diagonal
        if (src_a[k] != 0) { // don't bother copying non-zero values from the diagonal
          val = src_a[k];
          found = true;
        }
      } else {
        // copy one non-diagonal element
        for (size_t c = src_ija[k]; !found && c < src_ija[k+1]; ++c) {
          if (src_ija[c] == l) {
            val   = src_a[c];
            found = true;
          }
        }
      }

      if (found) {
        // Set value in destination matrix
        if (i == j) {
          dst_a[i] = val;
        } else {
          // copy non-diagonal element
          dst_ija[ija] = j;
          dst_a[ija]   = val;

          ++ija;
          for (size_t c = i + 1; c <= shape[0]; ++c) {
            dst_ija[c] = ija;
          }
        }
      }
    }
  }

  dst_ija[shape[0]] = ija; // indicate the end of the last row
  ns->ndnz = ndnz;
  return ns;
}
/*
 * Returns a pointer to the correct location in the A vector of a YALE_STORAGE object, given some set of coordinates
 * (the coordinates are stored in slice).
 */
template <typename DType,typename IType>
void* ref(YALE_STORAGE* storage, SLICE* slice) {
  size_t* coords = slice->coords;

  if (!slice->single) rb_raise(rb_eNotImpError, "This type slicing not supported yet.");

  DType* a = reinterpret_cast<DType*>(storage->a);
  IType* ija = reinterpret_cast<IType*>(storage->ija);

  if (coords[0] == coords[1])
    return &(a[ coords[0] ]); // return diagonal entry

  if (ija[coords[0]] == ija[coords[0]+1])
    return &(a[ storage->shape[0] ]); // return zero pointer

	// binary search for the column's location
  int pos = binary_search<IType>(storage,
                                          ija[coords[0]],
                                          ija[coords[0]+1]-1,
                                          coords[1]);

  if (pos != -1 && ija[pos] == coords[1])
    return &(a[pos]); // found exact value

  return &(a[ storage->shape[0] ]); // return a pointer that happens to be zero
}

/*
 * Attempt to set some cell in a YALE_STORAGE object. Must supply coordinates and a pointer to a value (which will be
 * copied into the storage object).
 */
template <typename DType, typename IType>
char set(YALE_STORAGE* storage, SLICE* slice, void* value) {
  DType* v = reinterpret_cast<DType*>(value);
  size_t* coords = slice->coords;

  bool found = false;
  char ins_type;

  if (coords[0] == coords[1]) {
    reinterpret_cast<DType*>(storage->a)[coords[0]] = *v; // set diagonal
    return 'r';
  }

  // Get IJA positions of the beginning and end of the row
  if (reinterpret_cast<IType*>(storage->ija)[coords[0]] == reinterpret_cast<IType*>(storage->ija)[coords[0]+1]) {
  	// empty row
    ins_type = vector_insert<DType,IType>(storage, reinterpret_cast<IType*>(storage->ija)[coords[0]], &(coords[1]), v, 1, false);
    increment_ia_after<IType>(storage, storage->shape[0], coords[0], 1);
    storage->ndnz++;
    
    return ins_type;
  }

  // non-empty row. search for coords[1] in the IJA array, between ija and ija_next
  // (including ija, not including ija_next)
  //ija_size = get_size<IType>(storage);

  // Do a binary search for the column
  size_t pos = insert_search<IType>(storage,
                                    reinterpret_cast<IType*>(storage->ija)[coords[0]],
                                    reinterpret_cast<IType*>(storage->ija)[coords[0]+1]-1,
                                    coords[1], &found);

  if (found) { // replace
    reinterpret_cast<IType*>(storage->ija)[pos] = coords[1];
    reinterpret_cast<DType*>(storage->a)[pos]   = *v;
  	return 'r';
  }

  ins_type = vector_insert<DType,IType>(storage, pos, &(coords[1]), v, 1, false);
  increment_ia_after<IType>(storage, storage->shape[0], coords[0], 1);
  storage->ndnz++;

  return ins_type;
}

///////////
// Tests //
///////////

/*
 * Yale eql? -- for whole-matrix comparison returning a single value.
 */
template <typename LDType, typename RDType, typename IType>
static bool eqeq(const YALE_STORAGE* left, const YALE_STORAGE* right) {
  LDType* la = reinterpret_cast<LDType*>(left->a);
  RDType* ra = reinterpret_cast<RDType*>(right->a);

  // Compare the diagonals first.
  for (size_t index = 0; index < left->shape[0]; ++index) {
    if (la[index] != ra[index]) return false;
  }

  IType* lij = reinterpret_cast<IType*>(left->ija);
  IType* rij = reinterpret_cast<IType*>(right->ija);

  for (IType i = 0; i < left->shape[0]; ++i) {

  // Get start and end positions of row
    IType l_ija = lij[i],
          l_ija_next = lij[i+1],
          r_ija = rij[i],
          r_ija_next = rij[i+1];

    // Check to see if one row is empty and the other isn't.
    if (ndrow_is_empty<LDType,IType>(left, l_ija, l_ija_next)) {
      if (!ndrow_is_empty<RDType,IType>(right, r_ija, r_ija_next)) {
      	return false;
      }

    } else if (ndrow_is_empty<RDType,IType>(right, r_ija, r_ija_next)) {
    	// one is empty but the other isn't
      return false;

    } else if (!ndrow_eqeq_ndrow<LDType,RDType,IType>(left, right, l_ija, l_ija_next, r_ija, r_ija_next)) {
    	// Neither row is empty. Must compare the rows directly.
      return false;
    }

  }

  return true;
}

/*
 * Are two non-diagonal rows the same? We already know.
 */
template <typename LDType, typename RDType, typename IType>
static bool ndrow_eqeq_ndrow(const YALE_STORAGE* l, const YALE_STORAGE* r, IType l_ija, const IType l_ija_next, IType r_ija, const IType r_ija_next) {
  bool l_no_more = false, r_no_more = false;

  IType *lij = reinterpret_cast<IType*>(l->ija),
        *rij = reinterpret_cast<IType*>(r->ija);

  LDType* la = reinterpret_cast<LDType*>(l->a);
  RDType* ra = reinterpret_cast<RDType*>(r->a);

  IType l_ja = lij[l_ija],
        r_ja = rij[r_ija];
        
  IType ja = std::min(l_ja, r_ja);

  while (!(l_no_more && r_no_more)) {
    if (l_ja == r_ja) {

      if (ra[r_ija] != la[l_ija]) return false; // Direct comparison

      ++l_ija;
      ++r_ija;

      if (l_ija < l_ija_next) {
      	l_ja = lij[l_ija];

      } else {
      	l_no_more = true;
      }

      if (r_ija < r_ija_next) {
      	r_ja = rij[r_ija];

      } else {
      	r_no_more = true;
      }

      ja = std::min(l_ja, r_ja);

    } else if (l_no_more || ja < l_ja) {

      if (ra[r_ija] != 0) return false;

      ++r_ija;
      if (r_ija < r_ija_next) {
      	// get next column
      	r_ja = rij[r_ija];
        ja = std::min(l_ja, r_ja);

      } else {
      	l_no_more = true;
      }

    } else if (r_no_more || ja < r_ja) {

      if (la[l_ija] != 0) return false;

      ++l_ija;
      if (l_ija < l_ija_next) {
      	// get next column
        l_ja = lij[l_ija];
        ja = std::min(l_ja, r_ja);
      } else {
      	l_no_more = true;
      }

    } else {
      std::fprintf(stderr, "Unhandled in eqeq: l_ja=%d, r_ja=%d\n", (int)l_ja, (int)r_ja);
    }
  }

	// every item matched
  return true;
}

/*
 * Is the non-diagonal portion of the row empty?
 */
template <typename DType, typename IType>
static bool ndrow_is_empty(const YALE_STORAGE* s, IType ija, const IType ija_next) {
  if (ija == ija_next) return true;

  DType* a = reinterpret_cast<DType*>(s->a);

	// do all the entries = zero?
  for (; ija < ija_next; ++ija) {
    if (a[ija] != 0) return false;
  }

  return true;
}

//////////
// Math //
//////////

#define YALE_IA(s) (reinterpret_cast<IType*>(s->ija))
#define YALE_IJ(s) (reinterpret_cast<IType*>(s->ija) + s->shape[0] + 1)
#define YALE_COUNT(yale) (yale->ndnz + yale->shape[0])

template <typename nm::ewop_t op, typename IType, typename DType>
YALE_STORAGE* ew_op(const YALE_STORAGE* left, const YALE_STORAGE* right, dtype_t dtype) {
	size_t  init_capacity;
	size_t* new_shape;
	
	unsigned int	da_index,
								la_index,
								ra_index,
								
								a_index_offset,
								
								la_row_max,
								ra_row_max,
								
								row_index;
	
	DType tmp_result;
	
	DType * la = reinterpret_cast<DType*> (left->a),
				* ra = reinterpret_cast<DType*>(right->a),
				* da;
	
	YALE_STORAGE* dest;
	
	new_shape			= reinterpret_cast<size_t*>(calloc(2, sizeof(size_t)));
	new_shape[0]	= left->shape[0];
	new_shape[1]	= left->shape[1];
	
	init_capacity = std::min(left->ndnz + right->ndnz + new_shape[0], new_shape[0] * new_shape[1]);
	
	dest	= nm_yale_storage_create(dtype, new_shape, 2, init_capacity, left->itype);
	da		= reinterpret_cast<DType*>(dest->a);
	
	// Calculate diagonal values.
	for (da_index = 0; da_index < dest->shape[0]; ++da_index) {
		da[da_index] = ew_op_switch<op, DType, DType>(la[da_index], ra[da_index]);
	}
	
	// Set the zero representation seperator.
	da[da_index] = typeid(DType) == typeid(RubyObject) ? INT2FIX(0) : 0;
	
	/*
	 * Calculate the offset between start of the A arrays and the non-diagonal
	 * entries.
	 */
	a_index_offset = dest->shape[0] + 1;
	
	// Re-base the A arrays.
	la = la + a_index_offset;
	ra = ra + a_index_offset;
	da = da + a_index_offset;
	
	// Initialize our A array indices.
	la_index = ra_index = da_index = 0;
	
	// Calculate the non-diagonal values.
	for (row_index = 0; row_index < dest->shape[0]; ++row_index) {
		/*
		 * Each row.
		 */
		
		printf("Row %d\n", row_index);
		
		// Get row bounds.
		la_row_max = YALE_IA( left)[row_index + 1] - a_index_offset;
		ra_row_max = YALE_IA(right)[row_index + 1] - a_index_offset;
		
		printf("Left  : Row Start: %d - Row End %d\n", la_index + a_index_offset, la_row_max + a_index_offset);
		printf("Right : Row Start: %d - Row End %d\n", ra_index + a_index_offset, ra_row_max + a_index_offset);
		
		/*
		 * Set this row's left bound (which is also the previous row's right
		 * bound).
		 */
		YALE_IA(dest)[row_index] = da_index + a_index_offset;
		
		printf("Left bound of row %d in destination: %d\n", (int)row_index, (int)YALE_IA(dest)[row_index]);
		
		// Iterate over non-diagonal entries in this row.
		while (la_index < la_row_max and ra_index < ra_row_max) {
			/*
			 * Elements are present on both the left- and right-hand side.
			 */
			
			printf("Marker 0\n");
			
			if (YALE_IJ(left)[la_index] == YALE_IJ(right)[ra_index]) {
				/*
				 * Current left- and right-hand values are in the same row and
				 * column.
				 */
				
				printf("Calculating value for [%d, %d].\n", (int)row_index, (int)YALE_IJ(left)[la_index]);
				
				tmp_result = ew_op_switch<op, DType, DType>(la[la_index], ra[ra_index]);
				
				if (tmp_result != 0) {
					printf("Setting value for [%d, %d] at index %d in destination's A array.\n", (int)row_index, (int)YALE_IJ(left)[la_index], (int)(da_index + a_index_offset));
					
					da[da_index]						= tmp_result;
					YALE_IJ(dest)[da_index] = YALE_IJ(left)[la_index];
					
					++da_index;
					
				} else {
					printf("Result was 0.  Skipping.\n");
				}
				
				++la_index;
				++ra_index;
				
			} else if (YALE_IJ(left)[la_index] < YALE_IJ(right)[ra_index]) {
				/*
				 * The right-hand index is ahead of the left-hand index.
				 */
				
				if (op != EW_MUL) {
					// If this is multiplion there is no point in doing the operation.
					
					tmp_result = ew_op_switch<op, DType, DType>(la[la_index], typeid(DType) == typeid(RubyObject) ? INT2FIX(0) : 0);
				
					printf("Setting value for [%d, %d].\n", (int)row_index, (int)YALE_IJ(left)[la_index]);
				
					if (tmp_result != 0) {
						da[da_index]						= tmp_result;
						YALE_IJ(dest)[da_index] = YALE_IJ(left)[la_index];
				
						++da_index;
					}
				}
				
				++la_index;
				
			} else {
				/*
				 * The left-hand index is ahead of the right-hand index.
				 */
				
				if (op != EW_MUL) {
					// If this is multiplion there is no point in doing the operation.
					
					tmp_result = ew_op_switch<op, DType, DType>(typeid(DType) == typeid(RubyObject) ? INT2FIX(0) : 0, ra[ra_index]);
				
					printf("Setting value for [%d, %d].\n", (int)row_index, (int)YALE_IJ(right)[ra_index]);
				
					if (tmp_result != 0) {
						da[da_index]						= tmp_result;
						YALE_IJ(dest)[da_index] = YALE_IJ(right)[ra_index];
				
						++da_index;
					}
				}
				
				++ra_index;
			}
		}
		
		if (op != EW_MUL) {
			/*
			 * Process the remaining elements on the left- or right-hand side.  One or
			 * the other, or neither, of the following loops may execute, but not
			 * both.
			 *
			 * If we are doing multiplication this is unnecessary as all remaining
			 * operations will produce a zero value.
			 */
		
			while (la_index < la_row_max) {
				/*
				 * Process the remaining elements on the left-hand side.
				 */
				
				printf("Marker 1\n");
				
				tmp_result = ew_op_switch<op, DType, DType>(la[la_index], typeid(DType) == typeid(RubyObject) ? INT2FIX(0) : 0);
				
				printf("Setting value for [%d, %d].\n", (int)row_index, (int)YALE_IJ(left)[la_index]);
				
				if (tmp_result != 0) {
					da[da_index]						= tmp_result;
					YALE_IJ(dest)[da_index] = YALE_IJ(left)[la_index];
					
					++da_index;
				}
				
				++la_index;
			}
		
			while (ra_index < ra_row_max) {
				/*
				 * Process the remaining elements on the right-hand side.
				 */
				
				printf("Marker 2\n");
				
				tmp_result = ew_op_switch<op, DType, DType>(typeid(DType) == typeid(RubyObject) ? INT2FIX(0) : 0, ra[ra_index]);
				
				printf("Setting value for [%d, %d].\n", (int)row_index, (int)YALE_IJ(right)[ra_index]);
				
				if (tmp_result != 0) {
					da[da_index]						= tmp_result;
					YALE_IJ(dest)[da_index] = YALE_IJ(right)[ra_index];
					
					++da_index;
				}
				
				++ra_index;
			}
		}
		
		// Advance the row indices.
		la_index = la_row_max;
		ra_index = ra_row_max;
		
		printf("End of row %d\n\n", row_index);
	}
	
	// Set the last row's right bound.
	YALE_IA(dest)[row_index] = da_index + a_index_offset;
	
	printf("Right bound of row %d in destination: %d\n", row_index - 1, da_index + a_index_offset);
	
	// Set the number of non-diagonal non-zero entries in the destination matrix.
	dest->ndnz = da_index;
	
	printf("Number of non-diagonal non-zero entries: %ld\n\n", (unsigned long)(dest->ndnz));
	
	// Set the capacity of the destination matrix.
	dest->capacity = dest->shape[0] + dest->ndnz + 1;
	
	// Resize the destination matrix.
	dest->a		= realloc(dest->a,   sizeof(DType) * dest->capacity);
	dest->ija = realloc(dest->ija, sizeof(IType) * dest->capacity);
	
	return dest;
}

/////////////
// Utility //
/////////////

/*
 * Binary search for returning stored values. Returns a non-negative position, or -1 for not found.
 */
template <typename IType>
int binary_search(YALE_STORAGE* s, IType left, IType right, IType key) {

  if (left > right) return -1;

  IType* ija = reinterpret_cast<IType*>(s->ija);

  IType mid = (left + right)/2;
  IType mid_j = ija[mid];

  if (mid_j == key)
  	return mid;

  else if (mid_j > key)
  	return binary_search<IType>(s, left, mid - 1, key);

  else
  	return binary_search<IType>(s, mid + 1, right, key);
}



/*
 * Resize yale storage vectors A and IJA in preparation for an insertion.
 */
template <typename DType, typename IType>
static char vector_insert_resize(YALE_STORAGE* s, size_t current_size, size_t pos, size_t* j, size_t n, bool struct_only) {
  // Determine the new capacity for the IJA and A vectors.
  size_t new_capacity = s->capacity * GROWTH_CONSTANT;
  size_t max_capacity = max_size(s);

  if (new_capacity > max_capacity) {
    new_capacity = max_capacity;

    if (current_size + n > max_capacity) rb_raise(rb_eNoMemError, "insertion size exceeded maximum yale matrix size");
  }

  if (new_capacity < current_size + n)
  	new_capacity = current_size + n;

  // Allocate the new vectors.
  IType* new_ija     = ALLOC_N( IType, new_capacity );
  NM_CHECK_ALLOC(new_ija);

  DType* new_a       = ALLOC_N( DType, new_capacity );
  NM_CHECK_ALLOC(new_a);

  IType* old_ija     = reinterpret_cast<IType*>(s->ija);
  DType* old_a       = reinterpret_cast<DType*>(s->a);

  // Copy all values prior to the insertion site to the new IJA and new A
  if (struct_only) {
    for (size_t i = 0; i < pos; ++i) {
      new_ija[i] = old_ija[i];
    }
  } else {
    for (size_t i = 0; i < pos; ++i) {
      new_ija[i] = old_ija[i];
      new_a[i]   = old_a[i];
    }
  }
  
  // Copy all values subsequent to the insertion site to the new IJA and new A, leaving room (size n) for insertion.
  if (struct_only) {
    for (size_t i = pos; i < current_size; ++i) {
      new_ija[i+n] = old_ija[i];
    }
  } else {
    for (size_t i = pos; i < current_size; ++i) {
      new_ija[i+n] = old_ija[i];
      new_a[i+n] = old_a[i];
    }
  }

  s->capacity = new_capacity;

  free(s->ija);
  free(s->a);

  s->ija = reinterpret_cast<void*>(new_ija);
  s->a   = reinterpret_cast<void*>(new_a);

  fprintf(stderr, "resize\n");

  return 'i';
}

/*
 * Insert a value or contiguous values in the ija and a vectors (after ja and
 * diag). Does not free anything; you are responsible!
 *
 * TODO: Improve this so it can handle non-contiguous element insertions
 *	efficiently. For now, we can just sort the elements in the row in
 *	question.)
 */
template <typename DType, typename IType>
static char vector_insert(YALE_STORAGE* s, size_t pos, size_t* j, void* val_, size_t n, bool struct_only) {
  if (pos < s->shape[0]) {
    rb_raise(rb_eArgError, "vector insert pos (%d) is before beginning of ja (%d); this should not happen", pos, s->shape[0]);
  }

  DType* val = reinterpret_cast<DType*>(val_);

  size_t size = get_size<IType>(s);

  IType* ija = reinterpret_cast<IType*>(s->ija);
  DType* a   = reinterpret_cast<DType*>(s->a);

  if (size + n > s->capacity) {
  	vector_insert_resize<DType,IType>(s, size, pos, j, n, struct_only);

    // Need to get the new locations for ija and a.
  	ija = reinterpret_cast<IType*>(s->ija);
    a   = reinterpret_cast<DType*>(s->a);

  } else {
    /*
     * No resize required:
     * easy (but somewhat slow), just copy elements to the tail, starting at
     * the end, one element at a time.
     *
     * TODO: This can be made slightly more efficient, but only after the tests
     *	are written.
     */

    if (struct_only) {
      for (size_t i = 0; i < size - pos; ++i) {
        ija[size+n-1-i] = ija[size-1-i];
      }
    } else {
      for (size_t i = 0; i < size - pos; ++i) {
        ija[size+n-1-i] = ija[size-1-i];
        a[size+n-1-i]   = a[size-1-i];
      }
    }
  }

  // Now insert the new values.
  if (struct_only) {
    for (size_t i = 0; i < n; ++i) {
      ija[pos+i]  = j[i];
    }
  } else {
    for (size_t i = 0; i < n; ++i) {
      ija[pos+i]  = j[i];
      a[pos+i]    = val[i];
    }
  }

  return 'i';
}

/*
 * If we add n items to row i, we need to increment ija[i+1] and onward.
 */
template <typename IType>
static void increment_ia_after(YALE_STORAGE* s, IType ija_size, IType i, IType n) {
  IType* ija = reinterpret_cast<IType*>(s->ija);

  ++i;
  for (; i <= ija_size; ++i) {
    ija[i] += n;
  }
}

/*
 * Binary search for returning insertion points.
 */
template <typename IType>
static IType insert_search(YALE_STORAGE* s, IType left, IType right, IType key, bool* found) {

  if (left > right) {
    *found = false;
    return left;
  }

  IType* ija = reinterpret_cast<IType*>(s->ija);
  IType mid = (left + right)/2;
  IType mid_j = ija[mid];

  if (mid_j == key) {
    *found = true;
    return mid;

  } else if (mid_j > key) {
  	return insert_search<IType>(s, left, mid-1, key, found);

  } else {
  	return insert_search<IType>(s, mid+1, right, key, found);
  }
}

/////////////////////////
// Copying and Casting //
/////////////////////////

/*
 * Templated copy constructor for changing dtypes.
 */
template <typename LDType, typename RDType, typename IType>
YALE_STORAGE* cast_copy(const YALE_STORAGE* rhs, dtype_t new_dtype) {

  // Allocate a new structure
  size_t size = get_size<IType>(rhs);
  YALE_STORAGE* lhs = copy_alloc_struct<IType>(rhs, new_dtype, rhs->capacity, size);

  if (rhs->dtype == new_dtype) {  // FIXME: Test if this condition is actually faster; second condition should work just as well.

    memcpy(lhs->a, rhs->a, size * DTYPE_SIZES[new_dtype]);

  } else {

    LDType* la = reinterpret_cast<LDType*>(lhs->a);
    RDType* ra = reinterpret_cast<RDType*>(rhs->a);

    for (size_t index = 0; index < size; ++index) {
      la[index] = ra[index];
    }

  }

  return lhs;
}

/*
 * Template access for getting the size of Yale storage.
 */
template <typename IType>
static inline size_t get_size(const YALE_STORAGE* storage) {
  return static_cast<size_t>(reinterpret_cast<IType*>(storage->ija)[ storage->shape[0] ]);
}


/*
 * Allocate for a copy or copy-cast operation, and copy the IJA portion of the
 * matrix (the structure).
 */
template <typename IType>
static YALE_STORAGE* copy_alloc_struct(const YALE_STORAGE* rhs, const dtype_t new_dtype, const size_t new_capacity, const size_t new_size) {
  YALE_STORAGE* lhs = ALLOC( YALE_STORAGE );
  lhs->dim          = rhs->dim;
  lhs->shape        = ALLOC_N( size_t, lhs->dim );
  memcpy(lhs->shape, rhs->shape, lhs->dim * sizeof(size_t));
  lhs->itype        = rhs->itype;
  lhs->capacity     = new_capacity;
  lhs->dtype        = new_dtype;
  lhs->ndnz         = rhs->ndnz;

  lhs->ija          = ALLOC_N( IType, lhs->capacity );
  lhs->a            = ALLOC_N( char, DTYPE_SIZES[new_dtype] * lhs->capacity );

  // Now copy the contents -- but only within the boundaries set by the size. Leave
  // the rest uninitialized.
  for (size_t i = 0; i < get_size<IType>(rhs); ++i)
    reinterpret_cast<IType*>(lhs->ija)[i] = reinterpret_cast<IType*>(rhs->ija)[i]; // copy indices

  return lhs;
}

template <typename DType, typename IType>
static STORAGE* matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector, nm::itype_t result_itype) {
  YALE_STORAGE *left  = (YALE_STORAGE*)(casted_storage.left),
               *right = (YALE_STORAGE*)(casted_storage.right);

  // We can safely get dtype from the casted matrices; post-condition of binary_storage_cast_alloc is that dtype is the
  // same for left and right.
  // int8_t dtype = left->dtype;

  // Create result storage.
  YALE_STORAGE* result = nm_yale_storage_create(left->dtype, resulting_shape, 2, left->capacity + right->capacity, result_itype);
  init<DType,IType>(result);

  // Massage the IType arrays into the correct form.

  IType* ijl;
  if (left->itype == result_itype) ijl = reinterpret_cast<IType*>(left->ija);
  else {  // make a temporary copy of the IJA vector for L with the correct itype
    std::cerr << "changing left itype from " << static_cast<uint8_t>(left->itype) << " to " << static_cast<int8_t>(result_itype) << std::endl;
    size_t length = nm_yale_storage_get_size(left);
    std::cerr << "length = " << length << std::endl;
    ijl = ALLOCA_N(IType, length);
    copy_recast_itype_vector(reinterpret_cast<void*>(left->ija), left->itype, reinterpret_cast<void*>(ijl), result_itype, length);
  }

  IType* ijr;
  if (right->itype == result_itype) ijr = reinterpret_cast<IType*>(right->ija);
  else {  // make a temporary copy of the IJA vector for R with the correct itype
    std::cerr << "changing right itype from " << static_cast<uint8_t>(right->itype) << " to " << static_cast<int8_t>(result_itype) << std::endl;
    size_t length = nm_yale_storage_get_size(right);
    std::cerr << "length = " << length << std::endl;
    ijr = ALLOCA_N(IType, length);
    copy_recast_itype_vector(reinterpret_cast<void*>(right->ija), right->itype, reinterpret_cast<void*>(ijr), result_itype, length);
  }

  IType* ija = reinterpret_cast<IType*>(result->ija);

  // Symbolic multiplication step (build the structure)
  nm::math::symbmm<IType>(result->shape[0], result->shape[1], ijl, ijl, true, ijr, ijr, true, ija, true);

  // Numeric multiplication step (fill in the elements)
  nm::math::numbmm<DType,IType>(result->shape[0], result->shape[1],
                                ijl, ijl, reinterpret_cast<DType*>(left->a), true,
                                ijr, ijr, reinterpret_cast<DType*>(right->a), true,
                                ija, ija, reinterpret_cast<DType*>(result->a), true);

  // Sort the columns
  nm::math::smmp_sort_columns<DType,IType>(result->shape[0], ija, ija, reinterpret_cast<DType*>(result->a));

  return reinterpret_cast<STORAGE*>(result);
}


} // end of namespace nm::yale_storage


// Helper function used only for the RETURN_SIZED_ENUMERATOR macro. Returns the length of
// the matrix's storage.
static VALUE nm_yale_enumerator_length(VALUE nmatrix) {
  long len = nm_yale_storage_get_size(NM_STORAGE_YALE(nmatrix));
  return LONG2NUM(len);
}


template <typename DType, typename IType>
struct yale_each_stored_with_indices_helper {
  static VALUE iterate(VALUE nm) {

    YALE_STORAGE* s = NM_STORAGE_YALE(nm);
    DType* a    = reinterpret_cast<DType*>(s->a);
    IType* ija  = reinterpret_cast<IType*>(s->ija);

    // If we don't have a block, return an enumerator.
    RETURN_SIZED_ENUMERATOR(nm, 0, 0, nm_yale_enumerator_length);

    // Iterate along diagonal
    for (size_t k = 0; k < s->shape[0]; ++k) {
      VALUE ii = LONG2NUM(k),
            jj = LONG2NUM(k);

      VALUE v = rubyobj_from_cval(&(a[k]), NM_DTYPE(nm)).rval;
      rb_yield_values(3, v, ii, jj );
    }

    // Iterate through non-diagonal elements, row by row
    for (long i = 0; i < s->shape[0]; ++i) {
      long p      = static_cast<long>( ija[i]   ),
           next_p = static_cast<long>( ija[i+1] );

      for (; p < next_p; ++p) {
        long j = static_cast<long>(ija[p]);
        VALUE ii = LONG2NUM(i),
              jj = LONG2NUM(j);

        VALUE v = rubyobj_from_cval(&(a[p]), NM_DTYPE(nm)).rval;
        rb_yield_values(3, v, ii, jj );
      }
    }

    return nm;
  }
};


template <typename IType>
struct yale_each_stored_with_indices_helper<RubyObject, IType> {
  static VALUE iterate(VALUE nm) {

    YALE_STORAGE* s = NM_STORAGE_YALE(nm);
    RubyObject* a   = reinterpret_cast<RubyObject*>(s->a);
    IType* ija      = reinterpret_cast<IType*>(s->ija);

    // If we don't have a block, return an enumerator.
    RETURN_SIZED_ENUMERATOR(nm, 0, 0, nm_yale_enumerator_length);

    // Iterate along diagonal
    for (size_t k = 0; k < s->shape[0]; ++k) {
      VALUE ii = LONG2NUM(k),
            jj = LONG2NUM(k);
      rb_yield_values(3, a[k].rval, ii, jj );  // yield element, i, j
    }

    // Iterate through non-diagonal elements, row by row
    for (long i = 0; i < s->shape[0]; ++i) {
      long p      = static_cast<long>( ija[i]   ),
           next_p = static_cast<long>( ija[i+1] );

      for (; p < next_p; ++p) {
        long j = static_cast<long>(ija[p]);
        VALUE ii = LONG2NUM(i),
              jj = LONG2NUM(j);

        rb_yield_values(3, a[p].rval, ii, jj );
      }
    }

    return nm;
  }
};


/*
 * This function and the two helper structs enable us to use partial template specialization.
 * See also: http://stackoverflow.com/questions/6623375/c-template-specialization-on-functions
 */
template <typename DType, typename IType>
static VALUE yale_each_stored_with_indices(VALUE nm) {
  return yale_each_stored_with_indices_helper<DType, IType>::iterate(nm);
}


} // end of namespace nm.

///////////////////
// Ruby Bindings //
///////////////////

/* These bindings are mostly only for debugging Yale. They are called from Init_nmatrix. */

extern "C" {

void nm_init_yale_functions() {
	/*
	 * This module stores methods that are useful for debugging Yale matrices,
	 * i.e. the ones with +:yale+ stype.	
	 */
  cNMatrix_YaleFunctions = rb_define_module_under(cNMatrix, "YaleFunctions");

  rb_define_method(cNMatrix_YaleFunctions, "yale_ija", (METHOD)nm_ija, -1);
  rb_define_method(cNMatrix_YaleFunctions, "yale_a", (METHOD)nm_a, -1);
  rb_define_method(cNMatrix_YaleFunctions, "yale_size", (METHOD)nm_size, 0);
  rb_define_method(cNMatrix_YaleFunctions, "yale_ia", (METHOD)nm_ia, 0);
  rb_define_method(cNMatrix_YaleFunctions, "yale_ja", (METHOD)nm_ja, 0);
  rb_define_method(cNMatrix_YaleFunctions, "yale_d", (METHOD)nm_d, -1);
  rb_define_method(cNMatrix_YaleFunctions, "yale_lu", (METHOD)nm_lu, 0);

  rb_define_method(cNMatrix_YaleFunctions, "yale_nd_row", (METHOD)nm_nd_row, -1);
  rb_define_method(cNMatrix_YaleFunctions, "yale_vector_insert", (METHOD)nm_vector_insert, -1);

  rb_define_const(cNMatrix_YaleFunctions, "YALE_GROWTH_CONSTANT", rb_float_new(nm::yale_storage::GROWTH_CONSTANT));
}


/////////////////
// C ACCESSORS //
/////////////////



VALUE nm_yale_each_stored_with_indices(VALUE nmatrix) {
  nm::dtype_t d = NM_DTYPE(nmatrix);
  nm::itype_t i = NM_ITYPE(nmatrix);

  NAMED_LI_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_each_stored_with_indices, VALUE, VALUE)

  return ttable[d][i](nmatrix);
}


/*
 * C accessor for inserting some value in a matrix (or replacing an existing cell).
 */
char nm_yale_storage_set(STORAGE* storage, SLICE* slice, void* v) {
  NAMED_LI_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::set, char, YALE_STORAGE* storage, SLICE* slice, void* value);

  YALE_STORAGE* casted_storage = (YALE_STORAGE*)storage;

  return ttable[casted_storage->dtype][casted_storage->itype](casted_storage, slice, v);
}

/*
 * C accessor for yale_storage::get, which returns a slice of YALE_STORAGE object by coppy
 *
 * Slicing-related.
 */
void* nm_yale_storage_get(STORAGE* storage, SLICE* slice) {
  NAMED_LI_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::get, void*, YALE_STORAGE* storage, SLICE* slice);
  YALE_STORAGE* s = (YALE_STORAGE*)storage;


  YALE_STORAGE* casted_storage = (YALE_STORAGE*)storage;
  return ttable[casted_storage->dtype][casted_storage->itype](casted_storage, slice);
}

/*
 * C accessor for yale_storage::vector_insert
 */
static char nm_yale_storage_vector_insert(YALE_STORAGE* s, size_t pos, size_t* js, void* vals, size_t n, bool struct_only, nm::dtype_t dtype, nm::itype_t itype) {
  NAMED_LI_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::vector_insert, char, YALE_STORAGE*, size_t, size_t*, void*, size_t, bool);

  return ttable[dtype][itype](s, pos, js, vals, n, struct_only);
}

/*
 * C accessor for yale_storage::increment_ia_after, typically called after ::vector_insert
 */
static void nm_yale_storage_increment_ia_after(YALE_STORAGE* s, size_t ija_size, size_t i, size_t n, nm::itype_t itype) {
  NAMED_ITYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::c_increment_ia_after, void, YALE_STORAGE*, size_t, size_t, size_t);

  ttable[itype](s, ija_size, i, n);
}


/*
 * C accessor for yale_storage::ref, which returns a pointer to the correct location in a YALE_STORAGE object
 * for some set of coordinates.
 */
void* nm_yale_storage_ref(STORAGE* storage, SLICE* slice) {
  NAMED_LI_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::ref, void*, YALE_STORAGE* storage, SLICE* slice);

  YALE_STORAGE* casted_storage = (YALE_STORAGE*)storage;
  return ttable[casted_storage->dtype][casted_storage->itype](casted_storage, slice);
}

/*
 * C accessor for determining whether two YALE_STORAGE objects have the same contents.
 *
 * FIXME: Is this for element-wise or whole-matrix equality?
 */
bool nm_yale_storage_eqeq(const STORAGE* left, const STORAGE* right) {
  NAMED_LRI_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::eqeq, bool, const YALE_STORAGE* left, const YALE_STORAGE* right);

  const YALE_STORAGE* casted_left = reinterpret_cast<const YALE_STORAGE*>(left);

  return ttable[casted_left->dtype][right->dtype][casted_left->itype](casted_left, (const YALE_STORAGE*)right);
}

/*
 * Copy constructor for changing dtypes. (C accessor)
 */
STORAGE* nm_yale_storage_cast_copy(const STORAGE* rhs, nm::dtype_t new_dtype) {
  NAMED_LRI_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::cast_copy, YALE_STORAGE*, const YALE_STORAGE* rhs, nm::dtype_t new_dtype);

  const YALE_STORAGE* casted_rhs = reinterpret_cast<const YALE_STORAGE*>(rhs);

  return (STORAGE*)ttable[new_dtype][casted_rhs->dtype][casted_rhs->itype](casted_rhs, new_dtype);
}

/*
 * Returns size of Yale storage as a size_t (no matter what the itype is). (C accessor)
 */
size_t nm_yale_storage_get_size(const YALE_STORAGE* storage) {
  NAMED_ITYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::get_size, size_t, const YALE_STORAGE* storage);

  return ttable[storage->itype](storage);
}

/*
 * C accessor for allocating a yale storage object for cast-copying. Copies the IJA vector, does not copy the A vector.
 */
static YALE_STORAGE* nm_copy_alloc_struct(const YALE_STORAGE* rhs, const nm::dtype_t new_dtype, const size_t new_capacity, const size_t new_size) {
  NAMED_ITYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::copy_alloc_struct, YALE_STORAGE*, const YALE_STORAGE* rhs, const nm::dtype_t new_dtype, const size_t new_capacity, const size_t new_size);

  return ttable[rhs->itype](rhs, new_dtype, new_capacity, new_size);
}

/*
 * Transposing copy constructor.
 */
STORAGE* nm_yale_storage_copy_transposed(const STORAGE* rhs_base) {
  YALE_STORAGE* rhs = (YALE_STORAGE*)rhs_base;

  size_t* shape = ALLOC_N(size_t, 2);
  shape[0] = rhs->shape[1];
  shape[1] = rhs->shape[0];

  size_t size   = nm_yale_storage_get_size(rhs);

  YALE_STORAGE* lhs = nm_yale_storage_create(rhs->dtype, shape, 2, size, nm::UINT8);
  nm_yale_storage_init(lhs);

  NAMED_LI_DTYPE_TEMPLATE_TABLE(transp, nm::math::transpose_yale, void, const size_t n, const size_t m, const void* ia_, const void* ja_, const void* a_, const bool diaga, void* ib_, void* jb_, void* b_, const bool move);

  transp[lhs->dtype][lhs->itype](rhs->shape[0], rhs->shape[1], rhs->ija, rhs->ija, rhs->a, true, lhs->ija, lhs->ija, lhs->a, true);

  return (STORAGE*)lhs;
}

/*
 * C accessor for multiplying two YALE_STORAGE matrices, which have already been casted to the same dtype.
 *
 * FIXME: There should be some mathematical way to determine the worst-case IType based on the input ITypes. Right now
 * it just uses the default.
 */
STORAGE* nm_yale_storage_matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector) {
  LI_DTYPE_TEMPLATE_TABLE(nm::yale_storage::matrix_multiply, STORAGE*, const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector, nm::itype_t resulting_itype);

  YALE_STORAGE* left = reinterpret_cast<YALE_STORAGE*>(casted_storage.left);
  YALE_STORAGE* right = reinterpret_cast<YALE_STORAGE*>(casted_storage.right);

  // Determine the itype for the matrix that will be returned.
  nm::itype_t itype = nm_yale_storage_itype_by_shape(resulting_shape),
              max_itype = NM_MAX_ITYPE(left->itype, right->itype);
  if (static_cast<int8_t>(itype) < static_cast<int8_t>(max_itype)) itype = max_itype;

  return ttable[left->dtype][itype](casted_storage, resulting_shape, vector, itype);
}

/*
 * Documentation goes here.
 */
STORAGE* nm_yale_storage_ew_op(nm::ewop_t op, const STORAGE* left, const STORAGE* right, VALUE scalar) {
	OP_ITYPE_DTYPE_TEMPLATE_TABLE(nm::yale_storage::ew_op, YALE_STORAGE*, const YALE_STORAGE*, const YALE_STORAGE*, nm::dtype_t);
	
	YALE_STORAGE* new_l = NULL, * new_r = NULL;
	YALE_STORAGE* result;
	
	const YALE_STORAGE* casted_l, * casted_r;
	
	nm::dtype_t new_dtype;
	
	if (left->dtype != right->dtype) {
		
		new_dtype = Upcast[left->dtype][right->dtype];
		
		if (left->dtype != new_dtype) {
			new_l = reinterpret_cast<YALE_STORAGE*>(nm_yale_storage_cast_copy( left, new_dtype));
		}
		
		if (right->dtype != new_dtype) {
			new_r = reinterpret_cast<YALE_STORAGE*>(nm_yale_storage_cast_copy(right, new_dtype));
		}
		
		if (static_cast<uint8_t>(op) < nm::NUM_NONCOMP_EWOPS) {
			result = ttable[op][new_l->itype][new_dtype](	left->dtype  == new_dtype ?
																											reinterpret_cast<const YALE_STORAGE*>( left) :
																											reinterpret_cast<const YALE_STORAGE*>(new_l),
																										
																										right->dtype == new_dtype ?
																											reinterpret_cast<const YALE_STORAGE*>(right) :
																											reinterpret_cast<const YALE_STORAGE*>(new_r),
																										
																										new_dtype);
			
		} else {
			rb_raise(rb_eNotImpError, "Elementwise comparison is not yet implemented for the Yale storage class.");
		}
		
		if (new_l != NULL) {
			nm_yale_storage_delete(new_l);
		}
		
		if (new_r != NULL) {
			nm_yale_storage_delete(new_r);
		}
		
		return result;
		
	} else {
		
		casted_l = reinterpret_cast<const YALE_STORAGE*>( left);
		casted_r = reinterpret_cast<const YALE_STORAGE*>(right);
		
		if (static_cast<uint8_t>(op) < nm::NUM_NONCOMP_EWOPS) {
			
			return ttable[op][casted_l->itype][casted_l->dtype](casted_l, casted_r, casted_l->dtype);
		
		} else {
			rb_raise(rb_eNotImpError, "Elementwise comparison is not yet implemented for the Yale storage class.");
		}
	}
}

///////////////
// Lifecycle //
///////////////

/*
 * C accessor function for creating a YALE_STORAGE object. Prior to calling this function, you MUST
 * allocate shape (should be size_t * 2) -- don't use use a regular size_t array!
 *
 * For this type, dim must always be 2. The final argument is the initial capacity with which to
 * create the storage.
 */

YALE_STORAGE* nm_yale_storage_create(nm::dtype_t dtype, size_t* shape, size_t dim, size_t init_capacity, nm::itype_t min_itype) {
  YALE_STORAGE* s;
  size_t max_capacity;

	// FIXME: This error should be handled in the nmatrix.c file.
  if (dim != 2) {
   	rb_raise(rb_eNotImpError, "Can only support 2D matrices");
  }

  s = alloc(dtype, shape, dim, min_itype);
  max_capacity = nm::yale_storage::max_size(s);

  // Set matrix capacity (and ensure its validity)
  if (init_capacity < NM_YALE_MINIMUM(s)) {
  	s->capacity = NM_YALE_MINIMUM(s);

  } else if (init_capacity > max_capacity) {
  	// Don't allow storage to be created larger than necessary
  	s->capacity = max_capacity;

  } else {
  	s->capacity = init_capacity;

  }

  s->ija = ALLOC_N( char, ITYPE_SIZES[s->itype] * s->capacity );
  s->a   = ALLOC_N( char, DTYPE_SIZES[s->dtype] * s->capacity );

  return s;
}

/*
 * Destructor for yale storage (C-accessible).
 */
void nm_yale_storage_delete(STORAGE* s) {
  if (s) {
    YALE_STORAGE* storage = (YALE_STORAGE*)s;
    free(storage->shape);
    free(storage->ija);
    free(storage->a);
    free(storage);
  }
}

/*
 * C accessor for yale_storage::init, a templated function.
 *
 * Initializes the IJA vector of the YALE_STORAGE matrix.
 */
void nm_yale_storage_init(YALE_STORAGE* s) {
  NAMED_LI_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::init, void, YALE_STORAGE* s);

  ttable[s->dtype][s->itype](s);
}


/*
 * Ruby GC mark function for YALE_STORAGE. C accessible.
 */
void nm_yale_storage_mark(void* storage_base) {
  YALE_STORAGE* storage = (YALE_STORAGE*)storage_base;
  size_t i;

  if (storage && storage->dtype == nm::RUBYOBJ) {
  	for (i = storage->capacity; i-- > 0;) {
      rb_gc_mark(*((VALUE*)((char*)(storage->a) + i*DTYPE_SIZES[nm::RUBYOBJ])));
    }
  }
}

/*
 * Allocates and initializes the basic struct (but not the IJA or A vectors).
 */
static YALE_STORAGE* alloc(nm::dtype_t dtype, size_t* shape, size_t dim, nm::itype_t min_itype) {
  YALE_STORAGE* s;

  s = ALLOC( YALE_STORAGE );

  s->ndnz        = 0;
  s->dtype       = dtype;
  s->shape       = shape;
  s->dim         = dim;
  s->itype       = nm_yale_storage_itype_by_shape(shape);

  // See if a higher itype has been requested.
  if (static_cast<int8_t>(s->itype) < static_cast<int8_t>(min_itype))
    s->itype    = min_itype;

  return s;
}

YALE_STORAGE* nm_yale_storage_create_from_old_yale(nm::dtype_t dtype, size_t* shape, void* ia, void* ja, void* a, nm::dtype_t from_dtype) {

  NAMED_LRI_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::create_from_old_yale, YALE_STORAGE*, nm::dtype_t dtype, size_t* shape, void* r_ia, void* r_ja, void* r_a);

  // With C++ templates, we don't want to have a 4-parameter template. That would be LDType, RDType, LIType, RIType.
  // We can prevent that by copying ia and ja into the correct itype (if necessary) before passing them to the yale
  // copy constructor.
  nm::itype_t to_itype = nm_yale_storage_itype_by_shape(shape);

  return ttable[dtype][from_dtype][to_itype](dtype, shape, ia, ja, a);

}

//////////////////////////////////////////////
// YALE-SPECIFIC FUNCTIONS (RUBY ACCESSORS) //
//////////////////////////////////////////////

/*
 * call-seq:
 *     yale_size -> Integer
 *
 * Get the size of a Yale matrix (the number of elements actually stored).
 *
 * For capacity (the maximum number of elements that can be stored without a resize), use capacity instead.
 */
static VALUE nm_size(VALUE self) {
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  return rubyobj_from_cval_by_itype((char*)(s->ija) + ITYPE_SIZES[s->itype]*(s->shape[0]), s->itype).rval;
}


/*
 * call-seq:
 *     yale_a -> Array
 *     yale_d(index) -> ...
 *
 * Get the A array of a Yale matrix (which stores the diagonal and the LU portions of the matrix).
 */
static VALUE nm_a(int argc, VALUE* argv, VALUE self) {
  VALUE idx;
  rb_scan_args(argc, argv, "01", &idx);

  YALE_STORAGE* s = NM_STORAGE_YALE(self);
  size_t size = nm_yale_storage_get_size(s);

  if (idx == Qnil) {
    VALUE* vals = ALLOCA_N(VALUE, size);

    for (size_t i = 0; i < size; ++i) {
      vals[i] = rubyobj_from_cval((char*)(s->a) + DTYPE_SIZES[s->dtype]*i, s->dtype).rval;
    }
    VALUE ary = rb_ary_new4(size, vals);

    for (size_t i = size; i < s->capacity; ++i)
      rb_ary_push(ary, Qnil);

    return ary;
  } else {
    size_t index = FIX2INT(idx);
    if (index >= size) rb_raise(rb_eRangeError, "out of range");

    return rubyobj_from_cval((char*)(s->a) + DTYPE_SIZES[s->dtype] * index, s->dtype).rval;
  }
}


/*
 * call-seq:
 *     yale_d -> Array
 *     yale_d(index) -> ...
 *
 * Get the diagonal ("D") portion of the A array of a Yale matrix.
 */
static VALUE nm_d(int argc, VALUE* argv, VALUE self) {
  VALUE idx;
  rb_scan_args(argc, argv, "01", &idx);

  YALE_STORAGE* s = NM_STORAGE_YALE(self);

  if (idx == Qnil) {
    VALUE* vals = ALLOCA_N(VALUE, s->shape[0]);

    for (size_t i = 0; i < s->shape[0]; ++i) {
      vals[i] = rubyobj_from_cval((char*)(s->a) + DTYPE_SIZES[s->dtype]*i, s->dtype).rval;
    }
    return rb_ary_new4(s->shape[0], vals);
  } else {
    size_t index = FIX2INT(idx);
    if (index >= s->shape[0]) rb_raise(rb_eRangeError, "out of range");

    return rubyobj_from_cval((char*)(s->a) + DTYPE_SIZES[s->dtype] * index, s->dtype).rval;
  }
}

/*
 * call-seq:
 *     yale_lu -> Array
 *
 * Get the non-diagonal ("LU") portion of the A array of a Yale matrix.
 */
static VALUE nm_lu(VALUE self) {
  YALE_STORAGE* s = NM_STORAGE_YALE(self);

  size_t size = nm_yale_storage_get_size(s);

  VALUE* vals = ALLOCA_N(VALUE, size - s->shape[0] - 1);

  for (size_t i = 0; i < size - s->shape[0] - 1; ++i) {
    vals[i] = rubyobj_from_cval((char*)(s->a) + DTYPE_SIZES[s->dtype]*(s->shape[0] + 1 + i), s->dtype).rval;
  }

  VALUE ary = rb_ary_new4(size - s->shape[0] - 1, vals);

  for (size_t i = size; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}

/*
 * call-seq:
 *     yale_ia -> Array
 *
 * Get the IA portion of the IJA array of a Yale matrix. This gives the start and end positions of rows in the
 * JA and LU portions of the IJA and A arrays, respectively.
 */
static VALUE nm_ia(VALUE self) {
  YALE_STORAGE* s = NM_STORAGE_YALE(self);

  VALUE* vals = ALLOCA_N(VALUE, s->shape[0] + 1);

  for (size_t i = 0; i < s->shape[0] + 1; ++i) {
    vals[i] = rubyobj_from_cval_by_itype((char*)(s->ija) + ITYPE_SIZES[s->itype]*i, s->itype).rval;
  }

  return rb_ary_new4(s->shape[0]+1, vals);
}

/*
 * call-seq:
 *     yale_ja -> Array
 *
 * Get the JA portion of the IJA array of a Yale matrix. This gives the column indices for entries in corresponding
 * positions in the LU portion of the A array.
 */
static VALUE nm_ja(VALUE self) {
  YALE_STORAGE* s = NM_STORAGE_YALE(self);

  size_t size = nm_yale_storage_get_size(s);

  VALUE* vals = ALLOCA_N(VALUE, size - s->shape[0] - 1);

  for (size_t i = 0; i < size - s->shape[0] - 1; ++i) {
    vals[i] = rubyobj_from_cval_by_itype((char*)(s->ija) + ITYPE_SIZES[s->itype]*(s->shape[0] + 1 + i), s->itype).rval;
  }

  VALUE ary = rb_ary_new4(size - s->shape[0] - 1, vals);

  for (size_t i = size; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}

/*
 * call-seq:
 *     yale_ija -> Array
 *     yale_ija(index) -> ...
 *
 * Get the IJA array of a Yale matrix (or a component of the IJA array).
 */
static VALUE nm_ija(int argc, VALUE* argv, VALUE self) {
  VALUE idx;
  rb_scan_args(argc, argv, "01", &idx);

  YALE_STORAGE* s = NM_STORAGE_YALE(self);
  size_t size = nm_yale_storage_get_size(s);

  if (idx == Qnil) {

    VALUE* vals = ALLOCA_N(VALUE, size);

    for (size_t i = 0; i < size; ++i) {
      vals[i] = rubyobj_from_cval_by_itype((char*)(s->ija) + ITYPE_SIZES[s->itype]*i, s->itype).rval;
    }

   VALUE ary = rb_ary_new4(size, vals);

    for (size_t i = size; i < s->capacity; ++i)
      rb_ary_push(ary, Qnil);

    return ary;

  } else {
    size_t index = FIX2INT(idx);
    if (index >= size) rb_raise(rb_eRangeError, "out of range");

    return rubyobj_from_cval_by_itype((char*)(s->ija) + ITYPE_SIZES[s->itype] * index, s->itype).rval;
  }
}


/*
 * call-seq:
 *     yale_nd_row -> ...
 *
 * This function gets the non-diagonal contents of a Yale matrix row.
 * The first argument should be the row index. The optional second argument may be :hash or :array, but defaults
 * to :hash. If :array is given, it will only return the Hash keys (the column indices).
 *
 * This function is meant to accomplish its purpose as efficiently as possible. It does not check for appropriate
 * range.
 *
 * FIXME: :array doesn't make sense. This should be :keys or :values to indicate which array we want.
 */
static VALUE nm_nd_row(int argc, VALUE* argv, VALUE self) {
  VALUE i_, as;
  rb_scan_args(argc, argv, "11", &i_, &as);

  bool array = false;
  if (as != Qnil && rb_to_id(as) != nm_rb_hash) array = true;

  size_t i = FIX2INT(i_);

  YALE_STORAGE* s   = NM_STORAGE_YALE(self);
  nm::dtype_t dtype = NM_DTYPE(self);
  nm::itype_t itype = NM_ITYPE(self);

  // get the position as a size_t
  // TODO: Come up with a faster way to get this than transforming to a Ruby object first.
  size_t pos = FIX2INT(rubyobj_from_cval_by_itype((char*)(s->ija) + ITYPE_SIZES[itype]*i, itype).rval);
  size_t nextpos = FIX2INT(rubyobj_from_cval_by_itype((char*)(s->ija) + ITYPE_SIZES[itype]*(i+1), itype).rval);
  size_t diff = nextpos - pos;

  //std::cerr << "diff = " << diff << "\tpos = " << pos << "\tnextpos = " << nextpos << std::endl;

  VALUE ret; // HERE
  if (array) {
    ret = rb_ary_new3(diff);

    for (size_t idx = pos; idx < nextpos; ++idx) {
      rb_ary_store(ret, idx - pos, rubyobj_from_cval_by_itype((char*)(s->ija) + ITYPE_SIZES[s->itype]*idx, s->itype).rval);
    }

  } else {
    ret = rb_hash_new();

    for (size_t idx = pos; idx < nextpos; ++idx) {
      rb_hash_aset(ret, rubyobj_from_cval_by_itype((char*)(s->ija) + ITYPE_SIZES[s->itype]*idx, s->itype).rval,
                        rubyobj_from_cval((char*)(s->a) + DTYPE_SIZES[s->dtype]*idx, s->dtype).rval);
    }
  }

  return ret;
}

/*
 * call-seq:
 *     yale_vector_insert -> Fixnum
 *
 * Insert at position pos an array of non-diagonal elements with column indices given. Note that the column indices and values
 * must be storage-contiguous -- that is, you can't insert them around existing elements in some row, only amid some
 * elements in some row. You *can* insert them around a diagonal element, since this is stored separately. This function
 * may not be used for the insertion of diagonal elements in most cases, as these are already present in the data
 * structure and are typically modified by replacement rather than insertion.
 *
 * The last argument, pos, may be nil if you want to insert at the beginning of a row. Otherwise it needs to be provided.
 * Don't expect this function to know the difference. It really does very little checking, because its goal is to make
 * multiple contiguous insertion as quick as possible.
 *
 * You should also not attempt to insert values which are the default (0). These are not supposed to be stored, and may
 * lead to undefined behavior.
 *
 * Example:
 *    m.yale_vector_insert(3, [0,3,4], [1,1,1], 15)
 *
 * The example above inserts the values 1, 1, and 1 in columns 0, 3, and 4, assumed to be located at position 15 (which
 * corresponds to row 3).
 *
 * Example:
 *    next = m.yale_vector_insert(3, [0,3,4], [1,1,1])
 *
 * This example determines that i=3 is at position 15 automatically. The value returned, next, is the position where the
 * next value(s) should be inserted.
 */
static VALUE nm_vector_insert(int argc, VALUE* argv, VALUE self) { //, VALUE i_, VALUE jv, VALUE vv, VALUE pos_) {

  // i, jv, vv are mandatory; pos is optional; thus "31"
  VALUE i_, jv, vv, pos_;
  rb_scan_args(argc, argv, "31", &i_, &jv, &vv, &pos_);

  size_t len   = RARRAY_LEN(jv); // need length in order to read the arrays in
  size_t vvlen = RARRAY_LEN(vv);
  if (len != vvlen)
    rb_raise(rb_eArgError, "lengths must match between j array (%d) and value array (%d)", len, vvlen);

  YALE_STORAGE* s   = NM_STORAGE_YALE(self);
  nm::dtype_t dtype = NM_DTYPE(self);
  nm::itype_t itype = NM_ITYPE(self);

  size_t i   = FIX2INT(i_);    // get the row

  // get the position as a size_t
  // TODO: Come up with a faster way to get this than transforming to a Ruby object first.
  if (pos_ == Qnil) pos_ = rubyobj_from_cval_by_itype((char*)(s->ija) + ITYPE_SIZES[itype]*i, itype).rval;
  size_t pos = FIX2INT(pos_);

  // Allocate the j array and the values array
  size_t* j  = ALLOCA_N(size_t, len);
  void* vals = ALLOCA_N(char, DTYPE_SIZES[dtype] * len);

  // Copy array contents
  for (size_t idx = 0; idx < len; ++idx) {
    j[idx] = FIX2INT(rb_ary_entry(jv, idx));
    rubyval_to_cval(rb_ary_entry(vv, idx), dtype, (char*)vals + idx * DTYPE_SIZES[dtype]);
  }

  char ins_type = nm_yale_storage_vector_insert(s, pos, j, vals, len, false, dtype, itype);
  nm_yale_storage_increment_ia_after(s, s->shape[0], i, len, itype);
  s->ndnz += len;

  // Return the updated position
  pos += len;
  return INT2FIX(pos);
}


} // end of extern "C" block
