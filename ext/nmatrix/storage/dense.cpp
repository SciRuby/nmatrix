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
// == dense.c
//
// Dense n-dimensional matrix storage.

/*
 * Standard Includes
 */

#include <ruby.h>

/*
 * Project Includes
 */
// #include "types.h"
#include "util/math.h"

#include "data/data.h"
#include "common.h"
#include "dense.h"

/*
 * Macros
 */

/*
 * Global Variables
 */

/*
 * Forward Declarations
 */

namespace nm { namespace dense_storage {

  template<typename LDType, typename RDType>
  void ref_slice_copy_transposed(const DENSE_STORAGE* rhs, DENSE_STORAGE* lhs);

  template <typename LDType, typename RDType>
  DENSE_STORAGE* cast_copy(const DENSE_STORAGE* rhs, nm::dtype_t new_dtype);
	
	template <typename LDType, typename RDType>
	bool eqeq(const DENSE_STORAGE* left, const DENSE_STORAGE* right);
	
	template <ewop_t op, typename LDType, typename RDType>
	static DENSE_STORAGE* ew_op(const DENSE_STORAGE* left, const DENSE_STORAGE* right, const void* rscalar);

  template <typename DType>
  static DENSE_STORAGE* matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector);

  template <typename DType>
  bool is_hermitian(const DENSE_STORAGE* mat, int lda);

  template <typename DType>
  bool is_symmetric(const DENSE_STORAGE* mat, int lda);

}} // end of namespace nm::dense_storage


extern "C" {

static size_t* stride(size_t* shape, size_t dim);
static void slice_copy(DENSE_STORAGE *dest, const DENSE_STORAGE *src, size_t* lengths, size_t pdest, size_t psrc, size_t n);

/*
 * Functions
 */

///////////////
// Lifecycle //
///////////////

/*
 * Note that elements and elements_length are for initial value(s) passed in.
 * If they are the correct length, they will be used directly. If not, they
 * will be concatenated over and over again into a new elements array. If
 * elements is NULL, the new elements array will not be initialized.
 */
DENSE_STORAGE* nm_dense_storage_create(nm::dtype_t dtype, size_t* shape, size_t dim, void* elements, size_t elements_length) {
  DENSE_STORAGE* s = ALLOC( DENSE_STORAGE );

  s->dim        = dim;
  s->shape      = shape;
  s->dtype      = dtype;

  s->offset     = ALLOC_N(size_t, dim);
  memset(s->offset, 0, sizeof(size_t)*dim);

  s->stride     = stride(shape, dim);
  s->count      = 1;
  s->src        = s;
	
	size_t count  = nm_storage_count_max_elements(s);

  if (elements_length == count) {
  	s->elements = elements;
    
  } else {
    s->elements = ALLOC_N(char, DTYPE_SIZES[dtype]*count);

    size_t copy_length = elements_length;

    if (elements_length > 0) {
      // Repeat elements over and over again until the end of the matrix.
      for (size_t i = 0; i < count; i += elements_length) {

        if (i + elements_length > count) {
        	copy_length = count - i;
        }
        
        memcpy((char*)(s->elements)+i*DTYPE_SIZES[dtype], (char*)(elements)+(i % elements_length)*DTYPE_SIZES[dtype], copy_length*DTYPE_SIZES[dtype]);
      }

      // Get rid of the init_val.
      free(elements);
    }
  }

  return s;
}

/*
 * Destructor for dense storage
 */
void nm_dense_storage_delete(STORAGE* s) {
  // Sometimes Ruby passes in NULL storage for some reason (probably on copy construction failure).
  if (s) {
    DENSE_STORAGE* storage = (DENSE_STORAGE*)s;
    if(storage->count-- == 1) {
      free(storage->shape);
      free(storage->offset);
      free(storage->stride);
      if (storage->elements != NULL)
        free(storage->elements);
      free(storage);
    }
  }
}

/*
 * Destructor for dense storage references (slicing).
 */
void nm_dense_storage_delete_ref(STORAGE* s) {
  // Sometimes Ruby passes in NULL storage for some reason (probably on copy construction failure).
  if (s) {
    DENSE_STORAGE* storage = (DENSE_STORAGE*)s;
    nm_dense_storage_delete( reinterpret_cast<STORAGE*>(storage->src) );
    free(storage->shape);
    free(storage->offset);
    free(storage);
  }
}

/*
 * Mark values in a dense matrix for garbage collection. This may not be necessary -- further testing required.
 */
void nm_dense_storage_mark(void* storage_base) {
  DENSE_STORAGE* storage = (DENSE_STORAGE*)storage_base;

  if (storage && storage->dtype == nm::RUBYOBJ) {
    VALUE* els = reinterpret_cast<VALUE*>(storage->elements);

  	for (size_t index = nm_storage_count_max_elements(storage); index-- > 0;) {
      rb_gc_mark(els[index]);
    }
  }
}

///////////////
// Accessors //
///////////////


VALUE nm_dense_each_with_indices(VALUE nmatrix) {
  volatile VALUE nm = nmatrix;

  RETURN_ENUMERATOR(nm, 0, 0);

  DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);

  // Create indices and initialize them to zero
  size_t* coords = ALLOCA_N(size_t, s->dim);
  memset(coords, 0, sizeof(size_t) * s->dim);

  size_t slice_index;
  size_t* shape_copy = ALLOC_N(size_t, s->dim);
  memcpy(shape_copy, s->shape, sizeof(size_t) * s->dim);

  DENSE_STORAGE* sliced_dummy = nm_dense_storage_create(s->dtype, shape_copy, s->dim, NULL, nm_storage_count_max_elements(s));
  

  for (size_t k = 0; k < nm_storage_count_max_elements(s); ++k) {
    nm_dense_storage_coords(sliced_dummy, k, coords);
    slice_index = nm_dense_storage_pos(s, coords);
    VALUE ary = rb_ary_new();
    if (NM_DTYPE(nm) == nm::RUBYOBJ) rb_ary_push(ary, reinterpret_cast<VALUE*>(s->elements)[slice_index]);
    else rb_ary_push(ary, rubyobj_from_cval((char*)(s->elements) + slice_index*DTYPE_SIZES[NM_DTYPE(nm)], NM_DTYPE(nm)).rval);

    for (size_t p = 0; p < s->dim; ++p) {
      rb_ary_push(ary, INT2FIX(coords[p]));
    }

    // yield the array which now consists of the value and the indices
    rb_yield(ary);

  }

  nm_dense_storage_delete(sliced_dummy);

}


/*
 * Borrowed this function from NArray. Handles 'each' iteration on a dense
 * matrix.
 *
 * Additionally, handles separately matrices containing VALUEs and matrices
 * containing other types of data.
 */
VALUE nm_dense_each(VALUE nmatrix) {
  volatile VALUE nm = nmatrix; // Not sure this actually does anything.
  DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);

  RETURN_ENUMERATOR(nm, 0, 0);

  size_t* temp_coords = ALLOCA_N(size_t, s->dim);
  size_t sliced_index;
  size_t* shape_copy = ALLOC_N(size_t, s->dim);
  memcpy(shape_copy, s->shape, sizeof(size_t) * s->dim);
  DENSE_STORAGE* sliced_dummy = nm_dense_storage_create(s->dtype, shape_copy, s->dim, NULL, nm_storage_count_max_elements(s));

  if (NM_DTYPE(nm) == nm::RUBYOBJ) {

    // matrix of Ruby objects -- yield those objects directly
    for (size_t i = 0; i < nm_storage_count_max_elements(s); ++i)
      nm_dense_storage_coords(sliced_dummy, i, temp_coords);
      sliced_index = nm_dense_storage_pos(s, temp_coords);
      rb_yield( reinterpret_cast<VALUE*>(s->elements)[sliced_index] );

  } else {

    // We're going to copy the matrix element into a Ruby VALUE and then operate on it. This way user can't accidentally
    // modify it and cause a seg fault.
    for (size_t i = 0; i < nm_storage_count_max_elements(s); ++i) {
      nm_dense_storage_coords(sliced_dummy, i, temp_coords);
      sliced_index = nm_dense_storage_pos(s, temp_coords);
      VALUE v = rubyobj_from_cval((char*)(s->elements) + sliced_index*DTYPE_SIZES[NM_DTYPE(nm)], NM_DTYPE(nm)).rval;
      rb_yield( v ); // yield to the copy we made
    }
  }

  nm_dense_storage_delete(sliced_dummy);
}


/*
 * Get a slice or one element, using copying.
 *
 * FIXME: Template the first condition.
 */
void* nm_dense_storage_get(STORAGE* storage, SLICE* slice) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)storage;
  DENSE_STORAGE* ns;

  if (slice->single)
    return (char*)(s->elements) + nm_dense_storage_pos(s, slice->coords) * DTYPE_SIZES[s->dtype];
  else { // Make references
    size_t *shape      = ALLOC_N(size_t, s->dim);
    for (size_t i = 0; i < s->dim; ++i) {
      shape[i]  = slice->lengths[i];
    }

    ns = nm_dense_storage_create(s->dtype, shape, s->dim, NULL, 0); 

    slice_copy(ns, 
        reinterpret_cast<const DENSE_STORAGE*>(s->src), 
        slice->lengths, 
        0, 
        nm_dense_storage_pos(s, slice->coords), 
        0);
    return ns;
  }
}

/*
 * Get a slice or one element by reference (no copy).
 *
 * FIXME: Template the first condition.
 */
void* nm_dense_storage_ref(STORAGE* storage, SLICE* slice) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)storage;

  if (slice->single)
    return (char*)(s->elements) + nm_dense_storage_pos(s, slice->coords) * DTYPE_SIZES[s->dtype];
    
  else {
    DENSE_STORAGE* ns = ALLOC( DENSE_STORAGE );
    ns->dim        = s->dim;
    ns->dtype      = s->dtype;
    ns->offset     = ALLOC_N(size_t, ns->dim);
    ns->shape      = ALLOC_N(size_t, ns->dim);

    for (size_t i = 0; i < ns->dim; ++i) {
      ns->offset[i] = slice->coords[i] + s->offset[i];
      ns->shape[i]  = slice->lengths[i];
    }

    ns->stride     = s->stride;
    ns->elements   = s->elements;
    
    s->src->count++;
    ns->src = s->src;

    return ns;
  }
}


/*
 * Does not free passed-in value! Different from list_storage_insert.
 */
void nm_dense_storage_set(STORAGE* storage, SLICE* slice, void* val) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)storage;
  memcpy((char*)(s->elements) + nm_dense_storage_pos(s, slice->coords) * DTYPE_SIZES[s->dtype], val, DTYPE_SIZES[s->dtype]);
}

///////////
// Tests //
///////////

/*
 * Do these two dense matrices have the same contents?
 *
 * TODO: Test the shape of the two matrices.
 * TODO: See if using memcmp is faster when the left- and right-hand matrices
 *				have the same dtype.
 */
bool nm_dense_storage_eqeq(const STORAGE* left, const STORAGE* right) {
	LR_DTYPE_TEMPLATE_TABLE(nm::dense_storage::eqeq, bool, const DENSE_STORAGE*, const DENSE_STORAGE*);
	
	return ttable[left->dtype][right->dtype]((const DENSE_STORAGE*)left, (const DENSE_STORAGE*)right);
}

/*
 * Test to see if the matrix is Hermitian.  If the matrix does not have a
 * dtype of Complex64 or Complex128 this is the same as testing for symmetry.
 */
bool nm_dense_storage_is_hermitian(const DENSE_STORAGE* mat, int lda) {
	if (mat->dtype == nm::COMPLEX64) {
		return nm::dense_storage::is_hermitian<nm::Complex64>(mat, lda);
		
	} else if (mat->dtype == nm::COMPLEX128) {
		return nm::dense_storage::is_hermitian<nm::Complex128>(mat, lda);
		
	} else {
		return nm_dense_storage_is_symmetric(mat, lda);
	}
}

/*
 * Is this dense matrix symmetric about the diagonal?
 */
bool nm_dense_storage_is_symmetric(const DENSE_STORAGE* mat, int lda) {
	DTYPE_TEMPLATE_TABLE(nm::dense_storage::is_symmetric, bool, const DENSE_STORAGE*, int);
	
	return ttable[mat->dtype](mat, lda);
}

//////////
// Math //
//////////

/*
 * Dense matrix-matrix and matrix-scalar element-wise operations.
 *
 * right or rscalar should be NULL; they should not both be initialized. If right is NULL, it'll use the scalar value instead.
 */
STORAGE* nm_dense_storage_ew_op(nm::ewop_t op, const STORAGE* left, const STORAGE* right, VALUE scalar) {
	OP_LR_DTYPE_TEMPLATE_TABLE(nm::dense_storage::ew_op, DENSE_STORAGE*, const DENSE_STORAGE* left, const DENSE_STORAGE* right, const void*);

	if (right)
	  return ttable[op][left->dtype][right->dtype](reinterpret_cast<const DENSE_STORAGE*>(left), reinterpret_cast<const DENSE_STORAGE*>(right), NULL);
	else {
	  nm::dtype_t r_dtype = nm_dtype_guess(scalar);
	  void* r_scalar  = ALLOCA_N(char, DTYPE_SIZES[r_dtype]);
	  rubyval_to_cval(scalar, r_dtype, r_scalar);

    return ttable[op][left->dtype][r_dtype](reinterpret_cast<const DENSE_STORAGE*>(left), NULL, r_scalar);
	}
}

/*
 * Dense matrix-matrix multiplication.
 */
STORAGE* nm_dense_storage_matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector) {
  DTYPE_TEMPLATE_TABLE(nm::dense_storage::matrix_multiply, DENSE_STORAGE*, const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector);

  return ttable[casted_storage.left->dtype](casted_storage, resulting_shape, vector);
}

/////////////
// Utility //
/////////////

/*
 * Determine the linear array position (in elements of s) of some set of coordinates
 * (given by slice).
 */
size_t nm_dense_storage_pos(const DENSE_STORAGE* s, const size_t* coords) {
  size_t pos = 0;

  for (size_t i = 0; i < s->dim; ++i)
    pos += (coords[i] + s->offset[i]) * s->stride[i];

  return pos;

}

/*
 * Determine the a set of slice coordinates from linear array position (in elements 
 * of s) of some set of coordinates (given by slice).  (Inverse of
 * nm_dense_storage_pos).  
 *
 * The parameter coords_out should be a pre-allocated array of size equal to s->dim.
 */
void nm_dense_storage_coords(const DENSE_STORAGE* s, const size_t slice_pos, size_t* coords_out) {

  size_t temp_pos = slice_pos;

  for (size_t i = 0; i < s->dim; ++i) {
    coords_out[i] = (temp_pos - temp_pos % s->stride[i])/s->stride[i] - s->offset[i];
    temp_pos = temp_pos % s->stride[i];
  }

}

/*
 * Calculate the stride length.
 */
static size_t* stride(size_t* shape, size_t dim) {
  size_t i, j;
  size_t* stride = ALLOC_N(size_t, dim);

  for (i = 0; i < dim; ++i) {
    stride[i] = 1;
    for (j = i+1; j < dim; ++j) {
      stride[i] *= shape[j];
    }
  }

  return stride;
}

/*
 * Recursive slicing for N-dimensional matrix.
 */
static void slice_copy(DENSE_STORAGE *dest, const DENSE_STORAGE *src, size_t* lengths, size_t pdest, size_t psrc, size_t n) {
  if (src->dim - n > 1) {
    for (size_t i = 0; i < lengths[n]; ++i) {
      slice_copy(dest, src, lengths,
                                    pdest + dest->stride[n]*i,
                                    psrc + src->stride[n]*i, 
                                    n + 1);
    }
  } else {
    memcpy((char*)dest->elements + pdest*DTYPE_SIZES[dest->dtype],
        (char*)src->elements + psrc*DTYPE_SIZES[src->dtype],
        dest->shape[n]*DTYPE_SIZES[dest->dtype]);
  }

}

/////////////////////////
// Copying and Casting //
/////////////////////////

/*
 * Copy dense storage, changing dtype if necessary.
 */
STORAGE* nm_dense_storage_cast_copy(const STORAGE* rhs, nm::dtype_t new_dtype) {
	NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::dense_storage::cast_copy, DENSE_STORAGE*, const DENSE_STORAGE* rhs, nm::dtype_t new_dtype);

	return (STORAGE*)ttable[new_dtype][rhs->dtype]((DENSE_STORAGE*)rhs, new_dtype);
}

/*
 * Copy dense storage without a change in dtype.
 */
DENSE_STORAGE* nm_dense_storage_copy(const DENSE_STORAGE* rhs) {
  size_t  count = 0;  
  size_t *shape  = ALLOC_N(size_t, rhs->dim);

  // copy shape and offset
  for (size_t i = 0; i < rhs->dim; ++i) {
    shape[i]  = rhs->shape[i];
  }

  DENSE_STORAGE* lhs = nm_dense_storage_create(rhs->dtype, shape, rhs->dim, NULL, 0);
  count = nm_storage_count_max_elements(lhs);


	// Ensure that allocation worked before copying.
  if (lhs && count) {
    if (rhs == rhs->src) // not a reference
      memcpy(lhs->elements, rhs->elements, DTYPE_SIZES[rhs->dtype] * count);
    else { // slice whole matrix
      size_t *offset = ALLOC_N(size_t, rhs->dim);
      memset(offset, 0, sizeof(size_t) * rhs->dim);

      slice_copy(lhs,
           reinterpret_cast<const DENSE_STORAGE*>(rhs->src),
           rhs->shape,
           0,
           nm_dense_storage_pos(rhs, offset),
           0);
    }
  }

  return lhs;
}


/*
 * Transpose dense storage into a new dense storage object. Basically a copy constructor.
 *
 * Not much point in templating this as it's pretty straight-forward.
 */
STORAGE* nm_dense_storage_copy_transposed(const STORAGE* rhs_base) {
  DENSE_STORAGE* rhs = (DENSE_STORAGE*)rhs_base;

  size_t *shape = ALLOC_N(size_t, rhs->dim);

  // swap shape and offset
  shape[0] = rhs->shape[1];
  shape[1] = rhs->shape[0];

  DENSE_STORAGE *lhs = nm_dense_storage_create(rhs->dtype, shape, rhs->dim, NULL, 0);
  lhs->offset[0] = rhs->offset[1];
  lhs->offset[1] = rhs->offset[0];

  if (rhs_base->src == rhs_base) {
    nm_math_transpose_generic(rhs->shape[0], rhs->shape[1], rhs->elements, rhs->shape[1], lhs->elements, lhs->shape[1], DTYPE_SIZES[rhs->dtype]);
  } else {
    NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::dense_storage::ref_slice_copy_transposed, void, const DENSE_STORAGE* rhs, DENSE_STORAGE* lhs);
    ttable[lhs->dtype][rhs->dtype](rhs, lhs);
  }

  return (STORAGE*)lhs;
}

} // end of extern "C" block

namespace nm { namespace dense_storage {

/////////////////////////
// Templated Functions //
/////////////////////////

template<typename LDType, typename RDType>
void ref_slice_copy_transposed(const DENSE_STORAGE* rhs, DENSE_STORAGE* lhs) {

  LDType* lhs_els = reinterpret_cast<LDType*>(lhs->elements);
  RDType* rhs_els = reinterpret_cast<RDType*>(rhs->elements);

  size_t count = nm_storage_count_max_elements(lhs);
  size_t* temp_coords = ALLOCA_N(size_t, lhs->dim);
  size_t coord_swap_temp;

  while (count-- > 0) {
    nm_dense_storage_coords(lhs, count, temp_coords);
    NM_SWAP(temp_coords[0], temp_coords[1], coord_swap_temp);
    size_t r_coord = nm_dense_storage_pos(rhs, temp_coords);
    lhs_els[count] = rhs_els[r_coord];
  }

}

template <typename LDType, typename RDType>
DENSE_STORAGE* cast_copy(const DENSE_STORAGE* rhs, dtype_t new_dtype) {
  size_t  count = nm_storage_count_max_elements(rhs);

  size_t *shape = ALLOC_N(size_t, rhs->dim);
  memcpy(shape, rhs->shape, sizeof(size_t) * rhs->dim);

  DENSE_STORAGE* lhs			= nm_dense_storage_create(new_dtype, shape, rhs->dim, NULL, 0);

  RDType*	rhs_els         = reinterpret_cast<RDType*>(rhs->elements);
  LDType* lhs_els	        = reinterpret_cast<LDType*>(lhs->elements);

	// Ensure that allocation worked before copying.
  if (lhs && count) {
    if (rhs->src != rhs) {
      /* Make a copy of a ref to a matrix. */

      DENSE_STORAGE* tmp = nm_dense_storage_copy(rhs);

      RDType* tmp_els    = reinterpret_cast<RDType*>(tmp->elements);
      while (count-- > 0)   {
        lhs_els[count] = tmp_els[count];
      }
      nm_dense_storage_delete(tmp);
    } else {
      /* Make a regular copy. */

    	while (count-- > 0)     		lhs_els[count] = rhs_els[count];
    }
  }
	
  return lhs;
}

template <typename LDType, typename RDType>
bool eqeq(const DENSE_STORAGE* left, const DENSE_STORAGE* right) {
  size_t index;
  DENSE_STORAGE *tmp1, *tmp2;
  tmp1 = NULL; tmp2 = NULL;
  bool result = true;
  /* FIXME: Very strange behavior! The GC calls the method directly with non-initialized data. */
  if (left->dim != right->dim) return false;


	LDType* left_elements	  = (LDType*)left->elements;
  RDType* right_elements	= (RDType*)right->elements;

  // Copy elements in temp matrix if you have refernce to the right.
  if (left->src != left) {
    tmp1 = nm_dense_storage_copy(left);
    left_elements = (LDType*)tmp1->elements;
  }
  if (right->src != right) {
    tmp2 = nm_dense_storage_copy(right);
    right_elements = (RDType*)tmp2->elements;
  }
  


	for (index = nm_storage_count_max_elements(left); index-- > 0;) {
		if (left_elements[index] != right_elements[index]) {
      result = false;
      break;
    }
	}

  if (tmp1)
    free(tmp1);
  if (tmp2)
    free(tmp2);

	return result;
}

template <typename DType>
bool is_hermitian(const DENSE_STORAGE* mat, int lda) {
	unsigned int i, j;
	register DType complex_conj;
	
	const DType* els = (DType*) mat->elements;
	
	for (i = mat->shape[0]; i-- > 0;) {
		for (j = i + 1; j < mat->shape[1]; ++j) {
			complex_conj		= els[j*lda + 1];
			complex_conj.i	= -complex_conj.i;
			
			if (els[i*lda+j] != complex_conj) {
	      return false;
	    }
		}
	}
	
	return true;
}

template <typename DType>
bool is_symmetric(const DENSE_STORAGE* mat, int lda) {
	unsigned int i, j;
	const DType* els = (DType*) mat->elements;
	
	for (i = mat->shape[0]; i-- > 0;) {
		for (j = i + 1; j < mat->shape[1]; ++j) {
			if (els[i*lda+j] != els[j*lda+i]) {
	      return false;
	    }
		}
	}
	
	return true;
}

/*
 * Templated dense storage element-wise operations which return the same DType.
 */
template <ewop_t op, typename LDType, typename RDType>
static DENSE_STORAGE* ew_op(const DENSE_STORAGE* left, const DENSE_STORAGE* right, const void* rscalar) {
	unsigned int count;
  size_t l_count;
  size_t r_count;

	size_t* temp_coords = ALLOCA_N(size_t, left->dim);

	size_t* new_shape = ALLOC_N(size_t, left->dim);
	memcpy(new_shape, left->shape, sizeof(size_t) * left->dim);

  // Determine the return dtype. This depends on the type of operation we're doing. Usually, it's going to be
  // set by the left matrix, but for comparisons, we'll use BYTE (in lieu of boolean).
  dtype_t new_dtype = static_cast<uint8_t>(op) < NUM_NONCOMP_EWOPS ? left->dtype : BYTE;

	DENSE_STORAGE* result = nm_dense_storage_create(new_dtype, new_shape, left->dim, NULL, 0);
	
	LDType* l_elems = reinterpret_cast<LDType*>(left->elements);

	if (right) { // matrix-matrix operation
	  RDType* r_elems = reinterpret_cast<RDType*>(right->elements);

    if (static_cast<uint8_t>(op) < NUM_NONCOMP_EWOPS) { // use left-dtype

      for (count = nm_storage_count_max_elements(result); count-- > 0;) {
        nm_dense_storage_coords(result, count, temp_coords);
        l_count = nm_dense_storage_pos(left, temp_coords);
        r_count = nm_dense_storage_pos(right, temp_coords);
        
        reinterpret_cast<LDType*>(result->elements)[count] = ew_op_switch<op,LDType,RDType>(l_elems[l_count], r_elems[r_count]);
      }

    } else { // new_dtype is BYTE: comparison operators
      uint8_t* res_elems = reinterpret_cast<uint8_t*>(result->elements);

      for (count = nm_storage_count_max_elements(result); count-- > 0;) {
        nm_dense_storage_coords(result, count, temp_coords);
        l_count = nm_dense_storage_pos(left, temp_coords);
        r_count = nm_dense_storage_pos(right, temp_coords);

        switch (op) {
          case EW_EQEQ:
            res_elems[count] = l_elems[l_count] == r_elems[r_count];
            break;

          case EW_NEQ:
            res_elems[count] = l_elems[l_count] != r_elems[r_count];
            break;

          case EW_LT:
            res_elems[count] = l_elems[l_count] < r_elems[r_count];
            break;

          case EW_GT:
            res_elems[count] = l_elems[l_count] > r_elems[r_count];
            break;

          case EW_LEQ:
            res_elems[count] = l_elems[l_count] <= r_elems[r_count];
            break;

          case EW_GEQ:
            res_elems[count] = l_elems[l_count] >= r_elems[r_count];
            break;

          default:
            rb_raise(rb_eStandardError, "this should not happen");
        }
      }
    }

  } else { // matrix-scalar operation
    const RDType* r_elem = reinterpret_cast<const RDType*>(rscalar);

    if (static_cast<uint8_t>(op) < NUM_NONCOMP_EWOPS) { // use left-dtype

      for (count = nm_storage_count_max_elements(result); count-- > 0;) {
        nm_dense_storage_coords(result, count, temp_coords);
        l_count = nm_dense_storage_pos(left, temp_coords);

        reinterpret_cast<LDType*>(result->elements)[count] = ew_op_switch<op,LDType,RDType>(l_elems[l_count], *r_elem);
      }

    } else {
      uint8_t* res_elems = reinterpret_cast<uint8_t*>(result->elements);

      for (count = nm_storage_count_max_elements(result); count-- > 0;) {
        nm_dense_storage_coords(result, count, temp_coords);
        l_count = nm_dense_storage_pos(left, temp_coords);
        
        switch (op) {
          case EW_EQEQ:
            res_elems[count] = l_elems[l_count] == *r_elem;
            break;

          case EW_NEQ:
            res_elems[count] = l_elems[l_count] != *r_elem;
            break;

          case EW_LT:
            res_elems[count] = l_elems[l_count] < *r_elem;
            break;

          case EW_GT:
            res_elems[count] = l_elems[l_count] > *r_elem;
            break;

          case EW_LEQ:
            res_elems[count] = l_elems[l_count] <= *r_elem;
            break;

          case EW_GEQ:
            res_elems[count] = l_elems[l_count] >= *r_elem;
            break;

          default:
            rb_raise(rb_eStandardError, "this should not happen");
        }
      }

    }
  }
	return result;
}


/*
 * DType-templated matrix-matrix multiplication for dense storage.
 */
template <typename DType>
static DENSE_STORAGE* matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector) {
  DENSE_STORAGE *left  = (DENSE_STORAGE*)(casted_storage.left),
                *right = (DENSE_STORAGE*)(casted_storage.right);

  // Create result storage.
  DENSE_STORAGE* result = nm_dense_storage_create(left->dtype, resulting_shape, 2, NULL, 0);

  DType *pAlpha = ALLOCA_N(DType, 1),
        *pBeta  = ALLOCA_N(DType, 1);

  *pAlpha = 1;
  *pBeta = 0;
  // Do the multiplication
  if (vector) nm::math::gemv<DType>(CblasNoTrans, left->shape[0], left->shape[1], pAlpha,
                                    reinterpret_cast<DType*>(left->elements), left->shape[1],
                                    reinterpret_cast<DType*>(right->elements), 1, pBeta,
                                    reinterpret_cast<DType*>(result->elements), 1);
  else        nm::math::gemm<DType>(CblasRowMajor, CblasNoTrans, CblasNoTrans, left->shape[0], right->shape[1], left->shape[1],
                                    pAlpha, reinterpret_cast<DType*>(left->elements), left->shape[1],
                                    reinterpret_cast<DType*>(right->elements), right->shape[1], pBeta,
                                    reinterpret_cast<DType*>(result->elements), result->shape[1]);

  return result;
}

}} // end of namespace nm::dense_storage
