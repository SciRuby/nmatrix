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
#ifndef NM_CHECK_ALLOC
#define NM_CHECK_ALLOC(x) if (!x) rb_raise(rb_eNoMemError, "insufficient memory");
#endif

/*
 * Global Variables
 */

/*
 * Forward Declarations
 */

namespace nm { namespace dense_storage {

  template <typename LDType, typename RDType>
  DENSE_STORAGE* cast_copy(const DENSE_STORAGE* rhs, dtype_t new_dtype);
	
	template <typename LDType, typename RDType>
	bool eqeq(const DENSE_STORAGE* left, const DENSE_STORAGE* right);
	
	template <ewop_t op, typename LDType, typename RDType>
	static DENSE_STORAGE* ew_op(const DENSE_STORAGE* left, const DENSE_STORAGE* right);

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
DENSE_STORAGE* nm_dense_storage_create(dtype_t dtype, size_t* shape, size_t dim, void* elements, size_t elements_length) {
  DENSE_STORAGE* s = ALLOC( DENSE_STORAGE );

  s->dim        = dim;
  s->shape      = shape;
  s->dtype      = dtype;

  s->offset     = ALLOC_N(size_t, dim);
  NM_CHECK_ALLOC(s->offset);

  memset(s->offset, 0, sizeof(size_t)*dim);

  s->stride     = stride(shape, dim);
  s->count      = 1;
  s->src        = s;
	
	size_t count  = nm_storage_count_max_elements(s);

  if (elements_length == count) {
  	s->elements = elements;
  	
  } else {
    s->elements = ALLOC_N(char, DTYPE_SIZES[dtype]*count);
    NM_CHECK_ALLOC(s->elements);

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

  if (storage && storage->dtype == RUBYOBJ) {
    VALUE* els = reinterpret_cast<VALUE*>(storage->elements);

  	for (size_t index = nm_storage_count_max_elements(storage); index-- > 0;) {
      rb_gc_mark(els[index]);
    }
  }
}

///////////////
// Accessors //
///////////////

/*
 * Get a slice or one element, using copying.
 *
 * FIXME: Template the first condition.
 */
void* nm_dense_storage_get(STORAGE* storage, SLICE* slice) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)storage;
  DENSE_STORAGE* ns;
  size_t count;

  if (slice->single)
    return (char*)(s->elements) + nm_dense_storage_pos(s, slice->coords) * DTYPE_SIZES[s->dtype];
  else { // Make references
    ns = ALLOC( DENSE_STORAGE );

    NM_CHECK_ALLOC(ns);

    ns->dim        = s->dim;
    ns->dtype      = s->dtype;

    ns->offset     = ALLOC_N(size_t, ns->dim);
    NM_CHECK_ALLOC(ns->offset);

    ns->shape      = ALLOC_N(size_t, ns->dim);
    NM_CHECK_ALLOC(ns->shape);

    for (size_t i = 0; i < ns->dim; ++i) {
      ns->offset[i] = 0;
      ns->shape[i]  = slice->lengths[i];
    }

    ns->stride     = stride(ns->shape, ns->dim);
    ns->count      = 1;
    ns->src        = ns;

    count          = nm_storage_count_max_elements(s);

    ns->elements   = ALLOC_N(char, DTYPE_SIZES[ns->dtype] * count);
    NM_CHECK_ALLOC(ns->elements);

    slice_copy(ns, s, slice->lengths, 0, nm_dense_storage_pos(s, slice->coords), 0);
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
    NM_CHECK_ALLOC(ns);

    ns->dim        = s->dim;
    ns->dtype      = s->dtype;

    ns->offset     = ALLOC_N(size_t, ns->dim);
    NM_CHECK_ALLOC(ns->offset);

    ns->shape      = ALLOC_N(size_t, ns->dim);
    NM_CHECK_ALLOC(ns->shape);

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
	if (mat->dtype == COMPLEX64) {
		return nm::dense_storage::is_hermitian<nm::Complex64>(mat, lda);
		
	} else if (mat->dtype == COMPLEX128) {
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
 * Dense element-wise operations.
 */
STORAGE* nm_dense_storage_ew_op(nm::ewop_t op, const STORAGE* left, const STORAGE* right) {
	OP_LR_DTYPE_TEMPLATE_TABLE(nm::dense_storage::ew_op, DENSE_STORAGE*, const DENSE_STORAGE*, const DENSE_STORAGE*);

	return ttable[op][left->dtype][right->dtype](reinterpret_cast<const DENSE_STORAGE*>(left), reinterpret_cast<const DENSE_STORAGE*>(right));
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
 * Calculate the stride length.
 */
static size_t* stride(size_t* shape, size_t dim) {
  size_t i, j;
  size_t* stride = ALLOC_N(size_t, dim);

  NM_CHECK_ALLOC(stride);

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
STORAGE* nm_dense_storage_cast_copy(const STORAGE* rhs, dtype_t new_dtype) {
	NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::dense_storage::cast_copy, DENSE_STORAGE*, const DENSE_STORAGE* rhs, dtype_t new_dtype);

	return (STORAGE*)ttable[new_dtype][rhs->dtype]((DENSE_STORAGE*)rhs, new_dtype);
}

/*
 * Copy dense storage without a change in dtype.
 */
DENSE_STORAGE* nm_dense_storage_copy(const DENSE_STORAGE* rhs) {
  size_t  count = nm_storage_count_max_elements(rhs);
  size_t *shape  = ALLOC_N(size_t, rhs->dim);

  NM_CHECK_ALLOC(shape);

  // copy shape and offset
  for (size_t i = 0; i < rhs->dim; ++i) {
    shape[i]  = rhs->shape[i];
  }

  DENSE_STORAGE* lhs = nm_dense_storage_create(rhs->dtype, shape, rhs->dim, NULL, 0);

	// Ensure that allocation worked before copying.
  if (lhs && count) {
    if (rhs == rhs->src) // not a reference
      memcpy(lhs->elements, rhs->elements, DTYPE_SIZES[rhs->dtype] * count);
    else // slice whole matrix
      slice_copy(lhs,
           reinterpret_cast<const DENSE_STORAGE*>(rhs->src),
           rhs->shape,
           0,
           nm_dense_storage_pos(reinterpret_cast<const DENSE_STORAGE*>(rhs->src), rhs->offset),
           0);
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

  NM_CHECK_ALLOC(shape);

  // swap shape and offset
  shape[0] = rhs->shape[1];
  shape[1] = rhs->shape[0];

  DENSE_STORAGE *lhs = nm_dense_storage_create(rhs->dtype, shape, rhs->dim, NULL, 0);
  lhs->offset[0] = rhs->offset[1];
  lhs->offset[1] = rhs->offset[0];

  nm_math_transpose_generic(rhs->shape[0], rhs->shape[1], rhs->elements, rhs->shape[1], lhs->elements, lhs->shape[1], DTYPE_SIZES[rhs->dtype]);

  return (STORAGE*)lhs;
}

} // end of extern "C" block

namespace nm { namespace dense_storage {

/////////////////////////
// Templated Functions //
/////////////////////////

template <typename LDType, typename RDType>
DENSE_STORAGE* cast_copy(const DENSE_STORAGE* rhs, dtype_t new_dtype) {
  size_t  count = nm_storage_count_max_elements(rhs);

  size_t *shape = ALLOC_N(size_t, rhs->dim);

  NM_CHECK_ALLOC(shape); // presumably we couldn't allocate offset if shape allocation failed.

  memcpy(shape, rhs->shape, sizeof(size_t) * rhs->dim);

  DENSE_STORAGE* lhs			= nm_dense_storage_create(new_dtype, shape, rhs->dim, NULL, 0);
  memcpy(lhs->offset, rhs->offset, sizeof(size_t) * rhs->dim);

  RDType*	rhs_els         = reinterpret_cast<RDType*>(rhs->elements);
  LDType* lhs_els	        = reinterpret_cast<LDType*>(lhs->elements);

	// Ensure that allocation worked before copying.
  if (lhs && count) {
    if (rhs->src != rhs) {
      /* Make a copy of a ref to a matrix. */

      DENSE_STORAGE* tmp = nm_dense_storage_copy(rhs);
      NM_CHECK_ALLOC(tmp);

      RDType* tmp_els    = reinterpret_cast<RDType*>(tmp->elements);
      while (count-- > 0)         lhs_els[count] = tmp_els[count];
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
  
  /* FIXME: Very strange behavior! The GC calls the method directly with non-initialized data. */
  if (left->dim != right->dim) return false;

	LDType* left_elements	  = (LDType*)left->elements;
	RDType* right_elements	= (RDType*)right->elements;
	
	for (index = nm_storage_count_max_elements(left); index-- > 0;) {
		if (left_elements[index] != right_elements[index]) return false;
	}

	return true;
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
 * Templated dense storage element-wise operations.
 */
template <ewop_t op, typename LDType, typename RDType>
static DENSE_STORAGE* ew_op(const DENSE_STORAGE* left, const DENSE_STORAGE* right) {
	unsigned int count;
	
	size_t* new_shape = (size_t*)calloc(left->dim, sizeof(size_t));
	memcpy(new_shape, left->shape, sizeof(size_t) * left->dim);
	
	DENSE_STORAGE* result = nm_dense_storage_create(left->dtype, new_shape, left->dim, NULL, 0);
	
	LDType* l_elems = reinterpret_cast<LDType*>(left->elements);
	RDType* r_elems = reinterpret_cast<RDType*>(right->elements);
	
	LDType* res_elems = reinterpret_cast<LDType*>(result->elements);

	for (count = nm_storage_count_max_elements(result); count-- > 0;) {
		switch (op) {
			case EW_ADD:
				res_elems[count] = l_elems[count] + r_elems[count];
				break;
				
			case EW_SUB:
				res_elems[count] = l_elems[count] - r_elems[count];
				break;
				
			case EW_MUL:
				res_elems[count] = l_elems[count] * r_elems[count];
				break;
				
			case EW_DIV:
				res_elems[count] = l_elems[count] / r_elems[count];
				break;
				
			case EW_MOD:
				rb_raise(rb_eNotImpError, "Element-wise modulo is currently not supported.");
				break;
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

  DType *pAlpha = new DType(1),
        *pBeta  = new DType(0);

  // Do the multiplication
  bool succ;

  if (vector) succ = nm::math::gemv<DType>(CblasNoTrans, left->shape[0], left->shape[1], pAlpha,
                                           reinterpret_cast<DType*>(left->elements), left->shape[1],
                                           reinterpret_cast<DType*>(right->elements), 1, pBeta,
                                           reinterpret_cast<DType*>(result->elements), 1);
  else        succ = nm::math::gemm<DType>(CblasNoTrans, CblasNoTrans, right->shape[1], left->shape[0], left->shape[1], pAlpha,
                                           reinterpret_cast<DType*>(right->elements), right->shape[1],
                                           reinterpret_cast<DType*>(left->elements), left->shape[1], pBeta,
                                           reinterpret_cast<DType*>(result->elements), result->shape[1]);

  delete pAlpha;
  delete pBeta;

  if (!succ)
    rb_raise(rb_eStandardError, "gemm/gemv failed for an unknown reason");

  return result;
}

}} // end of namespace nm::dense_storage
