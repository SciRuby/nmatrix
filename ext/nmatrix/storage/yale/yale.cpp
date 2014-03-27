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
#include <array>
#include <typeinfo>
#include <tuple>
#include <queue>

/*
 * Project Includes
 */

// #include "types.h"
#include "../../data/data.h"
#include "../../math/math.h"

#include "../common.h"

#include "../../nmatrix.h"
#include "../../data/meta.h"

#include "iterators/base.h"
#include "iterators/stored_diagonal.h"
#include "iterators/row_stored_nd.h"
#include "iterators/row_stored.h"
#include "iterators/row.h"
#include "iterators/iterator.h"
#include "class.h"
#include "yale.h"
#include "../../ruby_constants.h"

/*
 * Macros
 */

#ifndef NM_MAX
#define NM_MAX(a,b) (((a)>(b))?(a):(b))
#define NM_MIN(a,b) (((a)<(b))?(a):(b))
#endif

/*
 * Forward Declarations
 */

extern "C" {
  static YALE_STORAGE*	alloc(nm::dtype_t dtype, size_t* shape, size_t dim);

  static size_t yale_count_slice_copy_ndnz(const YALE_STORAGE* s, size_t*, size_t*);

  static void* default_value_ptr(const YALE_STORAGE* s);
  static VALUE default_value(const YALE_STORAGE* s);
  static VALUE obj_at(YALE_STORAGE* s, size_t k);

  /* Ruby-accessible functions */
  static VALUE nm_size(VALUE self);
  static VALUE nm_a(int argc, VALUE* argv, VALUE self);
  static VALUE nm_d(int argc, VALUE* argv, VALUE self);
  static VALUE nm_lu(VALUE self);
  static VALUE nm_ia(VALUE self);
  static VALUE nm_ja(VALUE self);
  static VALUE nm_ija(int argc, VALUE* argv, VALUE self);
  static VALUE nm_row_keys_intersection(VALUE m1, VALUE ii1, VALUE m2, VALUE ii2);

  static VALUE nm_nd_row(int argc, VALUE* argv, VALUE self);

  static inline size_t src_ndnz(const YALE_STORAGE* s) {
    return reinterpret_cast<YALE_STORAGE*>(s->src)->ndnz;
  }

} // end extern "C" block

namespace nm { namespace yale_storage {

template <typename LD, typename RD>
static VALUE map_merged_stored(VALUE left, VALUE right, VALUE init);

template <typename DType>
static bool						ndrow_is_empty(const YALE_STORAGE* s, IType ija, const IType ija_next);

template <typename LDType, typename RDType>
static bool						ndrow_eqeq_ndrow(const YALE_STORAGE* l, const YALE_STORAGE* r, IType l_ija, const IType l_ija_next, IType r_ija, const IType r_ija_next);

template <typename LDType, typename RDType>
static bool           eqeq(const YALE_STORAGE* left, const YALE_STORAGE* right);

template <typename LDType, typename RDType>
static bool eqeq_different_defaults(const YALE_STORAGE* s, const LDType& s_init, const YALE_STORAGE* t, const RDType& t_init);

static void						increment_ia_after(YALE_STORAGE* s, IType ija_size, IType i, long n);

static IType				  insert_search(YALE_STORAGE* s, IType left, IType right, IType key, bool& found);

template <typename DType>
static char           vector_insert(YALE_STORAGE* s, size_t pos, size_t* j, void* val_, size_t n, bool struct_only);

template <typename DType>
static char           vector_insert_resize(YALE_STORAGE* s, size_t current_size, size_t pos, size_t* j, size_t n, bool struct_only);

template <typename DType>
static std::tuple<long,bool,std::queue<std::tuple<IType,IType,int> > > count_slice_set_ndnz_change(YALE_STORAGE* s, size_t* coords, size_t* lengths, DType* v, size_t v_size);

static inline IType* IJA(const YALE_STORAGE* s) {
  return reinterpret_cast<YALE_STORAGE*>(s->src)->ija;
}

static inline IType IJA_SET(const YALE_STORAGE* s, size_t loc, IType val) {
  return IJA(s)[loc] = val;
}

template <typename DType>
static inline DType* A(const YALE_STORAGE* s) {
  return reinterpret_cast<DType*>(reinterpret_cast<YALE_STORAGE*>(s->src)->a);
}

template <typename DType>
static inline DType A_SET(const YALE_STORAGE* s, size_t loc, DType val) {
  return A<DType>(s)[loc] = val;
}


/*
 * Functions
 */

/*
 * Copy a vector from one DType to another.
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



/*
 * Create Yale storage from IA, JA, and A vectors given in Old Yale format (probably from a file, since NMatrix only uses
 * new Yale for its storage).
 *
 * This function is needed for Matlab .MAT v5 IO.
 */
template <typename LDType, typename RDType>
YALE_STORAGE* create_from_old_yale(dtype_t dtype, size_t* shape, char* r_ia, char* r_ja, char* r_a) {
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
  YALE_STORAGE* s = alloc(dtype, shape, 2);

  s->capacity = shape[0] + ndnz + 1;
  s->ndnz     = ndnz;

  // Setup IJA and A arrays
  s->ija = NM_ALLOC_N( IType, s->capacity );
  s->a   = NM_ALLOC_N( LDType, s->capacity );
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
 * Empty the matrix by initializing the IJA vector and setting the diagonal to 0.
 *
 * Called when most YALE_STORAGE objects are created.
 *
 * Can't go inside of class YaleStorage because YaleStorage creation requires that
 * IJA already be initialized.
 */
template <typename DType>
void init(YALE_STORAGE* s, void* init_val) {
  IType IA_INIT = s->shape[0] + 1;

  IType* ija = reinterpret_cast<IType*>(s->ija);
  // clear out IJA vector
  for (IType i = 0; i < IA_INIT; ++i) {
    ija[i] = IA_INIT; // set initial values for IJA
  }

  clear_diagonal_and_zero<DType>(s, init_val);
}


template <typename LDType, typename RDType>
static YALE_STORAGE* slice_copy(YALE_STORAGE* s) {
  YaleStorage<RDType> y(s);
  return y.template alloc_copy<LDType, false>();
}


/*
 * Template version of copy transposed. This could also, in theory, allow a map -- but transpose.h
 * would need to be updated.
 *
 * TODO: Update for slicing? Update for different dtype in and out? We can cast rather easily without
 * too much modification.
 */
template <typename D>
YALE_STORAGE* copy_transposed(YALE_STORAGE* rhs) {
  YaleStorage<D> y(rhs);
  return y.template alloc_copy_transposed<D, false>();
}


///////////////
// Accessors //
///////////////


/*
 * Determine the number of non-diagonal non-zeros in a not-yet-created copy of a slice or matrix.
 */
template <typename DType>
static size_t count_slice_copy_ndnz(const YALE_STORAGE* s, size_t* offset, size_t* shape) {
  IType* ija = s->ija;
  DType* a   = reinterpret_cast<DType*>(s->a);

  DType ZERO(*reinterpret_cast<DType*>(default_value_ptr(s)));

  // Calc ndnz for the destination
  size_t ndnz  = 0;
  size_t i, j; // indexes of destination matrix
  size_t k, l; // indexes of source matrix
  for (i = 0; i < shape[0]; i++) {
    k = i + offset[0];
    for (j = 0; j < shape[1]; j++) {
      l = j + offset[1];

      if (j == i)  continue;

      if (k == l) { // for diagonal element of source
        if (a[k] != ZERO) ++ndnz;
      } else { // for non-diagonal element
        for (size_t c = ija[k]; c < ija[k+1]; c++) {
          if (ija[c] == l) {
            ++ndnz;
            break;
          }
        }
      }
    }
  }

  return ndnz;
}



/*
 * Get a single element of a yale storage object
 */
template <typename DType>
static void* get_single(YALE_STORAGE* storage, SLICE* slice) {
  YaleStorage<DType> y(storage);
  return reinterpret_cast<void*>(y.get_single_p(slice));
}


/*
 * Returns a reference-slice of a matrix.
 */
template <typename DType>
YALE_STORAGE* ref(YALE_STORAGE* s, SLICE* slice) {
  return YaleStorage<DType>(s).alloc_ref(slice);
}


/*
 * Attempt to set a cell or cells in a Yale matrix.
 */
template <typename DType>
void set(VALUE left, SLICE* slice, VALUE right) {
  YALE_STORAGE* storage = NM_STORAGE_YALE(left);
  YaleStorage<DType> y(storage);
  y.insert(slice, right);
}

///////////
// Tests //
///////////

/*
 * Yale eql? -- for whole-matrix comparison returning a single value.
 */
template <typename LDType, typename RDType>
static bool eqeq(const YALE_STORAGE* left, const YALE_STORAGE* right) {
  return YaleStorage<LDType>(left) == YaleStorage<RDType>(right);
}


//////////
// Math //
//////////

#define YALE_IA(s) (reinterpret_cast<IType*>(s->ija))
#define YALE_IJ(s) (reinterpret_cast<IType*>(s->ija) + s->shape[0] + 1)
#define YALE_COUNT(yale) (yale->ndnz + yale->shape[0])

/////////////
// Utility //
/////////////


/*
 * Binary search for finding the beginning of a slice. Returns the position of the first element which is larger than
 * bound.
 */
IType binary_search_left_boundary(const YALE_STORAGE* s, IType left, IType right, IType bound) {
  if (left > right) return -1;

  IType* ija  = IJA(s);

  if (ija[left] >= bound) return left; // shortcut

  IType mid   = (left + right) / 2;
  IType mid_j = ija[mid];

  if (mid_j == bound)
    return mid;
  else if (mid_j > bound) { // eligible! don't exclude it.
    return binary_search_left_boundary(s, left, mid, bound);
  } else // (mid_j < bound)
    return binary_search_left_boundary(s, mid + 1, right, bound);
}


/*
 * Binary search for returning stored values. Returns a non-negative position, or -1 for not found.
 */
int binary_search(YALE_STORAGE* s, IType left, IType right, IType key) {
  if (s->src != s) throw; // need to fix this quickly

  if (left > right) return -1;

  IType* ija = s->ija;

  IType mid = (left + right)/2;
  IType mid_j = ija[mid];

  if (mid_j == key)
  	return mid;

  else if (mid_j > key)
  	return binary_search(s, left, mid - 1, key);

  else
  	return binary_search(s, mid + 1, right, key);
}


/*
 * Resize yale storage vectors A and IJA, copying values.
 */
static void vector_grow(YALE_STORAGE* s) {
  if (s != s->src) {
    throw; // need to correct this quickly.
  }
  nm_yale_storage_register(s);
  size_t new_capacity = s->capacity * GROWTH_CONSTANT;
  size_t max_capacity = YaleStorage<uint8_t>::max_size(s->shape);

  if (new_capacity > max_capacity) new_capacity = max_capacity;

  IType* new_ija      = NM_ALLOC_N(IType, new_capacity);
  void* new_a         = NM_ALLOC_N(char, DTYPE_SIZES[s->dtype] * new_capacity);

  IType* old_ija      = s->ija;
  void* old_a         = s->a;

  memcpy(new_ija, old_ija, s->capacity * sizeof(IType));
  memcpy(new_a,   old_a,   s->capacity * DTYPE_SIZES[s->dtype]);

  s->capacity         = new_capacity;

  if (s->dtype == nm::RUBYOBJ)
    nm_yale_storage_register_a(new_a, s->capacity * DTYPE_SIZES[s->dtype]);

  NM_FREE(old_ija);
  nm_yale_storage_unregister(s);
  NM_FREE(old_a);
  if (s->dtype == nm::RUBYOBJ)
    nm_yale_storage_unregister_a(new_a, s->capacity * DTYPE_SIZES[s->dtype]);

  s->ija         = new_ija;
  s->a           = new_a;

}


/*
 * Resize yale storage vectors A and IJA in preparation for an insertion.
 */
template <typename DType>
static char vector_insert_resize(YALE_STORAGE* s, size_t current_size, size_t pos, size_t* j, size_t n, bool struct_only) {
  if (s != s->src) throw;

  // Determine the new capacity for the IJA and A vectors.
  size_t new_capacity = s->capacity * GROWTH_CONSTANT;
  size_t max_capacity = YaleStorage<DType>::max_size(s->shape);

  if (new_capacity > max_capacity) {
    new_capacity = max_capacity;

    if (current_size + n > max_capacity) rb_raise(rb_eNoMemError, "insertion size exceeded maximum yale matrix size");
  }

  if (new_capacity < current_size + n)
  	new_capacity = current_size + n;

  nm_yale_storage_register(s);

  // Allocate the new vectors.
  IType* new_ija     = NM_ALLOC_N( IType, new_capacity );
  NM_CHECK_ALLOC(new_ija);

  DType* new_a       = NM_ALLOC_N( DType, new_capacity );
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
  if (s->dtype == nm::RUBYOBJ)
    nm_yale_storage_register_a(new_a, new_capacity);

  NM_FREE(s->ija);
  nm_yale_storage_unregister(s);
  NM_FREE(s->a);
  
  if (s->dtype == nm::RUBYOBJ)
    nm_yale_storage_unregister_a(new_a, new_capacity);

  s->ija = new_ija;
  s->a   = reinterpret_cast<void*>(new_a);

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
template <typename DType>
static char vector_insert(YALE_STORAGE* s, size_t pos, size_t* j, void* val_, size_t n, bool struct_only) {

  if (pos < s->shape[0]) {
    rb_raise(rb_eArgError, "vector insert pos (%lu) is before beginning of ja (%lu); this should not happen", pos, s->shape[0]);
  }

  DType* val = reinterpret_cast<DType*>(val_);

  size_t size = s->ija[s->shape[0]];

  IType* ija = s->ija;
  DType* a   = reinterpret_cast<DType*>(s->a);

  if (size + n > s->capacity) {
    vector_insert_resize<DType>(s, size, pos, j, n, struct_only);

    // Need to get the new locations for ija and a.
  	ija = s->ija;
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
static void increment_ia_after(YALE_STORAGE* s, IType ija_size, IType i, long n) {
  IType* ija = s->ija;

  ++i;
  for (; i <= ija_size; ++i) {
    ija[i] += n;
  }
}

/*
 * Binary search for returning insertion points.
 */
static IType insert_search(YALE_STORAGE* s, IType left, IType right, IType key, bool& found) {

  if (left > right) {
    found = false;
    return left;
  }

  IType* ija = s->ija;
  IType mid = (left + right)/2;
  IType mid_j = ija[mid];

  if (mid_j == key) {
    found = true;
    return mid;

  } else if (mid_j > key) {
  	return insert_search(s, left, mid-1, key, found);

  } else {
  	return insert_search(s, mid+1, right, key, found);
  }
}

/////////////////////////
// Copying and Casting //
/////////////////////////

/*
 * Templated copy constructor for changing dtypes.
 */
template <typename L, typename R>
YALE_STORAGE* cast_copy(const YALE_STORAGE* rhs) {
  YaleStorage<R> y(rhs);
  return y.template alloc_copy<L>();
}

/*
 * Template access for getting the size of Yale storage.
 */
size_t get_size(const YALE_STORAGE* storage) {
  return storage->ija[ storage->shape[0] ];
}


template <typename DType>
static STORAGE* matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector) {
  YALE_STORAGE *left  = (YALE_STORAGE*)(casted_storage.left),
               *right = (YALE_STORAGE*)(casted_storage.right);

  nm_yale_storage_register(left);
  nm_yale_storage_register(right);
  // We can safely get dtype from the casted matrices; post-condition of binary_storage_cast_alloc is that dtype is the
  // same for left and right.
  // int8_t dtype = left->dtype;

  IType* ijl = left->ija;
  IType* ijr = right->ija;

  // First, count the ndnz of the result.
  // TODO: This basically requires running symbmm twice to get the exact ndnz size. That's frustrating. Are there simple
  // cases where we can avoid running it?
  size_t result_ndnz = nm::math::symbmm(resulting_shape[0], left->shape[1], resulting_shape[1], ijl, ijl, true, ijr, ijr, true, NULL, true);

  // Create result storage.
  YALE_STORAGE* result = nm_yale_storage_create(left->dtype, resulting_shape, 2, result_ndnz);
  init<DType>(result, NULL);
  IType* ija = result->ija;

  // Symbolic multiplication step (build the structure)
  nm::math::symbmm(resulting_shape[0], left->shape[1], resulting_shape[1], ijl, ijl, true, ijr, ijr, true, ija, true);

  // Numeric multiplication step (fill in the elements)

  nm::math::numbmm<DType>(result->shape[0], left->shape[1], result->shape[1],
                                ijl, ijl, reinterpret_cast<DType*>(left->a), true,
                                ijr, ijr, reinterpret_cast<DType*>(right->a), true,
                                ija, ija, reinterpret_cast<DType*>(result->a), true);


  // Sort the columns
  nm::math::smmp_sort_columns<DType>(result->shape[0], ija, ija, reinterpret_cast<DType*>(result->a));

  nm_yale_storage_unregister(right);
  nm_yale_storage_unregister(left);
  return reinterpret_cast<STORAGE*>(result);
}


/*
 * Get the sum of offsets from the original matrix (for sliced iteration).
 */
static std::array<size_t,2> get_offsets(YALE_STORAGE* x) {
  std::array<size_t, 2> offsets{ {0,0} };
  while (x != x->src) {
    offsets[0] += x->offset[0];
    offsets[1] += x->offset[1];
    x = reinterpret_cast<YALE_STORAGE*>(x->src);
  }
  return offsets;
}


class RowIterator {
protected:
  YALE_STORAGE* s;
  IType* ija;
  void*  a;
  IType i, k, k_end;
  size_t j_offset, j_shape;
  bool diag, End;
  VALUE init;
public:
  RowIterator(YALE_STORAGE* s_, IType* ija_, IType i_, size_t j_shape_, size_t j_offset_ = 0)
    : s(s_),
      ija(ija_),
      a(reinterpret_cast<YALE_STORAGE*>(s->src)->a),
      i(i_),
      k(ija[i]),
      k_end(ija[i+1]),
      j_offset(j_offset_),
      j_shape(j_shape_),
      diag(row_has_no_nd() || diag_is_first()),
      End(false),
      init(default_value(s))
    { }

  RowIterator(YALE_STORAGE* s_, IType i_, size_t j_shape_, size_t j_offset_ = 0)
    : s(s_),
      ija(IJA(s)),
      a(reinterpret_cast<YALE_STORAGE*>(s->src)->a),
      i(i_),
      k(ija[i]),
      k_end(ija[i+1]),
      j_offset(j_offset_),
      j_shape(j_shape_),
      diag(row_has_no_nd() || diag_is_first()),
      End(false),
      init(default_value(s))
  { }

  RowIterator(const RowIterator& rhs) : s(rhs.s), ija(rhs.ija), a(reinterpret_cast<YALE_STORAGE*>(s->src)->a), i(rhs.i), k(rhs.k), k_end(rhs.k_end), j_offset(rhs.j_offset), j_shape(rhs.j_shape), diag(rhs.diag), End(rhs.End), init(rhs.init) { }

  VALUE obj() const {
    return diag ? obj_at(s, i) : obj_at(s, k);
  }

  template <typename T>
  T cobj() const {
    if (typeid(T) == typeid(RubyObject)) return obj();
    return A<T>(s)[diag ? i : k];
  }

  inline IType proper_j() const {
    return diag ? i : ija[k];
  }

  inline IType offset_j() const {
    return proper_j() - j_offset;
  }

  inline size_t capacity() const {
    return reinterpret_cast<YALE_STORAGE*>(s->src)->capacity;
  }

  inline void vector_grow() {
    YALE_STORAGE* src = reinterpret_cast<YALE_STORAGE*>(s->src);
    nm::yale_storage::vector_grow(src);
    ija = reinterpret_cast<IType*>(src->ija);
    a   = src->a;
  }

  /* Returns true if an additional value is inserted, false if it goes on the diagonal */
  bool insert(IType j, VALUE v) {
    if (j == i) { // insert regardless on diagonal
      reinterpret_cast<VALUE*>(a)[j] = v;
      return false;

    } else {
      if (rb_funcall(v, rb_intern("!="), 1, init) == Qtrue) {
        if (k >= capacity()) {
          vector_grow();
        }
        reinterpret_cast<VALUE*>(a)[k] = v;
        ija[k] = j;
        k++;
        return true;
      }
      return false;
    }
  }

  void update_row_end() {
    ija[i+1] = k;
    k_end    = k;
  }

  /* Past the j_shape? */
  inline bool end() const {
    if (End)  return true;
    //if (diag) return i - j_offset >= j_shape;
    //else return k >= s->capacity || ija[k] - j_offset >= j_shape;
    return (int)(diag ? i : ija[k]) - (int)(j_offset) >= (int)(j_shape);
  }

  inline bool row_has_no_nd() const { return ija[i] == k_end; /* k_start == k_end */  }
  inline bool diag_is_first() const { return i < ija[ija[i]];  }
  inline bool diag_is_last() const  { return i > ija[k_end-1]; } // only works if !row_has_no_nd()
  inline bool k_is_last_nd() const  { return k == k_end-1;     }
  inline bool k_is_last() const     { return k_is_last_nd() && !diag_is_last(); }
  inline bool diag_is_ahead() const { return i > ija[k]; }
  inline bool row_has_diag() const  { return i < s->shape[1];  }
  inline bool diag_is_next() const  { // assumes we've already tested for diag, row_has_no_nd(), diag_is_first()
    if (i == ija[k]+1) return true; // definite next
    else if (k+1 < k_end && i >= ija[k+1]+1) return false; // at least one item before it
    else return true;
  }

  RowIterator& operator++() {
    if (diag) {                                             // we're at the diagonal
      if (row_has_no_nd() || diag_is_last()) End = true;    //  and there are no non-diagonals (or none still to visit)
      diag = false;
    } else if (!row_has_diag()) {                           // row has no diagonal entries
      if (row_has_no_nd() || k_is_last_nd()) End = true;    // row is totally empty, or we're at last entry
      else k++;                                             // still entries to visit
    } else { // not at diag but it exists somewhere in the row, and row has at least one nd entry
      if (diag_is_ahead()) { // diag is ahead
        if (k_is_last_nd()) diag = true; // diag is next and last
        else if (diag_is_next()) {       // diag is next and not last
          diag = true;
          k++;
        } else k++;                      // diag is not next
      } else {                           // diag is past
        if (k_is_last_nd()) End = true;  //   and we're at the end
        else k++;                        //   and we're not at the end
      }
    }

    return *this;
  }


  RowIterator operator++(int unused) {
    RowIterator x(*this);
    ++(*this);
    return x;
  }
};


// Helper function used only for the RETURN_SIZED_ENUMERATOR macro. Returns the length of
// the matrix's storage.
static VALUE nm_yale_stored_enumerator_length(VALUE nmatrix) {
  NM_CONSERVATIVE(nm_register_value(nmatrix));
  YALE_STORAGE* s   = NM_STORAGE_YALE(nmatrix);
  YALE_STORAGE* src = s->src == s ? s : reinterpret_cast<YALE_STORAGE*>(s->src);
  size_t ia_size    = src->shape[0];
  // FIXME: This needs to be corrected for slicing.
  size_t len = std::min( s->shape[0] + s->offset[0], s->shape[1] + s->offset[1] ) + nm_yale_storage_get_size(src) -  ia_size;
  NM_CONSERVATIVE(nm_unregister_value(nmatrix));
  return INT2FIX(len);
}


// Helper function used only for the RETURN_SIZED_ENUMERATOR macro. Returns the length of
// the matrix's storage.
static VALUE nm_yale_stored_nondiagonal_enumerator_length(VALUE nmatrix) {
  NM_CONSERVATIVE(nm_register_value(nmatrix));
  YALE_STORAGE* s = NM_STORAGE_YALE(nmatrix);
  if (s->src != s) s = reinterpret_cast<YALE_STORAGE*>(s->src);  // need to get the original storage shape

  size_t ia_size = s->shape[0];
  size_t len     = nm_yale_storage_get_size(NM_STORAGE_YALE(nmatrix)) - ia_size;
  NM_CONSERVATIVE(nm_unregister_value(nmatrix));
  return INT2FIX(len);
}

// Helper function for diagonal length.
static VALUE nm_yale_stored_diagonal_enumerator_length(VALUE nmatrix) {
  NM_CONSERVATIVE(nm_register_value(nmatrix));
  YALE_STORAGE* s = NM_STORAGE_YALE(nmatrix);
  size_t len = std::min( s->shape[0] + s->offset[0], s->shape[1] + s->offset[1] );
  NM_CONSERVATIVE(nm_unregister_value(nmatrix));
  return INT2FIX(len);
}


// Helper function for full enumerator length.
static VALUE nm_yale_enumerator_length(VALUE nmatrix) {
  NM_CONSERVATIVE(nm_register_value(nmatrix));
  YALE_STORAGE* s = NM_STORAGE_YALE(nmatrix);
  size_t len = s->shape[0] * s->shape[1];
  NM_CONSERVATIVE(nm_unregister_value(nmatrix));
  return INT2FIX(len);
}


/*
 * Map the stored values of a matrix in storage order.
 */
template <typename D>
static VALUE map_stored(VALUE self) {
  NM_CONSERVATIVE(nm_register_value(self));
  YALE_STORAGE* s = NM_STORAGE_YALE(self);
  YaleStorage<D> y(s);
  
  RETURN_SIZED_ENUMERATOR_PRE
  NM_CONSERVATIVE(nm_unregister_value(self));
  RETURN_SIZED_ENUMERATOR(self, 0, 0, nm_yale_stored_enumerator_length);

  YALE_STORAGE* r = y.template alloc_copy<nm::RubyObject, true>();
  nm_yale_storage_register(r);
  NMATRIX* m      = nm_create(nm::YALE_STORE, reinterpret_cast<STORAGE*>(r));
  VALUE to_return = Data_Wrap_Struct(CLASS_OF(self), nm_mark, nm_delete, m);
  nm_yale_storage_unregister(r);
  NM_CONSERVATIVE(nm_unregister_value(self));
  return to_return;
}


/*
 * map_stored which visits the stored entries of two matrices in order.
 */
template <typename LD, typename RD>
static VALUE map_merged_stored(VALUE left, VALUE right, VALUE init) {
  nm::YaleStorage<LD> l(NM_STORAGE_YALE(left));
  nm::YaleStorage<RD> r(NM_STORAGE_YALE(right));
  VALUE to_return = l.map_merged_stored(CLASS_OF(left), r, init);
  return to_return;
}


/*
 * Iterate over the stored entries in Yale (diagonal and non-diagonal non-zeros)
 */
template <typename DType>
static VALUE each_stored_with_indices(VALUE nm) {
  NM_CONSERVATIVE(nm_register_value(nm));
  YALE_STORAGE* s = NM_STORAGE_YALE(nm);
  YaleStorage<DType> y(s);

  // If we don't have a block, return an enumerator.
  RETURN_SIZED_ENUMERATOR_PRE
  NM_CONSERVATIVE(nm_unregister_value(nm));
  RETURN_SIZED_ENUMERATOR(nm, 0, 0, nm_yale_stored_enumerator_length);

  for (typename YaleStorage<DType>::const_stored_diagonal_iterator d = y.csdbegin(); d != y.csdend(); ++d) {
    rb_yield_values(3, ~d, d.rb_i(), d.rb_j());
  }

  for (typename YaleStorage<DType>::const_row_iterator it = y.cribegin(); it != y.criend(); ++it) {
    for (auto jt = it.ndbegin(); jt != it.ndend(); ++jt) {
      rb_yield_values(3, ~jt, it.rb_i(), jt.rb_j());
    }
  }

  NM_CONSERVATIVE(nm_unregister_value(nm));

  return nm;
}


/*
 * Iterate over the stored diagonal entries in Yale.
 */
template <typename DType>
static VALUE stored_diagonal_each_with_indices(VALUE nm) {
  NM_CONSERVATIVE(nm_register_value(nm));

  YALE_STORAGE* s = NM_STORAGE_YALE(nm);
  YaleStorage<DType> y(s);

  // If we don't have a block, return an enumerator.
  RETURN_SIZED_ENUMERATOR_PRE
  NM_CONSERVATIVE(nm_unregister_value(nm));
  RETURN_SIZED_ENUMERATOR(nm, 0, 0, nm_yale_stored_diagonal_length); // FIXME: need diagonal length
  
  for (typename YaleStorage<DType>::const_stored_diagonal_iterator d = y.csdbegin(); d != y.csdend(); ++d) {
    rb_yield_values(3, ~d, d.rb_i(), d.rb_j());
  }

  NM_CONSERVATIVE(nm_unregister_value(nm));

  return nm;
}


/*
 * Iterate over the stored diagonal entries in Yale.
 */
template <typename DType>
static VALUE stored_nondiagonal_each_with_indices(VALUE nm) {
  NM_CONSERVATIVE(nm_register_value(nm));

  YALE_STORAGE* s = NM_STORAGE_YALE(nm);
  YaleStorage<DType> y(s);

  // If we don't have a block, return an enumerator.
  RETURN_SIZED_ENUMERATOR_PRE
  NM_CONSERVATIVE(nm_unregister_value(nm));
  RETURN_SIZED_ENUMERATOR(nm, 0, 0, 0); // FIXME: need diagonal length

  for (typename YaleStorage<DType>::const_row_iterator it = y.cribegin(); it != y.criend(); ++it) {
    for (auto jt = it.ndbegin(); jt != it.ndend(); ++jt) {
      rb_yield_values(3, ~jt, it.rb_i(), jt.rb_j());
    }
  }

  NM_CONSERVATIVE(nm_unregister_value(nm));

  return nm;
}


/*
 * Iterate over the stored entries in Yale in order of i,j. Visits every diagonal entry, even if it's the default.
 */
template <typename DType>
static VALUE each_ordered_stored_with_indices(VALUE nm) {
  NM_CONSERVATIVE(nm_register_value(nm));

  YALE_STORAGE* s = NM_STORAGE_YALE(nm);
  YaleStorage<DType> y(s);

  // If we don't have a block, return an enumerator.
  RETURN_SIZED_ENUMERATOR_PRE
  NM_CONSERVATIVE(nm_unregister_value(nm));
  RETURN_SIZED_ENUMERATOR(nm, 0, 0, nm_yale_stored_enumerator_length);

  for (typename YaleStorage<DType>::const_row_iterator it = y.cribegin(); it != y.criend(); ++it) {
    for (auto jt = it.begin(); jt != it.end(); ++jt) {
      rb_yield_values(3, ~jt, it.rb_i(), jt.rb_j());
    }
  }

  NM_CONSERVATIVE(nm_unregister_value(nm));

  return nm;
}


template <typename DType>
static VALUE each_with_indices(VALUE nm) {
  NM_CONSERVATIVE(nm_register_value(nm));

  YALE_STORAGE* s = NM_STORAGE_YALE(nm);
  YaleStorage<DType> y(s);

  // If we don't have a block, return an enumerator.
  RETURN_SIZED_ENUMERATOR_PRE
  NM_CONSERVATIVE(nm_unregister_value(nm));
  RETURN_SIZED_ENUMERATOR(nm, 0, 0, nm_yale_enumerator_length);

  for (typename YaleStorage<DType>::const_iterator iter = y.cbegin(); iter != y.cend(); ++iter) {
    rb_yield_values(3, ~iter, iter.rb_i(), iter.rb_j());
  }

  NM_CONSERVATIVE(nm_unregister_value(nm));

  return nm;
}

template <typename D>
static bool is_pos_default_value(YALE_STORAGE* s, size_t apos) {
  YaleStorage<D> y(s);
  return y.is_pos_default_value(apos);
}


} // end of namespace nm::yale_storage


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

  // Expert recommendation. Eventually this should go in a separate gem, or at least a separate module.
  rb_define_method(cNMatrix_YaleFunctions, "yale_row_keys_intersection", (METHOD)nm_row_keys_intersection, 3);

  // Debugging functions.
  rb_define_method(cNMatrix_YaleFunctions, "yale_ija", (METHOD)nm_ija, -1);
  rb_define_method(cNMatrix_YaleFunctions, "yale_a", (METHOD)nm_a, -1);
  rb_define_method(cNMatrix_YaleFunctions, "yale_size", (METHOD)nm_size, 0);
  rb_define_method(cNMatrix_YaleFunctions, "yale_ia", (METHOD)nm_ia, 0);
  rb_define_method(cNMatrix_YaleFunctions, "yale_ja", (METHOD)nm_ja, 0);
  rb_define_method(cNMatrix_YaleFunctions, "yale_d", (METHOD)nm_d, -1);
  rb_define_method(cNMatrix_YaleFunctions, "yale_lu", (METHOD)nm_lu, 0);

  rb_define_method(cNMatrix_YaleFunctions, "yale_nd_row", (METHOD)nm_nd_row, -1);

  rb_define_const(cNMatrix_YaleFunctions, "YALE_GROWTH_CONSTANT", rb_float_new(nm::yale_storage::GROWTH_CONSTANT));

  // This is so the user can easily check the IType size, mostly for debugging.
  size_t itype_size        = sizeof(IType);
  VALUE itype_dtype;
  if (itype_size == sizeof(uint64_t)) {
    itype_dtype = ID2SYM(rb_intern("int64"));
  } else if (itype_size == sizeof(uint32_t)) {
    itype_dtype = ID2SYM(rb_intern("int32"));
  } else if (itype_size == sizeof(uint16_t)) {
    itype_dtype = ID2SYM(rb_intern("int16"));
  } else {
    rb_raise(rb_eStandardError, "unhandled length for sizeof(IType): %lu; note that IType is probably defined as size_t", sizeof(IType));
  }
  rb_define_const(cNMatrix, "INDEX_DTYPE", itype_dtype);
}


/////////////////
// C ACCESSORS //
/////////////////


/* C interface for NMatrix#each_with_indices (Yale) */
VALUE nm_yale_each_with_indices(VALUE nmatrix) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::each_with_indices, VALUE, VALUE)

  return ttable[ NM_DTYPE(nmatrix) ](nmatrix);
}


/* C interface for NMatrix#each_stored_with_indices (Yale) */
VALUE nm_yale_each_stored_with_indices(VALUE nmatrix) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::each_stored_with_indices, VALUE, VALUE)

  return ttable[ NM_DTYPE(nmatrix) ](nmatrix);
}


/* Iterate along stored diagonal (not actual diagonal!) */
VALUE nm_yale_stored_diagonal_each_with_indices(VALUE nmatrix) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::stored_diagonal_each_with_indices, VALUE, VALUE)

  return ttable[ NM_DTYPE(nmatrix) ](nmatrix);
}

/* Iterate through stored nondiagonal (not actual diagonal!) */
VALUE nm_yale_stored_nondiagonal_each_with_indices(VALUE nmatrix) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::stored_nondiagonal_each_with_indices, VALUE, VALUE)

  return ttable[ NM_DTYPE(nmatrix) ](nmatrix);
}


/* C interface for NMatrix#each_ordered_stored_with_indices (Yale) */
VALUE nm_yale_each_ordered_stored_with_indices(VALUE nmatrix) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::each_ordered_stored_with_indices, VALUE, VALUE)

  return ttable[ NM_DTYPE(nmatrix) ](nmatrix);
}



/*
 * C accessor for inserting some value in a matrix (or replacing an existing cell).
 */
void nm_yale_storage_set(VALUE left, SLICE* slice, VALUE right) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::set, void, VALUE left, SLICE* slice, VALUE right);

  ttable[NM_DTYPE(left)](left, slice, right);
}


/*
 * Determine the number of non-diagonal non-zeros in a not-yet-created copy of a slice or matrix.
 */
static size_t yale_count_slice_copy_ndnz(const YALE_STORAGE* s, size_t* offset, size_t* shape) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::count_slice_copy_ndnz, size_t, const YALE_STORAGE*, size_t*, size_t*)

  return ttable[s->dtype](s, offset, shape);
}


/*
 * C accessor for yale_storage::get, which returns a slice of YALE_STORAGE object by copy
 *
 * Slicing-related.
 */
void* nm_yale_storage_get(const STORAGE* storage, SLICE* slice) {
  YALE_STORAGE* casted_storage = (YALE_STORAGE*)storage;

  if (slice->single) {
    NAMED_DTYPE_TEMPLATE_TABLE(elem_copy_table,  nm::yale_storage::get_single, void*, YALE_STORAGE*, SLICE*)

    return elem_copy_table[casted_storage->dtype](casted_storage, slice);
  } else {
    nm_yale_storage_register(casted_storage);
    //return reinterpret_cast<void*>(nm::YaleStorage<nm::dtype_enum_T<storage->dtype>::type>(casted_storage).alloc_ref(slice));
    NAMED_DTYPE_TEMPLATE_TABLE(ref_table, nm::yale_storage::ref, YALE_STORAGE*, YALE_STORAGE* storage, SLICE* slice)

    YALE_STORAGE* ref = ref_table[casted_storage->dtype](casted_storage, slice);

    NAMED_LR_DTYPE_TEMPLATE_TABLE(slice_copy_table, nm::yale_storage::slice_copy, YALE_STORAGE*, YALE_STORAGE*)

    YALE_STORAGE* ns = slice_copy_table[casted_storage->dtype][casted_storage->dtype](ref);

    NM_FREE(ref);

    nm_yale_storage_unregister(casted_storage);

    return ns;
  }
}

/*
 * C accessor for yale_storage::vector_insert
 */
static char nm_yale_storage_vector_insert(YALE_STORAGE* s, size_t pos, size_t* js, void* vals, size_t n, bool struct_only, nm::dtype_t dtype) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::vector_insert, char, YALE_STORAGE*, size_t, size_t*, void*, size_t, bool);

  return ttable[dtype](s, pos, js, vals, n, struct_only);
}

/*
 * C accessor for yale_storage::increment_ia_after, typically called after ::vector_insert
 */
static void nm_yale_storage_increment_ia_after(YALE_STORAGE* s, size_t ija_size, size_t i, long n) {
  nm::yale_storage::increment_ia_after(s, ija_size, i, n);
}


/*
 * C accessor for yale_storage::ref, which returns either a pointer to the correct location in a YALE_STORAGE object
 * for some set of coordinates, or a pointer to a single element.
 */
void* nm_yale_storage_ref(const STORAGE* storage, SLICE* slice) {
  YALE_STORAGE* casted_storage = (YALE_STORAGE*)storage;

  if (slice->single) {
    //return reinterpret_cast<void*>(nm::YaleStorage<nm::dtype_enum_T<storage->dtype>::type>(casted_storage).get_single_p(slice));
    NAMED_DTYPE_TEMPLATE_TABLE(elem_copy_table,  nm::yale_storage::get_single, void*, YALE_STORAGE*, SLICE*)
    return elem_copy_table[casted_storage->dtype](casted_storage, slice);
  } else {
    //return reinterpret_cast<void*>(nm::YaleStorage<nm::dtype_enum_T<storage->dtype>::type>(casted_storage).alloc_ref(slice));
    NAMED_DTYPE_TEMPLATE_TABLE(ref_table, nm::yale_storage::ref, YALE_STORAGE*, YALE_STORAGE* storage, SLICE* slice)
    return reinterpret_cast<void*>(ref_table[casted_storage->dtype](casted_storage, slice));

  }
}


/*
 * C accessor for determining whether two YALE_STORAGE objects have the same contents.
 */
bool nm_yale_storage_eqeq(const STORAGE* left, const STORAGE* right) {
  NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::eqeq, bool, const YALE_STORAGE* left, const YALE_STORAGE* right);

  const YALE_STORAGE* casted_left = reinterpret_cast<const YALE_STORAGE*>(left);

  return ttable[casted_left->dtype][right->dtype](casted_left, (const YALE_STORAGE*)right);
}


/*
 * Copy constructor for changing dtypes. (C accessor)
 */
STORAGE* nm_yale_storage_cast_copy(const STORAGE* rhs, nm::dtype_t new_dtype, void* dummy) {
  NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::cast_copy, YALE_STORAGE*, const YALE_STORAGE* rhs);

  const YALE_STORAGE* casted_rhs = reinterpret_cast<const YALE_STORAGE*>(rhs);
  //return reinterpret_cast<STORAGE*>(nm::YaleStorage<nm::dtype_enum_T< rhs->dtype >::type>(rhs).alloc_copy<nm::dtype_enum_T< new_dtype >::type>());
  return (STORAGE*)ttable[new_dtype][casted_rhs->dtype](casted_rhs);
}


/*
 * Returns size of Yale storage as a size_t (no matter what the itype is). (C accessor)
 */
size_t nm_yale_storage_get_size(const YALE_STORAGE* storage) {
  return nm::yale_storage::get_size(storage);
}



/*
 * Return a pointer to the matrix's default value entry.
 */
static void* default_value_ptr(const YALE_STORAGE* s) {
  return reinterpret_cast<void*>(reinterpret_cast<char*>(((YALE_STORAGE*)(s->src))->a) + (((YALE_STORAGE*)(s->src))->shape[0] * DTYPE_SIZES[s->dtype]));
}

/*
 * Return the Ruby object at a given location in storage.
 */
static VALUE obj_at(YALE_STORAGE* s, size_t k) {
  if (s->dtype == nm::RUBYOBJ)  return reinterpret_cast<VALUE*>(((YALE_STORAGE*)(s->src))->a)[k];
  else  return rubyobj_from_cval(reinterpret_cast<void*>(reinterpret_cast<char*>(((YALE_STORAGE*)(s->src))->a) + k * DTYPE_SIZES[s->dtype]), s->dtype).rval;
}


/*
 * Return the matrix's default value as a Ruby VALUE.
 */
static VALUE default_value(const YALE_STORAGE* s) {
  if (s->dtype == nm::RUBYOBJ) return *reinterpret_cast<VALUE*>(default_value_ptr(s));
  else return rubyobj_from_cval(default_value_ptr(s), s->dtype).rval;
}


/*
 * Check to see if a default value is some form of zero. Easy for non-Ruby object matrices, which should always be 0.
 */
static bool default_value_is_numeric_zero(const YALE_STORAGE* s) {
  return rb_funcall(default_value(s), rb_intern("=="), 1, INT2FIX(0)) == Qtrue;
}



/*
 * Transposing copy constructor.
 */
STORAGE* nm_yale_storage_copy_transposed(const STORAGE* rhs_base) {
  YALE_STORAGE* rhs = (YALE_STORAGE*)rhs_base;
  NAMED_DTYPE_TEMPLATE_TABLE(transp, nm::yale_storage::copy_transposed, YALE_STORAGE*, YALE_STORAGE*)
  return (STORAGE*)(transp[rhs->dtype](rhs));
}

/*
 * C accessor for multiplying two YALE_STORAGE matrices, which have already been casted to the same dtype.
 *
 * FIXME: There should be some mathematical way to determine the worst-case IType based on the input ITypes. Right now
 * it just uses the default.
 */
STORAGE* nm_yale_storage_matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector) {
  DTYPE_TEMPLATE_TABLE(nm::yale_storage::matrix_multiply, STORAGE*, const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector);

  YALE_STORAGE* left = reinterpret_cast<YALE_STORAGE*>(casted_storage.left);
  YALE_STORAGE* right = reinterpret_cast<YALE_STORAGE*>(casted_storage.right);

  if (!default_value_is_numeric_zero(left) || !default_value_is_numeric_zero(right)) {
    rb_raise(rb_eNotImpError, "matrix default value must be some form of zero (not false or nil) for multiplication");
    return NULL;
  }

  return ttable[left->dtype](casted_storage, resulting_shape, vector);
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

YALE_STORAGE* nm_yale_storage_create(nm::dtype_t dtype, size_t* shape, size_t dim, size_t init_capacity) {
  if (dim != 2) {
    rb_raise(nm_eStorageTypeError, "yale supports only 2-dimensional matrices");
  }
  DTYPE_OBJECT_STATIC_TABLE(nm::YaleStorage, create, YALE_STORAGE*, size_t* shape, size_t init_capacity)
  return ttable[dtype](shape, init_capacity);
}

/*
 * Destructor for yale storage (C-accessible).
 */
void nm_yale_storage_delete(STORAGE* s) {
  if (s) {
    YALE_STORAGE* storage = (YALE_STORAGE*)s;
    if (storage->count-- == 1) {
      NM_FREE(storage->shape);
      NM_FREE(storage->offset);
      NM_FREE(storage->ija);
      NM_FREE(storage->a);
      NM_FREE(storage);
    }
  }
}

/*
 * Destructor for the yale storage ref
 */
void nm_yale_storage_delete_ref(STORAGE* s) {
  if (s) {
    YALE_STORAGE* storage = (YALE_STORAGE*)s;
    nm_yale_storage_delete( reinterpret_cast<STORAGE*>(storage->src) );
    NM_FREE(storage->shape);
    NM_FREE(storage->offset);
    NM_FREE(s);
  }
}

/*
 * C accessor for yale_storage::init, a templated function.
 *
 * Initializes the IJA vector of the YALE_STORAGE matrix.
 */
void nm_yale_storage_init(YALE_STORAGE* s, void* init_val) {
  DTYPE_TEMPLATE_TABLE(nm::yale_storage::init, void, YALE_STORAGE*, void*);

  ttable[s->dtype](s, init_val);
}


/*
 * Ruby GC mark function for YALE_STORAGE. C accessible.
 */
void nm_yale_storage_mark(STORAGE* storage_base) {
  YALE_STORAGE* storage = (YALE_STORAGE*)storage_base;

  if (storage && storage->dtype == nm::RUBYOBJ) {

    VALUE* a = (VALUE*)(storage->a);
    rb_gc_mark_locations(a, &(a[storage->capacity-1]));
  }
}

void nm_yale_storage_register_a(void* a, size_t size) {
  nm_register_values(reinterpret_cast<VALUE*>(a), size);
}

void nm_yale_storage_unregister_a(void* a, size_t size) {
  nm_unregister_values(reinterpret_cast<VALUE*>(a), size);
}

void nm_yale_storage_register(const STORAGE* s) {
  const YALE_STORAGE* y = reinterpret_cast<const YALE_STORAGE*>(s);
  if (y->dtype == nm::RUBYOBJ) {
    nm_register_values(reinterpret_cast<VALUE*>(y->a), nm::yale_storage::get_size(y));
  }
}

void nm_yale_storage_unregister(const STORAGE* s) {
  const YALE_STORAGE* y = reinterpret_cast<const YALE_STORAGE*>(s);
  if (y->dtype == nm::RUBYOBJ) {
    nm_unregister_values(reinterpret_cast<VALUE*>(y->a), nm::yale_storage::get_size(y));
  }
}

/*
 * Allocates and initializes the basic struct (but not the IJA or A vectors).
 *
 * This function is ONLY used when creating from old yale.
 */
static YALE_STORAGE* alloc(nm::dtype_t dtype, size_t* shape, size_t dim) {
  YALE_STORAGE* s;

  s = NM_ALLOC( YALE_STORAGE );

  s->ndnz        = 0;
  s->dtype       = dtype;
  s->shape       = shape;
  s->offset      = NM_ALLOC_N(size_t, dim);
  for (size_t i = 0; i < dim; ++i)
    s->offset[i] = 0;
  s->dim         = dim;
  s->src         = reinterpret_cast<STORAGE*>(s);
  s->count       = 1;

  return s;
}

YALE_STORAGE* nm_yale_storage_create_from_old_yale(nm::dtype_t dtype, size_t* shape, char* ia, char* ja, char* a, nm::dtype_t from_dtype) {
  NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::create_from_old_yale, YALE_STORAGE*, nm::dtype_t dtype, size_t* shape, char* r_ia, char* r_ja, char* r_a);

  return ttable[dtype][from_dtype](dtype, shape, ia, ja, a);

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
  YALE_STORAGE* s = (YALE_STORAGE*)(NM_SRC(self));
  VALUE to_return = INT2FIX(nm::yale_storage::IJA(s)[s->shape[0]]);
  return to_return;
}


/*
 * Determine if some pos in the diagonal is the default. No bounds checking!
 */
static bool is_pos_default_value(YALE_STORAGE* s, size_t apos) {
  DTYPE_TEMPLATE_TABLE(nm::yale_storage::is_pos_default_value, bool, YALE_STORAGE*, size_t)
  return ttable[s->dtype](s, apos);
}


/*
 * call-seq:
 *     yale_row_keys_intersection(i, m2, i2) -> Array
 *
 * This function is experimental.
 *
 * It finds the intersection of row i of the current matrix with row i2 of matrix m2.
 * Both matrices must be Yale. They may not be slices.
 *
 * Only checks the stored indices; does not care about matrix default value.
 */
static VALUE nm_row_keys_intersection(VALUE m1, VALUE ii1, VALUE m2, VALUE ii2) {
  
  NM_CONSERVATIVE(nm_register_value(m1));
  NM_CONSERVATIVE(nm_register_value(m2));

  if (NM_SRC(m1) != NM_STORAGE(m1) || NM_SRC(m2) != NM_STORAGE(m2)) {
    NM_CONSERVATIVE(nm_unregister_value(m2));
    NM_CONSERVATIVE(nm_unregister_value(m1));
    rb_raise(rb_eNotImpError, "must be called on a real matrix and not a slice");
  }

  size_t i1 = FIX2INT(ii1),
         i2 = FIX2INT(ii2);

  YALE_STORAGE *s   = NM_STORAGE_YALE(m1),
               *t   = NM_STORAGE_YALE(m2);

  size_t pos1 = s->ija[i1],
         pos2 = t->ija[i2];

  size_t nextpos1 = s->ija[i1+1],
         nextpos2 = t->ija[i2+1];

  size_t diff1 = nextpos1 - pos1,
         diff2 = nextpos2 - pos2;

  // Does the diagonal have a nonzero in it?
  bool diag1 = i1 < s->shape[0] && !is_pos_default_value(s, i1),
       diag2 = i2 < t->shape[0] && !is_pos_default_value(t, i2);

  // Reserve max(diff1,diff2) space -- that's the max intersection possible.
  VALUE ret = rb_ary_new2(std::max(diff1,diff2)+1);
  nm_register_value(ret);

  // Handle once the special case where both have the diagonal in exactly
  // the same place.
  if (diag1 && diag2 && i1 == i2) {
    rb_ary_push(ret, INT2FIX(i1));
    diag1 = false; diag2 = false; // no need to deal with diagonals anymore.
  }

  // Now find the intersection.
  size_t idx1 = pos1, idx2 = pos2;
  while (idx1 < nextpos1 && idx2 < nextpos2) {
    if (s->ija[idx1] == t->ija[idx2]) {
      rb_ary_push(ret, INT2FIX(s->ija[idx1]));
      ++idx1; ++idx2;
    } else if (diag1 && i1 == t->ija[idx2]) {
      rb_ary_push(ret, INT2FIX(i1));
      diag1 = false;
      ++idx2;
    } else if (diag2 && i2 == s->ija[idx1]) {
      rb_ary_push(ret, INT2FIX(i2));
      diag2 = false;
      ++idx1;
    } else if (s->ija[idx1] < t->ija[idx2]) {
      ++idx1;
    } else { // s->ija[idx1] > t->ija[idx2]
      ++idx2;
    }
  }

  // Past the end of row i2's stored entries; need to try to find diagonal
  if (diag2 && idx1 < nextpos1) {
    idx1 = nm::yale_storage::binary_search_left_boundary(s, idx1, nextpos1, i2);
    if (s->ija[idx1] == i2) rb_ary_push(ret, INT2FIX(i2));
  }

  // Find the diagonal, if possible, in the other one.
  if (diag1 && idx2 < nextpos2) {
    idx2 = nm::yale_storage::binary_search_left_boundary(t, idx2, nextpos2, i1);
    if (t->ija[idx2] == i1) rb_ary_push(ret, INT2FIX(i1));
  }

  nm_unregister_value(ret);
  NM_CONSERVATIVE(nm_unregister_value(m1));
  NM_CONSERVATIVE(nm_unregister_value(m2));

  return ret;
}


/*
 * call-seq:
 *     yale_a -> Array
 *     yale_d(index) -> ...
 *
 * Get the A array of a Yale matrix (which stores the diagonal and the LU portions of the matrix).
 */
static VALUE nm_a(int argc, VALUE* argv, VALUE self) {
  NM_CONSERVATIVE(nm_register_value(self));

  VALUE idx;
  rb_scan_args(argc, argv, "01", &idx);
  NM_CONSERVATIVE(nm_register_value(idx));

  YALE_STORAGE* s = reinterpret_cast<YALE_STORAGE*>(NM_SRC(self));
  size_t size = nm_yale_storage_get_size(s);

  if (idx == Qnil) {

    VALUE* vals = NM_ALLOCA_N(VALUE, size);

    nm_register_values(vals, size);
    
    if (NM_DTYPE(self) == nm::RUBYOBJ) {
      for (size_t i = 0; i < size; ++i) {
        vals[i] = reinterpret_cast<VALUE*>(s->a)[i];
      }
    } else {
      for (size_t i = 0; i < size; ++i) {
        vals[i] = rubyobj_from_cval((char*)(s->a) + DTYPE_SIZES[s->dtype]*i, s->dtype).rval;
      }
    }
    VALUE ary = rb_ary_new4(size, vals);

    for (size_t i = size; i < s->capacity; ++i)
      rb_ary_push(ary, Qnil);

    nm_unregister_values(vals, size);
    NM_CONSERVATIVE(nm_unregister_value(idx));
    NM_CONSERVATIVE(nm_unregister_value(self));
    return ary;
  } else {
    size_t index = FIX2INT(idx);
    NM_CONSERVATIVE(nm_unregister_value(idx));
    NM_CONSERVATIVE(nm_unregister_value(self));
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
  NM_CONSERVATIVE(nm_register_value(self));
  VALUE idx;
  rb_scan_args(argc, argv, "01", &idx);
  NM_CONSERVATIVE(nm_register_value(idx));

  YALE_STORAGE* s = reinterpret_cast<YALE_STORAGE*>(NM_SRC(self));

  if (idx == Qnil) {
    VALUE* vals = NM_ALLOCA_N(VALUE, s->shape[0]);

    nm_register_values(vals, s->shape[0]);

    if (NM_DTYPE(self) == nm::RUBYOBJ) {
      for (size_t i = 0; i < s->shape[0]; ++i) {
        vals[i] = reinterpret_cast<VALUE*>(s->a)[i];
      }
    } else {
      for (size_t i = 0; i < s->shape[0]; ++i) {
        vals[i] = rubyobj_from_cval((char*)(s->a) + DTYPE_SIZES[s->dtype]*i, s->dtype).rval;
      }
    }
    nm_unregister_values(vals, s->shape[0]);
    NM_CONSERVATIVE(nm_unregister_value(idx));
    NM_CONSERVATIVE(nm_unregister_value(self));

    return rb_ary_new4(s->shape[0], vals);
  } else {
    size_t index = FIX2INT(idx);
    NM_CONSERVATIVE(nm_unregister_value(idx));
    NM_CONSERVATIVE(nm_unregister_value(self));
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
  NM_CONSERVATIVE(nm_register_value(self));

  YALE_STORAGE* s = reinterpret_cast<YALE_STORAGE*>(NM_SRC(self));

  size_t size = nm_yale_storage_get_size(s);

  VALUE* vals = NM_ALLOCA_N(VALUE, size - s->shape[0] - 1);

  nm_register_values(vals, size - s->shape[0] - 1);

  if (NM_DTYPE(self) == nm::RUBYOBJ) {
    for (size_t i = 0; i < size - s->shape[0] - 1; ++i) {
      vals[i] = reinterpret_cast<VALUE*>(s->a)[s->shape[0] + 1 + i];
    }
  } else {
    for (size_t i = 0; i < size - s->shape[0] - 1; ++i) {
      vals[i] = rubyobj_from_cval((char*)(s->a) + DTYPE_SIZES[s->dtype]*(s->shape[0] + 1 + i), s->dtype).rval;
    }
  }

  VALUE ary = rb_ary_new4(size - s->shape[0] - 1, vals);

  for (size_t i = size; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  nm_unregister_values(vals, size - s->shape[0] - 1);
  NM_CONSERVATIVE(nm_unregister_value(self));

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
  NM_CONSERVATIVE(nm_register_value(self));

  YALE_STORAGE* s = reinterpret_cast<YALE_STORAGE*>(NM_SRC(self));

  VALUE* vals = NM_ALLOCA_N(VALUE, s->shape[0] + 1);

  for (size_t i = 0; i < s->shape[0] + 1; ++i) {
    vals[i] = INT2FIX(s->ija[i]);
  }

  NM_CONSERVATIVE(nm_unregister_value(self)); 

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

  NM_CONSERVATIVE(nm_register_value(self));

  YALE_STORAGE* s = reinterpret_cast<YALE_STORAGE*>(NM_SRC(self));

  size_t size = nm_yale_storage_get_size(s);

  VALUE* vals = NM_ALLOCA_N(VALUE, size - s->shape[0] - 1);

  nm_register_values(vals, size - s->shape[0] - 1);

  for (size_t i = 0; i < size - s->shape[0] - 1; ++i) {
    vals[i] = INT2FIX(s->ija[s->shape[0] + 1 + i]);
  }

  VALUE ary = rb_ary_new4(size - s->shape[0] - 1, vals);

  for (size_t i = size; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  nm_unregister_values(vals, size - s->shape[0] - 1);
  NM_CONSERVATIVE(nm_unregister_value(self));

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
  NM_CONSERVATIVE(nm_register_value(self));

  VALUE idx;
  rb_scan_args(argc, argv, "01", &idx);
  NM_CONSERVATIVE(nm_register_value(idx));

  YALE_STORAGE* s = reinterpret_cast<YALE_STORAGE*>(NM_SRC(self));
  size_t size = nm_yale_storage_get_size(s);

  if (idx == Qnil) {

    VALUE* vals = NM_ALLOCA_N(VALUE, size);

    nm_register_values(vals, size);

    for (size_t i = 0; i < size; ++i) {
      vals[i] = INT2FIX(s->ija[i]);
    }

   VALUE ary = rb_ary_new4(size, vals);

    for (size_t i = size; i < s->capacity; ++i)
      rb_ary_push(ary, Qnil);

    nm_unregister_values(vals, size);
    NM_CONSERVATIVE(nm_unregister_value(idx));
    NM_CONSERVATIVE(nm_unregister_value(self));

    return ary;

  } else {
    size_t index = FIX2INT(idx);
    if (index >= size) rb_raise(rb_eRangeError, "out of range");
    NM_CONSERVATIVE(nm_unregister_value(self));
    NM_CONSERVATIVE(nm_unregister_value(idx));
    return INT2FIX(s->ija[index]);
  }
}


/*
 * call-seq:
 *     yale_nd_row -> ...
 *
 * This function gets the non-diagonal contents of a Yale matrix row.
 * The first argument should be the row index. The optional second argument may be :hash or :keys, but defaults
 * to :hash. If :keys is given, it will only return the Hash keys (the column indices).
 *
 * This function is meant to accomplish its purpose as efficiently as possible. It does not check for appropriate
 * range.
 */
static VALUE nm_nd_row(int argc, VALUE* argv, VALUE self) {

  NM_CONSERVATIVE(nm_register_value(self));
  
  if (NM_SRC(self) != NM_STORAGE(self)) {
    NM_CONSERVATIVE(nm_unregister_value(self));
    rb_raise(rb_eNotImpError, "must be called on a real matrix and not a slice");
  }  

  VALUE i_, as;
  rb_scan_args(argc, argv, "11", &i_, &as);
  NM_CONSERVATIVE(nm_register_value(as));
  NM_CONSERVATIVE(nm_register_value(i_));

  bool keys = false;
  if (as != Qnil && rb_to_id(as) != nm_rb_hash) keys = true;

  size_t i = FIX2INT(i_);

  YALE_STORAGE* s   = NM_STORAGE_YALE(self);
  //nm::dtype_t dtype = NM_DTYPE(self);

  if (i >= s->shape[0]) {
    NM_CONSERVATIVE(nm_unregister_value(self));
    NM_CONSERVATIVE(nm_unregister_value(as));
    NM_CONSERVATIVE(nm_unregister_value(i_));
    rb_raise(rb_eRangeError, "out of range (%lu >= %lu)", i, s->shape[0]);
  }

  size_t pos = s->ija[i];
  size_t nextpos = s->ija[i+1];
  size_t diff = nextpos - pos;

  VALUE ret;
  if (keys) {
    ret = rb_ary_new3(diff);

    for (size_t idx = pos; idx < nextpos; ++idx) {
      rb_ary_store(ret, idx - pos, INT2FIX(s->ija[idx]));
    }

  } else {
    ret = rb_hash_new();

    for (size_t idx = pos; idx < nextpos; ++idx) {
      rb_hash_aset(ret, INT2FIX(s->ija[idx]), rubyobj_from_cval((char*)(s->a) + DTYPE_SIZES[s->dtype]*idx, s->dtype).rval);
    }
  }
  NM_CONSERVATIVE(nm_unregister_value(as));
  NM_CONSERVATIVE(nm_unregister_value(i_));
  NM_CONSERVATIVE(nm_unregister_value(self));
  return ret;
}

/*
 * call-seq:
 *     yale_vector_set(i, column_index_array, cell_contents_array, pos) -> Fixnum
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
 *    m.yale_vector_set(3, [0,3,4], [1,1,1], 15)
 *
 * The example above inserts the values 1, 1, and 1 in columns 0, 3, and 4, assumed to be located at position 15 (which
 * corresponds to row 3).
 *
 * Example:
 *    next = m.yale_vector_set(3, [0,3,4], [1,1,1])
 *
 * This example determines that i=3 is at position 15 automatically. The value returned, next, is the position where the
 * next value(s) should be inserted.
 */
VALUE nm_vector_set(int argc, VALUE* argv, VALUE self) { //, VALUE i_, VALUE jv, VALUE vv, VALUE pos_) {

  NM_CONSERVATIVE(nm_register_value(self));

  if (NM_SRC(self) != NM_STORAGE(self)) {
    NM_CONSERVATIVE(nm_unregister_value(self));
    rb_raise(rb_eNotImpError, "must be called on a real matrix and not a slice");
  }

  // i, jv, vv are mandatory; pos is optional; thus "31"
  VALUE i_, jv, vv, pos_;
  rb_scan_args(argc, argv, "31", &i_, &jv, &vv, &pos_);
  NM_CONSERVATIVE(nm_register_value(i_));
  NM_CONSERVATIVE(nm_register_value(jv));
  NM_CONSERVATIVE(nm_register_value(vv));
  NM_CONSERVATIVE(nm_register_value(pos_));

  size_t len   = RARRAY_LEN(jv); // need length in order to read the arrays in
  size_t vvlen = RARRAY_LEN(vv);

  if (len != vvlen) {
    NM_CONSERVATIVE(nm_unregister_value(pos_));
    NM_CONSERVATIVE(nm_unregister_value(vv));
    NM_CONSERVATIVE(nm_unregister_value(jv));
    NM_CONSERVATIVE(nm_unregister_value(i_));
    NM_CONSERVATIVE(nm_unregister_value(self));
    rb_raise(rb_eArgError, "lengths must match between j array (%lu) and value array (%lu)", len, vvlen);
  }

  YALE_STORAGE* s   = NM_STORAGE_YALE(self);
  nm::dtype_t dtype = NM_DTYPE(self);

  size_t i   = FIX2INT(i_);    // get the row
  size_t pos = s->ija[i];

  // Allocate the j array and the values array
  size_t* j  = NM_ALLOCA_N(size_t, len);
  void* vals = NM_ALLOCA_N(char, DTYPE_SIZES[dtype] * len);
  if (dtype == nm::RUBYOBJ){
    nm_register_values(reinterpret_cast<VALUE*>(vals), len);
  }

  // Copy array contents
  for (size_t idx = 0; idx < len; ++idx) {
    j[idx] = FIX2INT(rb_ary_entry(jv, idx));
    rubyval_to_cval(rb_ary_entry(vv, idx), dtype, (char*)vals + idx * DTYPE_SIZES[dtype]);
  }

  nm_yale_storage_vector_insert(s, pos, j, vals, len, false, dtype);
  nm_yale_storage_increment_ia_after(s, s->shape[0], i, len);
  s->ndnz += len;

  if (dtype == nm::RUBYOBJ){
    nm_unregister_values(reinterpret_cast<VALUE*>(vals), len);
  }

  NM_CONSERVATIVE(nm_unregister_value(pos_));
  NM_CONSERVATIVE(nm_unregister_value(vv));
  NM_CONSERVATIVE(nm_unregister_value(jv));
  NM_CONSERVATIVE(nm_unregister_value(i_));
  NM_CONSERVATIVE(nm_unregister_value(self));

  // Return the updated position
  pos += len;
  return INT2FIX(pos);
}




/*
 * call-seq:
 *     __yale_default_value__ -> ...
 *
 * Get the default_value property from a yale matrix.
 */
VALUE nm_yale_default_value(VALUE self) {
  VALUE to_return = default_value(NM_STORAGE_YALE(self));
  return to_return;
}


/*
 * call-seq:
 *     __yale_map_merged_stored__(right) -> Enumerator
 *
 * A map operation on two Yale matrices which only iterates across the stored indices.
 */
VALUE nm_yale_map_merged_stored(VALUE left, VALUE right, VALUE init) {
  NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::map_merged_stored, VALUE, VALUE, VALUE, VALUE)
  return ttable[NM_DTYPE(left)][NM_DTYPE(right)](left, right, init);
  //return nm::yale_storage::map_merged_stored(left, right, init);
}


/*
 * call-seq:
 *     __yale_map_stored__ -> Enumerator
 *
 * A map operation on two Yale matrices which only iterates across the stored indices.
 */
VALUE nm_yale_map_stored(VALUE self) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::yale_storage::map_stored, VALUE, VALUE)
  return ttable[NM_DTYPE(self)](self);
}

} // end of extern "C" block
