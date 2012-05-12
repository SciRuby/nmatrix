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

#ifndef YALE_C
# define YALE_C

#include <ruby.h> // mostly for exceptions

#include "nmatrix.h"

extern nm_eqeq_t ElemEqEq;
extern VALUE nm_eStorageTypeError;


extern const char *nm_dtypestring[];


void print_vectors(YALE_STORAGE* s) {
  size_t i;
  fprintf(stderr, "------------------------------\n");
  fprintf(stderr, "dtype:%s\tshape:%dx%d\tndnz:%d\tcapacity:%d\tindex_dtype:%s\n", nm_dtypestring[s->dtype], s->shape[0], s->shape[1], s->ndnz, s->capacity, nm_dtypestring[s->index_dtype]);


  if (s->capacity > 60) rb_raise(rb_eArgError, "overflow in print_vectors; cannot handle that large of a vector");
  // print indices

  fprintf(stderr, "i:\t");
  for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5lu ", (unsigned long)i);

  fprintf(stderr, "\nija:\t");
  if (YALE_MAX_SIZE(s) < UINT8_MAX)
    for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5u ", *(u_int8_t*)YALE_IJA(s,nm_sizeof[s->index_dtype],i));
  else if (YALE_MAX_SIZE(s) < UINT16_MAX)
    for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5u ", *(u_int16_t*)YALE_IJA(s,nm_sizeof[s->index_dtype],i));
  else if (YALE_MAX_SIZE(s) < UINT32_MAX)
    for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5u ", *(u_int32_t*)YALE_IJA(s,nm_sizeof[s->index_dtype],i));
  else
    for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5llu ", *(u_int64_t*)YALE_IJA(s,nm_sizeof[s->index_dtype],i));
  fprintf(stderr, "\n");

  // print values
  fprintf(stderr, "a:\t");
  for (i = 0; i < s->capacity; ++i)
    fprintf(stderr, "%-*.3g ", 5, *(double*)((char*)(s->a) + nm_sizeof[NM_FLOAT64]*i));
  fprintf(stderr, "\n");

  fprintf(stderr, "------------------------------\n");
}


// Determine the index dtype (which will be used for the ija vector). This is determined by matrix shape, not IJA/A vector capacity.
// Note that it's MAX-2 because UINTX_MAX and UINTX_MAX-1 are both reserved for sparse matrix multiplication.
int8_t yale_index_dtype(YALE_STORAGE* s) {
  if (YALE_MAX_SIZE(s) < UINT8_MAX-2) return NM_INT8;
  else if (YALE_MAX_SIZE(s) < UINT16_MAX-2) return NM_INT16;
  else if (YALE_MAX_SIZE(s) < UINT32_MAX-2) return NM_INT32;
  else if (YALE_MAX_SIZE(s) >= UINT64_MAX-2)
    fprintf(stderr, "WARNING: Matrix can contain no more than %llu non-diagonal non-zero entries, or results may be unpredictable\n", UINT64_MAX - SMMP_MIN(s->shape[0],s->shape[1]) - 2);
    // TODO: Turn this into an exception somewhere else. It's pretty unlikely, but who knows.
  return NM_INT64;
}


// Is the non-diagonal portion of the row empty?
static bool ndrow_is_empty(const YALE_STORAGE* s, y_size_t ija, const y_size_t ija_next, const void* ZERO) {
  //fprintf(stderr, "ndrow_is_empty: ija=%d, ija_next=%d\n", (size_t)(ija), (size_t)(ija_next));
  if (ija == ija_next) return true;
  while (ija < ija_next) { // do all the entries = zero?
    if (!ElemEqEq[s->dtype][0]((char*)s->a + nm_sizeof[s->dtype]*ija, ZERO, 1, nm_sizeof[s->dtype])) return false;
    ++ija;
  }
  return true;
}


// Are two non-diagonal rows the same? We already know
static bool ndrow_eqeq_ndrow(const YALE_STORAGE* l, const YALE_STORAGE* r, y_size_t l_ija, const y_size_t l_ija_next, y_size_t r_ija, const y_size_t r_ija_next, const void* ZERO) {
  y_size_t l_ja, r_ja, ja;
  bool l_no_more = false, r_no_more = false;

  YaleGetIJA(l_ja,   l,   l_ija);
  YaleGetIJA(r_ja,   r,   r_ija);
  ja = SMMP_MIN(l_ja, r_ja);

  //fprintf(stderr, "ndrow_eqeq_ndrow\n");
  while (!(l_no_more && r_no_more)) {
    //fprintf(stderr, "ndrow_eqeq_ndrow(loop): l_ija=%d, l_ija_next=%d, r_ija=%d, r_ija_next=%d\n", (size_t)(l_ija), (size_t)(l_ija_next), (size_t)(r_ija), (size_t)(r_ija_next));
    if (l_ja == r_ja) {
      if (!ElemEqEq[r->dtype][0]((char*)r->a + nm_sizeof[r->dtype]*r_ija, (char*)l->a + nm_sizeof[l->dtype]*l_ija, 1, nm_sizeof[l->dtype])) return false;

      ++l_ija;
      ++r_ija;

      if (l_ija < l_ija_next) YaleGetIJA(l_ja, l, l_ija);
      else                    l_no_more = true;

      if (r_ija < r_ija_next) YaleGetIJA(r_ja, r, r_ija);
      else                    r_no_more = true;

      ja = SMMP_MIN(l_ja, r_ja);
    } else if (l_no_more || ja < l_ja) {
      if (!ElemEqEq[r->dtype][0]((char*)r->a + nm_sizeof[r->dtype]*r_ija, ZERO, 1, nm_sizeof[r->dtype])) return false;

      ++r_ija;
      if (r_ija < r_ija_next) {
        YaleGetIJA(r_ja,  r,  r_ija); // get next column
        ja = SMMP_MIN(l_ja, r_ja);
      } else l_no_more = true;

    } else if (r_no_more || ja < r_ja) {
      if (!ElemEqEq[l->dtype][0]((char*)l->a + nm_sizeof[l->dtype]*l_ija, ZERO, 1, nm_sizeof[l->dtype])) return false;

      ++l_ija;
      if (l_ija < l_ija_next) {
        YaleGetIJA(l_ja,  l,   l_ija); // get next column
        ja = SMMP_MIN(l_ja, r_ja);
      } else l_no_more = true;

    } else {
      fprintf(stderr, "Unhandled in eqeq: l_ja=%d, r_ja=%d\n", l_ja, r_ja);
    }
  }

  return true; // every item matched
}


bool yale_storage_eqeq(const YALE_STORAGE* left, const YALE_STORAGE* right) {
  y_size_t l_ija, l_ija_next, r_ija, r_ija_next;
  y_size_t i = 0;

  // Need to know zero.
  void* ZERO = alloca(nm_sizeof[left->dtype]);
  if (left->dtype == NM_ROBJ) *(VALUE*)ZERO = INT2FIX(0);
  else                        memset(ZERO, 0, nm_sizeof[left->dtype]);

  // Easy comparison: Do the IJA and A vectors match exactly?
  //if (!memcmp(left->ija, right->ija, nm_sizeof[left->index_dtype]) && !memcmp(left->a, right->a, nm_sizeof[left->dtype])) return true;
  //fprintf(stderr, "yale_storage_eqeq\n");

  // Compare the diagonals first.
  if (!ElemEqEq[left->dtype][0](left->a, right->a, left->shape[0], nm_sizeof[left->dtype])) return false;

  while (i < left->shape[0]) {
    // Get start and end positions of row
    YaleGetIJA(l_ija,      left,  i);
    YaleGetIJA(l_ija_next, left,  i+1);
    YaleGetIJA(r_ija,      right, i);
    YaleGetIJA(r_ija_next, right, i+1);

    //fprintf(stderr, "yale_storage_eqeq: l_ija=%d, l_ija_next=%d, r_ija=%d, r_ija_next=%d, i=%d\n", (size_t)(l_ija), (size_t)(l_ija_next), (size_t)(r_ija), (size_t)(r_ija_next), (size_t)(i));

    // Check to see if one row is empty and the other isn't.
    if (ndrow_is_empty(left, l_ija, l_ija_next, ZERO)) {
      if (!ndrow_is_empty(right, r_ija, r_ija_next, ZERO)) return false;

    } else if (ndrow_is_empty(right, r_ija, r_ija_next, ZERO)) { // one is empty but the other isn't
      return false;
    } else if (!ndrow_eqeq_ndrow(left, right, l_ija, l_ija_next, r_ija, r_ija_next, ZERO)) { // Neither row is empty. Must compare the rows directly.
      return false;
    }

    ++i;
  }

  return true;
}


// Just like yale_vector_replace, but doesn't replace the contents of the cell, only the column index.
static char yale_vector_replace_j(YALE_STORAGE* s, y_size_t pos, y_size_t* j) {
  SetFuncs[s->index_dtype][Y_SIZE_T](1,
                                     pos*nm_sizeof[s->index_dtype] + (char*)(s->ija),
                                     nm_sizeof[s->index_dtype],
                                     j,
                                     sizeof(y_size_t));

  return 'r';
}


char yale_vector_replace(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n) {

  // Now insert the new values
  SetFuncs[s->index_dtype][Y_SIZE_T](n,
                                     pos*nm_sizeof[s->index_dtype] + (char*)(s->ija),
                                     nm_sizeof[s->index_dtype],
                                     j,
                                     sizeof(y_size_t));
  SetFuncs[s->dtype][s->dtype](n,
                               pos*nm_sizeof[s->dtype] + (char*)(s->a),
                               nm_sizeof[s->dtype],
                               val,
                               nm_sizeof[s->dtype]);

  return 'r';
}


char yale_vector_insert_resize(YALE_STORAGE* s, y_size_t current_size, y_size_t pos, y_size_t* j, void* val, y_size_t n, bool struct_only) {
  void *new_ija, *new_a;
  // Determine the new capacity for the IJA and A vectors.
  size_t new_capacity = s->capacity * YALE_GROWTH_CONSTANT;
  if (new_capacity > YALE_MAX_SIZE(s)) {
    new_capacity = YALE_MAX_SIZE(s);
    if (current_size + n > YALE_MAX_SIZE(s)) rb_raise(rb_eNoMemError, "insertion size exceeded maximum yale matrix size");
  }
  if (new_capacity < current_size + n) new_capacity = current_size + n;

  // Allocate the new vectors.
  new_ija     = ALLOC_N( char, nm_sizeof[s->index_dtype] * new_capacity );
  new_a       = ALLOC_N( char, nm_sizeof[s->dtype]       * new_capacity );

  // Check that allocation succeeded.
  if (!new_ija || !new_a) {
    free(new_a); free(new_ija);
    rb_raise(rb_eNoMemError, "yale sparse vectors are full and there is insufficient memory for growing them");
    return (char)false;
  }

  // Copy all values prior to the insertion site to the new IJA and new A
  SetFuncs[s->index_dtype][s->index_dtype](pos, new_ija, nm_sizeof[s->index_dtype], s->ija, nm_sizeof[s->index_dtype]);
  if (!struct_only)
    SetFuncs[s->dtype    ][s->dtype      ](pos, new_a,   nm_sizeof[s->dtype],       s->a,   nm_sizeof[s->dtype]      );

  // Copy all values subsequent to the insertion site to the new IJA and new A, leaving room (size n) for insertion.
  SetFuncs[s->index_dtype][s->index_dtype](current_size-pos+n-1, (char*)new_ija + nm_sizeof[s->index_dtype]*(pos+n), nm_sizeof[s->index_dtype], (char*)(s->ija) + nm_sizeof[s->index_dtype]*pos, nm_sizeof[s->index_dtype]);
  if (!struct_only)
    SetFuncs[s->dtype    ][s->dtype      ](current_size-pos+n-1, (char*)new_a   + nm_sizeof[s->dtype]*(pos+n),       nm_sizeof[s->dtype],       (char*)(s->a)   + nm_sizeof[s->dtype      ]*pos, nm_sizeof[s->dtype      ]);

  s->capacity = new_capacity;

  free(s->ija);
  free(s->a);

  s->ija = new_ija;
  s->a   = new_a;

  return 'i';
}


// Insert a value or contiguous values in the ija and a vectors (after ja and diag). Does not free anything; you are responsible!
//
// TODO: Improve this so it can handle non-contiguous element insertions efficiently.
// (For now, we can just sort the elements in the row in question.)
char yale_vector_insert(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n, bool struct_only) {
  y_size_t sz, i;

  if (pos < s->shape[0])
    rb_raise(rb_eArgError, "vector insert pos is before beginning of ja; this should not happen");

  YaleGetSize(sz, s);

  if (sz + n > s->capacity) yale_vector_insert_resize(s, sz, pos, j, val, n, struct_only);
  else {

    // No resize required:
    // easy (but somewhat slow), just copy elements to the tail, starting at the end, one element at a time.
    // TODO: This can be made slightly more efficient, but only after the tests are written.
    for (i = 0; i < sz - pos; ++i) {
      SetFuncs[s->index_dtype][s->index_dtype](1, (char*)(s->ija) + (sz+n-1-i)*nm_sizeof[s->index_dtype], 0, (char*)(s->ija) + (sz-1-i)*nm_sizeof[s->index_dtype], 0);

      if (!struct_only)
        SetFuncs[s->dtype    ][s->dtype      ](1, (char*)(s->a)   + (sz+n-1-i)*nm_sizeof[s->dtype      ], 0, (char*)(s->a)   + (sz-1-i)*nm_sizeof[s->dtype      ], 0);
    }
  }

  // Now insert the new values.
  if (struct_only) yale_vector_replace_j(s, pos, j);
  else             yale_vector_replace(s, pos, j, val, n);

  return 'i';
}


void delete_yale_storage(YALE_STORAGE* s) {
  if (s) {
    free(s->shape);
    free(s->ija);
    free(s->a);
    free(s);
  }
}


void mark_yale_storage(void* m) {
  size_t i;
  YALE_STORAGE* storage;

  if (m) {
    storage = (YALE_STORAGE*)(((NMATRIX*)m)->storage);
    fprintf(stderr, "mark_yale_storage\n");
    if (storage && storage->dtype == NM_ROBJ)
      for (i = 0; i < storage->capacity; ++i)
        rb_gc_mark(*((VALUE*)((char*)(storage->a) + i*nm_sizeof[NM_ROBJ])));
  }
}


// Allocate for a copy or copy-cast operation, and copy the IJA portion of the matrix (the structure).
static YALE_STORAGE* copy_alloc_yale_storage_struct(const YALE_STORAGE* rhs, const int8_t new_dtype, const y_size_t new_capacity, const y_size_t new_size) {
  YALE_STORAGE* lhs = ALLOC( YALE_STORAGE );
  lhs->rank         = rhs->rank;
  lhs->shape        = ALLOC_N( size_t, lhs->rank );
  memcpy(lhs->shape, rhs->shape, lhs->rank * sizeof(size_t));
  lhs->index_dtype  = rhs->index_dtype;
  lhs->capacity     = new_capacity;
  lhs->dtype        = new_dtype;
  lhs->ndnz         = rhs->ndnz;

  lhs->ija          = ALLOC_N( char, nm_sizeof[lhs->index_dtype] * lhs->capacity );
  lhs->a            = ALLOC_N( char, nm_sizeof[lhs->dtype] * lhs->capacity );

  // Now copy the contents -- but only within the boundaries set by the size. Leave
  // the rest uninitialized.
  YaleGetSize(new_size, rhs);
  memcpy(lhs->ija, rhs->ija, new_size * nm_sizeof[lhs->index_dtype]); // indices

  return lhs;
}


// copy constructor
YALE_STORAGE* copy_yale_storage(YALE_STORAGE* rhs) {
  y_size_t size;
  YALE_STORAGE* lhs;

  YaleGetSize(size, rhs);
  lhs = copy_alloc_yale_storage_struct(rhs, rhs->dtype, rhs->capacity, size);

  // Now copy the contents -- but only within the boundaries set by the size. Leave
  // the rest uninitialized.
  memcpy(lhs->a, rhs->a, size * nm_sizeof[lhs->dtype]);

  return lhs;
}



// copy constructor
YALE_STORAGE* cast_copy_yale_storage(YALE_STORAGE* rhs, int8_t new_dtype) {
  y_size_t size;
  YALE_STORAGE* lhs;

  YaleGetSize(size, rhs);
  lhs = copy_alloc_yale_storage_struct(rhs, new_dtype, rhs->capacity, size);

  // Copy contents (not structure)
  if (lhs->dtype == rhs->dtype) memcpy(lhs->a, rhs->a, size * nm_sizeof[lhs->dtype]);
  else                          SetFuncs[new_dtype][rhs->dtype](size, lhs->a, nm_sizeof[lhs->dtype], rhs->a, nm_sizeof[rhs->dtype]);

  return lhs;
}


// Clear out the D portion of the A vector (clearing the diagonal and setting the zero value).
//
// Note: This sets a literal 0 value. If your dtype is NM_ROBJ (a Ruby object), it'll actually
// be INT2FIX(0) instead of a string of NULLs.
static void clear_diagonal_and_zero(YALE_STORAGE* s) {
  y_size_t i;
  // Clear out the diagonal + one extra entry
  if (s->dtype == NM_ROBJ) {
    for (i = 0; i < YALE_IA_SIZE(s)+1; ++i) // insert Ruby zeros
      *(VALUE*)( (char*)(s->a) + i*nm_sizeof[s->dtype] ) = INT2FIX(0);
  } else { // just insert regular zeros
    memset(s->a, 0, nm_sizeof[s->dtype] * (YALE_IA_SIZE(s)+1));
  }
}



YALE_STORAGE* scast_copy_yale_list(const LIST_STORAGE* rhs, int8_t l_dtype) {
  YALE_STORAGE* lhs;
  size_t* shape;
  NODE *i_curr, *j_curr;
  y_size_t ija;
  size_t ndnz = count_list_storage_nd_elements(rhs);

  if (rhs->rank != 2)
    rb_raise(nm_eStorageTypeError, "can only convert matrices of rank 2 to yale");

  if ((rhs->dtype == NM_ROBJ && *(VALUE*)(rhs->default_val) == INT2FIX(0)) || strncmp(rhs->default_val, "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", nm_sizeof[rhs->dtype]))
    rb_raise(nm_eStorageTypeError, "list matrix must have default value of 0 to convert to yale");

  // Copy shape for yale construction
  shape = ALLOC_N(size_t, 2);
  shape[0] = rhs->shape[0];
  shape[1] = rhs->shape[1];

  lhs = create_yale_storage(l_dtype, shape, 2, shape[0] + ndnz + 1);
  clear_diagonal_and_zero(lhs); // clear the diagonal and the zero location.

  ija = lhs->shape[0]+1;

  for (i_curr = rhs->rows->first; i_curr; i_curr = i_curr->next) {

    // indicate the beginning of a row in the IJA array
    YaleSetIJA(i_curr->key, lhs, ija);

    for (j_curr = ((LIST*)(i_curr->val))->first; j_curr; j_curr = j_curr->next) {
      if (i_curr->key == j_curr->key) {
        // set diagonal
        SetFuncs[l_dtype][rhs->dtype](1, (char*)(lhs->a) + (i_curr->key)*nm_sizeof[l_dtype], 0, j_curr->val, 0);
      } else {
        // set column value
        YaleSetIJA(ija, lhs, j_curr->key); // write column index

        // set cell value
        SetFuncs[l_dtype][rhs->dtype](1, (char*)(lhs->a) + ija*nm_sizeof[l_dtype], 0, j_curr->val, 0);

        ++ija;
      }
    }

    if (!i_curr->next) YaleSetIJA(i_curr->key, lhs, ija); // indicate the end of the last row
  }

  lhs->ndnz = ndnz;
  return lhs;
}


YALE_STORAGE* scast_copy_yale_dense(const DENSE_STORAGE* rhs, int8_t l_dtype) {
  YALE_STORAGE* lhs;
  size_t i, j;
  y_size_t pos = 0, ndnz = 0, ija;
  size_t* shape;

  // Figure out values to write for zero in yale and compare to zero in dense
  void *R_ZERO = ALLOCA_N(char, nm_sizeof[rhs->dtype]);
  if (rhs->dtype == NM_ROBJ) *(VALUE*)R_ZERO = INT2FIX(0);
  else                       memset(R_ZERO, 0, nm_sizeof[rhs->dtype]);

  if (rhs->rank != 2)
    rb_raise(nm_eStorageTypeError, "can only convert matrices of rank 2 to yale");

  // First, count the non-diagonal nonzeros
  for (i = 0; i < rhs->shape[0]; ++i) {
    for (j = 0; j < rhs->shape[1]; ++j) {
      if (i != j && !ElemEqEq[rhs->dtype][0]((char*)(rhs->elements) + pos*nm_sizeof[rhs->dtype], R_ZERO, 1, nm_sizeof[rhs->dtype])) ++ndnz;
      ++pos; // move forward 1 position in dense matrix elements array
    }
  }

  // Copy shape for yale construction
  shape = ALLOC_N(size_t, 2);
  shape[0] = rhs->shape[0];
  shape[1] = rhs->shape[1];

  // Create with minimum possible capacity -- just enough to hold all of the entries
  lhs = create_yale_storage(l_dtype, shape, 2, shape[0] + ndnz + 1);

  // Set the zero position in the yale matrix
  cast_copy_value_single((char*)(lhs->a) + shape[0]*nm_sizeof[l_dtype], R_ZERO, l_dtype, rhs->dtype);

  // Start just after the zero position.
  ija = lhs->shape[0]+1;
  pos = 0;

  // Copy contents
  for (i = 0; i < rhs->shape[0]; ++i) {
    // indicate the beginning of a row in the IJA array
    YaleSetIJA(i, lhs, ija);

    for (j = 0; j < rhs->shape[1]; ++j) {

      if (i == j) { // copy to diagonal
        cast_copy_value_single((char*)(lhs->a) + i*nm_sizeof[l_dtype], (char*)(rhs->elements) + pos*nm_sizeof[rhs->dtype], l_dtype, rhs->dtype);

      } else if (!ElemEqEq[rhs->dtype][0]((char*)(rhs->elements) + pos*nm_sizeof[rhs->dtype], R_ZERO, 1, nm_sizeof[rhs->dtype])) {      // copy nonzero to LU
        YaleSetIJA(ija, lhs, j); // write column index

        cast_copy_value_single((char*)(lhs->a) + ija*nm_sizeof[l_dtype], (char*)(rhs->elements) + pos*nm_sizeof[rhs->dtype], l_dtype, rhs->dtype);

        ++ija;
      }
      ++pos;
    }
  }
  YaleSetIJA(i, lhs, ija); // indicate the end of the last row

  lhs->ndnz = ndnz;

  return lhs;
}


// If we add n items to row i, we need to increment ija[i+1] and onward
// TODO: Add function pointer array for AddFuncs, the same way we have SetFuncs
static void yale_storage_increment_ia_after(YALE_STORAGE* s, y_size_t ija_size, y_size_t i, y_size_t n) {
  y_size_t val;

  ++i;
  for (; i <= ija_size; ++i) {
    YaleGetIJA(val, s, i);
    val += n;
    YaleSetIJA(i, s, val);
  }
}


// Binary search for returning insertion points
static y_size_t yale_storage_insert_search(YALE_STORAGE* s, y_size_t left, y_size_t right, y_size_t key, bool* found) {
  y_size_t mid = (left + right)/2, mid_j;

  if (left > right) {
    *found = false;
    return left;
  }

  YaleGetIJA(mid_j, s, mid);
  if (mid_j == key) {
    *found = true;
    return mid;
  }
  else if (mid_j > key) return yale_storage_insert_search(s, left, mid-1, key, found);
  else return yale_storage_insert_search(s, mid+1, right, key, found);
}



YALE_STORAGE* create_merged_yale_storage(const YALE_STORAGE* template, const YALE_STORAGE* other) {
  y_size_t ija, ija_next, o_ija, o_ija_next;
  y_size_t ja, o_ja, size, pos;
  YALE_STORAGE* s;
  bool found;
  char ins_type;
  size_t i;

  YaleGetSize(size, template);
  s = copy_alloc_yale_storage_struct(template, template->dtype, SMMP_MAX(template->capacity, other->capacity), size);

  // set the element between D and LU (the boundary in A), which should be 0.
  memcpy((char*)(s->a) + nm_sizeof[s->dtype] * s->shape[0], (char*)(template->a) + nm_sizeof[template->dtype] * template->shape[0], nm_sizeof[s->dtype]);

  if (other && other != template) { // some operations are unary and don't need this; others are x+x and don't need this
    for (i = 0; i < s->shape[0]; ++i) {
      YaleGetIJA(ija,        s, i);
      YaleGetIJA(ija_next,   s, i+1);

      YaleGetIJA(o_ija,      other, i);
      YaleGetIJA(o_ija_next, other, i+1);

      while (o_ija < o_ija_next) {
        YaleGetIJA(o_ja,     other, o_ija);
        YaleGetIJA(ja,       s,     ija);

        if (ija == ija_next) { // destination row is empty
          ins_type = yale_vector_insert(s, ija, &ja, NULL, 1, true);
          yale_storage_increment_ia_after(s, YALE_IA_SIZE(s), i, 1);
          ++(s->ndnz);
          ++ija;
          if (ins_type == 'i') ++ija_next;
        } else { // merge positions into destination row
          pos = yale_storage_insert_search(s, ija, ija_next-1, ja, &found);
          if (!found) {
            yale_vector_insert(s, pos, &ja, NULL, 1, true);
            yale_storage_increment_ia_after(s, YALE_IA_SIZE(s), i, 1);
            ++(s->ndnz);
            if (ins_type == 'i') ++ija_next;
          }
          ija = pos+1; // can now set a left boundary for the next search
        }

        ++o_ija;
      }
    }
  }

  return s;
}


// Allocates and initializes the basic struct (but not the IJA or A vectors).
static YALE_STORAGE* alloc_yale_storage(int8_t dtype, size_t* shape, size_t rank) {
  YALE_STORAGE* s;

  s = ALLOC( YALE_STORAGE );

  s->ndnz        = 0;
  s->dtype       = dtype;
  s->shape       = shape;
  s->rank        = rank;
  s->index_dtype = yale_index_dtype(s);

  return s;
}


YALE_STORAGE* create_yale_storage_from_old_yale(int8_t dtype, size_t* shape, char* ia, char* ja, char* a, int8_t from_dtype, int8_t from_index_dtype) {
  YALE_STORAGE* s;
  y_size_t i = 0, p, p_next, j, ndnz = 0, pp;

  // Read through ia and ja and figure out the ndnz count.

  // Walk down rows
  for (i = 0; i < shape[0]; ++i) {
    SetFuncs[Y_SIZE_T][from_index_dtype](1, &p, 0, ia + nm_sizeof[from_index_dtype]*i, 0);          // p = ia[i]
    SetFuncs[Y_SIZE_T][from_index_dtype](1, &p_next, 0, ia + nm_sizeof[from_index_dtype]*(i+1), 0); // p_next = ia[i+1]

    // Now walk through columns
    for (; p < p_next; ++p) {
      SetFuncs[Y_SIZE_T][from_index_dtype](1, &j, 0, ja + nm_sizeof[from_index_dtype]*p, 0);        // j = ja[p]

      if (i != j) ++ndnz; // entry is non-diagonal and probably nonzero
    }
  }

  // Having walked through the matrix, we now go about allocating the space for it.
  s = alloc_yale_storage(dtype, shape, 2);
  s->capacity = shape[0] + ndnz + 1;
  s->ndnz     = ndnz;

  // Setup IJA and A arrays
  s->ija = ALLOC_N( char, nm_sizeof[s->index_dtype] * s->capacity );
  s->a   = ALLOC_N( char, nm_sizeof[s->dtype] * s->capacity );

  // Figure out where to start writing JA in IJA:
  pp     = s->shape[0]+1;

  // Find beginning of first row
  SetFuncs[Y_SIZE_T][from_index_dtype](1, &p, 0, ia, 0);          // p = ia[i]

  // Now fill the arrays
  for (i = 0; i < s->shape[0]; ++i) {

    // Find end of row (of input)
    SetFuncs[Y_SIZE_T][from_index_dtype](1, &p_next, 0, ia + nm_sizeof[from_index_dtype]*(i+1), 0); // p_next = ia[i+1]

    // Set the beginning of the row (of output)
    SetFuncs[s->index_dtype][Y_SIZE_T](1, (char*)(s->ija) + nm_sizeof[s->index_dtype]*i, 0, &pp, 0);

    // Now walk through columns
    for (; p < p_next; ++p, ++pp) {
      SetFuncs[Y_SIZE_T][from_index_dtype](1, &j, 0, ja + nm_sizeof[from_index_dtype]*p, 0);        // j = ja[p]

      if (i == j) {
        SetFuncs[s->dtype][from_dtype](1, (char*)(s->a) + nm_sizeof[s->dtype]*i, 0, a + nm_sizeof[from_dtype]*p, 0);
        --pp;
      } else {
        SetFuncs[s->index_dtype][from_index_dtype](1, (char*)(s->ija) + nm_sizeof[s->index_dtype]*pp, 0, ja + nm_sizeof[from_index_dtype]*p, 0);
        SetFuncs[s->dtype][from_dtype](1,             (char*)(s->a)   + nm_sizeof[s->dtype]*pp,       0, a  + nm_sizeof[from_dtype]*p,       0);
      }
    }
  }

  // Set the end of the last row
  SetFuncs[s->index_dtype][Y_SIZE_T](1, (char*)(s->ija) + nm_sizeof[s->index_dtype]*i, 0, &pp, 0);

  // Set the zero position for our output matrix
  if (dtype == NM_ROBJ)     *(VALUE*)((char*)(s->a) + nm_sizeof[s->dtype]*i) = INT2FIX(0);
  else                      memset((char*)(s->a) + nm_sizeof[s->dtype]*i, 0, nm_sizeof[s->dtype]);

  return s;
}


YALE_STORAGE* create_yale_storage(int8_t dtype, size_t* shape, size_t rank, size_t init_capacity) {
  YALE_STORAGE* s;
  size_t max_capacity;

  if (rank != 2) rb_raise(rb_eNotImpError, "Can only support 2D matrices");

  s = alloc_yale_storage(dtype, shape, rank);
  max_capacity = count_storage_max_elements((STORAGE*)s) - s->shape[0] + 1;

  // Set matrix capacity (and ensure its validity)
  if (init_capacity < YALE_MINIMUM(s))   s->capacity = YALE_MINIMUM(s);
  else if (init_capacity > max_capacity) s->capacity = max_capacity; // Don't allow storage to be created larger than necessary
  else                                   s->capacity = init_capacity;

  s->ija = ALLOC_N( char, nm_sizeof[s->index_dtype] * s->capacity );
  s->a   = ALLOC_N( char, nm_sizeof[s->dtype] * s->capacity );

  return s;
}


// Empty the matrix
void init_yale_storage(YALE_STORAGE* s) {
  y_size_t IA_INIT = YALE_IA_SIZE(s)+1, i;

  // clear out IJA vector
  for (i = 0; i < YALE_IA_SIZE(s)+1; ++i)
    SetFuncs[s->index_dtype][Y_SIZE_T](1, (char*)(s->ija) + i*nm_sizeof[s->index_dtype], 0, &IA_INIT, 0); // set initial values for IJA

  clear_diagonal_and_zero(s);
}


char yale_storage_set_diagonal(YALE_STORAGE* s, y_size_t i, void* v) {
  memcpy(YALE_DIAG(s, nm_sizeof[s->dtype], i), v, nm_sizeof[s->dtype]); return 'r';
}


// Binary search for returning stored values
int yale_storage_binary_search(YALE_STORAGE* s, y_size_t left, y_size_t right, y_size_t key) {
  y_size_t mid = (left + right)/2, mid_j;

  if (left > right) return -1;

  YaleGetIJA(mid_j, s, mid);
  if (mid_j == key) return mid;
  else if (mid_j > key) return yale_storage_binary_search(s, left, mid-1, key);
  else return yale_storage_binary_search(s, mid+1, right, key);
}




char yale_storage_set(YALE_STORAGE* s, SLICE* slice, void* v) {
  size_t* coords = slice->coords;

  y_size_t i_next = coords[0] + 1;
  y_size_t ija, ija_next, ija_size;
  y_size_t pos;
  bool found = false;
  char ins_type;

  if (coords[0] == coords[1]) return yale_storage_set_diagonal(s, coords[0], v);

  // Get IJA positions of the beginning and end of the row
  YaleGetIJA(ija,      s, coords[0]);
  YaleGetIJA(ija_next, s, i_next);

  if (ija == ija_next) { // empty row
    ins_type = yale_vector_insert(s, ija, &(coords[1]), v, 1, false);
    yale_storage_increment_ia_after(s, YALE_IA_SIZE(s), coords[0], 1);
    s->ndnz++;
    return ins_type;
  }

  // non-empty row. search for coords[1] in the IJA array, between ija and ija_next
  // (including ija, not including ija_next)
  YaleGetSize(ija_size, s);
  //--ija_next;

  // Do a binary search for the column
  pos = yale_storage_insert_search(s, ija, ija_next-1, coords[1], &found);

  if (found) return yale_vector_replace(s, pos, &(coords[1]), v, 1);

  ins_type = yale_vector_insert(s, pos, &(coords[1]), v, 1, false);
  yale_storage_increment_ia_after(s, YALE_IA_SIZE(s), coords[0], 1);
  s->ndnz++;
  return ins_type;

}


void* yale_storage_ref(YALE_STORAGE* s, SLICE* slice) {
  size_t* coords = slice->coords;

  y_size_t l, r, i_plus_one = coords[0] + 1, test_j;
  int pos;

  if (coords[0] == coords[1]) return YALE_DIAG(s,nm_sizeof[s->dtype],coords[0]);

  YaleGetIJA(l, s, coords[0]);
  YaleGetIJA(r, s, i_plus_one);

  if (l == r) return YALE_A(s, nm_sizeof[s->dtype], s->shape[0]); // return zero pointer

  pos = yale_storage_binary_search(s, l, r-1, coords[1]); // binary search for the column's location
  if (pos != -1) {
    YaleGetIJA(test_j, s, pos);
    if (test_j == coords[1])
      return YALE_A(s, nm_sizeof[s->dtype], pos); // Found exact value.
  }

  // return a pointer that happens to be zero
  return YALE_A(s, nm_sizeof[s->dtype], s->shape[0]);
}

#endif
