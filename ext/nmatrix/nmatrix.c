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
// == nmatrix.c
//
#ifndef NMATRIX_C
# define NMATRIX_C

#include <ruby.h>

#include "nmatrix.h"
//#include "types.h"

VALUE cNMatrix, cNVector;
VALUE nm_eDataTypeError, nm_eStorageTypeError;

ID nm_id_real, nm_id_imag;
ID nm_id_numer, nm_id_denom;
ID nm_id_transpose, nm_id_no_transpose, nm_id_complex_conjugate; // cblas
ID nm_id_list, nm_id_dense;
ID nm_id_mult, nm_id_multeq;
ID nm_id_add;

#include "dtypes.c"

#ifdef BENCHMARK
double get_time() {
  struct timeval t;
  struct timezone tzp;
  gettimeofday(&t, &tzp);
  return t.tv_sec + t.tv_usec*1e-6;
}
#endif


const char *nm_stypestring[] = {
  "dense",
  "list",
  "yale",
  "stypes"
};


nm_delete_t DeleteFuncs = {
  delete_dense_storage,
  delete_list_storage,
  delete_yale_storage
};

nm_delete_t DeleteFuncsRef = {
  delete_dense_storage_ref,
  delete_list_storage,
  delete_yale_storage
};



nm_mark_t MarkFuncs = {
  mark_dense_storage,
  mark_list_storage,
  mark_yale_storage
};


nm_gemv_t GemvFuncs = {
  NULL,
  cblas_bgemv_,
  cblas_i8gemv_,
  cblas_i16gemv_,
  cblas_i32gemv_,
  cblas_i64gemv_,
  cblas_sgemv_,
  cblas_dgemv_,
  cblas_cgemv_,
  cblas_zgemv_,
  cblas_r32gemv_,
  cblas_r64gemv_,
  cblas_r128gemv_,
  NULL
};


nm_gemm_t GemmFuncs = { // by NM_TYPES
  NULL,
  cblas_bgemm_,
  cblas_i8gemm_,
  cblas_i16gemm_,
  cblas_i32gemm_,
  cblas_i64gemm_,
  cblas_sgemm_,
  cblas_dgemm_,
  cblas_cgemm_,
  cblas_zgemm_,
  cblas_r32gemm_,
  cblas_r64gemm_,
  cblas_r128gemm_,
  cblas_vgemm_
};

static void TransposeTypeErr(y_size_t n, y_size_t m, YALE_PARAM A, YALE_PARAM B, bool move) {
  rb_raise(nm_eDataTypeError, "illegal operation with this matrix type");
}


// First dimension is dtype, second dimension is index dtype (so lots of nulls)
nm_smmp_transpose_t SparseTransposeFuncs = {
  {TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_NONE
  {TransposeTypeErr, TransposeTypeErr, i8_b_transp_, i16_b_transp_, i32_b_transp_, i64_b_transp_}, // NM_BYTE
  {TransposeTypeErr, TransposeTypeErr, i8_i8_transp_, i16_i8_transp_, i32_i8_transp_, i64_i8_transp_}, // NM_INT8
  {TransposeTypeErr, TransposeTypeErr, i8_i16_transp_, i16_i16_transp_, i32_i16_transp_, i64_i16_transp_}, // NM_INT16
  {TransposeTypeErr, TransposeTypeErr, i8_i32_transp_, i16_i32_transp_, i32_i32_transp_, i64_i32_transp_}, // NM_INT32
  {TransposeTypeErr, TransposeTypeErr, i8_i64_transp_, i16_i64_transp_, i32_i64_transp_, i64_i64_transp_}, // NM_INT64
  {TransposeTypeErr, TransposeTypeErr, i8_f32_transp_, i16_f32_transp_, i32_f32_transp_, i64_f32_transp_}, // NM_FLOAT32
  {TransposeTypeErr, TransposeTypeErr, i8_f64_transp_, i16_f64_transp_, i32_f64_transp_, i64_f64_transp_}, // NM_FLOAT64
  {TransposeTypeErr, TransposeTypeErr, i8_c64_transp_, i16_c64_transp_, i32_c64_transp_, i64_c64_transp_}, // NM_COMPLEX64
  {TransposeTypeErr, TransposeTypeErr, i8_c128_transp_, i16_c128_transp_, i32_c128_transp_, i64_c128_transp_}, // NM_COMPLEX128
  {TransposeTypeErr, TransposeTypeErr, i8_r32_transp_, i16_r32_transp_, i32_r32_transp_, i64_r32_transp_}, // NM_RATIONAL32
  {TransposeTypeErr, TransposeTypeErr, i8_r64_transp_, i16_r64_transp_, i32_r64_transp_, i64_r64_transp_}, // NM_RATIONAL64
  {TransposeTypeErr, TransposeTypeErr, i8_r128_transp_, i16_r128_transp_, i32_r128_transp_, i64_r128_transp_}, // NM_RATIONAL128
  {TransposeTypeErr, TransposeTypeErr, i8_v_transp_, i16_v_transp_, i32_v_transp_, i64_v_transp_}  // NM_ROBJ
};

/*
// Currently commented out because dense_transpose_generic is about the same speed. Let's resurrect this when we write
// an in-place transpose (e.g., transpose!).

static void DenseTransTypeErr(int M, int N, void* A, int lda, void* B, int ldb, bool move) {
  rb_raise(nm_eDataTypeError, "illegal operation with this matrix type");
}

nm_dense_transpose_t DenseTransposeFuncs = {
  DenseTransTypeErr,
  btransp,
  i8transp,
  i16transp,
  i32transp,
  i64transp,
  f32transp,
  f64transp,
  c64transp,
  c128transp,
  r32transp,
  r64transp,
  r128transp,
  vtransp
}; */


static void SmmpTypeErr(y_size_t n, y_size_t m, YALE_PARAM A, YALE_PARAM B, YALE_PARAM C) {
  rb_raise(nm_eDataTypeError, "illegal operation with this matrix type");
}

// First dimension is dtype, second dimension is index dtype (so lots of nulls)
nm_smmp_t SmmpFuncs = {
  {SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_NONE
  {SmmpTypeErr, SmmpTypeErr, i8_b_smmp, i16_b_smmp, i32_b_smmp, i64_b_smmp}, // NM_BYTE
  {SmmpTypeErr, SmmpTypeErr, i8_i8_smmp, i16_i8_smmp, i32_i8_smmp, i64_i8_smmp}, // NM_INT8
  {SmmpTypeErr, SmmpTypeErr, i8_i16_smmp, i16_i16_smmp, i32_i16_smmp, i64_i16_smmp}, // NM_INT16
  {SmmpTypeErr, SmmpTypeErr, i8_i32_smmp, i16_i32_smmp, i32_i32_smmp, i64_i32_smmp}, // NM_INT32
  {SmmpTypeErr, SmmpTypeErr, i8_i64_smmp, i16_i64_smmp, i32_i64_smmp, i64_i64_smmp}, // NM_INT64
  {SmmpTypeErr, SmmpTypeErr, i8_f32_smmp, i16_f32_smmp, i32_f32_smmp, i64_f32_smmp}, // NM_FLOAT32
  {SmmpTypeErr, SmmpTypeErr, i8_f64_smmp, i16_f64_smmp, i32_f64_smmp, i64_f64_smmp}, // NM_FLOAT64
  {SmmpTypeErr, SmmpTypeErr, i8_c64_smmp, i16_c64_smmp, i32_c64_smmp, i64_c64_smmp}, // NM_COMPLEX64
  {SmmpTypeErr, SmmpTypeErr, i8_c128_smmp, i16_c128_smmp, i32_c128_smmp, i64_c128_smmp}, // NM_COMPLEX128
  {SmmpTypeErr, SmmpTypeErr, i8_r32_smmp, i16_r32_smmp, i32_r32_smmp, i64_r32_smmp}, // NM_RATIONAL32
  {SmmpTypeErr, SmmpTypeErr, i8_r64_smmp, i16_r64_smmp, i32_r64_smmp, i64_r64_smmp}, // NM_RATIONAL64
  {SmmpTypeErr, SmmpTypeErr, i8_r128_smmp, i16_r128_smmp, i32_r128_smmp, i64_r128_smmp}, // NM_RATIONAL128
  {SmmpTypeErr, SmmpTypeErr, i8_v_smmp, i16_v_smmp, i32_v_smmp, i64_v_smmp}  // NM_ROBJ
};


static inline DENSE_PARAM cblas_params_for_multiply(const DENSE_STORAGE* left, const DENSE_STORAGE* right, const DENSE_STORAGE* result, bool vector) {
  DENSE_PARAM p;

  p.A = left->elements;
  p.B = right->elements;    // for vector, this is actually x
  p.C = result->elements;   // vector Y

  p.M = left->shape[0];
  p.lda = left->shape[1];

  if (vector) {
    p.N = left->shape[1];

    p.ldb = 1;  // incX
    p.ldc = 1;  // incY
  } else {
    p.N = right->shape[1];
    p.K = left->shape[1];

    p.ldb = right->shape[1];
    p.ldc = result->shape[1];
  }

  switch(left->dtype) {
  case NM_FLOAT32:
  case NM_FLOAT64:
    p.alpha.d[0] = 1.0;
    p.beta.d[0] = 0.0;
    break;

  case NM_COMPLEX64:
    p.alpha.c[0].r = 1.0;
    p.alpha.c[0].i = 0.0;
    p.beta.c[0].r = 0.0;
    p.beta.c[0].i = 0.0;
    break;

  case NM_COMPLEX128:
    p.alpha.z.r = 1.0;
    p.alpha.z.i = 0.0;
    p.beta.z.r = 0.0;
    p.beta.z.i = 0.0;
    break;

  case NM_BYTE:
    p.alpha.b[0] = 1;
    p.beta.b[0] = 0;
    break;

  case NM_INT8:
  case NM_INT16:
  case NM_INT32:
  case NM_INT64:
    p.alpha.i[0] = 1;
    p.beta.i[0]  = 0;
    break;

  case NM_RATIONAL32:
    p.alpha.r[0].n = 1;
    p.alpha.r[0].d = 1;
    p.beta.r[0].n = 0;
    p.beta.r[0].d = 1;
    break;

  case NM_RATIONAL64:
    p.alpha.ra[0].n = 1;
    p.alpha.ra[0].d = 1;
    p.beta.ra[0].n = 0;
    p.beta.ra[0].d = 1;
    break;

  case NM_RATIONAL128:
    p.alpha.rat.n = 1;
    p.alpha.rat.d = 1;
    p.beta.rat.n = 0;
    p.beta.rat.d = 1;
    break;

  case NM_ROBJ:
    p.alpha.v[0] = INT2FIX(1);
    p.beta.v[0]  = RUBY_ZERO;
    break;

  default:
    rb_raise(nm_eDataTypeError, "unexpected dtype");

  }

  return p;
}


static NMATRIX* multiply_matrix_dense_casted(STORAGE_PAIR casted_storage, size_t* resulting_shape, bool vector) {
  DENSE_STORAGE *left  = (DENSE_STORAGE*)(casted_storage.left),
                *right = (DENSE_STORAGE*)(casted_storage.right),
                *result;

  // We can safely get dtype from the casted matrices; post-condition of binary_storage_cast_alloc is that dtype is the
  // same for left and right.
  int8_t dtype = left->dtype;

  // Create result storage.
  result = create_dense_storage(dtype, resulting_shape, 2, NULL, 0);

  // Do the multiplication
  if (vector) GemvFuncs[dtype](CblasRowMajor, CblasNoTrans, cblas_params_for_multiply(left, right, result, true));
  else        GemmFuncs[dtype](CblasRowMajor, CblasNoTrans, CblasNoTrans, cblas_params_for_multiply(left, right, result, false));

  return nm_create(S_DENSE, result);
}


static NMATRIX* multiply_matrix_yale_casted(STORAGE_PAIR casted_storage, size_t* resulting_shape, bool vector) {
  YALE_STORAGE *left  = (YALE_STORAGE*)(casted_storage.left),
               *right = (YALE_STORAGE*)(casted_storage.right),
               *result;
  YALE_PARAM A, B, C;

  // We can safely get dtype from the casted matrices; post-condition of binary_storage_cast_alloc is that dtype is the
  // same for left and right.
  int8_t dtype = left->dtype;

  // Create result storage.
  result = create_yale_storage(dtype, resulting_shape, 2, left->capacity + right->capacity);
  init_yale_storage(result);

  // Set multiplication parameters
  A.ia = A.ja = left->ija;
  A.a  = left->a;
  B.ia = B.ja = right->ija;
  B.a  = right->a;
  C.ia = C.ja = result->ija;
  C.a  = result->a;

  A.diag = B.diag = C.diag = true;

  // Do the multiplication
  SmmpFuncs[dtype][left->index_dtype](result->shape[0], result->shape[1], A, B, C);

  return nm_create(S_YALE, result);
}


static NMATRIX* multiply_matrix_list_casted(STORAGE_PAIR casted_storage, size_t* resulting_shape) {
  rb_raise(rb_eNotImpError, "multiplication not implemented for list-of-list matrices");
  free(resulting_shape);
  return NULL;
}


nm_matrix_multiply_op_t CastedMultiplyFuncs = {
  multiply_matrix_dense_casted,
  multiply_matrix_list_casted,
  multiply_matrix_yale_casted
};


nm_d_elementwise_binary_op_t DenseElementwiseFuncs = { // only for dense!
  NULL,
  nm_d_b_elementwise,
  nm_d_i8_elementwise,
  nm_d_i16_elementwise,
  nm_d_i32_elementwise,
  nm_d_i64_elementwise,
  nm_d_f32_elementwise,
  nm_d_f64_elementwise,
  nm_d_c64_elementwise,
  nm_d_c128_elementwise,
  nm_d_r32_elementwise,
  nm_d_r64_elementwise,
  nm_d_r128_elementwise,
  nm_d_v_elementwise,
  NULL
};

nm_det_t DenseDetExact = {
  NULL,
  bdet_exact,
  i8det_exact,
  i16det_exact,
  i32det_exact,
  i64det_exact,
  f32det_exact,
  f64det_exact,
  c64det_exact,
  c128det_exact,
  r32det_exact,
  r64det_exact,
  r128det_exact,
  vdet_exact,
  NULL
};

static void EwTypeErr(y_size_t n, enum NMatrix_Ops op, void* ija, void* ijb, void* ijc, void* a, void* b, void* c) {
  rb_raise(nm_eDataTypeError, "illegal operation with this matrix type");
}

// First dimension is dtype, second dimension is index dtype (so lots of nulls)
nm_y_elementwise_binary_op_t YaleElementwiseFuncs = { // only for yale!
  {EwTypeErr, EwTypeErr, EwTypeErr, EwTypeErr,  EwTypeErr,  EwTypeErr},
  {EwTypeErr, EwTypeErr, EwTypeErr, EwTypeErr,  EwTypeErr,  EwTypeErr},
  {EwTypeErr, EwTypeErr, i8_i8_ew,  i16_i8_ew,  i32_i8_ew,  i64_i8_ew},
  {EwTypeErr, EwTypeErr, i8_i16_ew, i16_i16_ew, i32_i16_ew, i64_i16_ew},
  {EwTypeErr, EwTypeErr, i8_i32_ew, i16_i32_ew, i32_i32_ew, i64_i32_ew},
  {EwTypeErr, EwTypeErr, i8_i64_ew, i16_i64_ew, i32_i64_ew, i64_i64_ew},
  {EwTypeErr, EwTypeErr, i8_f32_ew, i16_f32_ew, i32_f32_ew, i64_f32_ew},
  {EwTypeErr, EwTypeErr, i8_f64_ew, i16_f64_ew, i32_f64_ew, i64_f64_ew},
  {EwTypeErr, EwTypeErr, i8_c64_ew, i16_c64_ew, i32_c64_ew, i64_c64_ew},
  {EwTypeErr, EwTypeErr, i8_c128_ew,i16_c128_ew,i32_c128_ew,i64_c128_ew},
  {EwTypeErr, EwTypeErr, i8_r32_ew, i16_r32_ew, i32_r32_ew, i64_r32_ew},
  {EwTypeErr, EwTypeErr, i8_r64_ew, i16_r64_ew, i32_r64_ew, i64_r64_ew},
  {EwTypeErr, EwTypeErr, i8_r128_ew,i16_r128_ew,i32_r128_ew,i64_r128_ew},
  {EwTypeErr, EwTypeErr, i8_v_ew,   i16_v_ew,   i32_v_ew,   i64_v_ew}
};


static NMATRIX* elementwise_dense_casted(STORAGE_PAIR casted_storage, char op) {
  DENSE_STORAGE *left  = (DENSE_STORAGE*)(casted_storage.left),
                *right = (DENSE_STORAGE*)(casted_storage.right),
                *result;

  // We can safely get dtype from the casted matrices; post-condition of binary_storage_cast_alloc is that dtype is the
  // same for left and right.
  size_t i;
  int8_t dtype = left->dtype;

  // Setup matrix shape for result
  size_t* shape = ALLOC_N(size_t, left->rank);
  for (i = 0; i < left->rank; ++i) shape[i] = left->shape[i];

  // Create result storage.
  result = create_dense_storage(dtype, shape, left->rank, NULL, 0);

  // Do the operation
  DenseElementwiseFuncs[dtype](left->elements, right->elements, result->elements, count_dense_storage_elements(result), op);

  return nm_create(S_DENSE, result);
}


static NMATRIX* elementwise_list_casted(STORAGE_PAIR casted_storage, char op) {
  rb_raise(rb_eNotImpError, "elementwise operations not implemented for list-of-list matrices");
  return NULL;
}


static NMATRIX* elementwise_yale_casted(STORAGE_PAIR casted_storage, char op) {
  YALE_STORAGE *left  = (YALE_STORAGE*)(casted_storage.left),
               *right = (YALE_STORAGE*)(casted_storage.right);
  YALE_STORAGE *result = create_merged_yale_storage(left, right);

  fprintf(stderr, "result: %d, %d\n", result->dtype, result->index_dtype);

  //fprintf(stderr, "Remember to fix elementwise for yale!\n");
  YaleElementwiseFuncs[result->dtype][result->index_dtype](result->shape[0], result->shape[1], op, left->ija, right->ija, result->ija, left->a, right->a, result->a);

  return nm_create(S_YALE, result);
}


nm_elementwise_binary_op_casted_t CastedElementwiseFuncs = {
  elementwise_dense_casted,
  elementwise_list_casted,
  elementwise_yale_casted
};


nm_compare_t EqEqFuncs = {
  dense_storage_eqeq,
  list_storage_eqeq,
  yale_storage_eqeq
};


inline bool numeqeq(const void* x, const void* y, const size_t len, const size_t dtype_size) {
  return (!memcmp(x, y, len * dtype_size));
}


// element eqeq -- like memcmp but handles 0.0 == -0.0 for complex and floating points.
// Second dimension is for hermitians -- complex conjugate. Use 0 for regular equality and 1 for conjugate equality.
nm_eqeq_t ElemEqEq = {
  {NULL, NULL},
  {numeqeq, numeqeq}, // byte
  {numeqeq, numeqeq}, // int8
  {numeqeq, numeqeq}, // int16
  {numeqeq, numeqeq}, // int32
  {numeqeq, numeqeq}, // int64
  {f32eqeq, f32eqeq}, // float32
  {f64eqeq, f64eqeq}, // float64
  {c64eqeq, c64conjeq}, // complex64
  {c128eqeq, c128conjeq}, // complex128
  {numeqeq, numeqeq}, // rational32
  {numeqeq, numeqeq}, // rational64
  {numeqeq, numeqeq}, // rational128
  {numeqeq, numeqeq}  // Ruby object
};


static void nm_delete(NMATRIX* mat) {
  DeleteFuncs[mat->stype](mat->storage);
}

static void nm_delete_ref(NMATRIX* mat) {
  DeleteFuncsRef[mat->stype](mat->storage);
}


static STORAGE* nm_dense_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, size_t init_val_len, VALUE self) {
  return (STORAGE*)(create_dense_storage(dtype, shape, rank, init_val, init_val_len));
}


static STORAGE* nm_list_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, size_t init_val_len, VALUE self) {
  if (init_val_len > 1) {
    rb_raise(rb_eArgError, "list storage needs initial size, not initial value\n");
    return NULL;
  }
  return (STORAGE*)(create_list_storage(dtype, shape, rank, init_val));
}


static STORAGE* nm_yale_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, size_t init_val_len, VALUE self) {
  YALE_STORAGE* s;

  if (init_val_len > 1) {
    rb_raise(rb_eArgError, "list storage needs initial size, not initial value\n");
    return NULL;
  }

  s = create_yale_storage(dtype, shape, rank, *(size_t*)init_val);
  init_yale_storage(s);
  free(init_val);

  if (!s) rb_raise(rb_eNoMemError, "Yale allocation failed");

  return (STORAGE*)(s);
  //return Data_Wrap_Struct(self, NULL, nm_delete, matrix);
}


nm_create_storage_t CreateFuncs = {
  nm_dense_new,
  nm_list_new,
  nm_yale_new
};


nm_cast_copy_storage_t CastCopyFuncs = {
  cast_copy_dense_storage,
  cast_copy_list_storage,
  cast_copy_yale_storage
};



nm_scast_copy_storage_t ScastCopyFuncs = {
  {cast_copy_dense_storage, scast_copy_dense_list, scast_copy_dense_yale},
  {scast_copy_list_dense, cast_copy_list_storage, scast_copy_list_yale},
  {scast_copy_yale_dense, scast_copy_yale_list, cast_copy_yale_storage}
};


nm_stype_ref_t RefFuncs = {
  dense_storage_get,
  list_storage_get,
  yale_storage_ref
};


VALUE nm_dense_set(STORAGE* s, SLICE* slice, VALUE val) {
  void* v = ALLOCA_N(char, nm_sizeof[s->dtype]);
  SetFuncs[s->dtype][NM_ROBJ](1, v, 0, &val, 0);
  dense_storage_set( (DENSE_STORAGE*)s, slice, v );
  return val;
}


// Should work exactly the same as nm_dense_set.
VALUE nm_yale_set(STORAGE* s, SLICE* slice, VALUE val) {
  void* v = ALLOCA_N(char, nm_sizeof[s->dtype]);
  SetFuncs[s->dtype][NM_ROBJ](1, v, 0, &val, 0);
  yale_storage_set( (YALE_STORAGE*)s, slice, v );
  return val;
}


// TODO: Why can't you be more like your brothers, nm_dense_set and nm_yale_set?
VALUE nm_list_set(STORAGE* s, SLICE* slice, VALUE val) {
  void *v = ALLOC_N(char, nm_sizeof[s->dtype]), *rm;
  LIST_STORAGE* ls = (LIST_STORAGE*)s;

  //fprintf(stderr, "    create_val: %p\n", v);

  SetFuncs[s->dtype][NM_ROBJ](1, v, 0, &val, 0);

  if (!memcmp(ls->default_val, v, nm_sizeof[s->dtype])) {
    // User asked to insert default_value, which is actually node *removal*.
    // So let's do that instead.

    rm = list_storage_remove( ls, slice);

    //if (rm) fprintf(stderr, "    remove_val: %p\n", rm);

    if (rm) free(rm);
    return val;

  } else if (list_storage_insert( ls, slice, v ))    return val;
  return Qnil;
  // No need to free; the list keeps v.
}



nm_stype_ins_t InsFuncs = {
  nm_dense_set,
  nm_list_set,
  nm_yale_set,
};



// Converts a typestring to a typecode for storage. Only looks at the first three characters.
int8_t nm_stypestring_to_stype(VALUE str) {
  int8_t i;
  for (i = 0; i < S_TYPES; ++i)
    if ( !strncmp(RSTRING_PTR(str), nm_stypestring[i], 3) ) return i;
  return S_DENSE;
}

int8_t nm_stypesymbol_to_stype(VALUE sym) {
  int8_t i;
  for (i = 0; i < S_TYPES; ++i)
    if (SYM2ID(sym) == rb_intern(nm_stypestring[i])) return i;
  return S_DENSE;
}


int8_t nm_dtypestring_to_dtype(VALUE str) {
  int8_t i;
  for (i = 0; i < NM_TYPES; ++i)
    if ( !strncmp(RSTRING_PTR(str), nm_dtypestring[i], RSTRING_LEN(str)) ) return i;
  return NM_NONE;
}

int8_t nm_dtypesymbol_to_dtype(VALUE sym) {
  int8_t i;
  for (i = 0; i < NM_TYPES; ++i)
    if (SYM2ID(sym) == rb_intern(nm_dtypestring[i])) return i;
  return NM_NONE;
}


// TODO: Probably needs some work for Bignum.
int8_t nm_guess_dtype(VALUE v) {
  switch(TYPE(v)) {
  case T_TRUE:
  case T_FALSE:
    return NM_BYTE;
  case T_STRING:
    if (RSTRING_LEN(v) == 1) return NM_BYTE;
    else                     return NM_NONE;

#if SIZEOF_INT == 8
  case T_FIXNUM:
    return NM_INT64;
  case T_RATIONAL:
    return NM_RATIONAL128;
#else
# if SIZEOF_INT == 4
  case T_FIXNUM:
    return NM_INT32;
  case T_RATIONAL:
    return NM_RATIONAL64;
#else
  case T_FIXNUM:
    return NM_INT16;
  case T_RATIONAL:
    return NM_RATIONAL32;
# endif
#endif

  case T_BIGNUM:
    return NM_INT64;

#if SIZEOF_FLOAT == 4
  case T_COMPLEX:
    return NM_COMPLEX128;
  case T_FLOAT:
    return NM_FLOAT64;
#else
# if SIZEOF_FLOAT == 2
  case T_COMPLEX:
    return NM_COMPLEX64;
  case T_FLOAT:
    return NM_FLOAT32;
# endif
#endif

  case T_ARRAY: // may be passed for dense -- for now, just look at the first element.
    return nm_guess_dtype(RARRAY_PTR(v)[0]);
    // TODO: Look at entire array for most specific type.

  case T_NIL:
  default:
    return NM_NONE;
  }
}

// Used for scasting (changing stype)
inline void cast_copy_value_single(void* to, const void* from, int8_t l_dtype, int8_t r_dtype) {
  if (l_dtype == r_dtype) memcpy(to, from, nm_sizeof[l_dtype]);
  else                    SetFuncs[l_dtype][r_dtype](1, to, 0, from, 0);
}



// Read the shape argument to NMatrix.new, which may be either an array or a single number.
// Second argument is where the shape is stored at the end of the function; returns the rank.
// You are responsible for freeing shape!
size_t* nm_interpret_shape_arg(VALUE arg, size_t* rank) {
  size_t i;
  size_t* shape;

  if (TYPE(arg) == T_ARRAY) {
    *rank = RARRAY_LEN(arg);
    shape = ALLOC_N(size_t, *rank);
    for (i = 0; i < *rank; ++i)
      shape[i]  = (size_t)(FIX2UINT(RARRAY_PTR(arg)[i]));
  } else if (FIXNUM_P(arg)) {
    *rank = 2;
    shape = ALLOC_N(size_t, *rank);
    for (i = 0; i < *rank; ++i)
      shape[i]  = (size_t)(FIX2UINT(arg));
  } else {
    *rank  = 0;
    shape = NULL;
    rb_raise(rb_eArgError, "Expected an array of numbers or a single fixnum for matrix shape");
  }

  return shape;
}


// argv will be either 1 or 2 elements. If 1, could be either initial or dtype. If 2, is initial and dtype.
// This function returns the dtype.
int8_t nm_interpret_dtype(int argc, VALUE* argv, int8_t stype) {
  int offset = 0; // if argc == 1
  if (argc == 2) offset = 1;
  else if (argc != 1) rb_raise(rb_eArgError, "Need an initial value or a dtype");

  if (SYMBOL_P(argv[offset])) return nm_dtypesymbol_to_dtype(argv[offset]);
  else if (IS_STRING(argv[offset])) return nm_dtypestring_to_dtype(StringValue(argv[offset]));
  else if (stype == S_YALE) rb_raise(rb_eArgError, "yale requires dtype");
  else return nm_guess_dtype(argv[0]);

  return NM_NONE;
}

int8_t nm_interpret_stype(VALUE arg) {
  if (SYMBOL_P(arg)) return nm_stypesymbol_to_stype(arg);
  else if (IS_STRING(arg)) return nm_stypestring_to_stype(StringValue(arg));
  else rb_raise(rb_eArgError, "Expected storage type");
  return S_DENSE;
}


void* nm_interpret_initial_value(VALUE arg, int8_t dtype) {
  void* init_val;

  if (TYPE(arg) == T_ARRAY) { // array
    init_val = ALLOC_N(char, nm_sizeof[dtype] * RARRAY_LEN(arg));
    SetFuncs[dtype][NM_ROBJ](RARRAY_LEN(arg), init_val, nm_sizeof[dtype], RARRAY_PTR(arg), nm_sizeof[NM_ROBJ]);
  } else { // single value
    init_val = ALLOC_N(char, nm_sizeof[dtype]);
    SetFuncs[dtype][NM_ROBJ](1, init_val, 0, &arg, 0);
  }

  return init_val;
}


size_t* nm_interpret_initial_capacity(VALUE arg) {
  size_t* init_cap = ALLOC(size_t);
  *init_cap = FIX2UINT(arg);
  return init_cap;
}


/*
 * Create a new NMatrix helper for handling internal ia, ja, and a arguments.
 *
 * This constructor is only called by Ruby code, so we can skip most of the checks.
 */
static VALUE nm_init_yale_from_old_yale(VALUE shape, VALUE dtype, VALUE ia, VALUE ja, VALUE a, VALUE from_dtype, VALUE from_index_dtype, VALUE nm) {
  size_t rank     = 2;
  size_t* shape_  = nm_interpret_shape_arg(shape, &rank);
  int8_t dtype_   = nm_dtypesymbol_to_dtype(dtype);
  char *ia_       = RSTRING_PTR(ia),
       *ja_       = RSTRING_PTR(ja),
       *a_        = RSTRING_PTR(a);
  int8_t from_dtype_ = nm_dtypesymbol_to_dtype(from_dtype);
  int8_t from_index_dtype_ = nm_dtypesymbol_to_dtype(from_index_dtype);
  NMATRIX* nmatrix;

  UnwrapNMatrix( nm, nmatrix );

  nmatrix->stype   = S_YALE;
  nmatrix->storage = (STORAGE*)create_yale_storage_from_old_yale(dtype_, shape_, ia_, ja_, a_, from_dtype_, from_index_dtype_);

  return nm;
}


/*
 * Create a new NMatrix.
 *
 * There are several ways to do this. At a minimum, dimensions and either a dtype or initial values are needed, e.g.,
 *
 *     NMatrix.new(3, :int64)       # square 3x3 dense matrix
 *     NMatrix.new([3,4], :float32) # 3x4 matrix
 *     NMatrix.new(3, 0)            # 3x3 dense matrix initialized to all zeros
 *     NMatrix.new([3,3], [1,2,3])  # [[1,2,3],[1,2,3],[1,2,3]]
 *
 * NMatrix will try to guess the dtype from the first value in the initial values array.
 *
 * You can also provide the stype prior to the dimensions. However, non-dense matrices cannot take initial values, and
 * require a dtype (e.g., :int64):
 *
 *     NMatrix.new(:yale, [4,3], :int64)
 *     NMatrix.new(:list, 5, :rational128)
 *
 * For Yale, you can also give an initial size for the non-diagonal component of the matrix:
 *
 *     NMatrix.new(:yale, [4,3], 2, :int64)
 *
 * Finally, you can be extremely specific, and define a matrix very exactly:
 *
 *     NMatrix.new(:dense, [2,2,2], [0,1,2,3,4,5,6,7], :int8)
 *
 * There is one additional constructor for advanced users, which takes seven arguments and is only for creating Yale matrices
 * with known IA, JA, and A arrays. This is used primarily internally for IO, e.g., reading Matlab matrices, which are
 * stored in old Yale format.
 *
 * Just be careful! There are no overflow warnings in NMatrix.
 */
static VALUE nm_init(int argc, VALUE* argv, VALUE nm) {
  char    ZERO = 0;
  VALUE   QNIL = Qnil;
  int8_t  dtype, stype, offset = 0;
  size_t  rank;
  size_t* shape;
  size_t  init_val_len = 0;
  void*   init_val = NULL;
  NMATRIX* nmatrix;

  // READ ARGUMENTS

  //fprintf(stderr, "Called nmatrix new with %d arguments\n", argc);

  if (argc < 2) { rb_raise(rb_eArgError, "Expected 2-4 arguments (or 8 for internal Yale creation)"); return Qnil; }

  if (!SYMBOL_P(argv[0]) && !IS_STRING(argv[0])) {
    stype  = S_DENSE;
  } else {
    stype  = nm_interpret_stype(argv[0]);                        // 0: String or Symbol
    offset = 1;
  }

  // If there are 7 arguments and Yale, refer to a different init function with fewer sanity checks.
  if (argc == 8) {
    if (stype == S_YALE) return nm_init_yale_from_old_yale(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], nm);
    rb_raise(rb_eArgError, "Expected 2-4 arguments (or 7 for internal Yale creation)");
    return Qnil;
  }

  shape    = nm_interpret_shape_arg(argv[offset], &rank);        // 1: Array or Fixnum
  dtype    = nm_interpret_dtype(argc-1-offset, argv+offset+1, stype); // 2-3: dtype

  if (IS_NUMERIC(argv[1+offset]) || TYPE(argv[1+offset]) == T_ARRAY) { // initial value provided (could also be initial capacity, if yale)
    if (stype == S_YALE) {
      init_val     = nm_interpret_initial_capacity(argv[1+offset]);
      init_val_len = 1;
    } else {
      init_val = nm_interpret_initial_value(argv[1+offset], dtype);// 4: initial value / dtype
      if (TYPE(argv[1+offset]) == T_ARRAY) init_val_len = RARRAY_LEN(argv[1+offset]);
      else                                 init_val_len = 1;
    }
  } else {
    if (stype == S_DENSE) { // no need to initialize dense with any kind of default value unless it's an NM_ROBJ matrix
      if (dtype == NM_ROBJ) { // pretend [nil] was passed for ROBJ.
        init_val = ALLOC(VALUE);
        SetFuncs[NM_ROBJ][NM_ROBJ](1, init_val, 0, &QNIL, 0);
        init_val_len = 1;
      } else init_val = NULL;
    } else if (stype == S_YALE) { // if it's a list or compressed, we want to assume default of 0 even if none provided
      init_val = ALLOC(size_t);
      *(size_t*)init_val = 0;
    } else {
      init_val = ALLOC_N(char, nm_sizeof[dtype]);
      //memset(init_val, 0, nm_sizeof[dtype]); // TODO: See if this works instead of the next line (with NM_ROBJ matrix). Cleaner.
      SetFuncs[dtype][NM_BYTE](1, init_val, 0, &ZERO, 0);
    }
  }


  // TODO: Update to allow an array as the initial value.

  if (dtype == NM_NONE) {
    rb_raise(rb_eArgError, "Could not recognize dtype");
    free(init_val);
    free(shape);
    return nm;
  }

  if (stype < S_TYPES) {
    UnwrapNMatrix( nm, nmatrix );

    nmatrix->stype   = stype;
    nmatrix->storage = CreateFuncs[stype](shape, rank, dtype, init_val, init_val_len, nm);

    return nm;
  } else
    rb_raise(rb_eNotImpError, "Unrecognized storage type");

  free(shape);
  free(init_val);
  return nm;
}


static VALUE nm_alloc(VALUE klass) {
  NMATRIX* mat = ALLOC(NMATRIX);
  mat->storage = NULL;
  mat->stype   = S_TYPES;
  return Data_Wrap_Struct(klass, MarkFuncs[mat->stype], nm_delete, mat);
}


// This is the "back-door initializer," for when Ruby needs to create the object in an atypical way.
//
// Note that objects created this way will have NULL storage.
/*static VALUE nm_initialize(VALUE self, VALUE stype, VALUE dtype) {
  NMATRIX* matrix;
  UnwrapNMatrix(self, matrix);

  matrix->stype   = nm_interpret_stype(stype);
  matrix->dtype   = nm_interpret_dtype(1, &dtype, stype);
  matrix->storage = NULL;

  return self;
}*/


static VALUE nm_init_copy(VALUE copy, VALUE original) {
  NMATRIX *lhs, *rhs;

  CheckNMatrixType(original);

  if (copy == original) return copy;

  UnwrapNMatrix( original, rhs );
  UnwrapNMatrix( copy,     lhs );

  lhs->stype = rhs->stype;

  // Copy the storage
  lhs->storage = CastCopyFuncs[rhs->stype](rhs->storage, rhs->storage->dtype);

  return copy;
}


static VALUE nm_init_cast_copy(VALUE copy, VALUE original, VALUE new_dtype_symbol) {
  NMATRIX *lhs, *rhs;
  int8_t new_dtype = nm_dtypesymbol_to_dtype(new_dtype_symbol);
  //fprintf(stderr,"In copy constructor\n");

  CheckNMatrixType(original);

  if (copy == original) return copy;

  UnwrapNMatrix( original, rhs );
  UnwrapNMatrix( copy,     lhs );

  lhs->stype = rhs->stype;

  // Copy the storage
  lhs->storage = CastCopyFuncs[rhs->stype](rhs->storage, new_dtype);

  return copy;
}


/*
 * Create a copy of an NMatrix with a different dtype. See also cast.
 */
 // TODO: Deprecate this function and farm it out to scast_copy. as_dtype will still work, but it'll be in pure Ruby and
 // just use ::cast instead.
static VALUE nm_cast_copy(VALUE self, VALUE new_dtype_symbol) {
  NMATRIX *original, *copy;
  int8_t new_dtype = nm_dtypesymbol_to_dtype(new_dtype_symbol);

  CheckNMatrixType(self);

  UnwrapNMatrix(self, original);

  copy = ALLOC(NMATRIX);
  copy->stype = original->stype;
  copy->storage = CastCopyFuncs[original->stype](original->storage, new_dtype);

  return Data_Wrap_Struct(cNMatrix, MarkFuncs[copy->stype], nm_delete, copy);
}


/*
 * Create a copy of an NMatrix with a different stype and dtype. See also cast.
 *
 *     m.cast(:dense, :int64)
 *
 */
static VALUE nm_scast_copy(VALUE self, VALUE new_stype_symbol, VALUE new_dtype_symbol) {
  NMATRIX* original, *copy;
  int8_t new_dtype = nm_dtypesymbol_to_dtype(new_dtype_symbol);
  int8_t new_stype = nm_stypesymbol_to_stype(new_stype_symbol);

  CheckNMatrixType(self);

  UnwrapNMatrix(self, original);

  copy = ALLOC(NMATRIX);
  copy->stype = new_stype;

  // Copy and scast the storage.
  if (new_stype == original->stype) copy->storage = CastCopyFuncs[original->stype](original->storage, new_dtype);
  else                              copy->storage = ScastCopyFuncs[copy->stype][original->stype](original->storage, new_dtype);

  return Data_Wrap_Struct(cNMatrix, MarkFuncs[copy->stype], nm_delete, copy);
}



// Cast a single matrix to a new dtype (unless it's already casted, then just return it). Helper for binary_storage_cast_alloc.
static inline STORAGE* storage_cast_alloc(NMATRIX* matrix, int8_t new_dtype) {
  if (matrix->storage->dtype == new_dtype) return matrix->storage;
  else                                     return CastCopyFuncs[matrix->stype](matrix->storage, new_dtype);
}


// Cast a pair of matrices for a binary operation to a new dtype (which this function determines). Technically, only
// does an actual cast on matrices that are the wrong dtype; otherwise returns a reference to the original. Bear this in
// mind when freeing memory!
static inline STORAGE_PAIR binary_storage_cast_alloc(NMATRIX* left_matrix, NMATRIX* right_matrix) {
  STORAGE_PAIR casted;
  int8_t new_dtype = Upcast[left_matrix->storage->dtype][right_matrix->storage->dtype];

  casted.left  = storage_cast_alloc(left_matrix, new_dtype);
  casted.right = storage_cast_alloc(right_matrix, new_dtype);

  return casted;
}

/*
 * Equality operator. Returns a single true or false value indicating whether the matrices are equivalent.
 *
 * For elementwise, use == instead.
 *
 * This method will raise an exception if dimensions do not match.
 */
static VALUE nm_eqeq(VALUE left, VALUE right) {
  bool result;
  NMATRIX *l, *r;
  STORAGE_PAIR casted;

  CheckNMatrixType(left);
  CheckNMatrixType(right);

  UnwrapNMatrix(left, l);
  UnwrapNMatrix(right, r);

  if (l->stype != r->stype) //rb_raise(nm_eStorageTypeError, "wrong storage type");
    rb_raise(rb_eNotImpError, "comparison between different matrix stypes not yet implemented");

  casted = binary_storage_cast_alloc(l, r);

  result = EqEqFuncs[l->stype](casted.left, casted.right);

  // Free any casted-storage we created for the comparison.
  // TODO: Can we make the Ruby GC take care of this stuff now that we're using it?
  // If we did that, we night not have to re-create these every time, right? Or wrong? Need to do
  // more research.
  if (l->storage != casted.left)   DeleteFuncs[l->stype](casted.left);
  if (r->storage != casted.right)  DeleteFuncs[l->stype](casted.right);

  return result ? Qtrue : Qfalse;
}


static VALUE multiply_matrix(NMATRIX* left, NMATRIX* right) {
  ///TODO: multiplication for non-dense and/or non-decimal matrices
  size_t*  resulting_shape   = ALLOC_N(size_t, 2);
  NMATRIX* result;
  bool     vector = false;

  // Make sure both of our matrices are of the correct type.
  STORAGE_PAIR casted = binary_storage_cast_alloc(left, right);

  resulting_shape[0] = left->storage->shape[0];
  resulting_shape[1] = right->storage->shape[1];

  // Sometimes we only need to use matrix-vector multiplication (e.g., GEMM versus GEMV). Find out.
  if (resulting_shape[1] == 1) vector = true;

  result = CastedMultiplyFuncs[left->stype](casted, resulting_shape, vector);

  // Free any casted-storage we created for the multiplication.
  // TODO: Can we make the Ruby GC take care of this stuff now that we're using it?
  // If we did that, we night not have to re-create these every time, right? Or wrong? Need to do
  // more research.
  if (left->storage != casted.left)   DeleteFuncs[left->stype](casted.left);
  if (right->storage != casted.right) DeleteFuncs[left->stype](casted.right);

  if (result) return Data_Wrap_Struct(cNMatrix, MarkFuncs[result->stype], nm_delete, result);
  return Qnil; // Only if we try to multiply list matrices should we return Qnil.
}


static VALUE multiply_scalar(NMATRIX* left, VALUE scalar) {
  rb_raise(rb_eNotImpError, "matrix-scalar multiplication not implemented yet");
  return Qnil;
}


/*
 * Matrix multiply (dot product): against another matrix or a vector.
 *
 * For elementwise, use * instead.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 */
static VALUE nm_multiply(VALUE left_v, VALUE right_v) {
  NMATRIX *left, *right;

  // left has to be of type NMatrix.
  CheckNMatrixType(left_v);

  UnwrapNMatrix( left_v, left );

  if (IS_NUMERIC(right_v))
    return multiply_scalar(left, right_v);

  else if (TYPE(right_v) == T_ARRAY)
    rb_raise(rb_eNotImpError, "for matrix-vector multiplication, please use an NVector instead of an Array for now");

  //if (RDATA(right_v)->dfree != (RUBY_DATA_FUNC)nm_delete) {
  else if (TYPE(right_v) == T_DATA && RDATA(right_v)->dfree == (RUBY_DATA_FUNC)nm_delete) { // both are matrices
    UnwrapNMatrix( right_v, right );

    if (left->storage->shape[1] != right->storage->shape[0])
      rb_raise(rb_eArgError, "incompatible dimensions");

    if (left->stype != right->stype)
      rb_raise(rb_eNotImpError, "matrices must have same stype");

    return multiply_matrix(left, right);

  } else rb_raise(rb_eTypeError, "expected right operand to be NMatrix, NVector, or single numeric value");

  return Qnil;
}


static VALUE nm_elementwise(VALUE leftv, VALUE rightv, char op) {
  ///TODO: multiplication for non-dense and/or non-decimal matrices
  NMATRIX *result, *left, *right;
  STORAGE_PAIR casted;

  CheckNMatrixType(leftv);
  CheckNMatrixType(rightv);

  UnwrapNMatrix(rightv, right);
  UnwrapNMatrix(leftv, left);

  // Make sure both of our matrices are of the correct type.
  casted = binary_storage_cast_alloc(left, right);

  result = CastedElementwiseFuncs[left->stype](casted, op);

  // Free up temporary casted matrices
  if (left->storage != casted.left)   DeleteFuncs[left->stype](casted.left);
  if (right->storage != casted.right) DeleteFuncs[left->stype](casted.right);

  if (result) return Data_Wrap_Struct(cNMatrix, MarkFuncs[result->stype], nm_delete, result);
  return Qnil; // Only if we try to multiply list matrices should we return Qnil.
}


/*
 * Matrix element-wise addition.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 */
static VALUE nm_ew_add(VALUE left, VALUE right) {
  return nm_elementwise(left, right, '+');
}

/*
 * Matrix element-wise subtraction.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 */
static VALUE nm_ew_subtract(VALUE left, VALUE right) {
  return nm_elementwise(left, right, '-');
}

/*
 * Matrix element-wise multiplication.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * For dot product, use +dot+ instead.
 */
static VALUE nm_ew_multiply(VALUE left, VALUE right) {
  return nm_elementwise(left, right, '*');
}

/*
 * Matrix element-wise division.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 */
static VALUE nm_ew_divide(VALUE left, VALUE right) {
  return nm_elementwise(left, right, '/');
}


/*
 * Matrix element-wise comparison (equality) operator.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * Note that the matrix returned will be of the same dtype as the upcast of the input matrices. If that's not what you
 * want, use +cast+.
 */
static VALUE nm_ew_eqeq(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_OP_EQEQ);
}

/*
 * Matrix element-wise less-than-or-equals operator.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * Note that the matrix returned will be of the same dtype as the upcast of the input matrices. If that's not what you
 * want, use +cast+.
 */
static VALUE nm_ew_leq(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_OP_LTE);
}


/*
 * Matrix element-wise greater-than-or-equals operator.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * Note that the matrix returned will be of the same dtype as the upcast of the input matrices. If that's not what you
 * want, use +cast+.
 */
static VALUE nm_ew_geq(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_OP_GTE);
}


/*
 * Matrix element-wise strictly-less-than operator.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * Note that the matrix returned will be of the same dtype as the upcast of the input matrices. If that's not what you
 * want, use +cast+.
 */
static VALUE nm_ew_lt(VALUE left, VALUE right) {
  return nm_elementwise(left, right, '<');
}


/*
 * Matrix element-wise strictly-greater-than operator.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * Note that the matrix returned will be of the same dtype as the upcast of the input matrices. If that's not what you
 * want, use +cast+.
 */
static VALUE nm_ew_gt(VALUE left, VALUE right) {
  return nm_elementwise(left, right, '>');
}


/*
 * Matrix element-wise inequality operator.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * Note that the matrix returned will be of the same dtype as the upcast of the input matrices. If that's not what you
 * want, use +cast+.
 */
static VALUE nm_ew_neq(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_OP_NEQ);
}


// Borrowed this function from NArray. Handles 'each' iteration on a dense matrix.
//
// Additionally, handles separately matrices containing VALUEs and matrices containing
// other types of data.
static VALUE nm_dense_each(VALUE nmatrix) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)(NM_STORAGE(nmatrix));
  VALUE v;
  size_t i;

  void (*copy)();

  if (NM_DTYPE(nmatrix) == NM_ROBJ) {

    // matrix of Ruby objects -- yield directly
    for (i = 0; i < count_dense_storage_elements(s); ++i)
      rb_yield( *((VALUE*)((char*)(s->elements) + i*nm_sizeof[NM_DTYPE(nmatrix)])) );

  } else {
    // We're going to copy the matrix element into a Ruby VALUE and then operate on it.
    copy = SetFuncs[NM_ROBJ][NM_DTYPE(nmatrix)];

    for (i = 0; i < count_dense_storage_elements(s); ++i) {
      (*copy)(1, &v, 0, (char*)(s->elements) + i*nm_sizeof[NM_DTYPE(nmatrix)], 0);
      rb_yield(v); // yield to the copy we made
    }
  }

  return nmatrix;
}


/*
 * Iterate over the matrix as you would an Enumerable (e.g., Array).
 *
 * Currently only works for dense.
 */
static VALUE nm_each(VALUE nmatrix) {
  volatile VALUE nm = nmatrix; // not sure why we do this, but it gets done in ruby's array.c.

  switch(NM_STYPE(nm)) {
  case S_DENSE:
    return nm_dense_each(nm);
  default:
    rb_raise(rb_eNotImpError, "only dense each works right now");
  }
}


// Does not create storage, but does destroy it.
NMATRIX* nm_create(int8_t stype, void* storage) {
  NMATRIX* mat = ALLOC(NMATRIX);

  mat->stype   = stype;
  mat->storage = storage;

  return mat;
}


static SLICE* get_slice(size_t rank, VALUE* c, VALUE self) {
  size_t r;
  VALUE beg, end;
  int i;

  SLICE* slice = ALLOC(SLICE);
  slice->coords = ALLOC_N(size_t,rank);
  slice->lens = ALLOC_N(size_t, rank);
  slice->is_one_el = 1;

  for (r = 0; r < rank; ++r) {
    VALUE cl = CLASS_OF(c[r]);
    if (cl == rb_cFixnum) {
      slice->coords[r] = FIX2UINT(c[r]);
      slice->lens[r] = 1;
    }
    else if (cl == rb_cRange) {
      rb_range_values(c[r], &beg, &end, &i);
      slice->coords[r] = FIX2UINT(beg);
      slice->lens[r] = FIX2UINT(end) - slice->coords[r] + 1;
      slice->is_one_el = 0;
    }
    else rb_raise(rb_eArgError, "%s is not used for slicing", rb_class2name(cl));

    if (slice->coords[r] + slice->lens[r] > NM_SHAPE(self,r)) 
      rb_raise(rb_eArgError, "out of range");
  }

  return slice;
}


/*
 * Access the contents of an NMatrix at given coordinates.
 *
 *     n[3,3]  # => 5.0
 *
 */
VALUE nm_mref(int argc, VALUE* argv, VALUE self) {
  NMATRIX* mat;
  SLICE* slice;  
  void* v;

  if (NM_RANK(self) == (size_t)(argc)) {
    slice = get_slice((size_t)(argc), argv, self);
    // TODO: Slice for List, Yale types
    if (NM_STYPE(self) == S_DENSE && slice->is_one_el == 0) {

      mat = ALLOC(NMATRIX);
      mat->stype = S_DENSE; 
      mat->storage = RefFuncs[NM_STYPE(self)](NM_STORAGE(self), slice);
      return Data_Wrap_Struct(cNMatrix, MarkFuncs[mat->stype], nm_delete_ref, mat);
    } 
    else {
      SetFuncs[NM_ROBJ][NM_DTYPE(self)](1, &v, 0,
                RefFuncs[NM_STYPE(self)](NM_STORAGE(self), slice), 0);
      return v;
    }
                              
  } else if (NM_RANK(self) < (size_t)(argc)) {
    rb_raise(rb_eArgError, "Coordinates given exceed matrix rank");
  } else {
    rb_raise(rb_eNotImpError, "This type slicing not supported yet");
  }
  return Qnil;
}


/*
 * Modify the contents of an NMatrix in the given cell
 *
 *     n[3,3] = 5.0
 *
 * Also returns the new contents, so you can chain:
 *
 *     n[3,3] = n[2,3] = 5.0
 */
VALUE nm_mset(int argc, VALUE* argv, VALUE self) {
  size_t rank = argc - 1; // last arg is the value

  if (argc <= 1) {
    rb_raise(rb_eArgError, "Expected coordinates and r-value");

  } else if (NM_RANK(self) == rank) {
    return (*(InsFuncs[NM_STYPE(self)]))( NM_STORAGE(self),
                                         get_slice(rank, argv, self),
                                         argv[rank] );

  } else if (NM_RANK(self) < rank) {
    rb_raise(rb_eArgError, "Coordinates given exceed matrix rank");
  } else {
    rb_raise(rb_eNotImpError, "Slicing not supported yet");
  }
  return Qnil;
}


/*
 * Get the rank of an NMatrix (the number of dimensions).
 *
 * In other words, if you set your matrix to be 3x4, the rank is 2. If the matrix was initialized as 3x4x3, the rank
 * is 3.
 *
 * This function may lie slightly for NVectors, which are internally stored as rank 2 (and have an orientation), but
 * act as if they're rank 1.
 */
VALUE nm_rank(VALUE self) {
  VALUE ret;
  SetFuncs[NM_ROBJ][NM_INT64]( 1, &ret, 0, &(NM_STORAGE(self)->rank), 0 );
  return ret;
}


/*
 * Get the shape (dimensions) of a matrix.
 */
VALUE nm_shape(VALUE self) {
  STORAGE* s   = NM_STORAGE(self);

  // Copy elements into a VALUE array and then use those to create a Ruby array with rb_ary_new4.
  VALUE* shape = ALLOCA_N(VALUE, s->rank);
  SetFuncs[NM_ROBJ][NM_SIZE_T]( s->rank, shape, sizeof(VALUE), s->shape, sizeof(size_t));

  return rb_ary_new4(s->rank, shape);
}


/*
 * Get the storage type (stype) of a matrix, e.g., :yale, :dense, or :list.
 */
static VALUE nm_stype(VALUE self) {
  ID stype = rb_intern(nm_stypestring[NM_STYPE(self)]);
  return ID2SYM(stype);
}


/*
 * Get the data type (dtype) of a matrix, e.g., :byte, :int8, :int16, :int32, :int64, :float32, :float64, :complex64,
 * :complex128, :rational32, :rational64, :rational128, or :object (the last is a Ruby object).
 */
static VALUE nm_dtype(VALUE self) {
  ID dtype = rb_intern(nm_dtypestring[NM_DTYPE(self)]);
  return ID2SYM(dtype);
}


/* Interprets cblas argument which could be any of false/:no_transpose, :transpose, or :complex_conjugate,
 * into an enum recognized by cblas.
 *
 * Called by nm_cblas_gemm -- basically inline.
 *
 */
static char gemm_op_sym(VALUE op) {
  if (op == false || rb_to_id(op) == nm_id_no_transpose) return CblasNoTrans;
  else if (rb_to_id(op) == nm_id_transpose) return CblasTrans;
  else if (rb_to_id(op) == nm_id_complex_conjugate) return CblasConjTrans;
  else rb_raise(rb_eArgError, "Expected false, :transpose, or :complex_conjugate");
  return CblasNoTrans;
}


/* Call any of the cblas_xgemm functions as directly as possible.
 *
 * The cblas_xgemm functions (dgemm, sgemm, cgemm, and zgemm) define the following operation:
 *
 *    C = alpha*op(A)*op(B) + beta*C
 *
 * where op(X) is one of <tt>op(X) = X</tt>, <tt>op(X) = X**T</tt>, or the complex conjugate of X.
 *
 * Note that this will only work for dense matrices that are of types :float32, :float64, :complex64, and :complex128.
 * Other types are not implemented in BLAS, and while they exist in NMatrix, this method is intended only to
 * expose the ultra-optimized ATLAS versions.
 *
 * == Arguments
 * See: http://www.netlib.org/blas/dgemm.f
 *
 * You probably don't want to call this function. Instead, why don't you try cblas_gemm, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_cblas_gemm(VALUE self,
                           VALUE trans_a, VALUE trans_b,
                           VALUE m, VALUE n, VALUE k,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb,
                           VALUE beta,
                           VALUE c, VALUE ldc)
{
  struct cblas_param_t p = cblas_params_for_multiply(((DENSE_STORAGE*)(NM_STORAGE(a))), ((DENSE_STORAGE*)(NM_STORAGE(b))), ((DENSE_STORAGE*)(NM_STORAGE(c))), false);
  p.M = FIX2INT(m);
  p.N = FIX2INT(n);
  p.K = FIX2INT(k);

  p.lda = FIX2INT(lda);
  p.ldb = FIX2INT(ldb);
  p.ldc = FIX2INT(ldc);

  switch(NM_DTYPE(c)) {
  case NM_FLOAT32:
  case NM_FLOAT64:
    p.alpha.d[0] = NUM2DBL(alpha);
    p.beta.d[0] = NUM2DBL(beta);
    break;

  case NM_COMPLEX64:
    p.alpha.c[0].r = REAL2DBL(alpha);
    p.alpha.c[0].i = IMAG2DBL(alpha);
    p.beta.c[0].r = REAL2DBL(beta);
    p.beta.c[0].i = IMAG2DBL(beta);
    break;

  case NM_COMPLEX128:
    p.alpha.z.r = REAL2DBL(alpha);
    p.alpha.z.i = IMAG2DBL(alpha);
    p.beta.z.r = REAL2DBL(beta);
    p.beta.z.i = IMAG2DBL(beta);
    break;

  case NM_BYTE:
    p.alpha.b[0] = FIX2INT(alpha);
    p.beta.b[0] = FIX2INT(beta);
    break;

  case NM_INT8:
  case NM_INT16:
  case NM_INT32:
  case NM_INT64:
    p.alpha.i[0] = FIX2INT(alpha);
    p.beta.i[0]  = FIX2INT(beta);
    break;

  case NM_RATIONAL32:
    p.alpha.r[0].n = NUMER2INT(alpha);
    p.alpha.r[0].d = DENOM2INT(alpha);
    p.beta.r[0].n = NUMER2INT(beta);
    p.beta.r[0].d = DENOM2INT(beta);
    break;

  case NM_RATIONAL64:
    p.alpha.ra[0].n = NUMER2INT(alpha);
    p.alpha.ra[0].d = DENOM2INT(alpha);
    p.beta.ra[0].n = NUMER2INT(beta);
    p.beta.ra[0].d = DENOM2INT(beta);
    break;

  case NM_RATIONAL128:
    p.alpha.rat.n = NUMER2INT(alpha);
    p.alpha.rat.d = DENOM2INT(alpha);
    p.beta.rat.n = NUMER2INT(beta);
    p.beta.rat.d = DENOM2INT(beta);
    break;

  case NM_ROBJ:
    p.alpha.v[0] = alpha;
    p.beta.v[0]  = beta;
    break;

  default:
    rb_raise(nm_eDataTypeError, "unexpected dtype");

  }

  /* fprintf(stderr, "cblas_gemm: %d %d %d %d %d %f %d %d %f %d\n", trans_a_, trans_b_,
         m_, n_, k_, alpha_, lda_, ldb_, beta_, ldc_); */

  GemmFuncs[NM_DTYPE(c)](CblasRowMajor, gemm_op_sym(trans_a), gemm_op_sym(trans_b), p);

  return Qtrue;
}


/* Call any of the cblas_xgemv functions as directly as possible.
 *
 * The cblas_xgemv functions (dgemv, sgemv, cgemv, and zgemv) define the following operation:
 *
 *    y = alpha*op(A)*x + beta*y
 *
 * where op(A) is one of <tt>op(A) = A</tt>, <tt>op(A) = A**T</tt>, or the complex conjugate of A.
 *
 * Note that this will only work for dense matrices that are of types :float32, :float64, :complex64, and :complex128.
 * Other types are not implemented in BLAS, and while they exist in NMatrix, this method is intended only to
 * expose the ultra-optimized ATLAS versions.
 *
 * == Arguments
 * See: http://www.netlib.org/blas/dgemm.f
 *
 * You probably don't want to call this function. Instead, why don't you try cblas_gemv, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_cblas_gemv(VALUE self,
                           VALUE trans_a,
                           VALUE m, VALUE n,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE x, VALUE incx,
                           VALUE beta,
                           VALUE y, VALUE incy) {

  struct cblas_param_t p;
  p.M = FIX2INT(m);
  p.N = FIX2INT(n);
  p.A = ((DENSE_STORAGE*)(NM_STORAGE(a)))->elements;
  p.B = ((DENSE_STORAGE*)(NM_STORAGE(x)))->elements;
  p.C = ((DENSE_STORAGE*)(NM_STORAGE(y)))->elements;
  p.lda = FIX2INT(lda);
  p.ldb = FIX2INT(incx);
  p.ldc = FIX2INT(incy);

  switch(NM_DTYPE(y)) {
  case NM_FLOAT32:
  case NM_FLOAT64:
    p.alpha.d[0] = REAL2DBL(alpha);
    p.beta.d[0]  = REAL2DBL(beta);
    break;
  case NM_COMPLEX64:
    p.alpha.c[0].r = REAL2DBL(alpha);
    p.alpha.c[0].i = IMAG2DBL(alpha);
    p.beta.c[0].r = REAL2DBL(beta);
    p.beta.c[0].i = IMAG2DBL(beta);
    break;
  case NM_COMPLEX128:
    p.alpha.z.r = REAL2DBL(alpha);
    p.alpha.z.i = IMAG2DBL(alpha);
    p.beta.z.r = REAL2DBL(beta);
    p.beta.z.i = IMAG2DBL(beta);
    break;
  }

  /* fprintf(stderr, "cblas_gemm: %d %d %d %d %d %f %d %d %f %d\n", trans_a_, trans_b_,
         m_, n_, k_, alpha_, lda_, ldb_, beta_, ldc_); */

  GemvFuncs[NM_DTYPE(y)](CblasRowMajor, gemm_op_sym(trans_a), p);

  return Qtrue;
}


/*
 * Find the capacity of an NMatrix. The capacity only differs from the size for Yale matrices, which occasionally
 * allocate more space than they need. For list and dense, capacity gives the number of elements in the matrix.
 */
static VALUE nm_capacity(VALUE self) {
  VALUE cap;

  switch(NM_STYPE(self)) {
  case S_YALE:
    cap = UINT2NUM(((YALE_STORAGE*)(NM_STORAGE(self)))->capacity);
    break;

  case S_DENSE:
    cap = UINT2NUM(count_dense_storage_elements( (DENSE_STORAGE*)(NM_STORAGE(self)) ));
    break;

  case S_LIST:
    cap = UINT2NUM(count_list_storage_elements( (LIST_STORAGE*)(NM_STORAGE(self)) ));
    break;

  default:
    //rb_raise(rb_eNotImpError, "TODO: implement capacity/size on other storage types");
    rb_raise(nm_eStorageTypeError, "unrecognized stype");
  }

  return cap;
}


/*
 * Get the size of a Yale matrix (the number of elements actually stored).
 *
 * For capacity (the maximum number of elements that can be stored without a resize), use capacity instead.
 */
static VALUE nm_yale_size(VALUE self) {
  VALUE sz;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  SetFuncs[NM_ROBJ][s->index_dtype](1, &sz, 0, (YALE_SIZE_PTR((s), nm_sizeof[s->index_dtype])), 0);
  return sz;
}


/*
 * Get the A array of a Yale matrix (which stores the diagonal and the LU portions of the matrix).
 */
static VALUE nm_yale_a(VALUE self) {
  y_size_t sz, i;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*sz);

  SetFuncs[NM_ROBJ][s->dtype](sz, vals, nm_sizeof[NM_ROBJ], s->a, nm_sizeof[s->dtype]);
  ary = rb_ary_new4(sz, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}


/*
 * Get the diagonal ("D") portion of the A array of a Yale matrix.
 */
static VALUE nm_yale_d(VALUE self) {
  y_size_t sz;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*s->shape[0]);

  SetFuncs[NM_ROBJ][s->dtype](s->shape[0], vals, nm_sizeof[NM_ROBJ], s->a, nm_sizeof[s->dtype]);
  ary = rb_ary_new4(s->shape[0], vals);

  return ary;
}


/*
 * Get the non-diagonal ("LU") portion of the A array of a Yale matrix.
 */
static VALUE nm_yale_lu(VALUE self) {
  y_size_t sz, i;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*(s->capacity - s->shape[0]));

  SetFuncs[NM_ROBJ][s->dtype](sz - s->shape[0] - 1, vals, nm_sizeof[NM_ROBJ], (char*)(s->a) + (s->shape[0] + 1)*nm_sizeof[s->dtype], nm_sizeof[s->dtype]);
  ary = rb_ary_new4(sz - s->shape[0] - 1, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}


/*
 * Get the IA portion of the IJA array of a Yale matrix. This gives the start and end positions of rows in the
 * JA and LU portions of the IJA and A arrays, respectively.
 */
static VALUE nm_yale_ia(VALUE self) {
  y_size_t sz;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*(s->shape[0]+1));

  SetFuncs[NM_ROBJ][s->index_dtype](s->shape[0]+1, vals, nm_sizeof[NM_ROBJ], s->ija, nm_sizeof[s->index_dtype]);
  ary = rb_ary_new4(s->shape[0]+1, vals);

  return ary;
}


/*
 * Get the JA portion of the IJA array of a Yale matrix. This gives the column indices for entries in corresponding
 * positions in the LU portion of the A array.
 */
static VALUE nm_yale_ja(VALUE self) {
  y_size_t sz, i;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*(s->capacity - s->shape[0]));

  SetFuncs[NM_ROBJ][s->index_dtype](sz - s->shape[0] - 1, vals, nm_sizeof[NM_ROBJ], (char*)(s->ija) + (s->shape[0] + 1)*nm_sizeof[s->index_dtype], nm_sizeof[s->index_dtype]);
  ary = rb_ary_new4(sz - s->shape[0] - 1, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}


/*
 * Get the IJA array of a Yale matrix.
 */
static VALUE nm_yale_ija(VALUE self) {
  y_size_t sz, i;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*s->capacity);

  SetFuncs[NM_ROBJ][s->index_dtype](sz, vals, nm_sizeof[NM_ROBJ], s->ija, nm_sizeof[s->index_dtype]);
  ary = rb_ary_new4(sz, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}



/*
 * Calculate the exact determinant of a dense matrix.
 *
 * Returns nil for dense matrices which are not square or rank other than 2.
 *
 * Note: Currently only implemented for 2x2 and 3x3 matrices.
 */
static VALUE nm_det_exact(VALUE self) {
  NMATRIX* m;
  void* result = ALLOCA_N(char, nm_sizeof[NM_DTYPE(self)]);
  VALUE ret;

  if (NM_STYPE(self) != S_DENSE) rb_raise(nm_eStorageTypeError, "can only calculate exact determinant for dense matrices");

  UnwrapNMatrix(self, m);
  if (m->storage->rank != 2 || m->storage->shape[0] != m->storage->shape[1]) return Qnil;

  // Calculate the determinant and then assign it to the return value
  DenseDetExact[NM_DTYPE(self)](m->storage->shape[0], ((DENSE_STORAGE*)(m->storage))->elements, m->storage->shape[0], result);
  SetFuncs[NM_ROBJ][NM_DTYPE(self)](1, &ret, 0, result, 0);

  return ret;
}


// This is probably faster and smaller than writing an array of transpose functions. But if you want to see what it would look like,
// see transp.template.c (not the yale one).
//
// Note that this is a copy-transpose. In-place transpose is a whole different operation and bag of worms.
static void dense_transpose_generic(const unsigned int M, const unsigned int N, const char* A, const int lda, char* B, const int ldb, size_t dtype_size) {
  unsigned int i, j;

  for (i = 0; i < N; ++i) {
    for (j = 0; j < M; ++j) {
      memcpy(B + (i*ldb+j)*dtype_size, A + (j*lda+i)*dtype_size, dtype_size);
    }
  }
}


static NMATRIX* transpose_new_dense(NMATRIX* self_m, size_t* shape) {
  NMATRIX* result = nm_create(S_DENSE, create_dense_storage(self_m->storage->dtype, shape, 2, NULL, 0));

  dense_transpose_generic(self_m->storage->shape[0],
                          self_m->storage->shape[1],
                          ((DENSE_STORAGE*)(self_m->storage))->elements,
                          self_m->storage->shape[1],
                          ((DENSE_STORAGE*)(result->storage))->elements,
                          result->storage->shape[1],
                          nm_sizeof[self_m->storage->dtype]);
  return result;
}


static NMATRIX* transpose_new_yale(NMATRIX* self_m, size_t* shape) {
  YALE_PARAM A, B;
  NMATRIX* result;
  size_t sz;

  YaleGetSize(sz, (YALE_STORAGE*)(self_m->storage)); // size of new matrix is going to be size of old matrix
  result = nm_create(S_YALE, create_yale_storage(self_m->storage->dtype, shape, 2, sz));

  // TODO: Do we really need to initialize the whole thing? Or just the A portion?
  init_yale_storage((YALE_STORAGE*)(result->storage));

  A.ia = A.ja = ((YALE_STORAGE*)(self_m->storage))->ija;
  B.ia = B.ja = ((YALE_STORAGE*)(result->storage))->ija;
  A.a  = ((YALE_STORAGE*)(self_m->storage))->a;
  B.a  = ((YALE_STORAGE*)(result->storage))->a;
  A.diag = true;

  // call the appropriate function pointer
  SparseTransposeFuncs[ self_m->storage->dtype ][ ((YALE_STORAGE*)(self_m->storage))->index_dtype ](shape[0], shape[1], A, B, true);

  return result;
}


static NMATRIX* transpose_new_err(NMATRIX* self_m, size_t* shape) {
  free(shape);
  rb_raise(rb_eNotImpError, "no transpose written for this type");
  return self_m;
}

nm_transpose_t TransposeFuncs = {
  transpose_new_dense,
  transpose_new_err,
  transpose_new_yale
};

/*
 * Transform the matrix (in-place) to its complex conjugate. Only works on complex matrices.
 */
static VALUE nm_complex_conjugate_bang(VALUE self) {
  NMATRIX* m;
  void* elem;
  size_t sz, p;

  UnwrapNMatrix(self, m);

  if (m->stype == S_DENSE) sz = count_storage_max_elements(m->storage);
  else if (m->stype == S_YALE) YaleGetSize(sz, m->storage);
  else rb_raise(rb_eNotImpError, "please cast to yale or dense (complex) first");

  elem = m->storage->elements; // this gets A array or elements array from dense and yale

  // Walk through and negate the imaginary component
  if (NM_DTYPE(self) == NM_COMPLEX64) {
    for (p = 0; p < sz; ++p)      ((complex64*)elem)[p].i = -((complex64*)elem)[p].i;
  } else if (NM_DTYPE(self) == NM_COMPLEX128) {
    for (p = 0; p < sz; ++p)      ((complex128*)elem)[p].i = -((complex128*)elem)[p].i;
  } else {
    rb_raise(rb_eNotImpError, "can only calculate in-place complex conjugate on matrices of type :complex64 or :complex128");
  }

  return self;
}


/*
 * Create a transposed copy of this matrix.
 */
static VALUE nm_transpose_new(VALUE self) {
  NMATRIX *self_m, *result;
  size_t* shape   = ALLOC_N(size_t, 2);

  UnwrapNMatrix( self, self_m );

  // switch the dimensions
  shape[1] = self_m->storage->shape[0];
  shape[0] = self_m->storage->shape[1];

  result = TransposeFuncs[self_m->stype](self_m, shape);

  return Data_Wrap_Struct(cNMatrix, MarkFuncs[result->stype], nm_delete, result);
}


/*
 * Given a binary operation between types t1 and t2, what type will be returned?
 */
static VALUE nm_upcast(VALUE self, VALUE t1, VALUE t2) {
  int8_t dtype = Upcast[nm_dtypesymbol_to_dtype(t1)][nm_dtypesymbol_to_dtype(t2)];

  // The actual Upcast table returns NM_NONE if the types are unrecognized. If NM_NONE is the result, nil will be
  // returned instead.
  if (dtype == NM_NONE) return Qnil;

  return ID2SYM(rb_intern(nm_dtypestring[dtype]));
}


// Helper function for nm_symmetric and nm_hermitian.
static VALUE is_symmetric(VALUE self, bool hermitian) {
  NMATRIX* m;
  UnwrapNMatrix(self, m);

  if (m->storage->shape[0] == m->storage->shape[1] && m->storage->rank == 2) {

    if (NM_STYPE(self) == S_DENSE) {
      if (dense_is_symmetric((DENSE_STORAGE*)(m->storage), m->storage->shape[0], hermitian)) return Qtrue;
    } else {
      // TODO: Implement, at the very least, yale_is_symmetric. Model it after yale/transp.template.c.
      rb_raise(rb_eNotImpError, "symmetric? and hermitian? only implemented for dense currently");
    }

  }

  return Qfalse;
}


/*
 * Is this matrix symmetric?
 */
static VALUE nm_symmetric(VALUE self) {
  return is_symmetric(self, false);
}

/*
 * Is this matrix hermitian?
 *
 * Definition: http://en.wikipedia.org/wiki/Hermitian_matrix
 *
 * For non-complex matrices, this function should return the same result as symmetric?.
 */
static VALUE nm_hermitian(VALUE self) {
  return is_symmetric(self, true);
}



void Init_nmatrix() {
    /* Require Complex class */
    //rb_require("complex");
    //cComplex = rb_const_get( rb_cObject, rb_intern("Complex") );

    /* Define NMatrix class */
    cNMatrix = rb_define_class("NMatrix", rb_cObject);

    /* class methods */
    rb_define_singleton_method(cNMatrix, "__cblas_gemm__", nm_cblas_gemm, 13);
    rb_define_singleton_method(cNMatrix, "__cblas_gemv__", nm_cblas_gemv, 11);
    rb_define_singleton_method(cNMatrix, "upcast", nm_upcast, 2);

    rb_define_alloc_func(cNMatrix, nm_alloc);
    rb_define_method(cNMatrix, "initialize", nm_init, -1);
    // rb_define_singleton_method(cNMatrix, "new", nm_init, -1);


    rb_define_method(cNMatrix, "initialize_copy", nm_init_copy, 1);
    rb_define_method(cNMatrix, "initialize_cast_copy", nm_init_cast_copy, 2);
    rb_define_method(cNMatrix, "as_dtype", nm_cast_copy, 1);

    /* methods */
    rb_define_method(cNMatrix, "dtype", nm_dtype, 0);
    rb_define_method(cNMatrix, "stype", nm_stype, 0);
    rb_define_method(cNMatrix, "cast", nm_scast_copy, 2);

    rb_define_method(cNMatrix, "[]", nm_mref, -1);
    rb_define_method(cNMatrix, "[]=", nm_mset, -1);
    rb_define_method(cNMatrix, "rank", nm_rank, 0);
    rb_define_alias(cNMatrix, "dim", "rank");
    rb_define_method(cNMatrix, "shape", nm_shape, 0);
    rb_define_method(cNMatrix, "transpose", nm_transpose_new, 0);
    rb_define_method(cNMatrix, "det_exact", nm_det_exact, 0);
    //rb_define_method(cNMatrix, "transpose!", nm_transpose_auto, 0);
    rb_define_method(cNMatrix, "complex_conjugate!", nm_complex_conjugate_bang, 0);

    rb_define_method(cNMatrix, "each", nm_each, 0);

    rb_define_method(cNMatrix, "*", nm_ew_multiply, 1);
    rb_define_method(cNMatrix, "/", nm_ew_divide, 1);
    rb_define_method(cNMatrix, "+", nm_ew_add, 1);
    rb_define_method(cNMatrix, "-", nm_ew_subtract, 1);
    rb_define_method(cNMatrix, "==", nm_ew_eqeq, 1);
    rb_define_method(cNMatrix, "!=", nm_ew_neq, 1);
    rb_define_method(cNMatrix, "<=", nm_ew_leq, 1);
    rb_define_method(cNMatrix, ">=", nm_ew_geq, 1);
    rb_define_method(cNMatrix, "<", nm_ew_lt, 1);
    rb_define_method(cNMatrix, ">", nm_ew_gt, 1);
    rb_define_method(cNMatrix, "eql?", nm_eqeq, 1);
    rb_define_method(cNMatrix, "dot", nm_multiply, 1);
    //rb_define_alias(cNMatrix, "equal?", "eql?");

    rb_define_method(cNMatrix, "symmetric?", nm_symmetric, 0);
    rb_define_method(cNMatrix, "hermitian?", nm_hermitian, 0);


    rb_define_method(cNMatrix, "capacity", nm_capacity, 0);

    rb_define_method(cNMatrix, "__yale_ija__", nm_yale_ija, 0);
    rb_define_method(cNMatrix, "__yale_a__", nm_yale_a, 0);
    rb_define_method(cNMatrix, "__yale_size__", nm_yale_size, 0);
    rb_define_method(cNMatrix, "__yale_ia__", nm_yale_ia, 0);
    rb_define_method(cNMatrix, "__yale_ja__", nm_yale_ja, 0);
    rb_define_method(cNMatrix, "__yale_d__", nm_yale_d, 0);
    rb_define_method(cNMatrix, "__yale_lu__", nm_yale_lu, 0);
    rb_define_const(cNMatrix, "YALE_GROWTH_CONSTANT", rb_float_new(YALE_GROWTH_CONSTANT));


    cNVector = rb_define_class("NVector", cNMatrix);

    // Special exceptions
    nm_eDataTypeError    = rb_define_class("DataTypeError", rb_eStandardError);
    nm_eStorageTypeError = rb_define_class("StorageTypeError", rb_eStandardError);

    nm_id_real  = rb_intern("real");
    nm_id_imag  = rb_intern("imag");
    nm_id_numer = rb_intern("numerator");
    nm_id_denom = rb_intern("denominator");
    nm_id_mult  = rb_intern("*");
    nm_id_add   = rb_intern("+");
    nm_id_multeq= rb_intern("*=");

    nm_id_transpose = rb_intern("transpose");
    nm_id_no_transpose = rb_intern("no_transpose");
    nm_id_complex_conjugate = rb_intern("complex_conjugate");

    nm_id_dense = rb_intern("dense");
    nm_id_list = rb_intern("list");

}

#endif
