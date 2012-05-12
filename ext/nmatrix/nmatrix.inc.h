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
// == nmatrix.h
//

#ifndef NMATRIX_H
#define NMATRIX_H

#include "nmatrix_config.h"

#include <cblas.h>

#include <math.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ruby.h>
#define RUBY_ZERO INT2FIX(0)

#ifdef BENCHMARK // SOURCE: http://stackoverflow.com/questions/2349776/how-can-i-benchmark-a-c-program-easily
# include <sys/time.h>
# include <sys/resource.h>
#endif

#include "dtypes.h"

#include <stddef.h>
#ifdef HAVE_STDINT_H
# include <stdint.h>
#endif

/*
  Data types used in NArray / NMatrix :
  Please modify these types if your system has any different type.
*/


/* NM_BYTE : unsigned 8-bit integer */
#ifndef HAVE_U_INT8_T
# ifdef HAVE_UINT8_T
typedef uint8_t			u_int8_t;
# else
typedef unsigned char		u_int8_t;
# endif
#endif

//#ifndef HAVE_INT8_T
//typedef char                   int8_t;
//#endif

#ifndef HAVE_INT16_T
# if SIZEOF_SHORT == 2
typedef short                  int16_t;
# else
---->> Please define int16_t manually because sizeof(short) != 2. <<----
# endif
#endif /* HAVE_INT16_T */

#ifndef HAVE_INT32_T
# if SIZEOF_LONG == 4
typedef long                   int32_t;
# else
#  if SIZEOF_INT == 4
typedef int                    int32_t;
#  else
---->> Please define int32_t manually because sizeof(long) != 4. <<----
#  endif
# endif
#endif /* HAVE_INT32_T */

/* unsigned 32-bit integer */
#ifndef HAVE_U_INT32_T
# ifdef HAVE_UINT32_T
typedef uint32_t			u_int32_t;
# else
#  if SIZEOF_LONG == 4
typedef unsigned long                   u_int32_t;
#  else
#   if SIZEOF_INT == 4
typedef unsigned int                    u_int32_t;
#   else
---->> Please define u_int32_t manually because sizeof(long) != 4. <<----
#   endif
#  endif
# endif
#endif /* HAVE_U_INT32_T */

#ifndef HAVE_INT64_T
# if SIZEOF_QUAD == 8
typedef quad                   int64_t;
# else
#  if SIZEOF_LONG == 8
typedef long                   int64_t;
#  else
---->> Please define int64_t manually because sizeof(quad) != 8. <<----
#  endif
# endif
#endif /* HAVE_INT64_T */

/* unsigned 64-bit integer */
#ifndef HAVE_U_INT64_T
# ifdef HAVE_UINT64_T
typedef uint64_t            u_int64_t;
# else
#  if SIZEOF_QUAD == 8
typedef unsigned quad       u_int64_t;
#  else
#   if SIZEOF_LONG == 8
typedef unsigned long       u_int64_t;
#   else
---->> Please define u_int64_t manually because sizeof(quad) != 8. <<----
#   endif
#  endif
# endif
#endif /* HAVE_U_INT64_T */


#ifndef HAVE_SIZE_T /// If you modify this, make sure to modify the definition of y_size_t and Y_SIZE_T!
typedef u_int64_t    size_t;
# define NM_SIZE_T   NM_INT64
#else
# if SIZEOF_SIZE_T == 8
#  define NM_SIZE_T  NM_INT64
# else
#  if SIZEOF_SIZE_T == 4
#   define NM_SIZE_T NM_INT32
#  else
---->> Please define size_t and y_size_t manually because sizeof(size_t) is neither 8 nor 4. <<----
#  endif
# endif
#endif

// for when we need to return array indices.
// This must never be larger than size_t
typedef uint32_t    y_size_t;
#define Y_SIZE_T    NM_INT32


#ifdef HAVE_STDBOOL_H
# include <stdbool.h>
#else
typedef char    bool;
# define true    1;
# define false   0;
#endif


typedef struct { float r,i; } complex64;
typedef struct { double r,i; } complex128;
typedef struct { int16_t n,d; } rational32;
typedef struct { int32_t n,d; } rational64;
typedef struct { int64_t n,d; } rational128;


#if SIZEOF_INT == 8
# define DEFAULT_DTYPE  NM_INT64
#else
# if SIZEOF_INT == 4
#  define DEFAULT_DTYPE NM_INT32
# else
#  define DEFAULT_DTYPE NM_INT16
# endif
#endif


#define YALE_GROWTH_CONSTANT    1.5


enum NMatrix_STypes {
  S_DENSE,
  S_LIST,
  S_YALE,
  S_TYPES
};


// Element-wise operations (see blas/elementwise.template.c)
enum NMatrix_Ops {
  NM_OP_ADD = '+',
  NM_OP_SUB = '-',
  NM_OP_MUL = '*',
  NM_OP_DIV = '/',
  NM_OP_MOD = '%',
  NM_OP_BANG = '!',
  NM_OP_NEG, // unary minus
  NM_OP_EQEQ, // ==
  NM_OP_NEQ, // !=
  NM_OP_GT = '>', // >
  NM_OP_LT = '<', // <
  NM_OP_GTE = ',', // >=
  NM_OP_LTE = '.', // <=
  NM_OP_NOT = '~',
  NM_OP_AND = '&',
  NM_OP_OR = '|',
  NM_OP_XOR = '^',
  NM_OP_LSH, // <<
  NM_OP_RSH // >>
};


/* Singly-linked ordered list
 * - holds keys and values
 * - no duplicate keys
 * - keys are ordered
 * - values may be lists themselves
 */
typedef struct l_node { /* Linked list node */
  size_t key;
  void*  val;
  struct l_node * next; // next
} NODE;

typedef struct l_list {
  NODE* first;
} LIST;


// two vectors and a capacity
typedef struct y_vector {
  void*  ija;
  void*  a;
  size_t capacity;
} VECTOR;


typedef struct common_s { // Common elements found in all _s types.
  int8_t    dtype;
  size_t    rank;
  size_t*   shape;
  size_t*   offset;
  void*     elements;
} STORAGE;


typedef struct list_s {
  int8_t    dtype;
  size_t    rank;
  size_t*   shape;
  size_t*   offset;
  void*     default_val;
  LIST*     rows;
} LIST_STORAGE;


typedef struct dense_s {
  int8_t    dtype;
  size_t    rank;
  size_t*   shape;
  size_t*   offset;
  int       count;
  void*     src;
  void*     elements;
} DENSE_STORAGE;


typedef struct yale_s {
  int8_t    dtype;
  size_t    rank;
  size_t*   shape;
  size_t*   offset;
  void*     a;
  size_t    ndnz; // strictly non-diagonal non-zero count!
  size_t    capacity;
  int8_t    index_dtype;
  void*     ija;
} YALE_STORAGE;


typedef struct numeric_matrix {
  int8_t   stype;             /* method of storage (csc, dense, etc) */
  STORAGE* storage;           /* pointer to storage struct */
} NMATRIX;

typedef struct slice_s {
  size_t *coords;           // Coordinate of first element
  size_t *lens;             // Lenght of slice
  uint8_t is_one_el;        // 1 - if all lens eql 1
} SLICE;

/* Local */

typedef union {
  u_int8_t b[2];
  int16_t s;
} nm_size16_t;

typedef union {
  u_int8_t b[4];
  int32_t  i;
  float    f;
} nm_size32_t;

typedef union {
  u_int8_t  b[8];
  int64_t   q;
  float     f[2];
  double    d;
  complex64 c;
} nm_size64_t;

typedef union {
  u_int8_t   b[16];
  int64_t    i[2];
  double     d[2];
  float      f[4];
  complex64  c[2];
  complex128 z;
  rational32 r[4];
  rational64 ra[2];
  rational128 rat;
  VALUE      v[2];
} nm_size128_t;


// For calling cblas_gemm functions (see cblas.c)
typedef struct cblas_param_t {
  int M, N, K, lda, ldb, ldc;
  void *A, *B, *C;
  nm_size128_t alpha, beta;
} DENSE_PARAM;


// Formerly in smmp.h:
typedef struct smmp_param_t {
  void *ia, *ja, *a;
  bool diag;
} YALE_PARAM;

// Shouldn't be necessary, as they're defined in nmatrix.h:
// (Oddly, though, these fix the error.)
/*typedef uint8_t   u_int8_t;
typedef uint16_t  u_int16_t;
typedef uint32_t  u_int32_t;
typedef uint64_t  u_int64_t; */


// rational.c
int64_t nmrb_gcd(int64_t x, int64_t y);

// BLAS functions
#define SMMP_MAX_THREE(a,b,c) ((a)>(b) ? ( (a)>(c) ? (a) : (c) ) : ( (b)>(c) ? (b) : (c) ))
#define SMMP_MIN(a,b) ((a)>(b) ? (b) : (a))
#define SMMP_MAX(a,b) ((a)>(b) ? (a) : (b))

void transp(y_size_t n, y_size_t m, void* ia, void* ja, bool diaga, void* a, void* ib, void* jb, void* b, bool move, int8_t itype, int8_t dtype);


// For binary operations involving matrices that need to be casted.
typedef struct storage_pair_t {
  STORAGE* left;
  STORAGE* right;
} STORAGE_PAIR;


#ifndef NMATRIX_C
extern VALUE cNMatrix;

extern const int nm_sizeof[NM_TYPES+1];
#endif


#define NM_MAX_RANK 15

#define UnwrapNMatrix(obj,var)  Data_Get_Struct(obj, struct numeric_matrix, var)
#define IsNMatrix(obj)          (rb_obj_is_kind_of(obj, CNMatrix)==Qtrue)

#define NM_STORAGE(val)         (((struct numeric_matrix*)DATA_PTR(val))->storage)
//#define NM_PTR(a, p)            ((a)->ptr+(p)*nm_sizeof[(a)->type])
#define NM_STRUCT(val)          ((struct numeric_matrix*)DATA_PTR(val))
//#define NM_PTR_TYPE(val,type)   (type)(((struct numeric_matrix*)DATA_PTR(val))->ptr)
#define NM_RANK(val)            (((STORAGE*)(NM_STORAGE(val)))->rank)
#define NM_DTYPE(val)           (((STORAGE*)(NM_STORAGE(val)))->dtype)
#define NM_STYPE(val)           (((struct numeric_matrix*)DATA_PTR(val))->stype)
#define NM_SHAPE(val,i)         (((STORAGE*)(NM_STORAGE(val)))->shape[(i)])
#define NM_SHAPE0(val)          (((struct numeric_matrix*)DATA_PTR(val))->shape[0])
#define NM_SHAPE1(val)          (((struct numeric_matrix*)DATA_PTR(val))->shape[1])
#define NM_SIZEOF_DTYPE(val)    (nm_sizeof[NM_DTYPE(val)])
#define NM_REF(val,slice)      (RefFuncs[NM_STYPE(val)]( NM_STORAGE(val), slice, NM_SIZEOF_DTYPE(val) ))

#define NM_IsNMatrix(obj) (rb_obj_is_kind_of(obj, cNMatrix)==Qtrue)
#define NM_IsArray(obj)   (TYPE(obj)==T_ARRAY || rb_obj_is_kind_of(obj,cNMatrix)==Qtrue)
#define NM_IsROBJ(d) ((d)->dtype==NM_ROBJ)
#define NM_IsINTEGER(a) \
    (NM_DTYPE(a)==NM_BYTE || NM_DTYPE(a)==NM_INT8 || NM_DTYPE(a)==NM_INT16 || NM_DTYPE(a)==NM_INT32 || NM_DTYPE(a)==NM_INT64)
#define NM_IsCOMPLEX(a) \
    (NM_DTYPE(a)==NM_COMPLEX32 || NM_DTYPE(a)==NM_COMPLEX64)
#define NM_MAX(a,b) (((a)>(b))?(a):(b))
#define NM_SWAP(a,b,tmp) {(tmp)=(a);(a)=(b);(b)=(tmp);}

//#define NUM2REAL(v) NUM2DBL( rb_funcall((v),nm_id_real,0) ) // deprecated
#define REAL2DBL(v) NUM2DBL( rb_funcall((v),nm_id_real,0) )
//#define NUM2IMAG(v) NUM2DBL( rb_funcall((v),nm_id_imag,0) ) // deprecated
#define IMAG2DBL(v) NUM2DBL( rb_funcall((v),nm_id_imag,0) )

#define NUM2NUMER(v) NUM2INT( rb_funcall((v), nm_id_numer,0) ) // deprecated
#define NUMER2INT(v) NUM2INT( rb_funcall((v), nm_id_numer,0) )
#define NUM2DENOM(v) NUM2INT( rb_funcall((v), nm_id_denom,0) ) // deprecated
#define DENOM2INT(v) NUM2INT( rb_funcall((v), nm_id_denom,0) )

#define IS_NUMERIC(v)   (FIXNUM_P(v) || TYPE(v) == T_FLOAT || TYPE(v) == T_COMPLEX || TYPE(v) == T_RATIONAL)
#define IS_STRING(v)    (TYPE(v) == T_STRING)

#define CheckNMatrixType(v)   if (TYPE(v) != T_DATA || (RDATA(v)->dfree != (RUBY_DATA_FUNC)nm_delete && RDATA(v)->dfree != (RUBY_DATA_FUNC)nm_delete_ref)) rb_raise(rb_eTypeError, "expected NMatrix on left-hand side of operation");

//#define YALE_JA_START(sptr)             (((YALE_STORAGE*)(sptr))->shape[0]+1)
#define YALE_IJA(sptr,elem_size,i)          (void*)( (char*)(((YALE_STORAGE*)(sptr))->ija) + i * elem_size )
//#define YALE_JA(sptr,dtype,j)           ((((dtype)*)((YALE_STORAGE*)(sptr))->ija)[(YALE_JA_START(sptr))+j])
#define YALE_ROW_LENGTH(sptr,elem_size,i)   (*(size_t*)YALE_IA((sptr),(elem_size),(i)+1) - *(size_t*)YALE_IJA((sptr),(elem_size),(i)))
#define YALE_A(sptr,elem_size,i)            (void*)((char*)(((YALE_STORAGE*)(sptr))->a) + elem_size * i)
#define YALE_DIAG(sptr, elem_size, i)       ( YALE_A((sptr),(elem_size),(i)) )
//#define YALE_LU(sptr,dtype,i,j)             (((dtype)*)(((YALE_STORAGE*)(sptr))->a)[ YALE_JA_START(sptr) +  ])
#define YALE_MINIMUM(sptr)                  (((YALE_STORAGE*)(sptr))->shape[0]*2 + 1) // arbitrarily defined
#define YALE_SIZE_PTR(sptr,elem_size)       (void*)((char*)((YALE_STORAGE*)(sptr))->ija + ((YALE_STORAGE*)(sptr))->shape[0]*elem_size )
#define YALE_MAX_SIZE(sptr)                 (((YALE_STORAGE*)(sptr))->shape[0] * ((YALE_STORAGE*)(sptr))->shape[1] + 1)
#define YALE_IA_SIZE(sptr)                  ((YALE_STORAGE*)(sptr))->shape[0]

// None of these next three return anything. They set a reference directly.
#define YaleGetIJA(victim,s,i)              (SetFuncs[Y_SIZE_T][(s)->index_dtype](1, &(victim), 0, YALE_IJA((s), nm_sizeof[s->index_dtype], (i)), 0))
#define YaleSetIJA(i,s,from)                (SetFuncs[s->index_dtype][Y_SIZE_T](1, YALE_IJA((s), nm_sizeof[s->index_dtype], (i)), 0, &(from), 0))
#define YaleGetSize(sz,s)                   (SetFuncs[Y_SIZE_T][((YALE_STORAGE*)s)->index_dtype](1, &sz, 0, (YALE_SIZE_PTR(((YALE_STORAGE*)s), nm_sizeof[((YALE_STORAGE*)s)->index_dtype])), 0))
//#define YALE_FIRST_NZ_ROW_ENTRY(sptr,elem_size,i)


#if !defined RSTRING_LEN
#define RSTRING_LEN(a) RSTRING(a)->len
#endif
#if !defined RSTRING_PTR
#define RSTRING_PTR(a) RSTRING(a)->ptr
#endif
#if !defined RARRAY_LEN
#define RARRAY_LEN(a) RARRAY(a)->len
#endif
#if !defined RARRAY_PTR
#define RARRAY_PTR(a) RARRAY(a)->ptr
#endif

#define NM_INDEX_TYPES  NM_FLOAT32


typedef void     (*nm_setfunc_t[NM_TYPES][NM_TYPES])(); // copy functions
typedef void     (*nm_incfunc_t[NM_TYPES])();           // increment functions
typedef void*    (*nm_stype_ref_t[S_TYPES])(STORAGE*, SLICE*);        // get/ref
typedef VALUE    (*nm_stype_ins_t[S_TYPES])(STORAGE*, SLICE*, VALUE); // insert
typedef STORAGE* (*nm_create_storage_t[S_TYPES])();
typedef STORAGE* (*nm_cast_copy_storage_t[S_TYPES])();
typedef STORAGE* (*nm_scast_copy_storage_t[S_TYPES][S_TYPES])();
typedef NMATRIX* (*nm_matrix_multiply_op_t[S_TYPES])();
typedef NMATRIX* (*nm_elementwise_binary_op_casted_t[S_TYPES])();
typedef int      (*nm_d_elementwise_binary_op_t[NM_TYPES])();
typedef int      (*nm_y_elementwise_binary_op_t[NM_TYPES][NM_INDEX_TYPES])();
typedef bool     (*nm_compare_t[S_TYPES])();
typedef bool     (*nm_eqeq_t[NM_TYPES][2])(const void*, const void*, const size_t, const size_t);
typedef void     (*nm_delete_t[S_TYPES])();
typedef void     (*nm_mark_t[S_TYPES])(void*);
typedef void     (*nm_gemm_t[NM_TYPES])();           // general matrix/matrix multiply
typedef void     (*nm_det_t[NM_TYPES])(const int, const void*, const int, void*);            // determinant
typedef NMATRIX* (*nm_transpose_t[S_TYPES])();
typedef void     (*nm_dense_transpose_t[NM_TYPES])(); // dense transpose
typedef void     (*nm_gemv_t[NM_TYPES])();           // general matrix/vector multiply
typedef void     (*nm_smmp_t[NM_TYPES][NM_INDEX_TYPES])(); // sparse (yale) multiply
typedef void     (*nm_smmp_transpose_t[NM_TYPES][NM_INDEX_TYPES])(y_size_t, y_size_t, YALE_PARAM, YALE_PARAM, bool); // sparse (yale) transpose
//typedef void (*nm_setsf_t[S_TYPES][S_TYPES])();
//typedef void (*nm_setdf_t[NM_DTYPES][NM_DTYPES])();

extern nm_setfunc_t SetFuncs;
extern nm_incfunc_t Increment;
extern ID nm_id_real, nm_id_imag;
extern ID nm_id_denom, nm_id_numer;
extern ID nm_id_mult, nm_id_multeq, nm_id_add;

/**********************************************************************************************************************/
/* DO NOT MODIFY BETWEEN THESE LINES                                                                                  */
// %%INSERT_TEMPLATER_DECLARATIONS%%
/* DO NOT MODIFY BETWEEN THESE LINES                                                                                  */
/**********************************************************************************************************************/

/* cblas.c */
DENSE_PARAM init_cblas_params_for_nm_multiply_matrix(int8_t dtype);
void cblas_bgemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_bgemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_i8gemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_i8gemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_i16gemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_i16gemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_i32gemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_i32gemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_i64gemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_i64gemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_sgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_sgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_dgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_dgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_cgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_cgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_zgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_zgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_r32gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_r32gemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_r64gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_r64gemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_r128gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_r128gemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_vgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);


/* dense.c */
DENSE_STORAGE*  create_dense_storage(int8_t dtype, size_t* shape, size_t rank, void* elements, size_t elements_length);
void            delete_dense_storage_ref(DENSE_STORAGE* s);
void            delete_dense_storage(DENSE_STORAGE* s);
void            mark_dense_storage(void* s);
DENSE_STORAGE*  cast_copy_dense_storage(DENSE_STORAGE* rhs, int8_t new_dtype);

size_t          count_dense_storage_elements(const DENSE_STORAGE* s);
bool            dense_storage_eqeq(const DENSE_STORAGE*, const DENSE_STORAGE*);
bool            dense_is_symmetric(const DENSE_STORAGE*, int, bool);

size_t          dense_storage_pos(DENSE_STORAGE* s, SLICE* slice);
void*           dense_storage_get(DENSE_STORAGE* s, SLICE* slice);
void            dense_storage_set(DENSE_STORAGE* s, SLICE* slice, void* val);

/* list.c */
LIST_STORAGE*   create_list_storage(int8_t dtype, size_t* shape, size_t rank, void* init_val);
void            delete_list_storage(LIST_STORAGE* s);
void            mark_list_storage(void* s);
LIST_STORAGE*   cast_copy_list_storage(LIST_STORAGE* rhs, int8_t new_dtype);
size_t          count_storage_max_elements(const STORAGE*);

void*           list_storage_get(LIST_STORAGE* s, SLICE* slice);
void*           list_storage_insert(LIST_STORAGE* s, SLICE* slice, void* val);
void*           list_storage_remove(LIST_STORAGE* s, SLICE* slice);
bool            list_storage_eqeq(const LIST_STORAGE*, const LIST_STORAGE*);

/* yale.c */
void print_vectors(YALE_STORAGE* s);
YALE_STORAGE*   create_yale_storage(int8_t dtype, size_t* shape, size_t rank, size_t init_capacity);
YALE_STORAGE*   create_yale_storage_from_old_yale(int8_t dtype, size_t* shape, char* ia, char* ja, char* a, int8_t from_dtype, int8_t from_index_dtype);
void            init_yale_storage(YALE_STORAGE* s);
void            delete_yale_storage(YALE_STORAGE* s);
void            mark_yale_storage(void* s);
YALE_STORAGE*   cast_copy_yale_storage(YALE_STORAGE* rhs, int8_t new_dtype);
bool            yale_storage_eqeq(const YALE_STORAGE*, const YALE_STORAGE*);

void*           yale_storage_ref(YALE_STORAGE* s, SLICE* slice);
char            yale_storage_set(YALE_STORAGE* s, SLICE* slice, void* v);

YALE_STORAGE*   create_merged_yale_storage(const YALE_STORAGE*, const YALE_STORAGE*);

size_t          count_list_storage_nd_elements(const LIST_STORAGE*);
size_t          count_list_storage_elements(const LIST_STORAGE*);


/* stype casts */
DENSE_STORAGE* scast_copy_dense_yale(const YALE_STORAGE* rhs, int8_t l_dtype);
DENSE_STORAGE* scast_copy_dense_list(const LIST_STORAGE* rhs, int8_t l_dtype);
YALE_STORAGE* scast_copy_yale_dense(const DENSE_STORAGE* rhs, int8_t l_dtype);
YALE_STORAGE* scast_copy_yale_list(const LIST_STORAGE* rhs, int8_t l_dtype);
LIST_STORAGE* scast_copy_list_yale(const YALE_STORAGE* rhs, int8_t l_dtype);
LIST_STORAGE* scast_copy_list_dense(const DENSE_STORAGE* rhs, int8_t l_dtype);

/* nmatrix.c */
void cast_copy_value_single(void* to, const void* from, int8_t l_dtype, int8_t r_dtype);
int8_t nm_dtypestring_to_dtype(VALUE str);
int8_t nm_dtypesymbol_to_dtype(VALUE sym);
int8_t nm_stypestring_to_stype(VALUE str);
int8_t nm_stypesymbol_to_stype(VALUE sym);
int8_t nm_guess_dtype(VALUE v);
size_t* nm_interpret_shape_arg(VALUE arg, size_t* rank);
NMATRIX* nm_create(int8_t stype, void* storage);
void Init_nmatrix();

#endif
