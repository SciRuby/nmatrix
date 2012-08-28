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
// C and C++ API for NMatrix, and main header file.

#ifndef NMATRIX_H
#define NMATRIX_H

/*
 * Standard Includes
 */

#include <ruby.h>

#ifdef __cplusplus
  #include <cmath>
  #include <cstring>
#else
  #include <math.h>
  #include <string.h>
#endif



#ifdef BENCHMARK
	// SOURCE: http://stackoverflow.com/questions/2349776/how-can-i-benchmark-a-c-program-easily
  #ifdef __cplusplus
    #include <sys/ctime>
    #include <sys/cresource>
  #else
    #include <sys/time.h>
    #include <sys/resource.h>
	#endif
#endif

/*
 * Project Includes
 */

/*
 * Macros
 */


#define RUBY_ZERO INT2FIX(0)

#ifndef SIZEOF_INT
 #error SIZEOF_INT undefined
#else
 #if SIZEOF_INT == 8
	#define DEFAULT_DTYPE  INT64
	#define SIZE_T         INT64
 #else
	#if SIZEOF_INT == 4
		#define DEFAULT_DTYPE INT32
		#define SIZE_T        INT32
	#else
	 #if SIZEOF_INT == 2
		#define DEFAULT_DTYPE INT16
		#define SIZE_T        INT16
	 #else
	  #error Unhandled SIZEOF_INT -- please #define SIZE_T and DEFAULT_DTYPE manually.
	 #endif
	#endif
 #endif
#endif

// The definition above is partially overridden by the definition below when size_t is undefined.
#ifndef HAVE_SIZE_T
 typedef u_int64_t  size_t;
 #define SIZE_T     INT64
#endif


/*
 * == Macros for Concurrent C and C++ Header Maintenance
 *
 * These macros look complicated, but they're really not so bad. They're also important: they ensure that whether our
 * header file (nmatrix.h) is read by a C++ or a C compiler, all the same data structures and enumerators exist, albeit
 * with slightly different names.
 *
 * "But wait," you say, "You use structs. Structs exist in C and C++. Why use a macro to set them up?"
 *
 * Well, in C, you have to be explicit about what a struct is. You can actually get around that requirement by using a
 * typedef:
 *
 *   typedef struct STORAGE { ... } STORAGE;
 *
 * Also, we use C++ inheritance, which is obviously not allowed in C. So we have to ensure that the base class's members
 * are exposed properly to our child classes.
 *
 * The macros also allow us to put all of our C++ types into namespaces. For C, we prefix everything with either nm_ or
 * NM_ to distinguish our declarations from those in other libraries.
 */


#ifdef __cplusplus /* These are the C++ versions of the macros. */

  #define NM_DECL_ENUM(enum_type, name)   enum_type name;
  #define NM_DECL_STRUCT(type, name)      type      name;

  #define NM_DEF_STORAGE_ELEMENTS    \
    NM_DECL_ENUM(dtype_t, dtype);    \
    size_t  dim;                     \
    size_t* shape;                   \
    size_t* offset;                  \
	  int			count;                   \
	  STORAGE*		src;

  #define NM_DEF_STORAGE_CHILD_STRUCT_PRE(name)    struct name : STORAGE {
  #define NM_DEF_STORAGE_STRUCT_POST(name)         };

  #define NM_DEF_STORAGE_STRUCT      \
  struct STORAGE {                   \
    NM_DEF_STORAGE_ELEMENTS;         \
  };

  #define NM_DEF_STRUCT_PRE(name)  struct name {
  #define NM_DEF_STRUCT_POST(name) };

  #define NM_DEF_ENUM(name, ...)         \
    enum name {                           \
      __VA_ARGS__                         \
    };

#else   /* These are the C versions of the macros. */

  #define NM_DECL_ENUM(enum_type, name)   nm_ ## enum_type name;
  #define NM_DECL_STRUCT(type, name)      NM_ ## type      name;

  #define NM_DEF_STORAGE_ELEMENTS   \
    NM_DECL_ENUM(dtype_t, dtype);   \
    size_t      dim;                \
    size_t*     shape;              \
    size_t*     offset;             \
	  int			count;                  \ 
	  STORAGE*		src;              
  #define NM_DEF_STORAGE_CHILD_STRUCT_PRE(name)  typedef struct NM_ ## name { \
                                                    NM_DEF_STORAGE_ELEMENTS;

  #define NM_DEF_STORAGE_STRUCT_POST(name)       } NM_ ## name;

  #define NM_DEF_STORAGE_STRUCT      \
  typedef struct NM_STORAGE {        \
    NM_DEF_STORAGE_ELEMENTS;         \
  } NM_STORAGE;
>>>>>>> compile_point

  #define NM_DEF_STRUCT_PRE(name)                typedef struct NM_ ## name {
  #define NM_DEF_STRUCT_POST(name)               } NM_ ## name;

  #define NM_DEF_ENUM(name, ...)     \
    typedef enum nm_ ## name {       \
      __VA_ARGS__                    \
    } nm_ ## name;

#endif      /* End of C/C++ Parallel Header Macro Definitions */


/*
 * Types
 */

#define NM_NUM_DTYPES 13  // data/data.h
#define NM_NUM_ITYPES 4   // data/data.h
#define NM_NUM_STYPES 3   // storage/storage.h

#ifndef __cplusplus
//namespace nm {
#else

#endif

/* Storage Type -- Dense or Sparse */
NM_DEF_ENUM(stype_t,  DENSE_STORE = 0,
                      LIST_STORE = 1,
                      YALE_STORE = 2);

/* Data Type */
NM_DEF_ENUM(dtype_t,	BYTE				=  0,  // unsigned char
                    	INT8				=  1,  // char
                    	INT16				=  2,  // short
                    	INT32				=  3,  // int
                    	INT64				=  4,  // long
                    	FLOAT32			=  5,  // float
                    	FLOAT64			=  6,  // double
                    	COMPLEX64		=  7,  // Complex64 class
                    	COMPLEX128	=  8,  // Complex128 class
                    	RATIONAL32	=  9,  // Rational32 class
                    	RATIONAL64	= 10,  // Rational64 class
                    	RATIONAL128	= 11,  // Rational128 class
                    	RUBYOBJ			= 12);  // Ruby VALUE type

/* Index Type for Yale Matrices */
NM_DEF_ENUM(itype_t,  UINT8  = 0,
                      UINT16 = 1,
                      UINT32 = 2,
                      UINT64 = 3);

/* struct STORAGE */
NM_DEF_STORAGE_STRUCT;

/* Dense Storage */
NM_DEF_STORAGE_CHILD_STRUCT_PRE(DENSE_STORAGE); // struct DENSE_STORAGE : STORAGE {
	size_t*	stride;
	void*		elements;
NM_DEF_STORAGE_STRUCT_POST(DENSE_STORAGE);     // };

/* Yale Storage */
NM_DEF_STORAGE_CHILD_STRUCT_PRE(YALE_STORAGE);
	void* a;      // should go first
	size_t ndnz; // Strictly non-diagonal non-zero count!
	size_t	capacity;
	itype_t	itype;
	void*		ija;
NM_DEF_STORAGE_STRUCT_POST(YALE_STORAGE);

// FIXME: NODE and LIST should be put in some kind of namespace or something, at least in C++.
NM_DEF_STRUCT_PRE(NODE); // struct NODE {
  size_t key;
  void*  val;
  NM_DECL_STRUCT(NODE*, next);  // NODE* next;
NM_DEF_STRUCT_POST(NODE); // };

NM_DEF_STRUCT_PRE(LIST); // struct LIST {
  NM_DECL_STRUCT(NODE*, first); // NODE* first;
NM_DEF_STRUCT_POST(LIST); // };

/* List-of-Lists Storage */
NM_DEF_STORAGE_CHILD_STRUCT_PRE(LIST_STORAGE); // struct LIST_STORAGE : STORAGE {
	// List storage specific elements.
	void* default_val;
	NM_DECL_STRUCT(LIST*, rows); // LIST* rows;
NM_DEF_STORAGE_STRUCT_POST(LIST_STORAGE);      // };



/* NMATRIX Object */
NM_DEF_STRUCT_PRE(NMATRIX);   // struct NMATRIX {
  NM_DECL_ENUM(stype_t, stype);       // stype_t stype;     // Method of storage (csc, dense, etc).
  NM_DECL_STRUCT(STORAGE*, storage);  // STORAGE* storage;  // Pointer to storage struct.
NM_DEF_STRUCT_POST(NMATRIX);  // };

#define NM_MAX_RANK 15

#define UnwrapNMatrix(obj,var)  Data_Get_Struct(obj, NMATRIX, var)

#define NM_STORAGE(val)         (NM_STRUCT(val)->storage)
#ifdef __cplusplus
  #define NM_STRUCT(val)              ((NMATRIX*)(DATA_PTR(val)))
  #define NM_STORAGE_LIST(val)        ((LIST_STORAGE*)(NM_STORAGE(val)))
  #define NM_STORAGE_YALE(val)        ((YALE_STORAGE*)(NM_STORAGE(val)))
  #define NM_STORAGE_DENSE(val)       ((DENSE_STORAGE*)(NM_STORAGE(val)))
#else
  #define NM_STRUCT(val)              ((struct NM_NMATRIX*)(DATA_PTR(val)))
  #define NM_STORAGE_LIST(val)        ((struct NM_LIST_STORAGE*)(NM_STORAGE(val)))
  #define NM_STORAGE_YALE(val)        ((struct NM_YALE_STORAGE*)(NM_STORAGE(val)))
  #define NM_STORAGE_DENSE(val)       ((struct NM_DENSE_STORAGE*)(NM_STORAGE(val)))
#endif

#define NM_DENSE_SRC(val)       (NM_STORAGE_DENSE(val)->src)
#define NM_LIST_SRC(val)        (NM_STORAGE_LIST(val)->src)
#define NM_DIM(val)             (NM_STORAGE(val)->dim)
#define NM_DTYPE(val)           (NM_STORAGE(val)->dtype)
#define NM_ITYPE(val)           (NM_STORAGE_YALE(val)->itype)
#define NM_STYPE(val)           (NM_STRUCT(val)->stype)
#define NM_SHAPE(val,i)         (NM_STORAGE(val)->shape[(i)])
#define NM_SHAPE0(val)          (NM_STORAGE(val)->shape[0])
#define NM_SHAPE1(val)          (NM_STORAGE(val)->shape[1])

#define NM_DENSE_COUNT(val)     (storage_count_max_elements(NM_STORAGE_DENSE(val)))
#define NM_SIZEOF_DTYPE(val)    (DTYPE_SIZES[NM_DTYPE(val)])
#define NM_REF(val,slice)      (RefFuncs[NM_STYPE(val)]( NM_STORAGE(val), slice, NM_SIZEOF_DTYPE(val) ))
    
#define NM_MAX(a,b) (((a)>(b))?(a):(b))
#define NM_MIN(a,b) (((a)>(b))?(b):(a))
#define NM_SWAP(a,b,tmp) {(tmp)=(a);(a)=(b);(b)=(tmp);}

#define NM_CHECK_ALLOC(x) if (!x) rb_raise(rb_eNoMemError, "insufficient memory");

#define CheckNMatrixType(v)   if (TYPE(v) != T_DATA || (RDATA(v)->dfree != (RUBY_DATA_FUNC)nm_delete && RDATA(v)->dfree != (RUBY_DATA_FUNC)nm_delete_ref)) rb_raise(rb_eTypeError, "expected NMatrix on left-hand side of operation");

#define NM_IsNMatrix(obj) \
  (rb_obj_is_kind_of(obj, cNMatrix) == Qtrue)

#define NM_IsNVector(obj) \
  (rb_obj_is_kind_of(obj, cNVector) == Qtrue)


typedef VALUE (*METHOD)(...);

#ifdef __cplusplus
//}; // end of namespace nm
#endif

/*
 * Data
 */

/*
 * Functions
 */

extern "C" {
	void Init_nmatrix();
	
	// External API
	VALUE rb_nmatrix_dense_create(dtype_t dtype, size_t* shape, size_t dim, void* elements, size_t length);
	VALUE rb_nvector_dense_create(dtype_t dtype, void* elements, size_t length);
}

#endif // NMATRIX_H
