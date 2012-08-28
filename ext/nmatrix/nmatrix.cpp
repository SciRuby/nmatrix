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
// == nmatrix.cpp
//
// Main C++ source file for NMatrix. Contains Init_nmatrix and most Ruby instance and
// class methods for NMatrix. Also responsible for calling Init methods on related
// modules.

/*
 * Standard Includes
 */

#include <ruby.h>
#include <algorithm> // std::min

/*
 * Project Includes
 */
#include "types.h"
#include "data/data.h"
#include "util/math.h"
#include "util/io.h"
#include "storage/storage.h"

#include "nmatrix.h"

#include "ruby_constants.h"

/*
 * Macros
 */

/*
 * Global Variables
 */


extern "C" {

/*
 * Forward Declarations
 */

static VALUE nm_init(int argc, VALUE* argv, VALUE nm);
static VALUE nm_init_copy(VALUE copy, VALUE original);
static VALUE nm_init_transposed(VALUE self);
static VALUE nm_init_cast_copy(VALUE self, VALUE new_stype_symbol, VALUE new_dtype_symbol);
static VALUE nm_init_yale_from_old_yale(VALUE shape, VALUE dtype, VALUE ia, VALUE ja, VALUE a, VALUE from_dtype, VALUE nm);
static VALUE nm_alloc(VALUE klass);
static void  nm_delete(NMATRIX* mat);
static void  nm_delete_ref(NMATRIX* mat);
static VALUE nm_dtype(VALUE self);
static VALUE nm_itype(VALUE self);
static VALUE nm_stype(VALUE self);
static VALUE nm_dim(VALUE self);
static VALUE nm_shape(VALUE self);
static VALUE nm_capacity(VALUE self);
static VALUE nm_each(VALUE nmatrix);

static SLICE* get_slice(size_t dim, VALUE* c, VALUE self);
static VALUE nm_xslice(int argc, VALUE* argv, void* (*slice_func)(STORAGE*, SLICE*), void (*delete_func)(NMATRIX*), VALUE self);
static VALUE nm_mset(int argc, VALUE* argv, VALUE self);
static VALUE nm_mget(int argc, VALUE* argv, VALUE self);
static VALUE nm_mref(int argc, VALUE* argv, VALUE self);
static VALUE nm_is_ref(VALUE self);

static VALUE is_symmetric(VALUE self, bool hermitian);

static VALUE nm_ew_add(VALUE left_val, VALUE right_val);
static VALUE nm_ew_subtract(VALUE left_val, VALUE right_val);
static VALUE nm_ew_multiply(VALUE left_val, VALUE right_val);
static VALUE nm_ew_divide(VALUE left_val, VALUE right_val);
//static VALUE nm_ew_mod(VALUE left_val, VALUE right_val);

static VALUE elementwise_op(nm::ewop_t op, VALUE left_val, VALUE right_val);

static VALUE nm_symmetric(VALUE self);
static VALUE nm_hermitian(VALUE self);

static VALUE nm_eqeq(VALUE left, VALUE right);

static VALUE matrix_multiply_scalar(NMATRIX* left, VALUE scalar);
static VALUE matrix_multiply(NMATRIX* left, NMATRIX* right);
static VALUE nm_multiply(VALUE left_v, VALUE right_v);
static VALUE nm_factorize_lu_bang(VALUE self);
static VALUE nm_det_exact(VALUE self);
static VALUE nm_complex_conjugate_bang(VALUE self);

static dtype_t	dtype_guess(VALUE v);
static dtype_t	interpret_dtype(int argc, VALUE* argv, stype_t stype);
static void*		interpret_initial_value(VALUE arg, dtype_t dtype);
static size_t*	interpret_shape(VALUE arg, size_t* dim);
static stype_t	interpret_stype(VALUE arg);

/* Singleton methods */
static VALUE nm_itype_by_shape(VALUE self, VALUE shape_arg);
static VALUE nm_upcast(VALUE self, VALUE t1, VALUE t2);


#ifdef BENCHMARK
static double get_time(void);
#endif

/*
 * Functions
 */

///////////////////
// Ruby Bindings //
///////////////////

void Init_nmatrix() {
	
	///////////////////////
	// Class Definitions //
	///////////////////////
	
	cNMatrix = rb_define_class("NMatrix", rb_cObject);
	cNVector = rb_define_class("NVector", cNMatrix);
	
	// Special exceptions
	nm_eDataTypeError    = rb_define_class("DataTypeError",			rb_eStandardError);
	nm_eStorageTypeError = rb_define_class("StorageTypeError",	rb_eStandardError);

	///////////////////
	// Class Methods //
	///////////////////
	
	rb_define_alloc_func(cNMatrix, nm_alloc);

	///////////////////////
  // Singleton Methods //
  ///////////////////////

	rb_define_singleton_method(cNMatrix, "upcast", (METHOD)nm_upcast, 2);
	rb_define_singleton_method(cNMatrix, "itype_by_shape", (METHOD)nm_itype_by_shape, 1);
	
	//////////////////////
	// Instance Methods //
	//////////////////////

	rb_define_method(cNMatrix, "initialize", (METHOD)nm_init, -1);
	rb_define_method(cNMatrix, "initialize_copy", (METHOD)nm_init_copy, 1);

	// Technically, the following function is a copy constructor.
	rb_define_method(cNMatrix, "transpose", (METHOD)nm_init_transposed, 0);

	rb_define_method(cNMatrix, "dtype", (METHOD)nm_dtype, 0);
	rb_define_method(cNMatrix, "itype", (METHOD)nm_itype, 0);
	rb_define_method(cNMatrix, "stype", (METHOD)nm_stype, 0);
	rb_define_method(cNMatrix, "cast",  (METHOD)nm_init_cast_copy, 2);

	rb_define_method(cNMatrix, "[]", (METHOD)nm_mref, -1);
	rb_define_method(cNMatrix, "slice", (METHOD)nm_mget, -1);
	rb_define_method(cNMatrix, "[]=", (METHOD)nm_mset, -1);
	rb_define_method(cNMatrix, "is_ref?", (METHOD)nm_is_ref, 0);
	rb_define_method(cNMatrix, "dimensions", (METHOD)nm_dim, 0);

	rb_define_method(cNMatrix, "shape", (METHOD)nm_shape, 0);
	rb_define_method(cNMatrix, "det_exact", (METHOD)nm_det_exact, 0);
	//rb_define_method(cNMatrix, "transpose!", (METHOD)nm_transpose_self, 0);
	rb_define_method(cNMatrix, "complex_conjugate!", (METHOD)nm_complex_conjugate_bang, 0);

	rb_define_method(cNMatrix, "each", (METHOD)nm_each, 0);

	rb_define_method(cNMatrix, "+",			(METHOD)nm_ew_add,			1);
	rb_define_method(cNMatrix, "-",			(METHOD)nm_ew_subtract,	1);
  rb_define_method(cNMatrix, "*",			(METHOD)nm_ew_multiply,	1);
	rb_define_method(cNMatrix, "/",			(METHOD)nm_ew_divide,		1);
//rb_define_method(cNMatrix, "%",			(METHOD)nm_ew_mod,			1);
	rb_define_method(cNMatrix, "eql?",	(METHOD)nm_eqeq,				1);

	/////////////////////////
	// Matrix Math Methods //
	/////////////////////////
	rb_define_method(cNMatrix, "dot",		(METHOD)nm_multiply,		1);
	rb_define_method(cNMatrix, "factorize_lu!", (METHOD)nm_factorize_lu_bang, 0);
	
	/*
	 * TODO: Write new elementwise code for boolean operations
	rb_define_method(cNMatrix, "==", (METHOD)nm_ew_eqeq, 1);
	rb_define_method(cNMatrix, "!=", (METHOD)nm_ew_neq, 1);
	rb_define_method(cNMatrix, "<=", (METHOD)nm_ew_leq, 1);
	rb_define_method(cNMatrix, ">=", (METHOD)nm_ew_geq, 1);
	rb_define_method(cNMatrix, "<", (METHOD)nm_ew_lt, 1);
	rb_define_method(cNMatrix, ">", (METHOD)nm_ew_gt, 1);
	 */

	rb_define_method(cNMatrix, "symmetric?", (METHOD)nm_symmetric, 0);
	rb_define_method(cNMatrix, "hermitian?", (METHOD)nm_hermitian, 0);

	rb_define_method(cNMatrix, "capacity", (METHOD)nm_capacity, 0);
	
	/////////////
	// Aliases //
	/////////////
	
	rb_define_alias(cNMatrix, "dim", "dimensions");
	rb_define_alias(cNMatrix, "equal?", "eql?");
	
	///////////////////////
	// Symbol Generation //
	///////////////////////
	
	nm_init_ruby_constants();

	//////////////////////////
	// YaleFunctions module //
	//////////////////////////

	nm_init_yale_functions();

	/////////////////
	// BLAS module //
	/////////////////

	nm_math_init_blas();

	///////////////
	// IO module //
	///////////////
	nm_init_io();
}


//////////////////
// Ruby Methods //
//////////////////

/*
 * Allocator.
 */
static VALUE nm_alloc(VALUE klass) {
  NMATRIX* mat = ALLOC(NMATRIX);
  mat->storage = NULL;
  // FIXME: mark_table[mat->stype] should be passed to Data_Wrap_Struct, but can't be done without stype. Also, nm_delete depends on this.
  // mat->stype   = nm::NUM_STYPES;

  //STYPE_MARK_TABLE(mark_table);

  return Data_Wrap_Struct(klass, NULL, nm_delete, mat);
}

/*
 * Find the capacity of an NMatrix. The capacity only differs from the size for
 * Yale matrices, which occasionally allocate more space than they need. For
 * list and dense, capacity gives the number of elements in the matrix.
 */
static VALUE nm_capacity(VALUE self) {
  VALUE cap;

  switch(NM_STYPE(self)) {
  case YALE_STORE:
    cap = UINT2NUM(((YALE_STORAGE*)(NM_STORAGE(self)))->capacity);
    break;

  case DENSE_STORE:
    cap = UINT2NUM(nm_storage_count_max_elements( NM_STORAGE_DENSE(self) ));
    break;

  case LIST_STORE:
    cap = UINT2NUM(nm_list_storage_count_elements( NM_STORAGE_LIST(self) ));
    break;

  default:
    rb_raise(nm_eStorageTypeError, "unrecognized stype in nm_capacity()");
  }

  return cap;
}

/*
 * Destructor.
 */
static void nm_delete(NMATRIX* mat) {
  static void (*ttable[nm::NUM_STYPES])(STORAGE*) = {
    nm_dense_storage_delete,
    nm_list_storage_delete,
    nm_yale_storage_delete
  };
  ttable[mat->stype](mat->storage);
}

/*
 * Slicing destructor.
 */
static void nm_delete_ref(NMATRIX* mat) {
  static void (*ttable[nm::NUM_STYPES])(STORAGE*) = {
    nm_dense_storage_delete_ref,
    nm_list_storage_delete_ref, 
    nm_yale_storage_delete
  };
  ttable[mat->stype](mat->storage);
}

/*
 * Get the data type (dtype) of a matrix, e.g., :byte, :int8, :int16, :int32,
 * :int64, :float32, :float64, :complex64, :complex128, :rational32,
 * :rational64, :rational128, or :object (the last is a Ruby object).
 */
static VALUE nm_dtype(VALUE self) {
  ID dtype = rb_intern(DTYPE_NAMES[NM_DTYPE(self)]);
  return ID2SYM(dtype);
}


/*
 * Get the index data type (dtype) of a matrix. Defined only for yale; others return nil.
 */
static VALUE nm_itype(VALUE self) {
  if (NM_STYPE(self) == YALE_STORE) {
    ID itype = rb_intern(ITYPE_NAMES[NM_ITYPE(self)]);
    return ID2SYM(itype);
  }
  return Qnil;
}


/*
 * Get the index data type (dtype) of a matrix. Defined only for yale; others return nil.
 */
static VALUE nm_itype_by_shape(VALUE self, VALUE shape_arg) {

  size_t dim;
  size_t* shape = interpret_shape(shape_arg, &dim);

  itype_t itype = nm_yale_storage_itype_by_shape(shape);
  ID itype_id   = rb_intern(ITYPE_NAMES[itype]);

  return ID2SYM(itype_id);
}


/*
 * Given a binary operation between types t1 and t2, what type will be returned?
 *
 * This is a singleton method on NMatrix, e.g., NMatrix.upcast(:int32, :int64)
 */
static VALUE nm_upcast(VALUE self, VALUE t1, VALUE t2) {

  dtype_t d1    = nm_dtype_from_rbsymbol(t1),
          d2    = nm_dtype_from_rbsymbol(t2);

  return ID2SYM(rb_intern( DTYPE_NAMES[ Upcast[d1][d2] ] ));
}



/*
 * Borrowed this function from NArray. Handles 'each' iteration on a dense
 * matrix.
 *
 * Additionally, handles separately matrices containing VALUEs and matrices
 * containing other types of data.
 */
static VALUE nm_each_dense(VALUE nmatrix) {
  DENSE_STORAGE* s = NM_STORAGE_DENSE(nmatrix);
  VALUE v;
  size_t i;

  if (NM_DTYPE(nmatrix) == RUBYOBJ) {

    // matrix of Ruby objects -- yield those objects directly
    for (i = 0; i < nm_storage_count_max_elements(s); ++i)
      rb_yield( *((VALUE*)((char*)(s->elements) + i*DTYPE_SIZES[NM_DTYPE(nmatrix)])) );

  } else {
    // We're going to copy the matrix element into a Ruby VALUE and then operate on it. This way user can't accidentally
    // modify it and cause a seg fault.

    for (i = 0; i < nm_storage_count_max_elements(s); ++i) {
      v = rubyobj_from_cval((char*)(s->elements) + i*DTYPE_SIZES[NM_DTYPE(nmatrix)], NM_DTYPE(nmatrix)).rval;
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
  case DENSE_STORE:
    return nm_each_dense(nm);
  default:
    rb_raise(rb_eNotImpError, "only dense matrix's each method works right now");
  }
}



/*
 * Equality operator. Returns a single true or false value indicating whether
 * the matrices are equivalent.
 *
 * For elementwise, use == instead.
 *
 * This method will raise an exception if dimensions do not match.
 */
static VALUE nm_eqeq(VALUE left, VALUE right) {
  NMATRIX *l, *r;

  CheckNMatrixType(left);
  CheckNMatrixType(right);

  UnwrapNMatrix(left, l);
  UnwrapNMatrix(right, r);

  if (l->stype != r->stype)
    rb_raise(rb_eNotImpError, "comparison between different matrix stypes not yet implemented");

  bool result = false;

  switch(l->stype) {
  case DENSE_STORE:
    result = nm_dense_storage_eqeq(l->storage, r->storage);
    break;
  case LIST_STORE:
    result = nm_list_storage_eqeq(l->storage, r->storage);
    break;
  case YALE_STORE:
    result = nm_yale_storage_eqeq(l->storage, r->storage);
    break;
  }

  return result ? Qtrue : Qfalse;
}

/*
 * Simple n-dimensional matrix-matrix addition.
 */
static VALUE nm_ew_add(VALUE left_val, VALUE right_val) {
	return elementwise_op(nm::EW_ADD, left_val, right_val);
}

/*
 * Simple n-dimensional matrix-matrix subtraction.
 */
static VALUE nm_ew_subtract(VALUE left_val, VALUE right_val) {
	return elementwise_op(nm::EW_SUB, left_val, right_val);
}

/*
 * Simple n-dimensional matrix-matrix multiplication.
 */
static VALUE nm_ew_multiply(VALUE left_val, VALUE right_val) {
	return elementwise_op(nm::EW_MUL, left_val, right_val);
}

/*
 * Simple n-dimensional matrix-matrix division.
 */
static VALUE nm_ew_divide(VALUE left_val, VALUE right_val) {
	return elementwise_op(nm::EW_DIV, left_val, right_val);
}

/*
 * Simple n-dimensional matrix-matrix module.
 */
/*
static VALUE nm_ew_mod(VALUE left_val, VALUE right_val) {
	return elementwise_op(nm::EW_MOD, left_val, right_val);
}
*/

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

/*
 * Transform the matrix (in-place) to its complex conjugate. Only works on complex matrices.
 *
 * FIXME: For non-complex matrices, someone needs to implement a non-in-place complex conjugate (which doesn't use a bang).
 * Bang should imply that no copy is being made, even temporarily.
 */
static VALUE nm_complex_conjugate_bang(VALUE self) {
  NMATRIX* m;
  void* elem;
  size_t size, p;

  UnwrapNMatrix(self, m);

  if (m->stype == DENSE_STORE) {

    size = nm_storage_count_max_elements(NM_STORAGE(self));
    elem = NM_STORAGE_DENSE(self)->elements;

  } else if (m->stype == YALE_STORE) {

    size = nm_yale_storage_get_size(NM_STORAGE_YALE(self));
    elem = NM_STORAGE_YALE(self)->a;

  } else {
    rb_raise(rb_eNotImpError, "please cast to yale or dense (complex) first");
  }

  // Walk through and negate the imaginary component
  if (NM_DTYPE(self) == COMPLEX64) {

    for (p = 0; p < size; ++p) {
      reinterpret_cast<nm::Complex64*>(elem)[p].i = -reinterpret_cast<nm::Complex64*>(elem)[p].i;
    }

  } else if (NM_DTYPE(self) == COMPLEX128) {

    for (p = 0; p < size; ++p) {
      reinterpret_cast<nm::Complex128*>(elem)[p].i = -reinterpret_cast<nm::Complex128*>(elem)[p].i;
    }

  } else {
    rb_raise(nm_eDataTypeError, "can only calculate in-place complex conjugate on matrices of type :complex64 or :complex128");
  }

  return self;
}


/*
 * Helper function for creating a matrix. You have to create the storage and pass it in, but you don't
 * need to worry about deleting it.
 */
NMATRIX* nm_create(stype_t stype, STORAGE* storage) {
  NMATRIX* mat = ALLOC(NMATRIX);

  mat->stype   = stype;
  mat->storage = storage;

  return mat;
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

  if (argc < 2) {
  	rb_raise(rb_eArgError, "Expected 2-4 arguments (or 7 for internal Yale creation)");
  	return Qnil;
  }

  /* First, determine stype (dense by default) */
  stype_t stype;
  size_t  offset = 0;

  if (!SYMBOL_P(argv[0]) && TYPE(argv[0]) != T_STRING) {
    stype = DENSE_STORE;
    
  } else {
  	// 0: String or Symbol
    stype  = interpret_stype(argv[0]);
    offset = 1;
  }

  // If there are 7 arguments and Yale, refer to a different init function with fewer sanity checks.
  if (argc == 7) {
  	if (stype == YALE_STORE) {
			return nm_init_yale_from_old_yale(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], nm);
			
		} else {
			rb_raise(rb_eArgError, "Expected 2-4 arguments (or 7 for internal Yale creation)");
		}
  }
	
	// 1: Array or Fixnum
	size_t dim;
  size_t* shape = interpret_shape(argv[offset], &dim);

  // 2-3: dtype
  dtype_t dtype = interpret_dtype(argc-1-offset, argv+offset+1, stype);

  size_t init_cap = 0, init_val_len = 0;
  void* init_val  = NULL;
  if (NM_RUBYVAL_IS_NUMERIC(argv[1+offset]) || TYPE(argv[1+offset]) == T_ARRAY) {
  	// Initial value provided (could also be initial capacity, if yale).
  	
    if (stype == YALE_STORE) {
      init_cap = FIX2UINT(argv[1+offset]);
      
    } else {
    	// 4: initial value / dtype
      init_val = interpret_initial_value(argv[1+offset], dtype);
      
      if (TYPE(argv[1+offset]) == T_ARRAY) 	init_val_len = RARRAY_LEN(argv[1+offset]);
      else                                  init_val_len = 1;
    }
    
  } else {
  	// DType is RUBYOBJ.
  	
    if (stype == DENSE_STORE) {
    	/*
    	 * No need to initialize dense with any kind of default value unless it's
    	 * an RUBYOBJ matrix.
    	 */
      if (dtype == RUBYOBJ) {
      	// Pretend [nil] was passed for RUBYOBJ.
      	init_val = ALLOC(VALUE);
        *(VALUE*)init_val = Qnil;
        
        init_val_len = 1;
        
      } else {
      	init_val = NULL;
      }
    } else if (stype == LIST_STORE) {
    	init_val = ALLOC_N(char, DTYPE_SIZES[dtype]);
      std::memset(init_val, 0, DTYPE_SIZES[dtype]);
    }
  }
	
  // TODO: Update to allow an array as the initial value.
	NMATRIX* nmatrix;
  UnwrapNMatrix(nm, nmatrix);

  nmatrix->stype = stype;
  
  switch (stype) {
  	case DENSE_STORE:
  		nmatrix->storage = (STORAGE*)nm_dense_storage_create(dtype, shape, dim, init_val, init_val_len);
  		break;
  		
  	case LIST_STORE:
  		nmatrix->storage = (STORAGE*)nm_list_storage_create(dtype, shape, dim, init_val);
  		break;
  		
  	case YALE_STORE:
  		nmatrix->storage = (STORAGE*)nm_yale_storage_create(dtype, shape, dim, init_cap);
  		nm_yale_storage_init((YALE_STORAGE*)(nmatrix->storage));
  		break;
  }

  return nm;
}

/*
 * Copy constructor for changing dtypes and stypes.
 */
static VALUE nm_init_cast_copy(VALUE self, VALUE new_stype_symbol, VALUE new_dtype_symbol) {
  dtype_t new_dtype = nm_dtype_from_rbsymbol(new_dtype_symbol);
  stype_t new_stype = nm_stype_from_rbsymbol(new_stype_symbol);

  CheckNMatrixType(self);

  NMATRIX *lhs = ALLOC(NMATRIX),
          *rhs;
  lhs->stype = new_stype;

  UnwrapNMatrix( self, rhs );

  // Copy the storage
  STYPE_CAST_COPY_TABLE(cast_copy);
  lhs->storage = cast_copy[lhs->stype][rhs->stype](rhs->storage, new_dtype);

  STYPE_MARK_TABLE(mark);

  return Data_Wrap_Struct(cNMatrix, mark[lhs->stype], nm_delete, lhs);
}


/*
 * Copy constructor for transposing.
 */
static VALUE nm_init_transposed(VALUE self) {
  static STORAGE* (*storage_copy_transposed[nm::NUM_STYPES])(const STORAGE* rhs_base) = {
    nm_dense_storage_copy_transposed,
    nm_list_storage_copy_transposed,
    nm_yale_storage_copy_transposed
  };

  NMATRIX* lhs = nm_create( NM_STYPE(self),
                            storage_copy_transposed[NM_STYPE(self)]( NM_STORAGE(self) )
                          );

  STYPE_MARK_TABLE(mark);

  return Data_Wrap_Struct(cNMatrix, mark[lhs->stype], nm_delete, lhs);
}


/*
 * Copy constructor for no change of dtype or stype (used for #initialize_copy hook).
 */
static VALUE nm_init_copy(VALUE copy, VALUE original) {
  NMATRIX *lhs, *rhs;

  CheckNMatrixType(original);

  if (copy == original) return copy;

  UnwrapNMatrix( original, rhs );
  UnwrapNMatrix( copy,     lhs );

  lhs->stype = rhs->stype;

  // Copy the storage
  STYPE_CAST_COPY_TABLE(ttable);
  lhs->storage = ttable[lhs->stype][rhs->stype](rhs->storage, rhs->storage->dtype);

  return copy;
}


/*
 * Create a new NMatrix helper for handling internal ia, ja, and a arguments.
 *
 * This constructor is only called by Ruby code, so we can skip most of the
 * checks.
 */
static VALUE nm_init_yale_from_old_yale(VALUE shape, VALUE dtype, VALUE ia, VALUE ja, VALUE a, VALUE from_dtype, VALUE nm) {
  size_t dim     = 2;
  size_t* shape_  = interpret_shape(shape, &dim);
  dtype_t dtype_  = nm_dtype_from_rbsymbol(dtype);
  char *ia_       = RSTRING_PTR(ia),
       *ja_       = RSTRING_PTR(ja),
       *a_        = RSTRING_PTR(a);
  dtype_t from_dtype_ = nm_dtype_from_rbsymbol(from_dtype);
  NMATRIX* nmatrix;

  UnwrapNMatrix( nm, nmatrix );

  nmatrix->stype   = YALE_STORE;
  nmatrix->storage = (STORAGE*)nm_yale_storage_create_from_old_yale(dtype_, shape_, ia_, ja_, a_, from_dtype_);

  return nm;
}

/*
 * Check to determine whether matrix is a reference to another matrix.
 */
static VALUE nm_is_ref(VALUE self) {
	// Refs only allowed for dense and list matrices.
  if (NM_STYPE(self) == DENSE_STORE) {
    return (NM_DENSE_SRC(self) == NM_STORAGE(self)) ? Qfalse : Qtrue;
  }

  if (NM_STYPE(self) == LIST_STORE) {
    return (NM_LIST_SRC(self) == NM_STORAGE(self)) ? Qfalse : Qtrue;
  }

  return Qfalse;
}



/*
 * Access the contents of an NMatrix at given coordinates, using copying.
 *
 *     n.slice(3,3)  # => 5.0
 *     n.slice(0..1,0..1) #=> matrix [2,2]
 *
 */
static VALUE nm_mget(int argc, VALUE* argv, VALUE self) {
  static void* (*ttable[nm::NUM_STYPES])(STORAGE*, SLICE*) = {
    nm_dense_storage_get,
    nm_list_storage_get,
    nm_yale_storage_get
  };
  
  return nm_xslice(argc, argv, ttable[NM_STYPE(self)], nm_delete, self);
}

/*
 * Access the contents of an NMatrix at given coordinates by reference.
 *
 *     n[3,3]  # => 5.0
 *     n[0..1,0..1] #=> matrix [2,2]
 *
 */
static VALUE nm_mref(int argc, VALUE* argv, VALUE self) {
  static void* (*ttable[nm::NUM_STYPES])(STORAGE*, SLICE*) = {
    nm_dense_storage_ref,
    nm_list_storage_ref,
    nm_yale_storage_ref
  };
  return nm_xslice(argc, argv, ttable[NM_STYPE(self)], nm_delete_ref, self);
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
static VALUE nm_mset(int argc, VALUE* argv, VALUE self) {
  size_t dim = argc - 1; // last arg is the value

  if (argc <= 1) {
    rb_raise(rb_eArgError, "Expected coordinates and r-value");

  } else if (NM_DIM(self) == dim) {

    SLICE* slice = get_slice(dim, argv, self);

    void* value = rubyobj_to_cval(argv[dim], NM_DTYPE(self));

    // FIXME: Can't use a function pointer table here currently because these functions have different
    // signatures (namely the return type).
    switch(NM_STYPE(self)) {
    case DENSE_STORE:
      nm_dense_storage_set(NM_STORAGE(self), slice, value);
      break;
    case LIST_STORE:
      // Remove if it's a zero, insert otherwise
      if (!std::memcmp(value, NM_STORAGE_LIST(self)->default_val, DTYPE_SIZES[NM_DTYPE(self)])) {
        free(value);
        value = nm_list_storage_remove(NM_STORAGE(self), slice);
        free(value);
      } else {
        nm_list_storage_insert(NM_STORAGE(self), slice, value);
      }
      break;
    case YALE_STORE:
      nm_yale_storage_set(NM_STORAGE(self), slice, value);
      break;
    }

    return argv[dim];

  } else if (NM_DIM(self) < dim) {
    rb_raise(rb_eArgError, "Coordinates given exceed number of matrix dimensions");
  } else {
    rb_raise(rb_eNotImpError, "Slicing not supported yet");
  }
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

  if (NM_RUBYVAL_IS_NUMERIC(right_v))
    return matrix_multiply_scalar(left, right_v);

  else if (TYPE(right_v) == T_ARRAY)
    rb_raise(rb_eNotImpError, "for matrix-vector multiplication, please use an NVector instead of an Array for now");

  //if (RDATA(right_v)->dfree != (RUBY_DATA_FUNC)nm_delete) {
  else { // both are matrices
    CheckNMatrixType(right_v);
    UnwrapNMatrix( right_v, right );

    if (left->storage->shape[1] != right->storage->shape[0])
      rb_raise(rb_eArgError, "incompatible dimensions");

    if (left->stype != right->stype)
      rb_raise(rb_eNotImpError, "matrices must have same stype");

    return matrix_multiply(left, right);

  } 

  return Qnil;
}

/*
 * In-place LU factorization of a matrix.
 */
static VALUE nm_factorize_lu_bang(VALUE self) {
  if (NM_STYPE(self) != DENSE_STORE) {
    rb_raise(rb_eNotImpError, "only implemented for dense storage");
  }

  if (NM_DIM(self) != 2) {
    rb_raise(rb_eNotImpError, "matrix is not 2-dimensional");
  }

  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::clapack_getrf, bool, const int m, const int n, void* a, const int lda, int* ipiv);

  int* ipiv = ALLOCA_N(int, std::min(NM_SHAPE0(self), NM_SHAPE1(self)));

  bool result = ttable[NM_DTYPE(self)](NM_SHAPE0(self), NM_SHAPE1(self), NM_STORAGE_DENSE(self)->elements, NM_SHAPE1(self), ipiv) ? Qtrue : Qfalse;

  return self;
}

/*
 * Get the number of dimensions of a matrix.
 *
 * In other words, if you set your matrix to be 3x4, the dim is 2. If the
 * matrix was initialized as 3x4x3, the dim is 3.
 *
 * This function may lie slightly for NVectors, which are internally stored as
 * dim 2 (and have an orientation), but act as if they're dim 1.
 */
static VALUE nm_dim(VALUE self) {
  return INT2FIX(NM_STORAGE(self)->dim);
}

/*
 * Get the shape (dimensions) of a matrix.
 */
static VALUE nm_shape(VALUE self) {
  STORAGE* s   = NM_STORAGE(self);
  size_t index;

  // Copy elements into a VALUE array and then use those to create a Ruby array with rb_ary_new4.
  VALUE* shape = ALLOCA_N(VALUE, s->dim);
  for (index = 0; index < s->dim; ++index)
    shape[index] = INT2FIX(s->shape[index]);

  return rb_ary_new4(s->dim, shape);
}

/*
 * Get the storage type (stype) of a matrix, e.g., :yale, :dense, or :list.
 */
static VALUE nm_stype(VALUE self) {
  ID stype = rb_intern(STYPE_NAMES[NM_STYPE(self)]);
  return ID2SYM(stype);
}

/*
 * Is this matrix symmetric?
 */
static VALUE nm_symmetric(VALUE self) {
  return is_symmetric(self, false);
}

/*
 * Get a slice of an NMatrix.
 */
static VALUE nm_xslice(int argc, VALUE* argv, void* (*slice_func)(STORAGE*, SLICE*), void (*delete_func)(NMATRIX*), VALUE self) {
  VALUE result = Qnil;

  if (NM_DIM(self) == (size_t)(argc)) {
    SLICE* slice = get_slice((size_t)(argc), argv, self);

    // TODO: Slice for Yale types

    if (slice->single) {

      static void* (*ttable[nm::NUM_STYPES])(STORAGE*, SLICE*) = {
        nm_dense_storage_ref,
        nm_list_storage_ref,
        nm_yale_storage_ref
      };

      /* // Debugging for slice
      fprintf(stderr, "single: ");
      for (size_t i = 0; i < NM_DIM(self); ++i) {
        fprintf(stderr, "%u(%u) ", slice->coords[i], slice->lengths[i]);
      }
      fprintf(stderr, "\n");
      */

      //DENSE_STORAGE* s = NM_STORAGE_DENSE(self);

      if (NM_DTYPE(self) == RUBYOBJ)  result = *reinterpret_cast<VALUE*>( ttable[NM_STYPE(self)](NM_STORAGE(self), slice) );
      else                            result = rubyobj_from_cval( ttable[NM_STYPE(self)](NM_STORAGE(self), slice), NM_DTYPE(self) ).rval;

    } else {

      if (NM_STYPE(self) == DENSE_STORE 
          || NM_STYPE(self) == LIST_STORE) {
        STYPE_MARK_TABLE(mark_table);

        NMATRIX* mat = ALLOC(NMATRIX);
        mat->stype = NM_STYPE(self);
        mat->storage = (STORAGE*)((*slice_func)( NM_STORAGE(self), slice ));
        result = Data_Wrap_Struct(cNMatrix, mark_table[mat->stype], delete_func, mat);
      } else {
        rb_raise(rb_eNotImpError, "slicing only implemented for dense and list so far");
      }
    }

    free(slice);

  } else if (NM_DIM(self) < (size_t)(argc)) {
    rb_raise(rb_eArgError, "Coordinates given exceed number of matrix dimensions");
  } else {
    rb_raise(rb_eNotImpError, "This type slicing not supported yet");
  }

  return result;
}

//////////////////////
// Helper Functions //
//////////////////////

static VALUE elementwise_op(nm::ewop_t op, VALUE left_val, VALUE right_val) {
	STYPE_MARK_TABLE(mark);
	
	static STORAGE* (*ew_op[nm::NUM_STYPES])(nm::ewop_t, const STORAGE*, const STORAGE*) = {
		nm_dense_storage_ew_op,
		nm_list_storage_ew_op,
		NULL
	};
	
	NMATRIX* result = ALLOC(NMATRIX);
	
	CheckNMatrixType(left_val);
	CheckNMatrixType(right_val);

	// Check that the left- and right-hand sides have the same dimensionality.
	if (NM_DIM(left_val) != NM_DIM(right_val)) {
		rb_raise(rb_eArgError, "The left- and right-hand sides of the operation must have the same dimensionality.");
	}
	
	// Check that the left- and right-hand sides have the same shape.
	if (memcmp(&NM_SHAPE(left_val, 0), &NM_SHAPE(right_val, 0), sizeof(size_t) * NM_DIM(left_val)) != 0) {
		rb_raise(rb_eArgError, "The left- and right-hand sides of the operation must have the same shape.");
	}

	NMATRIX* left, * right;
	UnwrapNMatrix(left_val, left);
	UnwrapNMatrix(right_val, right);
	
	if (left->stype == right->stype) {
		
		if (ew_op[left->stype] == NULL) {
			rb_raise(rb_eArgError, "Element-wise operations are not supported for the given storage type.");
		}
		
		result->storage	= ew_op[left->stype](op, reinterpret_cast<STORAGE*>(left->storage), reinterpret_cast<STORAGE*>(right->storage));
		result->stype		= left->stype;
		
	} else {
		rb_raise(rb_eArgError, "Element-wise operations are not currently supported between matrices with differing stypes.");
	}

	return Data_Wrap_Struct(cNMatrix, mark[result->stype], nm_delete, result);
}

/*
 * Check to determine whether matrix is a reference to another matrix.
 */
bool is_ref(const NMATRIX* matrix) {
  // FIXME: Needs to work for other types
  if (matrix->stype != DENSE_STORE) {
    return false;
  }
  
  return ((DENSE_STORAGE*)(matrix->storage))->src != matrix->storage;
}

/*
 * Helper function for nm_symmetric and nm_hermitian.
 */
static VALUE is_symmetric(VALUE self, bool hermitian) {
  NMATRIX* m;
  UnwrapNMatrix(self, m);

  if (m->storage->shape[0] == m->storage->shape[1] and m->storage->dim == 2) {
		if (NM_STYPE(self) == DENSE_STORE) {
      if (hermitian) {
        nm_dense_storage_is_hermitian((DENSE_STORAGE*)(m->storage), m->storage->shape[0]);
        
      } else {
      	nm_dense_storage_is_symmetric((DENSE_STORAGE*)(m->storage), m->storage->shape[0]);
      }
      
    } else {
      // TODO: Implement, at the very least, yale_is_symmetric. Model it after yale/transp.template.c.
      rb_raise(rb_eNotImpError, "symmetric? and hermitian? only implemented for dense currently");
    }

  }

  return Qfalse;
}


///////////////////////
// Utility Functions //
///////////////////////

/*
 * Guess the data type given a value.
 *
 * TODO: Probably needs some work for Bignum.
 */
static dtype_t dtype_guess(VALUE v) {
  switch(TYPE(v)) {
  case T_TRUE:
  case T_FALSE:
    return BYTE;
    
  case T_STRING:
    if (RSTRING_LEN(v) == 1) {
    	return BYTE;
    	
    } else {
    	rb_raise(rb_eArgError, "Strings of length > 1 may not be stored in a matrix.");
    }

#if SIZEOF_INT == 8
  case T_FIXNUM:
    return INT64;
    
  case T_RATIONAL:
    return RATIONAL128;
    
#else
# if SIZEOF_INT == 4
  case T_FIXNUM:
    return INT32;
    
  case T_RATIONAL:
    return RATIONAL64;
    
#else
  case T_FIXNUM:
    return INT16;
    
  case T_RATIONAL:
    return RATIONAL32;
# endif
#endif

  case T_BIGNUM:
    return INT64;

#if SIZEOF_FLOAT == 4
  case T_COMPLEX:
    return COMPLEX128;
    
  case T_FLOAT:
    return FLOAT64;
    
#else
# if SIZEOF_FLOAT == 2
  case T_COMPLEX:
    return COMPLEX64;
    
  case T_FLOAT:
    return FLOAT32;
# endif
#endif

  case T_ARRAY:
  	/*
  	 * May be passed for dense -- for now, just look at the first element.
  	 * 
  	 * TODO: Look at entire array for most specific type.
  	 */
  	
    return dtype_guess(RARRAY_PTR(v)[0]);

  case T_NIL:
  default:
    rb_raise(rb_eArgError, "Unable to guess a data type from provided parameters; data type must be specified manually.");
  }
}

/*
 * Documentation goes here.
 */
static SLICE* get_slice(size_t dim, VALUE* c, VALUE self) {
  size_t r;
  VALUE beg, end;
  int exl;

  SLICE* slice = ALLOC(SLICE);
  slice->coords = ALLOC_N(size_t,dim);
  slice->lengths = ALLOC_N(size_t, dim);
  slice->single = true;

  for (r = 0; r < dim; ++r) {

    if (FIXNUM_P(c[r])) { // this used CLASS_OF before, which is inefficient for fixnum

      slice->coords[r]  = FIX2UINT(c[r]);
      slice->lengths[r] = 1;

    } else if (CLASS_OF(c[r]) == rb_cRange) {
        rb_range_values(c[r], &beg, &end, &exl);
        slice->coords[r]  = FIX2UINT(beg);
        slice->lengths[r] = FIX2UINT(end) - slice->coords[r] + 1;

        // Exclude last element for a...b range
        if (exl)
          slice->lengths[r] -= 1;

        slice->single     = false;

    } else {
      rb_raise(rb_eArgError, "cannot slice using class %s, needs a number or range or something", rb_obj_classname(c[r]));
    }

    if (slice->coords[r] + slice->lengths[r] > NM_SHAPE(self,r))
      rb_raise(rb_eArgError, "out of range");
  }

  return slice;
}

#ifdef BENCHMARK
/*
 * A simple function used when benchmarking NMatrix.
 */
static double get_time(void) {
  struct timeval t;
  struct timezone tzp;
  
  gettimeofday(&t, &tzp);
  
  return t.tv_sec + t.tv_usec*1e-6;
}
#endif

/*
 * The argv parameter will be either 1 or 2 elements.  If 1, could be either
 * initial or dtype.  If 2, is initial and dtype. This function returns the
 * dtype.
 */
static dtype_t interpret_dtype(int argc, VALUE* argv, stype_t stype) {
  int offset;
  
  switch (argc) {
  	case 1:
  		offset = 0;
  		break;
  	
  	case 2:
  		offset = 1;
  		break;
  		
  	default:
  		rb_raise(rb_eArgError, "Need an initial value or a dtype.");
  		break;
  }

  if (SYMBOL_P(argv[offset])) {
  	return nm_dtype_from_rbsymbol(argv[offset]);
  	
  } else if (TYPE(argv[offset]) == T_STRING) {
  	return nm_dtype_from_rbstring(StringValue(argv[offset]));
  	
  } else if (stype == YALE_STORE) {
  	rb_raise(rb_eArgError, "Yale storage class requires a dtype.");
  	
  } else {
  	return dtype_guess(argv[0]);
  }
}

/*
 * Convert an Ruby value or an array of Ruby values into initial C values.
 */
static void* interpret_initial_value(VALUE arg, dtype_t dtype) {
  unsigned int index;
  void* init_val;
  
  if (TYPE(arg) == T_ARRAY) {
  	// Array
    
    init_val = ALLOC_N(int8_t, DTYPE_SIZES[dtype] * RARRAY_LEN(arg));
    for (index = 0; index < RARRAY_LEN(arg); ++index) {
    	rubyval_to_cval(RARRAY_PTR(arg)[index], dtype, (char*)init_val + (index * DTYPE_SIZES[dtype]));
    }
    
  } else {
  	// Single value
  	
    init_val = rubyobj_to_cval(arg, dtype);
  }

  return init_val;
}

/*
 * Convert the shape argument, which may be either a Ruby value or an array of
 * Ruby values, into C values.  The second argument is where the dimensionality
 * of the matrix will be stored.  The function itself returns a pointer to the
 * array describing the shape, which must be freed manually.
 */
static size_t* interpret_shape(VALUE arg, size_t* dim) {
  size_t* shape;

  if (TYPE(arg) == T_ARRAY) {
    *dim = RARRAY_LEN(arg);
    shape = ALLOC_N(size_t, *dim);
    
    for (size_t index = 0; index < *dim; ++index) {
      shape[index] = FIX2UINT( RARRAY_PTR(arg)[index] );
    }
    
  } else if (FIXNUM_P(arg)) {
    *dim = 2;
    shape = ALLOC_N(size_t, *dim);
    
    shape[0] = FIX2UINT(arg);
    shape[1] = FIX2UINT(arg);
    
  } else {
    rb_raise(rb_eArgError, "Expected an array of numbers or a single Fixnum for matrix shape");
  }

  return shape;
}

/*
 * Convert a Ruby symbol or string into an storage type.
 */
static stype_t interpret_stype(VALUE arg) {
  if (SYMBOL_P(arg)) {
  	return nm_stype_from_rbsymbol(arg);
  	
  } else if (TYPE(arg) == T_STRING) {
  	return nm_stype_from_rbstring(StringValue(arg));
  	
  } else {
  	rb_raise(rb_eArgError, "Expected storage type");
  }
}


//////////////////
// Math Helpers //
//////////////////


STORAGE* matrix_storage_cast_alloc(NMATRIX* matrix, dtype_t new_dtype) {
  if (matrix->storage->dtype == new_dtype && !is_ref(matrix))
    return matrix->storage;

  STYPE_CAST_COPY_TABLE(cast_copy_storage);
  return cast_copy_storage[matrix->stype][matrix->stype](matrix->storage, new_dtype);
}


STORAGE_PAIR binary_storage_cast_alloc(NMATRIX* left_matrix, NMATRIX* right_matrix) {
  STORAGE_PAIR casted;
  dtype_t new_dtype = Upcast[left_matrix->storage->dtype][right_matrix->storage->dtype];

  casted.left  = matrix_storage_cast_alloc(left_matrix, new_dtype);
  casted.right = matrix_storage_cast_alloc(right_matrix, new_dtype);

  return casted;
}


static VALUE matrix_multiply_scalar(NMATRIX* left, VALUE scalar) {
  rb_raise(rb_eNotImpError, "matrix-scalar multiplication not implemented yet");
  return Qnil;
}

static VALUE matrix_multiply(NMATRIX* left, NMATRIX* right) {
  ///TODO: multiplication for non-dense and/or non-decimal matrices

  // Make sure both of our matrices are of the correct type.
  STORAGE_PAIR casted = binary_storage_cast_alloc(left, right);

  size_t*  resulting_shape   = ALLOC_N(size_t, 2);
  resulting_shape[0] = left->storage->shape[0];
  resulting_shape[1] = right->storage->shape[1];

  // Sometimes we only need to use matrix-vector multiplication (e.g., GEMM versus GEMV). Find out.
  bool vector = false;
  if (resulting_shape[1] == 1) vector = true;

  static STORAGE* (*storage_matrix_multiply[nm::NUM_STYPES])(const STORAGE_PAIR&, size_t*, bool) = {
    nm_dense_storage_matrix_multiply,
    nm_list_storage_matrix_multiply,
    nm_yale_storage_matrix_multiply
  };

  STORAGE* resulting_storage = storage_matrix_multiply[left->stype](casted, resulting_shape, vector);
  NMATRIX* result = nm_create(left->stype, resulting_storage);

  // Free any casted-storage we created for the multiplication.
  // TODO: Can we make the Ruby GC take care of this stuff now that we're using it?
  // If we did that, we night not have to re-create these every time, right? Or wrong? Need to do
  // more research.
  static void (*free_storage[nm::NUM_STYPES])(STORAGE*) = {
    nm_dense_storage_delete,
    nm_list_storage_delete,
    nm_yale_storage_delete
  };

  if (left->storage != casted.left)   free_storage[result->stype](casted.left);
  if (right->storage != casted.right) free_storage[result->stype](casted.right);


  STYPE_MARK_TABLE(mark_table);

  if (result) return Data_Wrap_Struct(cNMatrix, mark_table[result->stype], nm_delete, result);
  return Qnil; // Only if we try to multiply list matrices should we return Qnil.
}





/*
 * Calculate the exact determinant of a dense matrix.
 *
 * Returns nil for dense matrices which are not square or number of dimensions other than 2.
 *
 * Note: Currently only implemented for 2x2 and 3x3 matrices.
 */
static VALUE nm_det_exact(VALUE self) {
  if (NM_STYPE(self) != DENSE_STORE) rb_raise(nm_eStorageTypeError, "can only calculate exact determinant for dense matrices");

  if (NM_DIM(self) != 2 || NM_SHAPE0(self) != NM_SHAPE1(self)) return Qnil;

  // Calculate the determinant and then assign it to the return value
  void* result = ALLOCA_N(char, DTYPE_SIZES[NM_DTYPE(self)]);
  nm_math_det_exact(NM_SHAPE0(self), NM_STORAGE_DENSE(self)->elements, NM_SHAPE0(self), NM_DTYPE(self), result);

  return rubyobj_from_cval(result, NM_DTYPE(self)).rval;
}




/////////////////
// Exposed API //
/////////////////


/*
 * Create a dense matrix. Used by the NMatrix GSL fork. Unlike nm_create, this one copies all of the
 * arrays and such passed in -- so you don't have to allocate and pass a new shape object for every
 * matrix you want to create, for example. Same goes for elements.
 *
 * Returns a properly-wrapped Ruby object as a VALUE.
 *
 * TODO: Add a column-major option for libraries that use column-major matrices.
 */
VALUE rb_nmatrix_dense_create(dtype_t dtype, size_t* shape, size_t dim, void* elements, size_t length) {
  NMATRIX* nm;
  VALUE klass;
  size_t nm_dim;
  size_t* shape_copy;

  // Do not allow a dim of 1; if dim == 1, this should probably be an NVector instead, but that still has dim 2.
  if (dim == 1) {
    klass					= cNVector;
    nm_dim				= 2;
    shape_copy		= ALLOC_N(size_t, nm_dim);
    shape_copy[0]	= shape[0];
    shape_copy[1]	= 1;
    
  } else {
    klass				= cNMatrix;
    nm_dim			= dim;
    shape_copy	= ALLOC_N(size_t, nm_dim);
    memcpy(shape_copy, shape, sizeof(size_t)*nm_dim);
  }

  // Copy elements
  void* elements_copy = ALLOC_N(char, DTYPE_SIZES[dtype]*length);
  memcpy(elements_copy, elements, DTYPE_SIZES[dtype]*length);

  // allocate and create the matrix and its storage
  nm = nm_create(DENSE_STORE, nm_dense_storage_create(dtype, shape_copy, dim, elements_copy, length));

  // tell Ruby about the matrix and its storage, particularly how to garbage collect it.
  return Data_Wrap_Struct(klass, nm_dense_storage_mark, nm_dense_storage_delete, nm);
}


/*
 * Create a dense vector. Used by the NMatrix GSL fork.
 *
 * Basically just a convenience wrapper for rb_nmatrix_dense_create().
 *
 * Returns a properly-wrapped Ruby NVector object as a VALUE.
 *
 * TODO: Add a transpose option for setting the orientation of the vector?
 */
VALUE rb_nvector_dense_create(dtype_t dtype, void* elements, size_t length) {
  size_t dim = 1, shape = length;
  return rb_nmatrix_dense_create(dtype, &shape, dim, elements, length);
}

} // end of extern "C"

