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
// == io.h
//
// Header file for input/output support functions.

#ifndef NMATRIX_IO_H
#define NMATRIX_IO_H

/*
 * Project Includes
 */

#include "nmatrix.h"

#include "data/data.h"
#include "storage/storage.h"

/*
 * Extern Types
 */
extern const char* const DTYPE_NAMES[nm::NUM_DTYPES];

namespace nm { namespace io {
  /*
   * Types
   */
  enum matlab_dtype_t {
    miINT8 = 1,
    miUINT8 = 2,
    miINT16 = 3,
    miUINT16 = 4,
    miINT32 = 5,
    miUINT32 = 6,
    miSINGLE = 7,
    miDOUBLE = 9,
    miINT64 = 12,
    miUINT64 = 13,
    miMATRIX = 14
  };

  /*
   * Constants
   */

  const size_t NUM_MATLAB_DTYPES = 15;
}} // end of namespace nm::io

extern "C" {

  /*
   * C accessors.
   */
  nm::dtype_t nm_dtype_from_rbsymbol(VALUE sym);
  nm::dtype_t nm_dtype_from_rbstring(VALUE str);
  nm::stype_t nm_stype_from_rbsymbol(VALUE sym);
  nm::stype_t nm_stype_from_rbstring(VALUE str);

  void nm_init_io(void);


  /*
   * Macros.
   */

  /*
   * Macro for a function pointer table between NMatrix dtypes and MATLAB dtypes.
   *
   * You can't convert as freely between these two as you can between NMatrix dtypes, but there's no reason to. MATLAB
   * stores its complex numbers in two separate arrays, for example, not as a single unit of data. If you want to convert
   * to a VALUE, convert first to an appropriate integer or float type.
   *
   * FIXME: Maybe be a little more selective about which conversions we DO allow. This is really just for loading an
   * already-constructed MATLAB matrix into memory, and most of these functions will never get called.
   */
  #define NM_MATLAB_DTYPE_TEMPLATE_TABLE(name,fun,ret,...)    \
  static ret (*(name)[7][nm::io::NUM_MATLAB_DTYPES])(__VA_ARGS__) = {  \
      { NULL, fun<uint8_t,int8_t>, fun<uint8_t,uint8_t>, fun<uint8_t,int16_t>, fun<uint8_t,uint16_t>, fun<uint8_t,int32_t>, fun<uint8_t,uint32_t>, fun<uint8_t,float>, NULL, fun<uint8_t,double>, NULL, NULL, fun<uint8_t,int64_t>, fun<uint8_t,uint64_t>, NULL },  \
      { NULL, fun<int8_t,int8_t>, fun<int8_t,uint8_t>, fun<int8_t,int16_t>, fun<int8_t,uint16_t>, fun<int8_t,int32_t>, fun<int8_t,uint32_t>, fun<int8_t,float>, NULL, fun<int8_t,double>, NULL, NULL, fun<int8_t,int64_t>, fun<int8_t,uint64_t>, NULL },            \
      { NULL, fun<int16_t,int8_t>, fun<int16_t,uint8_t>, fun<int16_t,int16_t>, fun<int16_t,uint16_t>, fun<int16_t,int32_t>, fun<int16_t,uint32_t>, fun<int16_t,float>, NULL, fun<int16_t,double>, NULL, NULL, fun<int16_t,int64_t>, fun<int16_t,uint64_t>, NULL },  \
      { NULL, fun<int32_t,int8_t>, fun<int32_t,uint8_t>, fun<int32_t,int16_t>, fun<int32_t,uint16_t>, fun<int32_t,int32_t>, fun<int32_t,uint32_t>, fun<int32_t,float>, NULL, fun<int32_t,double>, NULL, NULL, fun<int32_t,int64_t>, fun<int32_t,uint64_t>, NULL },  \
      { NULL, fun<int64_t,int8_t>, fun<int64_t,uint8_t>, fun<int64_t,int16_t>, fun<int64_t,uint16_t>, fun<int64_t,int32_t>, fun<int64_t,uint32_t>, fun<int64_t,float>, NULL, fun<int64_t,double>, NULL, NULL, fun<int64_t,int64_t>, fun<int64_t,uint64_t>, NULL },  \
      { NULL, fun<float,int8_t>, fun<float,uint8_t>, fun<float,int16_t>, fun<float,uint16_t>, fun<float,int32_t>, fun<float,uint32_t>, fun<float,float>, NULL, fun<float,double>, NULL, NULL, fun<float,int64_t>, fun<float,uint64_t>, NULL },                      \
      { NULL, fun<double,int8_t>, fun<double,uint8_t>, fun<double,int16_t>, fun<double,uint16_t>, fun<double,int32_t>, fun<double,uint32_t>, fun<double,float>, NULL, fun<double,double>, NULL, NULL, fun<double,int64_t>, fun<double,uint64_t>, NULL }             \
    };

  /*
   * Hash#has_key? for symbols. Arguments are: hash (VALUE), string (char*).
   */
  #define RB_HASH_HAS_SYMBOL_KEY(hash, str)   (rb_funcall((hash), rb_intern("has_key?"), 1, ID2SYM(rb_intern(str))) == Qtrue)
}


#endif
