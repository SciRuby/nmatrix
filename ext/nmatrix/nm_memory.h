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
// == nm_memory.h
//
// Macros for memory allocation and freeing

/**
 * We define these macros, which just call the ruby ones, as this makes 
 * debugging memory issues (particularly those involving interaction with
 * the ruby GC) easier, as it's posssible to add debugging code temporarily.
 */
#ifndef __NM_MEMORY_H__
#define __NM_MEMORY_H__

#include <ruby.h>

#define NM_ALLOC(type) (ALLOC(type))

#define NM_ALLOC_N(type, n) (ALLOC_N(type, n))

#define NM_REALLOC_N(var, type, n) (REALLOC_N(var, type, n))

#define NM_ALLOCA_N(type, n) (ALLOCA_N(type, n))

#define NM_FREE(var) (xfree(var))

#define NM_ALLOC_NONRUBY(type) ((type*) malloc(sizeof(type)))

//Defines whether to do conservative gc registrations, i.e. those
//registrations that we're not that sure are necessary.
//#define NM_GC_CONSERVATIVE

#ifdef NM_GC_CONSERVATIVE
#define NM_CONSERVATIVE(statement) (statement)
#else
#define NM_CONSERVATIVE(statement)
#endif //NM_GC_CONSERVATIVE

#endif
