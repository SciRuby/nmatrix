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
// == hb_io.h
//
// Header file for Harwell Boeing input/output support functions.

#ifndef NMATRIX_HB_IO_H
#define NMATRIX_HB_IO_H

#include <cstdlib>
#include <fstream>
#include <ruby.h>

#include "nmatrix.h"

namespace nm {
	namespace hb_io {
		void hb_header_read ( std::ifstream &input, char **title, char **key, int *totcrd, 
		  int *ptrcrd, int *indcrd, int *valcrd, int *rhscrd, char **mxtype, int *nrow, 
		  int *ncol, int *nnzero, int *neltvl, char **ptrfmt, char **indfmt, char **valfmt, 
		  char **rhsfmt, char **rhstyp, int *nrhs, int *nrhsix );

		char *s_substring ( char *s, int a, int b );

		void s_trim ( char *s );

		void hb_structure_read ( std::ifstream &input, int ncol, char *mxtype, int nnzero, 
		  int neltvl, int ptrcrd, char *ptrfmt, int indcrd, char *indfmt, 
		  int colptr[], int rowind[] );

		void hb_values_read ( std::ifstream &input, int valcrd, char *mxtype, int nnzero,
		  int neltvl, char *valfmt, double values[] );

		void s_to_format ( char *s, int *r, char *code, int *w, int *m );

		int i4_min ( int i1, int i2 );

		int s_len_trim ( char *s );

		int ch_to_digit ( char c );

		bool ch_is_digit ( char c );

		bool ch_is_format_code ( char c );

		bool ch_eqi ( char c1, char c2 );
	}
}

#endif