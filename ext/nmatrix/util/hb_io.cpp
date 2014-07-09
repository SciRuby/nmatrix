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
// == util/hb_io.cpp
//
// Main file for definitions of functions required in Harwell Boeing format file IO. This file
// contains some functions taken from http://people.sc.fsu.edu/~jburkardt/c_src/hb_io/hb_io.c, which
// are useful for implementing only the functionality required in NMatrix.

/*
 * Standard Includes
 */

#include "hb_io.h"

namespace nm { 
	namespace hb_io {

		void hb_header_read ( std::ifstream &input, char **title, char **key, int *totcrd, 
		  int *ptrcrd, int *indcrd, int *valcrd, int *rhscrd, char **mxtype, int *nrow, 
		  int *ncol, int *nnzero, int *neltvl, char **ptrfmt, char **indfmt, char **valfmt, 
		  char **rhsfmt, char **rhstyp, int *nrhs, int *nrhsix )
		
		//****************************************************************************80
		//
		//  Purpose:
		//
		//    HB_HEADER_READ reads the header of an HB file.
		//
		//  Discussion:
		//
		//    The user should already have opened the file, and positioned it
		//    to the first record.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    09 April 2004
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Reference:
		//
		//    Iain Duff, Roger Grimes, John Lewis,
		//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
		//    October 1992.
		//
		//  Parameters:
		//
		//    Input, ifstream &INPUT, the unit from which data is read.
		//
		//    Output, char *TITLE, a 72 character title for the matrix.
		//
		//    Output, char *KEY, an 8 character identifier for the matrix.
		//
		//    Output, int *TOTCRD, the total number of lines of data.
		//
		//    Output, int *PTRCRD, the number of input lines for pointers.
		//
		//    Output, int *INDCRD, the number of input lines for row indices.
		//
		//    Output, int *VALCRD, the number of input lines for numerical values.
		//
		//    Output, int *RHSCRD, the number of input lines for right hand sides.
		//
		//    Output, char *MXTYPE, the 3 character matrix type.
		//    First character is R for Real, C for complex, P for pattern only.
		//    Second character is S for symmetric, U for unsymmetric, H for
		//      Hermitian, Z for skew symmetric, R for rectangular.
		//    Third character is A for assembled and E for unassembled
		//      finite element matrices.
		//
		//    Output, int *NROW, the number of rows or variables.
		//
		//    Output, int *NCOL, the number of columns or elements.
		//
		//    Output, int *NNZERO.  In the case of assembled sparse matrices,
		//    this is the number of nonzeroes.  In the case of unassembled finite
		//    element matrices, in which the right hand side vectors are also
		//    stored as unassembled finite element vectors, this is the total
		//    number of entries in a single unassembled right hand side vector.
		//
		//    Output, int *NELTVL, the number of finite element matrix entries,
		//    set to 0 in the case of assembled matrices.
		//
		//    Output, char *PTRFMT, the 16 character format for reading pointers.
		//
		//    Output, char *INDFMT, the 16 character format for reading indices.
		//
		//    Output, char *VALFMT, the 20 character format for reading values.
		//
		//    Output, char *RHSFMT, the 20 character format for reading values
		//    of the right hand side.
		//
		//    Output, char *RHSTYP, the 3 character right hand side type.
		//    First character is F for full storage or M for same as matrix.
		//    Second character is G if starting "guess" vectors are supplied.
		//    Third character is X if exact solution vectors are supplied.
		//
		//    Output, int *NRHS, the number of right hand sides.
		//
		//    Output, int *NRHSIX, the number of entries of storage for right
		//    hand side values, in the case where RHSTYP[0] = 'M' and
		//    MXTYPE[2] = 'A'.
		//
		{
		  char *field;
		  char line[255];
		//
		//  Read line 1.
		//
		  input.getline ( line, sizeof ( line ) );
		
		  if ( input.eof() ) {
 			  rb_raise(rb_eArgError, "I/O error reading header line 1.\n");
		    exit ( 1 );
		  }
		
		  *title = s_substring ( line, 1, 72 );
		  s_trim ( *title );
		
		  *key = s_substring ( line, 73, 80 );
		  s_trim ( *key );
		//
		//  Read line 2.
		//
		  input.getline ( line, sizeof ( line ) );
		
		  if ( input.eof() ) {
 			  rb_raise(rb_eArgError, "I/O error reading header line 2.\n");
		    exit ( 1 );
		  }
		
		  field = s_substring ( line,  1, 14 );
		  *totcrd = atoi ( field );
		
		  field = s_substring ( line, 15, 28 );
		  *ptrcrd = atoi ( field );
		
		  field = s_substring ( line, 29, 42 );
		  *indcrd = atoi ( field );
		
		  field = s_substring ( line, 43, 56 );
		  *valcrd = atoi ( field );
		
		  field = s_substring ( line, 57, 70 );
		  *rhscrd = atoi ( field );
		//
		//  Read line 3.
		//
		  input.getline ( line, sizeof ( line ) );
		  
		  if ( input.eof() ) {
 			  rb_raise(rb_eArgError, "I/O error reading header line 3.\n");
		    exit ( 1 );
		  }
		
		  *mxtype = s_substring ( line, 1, 3 );
		  s_trim ( *mxtype );
		
		  field = s_substring ( line, 15, 28 );
		  *nrow = atoi ( field );
		
		  field = s_substring ( line, 29, 42 );
		  *ncol = atoi ( field );
		
		  field = s_substring ( line, 43, 56 );
		  *nnzero = atoi ( field );
		
		  field = s_substring ( line, 57, 70 );
		  *neltvl = atoi ( field );
		//
		//  Read line 4.
		//
		  input.getline ( line, sizeof ( line ) );
		  
		  if ( input.eof() ) {
 			  rb_raise(rb_eArgError, "I/O error reading header line 4.\n");
		    exit ( 1 );
		  }
		
		  *ptrfmt = s_substring ( line,  1, 16 );
		  s_trim ( *ptrfmt );
		
		  *indfmt = s_substring ( line, 17, 32 );
		  s_trim ( *indfmt );
		
		  *valfmt = s_substring ( line, 33, 52 );
		  s_trim ( *valfmt );
		
		  *rhsfmt = s_substring ( line, 53, 72 );
		  s_trim ( *rhsfmt );
		//
		//  Read line 5.
		//
		  if ( 0 < rhscrd ) {
		
		    input.getline ( line, sizeof ( line ) );
		  
		    if ( input.eof() ) {
 			  	rb_raise(rb_eArgError, "I/O error reading header line 5.\n");
		      exit ( 1 );
		    }
		
		    *rhstyp = s_substring ( line, 1, 3 );
		    s_trim ( *rhstyp );
		
		    field = s_substring ( line, 15, 28 );
		    *nrhs = atoi ( field );
		
		    field = s_substring ( line, 29, 42 );
		    *nrhsix = atoi ( field );
		  }
		  else {
		    *rhstyp = NULL;
		    *nrhs = 0;
		    *nrhsix = 0;
		  }
		
		  return;
		}

		char *s_substring ( char *s, int a, int b )

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    S_SUBSTRING returns a substring of a given string.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    10 April 2004
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, char *S, a pointer to a string.
		//
		//    Input, int A, B, the indices of the first and last character of S to copy.
		//    These are 1-based indices!  B should be 
		//
		//    Output, char *S_SUBSTRING, a pointer to the substring.
		//
		{
		  int i;
		  int j;
		  char *t;

		  t = new char[b+2-a];

		  j = 0;
		  for ( i = a; i <= b; i++ ) {
		    t[j] = s[i-1];
		    j = j + 1;
		  }
		  t[j] = '\0';

		  return t;
		}

		void s_trim ( char *s )

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    S_TRIM promotes the final null forward through trailing blanks.
		//
		//  Discussion:
		//
		//    What we're trying to say is that we reposition the null character
		//    so that trailing blanks are no longer visible.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    10 April 2004
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input/output, char *S, the string to be trimmed.
		//
		{
		  char c;
		  int n;
		  char *t;

		  n = strlen ( s );
		  t = s + strlen ( s ) - 1;

		  while ( 0 < n )  {
		    if ( *t != ' ' ) {
		      return;
		    }
		    c      = *t;
		    *t     = *(t+1);
		    *(t+1) = c;
		    t--;
		    n--;
		  }

		  return;
		}

		void hb_structure_read ( std::ifstream &input, int ncol, char *mxtype, int nnzero, 
		  int neltvl, int ptrcrd, char *ptrfmt, int indcrd, char *indfmt, 
		  int colptr[], int rowind[] )

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    HB_STRUCTURE_READ reads the structure of an HB matrix.
		//
		//  Discussion:
		//
		//    The user should already have opened the file, and positioned it
		//    to just after the header records.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    09 April 2004
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Reference:
		//
		//    Iain Duff, Roger Grimes, John Lewis,
		//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
		//    October 1992.
		//
		//  Parameters:
		//
		//    Input, ifstream &INPUT, the unit from which data is read.
		//
		//    Input, int NCOL, the number of columns.
		//
		//    Input, char *MXTYPE, the 3 character matrix type.
		//    First character is R for Real, C for complex, P for pattern only.
		//    Second character is S for symmetric, U for unsymmetric, H for
		//      Hermitian, Z for skew symmetric, R for rectangular.
		//    Third character is A for assembled and E for unassembled
		//      finite element matrices.
		//
		//    Input, int NNZERO.  In the case of assembled sparse matrices,
		//    this is the number of nonzeroes.  In the case of unassembled finite
		//    element matrices, in which the right hand side vectors are also
		//    stored as unassembled finite element vectors, this is the total
		//    number of entries in a single unassembled right hand side vector.
		//
		//    Input, int NELTVL, the number of finite element matrix entries,
		//    set to 0 in the case of assembled matrices.
		//
		//    Input, int PTRCRD, the number of input lines for pointers.
		//
		//    Input, char *PTRFMT, the 16 character format for reading pointers.
		//
		//    Input, int INDCRD, the number of input lines for indices.
		//
		//    Input, char *INDFMT, the 16 character format for reading indices.
		//
		//    Output, int COLPTR[NCOL+1], COLPTR[I-1] points to the location of
		//    the first entry of column I in the sparse matrix structure.
		//
		//    If MXTYPE[2] == 'A':
		//
		//      Output, int ROWIND[NNZERO], the row index of each item.
		//
		//    If MXTYPE[2] == 'F':
		//
		//      Output, int ROWIND[NELTVL], the row index of each item.
		//
		{
		  char code;
		  int i;
		  int j;
		  int jhi;
		  int jlo;
		  int khi;
		  int klo;
		  char line[255];
		  int line_num;
		  int m;
		  int number;
		  int r;
		  char *s;
		  int w;

		  s_to_format ( ptrfmt, &r, &code, &w, &m );

		  if ( mxtype[2] == 'A' ) {
		    line_num = 1 + ( ( ncol + 1 ) - 1 ) / r;
		  }
		  else {
		    line_num = 1 + ( ( ncol     ) - 1 ) / r;
		  }

		  jhi = 0;
		  for ( i = 1; i <= line_num; i++ ) {
		    input.getline ( line, sizeof ( line ) );
		    jlo = jhi + 1;
		    if ( mxtype[2] == 'A' ) {
		      jhi = i4_min ( jlo + r - 1, ncol + 1 );
		    }
		    else {
		      jhi = i4_min ( jlo + r - 1, ncol     );
		    }
		    khi = 0;
		    for ( j = jlo; j <= jhi; j++ ) {
		      klo = khi + 1;
		      khi = i4_min ( klo + w - 1, strlen ( line ) );
		      s = s_substring ( line, klo, khi );
		      colptr[j-1] = atoi ( s );
		    }
		  }

		  if ( mxtype[2] == 'A' ) {
		    s_to_format ( indfmt, &r, &code, &w, &m );

		    line_num = 1 + ( nnzero - 1 ) / r;

		    jhi = 0;
		    for ( i = 1; i <= line_num; i++ ) {
		      input.getline ( line, sizeof ( line ) );
		      jlo = jhi + 1;
		      jhi = i4_min ( jlo + r - 1, nnzero );
		 
		      khi = 0;
		      for ( j = jlo; j <= jhi; j++ ) { 
		        klo = khi + 1;
		        khi = i4_min ( klo + w - 1, strlen ( line ) );
		        s = s_substring ( line, klo, khi );
		        rowind[j-1] = atoi ( s );
		      }
		    }
		  }
		  else if ( mxtype[2] == 'E' ) {
		    s_to_format ( indfmt, &r, &code, &w, &m );

		    number = colptr[ncol-1] - colptr[0];
		    line_num = 1 + ( number - 1 ) / r;

		    jhi = 0;
		    for ( i = 1; i <= line_num; i++ ) {
		      input.getline ( line, sizeof ( line ) );
		      jlo = jhi + 1;
		      jhi = i4_min ( jlo + r - 1, number );
		 
		      khi = 0;
		      for ( j = jlo; j <= jhi; j++ ) { 
		        klo = khi + 1;
		        khi = i4_min ( klo + w - 1, strlen ( line ) );
		        s = s_substring ( line, klo, khi );
		        rowind[j-1] = atoi ( s );
		      }
		    }
		  }
		  else {
		    rb_raise(rb_eStandardError, "Illegal value of MXTYPE character 3.\n" );
		    exit ( 1 );
		  }

		  return;
		}

		void hb_values_read ( std::ifstream &input, int valcrd, char *mxtype, int nnzero,
		  int neltvl, char *valfmt, double values[] )

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    HB_VALUES_READ reads the values of an HB matrix.
		//
		//  Discussion:
		//
		//    The user should already have opened the file, and positioned it
		//    to just after the header and structure records.
		//
		//    Values are contained in an HB file if the VALCRD parameter
		//    is nonzero.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    21 January 2014
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Reference:
		//
		//    Iain Duff, Roger Grimes, John Lewis,
		//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
		//    October 1992.
		//
		//  Parameters:
		//
		//    Input, ifstream &INPUT, the unit from which data is read.
		//
		//    Input, int VALCRD, the number of input lines for numerical values.
		//
		//    Input, char *MXTYPE, the 3 character matrix type.
		//    First character is R for Real, C for complex, P for pattern only.
		//    Second character is S for symmetric, U for unsymmetric, H for
		//      Hermitian, Z for skew symmetric, R for rectangular.
		//    Third character is A for assembled and E for unassembled
		//      finite element matrices.
		//
		//    Input, int NNZERO.  In the case of assembled sparse matrices,
		//    this is the number of nonzeroes.  In the case of unassembled finite
		//    element matrices, in which the right hand side vectors are also
		//    stored as unassembled finite element vectors, this is the total
		//    number of entries in a single unassembled right hand side vector.
		//
		//    Input, int NELTVL, the number of finite element matrix entries,
		//    set to 0 in the case of assembled matrices.
		//
		//    Input, char *VALFMT, the 20 character format for reading values.
		//
		//    If MXTYPE[2] == 'A':
		//
		//      Output, double VALUES[NNZERO], the nonzero values of the matrix.
		//
		//    If MXTYPE[2] == 'E':
		//
		//      Output, double VALUES[NELTVL], the nonzero values of the matrix.
		//
		{
		  char code;
		  int i;
		  int j;
		  int jhi;
		  int jlo;
		  int khi;
		  int klo;
		  char line[255];
		  int line_num;
		  int m;
		  int r;
		  char *s;
		  int w;

		  s_to_format ( valfmt, &r, &code, &w, &m );

		//
		//  Read the matrix values.
		//    case "A" = assembled;
		//    case "E" = unassembled finite element matrices.
		//
		  if ( 0 < valcrd ) {
		    if ( mxtype[2] == 'A' ) {
		      line_num = 1 + ( nnzero - 1 ) / r;
		    }
		    else if ( mxtype[2] == 'E' ) {
		      line_num = 1 + ( neltvl - 1 ) / r;
		    }
		    else {
		    	rb_raise(rb_eStandardError, "Illegal value of MXTYPE character 3.\n" );
		      exit ( 1 );
		    }

		    jhi = 0;
		    for ( i = 1; i <= line_num; i++ ) {

		      input.getline ( line, sizeof ( line ) );
		      jlo = jhi + 1;
		      if ( mxtype[2] == 'A' ) {
		        jhi = i4_min ( jlo + r - 1, nnzero );
		      }
		      else {
		        jhi = i4_min ( jlo + r - 1, neltvl );
		      }
		      khi = 0;
		      for ( j = jlo; j <= jhi; j++ ) {

		        klo = khi + 1;
		        khi = i4_min ( klo + w - 1, strlen ( line ) );
		        s = s_substring ( line, klo, khi );
		        values[j-1] = atof ( s );
		      }
		    }
		  }

		  return;
		}

		void s_to_format ( char *s, int *r, char *code, int *w, int *m )

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    S_TO_FORMAT reads a FORTRAN format from a string.
		//
		//  Discussion:
		//
		//    This routine will read as many characters as possible until it reaches
		//    the end of the string, or encounters a character which cannot be
		//    part of the format.  This routine is limited in its ability to
		//    recognize FORTRAN formats.  In particular, we are only expecting
		//    a single format specification, and cannot handle extra features
		//    such as 'ES' and 'EN' codes, '5X' spacing, and so on.
		//
		//    Legal input is:
		//
		//       0 nothing
		//       1 blanks
		//       2 optional '('
		//       3 blanks
		//       4 optional repeat factor R
		//       5 blanks
		//       6 CODE ( 'A', 'B', 'E', 'F', 'G', 'I', 'L', 'O', 'Z', '*' )
		//       7 blanks
		//       8 width W
		//       9 optional decimal point
		//      10 optional mantissa M
		//      11 blanks
		//      12 optional ')'
		//      13 blanks
		//
		//  Example:
		//
		//    S                 R   CODE   W    M
		//
		//    'I12              1   I      12   0
		//    'E8.0'            1   E       8   0
		//    'F10.5'           1   F      10   5
		//    '2G14.6'          2   G      14   6
		//    '*'               1   *      -1  -1
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    10 April 2004
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, char *S, the string containing the
		//    data to be read.  Reading will begin at position 1 and
		//    terminate at the end of the string, or when no more
		//    characters can be read.
		//
		//    Output, int *R, the repetition factor, which defaults to 1.
		//
		//    Output, char *CODE, the format code.
		//
		//    Output, int *W, the field width.
		//
		//    Output, int *M, the mantissa width.
		//
		{
		  char c;
		  int d;
		  bool debug = true;
		  int LEFT = 1;
		  int paren_sum;
		  int pos;
		  int RIGHT = -1;
		  int s_length;
		  int state;

		  state = 0;
		  paren_sum = 0;
		  pos = 0;
		  s_length = s_len_trim ( s );

		  *r = 0;
		  *w = 0;
		  *code = '?';
		  *m = 0;

		  while ( pos < s_length ) {
		    c = s[pos];
		    pos = pos + 1;
		//
		//  BLANK character:
		//
		    if ( c == ' ' ) {
		      if ( state == 4 ) {
		        state = 5;
		      }
		      else if ( state == 6 ) {
		        state = 7;
		      }
		      else if ( state == 10 ) {
		        state = 11;
		      }
		      else if ( state == 12 ) {
		        state = 13;
		      }
		    }
		//
		//  LEFT PAREN
		//
		    else if ( c == '(' ) {
		      if ( state < 2 ) {
		        paren_sum = paren_sum + LEFT;
		      }
		      else {
		        if ( debug ) {
		          rb_raise(rb_eStandardError, "S_TO_FORMAT - Fatal error!\nCurrent state = %d\nInput character = %c", state, c);
		        }
		        state = -1;
		        break;
		      }
		    }
		//
		//  DIGIT (R, F, or W)
		//
		    else if ( ch_is_digit ( c ) ) {
		      if ( state <= 3 ) {
		        state = 4;
		        *r = ch_to_digit ( c );
		      }
		      else if ( state == 4 ) {
		        d = ch_to_digit ( c );
		        *r = 10 * (*r) + d;
		      }
		      else if ( state == 6 || state == 7 ) {
		        if ( *code == '*' ) {
		          if ( debug ) {
		          	rb_raise(rb_eStandardError, "S_TO_FORMAT - Fatal error!\nCurrent state = %d \nCurrent code = %s\nInput character = %c", state, code, c);
		          }
		          state = -1;
		          break;
		        }
		        state = 8;
		        *w = ch_to_digit ( c );
		      }
		      else if ( state == 8 ) {
		        d = ch_to_digit ( c );
		        *w = 10 * (*w) + d;
		      }
		      else if ( state == 9 ) {		     
		        state = 10;
		        *m = ch_to_digit ( c );
		      }
		      else if ( state == 10 ) {
		        d = ch_to_digit ( c );
		        *m = 10 * (*m) + d;
		      }
		      else {
		        if ( debug ) {
		        	rb_raise(rb_eStandardError, "S_TO_FORMAT - Fatal error!\nCurrent state = %d\nInput character = %c", state, c);
		        }
		        state = -1;
		        break;
		      }
		    }
		//
		//  DECIMAL POINT
		//
		    else if ( c == '.' ) {
		      if ( state == 8 ) {
		        state = 9;
		      }
		      else {
		        if ( debug ) {
		        	rb_raise(rb_eStandardError, "S_TO_FORMAT - Fatal error!\nCurrent state = %d\nInput character = %c", state, c);
		        }
		        state = -1;
		        break;
		      }
		    }
		//
		//  RIGHT PAREN
		//
		    else if ( c == ')' )
		    {
		      paren_sum = paren_sum + RIGHT;

		      if ( paren_sum != 0 )
		      {
		        if ( debug ) {
		        	rb_raise(rb_eStandardError, "S_TO_FORMAT - Fatal error!\nCurrent paren sum = %d\nInput character = %c", paren_sum, c);
		        }
		        state = -1;
		        break;
		      }

		      if ( state == 6 && *code == '*' ) {
		        state = 12;
		      }
		      else if ( 6 <= state ) {
		        state = 12;
		      }
		      else {
		        if ( debug ) {
		        	rb_raise(rb_eStandardError, "S_TO_FORMAT - Fatal error!\nCurrent state = %d\nInput character = %c", state, c);
		        }
		        state = -1;
		        break;
		      }
		    }
		//
		//  Code
		//
		    else if ( ch_is_format_code ( c ) ) {
		      if ( state < 6 ) {
		        state = 6;
		        *code = c;
		      }
		      else {
		        if ( debug ) {
		        	rb_raise(rb_eStandardError, "S_TO_FORMAT - Fatal error!\nCurrent state = %d\nInput character = %c", state, c);
		        }
		        state = -1;
		        break;
		      }
		    }
		//
		//  Unexpected character
		//
		    else {
		      if ( debug ) {
		      	rb_raise(rb_eStandardError, "S_TO_FORMAT - Fatal error!\nCurrent state = %d\nInput character = %c", state, c);
		      }
		      state = -1;
		      break;
		    }
		  }

		  if ( paren_sum != 0 ) {
		  	rb_raise(rb_eStandardError, "S_TO_FORMAT - Fatal error!\nParentheses mismatch.\n");
		    exit ( 1 );
		  }

		  if ( state < 0 ) {
		  	rb_raise(rb_eStandardError, "S_TO_FORMAT - Fatal error!\nParsing error.\n");
		    exit ( 1 );
		  }

		  if ( *r == 0 ) {
		    *r = 1;
		  }

		  return;
		}

		int i4_min ( int i1, int i2 )

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    I4_MIN returns the smaller of two I4's.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    13 October 1998
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, int I1, I2, two integers to be compared.
		//
		//    Output, int I4_MIN, the smaller of I1 and I2.
		//
		//
		{
		  if ( i1 < i2 ) {
		    return i1;
		  }
		  else {
		    return i2;
		  }
		}

		int s_len_trim ( char *s )

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    S_LEN_TRIM returns the length of a string to the last nonblank.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    26 April 2003
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, char *S, a pointer to a string.
		//
		//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
		//    If S_LEN_TRIM is 0, then the string is entirely blank.
		//
		{
		  int n;
		  char* t;

		  n = strlen ( s );
		  t = s + strlen ( s ) - 1;

		  while ( 0 < n ) {
		    if ( *t != ' ' ) {
		      return n;
		    }
		    t--;
		    n--;
		  }

		  return n;
		}

		int ch_to_digit ( char c )

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    CH_TO_DIGIT returns the integer value of a base 10 digit.
		//
		//  Example:
		//
		//     C   DIGIT
		//    ---  -----
		//    '0'    0
		//    '1'    1
		//    ...  ...
		//    '9'    9
		//    ' '    0
		//    'X'   -1
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    13 June 2003
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, char C, the decimal digit, '0' through '9' or blank are legal.
		//
		//    Output, int CH_TO_DIGIT, the corresponding integer value.  If C was
		//    'illegal', then DIGIT is -1.
		//
		{
		  int digit;

		  if ( '0' <= c && c <= '9' ) {
		    digit = c - '0';
		  }
		  else if ( c == ' ' ) {
		    digit = 0;
		  }
		  else {
		    digit = -1;
		  }

		  return digit;
		}

		bool ch_is_digit ( char c )

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    CH_IS_DIGIT returns TRUE if a character is a decimal digit.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    05 December 2003
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, char C, the character to be analyzed.
		//
		//    Output, bool CH_IS_DIGIT, is TRUE if C is a digit.
		//
		{
		  if ( '0' <= c && c <= '9' ) {
		    return true;
		  }
		  else {
		    return false;
		  }
		}

		bool ch_is_format_code ( char c )

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    CH_IS_FORMAT_CODE returns TRUE if a character is a FORTRAN format code.
		//
		//  Discussion:
		//
		//    The format codes accepted here are not the only legal format
		//    codes in FORTRAN90.  However, they are more than sufficient
		//    for my needs!
		//
		//  Table:
		//
		//    A  Character
		//    B  Binary digits
		//    D  Real number, exponential representation
		//    E  Real number, exponential representation
		//    F  Real number, fixed point
		//    G  General format
		//    I  Integer
		//    L  Logical variable
		//    O  Octal digits
		//    Z  Hexadecimal digits
		//    *  Free format
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    21 November 2003
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, char C, the character to be analyzed.
		//
		//    Output, bool CH_IS_FORMAT_CODE, is TRUE if C is a FORTRAN format code.
		//
		{
		  if ( ch_eqi ( c, 'A' ) )  {
		    return true;
		  }
		  else if ( ch_eqi ( c, 'B' ) ) {
		    return true;
		  }
		  else if ( ch_eqi ( c, 'D' ) ) {
		    return true;
		  }
		  else if ( ch_eqi ( c, 'E' ) ) {
		    return true;
		  }
		  else if ( ch_eqi ( c, 'F' ) ) {
		    return true;
		  }
		  else if ( ch_eqi ( c, 'G' ) ) {
		    return true;
		  }
		  else if ( ch_eqi ( c, 'I' ) ) {
		    return true;
		  }
		  else if ( ch_eqi ( c, 'L' ) ) {
		    return true;
		  }
		  else if ( ch_eqi ( c, 'O' ) ) {
		    return true;
		  }
		  else if ( ch_eqi ( c, 'Z' ) ) {
		    return true;
		  }
		  else if ( c == '*' ) {
		    return true;
		  }
		  else {
		    return false;
		  }
		}

		bool ch_eqi ( char c1, char c2 )

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    CH_EQI is true if two characters are equal, disregarding case.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    13 June 2003
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, char C1, C2, the characters to compare.
		//
		//    Output, bool CH_EQI, is true if the two characters are equal,
		//    disregarding case.
		//
		{
		  if ( 97 <= c1 && c1 <= 122 ) {
		    c1 = c1 - 32;
		  } 
		  if ( 97 <= c2 && c2 <= 122 ) {
		    c2 = c2 - 32;
		  }     

		  return ( c1 == c2 );
		}
	}
}