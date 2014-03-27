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
// == getf2.h
//
// LAPACK getf2 function in native C++.
//

#ifndef GETF2_H
#define GETF2_H

namespace nm { namespace math {

template <typename DType>
inline int getf2(const int m, const int n, DType* a, const int lda, int *ipiv) {

  /* Function Body */
  if (m < 0)                      return -1; // error
  else if (n < 0)                 return -2; // error
  else if (lda < std::max(1,m))   return -4; // error


  if (m == 0 || n == 0)     return 0;   /* Quick return if possible */

  for (size_t j = 0; j < std::min(m,n); ++j) { // changed

    /* Find pivot and test for singularity. */

    int jp = j - 1 + idamax<DType>(m-j+1, &a[j + j * lda], 1);

    ipiv[j] = jp;


    if (a[jp + j*lda] != 0) {

      /* Apply the interchange to columns 1:N. */
      // (Don't swap two columns that are the same.)
      if (jp != j) swap<DType>(n, &a[j], lda, &a[jp], lda);

      /* Compute elements J+1:M of J-th column. */

	    if (j < m-1) {
        if (std::abs(a[j+j*lda]) >= std::numeric_limits<DType>::min()) {
          scal<DType>(m-j, 1.0 / a[j+j*lda], &a[j+1+j*lda], 1);
		    } else {
		      for (size_t i = 0; i < m-j; ++i) { // changed
			      a[j+i+j*lda] /= a[j+j*lda];
		      }
		    }
	    }

    } else { // singular matrix
      return j; // U(j,j) is exactly zero, div by zero if answer is used to solve a system of equations.
    }

    if (j < std::min(m,n)-1) /*           Update trailing submatrix. */
      ger<DType>(m-j, n-j, -1.0, &a[j+1+j*lda], 1, &a[j+(j+1)*lda], lda, &a[j+1+(j+1)*lda], lda);

  }
  return 0;
} /* dgetf2_ */


}} // end of namespace nm::math

#endif // GETF2