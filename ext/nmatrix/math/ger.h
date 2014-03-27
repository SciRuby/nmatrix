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
// == ger.h
//
// BLAS level-2 ger function in native C++.
//

#ifndef GER_H
#define GER_H

namespace nm { namespace math {

template <typename DType>
inline int ger(int m, int n, DType alpha, DType* x, int incx, DType* y, int incy, DType* a, int lda) {

  // FIXME: Call BLAS ger if available

  if (m < 0) {
    return 1;
  } else if (n < 0) {
    return 2;
  } else if (incx == 0) {
    return 5;
  } else if (incy == 0) {
    return 7;
  } else if (lda < std::max(1,m)) {
    return 9;
  }

  if (m == 0 || n == 0 || alpha == 0) return 0; /* Quick return if possible. */

  /*     Start the operations. In this version the elements of A are */
  /*     accessed sequentially with one pass through A. */

  // FIXME: These have been unrolled in a way that the compiler can handle. Collapse into a single case, or optimize
  // FIXME: in a more modern way.

  int jy = incy > 0 ? 0 :  -(n-1) * incy;

  if (incx == 1) {

	  for (size_t j = 0; j < n; ++j, jy += incy) {
	    if (y[jy] != 0) {
		    DType temp = alpha * y[jy];
		    for (size_t i = 0; i < m; ++i) {
		      a[i + j * lda] += x[i] * temp;
		    }
	    }
	  }

  } else {

    int kx = incx > 0 ? 0 : -(m-1) * incx;

	  for (size_t j = 0; j < n; ++j, jy += incy) {
	    if (y[jy] != 0) {
    		DType temp = alpha * y[jy];

    		for (size_t i = 0, ix = kx; i < m; ++i, ix += incx) {
          a[i + j * lda] += x[ix] * temp;
    		}
	    }
	  }

  }

  return 0;

/*     End of DGER  . */

} /* dger_ */

}} // end nm::math

#endif // GER_H