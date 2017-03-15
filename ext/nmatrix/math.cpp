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
// SciRuby is Copyright (c) 2010 - present, Ruby Science Foundation
// NMatrix is Copyright (c) 2012 - present, John Woods and the Ruby Science Foundation
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
// == math.cpp
//
// Ruby-exposed CBLAS and LAPACK functions that are available without
// an external library.
//
// === Procedure for adding CBLAS functions to math.cpp/math.h:
//
// This procedure is written as if for a fictional function with double
// version dbacon, which we'll say is from CBLAS.
//
// 1. Write a default templated version which probably returns a boolean.
//    Call it bacon, and put it in math.h.
//
//    template <typename DType>
//    bool bacon(const CBLAS_TRANSPOSE trans, const int M, const int N, DType* A, ...) {
//      rb_raise(rb_eNotImpError, "only implemented for ATLAS types (float32, float64, complex64, complex128)");
//    }
//
//    Make sure this is in namespace nm::math
//
// 2. In math.cpp, add a templated inline static version of the function which takes
//    only void* pointers and uses static_cast to convert them to the
//    proper dtype. This should also be in namespace nm::math
//
//    This function may also need to switch m and n if these arguments are given.
//
//    For an example, see cblas_gemm. This function should do nothing other than cast
//    appropriately. If cblas_dbacon, cblas_sbacon, cblas_cbacon, and cblas_zbacon
//    all take void* only, and no other pointers that vary between functions, you can skip
//    this particular step -- as we can call them directly using a custom function pointer
//    array (same function signature!).
//
//    This version of the function will be the one exposed through NMatrix::BLAS. We
//    want it to be as close to the actual BLAS version of the function as possible,
//    and with as few checks as possible.
//
//    You will probably need a forward declaration in the extern "C" block.
//
//    Note: In that case, the function you wrote in Step 1 should also take exactly the
//    same arguments as cblas_xbacon. Otherwise Bad Things will happen.
//
// 3. In cblas_templates_core.h, add a default template like in step 1 (which will just
//    call nm::math::bacon()) and also
//    inline specialized versions of bacon for the different BLAS types.
//    This will allow both nmatrix-atlas and nmatrix-lapacke to use the optimized version
//    of bacon from whatever external library is available, as well as the internal version
//    if an external version is not available. These functions will end up in a namsespace
//    like nm::math::atlas, but don't explicitly put them in a namespace, they will get
//    put in the appropriate namespace when cblas_templates_core.h is included.
//
//    template <typename DType>
//    inline bool bacon(const CBLAS_TRANSPOSE trans, const int M, const int N, DType* A, ...) {
//      nm::math::bacon(trans, M, N, A, ...);
//    }
//
//    template <>
//    inline bool bacon(const CBLAS_TRANSPOSE trans, const int M, const int N, float* A, ...) {
//      cblas_sbacon(trans, M, N, A, ...);
//      return true;
//    }
//
//    Note that you should do everything in your power here to parse any return values
//    cblas_sbacon may give you. We're not trying very hard in this example, but you might
//    look at getrf to see how it might be done.
//
// 4. Write the C function nm_cblas_bacon, which is what Ruby will call. Use the example
//    of nm_cblas_gemm below. Also you must add a similar function in math_atlas.cpp
//    and math_lapacke.cpp
//
// 5. Expose the function in nm_math_init_blas(), in math.cpp:
//
//    rb_define_singleton_method(cNMatrix_Internal_BLAS, "cblas_bacon", (METHOD)nm_cblas_bacon, 5);
//
//    Do something similar in math_atlas.cpp and math_lapacke.cpp to add the function
//    to the plugin gems.
//
//    Here, we're telling Ruby that nm_cblas_bacon takes five arguments as a Ruby function.
//
// 6. In blas.rb, write a bacon function which accesses cblas_bacon, but does all the
//    sanity checks we left out in step 2.
//
// 7. Write tests for NMatrix::BLAS::bacon, confirming that it works for the ATLAS dtypes.
//
// 8. After you get it working properly with CBLAS, download dbacon.f from NETLIB, and use
//    f2c to convert it to C. Clean it up so it's readable. Remove the extra indices -- f2c
//    inserts a lot of unnecessary stuff.
//
//    Copy and paste the output into the default templated function you wrote in Step 1.
//    Fix it so it works as a template instead of just for doubles.
//
//    Because of step 3, this will automatically also work for the nmatrix-atlas
//    and nmatrix-lapacke implementations.
//
// 9. Write tests to confirm that it works for all data types.
//
// 10. See about adding a Ruby-like interface, such as matrix_matrix_multiply for cblas_gemm,
//    or matrix_vector_multiply for cblas_gemv. This step is not mandatory.
//
// 11. Pull request!

/*
 * Project Includes
 */


#include <ruby.h>
#include <algorithm>
#include <limits>
#include <cmath>

#include "math/cblas_enums.h"

#include "data/data.h"
#include "math/magnitude.h"
#include "math/imax.h"
#include "math/scal.h"
#include "math/laswp.h"
#include "math/trsm.h"
#include "math/gemm.h"
#include "math/gemv.h"
#include "math/asum.h"
#include "math/nrm2.h"
#include "math/getrf.h"
#include "math/getrs.h"
#include "math/rot.h"
#include "math/rotg.h"
#include "math/math.h"
#include "math/util.h"
#include "storage/dense/dense.h"

#include "nmatrix.h"
#include "ruby_constants.h"

/*
 * Forward Declarations
 */

extern "C" {
  /* BLAS Level 1. */
  static VALUE nm_cblas_scal(VALUE self, VALUE n, VALUE scale, VALUE vector, VALUE incx);
  static VALUE nm_cblas_nrm2(VALUE self, VALUE n, VALUE x, VALUE incx);
  static VALUE nm_cblas_asum(VALUE self, VALUE n, VALUE x, VALUE incx);
  static VALUE nm_cblas_rot(VALUE self, VALUE n, VALUE x, VALUE incx, VALUE y, VALUE incy, VALUE c, VALUE s);
  static VALUE nm_cblas_rotg(VALUE self, VALUE ab);
  static VALUE nm_cblas_imax(VALUE self, VALUE n, VALUE x, VALUE incx);

  /* BLAS Level 2. */
  static VALUE nm_cblas_gemv(VALUE self, VALUE trans_a, VALUE m, VALUE n, VALUE vAlpha, VALUE a, VALUE lda,
                             VALUE x, VALUE incx, VALUE vBeta, VALUE y, VALUE incy);

  /* BLAS Level 3. */
  static VALUE nm_cblas_gemm(VALUE self, VALUE order, VALUE trans_a, VALUE trans_b, VALUE m, VALUE n, VALUE k, VALUE vAlpha,
                             VALUE a, VALUE lda, VALUE b, VALUE ldb, VALUE vBeta, VALUE c, VALUE ldc);
  static VALUE nm_cblas_trsm(VALUE self, VALUE order, VALUE side, VALUE uplo, VALUE trans_a, VALUE diag, VALUE m, VALUE n,
                             VALUE vAlpha, VALUE a, VALUE lda, VALUE b, VALUE ldb);

  /* LAPACK. */
  static VALUE nm_has_clapack(VALUE self);
  static VALUE nm_clapack_getrf(VALUE self, VALUE order, VALUE m, VALUE n, VALUE a, VALUE lda);
  static VALUE nm_clapack_getrs(VALUE self, VALUE order, VALUE trans, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE ipiv, VALUE b, VALUE ldb);
  static VALUE nm_clapack_laswp(VALUE self, VALUE n, VALUE a, VALUE lda, VALUE k1, VALUE k2, VALUE ipiv, VALUE incx);
} // end of extern "C" block

////////////////////
// Math Functions //
////////////////////

namespace nm {
  namespace math {

    /*
     * Calculate the determinant for a dense matrix (A [elements]) of size 2 or 3. Return the result.
     */
    template <typename DType>
    void det_exact_from_dense(const int M, const void* A_elements, const int lda, void* result_arg) {
      DType* result  = reinterpret_cast<DType*>(result_arg);
      const DType* A = reinterpret_cast<const DType*>(A_elements);

      typename LongDType<DType>::type x, y;

      if (M == 2) {
        *result = A[0] * A[lda+1] - A[1] * A[lda];
      } else if (M == 3) {
        x = A[lda+1] * A[2*lda+2] - A[lda+2] * A[2*lda+1]; // ei - fh
        y = A[lda] * A[2*lda+2] -   A[lda+2] * A[2*lda];   // fg - di
        x = A[0]*x - A[1]*y ; // a*(ei-fh) - b*(fg-di)

        y = A[lda] * A[2*lda+1] - A[lda+1] * A[2*lda];    // dh - eg
        *result = A[2]*y + x; // c*(dh-eg) + _
      } else if (M < 2) {
        rb_raise(rb_eArgError, "can only calculate exact determinant of a square matrix of size 2 or larger");
      } else {
        rb_raise(rb_eNotImpError, "exact determinant calculation needed for matrices larger than 3x3");
      }
    }

    //we can't do det_exact on byte, because it will want to return a byte (unsigned), but determinants can be negative, even if all elements of the matrix are positive
    template <>
    void det_exact_from_dense<uint8_t>(const int M, const void* A_elements, const int lda, void* result_arg) {
      rb_raise(nm_eDataTypeError, "cannot call det_exact on unsigned type");
    }
    /*
     * Calculate the determinant for a yale matrix (storage) of size 2 or 3. Return the result.
     */
    template <typename DType>
    void det_exact_from_yale(const int M, const YALE_STORAGE* storage, const int lda, void* result_arg) {
      DType* result  = reinterpret_cast<DType*>(result_arg);
      IType* ija = reinterpret_cast<IType *>(storage->ija);
      DType* a = reinterpret_cast<DType*>(storage->a);
      IType col_pos = storage->shape[0] + 1;
      if (M == 2) {
        if (ija[2] - ija[0] == 2) {
          *result = a[0] * a[1] - a[col_pos] * a[col_pos+1];
        }
        else { *result = a[0] * a[1]; }
      } else if (M == 3) {
        DType m[3][3];
        for (int i = 0; i < 3; ++i) {
          m[i][i] = a[i];
          switch(ija[i+1] - ija[i]) {
          case 2:
            m[i][ija[col_pos]] = a[col_pos];
            m[i][ija[col_pos+1]] = a[col_pos+1];
            col_pos += 2;
            break;
          case 1:
            m[i][(i+1)%3] = m[i][(i+2)%3] = 0;
            m[i][ija[col_pos]] = a[col_pos];
            ++col_pos;
            break;
          case 0:
            m[i][(i+1)%3] = m[i][(i+2)%3] = 0;
            break;
          default:
            rb_raise(rb_eArgError, "some value in IJA is incorrect!");
          }
        }
        *result =
          m[0][0] * m[1][1] * m[2][2] + m[0][1] * m[1][2] * m[2][0] + m[0][2] * m[1][0] * m[2][1]
        - m[0][0] * m[1][2] * m[2][1] - m[0][1] * m[1][0] * m[2][2] - m[0][2] * m[1][1] * m[2][0];

      } else if (M < 2) {
        rb_raise(rb_eArgError, "can only calculate exact determinant of a square matrix of size 2 or larger");
      } else {
        rb_raise(rb_eNotImpError, "exact determinant calculation needed for matrices larger than 3x3");
      }
    }

    /*
     * Solve a system of linear equations using forward-substution followed by
     * back substution from the LU factorization of the matrix of co-efficients.
     * Replaces x_elements with the result. Works only with non-integer, non-object
     * data types.
     *
     * args - r           -> The number of rows of the matrix.
     *        lu_elements -> Elements of the LU decomposition of the co-efficients
     *                       matrix, as a contiguos array.
     *        b_elements  -> Elements of the the right hand sides, as a contiguous array.
     *        x_elements  -> The array that will contain the results of the computation.
     *        pivot       -> Positions of permuted rows.
     */
    template <typename DType>
    void solve(const int r, const void* lu_elements, const void* b_elements, void* x_elements, const int* pivot) {
      int ii = 0, ip;
      DType sum;

      const DType* matrix = reinterpret_cast<const DType*>(lu_elements);
      const DType* b      = reinterpret_cast<const DType*>(b_elements);
      DType* x            = reinterpret_cast<DType*>(x_elements);

      for (int i = 0; i < r; ++i) { x[i] = b[i]; }
      for (int i = 0; i < r; ++i) { // forward substitution loop
        ip = pivot[i];
        sum = x[ip];
        x[ip] = x[i];

        if (ii != 0) {
          for (int j = ii - 1;j < i; ++j) { sum = sum - matrix[i * r + j] * x[j]; }
        }
        else if (sum != 0.0) {
          ii = i + 1;
        }
        x[i] = sum;
      }

      for (int i = r - 1; i >= 0; --i) { // back substitution loop
        sum = x[i];
        for (int j = i + 1; j < r; j++) { sum = sum - matrix[i * r + j] * x[j]; }
        x[i] = sum/matrix[i * r + i];
      }
    }

    /*
     * Calculates in-place inverse of A_elements. Uses Gauss-Jordan elimination technique.
     * In-place inversion of the matrix saves on memory and time.
     *
     * args - M - Shape of the matrix
     *        a_elements - A duplicate of the original expressed as a contiguos array
     */
    template <typename DType>
    void inverse(const int M, void* a_elements) {
      DType* matrix   = reinterpret_cast<DType*>(a_elements);
      int row_index[M]; // arrays for keeping track of column scrambling
      int col_index[M];

      for (int k = 0;k < M; ++k) {
        typename MagnitudeDType<DType>::type akk;
        akk = magnitude( matrix[k * (M + 1)] ); // diagonal element

        int interchange = k;

        for (int row = k + 1; row < M; ++row) {
          typename MagnitudeDType<DType>::type big;
          big = magnitude( matrix[M*row + k] ); // element below the temp pivot

          if ( big > akk ) {
            interchange = row;
            akk = big;
          }
        }

        if (interchange != k) { // check if rows need flipping
          DType temp;

          for (int col = 0; col < M; ++col) {
            NM_SWAP(matrix[interchange*M + col], matrix[k*M + col], temp);
          }
        }

        row_index[k] = interchange;
        col_index[k] = k;

        if (matrix[k * (M + 1)] == (DType)(0)) {
          rb_raise(rb_eZeroDivError, "Expected Non-Singular Matrix.");
        }

        DType pivot = matrix[k * (M + 1)];
        matrix[k * (M + 1)] = (DType)(1); // set diagonal as 1 for in-place inversion

        for (int col = 0; col < M; ++col) {
          // divide each element in the kth row with the pivot
          matrix[k*M + col] = matrix[k*M + col] / pivot;
        }

        for (int kk = 0; kk < M; ++kk) { // iterate and reduce all rows
          if (kk == k) continue;

          DType dum = matrix[k + M*kk];
          matrix[k + M*kk] = (DType)(0); // prepare for inplace inversion
          for (int col = 0; col < M; ++col) {
            matrix[M*kk + col] = matrix[M*kk + col] - matrix[M*k + col] * dum;
          }
        }
      }

      // Unscramble columns
      DType temp;

      for (int k = M - 1; k >= 0; --k) {
        if (row_index[k] != col_index[k]) {

          for (int row = 0; row < M; ++row) {
            NM_SWAP(matrix[row * M + row_index[k]], matrix[row * M + col_index[k]],
              temp);
          }
        }
      }
    }

    /*
     * Reduce a square matrix to hessenberg form with householder transforms
     *
     * == Arguments
     *
     * nrows - The number of rows present in matrix a.
     * a_elements - Elements of the matrix to be reduced in 1D array form.
     *
     * == References
     *
     * http://www.mymathlib.com/c_source/matrices/eigen/hessenberg_orthog.c
     * This code has been included by permission of the author.
     */
    template <typename DType>
    void hessenberg(const int nrows, void* a_elements) {
      DType* a = reinterpret_cast<DType*>(a_elements);
      DType* u = new DType[nrows]; // auxillary storage for the chosen vector
      DType sum_of_squares, *p_row, *psubdiag, *p_a, scale, innerproduct;
      int i, k, col;

      // For each column use a Householder transformation to zero all entries
      // below the subdiagonal.
      for (psubdiag = a + nrows, col = 0; col < nrows - 2; psubdiag += nrows + 1,
        col++) {
        // Calculate the signed square root of the sum of squares of the
        // elements below the diagonal.

        for (p_a = psubdiag, sum_of_squares = 0.0, i = col + 1; i < nrows;
          p_a += nrows, i++) {
          sum_of_squares += *p_a * *p_a;
        }
        if (sum_of_squares == 0.0) { continue; }
        sum_of_squares = std::sqrt(sum_of_squares);

        if ( *psubdiag >= 0.0 ) { sum_of_squares = -sum_of_squares; }

        // Calculate the Householder transformation Q = I - 2uu'/u'u.
        u[col + 1] = *psubdiag - sum_of_squares;
        *psubdiag = sum_of_squares;

        for (p_a = psubdiag + nrows, i = col + 2; i < nrows; p_a += nrows, i++) {
          u[i] = *p_a;
          *p_a = 0.0;
        }

        // Premultiply A by Q
        scale = -1.0 / (sum_of_squares * u[col+1]);
        for (p_row = psubdiag - col, i = col + 1; i < nrows; i++) {
          p_a = a + nrows * (col + 1) + i;
          for (innerproduct = 0.0, k = col + 1; k < nrows; p_a += nrows, k++) {
            innerproduct += u[k] * *p_a;
          }
          innerproduct *= scale;
          for (p_a = p_row + i, k = col + 1; k < nrows; p_a += nrows, k++) {
            *p_a -= u[k] * innerproduct;
          }
        }

        // Postmultiply QA by Q
        for (p_row = a, i = 0; i < nrows; p_row += nrows, i++) {
          for (innerproduct = 0.0, k = col + 1; k < nrows; k++) {
            innerproduct += u[k] * *(p_row + k);
          }
          innerproduct *= scale;

          for (k = col + 1; k < nrows; k++) {
            *(p_row + k) -= u[k] * innerproduct;
          }
        }
      }

      delete[] u;
    }

    void raise_not_invertible_error() {
        rb_raise(nm_eNotInvertibleError,
            "matrix must have non-zero determinant to be invertible (not getting this error does not mean matrix is invertible if you're dealing with floating points)");
    }

    /*
     * Calculate the exact inverse for a dense matrix (A [elements]) of size 2 or 3. Places the result in B_elements.
     */
    template <typename DType>
    void inverse_exact_from_dense(const int M, const void* A_elements,
        const int lda, void* B_elements, const int ldb) {

      const DType* A = reinterpret_cast<const DType*>(A_elements);
      DType* B       = reinterpret_cast<DType*>(B_elements);

      if (M == 2) {
        DType det = A[0] * A[lda+1] - A[1] * A[lda];
        if (det == 0) { raise_not_invertible_error(); }
        B[0] = A[lda+1] / det;
        B[1] = -A[1] / det;
        B[ldb] = -A[lda] / det;
        B[ldb+1] = A[0] / det;

      } else if (M == 3) {
        // Calculate the exact determinant.
        DType det;
        det_exact_from_dense<DType>(M, A_elements, lda, reinterpret_cast<void*>(&det));
        if (det == 0) { raise_not_invertible_error(); }

        B[0]      = (  A[lda+1] * A[2*lda+2] - A[lda+2] * A[2*lda+1]) / det; // A = ei - fh
        B[1]      = (- A[1]     * A[2*lda+2] + A[2]     * A[2*lda+1]) / det; // D = -bi + ch
        B[2]      = (  A[1]     * A[lda+2]   - A[2]     * A[lda+1])   / det; // G = bf - ce
        B[ldb]    = (- A[lda]   * A[2*lda+2] + A[lda+2] * A[2*lda])   / det; // B = -di + fg
        B[ldb+1]  = (  A[0]     * A[2*lda+2] - A[2]     * A[2*lda])   / det; // E = ai - cg
        B[ldb+2]  = (- A[0]     * A[lda+2]   + A[2]     * A[lda])     / det; // H = -af + cd
        B[2*ldb]  = (  A[lda]   * A[2*lda+1] - A[lda+1] * A[2*lda])   / det; // C = dh - eg
        B[2*ldb+1]= ( -A[0]     * A[2*lda+1] + A[1]     * A[2*lda])   / det; // F = -ah + bg
        B[2*ldb+2]= (  A[0]     * A[lda+1]   - A[1]     * A[lda])     / det; // I = ae - bd
      } else if (M == 1) {
        B[0] = 1 / A[0];
      } else {
        rb_raise(rb_eNotImpError, "exact inverse calculation needed for matrices larger than 3x3");
      }
    }

    template <typename DType>
    void inverse_exact_from_yale(const int M, const YALE_STORAGE* storage,
        const int lda, YALE_STORAGE* inverse, const int ldb) {

      // inverse is a clone of storage
      const DType* a = reinterpret_cast<const DType*>(storage->a);
      const IType* ija = reinterpret_cast<const IType *>(storage->ija);
      DType* b       = reinterpret_cast<DType*>(inverse->a);
      IType* ijb = reinterpret_cast<IType *>(inverse->ija);
      IType col_pos = storage->shape[0] + 1;
      // Calculate the exact determinant.
      DType det;

      if (M == 2) {
        IType ndnz = ija[2] - ija[0];
        if (ndnz == 2) {
          det = a[0] * a[1] - a[col_pos] * a[col_pos+1];
        }
        else { det = a[0] * a[1]; }
        if (det == 0) { raise_not_invertible_error(); }
        b[0] = a[1] / det;
        b[1] = a[0] / det;
        if (ndnz == 2) {
          b[col_pos] = -a[col_pos] / det;
          b[col_pos+1] = -a[col_pos+1] / det;
        }
        else if (ndnz == 1) {
          b[col_pos] = -a[col_pos] / det;
        }

      } else if (M == 3) {
        DType *A = new DType[lda*3];
        for (int i = 0; i < lda; ++i) {
          A[i*3+i] = a[i];
          switch (ija[i+1] - ija[i]) {
          case 2:
            A[i*3 + ija[col_pos]] = a[col_pos];
            A[i*3 + ija[col_pos+1]] = a[col_pos+1];
            col_pos += 2;
            break;
          case 1:
            A[i*3 + (i+1)%3] = A[i*3 + (i+2)%3] = 0;
            A[i*3 + ija[col_pos]] = a[col_pos];
            col_pos += 1;
            break;
          case 0:
            A[i*3 + (i+1)%3] = A[i*3 + (i+2)%3] = 0;
            break;
          default:
            rb_raise(rb_eArgError, "some value in IJA is incorrect!");
          }
        }
        det =
          A[0] * A[lda+1] * A[2*lda+2] + A[1] * A[lda+2] * A[2*lda] + A[2] * A[lda] * A[2*lda+1]
        - A[0] * A[lda+2] * A[2*lda+1] - A[1] * A[lda] * A[2*lda+2] - A[2] * A[lda+1] * A[2*lda];
        if (det == 0) { raise_not_invertible_error(); }

        DType *B = new DType[3*ldb];
        B[0]      = (  A[lda+1] * A[2*lda+2] - A[lda+2] * A[2*lda+1]) / det; // A = ei - fh
        B[1]      = (- A[1]     * A[2*lda+2] + A[2]     * A[2*lda+1]) / det; // D = -bi + ch
        B[2]      = (  A[1]     * A[lda+2]   - A[2]     * A[lda+1])   / det; // G = bf - ce
        B[ldb]    = (- A[lda]   * A[2*lda+2] + A[lda+2] * A[2*lda])   / det; // B = -di + fg
        B[ldb+1]  = (  A[0]     * A[2*lda+2] - A[2]     * A[2*lda])   / det; // E = ai - cg
        B[ldb+2]  = (- A[0]     * A[lda+2]   + A[2]     * A[lda])     / det; // H = -af + cd
        B[2*ldb]  = (  A[lda]   * A[2*lda+1] - A[lda+1] * A[2*lda])   / det; // C = dh - eg
        B[2*ldb+1]= ( -A[0]     * A[2*lda+1] + A[1]     * A[2*lda])   / det; // F = -ah + bg
        B[2*ldb+2]= (  A[0]     * A[lda+1]   - A[1]     * A[lda])     / det; // I = ae - bd

        // Calculate the size of ijb and b, then reallocate them.
        IType ndnz = 0;
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            if (j != i && B[i*ldb + j] != 0) { ++ndnz; }
          }
        }
        inverse->ndnz = ndnz;
        col_pos = 4; // shape[0] + 1
        inverse->capacity = 4 + ndnz;
        NM_REALLOC_N(inverse->a, DType, 4 + ndnz);
        NM_REALLOC_N(inverse->ija, IType, 4 + ndnz);
        b = reinterpret_cast<DType*>(inverse->a);
        ijb = reinterpret_cast<IType *>(inverse->ija);

        for (int i = 0; i < 3; ++i) {
          ijb[i] = col_pos;
          for (int j = 0; j < 3; ++j) {
            if (j == i) {
              b[i] = B[i*ldb + j];
            }
            else if (B[i*ldb + j] != 0) {
              b[col_pos] = B[i*ldb + j];
              ijb[col_pos] = j;
              ++col_pos;
            }
          }
        }
        b[3] = 0;
        ijb[3] = col_pos;
        delete [] B;
        delete [] A;
      } else if (M == 1) {
        b[0] = 1 / a[0];
      } else {
        rb_raise(rb_eNotImpError, "exact inverse calculation needed for matrices larger than 3x3");
      }
    }

    /*
     * Function signature conversion for calling CBLAS' gemm functions as directly as possible.
     *
     * For documentation: http://www.netlib.org/blas/dgemm.f
     */
    template <typename DType>
    inline static void cblas_gemm(const enum CBLAS_ORDER order,
                                  const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_TRANSPOSE trans_b,
                                  int m, int n, int k,
                                  void* alpha,
                                  void* a, int lda,
                                  void* b, int ldb,
                                  void* beta,
                                  void* c, int ldc)
    {
      gemm<DType>(order, trans_a, trans_b, m, n, k, reinterpret_cast<DType*>(alpha),
                  reinterpret_cast<DType*>(a), lda,
                  reinterpret_cast<DType*>(b), ldb, reinterpret_cast<DType*>(beta),
                  reinterpret_cast<DType*>(c), ldc);
    }


    /*
     * Function signature conversion for calling CBLAS's gemv functions as directly as possible.
     *
     * For documentation: http://www.netlib.org/lapack/double/dgetrf.f
     */
    template <typename DType>
    inline static bool cblas_gemv(const enum CBLAS_TRANSPOSE trans,
                                  const int m, const int n,
                                  const void* alpha,
                                  const void* a, const int lda,
                                  const void* x, const int incx,
                                  const void* beta,
                                  void* y, const int incy)
    {
      return gemv<DType>(trans,
                         m, n, reinterpret_cast<const DType*>(alpha),
                         reinterpret_cast<const DType*>(a), lda,
                         reinterpret_cast<const DType*>(x), incx, reinterpret_cast<const DType*>(beta),
                         reinterpret_cast<DType*>(y), incy);
    }


    /*
     * Function signature conversion for calling CBLAS' trsm functions as directly as possible.
     *
     * For documentation: http://www.netlib.org/blas/dtrsm.f
     */
    template <typename DType>
    inline static void cblas_trsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                                   const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_DIAG diag,
                                   const int m, const int n, const void* alpha, const void* a,
                                   const int lda, void* b, const int ldb)
    {
      trsm<DType>(order, side, uplo, trans_a, diag, m, n, *reinterpret_cast<const DType*>(alpha),
                  reinterpret_cast<const DType*>(a), lda, reinterpret_cast<DType*>(b), ldb);
    }

  }
} // end of namespace nm::math


extern "C" {

///////////////////
// Ruby Bindings //
///////////////////

void nm_math_init_blas() {
  VALUE cNMatrix_Internal = rb_define_module_under(cNMatrix, "Internal");

  rb_define_singleton_method(cNMatrix, "has_clapack?", (METHOD)nm_has_clapack, 0);

  VALUE cNMatrix_Internal_LAPACK = rb_define_module_under(cNMatrix_Internal, "LAPACK");

  /* ATLAS-CLAPACK Functions that are implemented internally */
  rb_define_singleton_method(cNMatrix_Internal_LAPACK, "clapack_getrf", (METHOD)nm_clapack_getrf, 5);
  rb_define_singleton_method(cNMatrix_Internal_LAPACK, "clapack_getrs", (METHOD)nm_clapack_getrs, 9);
  rb_define_singleton_method(cNMatrix_Internal_LAPACK, "clapack_laswp", (METHOD)nm_clapack_laswp, 7);

  VALUE cNMatrix_Internal_BLAS = rb_define_module_under(cNMatrix_Internal, "BLAS");

  rb_define_singleton_method(cNMatrix_Internal_BLAS, "cblas_scal", (METHOD)nm_cblas_scal, 4);
  rb_define_singleton_method(cNMatrix_Internal_BLAS, "cblas_nrm2", (METHOD)nm_cblas_nrm2, 3);
  rb_define_singleton_method(cNMatrix_Internal_BLAS, "cblas_asum", (METHOD)nm_cblas_asum, 3);
  rb_define_singleton_method(cNMatrix_Internal_BLAS, "cblas_rot",  (METHOD)nm_cblas_rot,  7);
  rb_define_singleton_method(cNMatrix_Internal_BLAS, "cblas_rotg", (METHOD)nm_cblas_rotg, 1);
  rb_define_singleton_method(cNMatrix_Internal_BLAS, "cblas_imax", (METHOD)nm_cblas_imax, 3);

  rb_define_singleton_method(cNMatrix_Internal_BLAS, "cblas_gemm", (METHOD)nm_cblas_gemm, 14);
  rb_define_singleton_method(cNMatrix_Internal_BLAS, "cblas_gemv", (METHOD)nm_cblas_gemv, 11);
  rb_define_singleton_method(cNMatrix_Internal_BLAS, "cblas_trsm", (METHOD)nm_cblas_trsm, 12);
}

/*
 * call-seq:
 *     NMatrix::BLAS.cblas_scal(n, alpha, vector, inc) -> NMatrix
 *
 * BLAS level 1 function +scal+. Works with all dtypes.
 *
 * Scale +vector+ in-place by +alpha+ and also return it. The operation is as
 * follows:
 *  x <- alpha * x
 *
 * - +n+ -> Number of elements of +vector+.
 * - +alpha+ -> Scalar value used in the operation.
 * - +vector+ -> NMatrix of shape [n,1] or [1,n]. Modified in-place.
 * - +inc+ -> Increment used in the scaling function. Should generally be 1.
 */
static VALUE nm_cblas_scal(VALUE self, VALUE n, VALUE alpha, VALUE vector, VALUE incx) {
  nm::dtype_t dtype = NM_DTYPE(vector);

  void* scalar = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, scalar);

  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::cblas_scal, void, const int n,
      const void* scalar, void* x, const int incx);

  ttable[dtype](FIX2INT(n), scalar, NM_STORAGE_DENSE(vector)->elements,
      FIX2INT(incx));

  return vector;
}

/*
 * Call any of the cblas_xrotg functions as directly as possible.
 *
 * xROTG computes the elements of a Givens plane rotation matrix such that:
 *
 *  |  c s |   | a |   | r |
 *  | -s c | * | b | = | 0 |
 *
 * where r = +- sqrt( a**2 + b**2 ) and c**2 + s**2 = 1.
 *
 * The Givens plane rotation can be used to introduce zero elements into a matrix selectively.
 *
 * This function differs from most of the other raw BLAS accessors. Instead of
 * providing a, b, c, s as arguments, you should only provide a and b (the
 * inputs), and you should provide them as the first two elements of any dense
 * NMatrix type.
 *
 * The outputs [c,s] will be returned in a Ruby Array at the end; the input
 * NMatrix will also be modified in-place.
 *
 * This function, like the other cblas_ functions, does minimal type-checking.
 */
static VALUE nm_cblas_rotg(VALUE self, VALUE ab) {
  static void (*ttable[nm::NUM_DTYPES])(void* a, void* b, void* c, void* s) = {
      NULL, NULL, NULL, NULL, NULL, // can't represent c and s as integers, so no point in having integer operations.
      nm::math::cblas_rotg<float>,
      nm::math::cblas_rotg<double>,
      nm::math::cblas_rotg<nm::Complex64>,
      nm::math::cblas_rotg<nm::Complex128>,
      NULL //nm::math::cblas_rotg<nm::RubyObject>
  };

  nm::dtype_t dtype = NM_DTYPE(ab);

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this operation undefined for integer vectors");
    return Qnil;

  } else {
    NM_CONSERVATIVE(nm_register_value(&self));
    NM_CONSERVATIVE(nm_register_value(&ab));
    void *pC = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]),
         *pS = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);

    // extract A and B from the NVector (first two elements)
    void* pA = NM_STORAGE_DENSE(ab)->elements;
    void* pB = (char*)(NM_STORAGE_DENSE(ab)->elements) + DTYPE_SIZES[dtype];
    // c and s are output

    ttable[dtype](pA, pB, pC, pS);

    VALUE result = rb_ary_new2(2);

    if (dtype == nm::RUBYOBJ) {
      rb_ary_store(result, 0, *reinterpret_cast<VALUE*>(pC));
      rb_ary_store(result, 1, *reinterpret_cast<VALUE*>(pS));
    } else {
      rb_ary_store(result, 0, nm::rubyobj_from_cval(pC, dtype).rval);
      rb_ary_store(result, 1, nm::rubyobj_from_cval(pS, dtype).rval);
    }
    NM_CONSERVATIVE(nm_unregister_value(&ab));
    NM_CONSERVATIVE(nm_unregister_value(&self));
    return result;
  }
}


/*
 * Call any of the cblas_xrot functions as directly as possible.
 *
 * xROT is a BLAS level 1 routine (taking two vectors) which applies a plane rotation.
 *
 * It's tough to find documentation on xROT. Here are what we think the arguments are for:
 *  * n     :: number of elements to consider in x and y
 *  * x     :: a vector (expects an NVector)
 *  * incx  :: stride of x
 *  * y     :: a vector (expects an NVector)
 *  * incy  :: stride of y
 *  * c     :: cosine of the angle of rotation
 *  * s     :: sine of the angle of rotation
 *
 * Note that c and s will be the same dtype as x and y, except when x and y are complex. If x and y are complex, c and s
 * will be float for Complex64 or double for Complex128.
 *
 * You probably don't want to call this function. Instead, why don't you try rot, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_cblas_rot(VALUE self, VALUE n, VALUE x, VALUE incx, VALUE y, VALUE incy, VALUE c, VALUE s) {
  static void (*ttable[nm::NUM_DTYPES])(const int N, void*, const int, void*, const int, const void*, const void*) = {
      NULL, NULL, NULL, NULL, NULL, // can't represent c and s as integers, so no point in having integer operations.
      nm::math::cblas_rot<float,float>,
      nm::math::cblas_rot<double,double>,
      nm::math::cblas_rot<nm::Complex64,float>,
      nm::math::cblas_rot<nm::Complex128,double>,
      nm::math::cblas_rot<nm::RubyObject,nm::RubyObject>
  };

  nm::dtype_t dtype = NM_DTYPE(x);


  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this operation undefined for integer vectors");
    return Qfalse;
  } else {
    void *pC, *pS;

    // We need to ensure the cosine and sine arguments are the correct dtype -- which may differ from the actual dtype.
    if (dtype == nm::COMPLEX64) {
      pC = NM_ALLOCA_N(float,1);
      pS = NM_ALLOCA_N(float,1);
      rubyval_to_cval(c, nm::FLOAT32, pC);
      rubyval_to_cval(s, nm::FLOAT32, pS);
    } else if (dtype == nm::COMPLEX128) {
      pC = NM_ALLOCA_N(double,1);
      pS = NM_ALLOCA_N(double,1);
      rubyval_to_cval(c, nm::FLOAT64, pC);
      rubyval_to_cval(s, nm::FLOAT64, pS);
    } else {
      pC = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
      pS = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
      rubyval_to_cval(c, dtype, pC);
      rubyval_to_cval(s, dtype, pS);
    }


    ttable[dtype](FIX2INT(n), NM_STORAGE_DENSE(x)->elements, FIX2INT(incx), NM_STORAGE_DENSE(y)->elements, FIX2INT(incy), pC, pS);

    return Qtrue;
  }
}


/*
 * Call any of the cblas_xnrm2 functions as directly as possible.
 *
 * xNRM2 is a BLAS level 1 routine which calculates the 2-norm of an n-vector x.
 *
 * Arguments:
 *  * n     :: length of x, must be at least 0
 *  * x     :: pointer to first entry of input vector
 *  * incx  :: stride of x, must be POSITIVE (ATLAS says non-zero, but 3.8.4 code only allows positive)
 *
 * You probably don't want to call this function. Instead, why don't you try nrm2, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_cblas_nrm2(VALUE self, VALUE n, VALUE x, VALUE incx) {

  static void (*ttable[nm::NUM_DTYPES])(const int N, const void* X, const int incX, void* sum) = {
      NULL, NULL, NULL, NULL, NULL, // no help for integers
      nm::math::cblas_nrm2<float32_t>,
      nm::math::cblas_nrm2<float64_t>,
      nm::math::cblas_nrm2<nm::Complex64>,
      nm::math::cblas_nrm2<nm::Complex128>,
      nm::math::cblas_nrm2<nm::RubyObject>
  };

  nm::dtype_t dtype  = NM_DTYPE(x);

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this operation undefined for integer vectors");
    return Qnil;

  } else {
    // Determine the return dtype and allocate it
    nm::dtype_t rdtype = dtype;
    if      (dtype == nm::COMPLEX64)  rdtype = nm::FLOAT32;
    else if (dtype == nm::COMPLEX128) rdtype = nm::FLOAT64;

    void *Result = NM_ALLOCA_N(char, DTYPE_SIZES[rdtype]);

    ttable[dtype](FIX2INT(n), NM_STORAGE_DENSE(x)->elements, FIX2INT(incx), Result);

    return nm::rubyobj_from_cval(Result, rdtype).rval;
  }
}



/*
 * Call any of the cblas_xasum functions as directly as possible.
 *
 * xASUM is a BLAS level 1 routine which calculates the sum of absolute values of the entries
 * of a vector x.
 *
 * Arguments:
 *  * n     :: length of x, must be at least 0
 *  * x     :: pointer to first entry of input vector
 *  * incx  :: stride of x, must be POSITIVE (ATLAS says non-zero, but 3.8.4 code only allows positive)
 *
 * You probably don't want to call this function. Instead, why don't you try asum, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_cblas_asum(VALUE self, VALUE n, VALUE x, VALUE incx) {

  static void (*ttable[nm::NUM_DTYPES])(const int N, const void* X, const int incX, void* sum) = {
      nm::math::cblas_asum<uint8_t>,
      nm::math::cblas_asum<int8_t>,
      nm::math::cblas_asum<int16_t>,
      nm::math::cblas_asum<int32_t>,
      nm::math::cblas_asum<int64_t>,
      nm::math::cblas_asum<float32_t>,
      nm::math::cblas_asum<float64_t>,
      nm::math::cblas_asum<nm::Complex64>,
      nm::math::cblas_asum<nm::Complex128>,
      nm::math::cblas_asum<nm::RubyObject>
  };

  nm::dtype_t dtype  = NM_DTYPE(x);

  // Determine the return dtype and allocate it
  nm::dtype_t rdtype = dtype;
  if      (dtype == nm::COMPLEX64)  rdtype = nm::FLOAT32;
  else if (dtype == nm::COMPLEX128) rdtype = nm::FLOAT64;

  void *Result = NM_ALLOCA_N(char, DTYPE_SIZES[rdtype]);

  ttable[dtype](FIX2INT(n), NM_STORAGE_DENSE(x)->elements, FIX2INT(incx), Result);

  return nm::rubyobj_from_cval(Result, rdtype).rval;
}

/*
 * call-seq:
 *    NMatrix::BLAS.cblas_imax(n, vector, inc) -> Fixnum
 *
 * BLAS level 1 routine.
 *
 * Return the index of the largest element of +vector+.
 *
 * - +n+ -> Vector's size. Generally, you can use NMatrix#rows or NMatrix#cols.
 * - +vector+ -> A NMatrix of shape [n,1] or [1,n] with any dtype.
 * - +inc+ -> It's the increment used when searching. Use 1 except if you know
 *   what you're doing.
 */
static VALUE nm_cblas_imax(VALUE self, VALUE n, VALUE x, VALUE incx) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::cblas_imax, int, const int n, const void* x, const int incx);

  nm::dtype_t dtype = NM_DTYPE(x);

  int index = ttable[dtype](FIX2INT(n), NM_STORAGE_DENSE(x)->elements, FIX2INT(incx));

  // Convert to Ruby's Int value.
  return INT2FIX(index);
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
 * You probably don't want to call this function. Instead, why don't you try gemm, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_cblas_gemm(VALUE self,
                           VALUE order,
                           VALUE trans_a, VALUE trans_b,
                           VALUE m, VALUE n, VALUE k,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb,
                           VALUE beta,
                           VALUE c, VALUE ldc)
{
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::cblas_gemm, void, const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_TRANSPOSE trans_b, int m, int n, int k, void* alpha, void* a, int lda, void* b, int ldb, void* beta, void* c, int ldc);

  nm::dtype_t dtype = NM_DTYPE(a);

  void *pAlpha = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]),
       *pBeta  = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, pAlpha);
  rubyval_to_cval(beta, dtype, pBeta);

  ttable[dtype](blas_order_sym(order), blas_transpose_sym(trans_a), blas_transpose_sym(trans_b), FIX2INT(m), FIX2INT(n), FIX2INT(k), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb), pBeta, NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc));

  return c;
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
                           VALUE y, VALUE incy)
{
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::cblas_gemv, bool, const enum CBLAS_TRANSPOSE, const int, const int, const void*, const void*, const int, const void*, const int, const void*, void*, const int)

  nm::dtype_t dtype = NM_DTYPE(a);

  void *pAlpha = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]),
       *pBeta  = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, pAlpha);
  rubyval_to_cval(beta, dtype, pBeta);

  return ttable[dtype](blas_transpose_sym(trans_a), FIX2INT(m), FIX2INT(n), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(x)->elements, FIX2INT(incx), pBeta, NM_STORAGE_DENSE(y)->elements, FIX2INT(incy)) ? Qtrue : Qfalse;
}


static VALUE nm_cblas_trsm(VALUE self,
                           VALUE order,
                           VALUE side, VALUE uplo,
                           VALUE trans_a, VALUE diag,
                           VALUE m, VALUE n,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb)
{
  static void (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const enum CBLAS_SIDE, const enum CBLAS_UPLO,
                                        const enum CBLAS_TRANSPOSE, const enum CBLAS_DIAG,
                                        const int m, const int n, const void* alpha, const void* a,
                                        const int lda, void* b, const int ldb) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::cblas_trsm<float>,
      nm::math::cblas_trsm<double>,
      nm::math::cblas_trsm<nm::Complex64>,
      nm::math::cblas_trsm<nm::Complex128>,
      nm::math::cblas_trsm<nm::RubyObject>
  };

  nm::dtype_t dtype = NM_DTYPE(a);

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    void *pAlpha = NM_ALLOCA_N(char, DTYPE_SIZES[dtype]);
    rubyval_to_cval(alpha, dtype, pAlpha);

    ttable[dtype](blas_order_sym(order), blas_side_sym(side), blas_uplo_sym(uplo), blas_transpose_sym(trans_a), blas_diag_sym(diag), FIX2INT(m), FIX2INT(n), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb));
  }

  return Qtrue;
}

/* Call any of the clapack_xgetrf functions as directly as possible.
 *
 * The clapack_getrf functions (dgetrf, sgetrf, cgetrf, and zgetrf) compute an LU factorization of a general M-by-N
 * matrix A using partial pivoting with row interchanges.
 *
 * The factorization has the form:
 *    A = P * L * U
 * where P is a permutation matrix, L is lower triangular with unit diagonal elements (lower trapezoidal if m > n),
 * and U is upper triangular (upper trapezoidal if m < n).
 *
 * This is the right-looking level 3 BLAS version of the algorithm.
 *
 * == Arguments
 * See: http://www.netlib.org/lapack/double/dgetrf.f
 * (You don't need argument 5; this is the value returned by this function.)
 *
 * You probably don't want to call this function. Instead, why don't you try clapack_getrf, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 *
 * Returns an array giving the pivot indices (normally these are argument #5).
 */
static VALUE nm_clapack_getrf(VALUE self, VALUE order, VALUE m, VALUE n, VALUE a, VALUE lda) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const int m, const int n, void* a, const int lda, int* ipiv) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::clapack_getrf<float>,
      nm::math::clapack_getrf<double>,
      nm::math::clapack_getrf<nm::Complex64>,
      nm::math::clapack_getrf<nm::Complex128>,
      nm::math::clapack_getrf<nm::RubyObject>
  };

  int M = FIX2INT(m),
      N = FIX2INT(n);

  // Allocate the pivot index array, which is of size MIN(M, N).
  size_t ipiv_size = std::min(M,N);
  int* ipiv = NM_ALLOCA_N(int, ipiv_size);

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    // Call either our version of getrf or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), M, N, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), ipiv);
  }

  // Result will be stored in a. We return ipiv as an array.
  VALUE ipiv_array = rb_ary_new2(ipiv_size);
  for (size_t i = 0; i < ipiv_size; ++i) {
    rb_ary_store(ipiv_array, i, INT2FIX(ipiv[i]));
  }

  return ipiv_array;
}


/*
 * Call any of the clapack_xgetrs functions as directly as possible.
 */
static VALUE nm_clapack_getrs(VALUE self, VALUE order, VALUE trans, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE ipiv, VALUE b, VALUE ldb) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans, const int N,
                                       const int NRHS, const void* A, const int lda, const int* ipiv, void* B,
                                       const int ldb) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::clapack_getrs<float>,
      nm::math::clapack_getrs<double>,
      nm::math::clapack_getrs<nm::Complex64>,
      nm::math::clapack_getrs<nm::Complex128>,
      nm::math::clapack_getrs<nm::RubyObject>
  };

  // Allocate the C version of the pivot index array
  int* ipiv_;
  if (!RB_TYPE_P(ipiv, T_ARRAY)) {
    rb_raise(rb_eArgError, "ipiv must be of type Array");
  } else {
    ipiv_ = NM_ALLOCA_N(int, RARRAY_LEN(ipiv));
    for (int index = 0; index < RARRAY_LEN(ipiv); ++index) {
      ipiv_[index] = FIX2INT( RARRAY_AREF(ipiv, index) );
    }
  }

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {

    // Call either our version of getrs or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), blas_transpose_sym(trans), FIX2INT(n), FIX2INT(nrhs), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda),
                        ipiv_, NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb));
  }

  // b is both returned and modified directly in the argument list.
  return b;
}

/*
 * Simple way to check from within Ruby code if clapack functions are available, without
 * having to wait around for an exception to be thrown.
 */
static VALUE nm_has_clapack(VALUE self) {
  return Qfalse;
}

/*
 * Call any of the clapack_xlaswp functions as directly as possible.
 *
 * Note that LAPACK's xlaswp functions accept a column-order matrix, but NMatrix uses row-order. Thus, n should be the
 * number of rows and lda should be the number of columns, no matter what it says in the documentation for dlaswp.f.
 */
static VALUE nm_clapack_laswp(VALUE self, VALUE n, VALUE a, VALUE lda, VALUE k1, VALUE k2, VALUE ipiv, VALUE incx) {
  static void (*ttable[nm::NUM_DTYPES])(const int n, void* a, const int lda, const int k1, const int k2, const int* ipiv, const int incx) = {
      nm::math::clapack_laswp<uint8_t>,
      nm::math::clapack_laswp<int8_t>,
      nm::math::clapack_laswp<int16_t>,
      nm::math::clapack_laswp<int32_t>,
      nm::math::clapack_laswp<int64_t>,
      nm::math::clapack_laswp<float>,
      nm::math::clapack_laswp<double>,
      nm::math::clapack_laswp<nm::Complex64>,
      nm::math::clapack_laswp<nm::Complex128>,
      nm::math::clapack_laswp<nm::RubyObject>
  };

  // Allocate the C version of the pivot index array
  int* ipiv_;
  if (!RB_TYPE_P(ipiv, T_ARRAY)) {
    rb_raise(rb_eArgError, "ipiv must be of type Array");
  } else {
    ipiv_ = NM_ALLOCA_N(int, RARRAY_LEN(ipiv));
    for (int index = 0; index < RARRAY_LEN(ipiv); ++index) {
      ipiv_[index] = FIX2INT( RARRAY_AREF(ipiv, index) );
    }
  }

  // Call either our version of laswp or the LAPACK version.
  ttable[NM_DTYPE(a)](FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), FIX2INT(k1), FIX2INT(k2), ipiv_, FIX2INT(incx));

  // a is both returned and modified directly in the argument list.
  return a;
}


/*
 * C accessor for calculating an exact determinant. Dense matrix version.
 */
void nm_math_det_exact_from_dense(const int M, const void* elements, const int lda,
        nm::dtype_t dtype, void* result) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::det_exact_from_dense, void, const int M,
          const void* A_elements, const int lda, void* result_arg);

  ttable[dtype](M, elements, lda, result);
}

/*
 * C accessor for calculating an exact determinant. Yale matrix version.
 */
void nm_math_det_exact_from_yale(const int M, const YALE_STORAGE* storage, const int lda,
        nm::dtype_t dtype, void* result) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::det_exact_from_yale, void, const int M,
          const YALE_STORAGE* storage, const int lda, void* result_arg);

  ttable[dtype](M, storage, lda, result);
}

/*
 * C accessor for solving a system of linear equations.
 */
void nm_math_solve(VALUE lu, VALUE b, VALUE x, VALUE ipiv) {
  int* pivot = new int[RARRAY_LEN(ipiv)];

  for (int i = 0; i < RARRAY_LEN(ipiv); ++i) {
    pivot[i] = FIX2INT(rb_ary_entry(ipiv, i));
  }

  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::solve, void, const int, const void*, const void*, void*, const int*);

  ttable[NM_DTYPE(x)](NM_SHAPE0(b), NM_STORAGE_DENSE(lu)->elements,
    NM_STORAGE_DENSE(b)->elements, NM_STORAGE_DENSE(x)->elements, pivot);
}

/*
 * C accessor for reducing a matrix to hessenberg form.
 */
void nm_math_hessenberg(VALUE a) {
  static void (*ttable[nm::NUM_DTYPES])(const int, void*) = {
      NULL, NULL, NULL, NULL, NULL, // does not support ints
      nm::math::hessenberg<float>,
      nm::math::hessenberg<double>,
      NULL, NULL, // does not support Complex
      NULL // no support for Ruby Object
  };

  ttable[NM_DTYPE(a)](NM_SHAPE0(a), NM_STORAGE_DENSE(a)->elements);
}
/*
 * C accessor for calculating an in-place inverse.
 */
void nm_math_inverse(const int M, void* a_elements, nm::dtype_t dtype) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::inverse, void, const int, void*);

  ttable[dtype](M, a_elements);
}

/*
 * C accessor for calculating an exact inverse. Dense matrix version.
 */
void nm_math_inverse_exact_from_dense(const int M, const void* A_elements,
    const int lda, void* B_elements, const int ldb, nm::dtype_t dtype) {

  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::inverse_exact_from_dense, void,
      const int, const void*, const int, void*, const int);

  ttable[dtype](M, A_elements, lda, B_elements, ldb);
}

/*
 * C accessor for calculating an exact inverse. Yale matrix version.
 */
void nm_math_inverse_exact_from_yale(const int M, const YALE_STORAGE* storage,
    const int lda, YALE_STORAGE* inverse, const int ldb, nm::dtype_t dtype) {

  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::inverse_exact_from_yale, void,
      const int, const YALE_STORAGE*, const int, YALE_STORAGE*, const int);

  ttable[dtype](M, storage, lda, inverse, ldb);
}

/*
 * Transpose an array of elements that represent a row-major dense matrix. Does not allocate anything, only does an memcpy.
 */
void nm_math_transpose_generic(const size_t M, const size_t N, const void* A, const int lda, void* B, const int ldb, size_t element_size) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {

      memcpy(reinterpret_cast<char*>(B) + (i*ldb+j)*element_size,
             reinterpret_cast<const char*>(A) + (j*lda+i)*element_size,
             element_size);

    }
  }
}


} // end of extern "C" block
