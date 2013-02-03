#--
# = NMatrix
#
# A linear algebra library for scientific computation in Ruby.
# NMatrix is part of SciRuby.
#
# NMatrix was originally inspired by and derived from NArray, by
# Masahiro Tanaka: http://narray.rubyforge.org
#
# == Copyright Information
#
# SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
# NMatrix is Copyright (c) 2012, Ruby Science Foundation
#
# Please see LICENSE.txt for additional copyright notices.
#
# == Contributing
#
# By contributing source code to SciRuby, you agree to be bound by
# our Contributor Agreement:
#
# * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
#
# == lapack.rb
#
# This file contains LAPACK functions accessible in their C versions,
# e.g., NMatrix::LAPACK::clapack_func. There are some exceptions,
# such as clapack_gesv, which is implemented in Ruby but calls
# clapack_getrf and clapack_getrs.
#
# Note: most of these functions are borrowed from ATLAS, which is available under a BSD-
# style license.
#++

class NMatrix
  module LAPACK
    class << self
      #
      # call-seq:
      #     clapack_gesv(order, n, nrhs, a, lda, ipiv, b, ldb) -> NMatrix
      #
      # Computes the solution to a system of linear equations
      #   A * X = B,
      # where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
      #
      # The LU factorization used to factor A is dependent on the +order+
      # parameter, as detailed in the leading comments of clapack_getrf.
      #
      # The factored form of A is then used solve the system of equations
      # A * X = B.
      #
      # A is overwritten with the appropriate LU factorization, and B, which
      # contains B on input, is overwritten with the solution X on output.
      #
      # From ATLAS 3.8.0.
      #
      # Note: Because this function is implemented in Ruby, the ATLAS lib 
      # version is never called! For float32, float64, complex64, and 
      # complex128, the ATLAS lib versions of getrf and getrs *will* be called.
      #
      # * *Arguments* :
      #   - +order+ ->
      #   - +n+ ->
      #   - +nrhs+ ->
      #   - +a+ ->
      #   - +lda+ ->
      #   - +ipiv+ ->
      #   - +b+ ->
      #   - +ldb+ ->
      # * *Returns* :
      #   -
      # * *Raises* :
      #   - ++ ->
      #
      def clapack_gesv(order, n, nrhs, a, lda, ipiv, b, ldb)
        clapack_getrf(order, n, n, a, lda, ipiv)
        clapack_getrs(order, :no_transpose, n, nrhs, a, lda, ipiv, b, ldb)
      end

      #
      # call-seq:
      #     clapack_posv(order, uplo, n ,nrhs, a, lda, b, ldb) -> ...
      #
      # TODO Complete this description.
      #
      # Computes the solution to a real system of linear equations
      #   A * X = B,
      # where A is an N-by-N symmetric positive definite matrix and X and B
      # are N-by-NRHS matrices.
      #
      # The Cholesky decomposition is used to factor A as
      #   A = U**T* U,  if UPLO = 'U', or
      #   A = L * L**T,  if UPLO = 'L',
      # where U is an upper triangular matrix and L is a lower triangular
      # matrix.  The factored form of A is then used to solve the system of
      # equations A * X = B.
      #
      # From ATLAS 3.8.0.
      #
      # Note: Because this function is implemented in Ruby, the ATLAS lib
      # version is never called! For float32, float64, complex64, and 
      # complex128, the ATLAS lib versions of potrf and potrs *will* be called.
      #
      # * *Arguments* :
      #   - +order+ ->
      #   - +uplo+ ->
      #   - +n+ ->
      #   - +nrhs+ ->
      #   - +a+ ->
      #   - +lda+ ->
      #   - +b+ ->
      #   - +ldb+ ->
      # * *Returns* :
      #   -
      # * *Raises* :
      #   - ++ ->
      #
      def clapack_posv(order, uplo, n, nrhs, a, lda, b, ldb)
        clapack_potrf(order, uplo, n, a, lda)
        clapack_potrs(order, uplo, n, nrhs, a, lda, b, ldb)
      end

    end
  end
end