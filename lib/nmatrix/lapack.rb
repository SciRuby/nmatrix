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
# SciRuby is Copyright (c) 2010 - 2014, Ruby Science Foundation
# NMatrix is Copyright (c) 2012 - 2014, John Woods and the Ruby Science Foundation
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
      #   - +b+ ->
      #   - +ldb+ ->
      #   - +ipiv+ -> A pivot array (if nil, one will be generated with +clapack_getrf+)
      # * *Returns* :
      #   -
      # * *Raises* :
      #   - ++ ->
      #
      def clapack_gesv(order, n, nrhs, a, lda, b, ldb, ipiv=nil)
        ipiv ||= clapack_getrf(order, n, n, a, lda)
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

      #     laswp(matrix, ipiv) -> NMatrix
      #
      # Permute the columns of a matrix (in-place) according to the Array +ipiv+.
      #
      def laswp(matrix, ipiv)
        raise(ArgumentError, "expected NMatrix for argument 0") unless matrix.is_a?(NMatrix)
        raise(StorageTypeError, "LAPACK functions only work on :dense NMatrix instances") unless matrix.stype == :dense
        raise(ArgumentError, "expected Array ipiv to have no more entries than NMatrix a has columns") if ipiv.size > matrix.shape[1]

        clapack_laswp(matrix.shape[0], matrix, matrix.shape[1], 0, ipiv.size-1, ipiv, 1)
      end

      def alloc_svd_result(matrix)
        [
          NMatrix.new(matrix.shape[0], dtype: matrix.dtype),
          NMatrix.new([matrix.shape[0],1], dtype: matrix.dtype),
          NMatrix.new(matrix.shape[1], dtype: matrix.dtype)
        ]
      end


      #
      # call-seq:
      #     gesvd(matrix) -> [u, sigma, v_transpose]
      #     gesvd(matrix) -> [u, sigma, v_conjugate_transpose] # complex
      #
      # Compute the singular value decomposition of a matrix using LAPACK's GESVD function.
      #
      # Optionally accepts a +workspace_size+ parameter, which will be honored only if it is larger than what LAPACK
      # requires.
      #
      def gesvd(matrix, workspace_size=1)
        result = alloc_svd_result(matrix)
        NMatrix::LAPACK::lapack_gesvd(:a, :a, matrix.shape[0], matrix.shape[1], matrix, matrix.shape[0], result[1], result[0], matrix.shape[0], result[2], matrix.shape[1], workspace_size)
        result
      end

      #
      # call-seq:
      #     gesdd(matrix) -> [u, sigma, v_transpose]
      #     gesdd(matrix) -> [u, sigma, v_conjugate_transpose] # complex
      #
      # Compute the singular value decomposition of a matrix using LAPACK's GESDD function. This uses a divide-and-conquer
      # strategy. See also #gesvd.
      #
      # Optionally accepts a +workspace_size+ parameter, which will be honored only if it is larger than what LAPACK
      # requires.
      #
      def gesdd(matrix, workspace_size=100000)
        result = alloc_svd_result(matrix)
        NMatrix::LAPACK::lapack_gesdd(:a, matrix.shape[0], matrix.shape[1], matrix, matrix.shape[0], result[1], result[0], matrix.shape[0], result[2], matrix.shape[1], workspace_size)
        result
      end

      def alloc_evd_result(matrix)
        [
          NMatrix.new(matrix.shape[0], dtype: matrix.dtype),
          NMatrix.new(matrix.shape[0], dtype: matrix.dtype),
          NMatrix.new([matrix.shape[0],1], dtype: matrix.dtype),
          NMatrix.new([matrix.shape[0],1], dtype: matrix.dtype),
        ]
      end


      #
      # call-seq:
      #     geev(matrix) -> [eigenvalues, left_eigenvectors, right_eigenvectors]
      #     geev(matrix, :left) -> [eigenvalues, left_eigenvectors]
      #     geev(matrix, :right) -> [eigenvalues, right_eigenvectors]
      #
      # Perform eigenvalue decomposition on a matrix using LAPACK's xGEEV function.
      #
      def geev(matrix, which=:both)
        jobvl = (which == :both || which == :left) ? :left : false
        jobvr = (which == :both || which == :right) ? :right : false

        # Copy the matrix so it doesn't get overwritten.
        temporary_matrix = matrix.clone

        # Outputs
        real_eigenvalues = NMatrix.new([matrix.shape[0], 1], dtype: matrix.dtype)
        imag_eigenvalues = NMatrix.new([matrix.shape[0], 1], dtype: matrix.dtype)

        left_output      = jobvl == :left ? matrix.clone_structure : NMatrix.new(1, dtype: matrix.dtype)
        right_output     = jobvr == :right ? matrix.clone_structure : NMatrix.new(1, dtype: matrix.dtype)

        NMatrix::LAPACK::lapack_geev(jobvl, # compute left eigenvectors of A?
                                     jobvr, # compute right eigenvectors of A? (left eigenvectors of A**T)
                                     matrix.shape[0], # order of the matrix
                                     temporary_matrix,# input matrix (used as work)
                                     matrix.shape[0], # leading dimension of matrix
                                     real_eigenvalues,# real part of computed eigenvalues
                                     imag_eigenvalues,# imag part of computed eigenvalues
                                     left_output,     # left eigenvectors, if applicable
                                     left_output.shape[0], # leading dimension of left_output
                                     right_output,    # right eigenvectors, if applicable
                                     right_output.shape[0], # leading dimension of right_output
                                     2*matrix.shape[0])

        # Put the real and imaginary parts together
        eigenvalues = real_eigenvalues.to_a.flatten.map.with_index do |real,i|
          imag_eigenvalues[i] != 0 ? Complex(real, imag_eigenvalues[i]) : real
        end

        if which == :both
          return [eigenvalues, left_output.transpose, right_output.transpose]
        elsif which == :left
          return [eigenvalues, left_output.transpose]
        else
          return [eigenvalues, right_output]
        end
      end

    end
  end
end
