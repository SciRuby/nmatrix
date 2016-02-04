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
# == lapacke.rb
#
# ruby file for the nmatrix-lapacke gem. Loads the C extension and defines
# nice ruby interfaces for LAPACK functions.
#++

require 'nmatrix/nmatrix.rb' #need to have nmatrix required first or else bad things will happen
require_relative 'lapack_ext_common'

NMatrix.register_lapack_extension("nmatrix-lapacke")

require "nmatrix_lapacke.so"

class NMatrix
  #Add functions from the LAPACKE C extension to the main LAPACK and BLAS modules.
  #This will overwrite the original functions where applicable.
  module LAPACK
    class << self
      NMatrix::LAPACKE::LAPACK.singleton_methods.each do |m|
        define_method m, NMatrix::LAPACKE::LAPACK.method(m).to_proc
      end
    end
  end

  module BLAS
    class << self
      NMatrix::LAPACKE::BLAS.singleton_methods.each do |m|
        define_method m, NMatrix::LAPACKE::BLAS.method(m).to_proc
      end
    end
  end

  module LAPACK
    class << self
      def posv(uplo, a, b)
        raise(ShapeError, "a must be square") unless a.dim == 2 && a.shape[0] == a.shape[1]
        raise(ShapeError, "number of rows of b must equal number of cols of a") unless a.shape[1] == b.shape[0]
        raise(StorageTypeError, "only works with dense matrices") unless a.stype == :dense && b.stype == :dense
        raise(DataTypeError, "only works for non-integer, non-object dtypes") if 
          a.integer_dtype? || a.object_dtype? || b.integer_dtype? || b.object_dtype?

        x     = b.clone
        clone = a.clone
        n = a.shape[0]
        nrhs = b.shape[1]
        lapacke_potrf(:row, uplo, n, clone, n)
        lapacke_potrs(:row, uplo, n, nrhs, clone, n, x, b.shape[1])
        x
      end

      def geev(matrix, which=:both)
        raise(StorageTypeError, "LAPACK functions only work on dense matrices") unless matrix.dense?
        raise(ShapeError, "eigenvalues can only be computed for square matrices") unless matrix.dim == 2 && matrix.shape[0] == matrix.shape[1]

        jobvl = (which == :both || which == :left) ? :t : false
        jobvr = (which == :both || which == :right) ? :t : false

        # Copy the matrix so it doesn't get overwritten.
        temporary_matrix = matrix.clone
        n = matrix.shape[0]

        # Outputs
        eigenvalues = NMatrix.new([n, 1], dtype: matrix.dtype) # For real dtypes this holds only the real part of the eigenvalues.
        imag_eigenvalues = matrix.complex_dtype? ? nil : NMatrix.new([n, 1], dtype: matrix.dtype) # For complex dtypes, this is unused.
        left_output      = jobvl ? matrix.clone_structure : nil
        right_output     = jobvr ? matrix.clone_structure : nil

        NMatrix::LAPACK::lapacke_geev(:row,
                                      jobvl, # compute left eigenvectors of A?
                                      jobvr, # compute right eigenvectors of A? (left eigenvectors of A**T)
                                      n, # order of the matrix
                                      temporary_matrix,# input matrix (used as work)
                                      n, # leading dimension of matrix
                                      eigenvalues,# real part of computed eigenvalues
                                      imag_eigenvalues,# imag part of computed eigenvalues
                                      left_output,     # left eigenvectors, if applicable
                                      n, # leading dimension of left_output
                                      right_output,    # right eigenvectors, if applicable
                                      n) # leading dimension of right_output


        # For real dtypes, transform left_output and right_output into correct forms.
        # If the j'th and the (j+1)'th eigenvalues form a complex conjugate
        # pair, then the j'th and (j+1)'th columns of the matrix are
        # the real and imag parts of the eigenvector corresponding
        # to the j'th eigenvalue.
        if !matrix.complex_dtype?
          complex_indices = []
          n.times do |i|
            complex_indices << i if imag_eigenvalues[i] != 0.0
          end

          if !complex_indices.empty?
            # For real dtypes, put the real and imaginary parts together
            eigenvalues = eigenvalues + imag_eigenvalues*Complex(0.0,1.0)
            left_output = left_output.cast(dtype: NMatrix.upcast(:complex64, matrix.dtype)) if left_output
            right_output = right_output.cast(dtype: NMatrix.upcast(:complex64, matrix.dtype)) if right_output
          end

          complex_indices.each_slice(2) do |i, _|
            if right_output
              right_output[0...n,i] = right_output[0...n,i] + right_output[0...n,i+1]*Complex(0.0,1.0)
              right_output[0...n,i+1] = right_output[0...n,i].complex_conjugate
            end

            if left_output
              left_output[0...n,i] = left_output[0...n,i] + left_output[0...n,i+1]*Complex(0.0,1.0)
              left_output[0...n,i+1] = left_output[0...n,i].complex_conjugate
            end
          end
        end

        if which == :both
          return [eigenvalues, left_output, right_output]
        elsif which == :left
          return [eigenvalues, left_output]
        else
          return [eigenvalues, right_output]
        end
      end

      def gesvd(matrix, workspace_size=1)
        result = alloc_svd_result(matrix)

        m = matrix.shape[0]
        n = matrix.shape[1]

        superb = NMatrix.new([[m,n].min], dtype: matrix.abs_dtype)

        NMatrix::LAPACK::lapacke_gesvd(:row, :a, :a, m, n, matrix, n, result[1], result[0], m, result[2], n, superb)
        result
      end

      def gesdd(matrix, workspace_size=nil)
        result = alloc_svd_result(matrix)

        m = matrix.shape[0]
        n = matrix.shape[1]

        NMatrix::LAPACK::lapacke_gesdd(:row, :a, m, n, matrix, n, result[1], result[0], m, result[2], n)
        result
      end
    end
  end

  def getrf!
    raise(StorageTypeError, "LAPACK functions only work on dense matrices") unless self.dense?

    ipiv = NMatrix::LAPACK::lapacke_getrf(:row, self.shape[0], self.shape[1], self, self.shape[1])

    return ipiv
  end

  def invert!
    raise(StorageTypeError, "invert only works on dense matrices currently") unless self.dense?
    raise(ShapeError, "Cannot invert non-square matrix") unless shape[0] == shape[1]
    raise(DataTypeError, "Cannot invert an integer matrix in-place") if self.integer_dtype?

    # Get the pivot array; factor the matrix
    n = self.shape[0]
    pivot = NMatrix::LAPACK::lapacke_getrf(:row, n, n, self, n)
    # Now calculate the inverse using the pivot array
    NMatrix::LAPACK::lapacke_getri(:row, n, self, n, pivot)

    self
  end

  def potrf!(which)
    raise(StorageTypeError, "LAPACK functions only work on dense matrices") unless self.dense?
    raise(ShapeError, "Cholesky decomposition only valid for square matrices") unless self.dim == 2 && self.shape[0] == self.shape[1]

    NMatrix::LAPACK::lapacke_potrf(:row, which, self.shape[0], self, self.shape[1])
  end

  def solve(b, opts = {})
    raise(ShapeError, "Must be called on square matrix") unless self.dim == 2 && self.shape[0] == self.shape[1]
    raise(ShapeError, "number of rows of b must equal number of cols of self") if 
      self.shape[1] != b.shape[0]
    raise(ArgumentError, "only works with dense matrices") if self.stype != :dense
    raise(ArgumentError, "only works for non-integer, non-object dtypes") if 
      integer_dtype? or object_dtype? or b.integer_dtype? or b.object_dtype?

    opts = { form: :general }.merge(opts)
    x    = b.clone
    n    = self.shape[0]
    nrhs = b.shape[1]

    case opts[:form] 
    when :general
      clone = self.clone
      ipiv = NMatrix::LAPACK.lapacke_getrf(:row, n, n, clone, n)
      NMatrix::LAPACK.lapacke_getrs(:row, :no_transpose, n, nrhs, clone, n, ipiv, x, nrhs)
      x
    when :upper_tri, :upper_triangular
      raise(ArgumentError, "upper triangular solver does not work with complex dtypes") if
        complex_dtype? or b.complex_dtype?
      NMatrix::BLAS::cblas_trsm(:row, :left, :upper, false, :nounit, n, nrhs, 1.0, self, n, x, nrhs)
      x
    when :lower_tri, :lower_triangular
      raise(ArgumentError, "lower triangular solver does not work with complex dtypes") if
        complex_dtype? or b.complex_dtype?
      NMatrix::BLAS::cblas_trsm(:row, :left, :lower, false, :nounit, n, nrhs, 1.0, self, n, x, nrhs)
      x
    when :pos_def, :positive_definite
      u, l = self.factorize_cholesky
      z = l.solve(b, form: :lower_tri)
      u.solve(z, form: :upper_tri)
    else
      raise(ArgumentError, "#{opts[:form]} is not a valid form option")
    end
  end
end
