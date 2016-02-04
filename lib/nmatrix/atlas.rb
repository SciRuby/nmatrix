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
# == atlas.rb
#
# ruby file for the nmatrix-atlas gem. Loads the C extension and defines
# nice ruby interfaces for ATLAS functions.
#++

require 'nmatrix/nmatrix.rb' #need to have nmatrix required first or else bad things will happen
require_relative 'lapack_ext_common'

NMatrix.register_lapack_extension("nmatrix-atlas")

require "nmatrix_atlas.so"

class NMatrix

  #Add functions from the ATLAS C extension to the main LAPACK and BLAS modules.
  #This will overwrite the original functions where applicable.
  module LAPACK
    class << self
      NMatrix::ATLAS::LAPACK.singleton_methods.each do |m|
        define_method m, NMatrix::ATLAS::LAPACK.method(m).to_proc
      end
    end
  end

  module BLAS
    class << self
      NMatrix::ATLAS::BLAS.singleton_methods.each do |m|
        define_method m, NMatrix::ATLAS::BLAS.method(m).to_proc
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
        clapack_potrf(:row, uplo, n, clone, n)
        # Must transpose b before and after: http://math-atlas.sourceforge.net/faq.html#RowSolve
        x = x.transpose
        clapack_potrs(:row, uplo, n, nrhs, clone, n, x, n)
        x.transpose
      end

      def geev(matrix, which=:both)
        raise(StorageTypeError, "LAPACK functions only work on dense matrices") unless matrix.dense?
        raise(ShapeError, "eigenvalues can only be computed for square matrices") unless matrix.dim == 2 && matrix.shape[0] == matrix.shape[1]

        jobvl = (which == :both || which == :left) ? :t : false
        jobvr = (which == :both || which == :right) ? :t : false

        n = matrix.shape[0]

        # Outputs
        eigenvalues = NMatrix.new([n, 1], dtype: matrix.dtype) # For real dtypes this holds only the real part of the eigenvalues.
        imag_eigenvalues = matrix.complex_dtype? ? nil : NMatrix.new([n, 1], dtype: matrix.dtype) # For complex dtypes, this is unused.
        left_output      = jobvl ? matrix.clone_structure : nil
        right_output     = jobvr ? matrix.clone_structure : nil

        # lapack_geev is a pure LAPACK routine so it expects column-major matrices,
        # so we need to transpose the input as well as the output.
        temporary_matrix = matrix.transpose
        NMatrix::LAPACK::lapack_geev(jobvl, # compute left eigenvectors of A?
                                     jobvr, # compute right eigenvectors of A? (left eigenvectors of A**T)
                                     n, # order of the matrix
                                     temporary_matrix,# input matrix (used as work)
                                     n, # leading dimension of matrix
                                     eigenvalues,# real part of computed eigenvalues
                                     imag_eigenvalues,# imag part of computed eigenvalues
                                     left_output,     # left eigenvectors, if applicable
                                     n, # leading dimension of left_output
                                     right_output,    # right eigenvectors, if applicable
                                     n, # leading dimension of right_output
                                     2*n)
        left_output = left_output.transpose if jobvl
        right_output = right_output.transpose if jobvr


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

        # This is a pure LAPACK function so it expects column-major functions.
        # So we need to transpose the input as well as the output.
        matrix = matrix.transpose
        NMatrix::LAPACK::lapack_gesvd(:a, :a, m, n, matrix, m, result[1], result[0], m, result[2], n, workspace_size)
        result[0] = result[0].transpose
        result[2] = result[2].transpose
        result
      end

      def gesdd(matrix, workspace_size=nil)
        min_workspace_size = matrix.shape.min * (6 + 4 * matrix.shape.min) + matrix.shape.max
        workspace_size = min_workspace_size if workspace_size.nil? || workspace_size < min_workspace_size

        result = alloc_svd_result(matrix)

        m = matrix.shape[0]
        n = matrix.shape[1]

        # This is a pure LAPACK function so it expects column-major functions.
        # So we need to transpose the input as well as the output.
        matrix = matrix.transpose
        NMatrix::LAPACK::lapack_gesdd(:a, m, n, matrix, m, result[1], result[0], m, result[2], n, workspace_size)
        result[0] = result[0].transpose
        result[2] = result[2].transpose
        result
      end
    end
  end

  def invert!
    raise(StorageTypeError, "invert only works on dense matrices currently") unless self.dense?
    raise(ShapeError, "Cannot invert non-square matrix") unless shape[0] == shape[1]
    raise(DataTypeError, "Cannot invert an integer matrix in-place") if self.integer_dtype?

    # Even though we are using the ATLAS plugin, we still might be missing
    # CLAPACK (and thus clapack_getri) if we are on OS X.
    if NMatrix.has_clapack?
      # Get the pivot array; factor the matrix
      # We can't used getrf! here since it doesn't have the clapack behavior,
      # so it doesn't play nicely with clapack_getri
      n = self.shape[0]
      pivot = NMatrix::LAPACK::clapack_getrf(:row, n, n, self, n)
      # Now calculate the inverse using the pivot array
      NMatrix::LAPACK::clapack_getri(:row, n, self, n, pivot)
      self
    else
      __inverse__(self,true)
    end
  end

  def potrf!(which)
    raise(StorageTypeError, "ATLAS functions only work on dense matrices") unless self.dense?
    raise(ShapeError, "Cholesky decomposition only valid for square matrices") unless self.dim == 2 && self.shape[0] == self.shape[1]

    NMatrix::LAPACK::clapack_potrf(:row, which, self.shape[0], self, self.shape[1])
  end
end
