# This is the main ruby file for the nmatrix-lapacke gem
require 'nmatrix' #need to have nmatrix required first or else bad things will happen
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

  def getrf!
    raise(StorageTypeError, "ATLAS functions only work on dense matrices") unless self.dense?

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
    raise(StorageTypeError, "ATLAS functions only work on dense matrices") unless self.dense?
    raise(ShapeError, "Cholesky decomposition only valid for square matrices") unless self.dim == 2 && self.shape[0] == self.shape[1]

    NMatrix::LAPACK::lapacke_potrf(:row, which, self.shape[0], self, self.shape[1])
  end

  def solve b
    raise(ShapeError, "b must be a column vector") unless b.dim == 2 && b.shape[1] == 1
    raise(ShapeError, "Must be called on square matrix") unless self.dim == 2 && self.shape[0] == self.shape[1]
    raise(ShapeError, "number of rows of b must equal number of cols of self") if 
      self.shape[1] != b.shape[0]
    raise ArgumentError, "only works with dense matrices" if self.stype != :dense
    raise ArgumentError, "only works for non-integer, non-object dtypes" if 
      integer_dtype? or object_dtype? or b.integer_dtype? or b.object_dtype?

    x     = b.clone
    clone = self.clone
    n = self.shape[0]
    ipiv = NMatrix::LAPACK.lapacke_getrf(:row, n, n, clone, n)
    NMatrix::LAPACK.lapacke_getrs(:row, :no_transpose, n, b.shape[1], clone, n, ipiv, x, b.shape[1])
    x
  end
end