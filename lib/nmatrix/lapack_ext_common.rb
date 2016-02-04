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
# == lapack_ext_common.rb
#
# Contains functions shared by nmatrix-atlas and nmatrix-lapacke gems.
#++

class NMatrix
  def NMatrix.register_lapack_extension(name)
    if (defined? @@lapack_extension)
      raise "Attempting to load #{name} when #{@@lapack_extension} is already loaded. You can only load one LAPACK extension."
    end

    @@lapack_extension = name
  end

  alias_method :internal_dot, :dot

  def dot(right_v)
    if (right_v.is_a?(NMatrix) && self.stype == :dense && right_v.stype == :dense &&
        self.dim == 2 && right_v.dim == 2 && self.shape[1] == right_v.shape[0])

      result_dtype = NMatrix.upcast(self.dtype,right_v.dtype)
      left = self.dtype == result_dtype ? self : self.cast(dtype: result_dtype)
      right = right_v.dtype == result_dtype ? right_v : right_v.cast(dtype: result_dtype)

      left = left.clone if left.is_ref?
      right = right.clone if right.is_ref?

      result_m = left.shape[0]
      result_n = right.shape[1]
      left_n = left.shape[1]
      vector = result_n == 1
      result = NMatrix.new([result_m,result_n], dtype: result_dtype)

      if vector
        NMatrix::BLAS.cblas_gemv(false, result_m, left_n, 1, left, left_n, right, 1, 0, result, 1)
      else
        NMatrix::BLAS.cblas_gemm(:row, false, false, result_m, result_n, left_n, 1, left, left_n, right, result_n, 0, result, result_n)
      end
      return result
    else
      #internal_dot will handle non-dense matrices (and also dot-products for NMatrix's with dim=1),
      #and also all error-handling if the input is not valid
      self.internal_dot(right_v)
    end
  end
end
