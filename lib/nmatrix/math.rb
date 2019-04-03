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
# == math.rb
#
# Math functionality for NMatrix, along with any NMatrix instance
# methods that correspond to ATLAS/BLAS/LAPACK functions (e.g.,
# laswp).
#++

class NMatrix
  module NMMath #:nodoc:
    METHODS_ARITY_2 = [:atan2, :ldexp, :hypot]
    METHODS_ARITY_1 = [:cos, :sin, :tan, :acos, :asin, :atan, :cosh, :sinh, :tanh, :acosh,
                       :asinh, :atanh, :exp, :log2, :log10, :sqrt, :cbrt, :erf, :erfc, :gamma, :-@,]
  end

  # Methods for generating permutation matrix from LU factorization results.
  module FactorizeLUMethods
    class << self
      def permutation_matrix_from(pivot_array)
        perm_arry = permutation_array_for(pivot_array)
        n         = NMatrix.zeros(perm_arry.size, dtype: :byte)

        perm_arry.each_with_index { |e, i| n[e, i] = 1 }

        n
      end

      def permutation_array_for(pivot_array)
        perm_arry = Array.new(pivot_array.size) { |i| i }
        perm_arry.each_index do |i|
          # the pivot indices returned by LAPACK getrf are indexed starting
          # from 1, so we need to subtract 1 here
          perm_arry[i], perm_arry[pivot_array[i] - 1] = perm_arry[pivot_array[i] - 1], perm_arry[i]
        end

        perm_arry
      end
    end
  end

  #
  # call-seq:
  #     invert! -> NMatrix
  #
  # Use LAPACK to calculate the inverse of the matrix (in-place) if available.
  # Only works on dense matrices. Alternatively uses in-place Gauss-Jordan
  # elimination.
  #
  # * *Raises* :
  #   - +StorageTypeError+ -> only implemented on dense matrices.
  #   - +ShapeError+ -> matrix must be square.
  #   - +DataTypeError+ -> cannot invert an integer matrix in-place.
  #
  def invert!
    raise(StorageTypeError, "invert only works on dense matrices currently") unless dense?
    raise(ShapeError, "Cannot invert non-square matrix") unless dim == 2 && shape[0] == shape[1]
    raise(DataTypeError, "Cannot invert an integer matrix in-place") if integer_dtype?

    # No internal implementation of getri, so use this other function
    __inverse__(self, true)
  end

  #
  # call-seq:
  #     invert -> NMatrix
  #
  # Make a copy of the matrix, then invert using Gauss-Jordan elimination.
  # Works without LAPACK.
  #
  # * *Returns* :
  #   - A dense NMatrix. Will be the same type as the input NMatrix,
  #   except if the input is an integral dtype, in which case it will be a
  #   :float64 NMatrix.
  #
  # * *Raises* :
  #   - +StorageTypeError+ -> only implemented on dense matrices.
  #   - +ShapeError+ -> matrix must be square.
  #
  def invert
    # write this in terms of invert! so plugins will only have to overwrite
    # invert! and not invert
    if integer_dtype?
      cloned = cast(dtype: :float64)
      cloned.invert!
    else
      cloned = clone
      cloned.invert!
    end
  end
  alias inverse invert

  # call-seq:
  #     exact_inverse! -> NMatrix
  #
  # Calulates inverse_exact of a matrix of size 2 or 3.
  # Only works on dense matrices.
  #
  # * *Raises* :
  #   - +DataTypeError+ -> cannot invert an integer matrix in-place.
  #   - +NotImplementedError+ -> cannot find exact inverse of matrix with size greater than 3  #
  def exact_inverse!
    raise(ShapeError, "Cannot invert non-square matrix") unless dim == 2 && shape[0] == shape[1]
    raise(DataTypeError, "Cannot invert an integer matrix in-place") if integer_dtype?
    # No internal implementation of getri, so use this other function
    n = shape[0]
    if n > 3
      raise(NotImplementedError, "Cannot find exact inverse of matrix of size greater than 3")
    else
      clond = clone
      __inverse_exact__(clond, n, n)
    end
  end

  #
  # call-seq:
  #     exact_inverse -> NMatrix
  #
  # Make a copy of the matrix, then invert using exact_inverse
  #
  # * *Returns* :
  #   - A dense NMatrix. Will be the same type as the input NMatrix,
  #   except if the input is an integral dtype, in which case it will be a
  #   :float64 NMatrix.
  #
  # * *Raises* :
  #   - +StorageTypeError+ -> only implemented on dense matrices.
  #   - +ShapeError+ -> matrix must be square.
  #   - +NotImplementedError+ -> cannot find exact inverse of matrix with size greater than 3
  #
  def exact_inverse
    # write this in terms of exact_inverse! so plugins will only have to overwrite
    # exact_inverse! and not exact_inverse
    if integer_dtype?
      cloned = cast(dtype: :float64)
      cloned.exact_inverse!
    else
      cloned = clone
      cloned.exact_inverse!
    end
  end
  alias invert_exactly exact_inverse

  #
  # call-seq:
  #     pinv -> NMatrix
  #
  # Compute the Moore-Penrose pseudo-inverse of a matrix using its
  # singular value decomposition (SVD).
  #
  # This function requires the nmatrix-atlas gem installed.
  #
  # * *Arguments* :
  #  - +tolerance(optional)+ -> Cutoff for small singular values.
  #
  # * *Returns* :
  #   -  Pseudo-inverse matrix.
  #
  # * *Raises* :
  #   - +NotImplementedError+ -> If called without nmatrix-atlas or nmatrix-lapacke gem.
  #   - +TypeError+ -> If called without float or complex data type.
  #
  # * *Examples* :
  #
  #  a = NMatrix.new([2,2],[1,2,
  #                         3,4], dtype: :float64)
  #  a.pinv # => [ [-2.0000000000000018, 1.0000000000000007]
  #                [1.5000000000000016, -0.5000000000000008] ]
  #
  #  b = NMatrix.new([4,1],[1,2,3,4], dtype: :float64)
  #  b.pinv # => [ [ 0.03333333, 0.06666667, 0.99999999, 0.13333333] ]
  #
  # == References
  #
  # * https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse
  # * G. Strang, Linear Algebra and Its Applications, 2nd Ed., Orlando, FL, Academic Press
  #
  def pinv(tolerance = 1e-15)
    raise DataTypeError, "pinv works only with matrices of float or complex data type" unless
      [:float32, :float64, :complex64, :complex128].include?(dtype)
    if complex_dtype?
      u, s, vt = complex_conjugate.gesvd # singular value decomposition
    else
      u, s, vt = gesvd
    end
    rows = shape[0]
    cols = shape[1]
    if rows < cols
      u_reduced = u
      vt_reduced = vt[0..rows - 1, 0..cols - 1].transpose
    else
      u_reduced = u[0..rows - 1, 0..cols - 1]
      vt_reduced = vt.transpose
    end
    largest_singular_value = s.max.to_f
    cutoff = tolerance * largest_singular_value
    (0...[rows, cols].min).each do |i|
      s[i] = 1 / s[i] if s[i] > cutoff
      s[i] = 0        if s[i] <= cutoff
    end
    multiplier = u_reduced.dot(NMatrix.diagonal(s.to_a)).transpose
    vt_reduced.dot(multiplier)
  end
  alias pseudo_inverse pinv
  alias pseudoinverse pinv

  #
  # call-seq:
  #     adjugate! -> NMatrix
  #
  # Calculate the adjugate of the matrix (in-place).
  # Only works on dense matrices.
  #
  # * *Raises* :
  #   - +StorageTypeError+ -> only implemented on dense matrices.
  #   - +ShapeError+ -> matrix must be square.
  #   - +DataTypeError+ -> cannot calculate adjugate of an integer matrix in-place.
  #
  def adjugate!
    raise(StorageTypeError, "adjugate only works on dense matrices currently") unless dense?
    raise(ShapeError, "Cannot calculate adjugate of a non-square matrix") unless dim == 2 && shape[0] == shape[1]
    raise(DataTypeError, "Cannot calculate adjugate of an integer matrix in-place") if integer_dtype?
    d = det
    invert!
    map! { |e| e * d }
    self
  end
  alias adjoint! adjugate!

  #
  # call-seq:
  #     adjugate -> NMatrix
  #
  # Make a copy of the matrix and calculate the adjugate of the matrix.
  # Only works on dense matrices.
  #
  # * *Returns* :
  #   - A dense NMatrix. Will be the same type as the input NMatrix,
  #   except if the input is an integral dtype, in which case it will be a
  #   :float64 NMatrix.
  #
  # * *Raises* :
  #   - +StorageTypeError+ -> only implemented on dense matrices.
  #   - +ShapeError+ -> matrix must be square.
  #
  def adjugate
    raise(StorageTypeError, "adjugate only works on dense matrices currently") unless dense?
    raise(ShapeError, "Cannot calculate adjugate of a non-square matrix") unless dim == 2 && shape[0] == shape[1]
    d = det
    mat = invert
    mat.map! { |e| e * d }
    mat
  end
  alias adjoint adjugate

  # Reduce self to upper hessenberg form using householder transforms.
  #
  # == References
  #
  # * http://en.wikipedia.org/wiki/Hessenberg_matrix
  # * http://www.mymathlib.com/c_source/matrices/eigen/hessenberg_orthog.c
  def hessenberg
    clone.hessenberg!
  end

  # Destructive version of #hessenberg
  def hessenberg!
    raise ShapeError, "Trying to reduce non 2D matrix to hessenberg form" if
      shape.size != 2
    raise ShapeError, "Trying to reduce non-square matrix to hessenberg form" if
      shape[0] != shape[1]
    raise StorageTypeError, "Matrix must be dense" if stype != :dense
    raise TypeError, "Works with float matrices only" unless
      [:float64, :float32].include?(dtype)

    __hessenberg__(self)
    self
  end

  # call-seq:
  #   matrix_norm -> Numeric
  #
  #  Calculates the selected norm (defaults to 2-norm) of a 2D matrix.
  #
  #  This should be used for small or medium sized matrices.
  #  For greater matrices, there should be a separate implementation where
  #  the norm is estimated rather than computed, for the sake of computation speed.
  #
  #  Currently implemented norms are 1-norm, 2-norm, Frobenius, Infinity.
  #  A minus on the 1, 2 and inf norms returns the minimum instead of the maximum value.
  #
  #  Tested mainly with dense matrices. Further checks and modifications might
  #  be necessary for sparse matrices.
  #
  # * *Returns* :
  # - The selected norm of the matrix.
  # * *Raises* :
  # - +NotImplementedError+ -> norm can be calculated only for 2D matrices
  # - +ArgumentError+ -> unrecognized norm
  #
  def matrix_norm type = 2
    raise(NotImplementedError, "norm can be calculated only for 2D matrices") unless dim == 2
    raise(NotImplementedError, "norm only implemented for dense storage") unless stype == :dense
    raise(ArgumentError, "norm not defined for byte dtype") if dtype == :byte
    case type
    when nil, 2, -2
      return two_matrix_norm (type == -2)
    when 1, -1
      return one_matrix_norm (type == -1)
    when :frobenius, :fro
      return fro_matrix_norm
    when :infinity, :inf, :'-inf', :'-infinity'
      return inf_matrix_norm (type == :'-inf' || type == :'-infinity')
    else
      raise ArgumentError.new("argument must be a valid integer or symbol")
    end
  end

  # Calculate the variance co-variance matrix
  #
  # == Options
  #
  # * +:for_sample_data+ - Default true. If set to false will consider the denominator for
  #   population data (i.e. N, as opposed to N-1 for sample data).
  #
  # == References
  #
  # * http://stattrek.com/matrix-algebra/covariance-matrix.aspx
  def cov(opts = {})
    raise TypeError, "Only works for non-integer dtypes" if integer_dtype?
    opts = {
      for_sample_data: true,
    }.merge(opts)

    denominator      = opts[:for_sample_data] ? rows - 1 : rows
    ones             = NMatrix.ones [rows, 1]
    deviation_scores = self - ones.dot(ones.transpose).dot(self) / rows
    deviation_scores.transpose.dot(deviation_scores) / denominator
  end

  # Calculate the correlation matrix.
  def corr
    raise NotImplementedError, "Does not work for complex dtypes" if complex_dtype?
    standard_deviation = std
    cov / standard_deviation.transpose.dot(standard_deviation)
  end

  # Raise a square matrix to a power. Be careful of numeric overflows!
  # In case *n* is 0, an identity matrix of the same dimension is returned. In case
  # of negative *n*, the matrix is inverted and the absolute value of *n* taken
  # for computing the power.
  #
  # == Arguments
  #
  # * +n+ - Integer to which self is to be raised.
  #
  # == References
  #
  # * R.G Dromey - How to Solve it by Computer. Link -
  #     http://www.amazon.com/Solve-Computer-Prentice-Hall-International-Science/dp/0134340019/ref=sr_1_1?ie=UTF8&qid=1422605572&sr=8-1&keywords=how+to+solve+it+by+computer
  def pow n
    raise ShapeError, "Only works with 2D square matrices." if
      (shape[0] != shape[1]) || (shape.size != 2)
    raise TypeError, "Only works with integer powers" unless n.is_a?(Integer)

    sequence = (integer_dtype? ? cast(dtype: :int64) : self).clone
    product  = NMatrix.eye shape[0], dtype: sequence.dtype, stype: sequence.stype

    if n == 0
      return NMatrix.eye(shape, dtype: dtype, stype: stype)
    elsif n == 1
      return sequence
    elsif n < 0
      n = n.abs
      sequence.invert!
      product = NMatrix.eye shape[0], dtype: sequence.dtype, stype: sequence.stype
    end

    # Decompose n to reduce the number of multiplications.
    while n > 0
      product = product.dot(sequence) if n % 2 == 1
      n /= 2
      sequence = sequence.dot(sequence)
    end

    product
  end

  # Compute the Kronecker product of +self+ and other NMatrix
  #
  # === Arguments
  #
  #   * +mat+ - A 2D NMatrix object
  #
  # === Usage
  #
  #  a = NMatrix.new([2,2],[1,2,
  #                         3,4])
  #  b = NMatrix.new([2,3],[1,1,1,
  #                         1,1,1], dtype: :float64)
  #  a.kron_prod(b) # => [ [1.0, 1.0, 1.0, 2.0, 2.0, 2.0]
  #                        [1.0, 1.0, 1.0, 2.0, 2.0, 2.0]
  #                        [3.0, 3.0, 3.0, 4.0, 4.0, 4.0]
  #                        [3.0, 3.0, 3.0, 4.0, 4.0, 4.0] ]
  #
  def kron_prod(mat)
    unless (dimensions == 2) && (mat.dimensions == 2)
      raise ShapeError, "Implemented for 2D NMatrix objects only."
    end

    # compute the shape [n,m] of the product matrix
    n, m = shape[0] * mat.shape[0], shape[1] * mat.shape[1]
    # compute the entries of the product matrix
    kron_prod_array = []
    if yale?
      # +:yale+ requires to get the row by copy in order to apply +#transpose+ to it
      each_row(getby = :copy) do |selfr|
        mat.each_row do |matr|
          kron_prod_array += (selfr.transpose.dot matr).to_flat_a
        end
      end
    else
      each_row do |selfr|
        mat.each_row do |matr|
          kron_prod_array += (selfr.transpose.dot matr).to_flat_a
        end
      end
    end

    NMatrix.new([n, m], kron_prod_array)
  end

  #
  # call-seq:
  #     trace -> Numeric
  #
  # Calculates the trace of an nxn matrix.
  #
  # * *Raises* :
  #   - +ShapeError+ -> Expected square matrix
  #
  # * *Returns* :
  #   - The trace of the matrix (a numeric value)
  #
  def trace
    raise(ShapeError, "Expected square matrix") unless shape[0] == shape[1] && dim == 2

    (0...shape[0]).inject(0) do |total, i|
      total + self[i, i]
    end
  end

  ##
  # call-seq:
  #   mean() -> NMatrix
  #   mean(dimen) -> NMatrix
  #
  # Calculates the mean along the specified dimension.
  #
  # This will force integer types to float64 dtype.
  #
  # @see #inject_rank
  #
  def mean(dimen = 0)
    reduce_dtype = nil
    if integer_dtype?
      reduce_dtype = :float64
    end
    inject_rank(dimen, 0.0, reduce_dtype) { |mean, sub_mat|
      mean + sub_mat
    } / shape[dimen]
  end

  ##
  # call-seq:
  #   sum() -> NMatrix
  #   cumsum() -> NMatrix
  #   sum(dimen) -> NMatrix
  #   cumsum(dimen) -> NMatrix
  #
  # Calculates the sum along the specified dimension.
  #
  # @see #inject_rank
  def sum(dimen = 0)
    inject_rank(dimen, 0.0) do |sum, sub_mat|
      sum + sub_mat
    end
  end
  alias cumsum sum

  ##
  # call-seq:
  #   min() -> NMatrix
  #   min(dimen) -> NMatrix
  #
  # Calculates the minimum along the specified dimension.
  #
  # @see #inject_rank
  #
  def min(dimen = 0)
    inject_rank(dimen) do |min, sub_mat|
      if min.is_a? NMatrix
        min * (min <= sub_mat).cast(stype, dtype) + (min * 0.0 + (min > sub_mat).cast(stype, dtype)) * sub_mat
      else
        min <= sub_mat ? min : sub_mat
      end
    end
  end

  ##
  # call-seq:
  #   max() -> NMatrix
  #   max(dimen) -> NMatrix
  #
  # Calculates the maximum along the specified dimension.
  #
  # @see #inject_rank
  #
  def max(dimen = 0)
    inject_rank(dimen) do |max, sub_mat|
      if max.is_a? NMatrix
        max * (max >= sub_mat).cast(stype, dtype) + (max * 0.0 + (max < sub_mat).cast(stype, dtype)) * sub_mat
      else
        max >= sub_mat ? max : sub_mat
      end
    end
  end

  ##
  # call-seq:
  #   variance() -> NMatrix
  #   variance(dimen) -> NMatrix
  #
  # Calculates the sample variance along the specified dimension.
  #
  # This will force integer types to float64 dtype.
  #
  # @see #inject_rank
  #
  def variance(dimen = 0)
    reduce_dtype = nil
    if integer_dtype?
      reduce_dtype = :float64
    end
    m = mean(dimen)
    inject_rank(dimen, 0.0, reduce_dtype) do |var, sub_mat|
      var + (m - sub_mat) * (m - sub_mat) / (shape[dimen] - 1)
    end
  end

  ##
  # call-seq:
  #   std() -> NMatrix
  #   std(dimen) -> NMatrix
  #
  #
  # Calculates the sample standard deviation along the specified dimension.
  #
  # This will force integer types to float64 dtype.
  #
  # @see #inject_rank
  #
  def std(dimen = 0)
    variance(dimen).sqrt
  end

  #
  # call-seq:
  #     abs_dtype -> Symbol
  #
  # Returns the dtype of the result of a call to #abs. In most cases, this is the same as dtype; it should only differ
  # for :complex64 (where it's :float32) and :complex128 (:float64).
  def abs_dtype
    if dtype == :complex64
      :float32
    elsif dtype == :complex128
      :float64
    else
      dtype
    end
  end

  #
  # call-seq:
  #     abs -> NMatrix
  #
  # Maps all values in a matrix to their absolute values.
  def abs
    if stype == :dense
      __dense_map__ { |v| v.abs }
    elsif stype == :list
      # FIXME: Need __list_map_stored__, but this will do for now.
      __list_map_merged_stored__(nil, nil) { |v, dummy| v.abs }
    else
      __yale_map_stored__ { |v| v.abs }
    end.cast(stype, abs_dtype)
  end

  # Norm calculation methods
  # Frobenius norm: the Euclidean norm of the matrix, treated as if it were a vector
  def fro_matrix_norm
    # float64 has to be used in any case, since nrm2 will not yield correct result for float32
    self_cast = cast(dtype: :float64)

    column_vector = self_cast.reshape([size, 1])

    column_vector.nrm2
  end

  # 2-norm: the largest/smallest singular value of the matrix
  def two_matrix_norm minus = false
    self_cast = cast(dtype: :float64)

    # TODO: confirm if this is the desired svd calculation
    svd = self_cast.gesvd
    return svd[1][0, 0] unless minus
    svd[1][svd[1].rows - 1, svd[1].cols - 1]
  end

  # 1-norm: the maximum/minimum absolute column sum of the matrix
  def one_matrix_norm minus = false
    # TODO: change traversing method for sparse matrices
    number_of_columns = cols
    col_sums = []

    number_of_columns.times do |i|
      col_sums << col(i).inject(0) { |sum, number| sum += number.abs}
    end

    return col_sums.max unless minus
    col_sums.min
  end

  # Infinity norm: the maximum/minimum absolute row sum of the matrix
  def inf_matrix_norm minus = false
    number_of_rows = rows
    row_sums = []

    number_of_rows.times do |i|
      row_sums << row(i).inject(0) { |sum, number| sum += number.abs}
    end

    return row_sums.max unless minus
    row_sums.min
  end

  #
  # call-seq:
  #     positive_definite? -> boolean
  #
  # A matrix is positive definite if it’s symmetric and all its eigenvalues are positive
  #
  # * *Returns* :
  #   - A boolean value telling if the NMatrix is positive definite or not.
  # * *Raises* :
  #   - +ShapeError+ -> Must be used on square matrices.
  #
  def positive_definite?
    raise(ShapeError, "positive definite calculated only for square matrices") unless
      dim == 2 && shape[0] == shape[1]
    cond = 0
    while cond != cols
      if self[0..cond, 0..cond].det <= 0
        return false
      end
      cond += 1
    end
    true
  end

  #
  # call-seq:
  #   svd_rank() -> int
  #   svd_rank(tolerence) ->int
  # Gives rank of the matrix based on the singular value decomposition.
  # The rank of a matrix  is computed as the number of diagonal elements in Sigma that are larger than a tolerance
  #
  # * *Returns* :
  # - An integer equal to the rank of the matrix
  # * *Raises* :
  #  - +ShapeError+ -> Is only computable on 2-D matrices
  #
  def svd_rank(tolerence = "default")
    raise(ShapeError, "rank calculated only for 2-D matrices") unless
      dim == 2

    sigmas = gesvd[1].to_a.flatten
    eps = NMatrix::FLOAT64_EPSILON

    # epsilon depends on the width of the number
    if dtype == :float32 || dtype == :complex64
      eps = NMatrix::FLOAT32_EPSILON
    end
    case tolerence
      when "default"
        tolerence = shape.max * sigmas.max * eps # tolerence of a Matrix A is max(size(A))*eps(norm(A)). norm(A) is nearly equal to max(sigma of A)
    end
    sigmas.map { |x| x > tolerence ? 1 : 0 }.reduce(:+)
  end

  protected

  # Define the element-wise operations for lists. Note that the __list_map_merged_stored__ iterator returns a Ruby Object
  # matrix, which we then cast back to the appropriate type. If you don't want that, you can redefine these functions in
  # your own code.
  {add: :+, sub: :-, mul: :*, div: :/, pow: :**, mod: :%}.each_pair do |ewop, op|
    define_method("__list_elementwise_#{ewop}__") do |rhs|
      __list_map_merged_stored__(rhs, nil) { |l, r| l.send(op, r) }.cast(stype, NMatrix.upcast(dtype, rhs.dtype))
    end
    define_method("__dense_elementwise_#{ewop}__") do |rhs|
      __dense_map_pair__(rhs) { |l, r| l.send(op, r) }.cast(stype, NMatrix.upcast(dtype, rhs.dtype))
    end
    define_method("__yale_elementwise_#{ewop}__") do |rhs|
      __yale_map_merged_stored__(rhs, nil) { |l, r| l.send(op, r) }.cast(stype, NMatrix.upcast(dtype, rhs.dtype))
    end
    define_method("__list_scalar_#{ewop}__") do |rhs|
      __list_map_merged_stored__(rhs, nil) { |l, r| l.send(op, r) }.cast(stype, NMatrix.upcast(dtype, NMatrix.min_dtype(rhs)))
    end
    define_method("__yale_scalar_#{ewop}__") do |rhs|
      __yale_map_stored__ { |l| l.send(op, rhs) }.cast(stype, NMatrix.upcast(dtype, NMatrix.min_dtype(rhs)))
    end
    define_method("__dense_scalar_#{ewop}__") do |rhs|
      __dense_map__ { |l| l.send(op, rhs) }.cast(stype, NMatrix.upcast(dtype, NMatrix.min_dtype(rhs)))
    end
  end

  # These don't actually take an argument -- they're called reverse-polish style on the matrix.
  # This group always gets casted to float64.
  [:log, :log2, :log10, :sqrt, :sin, :cos, :tan, :acos, :asin, :atan, :cosh, :sinh, :tanh, :acosh,
   :asinh, :atanh, :exp, :erf, :erfc, :gamma, :cbrt, :round,].each do |ewop|
    define_method("__list_unary_#{ewop}__") do
      __list_map_stored__(nil) { |l| Math.send(ewop, l) }.cast(stype, NMatrix.upcast(dtype, :float64))
    end
    define_method("__yale_unary_#{ewop}__") do
      __yale_map_stored__ { |l| Math.send(ewop, l) }.cast(stype, NMatrix.upcast(dtype, :float64))
    end
    define_method("__dense_unary_#{ewop}__") do
      __dense_map__ { |l| Math.send(ewop, l) }.cast(stype, NMatrix.upcast(dtype, :float64))
    end
  end

  #:stopdoc:
  # log takes an optional single argument, the base. Default to natural log.
  def __list_unary_log__(base)
    __list_map_stored__(nil) { |l| Math.log(l, base) }.cast(stype, NMatrix.upcast(dtype, :float64))
  end

  def __yale_unary_log__(base)
    __yale_map_stored__ { |l| Math.log(l, base) }.cast(stype, NMatrix.upcast(dtype, :float64))
  end

  def __dense_unary_log__(base)
    __dense_map__ { |l| Math.log(l, base) }.cast(stype, NMatrix.upcast(dtype, :float64))
  end

  # These are for negating matrix contents using -@
  def __list_unary_negate__
    __list_map_stored__(nil) { |l| -l }.cast(stype, dtype)
  end

  def __yale_unary_negate__
    __yale_map_stored__ { |l| -l }.cast(stype, dtype)
  end

  def __dense_unary_negate__
    __dense_map__ { |l| -l }.cast(stype, dtype)
  end
  #:startdoc:

  # These are for rounding each value of a matrix. Takes an optional argument
  def __list_unary_round__(precision)
    if complex_dtype?
      __list_map_stored__(nil) { |l| Complex(l.real.round(precision), l.imag.round(precision)) }
        .cast(stype, dtype)
    else
      __list_map_stored__(nil) { |l| l.round(precision) }.cast(stype, dtype)
    end
  end

  def __yale_unary_round__(precision)
    if complex_dtype?
      __yale_map_stored__ { |l| Complex(l.real.round(precision), l.imag.round(precision)) }
        .cast(stype, dtype)
    else
      __yale_map_stored__ { |l| l.round(precision) }.cast(stype, dtype)
    end
  end

  def __dense_unary_round__(precision)
    if complex_dtype?
      __dense_map__ { |l| Complex(l.real.round(precision), l.imag.round(precision)) }
        .cast(stype, dtype)
    else
      __dense_map__ { |l| l.round(precision) }.cast(stype, dtype)
    end
  end

  # These are for calculating the floor or ceil of matrix
  def dtype_for_floor_or_ceil
    if integer_dtype? || [:complex64, :complex128, :object].include?(dtype)
      return_dtype = dtype
    elsif [:float32, :float64].include?(dtype)
      return_dtype = :int64
    end

    return_dtype
  end

  [:floor, :ceil].each do |meth|
    define_method("__list_unary_#{meth}__") do
      return_dtype = dtype_for_floor_or_ceil

      if [:complex64, :complex128].include?(dtype)
        __list_map_stored__(nil) { |l| Complex(l.real.send(meth), l.imag.send(meth)) }.cast(stype, return_dtype)
      else
        __list_map_stored__(nil) { |l| l.send(meth) }.cast(stype, return_dtype)
      end
    end

    define_method("__yale_unary_#{meth}__") do
      return_dtype = dtype_for_floor_or_ceil

      if [:complex64, :complex128].include?(dtype)
        __yale_map_stored__ { |l| Complex(l.real.send(meth), l.imag.send(meth)) }.cast(stype, return_dtype)
      else
        __yale_map_stored__ { |l| l.send(meth) }.cast(stype, return_dtype)
      end
    end

    define_method("__dense_unary_#{meth}__") do
      return_dtype = dtype_for_floor_or_ceil

      if [:complex64, :complex128].include?(dtype)
        __dense_map__ { |l| Complex(l.real.send(meth), l.imag.send(meth)) }.cast(stype, return_dtype)
      else
        __dense_map__ { |l| l.send(meth) }.cast(stype, return_dtype)
      end
    end
  end

  # These take two arguments. One might be a matrix, and one might be a scalar.
  # See also monkeys.rb, which contains Math module patches to let the first
  # arg be a scalar
  [:atan2, :ldexp, :hypot].each do |ewop|
    define_method("__list_elementwise_#{ewop}__") do |rhs, order|
      if order
        __list_map_merged_stored__(rhs, nil) { |r, l| Math.send(ewop, l, r) }
      else
        __list_map_merged_stored__(rhs, nil) { |l, r| Math.send(ewop, l, r) }
      end.cast(stype, NMatrix.upcast(dtype, :float64))
    end

    define_method("__dense_elementwise_#{ewop}__") do |rhs, order|
      if order
        __dense_map_pair__(rhs) { |r, l| Math.send(ewop, l, r) }
      else
        __dense_map_pair__(rhs) { |l, r| Math.send(ewop, l, r) }
      end.cast(stype, NMatrix.upcast(dtype, :float64))
    end

    define_method("__yale_elementwise_#{ewop}__") do |rhs, order|
      if order
        __yale_map_merged_stored__(rhs, nil) { |r, l| Math.send(ewop, l, r) }
      else
        __yale_map_merged_stored__(rhs, nil) { |l, r| Math.send(ewop, l, r) }
      end.cast(stype, NMatrix.upcast(dtype, :float64))
    end

    define_method("__list_scalar_#{ewop}__") do |rhs, order|
      if order
        __list_map_stored__(nil) { |l| Math.send(ewop, rhs, l) }
      else
        __list_map_stored__(nil) { |l| Math.send(ewop, l, rhs) }
      end.cast(stype, NMatrix.upcast(dtype, :float64))
    end

    define_method("__yale_scalar_#{ewop}__") do |rhs, order|
      if order
        __yale_map_stored__ { |l| Math.send(ewop, rhs, l) }
      else
        __yale_map_stored__ { |l| Math.send(ewop, l, rhs) }
      end.cast(stype, NMatrix.upcast(dtype, :float64))
    end

    define_method("__dense_scalar_#{ewop}__") do |rhs, order|
      if order
        __dense_map__ { |l| Math.send(ewop, rhs, l) }
      else
        __dense_map__ { |l| Math.send(ewop, l, rhs) }
      end.cast(stype, NMatrix.upcast(dtype, :float64))
    end
  end

  # Equality operators do not involve a cast. We want to get back matrices of TrueClass and FalseClass.
  {eqeq: :==, neq: :!=, lt: :<, gt: :>, leq: :<=, geq: :>=}.each_pair do |ewop, op|
    define_method("__list_elementwise_#{ewop}__") do |rhs|
      __list_map_merged_stored__(rhs, nil) { |l, r| l.send(op, r) }
    end
    define_method("__dense_elementwise_#{ewop}__") do |rhs|
      __dense_map_pair__(rhs) { |l, r| l.send(op, r) }
    end
    define_method("__yale_elementwise_#{ewop}__") do |rhs|
      __yale_map_merged_stored__(rhs, nil) { |l, r| l.send(op, r) }
    end

    define_method("__list_scalar_#{ewop}__") do |rhs|
      __list_map_merged_stored__(rhs, nil) { |l, r| l.send(op, r) }
    end
    define_method("__yale_scalar_#{ewop}__") do |rhs|
      __yale_map_stored__ { |l| l.send(op, rhs) }
    end
    define_method("__dense_scalar_#{ewop}__") do |rhs|
      __dense_map__ { |l| l.send(op, rhs) }
    end
  end
end

if jruby?
  require_relative "./jruby/math.rb"
else
  require_relative "./cruby/math.rb"
end
