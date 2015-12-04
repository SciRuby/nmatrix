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
      :asinh, :atanh, :exp, :log2, :log10, :sqrt, :cbrt, :erf, :erfc, :gamma, :-@]
  end

  # Methods for generating permutation matrix from LU factorization results.
  module FactorizeLUMethods
    class << self
      def permutation_matrix_from pivot_array
        perm_arry = permutation_array_for pivot_array
        n         = NMatrix.zeros perm_arry.size, dtype: :byte

        perm_arry.each_with_index { |e, i| n[i,e] = 1 }

        n
      end

      def permutation_array_for pivot_array
        perm_arry = Array.new(pivot_array.size) { |i| i }
        perm_arry.each_index do |i|
          perm_arry[i], perm_arry[pivot_array[i]] = perm_arry[pivot_array[i]], perm_arry[i]
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
  def invert!
    if NMatrix.has_clapack?
      # Get the pivot array; factor the matrix
      pivot = self.getrf!

      # Now calculate the inverse using the pivot array
      NMatrix::LAPACK::clapack_getri(:row, self.shape[1], self, self.shape[1], pivot)

      self
    else
      if self.integer_dtype?
        __inverse__(self.cast(dtype: :float64), true)
      else
        dtype = self.dtype
        __inverse__(self, true)
      end
    end
  end

  #
  # call-seq:
  #     invert -> NMatrix
  #
  # Make a copy of the matrix, then invert using Gauss-Jordan elimination.
  # Works without LAPACK.
  #
  #
  # * *Returns* :
  #   - A dense NMatrix.
  #
  def invert lda=nil, ldb=nil
    if lda.nil? and ldb.nil?
      if NMatrix.has_clapack?
        begin
          self.cast(:dense, self.dtype).invert! # call CLAPACK version
        rescue NotImplementedError
          inverse = self.clone
          __inverse__(inverse, false)
        end
      elsif self.integer_dtype? # FIXME: This check is probably too slow.
        casted = self.cast(dtype: :float64)
        inverse       = casted.clone
        casted.__inverse__(inverse, false)
      else
        inverse       = self.clone
        __inverse__(inverse, false)
      end
    else
      inverse = self.clone_structure
      if self.integer_dtype?
        __inverse_exact__(inverse.cast(dtype: :float64), lda, ldb)
      else
        __inverse_exact__(inverse, lda, ldb)
      end
    end
  end
  alias :inverse :invert

  #
  # call-seq:
  #     getrf! -> NMatrix
  #
  # LU factorization of a general M-by-N matrix +A+ using partial pivoting with
  # row interchanges. Only works in dense matrices.
  #
  # * *Returns* :
  #   - The IPIV vector. The L and U matrices are stored in A.
  # * *Raises* :
  #   - +StorageTypeError+ -> ATLAS functions only work on dense matrices.
  #
  def getrf!
    raise(StorageTypeError, "ATLAS functions only work on dense matrices") unless self.dense?
    NMatrix::LAPACK::clapack_getrf(:row, self.shape[0], self.shape[1], self, self.shape[1])
  end

  alias :lu_decomposition! :getrf!

  #
  # call-seq:
  #     getrf -> NMatrix
  #
  # In-place version of #getrf!. Returns the new matrix, which contains L and U matrices.
  #
  # * *Raises* :
  #   - +StorageTypeError+ -> ATLAS functions only work on dense matrices.
  #
  def getrf
    a = self.clone
    a.getrf!
    return a
  end


  #
  # call-seq:
  #     potrf!(upper_or_lower) -> NMatrix
  #
  # Cholesky factorization of a symmetric positive-definite matrix -- or, if complex,
  # a Hermitian positive-definite matrix +A+. This uses the ATLAS function clapack_potrf,
  # so the result will be written in either the upper or lower triangular portion of the
  # matrix upon which it is called.
  #
  # * *Returns* :
  #   the triangular portion specified by the parameter
  # * *Raises* :
  #   - +StorageTypeError+ -> ATLAS functions only work on dense matrices.
  #
  def potrf!(which)
    raise(StorageTypeError, "ATLAS functions only work on dense matrices") unless self.dense?
    # FIXME: Surely there's an easy way to calculate one of these from the other. Do we really need to run twice?
    NMatrix::LAPACK::clapack_potrf(:row, which, self.shape[0], self, self.shape[1])
  end

  def potrf_upper!
    potrf! :upper
  end

  def potrf_lower!
    potrf! :lower
  end


  #
  # call-seq:
  #     factorize_cholesky -> ...
  #
  # Cholesky factorization of a matrix.
  def factorize_cholesky
    [self.clone.potrf_upper!.triu!,
    self.clone.potrf_lower!.tril!]
  end

  #
  # call-seq:
  #     factorize_lu -> ...
  #
  # LU factorization of a matrix. Optionally return the permutation matrix.
  #   Note that computing the permutation matrix will introduce a slight memory
  #   and time overhead. 
  # 
  # == Arguments
  # 
  # +with_permutation_matrix+ - If set to *true* will return the permutation 
  #   matrix alongwith the LU factorization as a second return value.
  # 
  # FIXME: For some reason, getrf seems to require that the matrix be transposed first -- and then you have to transpose the
  # FIXME: result again. Ideally, this would be an in-place factorize instead, and would be called nm_factorize_lu_bang.
  #
  def factorize_lu with_permutation_matrix=nil
    raise(NotImplementedError, "only implemented for dense storage") unless self.stype == :dense
    raise(NotImplementedError, "matrix is not 2-dimensional") unless self.dimensions == 2

    t     = self.transpose
    pivot = NMatrix::LAPACK::clapack_getrf(:row, t.shape[0], t.shape[1], t, t.shape[1])
    return t.transpose unless with_permutation_matrix

    [t.transpose, FactorizeLUMethods.permutation_matrix_from(pivot)]
  end

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
      [:float64,:float32].include?(dtype)

    __hessenberg__(self)
    self
  end

  # Solve a system of linear equations where *self* is the matrix of co-efficients
  # and *b* is the vertical vector of right hand sides. Only works with dense
  # matrices and non-integer, non-object data types.
  # 
  # == Arguments
  # 
  # +b+ - Vector of Right Hand Sides.
  # 
  # == Usage
  # 
  #   a = NMatrix.new [2,2], [3,1,1,2], dtype: dtype
  #   b = NMatrix.new [2,1], [9,8], dtype: dtype
  #   a.solve(b)
  def solve b
    raise ArgumentError, "b must be a column vector" if b.shape[1] != 1
    raise ArgumentError, "number of rows of b must equal number of cols of self" if 
      self.shape[1] != b.shape[0]
    raise ArgumentError, "only works with dense matrices" if self.stype != :dense
    raise ArgumentError, "only works for non-integer, non-object dtypes" if 
      integer_dtype? or object_dtype? or b.integer_dtype? or b.object_dtype?

    x     = b.clone_structure
    clone = self.clone
    t     = clone.transpose # transpose because of the getrf anomaly described above.
    pivot = t.lu_decomposition!
    t     = t.transpose
    
    __solve__(t, b, x, pivot)
    x
  end

  #
  # call-seq:
  #     gesvd! -> [u, sigma, v_transpose]
  #     gesvd! -> [u, sigma, v_conjugate_transpose] # complex
  #
  # Compute the singular value decomposition of a matrix using LAPACK's GESVD function.
  # This is destructive, modifying the source NMatrix.  See also #gesdd.
  #
  # Optionally accepts a +workspace_size+ parameter, which will be honored only if it is larger than what LAPACK
  # requires.
  #
  def gesvd!(workspace_size=1)
    NMatrix::LAPACK::gesvd(self, workspace_size)
  end

  #
  # call-seq:
  #     gesvd -> [u, sigma, v_transpose]
  #     gesvd -> [u, sigma, v_conjugate_transpose] # complex
  #
  # Compute the singular value decomposition of a matrix using LAPACK's GESVD function.
  #
  # Optionally accepts a +workspace_size+ parameter, which will be honored only if it is larger than what LAPACK
  # requires.
  #
  def gesvd(workspace_size=1)
    self.clone.gesvd!(workspace_size)
  end



  #
  # call-seq:
  #     gesdd! -> [u, sigma, v_transpose]
  #     gesdd! -> [u, sigma, v_conjugate_transpose] # complex
  #
  # Compute the singular value decomposition of a matrix using LAPACK's GESDD function. This uses a divide-and-conquer
  # strategy. This is destructive, modifying the source NMatrix.  See also #gesvd.
  #
  # Optionally accepts a +workspace_size+ parameter, which will be honored only if it is larger than what LAPACK
  # requires.
  #
  def gesdd!(workspace_size=nil)
    NMatrix::LAPACK::gesdd(self, workspace_size)
  end

  #
  # call-seq:
  #     gesdd -> [u, sigma, v_transpose]
  #     gesdd -> [u, sigma, v_conjugate_transpose] # complex
  #
  # Compute the singular value decomposition of a matrix using LAPACK's GESDD function. This uses a divide-and-conquer
  # strategy. See also #gesvd.
  #
  # Optionally accepts a +workspace_size+ parameter, which will be honored only if it is larger than what LAPACK
  # requires.
  #
  def gesdd(workspace_size=nil)
    self.clone.gesdd!(workspace_size)
  end

  #
  # call-seq:
  #     laswp!(ary) -> NMatrix
  #
  # In-place permute the columns of a dense matrix using LASWP according to the order given as an array +ary+.
  #
  # If +:convention+ is +:lapack+, then +ary+ represents a sequence of pair-wise permutations which are 
  # performed successively. That is, the i'th entry of +ary+ is the index of the column to swap 
  # the i'th column with, having already applied all earlier swaps. 
  #
  # If +:convention+ is +:intuitive+, then +ary+ represents the order of columns after the permutation. 
  # That is, the i'th entry of +ary+ is the index of the column that will be in position i after the 
  # reordering (Matlab-like behaviour). This is the default.
  #
  # Not yet implemented for yale or list. 
  #
  # == Arguments
  #
  # * +ary+ - An Array specifying the order of the columns. See above for details.
  # 
  # == Options
  # 
  # * +:covention+ - Possible values are +:lapack+ and +:intuitive+. Default is +:intuitive+. See above for details.
  #
  def laswp!(ary, opts={})
    raise(StorageTypeError, "ATLAS functions only work on dense matrices") unless self.dense?
    opts = { convention: :intuitive }.merge(opts)
    
    if opts[:convention] == :intuitive
      if ary.length != ary.uniq.length
        raise(ArgumentError, "No duplicated entries in the order array are allowed under convention :intuitive")
      end
      n = self.shape[1]
      p = []
      order = (0...n).to_a
      0.upto(n-2) do |i|
        p[i] = order.index(ary[i])
        order[i], order[p[i]] = order[p[i]], order[i]
      end
      p[n-1] = n-1
    else
      p = ary
    end

    NMatrix::LAPACK::laswp(self, p)
  end

  #
  # call-seq:
  #     laswp(ary) -> NMatrix
  #
  # Permute the columns of a dense matrix using LASWP according to the order given in an array +ary+.
  #
  # If +:convention+ is +:lapack+, then +ary+ represents a sequence of pair-wise permutations which are 
  # performed successively. That is, the i'th entry of +ary+ is the index of the column to swap 
  # the i'th column with, having already applied all earlier swaps. This is the default.
  #
  # If +:convention+ is +:intuitive+, then +ary+ represents the order of columns after the permutation. 
  # That is, the i'th entry of +ary+ is the index of the column that will be in position i after the 
  # reordering (Matlab-like behaviour). 
  #
  # Not yet implemented for yale or list. 
  #
  # == Arguments
  #
  # * +ary+ - An Array specifying the order of the columns. See above for details.
  # 
  # == Options
  # 
  # * +:covention+ - Possible values are +:lapack+ and +:intuitive+. Default is +:lapack+. See above for details.
  #
  def laswp(ary, opts={})
    self.clone.laswp!(ary, opts)
  end

  #
  # call-seq:
  #     det -> determinant
  #
  # Calculate the determinant by way of LU decomposition. This is accomplished
  # using clapack_getrf, and then by summing the diagonal elements. There is a
  # risk of underflow/overflow.
  #
  # There are probably also more efficient ways to calculate the determinant.
  # This method requires making a copy of the matrix, since clapack_getrf
  # modifies its input.
  #
  # For smaller matrices, you may be able to use +#det_exact+.
  #
  # This function is guaranteed to return the same type of data in the matrix
  # upon which it is called.
  #
  # Integer matrices are converted to floating point matrices for the purposes of
  # performing the calculation, as xGETRF can't work on integer matrices.
  #
  # * *Returns* :
  #   - The determinant of the matrix. It's the same type as the matrix's dtype.
  # * *Raises* :
  #   - +NotImplementedError+ -> Must be used in 2D matrices.
  #
  def det
    raise(NotImplementedError, "determinant can be calculated only for 2D matrices") unless self.dim == 2

    # Cast to a dtype for which getrf is implemented
    new_dtype = [:byte,:int8,:int16,:int32,:int64].include?(self.dtype) ? :float64 : self.dtype
    copy = self.cast(:dense, new_dtype)

    # Need to know the number of permutations. We'll add up the diagonals of
    # the factorized matrix.
    pivot = copy.getrf!

    num_perm = 0 #number of permutations
    pivot.each_with_index do |swap, i|
      num_perm += 1 if swap != i
    end
    prod = num_perm % 2 == 1 ? -1 : 1 # odd permutations => negative
    [shape[0],shape[1]].min.times do |i|
      prod *= copy[i,i]
    end

    # Convert back to an integer if necessary
    new_dtype != self.dtype ? prod.round : prod #prevent rounding errors
  end

  #
  # call-seq:
  #     complex_conjugate -> NMatrix
  #     complex_conjugate(new_stype) -> NMatrix
  #
  # Get the complex conjugate of this matrix. See also complex_conjugate! for
  # an in-place operation (provided the dtype is already +:complex64+ or
  # +:complex128+).
  #
  # Doesn't work on list matrices, but you can optionally pass in the stype you
  # want to cast to if you're dealing with a list matrix.
  #
  # * *Arguments* :
  #   - +new_stype+ -> stype for the new matrix.
  # * *Returns* :
  #   - If the original NMatrix isn't complex, the result is a +:complex128+ NMatrix. Otherwise, it's the original dtype.
  #
  def complex_conjugate(new_stype = self.stype)
    self.cast(new_stype, NMatrix::upcast(dtype, :complex64)).complex_conjugate!
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
  def cov(opts={})
    raise TypeError, "Only works for non-integer dtypes" if integer_dtype?
     opts = {
      for_sample_data: true
    }.merge(opts)
    
    denominator      = opts[:for_sample_data] ? rows - 1 : rows
    ones             = NMatrix.ones [rows,1] 
    deviation_scores = self - ones.dot(ones.transpose).dot(self) / rows
    deviation_scores.transpose.dot(deviation_scores) / denominator
  end

  # Calculate the correlation matrix.
  def corr
    raise NotImplementedError, "Does not work for complex dtypes" if complex_dtype?
    standard_deviation = std
    cov / (standard_deviation.transpose.dot(standard_deviation))
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
      shape[0] != shape[1] or shape.size != 2
    raise TypeError, "Only works with integer powers" unless n.is_a?(Integer)

    sequence = (integer_dtype? ? self.cast(dtype: :int64) : self).clone
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
      n = n / 2
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
    unless self.dimensions==2 and mat.dimensions==2
      raise ShapeError, "Implemented for 2D NMatrix objects only."
    end

    # compute the shape [n,m] of the product matrix
    n, m = self.shape[0]*mat.shape[0], self.shape[1]*mat.shape[1]
    # compute the entries of the product matrix
    kron_prod_array = []
    if self.yale?
      # +:yale+ requires to get the row by copy in order to apply +#transpose+ to it
      self.each_row(getby=:copy) do |selfr|
        mat.each_row do |matr|
          kron_prod_array += (selfr.transpose.dot matr).to_flat_a
        end
      end
    else
      self.each_row do |selfr|
        mat.each_row do |matr|
          kron_prod_array += (selfr.transpose.dot matr).to_flat_a
        end
      end
    end

    NMatrix.new([n,m], kron_prod_array) 
  end

  #
  # call-seq:
  #     conjugate_transpose -> NMatrix
  #
  # Calculate the conjugate transpose of a matrix. If your dtype is already
  # complex, this should only require one copy (for the transpose).
  #
  # * *Returns* :
  #   - The conjugate transpose of the matrix as a copy.
  #
  def conjugate_transpose
    self.transpose.complex_conjugate!
  end

  def get_rank(workspace_size=1)
    sigmas = self.gesvd[1].to_a.flatten
    eps = NMatrix::FLOAT64_EPSILON

    # epsilon depends on the width of the number
    if (self.dtype == :float32 || self.dtype == :complex64) 
      eps = NMatrix::FLOAT32_EPSILON
    end

    tol = self.shape.max * sigmas.max * eps

    sigmas.map { |x| x > tol ? 1 : 0 }.reduce(:+)
  end

  def full_rank?
    # if it's a square matrix, save some computation time
    # and just check if determinant is > 0
    return self.det != 0 if self.shape[0] == self.shape[1]
    self.get_rank == self.shape.min
  end

  def rank_deficient?
    !full_rank?
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
    raise(ShapeError, "Expected square matrix") unless self.shape[0] == self.shape[1] && self.dim == 2

    (0...self.shape[0]).inject(0) do |total,i|
      total + self[i,i]
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
  def mean(dimen=0)
    reduce_dtype = nil
    if integer_dtype? then
      reduce_dtype = :float64
    end
    inject_rank(dimen, 0.0, reduce_dtype) do |mean, sub_mat|
      mean + sub_mat
    end / shape[dimen]
  end

  ##
  # call-seq:
  #   sum() -> NMatrix
  #   sum(dimen) -> NMatrix
  #
  # Calculates the sum along the specified dimension.
  #
  # @see #inject_rank
  def sum(dimen=0)
    inject_rank(dimen, 0.0) do |sum, sub_mat|
      sum + sub_mat
    end
  end


  ##
  # call-seq:
  #   min() -> NMatrix
  #   min(dimen) -> NMatrix
  #
  # Calculates the minimum along the specified dimension.
  #
  # @see #inject_rank
  #
  def min(dimen=0)
    inject_rank(dimen) do |min, sub_mat|
      if min.is_a? NMatrix then
        min * (min <= sub_mat).cast(self.stype, self.dtype) + ((min)*0.0 + (min > sub_mat).cast(self.stype, self.dtype)) * sub_mat
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
  def max(dimen=0)
    inject_rank(dimen) do |max, sub_mat|
      if max.is_a? NMatrix then
        max * (max >= sub_mat).cast(self.stype, self.dtype) + ((max)*0.0 + (max < sub_mat).cast(self.stype, self.dtype)) * sub_mat
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
  def variance(dimen=0)
    reduce_dtype = nil
    if integer_dtype? then
      reduce_dtype = :float64
    end
    m = mean(dimen)
    inject_rank(dimen, 0.0, reduce_dtype) do |var, sub_mat|
      var + (m - sub_mat)*(m - sub_mat)/(shape[dimen]-1)
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
  def std(dimen=0)
    variance(dimen).sqrt
  end


  #
  # call-seq:
  #     abs_dtype -> Symbol
  #
  # Returns the dtype of the result of a call to #abs. In most cases, this is the same as dtype; it should only differ
  # for :complex64 (where it's :float32) and :complex128 (:float64).
  def abs_dtype
    if self.dtype == :complex64
      :float32
    elsif self.dtype == :complex128
      :float64
    else
      self.dtype
    end
  end


  #
  # call-seq:
  #     abs -> NMatrix
  #
  # Maps all values in a matrix to their absolute values.
  def abs
    if stype == :dense
      self.__dense_map__ { |v| v.abs }
    elsif stype == :list
      # FIXME: Need __list_map_stored__, but this will do for now.
      self.__list_map_merged_stored__(nil, nil) { |v,dummy| v.abs }
    else
      self.__yale_map_stored__ { |v| v.abs }
    end.cast(self.stype, abs_dtype)
  end


  #
  # call-seq:
  #     absolute_sum -> Numeric
  #
  # == Arguments
  #   - +incx+ -> the skip size (defaults to 1, no skip)
  #   - +n+ -> the number of elements to include
  #
  # Return the sum of the contents of the vector. This is the BLAS asum routine.
  def asum incx=1, n=nil
    if self.shape == [1]
      return self[0].abs unless self.complex_dtype?
      return self[0].real.abs + self[0].imag.abs
    end
    return method_missing(:asum, incx, n) unless vector?
    NMatrix::BLAS::asum(self, incx, self.size / incx)
  end
  alias :absolute_sum :asum

  #
  # call-seq:
  #     norm2 -> Numeric
  #
  # == Arguments
  #   - +incx+ -> the skip size (defaults to 1, no skip)
  #   - +n+ -> the number of elements to include
  #
  # Return the 2-norm of the vector. This is the BLAS nrm2 routine.
  def nrm2 incx=1, n=nil
    return method_missing(:nrm2, incx, n) unless vector?
    NMatrix::BLAS::nrm2(self, incx, self.size / incx)
  end
  alias :norm2 :nrm2


  alias :permute_columns  :laswp
  alias :permute_columns! :laswp!

protected
  # Define the element-wise operations for lists. Note that the __list_map_merged_stored__ iterator returns a Ruby Object
  # matrix, which we then cast back to the appropriate type. If you don't want that, you can redefine these functions in
  # your own code.
  {add: :+, sub: :-, mul: :*, div: :/, pow: :**, mod: :%}.each_pair do |ewop, op|
    define_method("__list_elementwise_#{ewop}__") do |rhs|
      self.__list_map_merged_stored__(rhs, nil) { |l,r| l.send(op,r) }.cast(stype, NMatrix.upcast(dtype, rhs.dtype))
    end
    define_method("__dense_elementwise_#{ewop}__") do |rhs|
      self.__dense_map_pair__(rhs) { |l,r| l.send(op,r) }.cast(stype, NMatrix.upcast(dtype, rhs.dtype))
    end
    define_method("__yale_elementwise_#{ewop}__") do |rhs|
      self.__yale_map_merged_stored__(rhs, nil) { |l,r| l.send(op,r) }.cast(stype, NMatrix.upcast(dtype, rhs.dtype))
    end
    define_method("__list_scalar_#{ewop}__") do |rhs|
      self.__list_map_merged_stored__(rhs, nil) { |l,r| l.send(op,r) }.cast(stype, NMatrix.upcast(dtype, NMatrix.min_dtype(rhs)))
    end
    define_method("__yale_scalar_#{ewop}__") do |rhs|
      self.__yale_map_stored__ { |l| l.send(op,rhs) }.cast(stype, NMatrix.upcast(dtype, NMatrix.min_dtype(rhs)))
    end
    define_method("__dense_scalar_#{ewop}__") do |rhs|
      self.__dense_map__ { |l| l.send(op,rhs) }.cast(stype, NMatrix.upcast(dtype, NMatrix.min_dtype(rhs)))
    end
  end

  # These don't actually take an argument -- they're called reverse-polish style on the matrix.
  # This group always gets casted to float64.
  [:log, :log2, :log10, :sqrt, :sin, :cos, :tan, :acos, :asin, :atan, :cosh, :sinh, :tanh, :acosh,
   :asinh, :atanh, :exp, :erf, :erfc, :gamma, :cbrt, :round].each do |ewop|
    define_method("__list_unary_#{ewop}__") do
      self.__list_map_stored__(nil) { |l| Math.send(ewop, l) }.cast(stype, NMatrix.upcast(dtype, :float64))
    end
    define_method("__yale_unary_#{ewop}__") do
      self.__yale_map_stored__ { |l| Math.send(ewop, l) }.cast(stype, NMatrix.upcast(dtype, :float64))
    end
    define_method("__dense_unary_#{ewop}__") do
      self.__dense_map__ { |l| Math.send(ewop, l) }.cast(stype, NMatrix.upcast(dtype, :float64))
    end
  end

  #:stopdoc:
  # log takes an optional single argument, the base. Default to natural log.
  def __list_unary_log__(base)
    self.__list_map_stored__(nil) { |l| Math.log(l, base) }.cast(stype, NMatrix.upcast(dtype, :float64))
  end

  def __yale_unary_log__(base)
    self.__yale_map_stored__ { |l| Math.log(l, base) }.cast(stype, NMatrix.upcast(dtype, :float64))
  end

  def __dense_unary_log__(base)
    self.__dense_map__ { |l| Math.log(l, base) }.cast(stype, NMatrix.upcast(dtype, :float64))
  end

  # These are for negating matrix contents using -@
  def __list_unary_negate__
    self.__list_map_stored__(nil) { |l| -l }.cast(stype, dtype)
  end

  def __yale_unary_negate__
    self.__yale_map_stored__ { |l| -l }.cast(stype, dtype)
  end

  def __dense_unary_negate__
    self.__dense_map__ { |l| -l }.cast(stype, dtype)
  end
  #:startdoc:

  # These are for rounding each value of a matrix. Takes an optional argument
  def __list_unary_round__(precision)
    if self.complex_dtype?
      self.__list_map_stored__(nil) { |l| Complex(l.real.round(precision), l.imag.round(precision)) }
                                    .cast(stype, dtype)
    else
      self.__list_map_stored__(nil) { |l| l.round(precision) }.cast(stype, dtype)
    end
  end

  def __yale_unary_round__(precision)
    if self.complex_dtype?
      self.__yale_map_stored__ { |l| Complex(l.real.round(precision), l.imag.round(precision)) }
                                    .cast(stype, dtype)
    else
      self.__yale_map_stored__ { |l| l.round(precision) }.cast(stype, dtype)
    end
  end

  def __dense_unary_round__(precision)
    if self.complex_dtype?
      self.__dense_map__ { |l| Complex(l.real.round(precision), l.imag.round(precision)) }
                                    .cast(stype, dtype)
    else
      self.__dense_map__ { |l| l.round(precision) }.cast(stype, dtype)
    end
  end

  # These are for calculating the floor or ceil of matrix
  def dtype_for_floor_or_ceil
    if self.integer_dtype? or [:complex64, :complex128, :object].include?(self.dtype)
      return_dtype = dtype
    elsif [:float32, :float64].include?(self.dtype)
      return_dtype = :int64
    end

    return_dtype
  end

  [:floor, :ceil].each do |meth|
    define_method("__list_unary_#{meth}__") do
      return_dtype = dtype_for_floor_or_ceil

      if [:complex64, :complex128].include?(self.dtype)
        self.__list_map_stored__(nil) { |l| Complex(l.real.send(meth), l.imag.send(meth)) }.cast(stype, return_dtype)
      else
        self.__list_map_stored__(nil) { |l| l.send(meth) }.cast(stype, return_dtype)
      end
    end

    define_method("__yale_unary_#{meth}__") do
      return_dtype = dtype_for_floor_or_ceil

      if [:complex64, :complex128].include?(self.dtype)
        self.__yale_map_stored__ { |l| Complex(l.real.send(meth), l.imag.send(meth)) }.cast(stype, return_dtype)
      else
        self.__yale_map_stored__ { |l| l.send(meth) }.cast(stype, return_dtype)
      end
    end

    define_method("__dense_unary_#{meth}__") do
      return_dtype = dtype_for_floor_or_ceil

      if [:complex64, :complex128].include?(self.dtype)
        self.__dense_map__ { |l| Complex(l.real.send(meth), l.imag.send(meth)) }.cast(stype, return_dtype)
      else
        self.__dense_map__ { |l| l.send(meth) }.cast(stype, return_dtype)
      end
    end
  end

  # These take two arguments. One might be a matrix, and one might be a scalar.
  # See also monkeys.rb, which contains Math module patches to let the first
  # arg be a scalar
  [:atan2, :ldexp, :hypot].each do |ewop|
    define_method("__list_elementwise_#{ewop}__") do |rhs,order|
      if order then
        self.__list_map_merged_stored__(rhs, nil) { |r,l| Math.send(ewop,l,r) }
      else
        self.__list_map_merged_stored__(rhs, nil) { |l,r| Math.send(ewop,l,r) }
      end.cast(stype, NMatrix.upcast(dtype, :float64))
    end

    define_method("__dense_elementwise_#{ewop}__") do |rhs, order|
      if order then
        self.__dense_map_pair__(rhs) { |r,l| Math.send(ewop,l,r) }
      else
        self.__dense_map_pair__(rhs) { |l,r| Math.send(ewop,l,r) }
      end.cast(stype, NMatrix.upcast(dtype, :float64))
    end

    define_method("__yale_elementwise_#{ewop}__") do |rhs, order|
      if order then
        self.__yale_map_merged_stored__(rhs, nil) { |r,l| Math.send(ewop,l,r) }
      else
        self.__yale_map_merged_stored__(rhs, nil) { |l,r| Math.send(ewop,l,r) }
      end.cast(stype, NMatrix.upcast(dtype, :float64))
    end

    define_method("__list_scalar_#{ewop}__") do |rhs,order|
      if order then
        self.__list_map_stored__(nil) { |l| Math.send(ewop, rhs, l) }
      else
        self.__list_map_stored__(nil) { |l| Math.send(ewop, l, rhs) }
      end.cast(stype, NMatrix.upcast(dtype, :float64))
    end

    define_method("__yale_scalar_#{ewop}__") do |rhs,order|
      if order then
        self.__yale_map_stored__ { |l| Math.send(ewop, rhs, l) }
      else
        self.__yale_map_stored__ { |l| Math.send(ewop, l, rhs) }
      end.cast(stype, NMatrix.upcast(dtype, :float64))
    end

    define_method("__dense_scalar_#{ewop}__") do |rhs,order|
      if order
        self.__dense_map__ { |l| Math.send(ewop, rhs, l) }
      else
        self.__dense_map__ { |l| Math.send(ewop, l, rhs) }
      end.cast(stype, NMatrix.upcast(dtype, :float64))
    end
  end

  # Equality operators do not involve a cast. We want to get back matrices of TrueClass and FalseClass.
  {eqeq: :==, neq: :!=, lt: :<, gt: :>, leq: :<=, geq: :>=}.each_pair do |ewop, op|
    define_method("__list_elementwise_#{ewop}__") do |rhs|
      self.__list_map_merged_stored__(rhs, nil) { |l,r| l.send(op,r) }
    end
    define_method("__dense_elementwise_#{ewop}__") do |rhs|
      self.__dense_map_pair__(rhs) { |l,r| l.send(op,r) }
    end
    define_method("__yale_elementwise_#{ewop}__") do |rhs|
      self.__yale_map_merged_stored__(rhs, nil) { |l,r| l.send(op,r) }
    end

    define_method("__list_scalar_#{ewop}__") do |rhs|
      self.__list_map_merged_stored__(rhs, nil) { |l,r| l.send(op,r) }
    end
    define_method("__yale_scalar_#{ewop}__") do |rhs|
      self.__yale_map_stored__ { |l| l.send(op,rhs) }
    end
    define_method("__dense_scalar_#{ewop}__") do |rhs|
      self.__dense_map__ { |l| l.send(op,rhs) }
    end
  end
end
