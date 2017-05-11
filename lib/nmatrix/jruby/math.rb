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
      def permutation_matrix_from(pivot_array)
        perm_arry = permutation_array_for(pivot_array)
        n         = NMatrix.zeros(perm_arry.size, dtype: :byte)

        perm_arry.each_with_index { |e, i| n[e,i] = 1 }

        n
      end

      def permutation_array_for(pivot_array)
        perm_arry = Array.new(pivot_array.size) { |i| i }
        perm_arry.each_index do |i|
          #the pivot indices returned by LAPACK getrf are indexed starting
          #from 1, so we need to subtract 1 here
          perm_arry[i], perm_arry[pivot_array[i]-1] = perm_arry[pivot_array[i]-1], perm_arry[i]
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
    raise(StorageTypeError, "invert only works on dense matrices currently") unless self.dense?
    raise(ShapeError, "Cannot invert non-square matrix") unless self.dim == 2 && self.shape[0] == self.shape[1]
    raise(DataTypeError, "Cannot invert an integer matrix in-place") if self.integer_dtype?

    #No internal implementation of getri, so use this other function
    __inverse__!
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
    #write this in terms of invert! so plugins will only have to overwrite
    #invert! and not invert
    __inverse__
  end
  alias :inverse :invert

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
    raise(StorageTypeError, "adjugate only works on dense matrices currently") unless self.dense?
    raise(ShapeError, "Cannot calculate adjugate of a non-square matrix") unless self.dim == 2 && self.shape[0] == self.shape[1]
    raise(DataTypeError, "Cannot calculate adjugate of an integer matrix in-place") if self.integer_dtype?
    d = self.det
    self.invert!
    self.map! { |e| e * d }
    self
  end
  alias :adjoint! :adjugate!

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
    raise(StorageTypeError, "adjugate only works on dense matrices currently") unless self.dense?
    raise(ShapeError, "Cannot calculate adjugate of a non-square matrix") unless self.dim == 2 && self.shape[0] == self.shape[1]
    d = self.det
    mat = self.invert
    mat.map! { |e| e * d }
    mat
  end
  alias :adjoint :adjugate

  #
  # call-seq:
  #     getrf! -> Array
  #
  # LU factorization of a general M-by-N matrix +A+ using partial pivoting with
  # row interchanges. The LU factorization is A = PLU, where P is a row permutation
  # matrix, L is a lower triangular matrix with unit diagonals, and U is an upper
  # triangular matrix (note that this convention is different from the
  # clapack_getrf behavior, but matches the standard LAPACK getrf).
  # +A+ is overwritten with the elements of L and U (the unit
  # diagonal elements of L are not saved). P is not returned directly and must be
  # constructed from the pivot array ipiv. The row indices in ipiv are indexed
  # starting from 1.
  # Only works for dense matrices.
  #
  # * *Returns* :
  #   - The IPIV vector. The L and U matrices are stored in A.
  # * *Raises* :
  #   - +StorageTypeError+ -> ATLAS functions only work on dense matrices.
  #
  def getrf!
    ipiv = LUDecomposition.new(self.twoDMat).getPivot.to_a
    return ipiv
  end

  #
  # call-seq:
  #     geqrf! -> shape.min x 1 NMatrix
  #
  # QR factorization of a general M-by-N matrix +A+.
  #
  # The QR factorization is A = QR, where Q is orthogonal and R is Upper Triangular
  # +A+ is overwritten with the elements of R and Q with Q being represented by the
  # elements below A's diagonal and an array of scalar factors in the output NMatrix.
  #
  # The matrix Q is represented as a product of elementary reflectors
  #     Q = H(1) H(2) . . . H(k), where k = min(m,n).
  #
  # Each H(i) has the form
  #
  #     H(i) = I - tau * v * v'
  #
  # http://www.netlib.org/lapack/explore-html/d3/d69/dgeqrf_8f.html
  #
  # Only works for dense matrices.
  #
  # * *Returns* :
  #   - Vector TAU. Q and R are stored in A. Q is represented by TAU and A
  # * *Raises* :
  #   - +StorageTypeError+ -> LAPACK functions only work on dense matrices.
  #
  def geqrf!
    # The real implementation is in lib/nmatrix/lapacke.rb
    raise(NotImplementedError, "geqrf! requires the nmatrix-lapacke gem")
  end

  #
  # call-seq:
  #     ormqr(tau) -> NMatrix
  #     ormqr(tau, side, transpose, c) -> NMatrix
  #
  # Returns the product Q * c or c * Q after a call to geqrf! used in QR factorization.
  # +c+ is overwritten with the elements of the result NMatrix if supplied. Q is the orthogonal matrix
  # represented by tau and the calling NMatrix
  #
  # Only works on float types, use unmqr for complex types.
  #
  # == Arguments
  #
  # * +tau+ - vector containing scalar factors of elementary reflectors
  # * +side+ - direction of multiplication [:left, :right]
  # * +transpose+ - apply Q with or without transpose [false, :transpose]
  # * +c+ - NMatrix multplication argument that is overwritten, no argument assumes c = identity
  #
  # * *Returns* :
  #
  #   - Q * c or c * Q Where Q may be transposed before multiplication.
  #
  #
  # * *Raises* :
  #   - +StorageTypeError+ -> LAPACK functions only work on dense matrices.
  #   - +TypeError+ -> Works only on floating point matrices, use unmqr for complex types
  #   - +TypeError+ -> c must have the same dtype as the calling NMatrix
  #
  def ormqr(tau, side=:left, transpose=false, c=nil)
    # The real implementation is in lib/nmatrix/lapacke.rb
    raise(NotImplementedError, "ormqr requires the nmatrix-lapacke gem")

  end

  #
  # call-seq:
  #     unmqr(tau) -> NMatrix
  #     unmqr(tau, side, transpose, c) -> NMatrix
  #
  # Returns the product Q * c or c * Q after a call to geqrf! used in QR factorization.
  # +c+ is overwritten with the elements of the result NMatrix if it is supplied. Q is the orthogonal matrix
  # represented by tau and the calling NMatrix
  #
  # Only works on complex types, use ormqr for float types.
  #
  # == Arguments
  #
  # * +tau+ - vector containing scalar factors of elementary reflectors
  # * +side+ - direction of multiplication [:left, :right]
  # * +transpose+ - apply Q as Q or its complex conjugate [false, :complex_conjugate]
  # * +c+ - NMatrix multplication argument that is overwritten, no argument assumes c = identity
  #
  # * *Returns* :
  #
  #   - Q * c or c * Q Where Q may be transformed to its complex conjugate before multiplication.
  #
  #
  # * *Raises* :
  #   - +StorageTypeError+ -> LAPACK functions only work on dense matrices.
  #   - +TypeError+ -> Works only on floating point matrices, use unmqr for complex types
  #   - +TypeError+ -> c must have the same dtype as the calling NMatrix
  #
  def unmqr(tau, side=:left, transpose=false, c=nil)
    # The real implementation is in lib/nmatrix/lapacke.rb
    raise(NotImplementedError, "unmqr requires the nmatrix-lapacke gem")
  end

  #
  # call-seq:
  #     potrf!(upper_or_lower) -> NMatrix
  #
  # Cholesky factorization of a symmetric positive-definite matrix -- or, if complex,
  # a Hermitian positive-definite matrix +A+.
  # The result will be written in either the upper or lower triangular portion of the
  # matrix, depending on whether the argument is +:upper+ or +:lower+.
  # Also the function only reads in the upper or lower part of the matrix,
  # so it doesn't actually have to be symmetric/Hermitian.
  # However, if the matrix (i.e. the symmetric matrix implied by the lower/upper
  # half) is not positive-definite, the function will return nonsense.
  #
  # This functions requires either the nmatrix-atlas or nmatrix-lapacke gem
  # installed.
  #
  # * *Returns* :
  #   the triangular portion specified by the parameter
  # * *Raises* :
  #   - +StorageTypeError+ -> ATLAS functions only work on dense matrices.
  #   - +ShapeError+ -> Must be square.
  #   - +NotImplementedError+ -> If called without nmatrix-atlas or nmatrix-lapacke gem
  #
  def potrf!(which)
    # The real implementation is in the plugin files.
    cholesky = CholeskyDecomposition.new(self.twoDMat)
    if which == :upper
      u = create_dummy_nmatrix
      twoDMat = cholesky.getLT
      u.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(twoDMat.getData, @shape[0], @shape[1]))
      return u
    else
      l = create_dummy_nmatrix
      twoDMat = cholesky.getL
      l.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(twoDMat.getData, @shape[0], @shape[1]))
      return l
    end
  end

  def potrf_upper!
    potrf! :upper
  end

  def potrf_lower!
    potrf! :lower
  end


  #
  # call-seq:
  #     factorize_cholesky -> [upper NMatrix, lower NMatrix]
  #
  # Calculates the Cholesky factorization of a matrix and returns the
  # upper and lower matrices such that A=LU and L=U*, where * is
  # either the transpose or conjugate transpose.
  #
  # Unlike potrf!, this makes method requires that the original is matrix is
  # symmetric or Hermitian. However, it is still your responsibility to make
  # sure it is positive-definite.
  def factorize_cholesky
    # raise "Matrix must be symmetric/Hermitian for Cholesky factorization" unless self.hermitian?
    cholesky = CholeskyDecomposition.new(self.twoDMat)
    l = create_dummy_nmatrix
    twoDMat = cholesky.getL
    l.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(twoDMat.getData, @shape[0], @shape[1]))
    u = create_dummy_nmatrix
    twoDMat = cholesky.getLT
    u.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(twoDMat.getData, @shape[0], @shape[1]))
    return [u,l]
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
  def factorize_lu with_permutation_matrix=nil
    raise(NotImplementedError, "only implemented for dense storage") unless self.stype == :dense
    raise(NotImplementedError, "matrix is not 2-dimensional") unless self.dimensions == 2
    t = self.clone
    pivot = create_dummy_nmatrix
    twoDMat = LUDecomposition.new(self.twoDMat).getP
    pivot.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(twoDMat.getData, @shape[0], @shape[1]))
    return [t,pivot]
  end

  #
  # call-seq:
  #     factorize_qr -> [Q,R]
  #
  # QR factorization of a matrix without column pivoting.
  # Q is orthogonal and R is upper triangular if input is square or upper trapezoidal if
  # input is rectangular.
  #
  # Only works for dense matrices.
  #
  # * *Returns* :
  #   - Array containing Q and R matrices
  #
  # * *Raises* :
  #   - +StorageTypeError+ -> only implemented for desnse storage.
  #   - +ShapeError+ -> Input must be a 2-dimensional matrix to have a QR decomposition.
  #
  def factorize_qr

    raise(NotImplementedError, "only implemented for dense storage") unless self.stype == :dense
    raise(ShapeError, "Input must be a 2-dimensional matrix to have a QR decomposition") unless self.dim == 2
    qrdecomp = QRDecomposition.new(self.twoDMat)

    qmat = create_dummy_nmatrix
    qtwoDMat = qrdecomp.getQ
    qmat.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(qtwoDMat.getData, @shape[0], @shape[1]))

    rmat = create_dummy_nmatrix
    rtwoDMat = qrdecomp.getR
    rmat.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(rtwoDMat.getData, @shape[0], @shape[1]))
    return [qmat,rmat]

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

  # Solve the matrix equation AX = B, where A is +self+, B is the first
  # argument, and X is returned. A must be a nxn square matrix, while B must be
  # nxm. Only works with dense matrices and non-integer, non-object data types.
  #
  # == Arguments
  #
  # * +b+ - the right hand side
  #
  # == Options
  #
  # * +form+ - Signifies the form of the matrix A in the linear system AX=B.
  #   If not set then it defaults to +:general+, which uses an LU solver.
  #   Other possible values are +:lower_tri+, +:upper_tri+ and +:pos_def+ (alternatively,
  #   non-abbreviated symbols +:lower_triangular+, +:upper_triangular+,
  #   and +:positive_definite+ can be used.
  #   If +:lower_tri+ or +:upper_tri+ is set, then a specialized linear solver for linear
  #   systems AX=B with a lower or upper triangular matrix A is used. If +:pos_def+ is chosen,
  #   then the linear system is solved via the Cholesky factorization.
  #   Note that when +:lower_tri+ or +:upper_tri+ is used, then the algorithm just assumes that
  #   all entries in the lower/upper triangle of the matrix are zeros without checking (which
  #   can be useful in certain applications).
  #
  #
  # == Usage
  #
  #   a = NMatrix.new [2,2], [3,1,1,2], dtype: dtype
  #   b = NMatrix.new [2,1], [9,8], dtype: dtype
  #   a.solve(b)
  #
  #   # solve an upper triangular linear system more efficiently:
  #   require 'benchmark'
  #   require 'nmatrix/lapacke'
  #   rand_mat = NMatrix.random([10000, 10000], dtype: :float64)
  #   a = rand_mat.triu
  #   b = NMatrix.random([10000, 10], dtype: :float64)
  #   Benchmark.bm(10) do |bm|
  #     bm.report('general') { a.solve(b) }
  #     bm.report('upper_tri') { a.solve(b, form: :upper_tri) }
  #   end
  #   #                   user     system      total        real
  #   #  general     73.170000   0.670000  73.840000 ( 73.810086)
  #   #  upper_tri    0.180000   0.000000   0.180000 (  0.182491)
  #
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

    nmatrix = create_dummy_nmatrix
    case opts[:form]
    when :general, :upper_tri, :upper_triangular, :lower_tri, :lower_triangular
      #LU solver
      solver = LUDecomposition.new(self.twoDMat).getSolver
      nmatrix.s = solver.solve(b.s)
      return nmatrix
    when :pos_def, :positive_definite
      solver = CholeskyDecomposition.new(self.twoDMat).getSolver
      nmatrix.s = solver.solve(b.s)
      return nmatrix
    else
      raise(ArgumentError, "#{opts[:form]} is not a valid form option")
    end

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
  # using clapack_getrf, and then by taking the product of the diagonal elements. There is a
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
  #   - +ShapeError+ -> Must be used on square matrices.
  #
  def det
    raise(ShapeError, "determinant can be calculated only for square matrices") unless self.dim == 2 && self.shape[0] == self.shape[1]
    self.det_exact2
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
    self.twoDMat.getFrobeniusNorm()
  end
  alias :norm2 :nrm2

  #
  # call-seq:
  #     scale! -> NMatrix
  #
  # == Arguments
  #   - +alpha+ -> Scalar value used in the operation.
  #   - +inc+ -> Increment used in the scaling function. Should generally be 1.
  #   - +n+ -> Number of elements of +vector+.
  #
  # This is a destructive method, modifying the source NMatrix.  See also #scale.
  # Return the scaling result of the matrix. BLAS scal will be invoked if provided.

  def scale!(alpha, incx=1, n=nil)
    #FIXME
    # raise(DataTypeError, "Incompatible data type for the scaling factor") unless
    #     NMatrix::upcast(self.dtype, NMatrix::min_dtype(alpha)) == self.dtype
    raise(DataTypeError, "Incompatible data type for the scaling factor") if
        self.dtype == :int8
    @s.mapMultiplyToSelf(alpha)
    return self
  end

  #
  # call-seq:
  #     scale -> NMatrix
  #
  # == Arguments
  #   - +alpha+ -> Scalar value used in the operation.
  #   - +inc+ -> Increment used in the scaling function. Should generally be 1.
  #   - +n+ -> Number of elements of +vector+.
  #
  # Return the scaling result of the matrix. BLAS scal will be invoked if provided.

  def scale(alpha, incx=1, n=nil)
    # FIXME
    # raise(DataTypeError, "Incompatible data type for the scaling factor") unless
    #     NMatrix::upcast(self.dtype, NMatrix::min_dtype(alpha)) == self.dtype
    raise(DataTypeError, "Incompatible data type for the scaling factor") if
        self.dtype == :byte || self.dtype == :int8 || self.dtype == :int16 ||
        self.dtype == :int32 || self.dtype == :int64
    nmatrix = NMatrix.new :copy
    nmatrix.shape = @shape.clone
    nmatrix.s = ArrayRealVector.new(@s.toArray.clone).mapMultiplyToSelf(alpha)
    return nmatrix
  end

  alias :permute_columns  :laswp
  alias :permute_columns! :laswp!

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
      self.dim == 2 && self.shape[0] == self.shape[1]
    cond = 0
    while cond != self.cols
      if self[0..cond, 0..cond].det <= 0
        return false
      end
      cond += 1
    end
    return true
  end

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

