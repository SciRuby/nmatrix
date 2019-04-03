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
    ipiv = LUDecomposition.new(twoDMat).getPivot.to_a
    ipiv
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
  def ormqr(tau, side = :left, transpose = false, c = nil)
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
  def unmqr(tau, side = :left, transpose = false, c = nil)
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
    cholesky = CholeskyDecomposition.new(twoDMat)
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
    cholesky = CholeskyDecomposition.new(twoDMat)
    l = create_dummy_nmatrix
    twoDMat = cholesky.getL
    l.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(twoDMat.getData, @shape[0], @shape[1]))
    u = create_dummy_nmatrix
    twoDMat = cholesky.getLT
    u.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(twoDMat.getData, @shape[0], @shape[1]))
    [u, l]
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
  def factorize_lu with_permutation_matrix = nil
    raise(NotImplementedError, "only implemented for dense storage") unless stype == :dense
    raise(NotImplementedError, "matrix is not 2-dimensional") unless dimensions == 2
    t = clone
    pivot = create_dummy_nmatrix
    twoDMat = LUDecomposition.new(self.twoDMat).getP
    pivot.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(twoDMat.getData, @shape[0], @shape[1]))
    [t, pivot]
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
    raise(NotImplementedError, "only implemented for dense storage") unless stype == :dense
    raise(ShapeError, "Input must be a 2-dimensional matrix to have a QR decomposition") unless dim == 2
    qrdecomp = QRDecomposition.new(twoDMat)

    qmat = create_dummy_nmatrix
    qtwoDMat = qrdecomp.getQ
    qmat.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(qtwoDMat.getData, @shape[0], @shape[1]))

    rmat = create_dummy_nmatrix
    rtwoDMat = qrdecomp.getR
    rmat.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(rtwoDMat.getData, @shape[0], @shape[1]))
    [qmat, rmat]
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
    raise(ShapeError, "Must be called on square matrix") unless dim == 2 && shape[0] == shape[1]
    raise(ShapeError, "number of rows of b must equal number of cols of self") if
      shape[1] != b.shape[0]
    raise(ArgumentError, "only works with dense matrices") if stype != :dense
    raise(ArgumentError, "only works for non-integer, non-object dtypes") if
      integer_dtype? || object_dtype? || b.integer_dtype? || b.object_dtype?

    opts = {form: :general}.merge(opts)
    x    = b.clone
    n    = shape[0]
    nrhs = b.shape[1]

    nmatrix = create_dummy_nmatrix
    case opts[:form]
    when :general, :upper_tri, :upper_triangular, :lower_tri, :lower_triangular
      # LU solver
      solver = LUDecomposition.new(twoDMat).getSolver
      nmatrix.s = solver.solve(b.s)
      return nmatrix
    when :pos_def, :positive_definite
      solver = CholeskyDecomposition.new(twoDMat).getSolver
      nmatrix.s = solver.solve(b.s)
      return nmatrix
    else
      raise(ArgumentError, "#{opts[:form]} is not a valid form option")
    end
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
    raise(ShapeError, "determinant can be calculated only for square matrices") unless dim == 2 && shape[0] == shape[1]
    det_exact2
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
  def complex_conjugate(new_stype = stype)
    cast(new_stype, NMatrix.upcast(dtype, :complex64)).complex_conjugate!
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
    transpose.complex_conjugate!
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
  def asum incx = 1, n = nil
    if shape == [1]
      return self[0].abs unless complex_dtype?
      return self[0].real.abs + self[0].imag.abs
    end
    return method_missing(:asum, incx, n) unless vector?
    NMatrix::BLAS.asum(self, incx, size / incx)
  end
  alias absolute_sum asum

  #
  # call-seq:
  #     norm2 -> Numeric
  #
  # == Arguments
  #   - +incx+ -> the skip size (defaults to 1, no skip)
  #   - +n+ -> the number of elements to include
  #
  # Return the 2-norm of the vector. This is the BLAS nrm2 routine.
  def nrm2 incx = 1, n = nil
    twoDMat.getFrobeniusNorm
  end
  alias norm2 nrm2

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

  def scale!(alpha, incx = 1, n = nil)
    # FIXME
    # raise(DataTypeError, "Incompatible data type for the scaling factor") unless
    #     NMatrix::upcast(self.dtype, NMatrix::min_dtype(alpha)) == self.dtype
    raise(DataTypeError, "Incompatible data type for the scaling factor") if
        dtype == :int8
    @s.mapMultiplyToSelf(alpha)
    self
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

  def scale(alpha, incx = 1, n = nil)
    # FIXME
    # raise(DataTypeError, "Incompatible data type for the scaling factor") unless
    #     NMatrix::upcast(self.dtype, NMatrix::min_dtype(alpha)) == self.dtype
    raise(DataTypeError, "Incompatible data type for the scaling factor") if
        dtype == :byte || dtype == :int8 || dtype == :int16 ||
          dtype == :int32 || dtype == :int64
    nmatrix = NMatrix.new :copy
    nmatrix.shape = @shape.clone
    nmatrix.s = ArrayRealVector.new(@s.toArray.clone).mapMultiplyToSelf(alpha)
    nmatrix
  end
end
