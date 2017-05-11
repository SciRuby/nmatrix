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
    #write this in terms of invert! so plugins will only have to overwrite
    #invert! and not invert
    if self.integer_dtype?
      cloned = self.cast(dtype: :float64)
      cloned.invert!
    else
      cloned = self.clone
      cloned.invert!
    end
  end
  alias :inverse :invert

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
    raise(ShapeError, "Cannot invert non-square matrix") unless self.dim == 2 && self.shape[0] == self.shape[1]
    raise(DataTypeError, "Cannot invert an integer matrix in-place") if self.integer_dtype?
    #No internal implementation of getri, so use this other function
    n = self.shape[0]
    if n>3
      raise(NotImplementedError, "Cannot find exact inverse of matrix of size greater than 3")
    else
      clond=self.clone
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
    #write this in terms of exact_inverse! so plugins will only have to overwrite
    #exact_inverse! and not exact_inverse
    if self.integer_dtype?
      cloned = self.cast(dtype: :float64)
      cloned.exact_inverse!
    else
      cloned = self.clone
      cloned.exact_inverse!
    end
  end
  alias :invert_exactly :exact_inverse



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
    if self.complex_dtype?
      u, s, vt = self.complex_conjugate.gesvd # singular value decomposition
    else
      u, s, vt = self.gesvd
    end
    rows = self.shape[0]
    cols = self.shape[1]
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
  alias :pseudo_inverse :pinv
  alias :pseudoinverse :pinv


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
    raise(StorageTypeError, "ATLAS functions only work on dense matrices") unless self.dense?

    #For row-major matrices, clapack_getrf uses a different convention than
    #described above (U has unit diagonal elements instead of L and columns
    #are interchanged rather than rows). For column-major matrices, clapack
    #uses the stanard conventions. So we just transpose the matrix before
    #and after calling clapack_getrf.
    #Unfortunately, this is not a very good way, uses a lot of memory.
    temp = self.transpose
    ipiv = NMatrix::LAPACK::clapack_getrf(:col, self.shape[0], self.shape[1], temp, self.shape[0])
    temp = temp.transpose
    self[0...self.shape[0], 0...self.shape[1]] = temp

    #for some reason, in clapack_getrf, the indices in ipiv start from 0
    #instead of 1 as in LAPACK.
    ipiv.each_index { |i| ipiv[i]+=1 }

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
    raise(NotImplementedError, "potrf! requires either the nmatrix-atlas or nmatrix-lapacke gem")
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
    raise "Matrix must be symmetric/Hermitian for Cholesky factorization" unless self.hermitian?
    l = self.clone.potrf_lower!.tril!
    u = l.conjugate_transpose
    [u,l]
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

    t     = self.clone
    pivot = t.getrf!
    return t unless with_permutation_matrix

    [t, FactorizeLUMethods.permutation_matrix_from(pivot)]
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

    rows, columns = self.shape
    r = self.clone
    tau =  r.geqrf!

    #Obtain Q
    q = self.complex_dtype? ? r.unmqr(tau) : r.ormqr(tau)

    #Obtain R
    if rows <= columns
      r.upper_triangle!
    #Need to account for upper trapezoidal structure if R is a tall rectangle (rows > columns)
    else
      r[0...columns, 0...columns].upper_triangle!
      r[columns...rows, 0...columns] = 0
    end

    [q,r]
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

    case opts[:form]
    when :general
      clone = self.clone
      ipiv = NMatrix::LAPACK.clapack_getrf(:row, n, n, clone, n)
      # When we call clapack_getrs with :row, actually only the first matrix
      # (i.e. clone) is interpreted as row-major, while the other matrix (x)
      # is interpreted as column-major. See here: http://math-atlas.sourceforge.net/faq.html#RowSolve
      # So we must transpose x before and after
      # calling it.
      x = x.transpose
      NMatrix::LAPACK.clapack_getrs(:row, :no_transpose, n, nrhs, clone, n, ipiv, x, n)
      x.transpose
    when :upper_tri, :upper_triangular
      raise(ArgumentError, "upper triangular solver does not work with complex dtypes") if
        complex_dtype? or b.complex_dtype?
      # this is the correct function call; see https://github.com/SciRuby/nmatrix/issues/374
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

  #
  # call-seq:
  #     least_squares(b) -> NMatrix
  #     least_squares(b, tolerance: 10e-10) -> NMatrix
  #
  # Provides the linear least squares approximation of an under-determined system
  # using QR factorization provided that the matrix is not rank-deficient.
  #
  # Only works for dense matrices.
  #
  # * *Arguments* :
  #   - +b+ -> The solution column vector NMatrix of A * X = b.
  #   - +tolerance:+ -> Absolute tolerance to check if a diagonal element in A = QR is near 0
  #
  # * *Returns* :
  #   - NMatrix that is a column vector with the LLS solution
  #
  # * *Raises* :
  #   - +ArgumentError+ -> least squares approximation only works for non-complex types
  #   - +ShapeError+ -> system must be under-determined ( rows > columns )
  #
  # Examples :-
  #
  #   a = NMatrix.new([3,2], [2.0, 0, -1, 1, 0, 2])
  #
  #   b = NMatrix.new([3,1], [1.0, 0, -1])
  #
  #   a.least_squares(b)
  #     =>[
  #         [ 0.33333333333333326 ]
  #         [ -0.3333333333333334 ]
  #       ]
  #
  def least_squares(b, tolerance: 10e-6)
    raise(ArgumentError, "least squares approximation only works for non-complex types") if
      self.complex_dtype?

    rows, columns = self.shape

    raise(ShapeError, "system must be under-determined ( rows > columns )") unless
      rows > columns

    #Perform economical QR factorization
    r = self.clone
    tau = r.geqrf!
    q_transpose_b = r.ormqr(tau, :left, :transpose, b)

    #Obtain R from geqrf! intermediate
    r[0...columns, 0...columns].upper_triangle!
    r[columns...rows, 0...columns] = 0

    diagonal = r.diagonal

    raise(ArgumentError, "rank deficient matrix") if diagonal.any? { |x| x == 0 }

    if diagonal.any? { |x| x.abs < tolerance }
      warn "warning: A diagonal element of R in A = QR is close to zero ;" <<
           " indicates a possible loss of precision"
    end

    # Transform the system A * X = B to R1 * X = B2 where B2 = Q1_t * B
    r1 = r[0...columns, 0...columns]
    b2 = q_transpose_b[0...columns]

    nrhs = b2.shape[1]

    #Solve the upper triangular system
    NMatrix::BLAS::cblas_trsm(:row, :left, :upper, false, :nounit, r1.shape[0], nrhs, 1.0, r1, r1.shape[0], b2, nrhs)
    b2
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
    raise(NotImplementedError, "norm can be calculated only for 2D matrices") unless self.dim == 2
    raise(NotImplementedError, "norm only implemented for dense storage") unless self.stype == :dense
    raise(ArgumentError, "norm not defined for byte dtype")if self.dtype == :byte
    case type
    when nil, 2, -2
      return self.two_matrix_norm (type == -2)
    when 1, -1
      return self.one_matrix_norm (type == -1)
    when :frobenius, :fro
      return self.fro_matrix_norm
    when :infinity, :inf, :'-inf', :'-infinity'
      return self.inf_matrix_norm  (type == :'-inf' || type == :'-infinity')
    else
      raise ArgumentError.new("argument must be a valid integer or symbol")
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
    raise(ShapeError, "determinant can be calculated only for square matrices") unless self.dim == 2 && self.shape[0] == self.shape[1]

    # Cast to a dtype for which getrf is implemented
    new_dtype = self.integer_dtype? ? :float64 : self.dtype
    copy = self.cast(:dense, new_dtype)

    # Need to know the number of permutations. We'll add up the diagonals of
    # the factorized matrix.
    pivot = copy.getrf!

    num_perm = 0 #number of permutations
    pivot.each_with_index do |swap, i|
      #pivot indexes rows starting from 1, instead of 0, so need to subtract 1 here
      num_perm += 1 if swap-1 != i
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
  #   cumsum() -> NMatrix
  #   sum(dimen) -> NMatrix
  #   cumsum(dimen) -> NMatrix
  #
  # Calculates the sum along the specified dimension.
  #
  # @see #inject_rank
  def sum(dimen=0)
    inject_rank(dimen, 0.0) do |sum, sub_mat|
      sum + sub_mat
    end
  end
  alias :cumsum :sum

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

  # Norm calculation methods
  # Frobenius norm: the Euclidean norm of the matrix, treated as if it were a vector
  def fro_matrix_norm
    #float64 has to be used in any case, since nrm2 will not yield correct result for float32
    self_cast = self.cast(:dtype => :float64)

    column_vector = self_cast.reshape([self.size, 1])

    return column_vector.nrm2
  end

  # 2-norm: the largest/smallest singular value of the matrix
  def two_matrix_norm minus = false

    self_cast = self.cast(:dtype => :float64)

    #TODO: confirm if this is the desired svd calculation
    svd = self_cast.gesvd
    return svd[1][0, 0] unless minus
    return svd[1][svd[1].rows-1, svd[1].cols-1]
  end

  # 1-norm: the maximum/minimum absolute column sum of the matrix
  def one_matrix_norm minus = false
    #TODO: change traversing method for sparse matrices
    number_of_columns = self.cols
    col_sums = []

    number_of_columns.times do |i|
      col_sums << self.col(i).inject(0) { |sum, number| sum += number.abs}
    end

    return col_sums.max unless minus
    return col_sums.min
  end

  # Infinity norm: the maximum/minimum absolute row sum of the matrix
  def inf_matrix_norm minus = false
    number_of_rows = self.rows
    row_sums = []

    number_of_rows.times do |i|
      row_sums << self.row(i).inject(0) { |sum, number| sum += number.abs}
    end

    return row_sums.max unless minus
    return row_sums.min
  end

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
    raise(DataTypeError, "Incompatible data type for the scaling factor") unless
        NMatrix::upcast(self.dtype, NMatrix::min_dtype(alpha)) == self.dtype
    return NMatrix::BLAS::scal(alpha, self, incx, self.size / incx) if NMatrix::BLAS.method_defined? :scal
    self.each_stored_with_indices do |e, *i|
      self[*i] = e*alpha
    end
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
    return self.clone.scale!(alpha, incx, n)
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
