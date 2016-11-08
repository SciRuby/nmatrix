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
# SciRuby is Copyright (c) 2010 - 2016, Ruby Science Foundation
# NMatrix is Copyright (c) 2012 - 2016, John Woods and the Ruby Science Foundation
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
# == blas.rb
#
# This file contains the safer accessors for the BLAS functions
# supported by NMatrix.
#++

module NMatrix::BLAS

  #Add functions from C extension to main BLAS module
  class << self
    if jruby?
      # BLAS functionalities for JRuby need to be implemented
    else
      NMatrix::Internal::BLAS.singleton_methods.each do |m|
        define_method m, NMatrix::Internal::BLAS.method(m).to_proc
      end
    end
  end

  class << self
    #
    # call-seq:
    #     gemm(a, b) -> NMatrix
    #     gemm(a, b, c) -> NMatrix
    #     gemm(a, b, c, alpha, beta) -> NMatrix
    #
    # Updates the value of C via the matrix multiplication
    #   C = (alpha * A * B) + (beta * C)
    # where +alpha+ and +beta+ are scalar values.
    #
    # * *Arguments* :
    #   - +a+ -> Matrix A.
    #   - +b+ -> Matrix B.
    #   - +c+ -> Matrix C.
    #   - +alpha+ -> A scalar value that multiplies A * B.
    #   - +beta+ -> A scalar value that multiplies C.
    #   - +transpose_a+ ->
    #   - +transpose_b+ ->
    #   - +m+ ->
    #   - +n+ ->
    #   - +k+ ->
    #   - +lda+ ->
    #   - +ldb+ ->
    #   - +ldc+ ->
    # * *Returns* :
    #   - A NMatrix equal to (alpha * A * B) + (beta * C).
    # * *Raises* :
    #   - +ArgumentError+ -> +a+ and +b+ must be dense matrices.
    #   - +ArgumentError+ -> +c+ must be +nil+ or a dense matrix.
    #   - +ArgumentError+ -> The dtype of the matrices must be equal.
    #
    def gemm(a, b, c = nil, alpha = 1.0, beta = 0.0,
             transpose_a = false, transpose_b = false, m = nil,
             n = nil, k = nil, lda = nil, ldb = nil, ldc = nil)

      raise(ArgumentError, 'Expected dense NMatrices as first two arguments.') \
            unless a.is_a?(NMatrix) and b.is_a? \
            (NMatrix) and a.stype == :dense and b.stype == :dense

      raise(ArgumentError, 'Expected nil or dense NMatrix as third argument.') \
            unless c.nil? or (c.is_a?(NMatrix)  \
            and c.stype == :dense)
      raise(ArgumentError, 'NMatrix dtype mismatch.') \
            unless a.dtype == b.dtype and (c ? a.dtype == c.dtype : true)

      # First, set m, n, and k, which depend on whether we're taking the
      # transpose of a and b.
      if c
        m ||= c.shape[0]
        n ||= c.shape[1]
        k ||= transpose_a ? a.shape[0] : a.shape[1]

      else
        if transpose_a
          # Either :transpose or :complex_conjugate.
          m ||= a.shape[1]
          k ||= a.shape[0]

        else
          # No transpose.
          m ||= a.shape[0]
          k ||= a.shape[1]
        end

        n ||= transpose_b ? b.shape[0] : b.shape[1]
        c  = NMatrix.new([m, n], dtype: a.dtype)
      end

      # I think these are independent of whether or not a transpose occurs.
      lda ||= a.shape[1]
      ldb ||= b.shape[1]
      ldc ||= c.shape[1]

      # NM_COMPLEX64 and NM_COMPLEX128 both require complex alpha and beta.
      if a.dtype == :complex64 or a.dtype == :complex128
        alpha = Complex(1.0, 0.0) if alpha == 1.0
        beta  = Complex(0.0, 0.0) if beta  == 0.0
      end

      # For argument descriptions, see: http://www.netlib.org/blas/dgemm.f
      ::NMatrix::BLAS.cblas_gemm(:row, transpose_a, transpose_b,
       m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)

      return c
    end

    #
    # call-seq:
    #     gemv(a, x) -> NMatrix
    #     gemv(a, x, y) -> NMatrix
    #     gemv(a, x, y, alpha, beta) -> NMatrix
    #
    # Implements matrix-vector product via
    #   y = (alpha * A * x) + (beta * y)
    # where +alpha+ and +beta+ are scalar values.
    #
    # * *Arguments* :
    #   - +a+ -> Matrix A.
    #   - +x+ -> Vector x.
    #   - +y+ -> Vector y.
    #   - +alpha+ -> A scalar value that multiplies A * x.
    #   - +beta+ -> A scalar value that multiplies y.
    #   - +transpose_a+ ->
    #   - +m+ ->
    #   - +n+ ->
    #   - +lda+ ->
    #   - +incx+ ->
    #   - +incy+ ->
    # * *Returns* :
    #   -
    # * *Raises* :
    #   - ++ ->
    #
    def gemv(a, x, y = nil, alpha = 1.0, beta = 0.0,
             transpose_a = false, m = nil, n = nil, lda = nil,
             incx = nil, incy = nil)
      raise(ArgumentError, 'Expected dense NMatrices as first two arguments.') \
       unless a.is_a?(NMatrix) and x.is_a?(NMatrix) and \
       a.stype == :dense and x.stype == :dense

      raise(ArgumentError, 'Expected nil or dense NMatrix as third argument.') \
       unless y.nil? or (y.is_a?(NMatrix) and y.stype == :dense)

      raise(ArgumentError, 'NMatrix dtype mismatch.') \
       unless a.dtype == x.dtype and (y ? a.dtype == y.dtype : true)

      m ||= transpose_a == :transpose ? a.shape[1] : a.shape[0]
      n ||= transpose_a == :transpose ? a.shape[0] : a.shape[1]
      raise(ArgumentError, "dimensions don't match") \
       unless x.shape[0] == n && x.shape[1] == 1

      if y
        raise(ArgumentError, "dimensions don't match") \
         unless y.shape[0] == m && y.shape[1] == 1
      else
        y = NMatrix.new([m,1], dtype: a.dtype)
      end

      lda  ||= a.shape[1]
      incx ||= 1
      incy ||= 1

      ::NMatrix::BLAS.cblas_gemv(transpose_a, m, n,
       alpha, a, lda, x, incx, beta, y, incy)

      return y
    end

    #
    # call-seq:
    #     rot(x, y, c, s) -> [NMatrix, NMatrix]
    #
    # Apply plane rotation.
    #
    # * *Arguments* :
    #   - +x+ -> NMatrix
    #   - +y+ -> NMatrix
    #   - +c+ -> cosine of the angle of rotation
    #   - +s+ -> sine of the angle of rotation
    #   - +incx+ -> stride of NMatrix +x+
    #   - +incy+ -> stride of NMatrix +y+
    #   - +n+ -> number of elements to consider in x and y
    #   - +in_place+ -> true   if it's okay to modify the supplied
    #                           +x+ and +y+ parameters directly;
    #                   false if not. Default is false.
    # * *Returns* :
    #   - Array with the results, in the format [xx, yy]
    # * *Raises* :
    #   - +ArgumentError+ -> Expected dense NMatrices as first two arguments.
    #   - +ArgumentError+ -> NMatrix dtype mismatch.
    #   - +ArgumentError+ -> Need to supply n for non-standard incx,
    #                         incy values.
    #
    def rot(x, y, c, s, incx = 1, incy = 1, n = nil, in_place=false)
      raise(ArgumentError, 'Expected dense NMatrices as first two arguments.') \
       unless x.is_a?(NMatrix) and y.is_a?(NMatrix) \
       and x.stype == :dense and y.stype == :dense

      raise(ArgumentError, 'NMatrix dtype mismatch.') \
       unless x.dtype == y.dtype

      raise(ArgumentError, 'Need to supply n for non-standard incx, incy values') \
       if n.nil? && incx != 1 && incx != -1 && incy != 1 && incy != -1

      n ||= [x.size/incx.abs, y.size/incy.abs].min

      if in_place
        ::NMatrix::BLAS.cblas_rot(n, x, incx, y, incy, c, s)
        return [x,y]
      else
        xx = x.clone
        yy = y.clone

        ::NMatrix::BLAS.cblas_rot(n, xx, incx, yy, incy, c, s)

        return [xx,yy]
      end
    end


    #
    # call-seq:
    #     rot!(x, y, c, s) -> [NMatrix, NMatrix]
    #
    # Apply plane rotation directly to +x+ and +y+.
    #
    # See rot for arguments.
    def rot!(x, y, c, s, incx = 1, incy = 1, n = nil)
      rot(x,y,c,s,incx,incy,n,true)
    end


    #
    # call-seq:
    #     rotg(ab) -> [Numeric, Numeric]
    #
    # Apply givens plane rotation to the coordinates (a,b),
    #  returning the cosine and sine of the angle theta.
    #
    # Since the givens rotation includes a square root,
    #  integers are disallowed.
    #
    # * *Arguments* :
    #   - +ab+ -> NMatrix with two elements
    # * *Returns* :
    #   - Array with the results, in the format [cos(theta), sin(theta)]
    # * *Raises* :
    #   - +ArgumentError+ -> Expected dense NMatrix of size 2
    #
    def rotg(ab)
      raise(ArgumentError, "Expected dense NMatrix of shape [2,1] or [1,2]") \
       unless ab.is_a?(NMatrix) && ab.stype == :dense && ab.size == 2

      ::NMatrix::BLAS.cblas_rotg(ab)
    end


    #
    # call-seq:
    #     asum(x, incx, n) -> Numeric
    #
    # Calculate the sum of absolute values of the entries of a
    #  vector +x+ of size +n+
    #
    # * *Arguments* :
    #   - +x+ -> an NMatrix (will also allow an NMatrix,
    #             but will treat it as if it's a vector )
    #   - +incx+ -> the skip size (defaults to 1)
    #   - +n+ -> the size of +x+ (defaults to +x.size / incx+)
    # * *Returns* :
    #   - The sum
    # * *Raises* :
    #   - +ArgumentError+ -> Expected dense NMatrix for arg 0
    #   - +RangeError+ -> n out of range
    #
    def asum(x, incx = 1, n = nil)
      n ||= x.size / incx
      raise(ArgumentError, "Expected dense NMatrix for arg 0") \
       unless x.is_a?(NMatrix)

      raise(RangeError, "n out of range") \
       if n*incx > x.size || n*incx <= 0 || n <= 0
       ::NMatrix::BLAS.cblas_asum(n, x, incx)
    end

    #
    # call-seq:
    #     nrm2(x, incx, n)
    #
    # Calculate the 2-norm of a vector +x+ of size +n+
    #
    # * *Arguments* :
    #   - +x+ -> an NMatrix (will also allow an
    #             NMatrix, but will treat it as if it's a vector )
    #   - +incx+ -> the skip size (defaults to 1)
    #   - +n+ -> the size of +x+ (defaults to +x.size / incx+)
    # * *Returns* :
    #   - The 2-norm
    # * *Raises* :
    #   - +ArgumentError+ -> Expected dense NMatrix for arg 0
    #   - +RangeError+ -> n out of range
    #
    def nrm2(x, incx = 1, n = nil)
      n ||= x.size / incx
      raise(ArgumentError, "Expected dense NMatrix for arg 0") \
       unless x.is_a?(NMatrix)

      raise(RangeError, "n out of range") \
       if n*incx > x.size || n*incx <= 0 || n <= 0
       ::NMatrix::BLAS.cblas_nrm2(n, x, incx)
    end

    #
    # call-seq:
    #     scal(alpha, vector, incx, n)
    #
    # Scale a matrix by a given scaling factor
    #
    # * *Arguments* :
    #   - +alpha+ -> a scaling factor
    #   - +vector+ -> an NMatrix
    #   - +incx+ -> the skip size (defaults to 1)
    #   - +n+ -> the size of +x+ (defaults to +x.size / incx+)
    # * *Returns* :
    #   - The scaling result
    # * *Raises* :
    #   - +ArgumentError+ -> Expected dense NMatrix for arg 0
    #   - +RangeError+ -> n out of range
    #
    def scal(alpha, vector, incx=1, n=nil)
      n ||= vector.size / incx
      raise(ArgumentError, "Expected dense NMatrix for arg 0") unless vector.is_a?(NMatrix)
      raise(RangeError, "n out of range") if n*incx > vector.size || n*incx <= 0 || n <= 0
      ::NMatrix::BLAS.cblas_scal(n, alpha, vector, incx)
    end

    # The following are functions that used to be implemented in C, but
    # now require nmatrix-atlas or nmatrix-lapcke to run properly, so we can just
    # implemented their stubs in Ruby.
    def cblas_trmm(order, side, uplo, trans_a, diag, m, n, alpha, a, lda, b, ldb)
      raise(NotImplementedError,"cblas_trmm requires either the
       nmatrix-lapacke or nmatrix-atlas gem")
    end

    def cblas_syrk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      raise(NotImplementedError,"cblas_syrk requires either the
       nmatrix-lapacke or nmatrix-atlas gem")
    end

    def cblas_herk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      raise(NotImplementedError,"cblas_herk requires either the
       nmatrix-lapacke or nmatrix-atlas gem")
    end
  end
end
