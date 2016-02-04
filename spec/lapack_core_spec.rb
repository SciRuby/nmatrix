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
# == lapack_core_spec.rb
#
# Tests for LAPACK functions that have internal implementations (i.e. they
# don't rely on external libraries) and also functions that are implemented
# by both nmatrix-atlas and nmatrix-lapacke. These tests will also be run for the
# plugins that do use external libraries, since they will override the
# internal implmentations.
#

require 'spec_helper'

describe "NMatrix::LAPACK functions with internal implementations" do
  # where integer math is allowed
  [:byte, :int8, :int16, :int32, :int64, :float32, :float64, :complex64, :complex128].each do |dtype|
    context dtype do
      # This spec seems a little weird. It looks like laswp ignores the last
      # element of piv, though maybe I misunderstand smth. It would make
      # more sense if piv were [2,1,3,3]
      it "exposes clapack laswp" do
        a = NMatrix.new(:dense, [3,4], [1,2,3,4,5,6,7,8,9,10,11,12], dtype)
        NMatrix::LAPACK::clapack_laswp(3, a, 4, 0, 3, [2,1,3,0], 1)
        b = NMatrix.new(:dense, [3,4], [3,2,4,1,7,6,8,5,11,10,12,9], dtype)
        expect(a).to eq(b)
      end

      # This spec is OK, because the default behavior for permute_columns
      # is :intuitive, which is different from :lapack (default laswp behavior)
      it "exposes NMatrix#permute_columns and #permute_columns! (user-friendly laswp)" do
        a = NMatrix.new(:dense, [3,4], [1,2,3,4,5,6,7,8,9,10,11,12], dtype)
        b = NMatrix.new(:dense, [3,4], [3,2,4,1,7,6,8,5,11,10,12,9], dtype)
        piv = [2,1,3,0]
        r = a.permute_columns(piv)
        expect(r).not_to eq(a)
        expect(r).to eq(b)
        a.permute_columns!(piv)
        expect(a).to eq(b)
      end
    end
  end

  # where integer math is not allowed
  [:float32, :float64, :complex64, :complex128].each do |dtype|
    context dtype do

      # clapack_getrf performs a LU decomposition, but unlike the
      # standard LAPACK getrf, it's the upper matrix that has unit diagonals
      # and the permutation is done in columns not rows. See the code for
      # details.
      # Also the rows in the pivot vector are indexed starting from 0,
      # rather than 1 as in LAPACK
      it "calculates LU decomposition using clapack_getrf (row-major, square)" do
        a = NMatrix.new(3, [4,9,2,3,5,7,8,1,6], dtype: dtype)
        ipiv = NMatrix::LAPACK::clapack_getrf(:row, a.shape[0], a.shape[1], a, a.shape[1])
        b = NMatrix.new(3,[9, 2.0/9, 4.0/9,
                           5, 53.0/9, 7.0/53,
                           1, 52.0/9, 360.0/53], dtype: dtype)
        ipiv_true = [1,2,2]

        # delta varies for different dtypes
        err = case dtype
                when :float32, :complex64
                  1e-6
                when :float64, :complex128
                  1e-15
              end

        expect(a).to be_within(err).of(b)
        expect(ipiv).to eq(ipiv_true)
      end

      it "calculates LU decomposition using clapack_getrf (row-major, rectangular)" do
        a = NMatrix.new([3,4], GETRF_EXAMPLE_ARRAY, dtype: dtype)
        ipiv = NMatrix::LAPACK::clapack_getrf(:row, a.shape[0], a.shape[1], a, a.shape[1])
        #we can't use GETRF_SOLUTION_ARRAY here, because of the different
        #conventions of clapack_getrf
        b = NMatrix.new([3,4],[10.0, -0.1,      0.0,       0.4,
                               3.0,   9.3,  20.0/93,   38.0/93,
                               1.0,   7.1, 602.0/93, 251.0/602], dtype: dtype)
        ipiv_true = [2,2,2]

        # delta varies for different dtypes
        err = case dtype
                when :float32, :complex64
                  1e-6
                when :float64, :complex128
                  1e-15
              end

        expect(a).to be_within(err).of(b)
        expect(ipiv).to eq(ipiv_true)
      end

      #Normally we wouldn't check column-major routines, since all our matrices
      #are row-major, but we use the column-major version in #getrf!, so we
      #want to test it here.
      it "calculates LU decomposition using clapack_getrf (col-major, rectangular)" do
        #this is supposed to represent the 3x2 matrix
        # -1  2
        #  0  3
        #  1 -2
        a = NMatrix.new([1,6], [-1,0,1,2,3,-2], dtype: dtype)
        ipiv = NMatrix::LAPACK::clapack_getrf(:col, 3, 2, a, 3)
        b = NMatrix.new([1,6], [-1,0,-1,2,3,0], dtype: dtype)
        ipiv_true = [0,1]

        # delta varies for different dtypes
        err = case dtype
                when :float32, :complex64
                  1e-6
                when :float64, :complex128
                  1e-15
              end

        expect(a).to be_within(err).of(b)
        expect(ipiv).to eq(ipiv_true)
      end

      it "calculates LU decomposition using #getrf! (rectangular)" do
        a = NMatrix.new([3,4], GETRF_EXAMPLE_ARRAY, dtype: dtype)
        ipiv = a.getrf!
        b = NMatrix.new([3,4], GETRF_SOLUTION_ARRAY, dtype: dtype)
        ipiv_true = [2,3,3]

        # delta varies for different dtypes
        err = case dtype
                when :float32, :complex64
                  1e-6
                when :float64, :complex128
                  1e-14
              end

        expect(a).to be_within(err).of(b)
        expect(ipiv).to eq(ipiv_true)
      end

      it "calculates LU decomposition using #getrf! (square)" do
        a = NMatrix.new([4,4], [0,1,2,3, 1,1,1,1, 0,-1,-2,0, 0,2,0,2], dtype: dtype)
        ipiv = a.getrf!

        b = NMatrix.new([4,4], [1,1,1,1, 0,2,0,2, 0,-0.5,-2,1, 0,0.5,-1,3], dtype: dtype)
        ipiv_true = [2,4,3,4]

        expect(a).to eq(b)
        expect(ipiv).to eq(ipiv_true)
      end

      # Together, these calls are basically xGESV from LAPACK: http://www.netlib.org/lapack/double/dgesv.f
      it "exposes clapack_getrs" do
        a     = NMatrix.new(3, [-2,4,-3, 3,-2,1, 0,-4,3], dtype: dtype)
        ipiv  = NMatrix::LAPACK::clapack_getrf(:row, 3, 3, a, 3)
        b     = NMatrix.new([3,1], [-1, 17, -9], dtype: dtype)

        NMatrix::LAPACK::clapack_getrs(:row, false, 3, 1, a, 3, ipiv, b, 3)

        expect(b[0]).to eq(5)
        expect(b[1]).to eq(-15.0/2)
        expect(b[2]).to eq(-13)
      end

      it "solves matrix equation (non-vector rhs) using clapack_getrs" do
        a     = NMatrix.new(3, [-2,4,-3, 3,-2,1, 0,-4,3], dtype: dtype)
        b     = NMatrix.new([3,2], [-1,2, 17,1, -9,-4], dtype: dtype)

        n = a.shape[0]
        nrhs = b.shape[1]

        ipiv  = NMatrix::LAPACK::clapack_getrf(:row, n, n, a, n)
        # Even though we pass :row to clapack_getrs, it still interprets b as
        # column-major, so need to transpose b before and after:
        b = b.transpose
        NMatrix::LAPACK::clapack_getrs(:row, false, n, nrhs, a, n, ipiv, b, n)
        b = b.transpose

        b_true = NMatrix.new([3,2], [5,1, -7.5,1, -13,0], dtype: dtype)
        expect(b).to eq(b_true)
      end

      #posv is like potrf+potrs
      #posv is implemented in both nmatrix-atlas and nmatrix-lapacke, so the spec
      #needs to be shared here
      it "solves a (symmetric positive-definite) matrix equation using posv (vector rhs)" do
        a = NMatrix.new(3, [4, 0,-1,
                            0, 2, 1,
                            0, 0, 1], dtype: dtype)
        b = NMatrix.new([3,1], [4,2,0], dtype: dtype)

        begin
          x = NMatrix::LAPACK::posv(:upper, a, b)
        rescue NotImplementedError => e
          pending e.to_s
        end

        x_true = NMatrix.new([3,1], [1, 1, 0], dtype: dtype)

        err = case dtype
                when :float32, :complex64
                  1e-5
                when :float64, :complex128
                  1e-14
              end

        expect(x).to be_within(err).of(x_true)
      end

      it "solves a (symmetric positive-definite) matrix equation using posv (non-vector rhs)" do
        a = NMatrix.new(3, [4, 0,-1,
                            0, 2, 1,
                            0, 0, 1], dtype: dtype)
        b = NMatrix.new([3,2], [4,-1, 2,-1, 0,0], dtype: dtype)

        begin
          x = NMatrix::LAPACK::posv(:upper, a, b)
        rescue NotImplementedError => e
          pending e.to_s
        end

        x_true = NMatrix.new([3,2], [1,0, 1,-1, 0,1], dtype: dtype)

        err = case dtype
                when :float32, :complex64
                  1e-5
                when :float64, :complex128
                  1e-14
              end

        expect(x).to be_within(err).of(x_true)
      end

      it "calculates the singular value decomposition with NMatrix#gesvd" do
        #example from Wikipedia
        m = 4
        n = 5
        mn_min = [m,n].min
        a = NMatrix.new([m,n],[1,0,0,0,2, 0,0,3,0,0, 0,0,0,0,0, 0,4,0,0,0], dtype: dtype)

        begin
          u, s, vt = a.gesvd
        rescue NotImplementedError => e
          pending e.to_s
        end

        s_true = NMatrix.new([mn_min,1], [4,3,Math.sqrt(5),0], dtype: a.abs_dtype)
        u_true = NMatrix.new([m,m], [0,0,1,0, 0,1,0,0, 0,0,0,-1, 1,0,0,0], dtype: dtype)
        vt_true = NMatrix.new([n,n], [0,1,0,0,0, 0,0,1,0,0, Math.sqrt(0.2),0,0,0,Math.sqrt(0.8), 0,0,0,1,0, -Math.sqrt(0.8),0,0,0,Math.sqrt(0.2)], dtype: dtype)

        err = case dtype
                when :float32, :complex64
                  1e-5
                when :float64, :complex128
                  1e-14
              end

        expect(s).to be_within(err).of(s_true)
        expect(u).to be_within(err).of(u_true)
        expect(vt).to be_within(err).of(vt_true)

        expect(s.dtype).to eq(a.abs_dtype)
        expect(u.dtype).to eq(dtype)
        expect(vt.dtype).to eq(dtype)
      end

      it "calculates the singular value decomposition with NMatrix#gesdd" do
        #example from Wikipedia
        m = 4
        n = 5
        mn_min = [m,n].min
        a = NMatrix.new([m,n],[1,0,0,0,2, 0,0,3,0,0, 0,0,0,0,0, 0,4,0,0,0], dtype: dtype)

        begin
          u, s, vt = a.gesdd
        rescue NotImplementedError => e
          pending e.to_s
        end

        s_true = NMatrix.new([mn_min,1], [4,3,Math.sqrt(5),0], dtype: a.abs_dtype)
        u_true = NMatrix.new([m,m], [0,0,1,0, 0,1,0,0, 0,0,0,-1, 1,0,0,0], dtype: dtype)
        vt_true = NMatrix.new([n,n], [0,1,0,0,0, 0,0,1,0,0, Math.sqrt(0.2),0,0,0,Math.sqrt(0.8), 0,0,0,1,0, -Math.sqrt(0.8),0,0,0,Math.sqrt(0.2)], dtype: dtype)

        err = case dtype
                when :float32, :complex64
                  1e-5
                when :float64, :complex128
                  1e-14
              end

        expect(s).to be_within(err).of(s_true)
        expect(u).to be_within(err).of(u_true)
        expect(vt).to be_within(err).of(vt_true)
      end


      it "calculates eigenvalues and eigenvectors NMatrix::LAPACK.geev (real matrix, complex eigenvalues)" do
        n = 3
        a = NMatrix.new([n,n], [-1,0,0, 0,1,-2, 0,1,-1], dtype: dtype)

        begin
          eigenvalues, vl, vr = NMatrix::LAPACK.geev(a)
        rescue NotImplementedError => e
          pending e.to_s
        end

        eigenvalues_true = NMatrix.new([n,1], [Complex(0,1), -Complex(0,1), -1], dtype: NMatrix.upcast(dtype, :complex64))
        vr_true = NMatrix.new([n,n],[0,0,1,
                                     2/Math.sqrt(6),2/Math.sqrt(6),0,
                                     Complex(1,-1)/Math.sqrt(6),Complex(1,1)/Math.sqrt(6),0], dtype: NMatrix.upcast(dtype, :complex64))
        vl_true = NMatrix.new([n,n],[0,0,1,
                                     Complex(-1,1)/Math.sqrt(6),Complex(-1,-1)/Math.sqrt(6),0,
                                     2/Math.sqrt(6),2/Math.sqrt(6),0], dtype: NMatrix.upcast(dtype, :complex64))

        err = case dtype
                when :float32, :complex64
                  1e-6
                when :float64, :complex128
                  1e-15
              end

        expect(eigenvalues).to be_within(err).of(eigenvalues_true)
        expect(vr).to be_within(err).of(vr_true)
        expect(vl).to be_within(err).of(vl_true)

        expect(eigenvalues.dtype).to eq(NMatrix.upcast(dtype, :complex64))
        expect(vr.dtype).to eq(NMatrix.upcast(dtype, :complex64))
        expect(vl.dtype).to eq(NMatrix.upcast(dtype, :complex64))
      end

      it "calculates eigenvalues and eigenvectors NMatrix::LAPACK.geev (real matrix, real eigenvalues)" do
        n = 3
        a = NMatrix.new([n,n], [2,0,0, 0,3,2, 0,1,2], dtype: dtype)

        begin
          eigenvalues, vl, vr = NMatrix::LAPACK.geev(a)
        rescue NotImplementedError => e
          pending e.to_s
        end

        eigenvalues_true = NMatrix.new([n,1], [1, 4, 2], dtype: dtype)

        # For some reason, some of the eigenvectors have different signs
        # when we use the complex versions of geev. This is totally fine, since
        # they are still normalized eigenvectors even with the sign flipped.
        if a.complex_dtype?
          vr_true = NMatrix.new([n,n],[0,0,1,
                                       1/Math.sqrt(2),2/Math.sqrt(5),0,
                                       -1/Math.sqrt(2),1/Math.sqrt(5),0], dtype: dtype)
          vl_true = NMatrix.new([n,n],[0,0,1,
                                       -1/Math.sqrt(5),1/Math.sqrt(2),0,
                                       2/Math.sqrt(5),1/Math.sqrt(2),0], dtype: dtype)
        else
          vr_true = NMatrix.new([n,n],[0,0,1,
                                       1/Math.sqrt(2),-2/Math.sqrt(5),0,
                                       -1/Math.sqrt(2),-1/Math.sqrt(5),0], dtype: dtype)
          vl_true = NMatrix.new([n,n],[0,0,1,
                                       1/Math.sqrt(5),-1/Math.sqrt(2),0,
                                       -2/Math.sqrt(5),-1/Math.sqrt(2),0], dtype: dtype)
        end

        err = case dtype
                when :float32, :complex64
                  1e-6
                when :float64, :complex128
                  1e-15
              end

        expect(eigenvalues).to be_within(err).of(eigenvalues_true)
        expect(vr).to be_within(err).of(vr_true)
        expect(vl).to be_within(err).of(vl_true)

        expect(eigenvalues.dtype).to eq(dtype)
        expect(vr.dtype).to eq(dtype)
        expect(vl.dtype).to eq(dtype)
      end

      it "calculates eigenvalues and eigenvectors NMatrix::LAPACK.geev (left eigenvectors only)" do
        n = 3
        a = NMatrix.new([n,n], [-1,0,0, 0,1,-2, 0,1,-1], dtype: dtype)

        begin
          eigenvalues, vl = NMatrix::LAPACK.geev(a, :left)
        rescue NotImplementedError => e
          pending e.to_s
        end

        eigenvalues_true = NMatrix.new([n,1], [Complex(0,1), -Complex(0,1), -1], dtype: NMatrix.upcast(dtype, :complex64))
        vl_true = NMatrix.new([n,n],[0,0,1,
                                     Complex(-1,1)/Math.sqrt(6),Complex(-1,-1)/Math.sqrt(6),0,
                                     2/Math.sqrt(6),2/Math.sqrt(6),0], dtype: NMatrix.upcast(dtype, :complex64))

        err = case dtype
                when :float32, :complex64
                  1e-6
                when :float64, :complex128
                  1e-15
              end

        expect(eigenvalues).to be_within(err).of(eigenvalues_true)
        expect(vl).to be_within(err).of(vl_true)
      end

      it "calculates eigenvalues and eigenvectors NMatrix::LAPACK.geev (right eigenvectors only)" do
        n = 3
        a = NMatrix.new([n,n], [-1,0,0, 0,1,-2, 0,1,-1], dtype: dtype)

        begin
          eigenvalues, vr = NMatrix::LAPACK.geev(a, :right)
        rescue NotImplementedError => e
          pending e.to_s
        end

        eigenvalues_true = NMatrix.new([n,1], [Complex(0,1), -Complex(0,1), -1], dtype: NMatrix.upcast(dtype, :complex64))
        vr_true = NMatrix.new([n,n],[0,0,1,
                                     2/Math.sqrt(6),2/Math.sqrt(6),0,
                                     Complex(1,-1)/Math.sqrt(6),Complex(1,1)/Math.sqrt(6),0], dtype: NMatrix.upcast(dtype, :complex64))

        err = case dtype
                when :float32, :complex64
                  1e-6
                when :float64, :complex128
                  1e-15
              end

        expect(eigenvalues).to be_within(err).of(eigenvalues_true)
        expect(vr).to be_within(err).of(vr_true)
      end
    end
  end

  [:complex64, :complex128].each do |dtype|
    context dtype do
      it "calculates eigenvalues and eigenvectors NMatrix::LAPACK.geev (complex matrix)" do
        n = 3
        a = NMatrix.new([n,n], [Complex(0,1),0,0, 0,3,2, 0,1,2], dtype: dtype)

        begin
          eigenvalues, vl, vr = NMatrix::LAPACK.geev(a)
        rescue NotImplementedError => e
          pending e.to_s
        end

        eigenvalues_true = NMatrix.new([n,1], [1, 4, Complex(0,1)], dtype: dtype)
        vr_true = NMatrix.new([n,n],[0,0,1,
                                     1/Math.sqrt(2),2/Math.sqrt(5),0,
                                     -1/Math.sqrt(2),1/Math.sqrt(5),0], dtype: dtype)
        vl_true = NMatrix.new([n,n],[0,0,1,
                                     -1/Math.sqrt(5),1/Math.sqrt(2),0,
                                     2/Math.sqrt(5),1/Math.sqrt(2),0], dtype: dtype)

        err = case dtype
                when :float32, :complex64
                  1e-6
                when :float64, :complex128
                  1e-15
              end

        expect(eigenvalues).to be_within(err).of(eigenvalues_true)
        expect(vr).to be_within(err).of(vr_true)
        expect(vl).to be_within(err).of(vl_true)
      end
    end
  end
end
