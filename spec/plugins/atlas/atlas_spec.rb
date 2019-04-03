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
# == atlas_spec.rb
#
# Tests for interfaces that are only exposed by nmatrix-atlas
#

require "spec_helper"
require "./lib/nmatrix/atlas"

describe "NMatrix::LAPACK implementation from nmatrix-atlas plugin" do
  [:float32, :float64, :complex64, :complex128].each do |dtype|
    context dtype do
      it "exposes clapack_getri" do
        a = NMatrix.new(:dense, 3, [1, 0, 4, 1, 1, 6, -3, 0, -10], dtype)
        ipiv = NMatrix::LAPACK.clapack_getrf(:row, 3, 3, a, 3) # get pivot from getrf, use for getri

        begin
          NMatrix::LAPACK.clapack_getri(:row, 3, a, 3, ipiv)

          b = NMatrix.new(:dense, 3, [-5, 0, -2, -4, 1, -1, 1.5, 0, 0.5], dtype)
          expect(a).to eq(b)
        rescue NotImplementedError => e
          pending e.to_s
        end
      end

      # potrf decomposes a symmetric (or Hermitian)
      # positive-definite matrix. The matrix tested below isn't symmetric.
      # But this is okay since potrf just examines the upper/lower half
      # (as requested) of the matrix and assumes that the rest is symmetric,
      # so we just set the other part of the matrix to zero.
      it "exposes clapack_potrf upper" do
        pending "potrf requires clapack" unless NMatrix.has_clapack?

        a = NMatrix.new(:dense, 3, [25, 15, -5, 0, 18, 0, 0, 0, 11], dtype)
        NMatrix::LAPACK.clapack_potrf(:row, :upper, 3, a, 3)
        b = NMatrix.new(:dense, 3, [5, 3, -1, 0, 3, 1, 0, 0, 3], dtype)
        expect(a).to eq(b)
      end

      it "exposes clapack_potrf lower" do
        pending "potrf requires clapack" unless NMatrix.has_clapack?

        a = NMatrix.new(:dense, 3, [25, 0, 0, 15, 18, 0, -5, 0, 11], dtype)
        NMatrix::LAPACK.clapack_potrf(:row, :lower, 3, a, 3)
        b = NMatrix.new(:dense, 3, [5, 0, 0, 3, 3, 0, -1, 1, 3], dtype)
        expect(a).to eq(b)
      end

      it "exposes clapack_potri" do
        pending "potri requires clapack" unless NMatrix.has_clapack?

        a = NMatrix.new(3, [4, 0, -1,
                            0, 2, 1,
                            0, 0, 1,], dtype: dtype)
        NMatrix::LAPACK.clapack_potrf(:row, :upper, 3, a, 3)
        NMatrix::LAPACK.clapack_potri(:row, :upper, 3, a, 3)
        b = NMatrix.new(3, [0.5, -0.5, 1, 0, 1.5, -2,  0, 0, 4], dtype: dtype)
        err = case dtype
                when :float32, :complex64
                  1e-6
                when :float64, :complex128
                  1e-14
        end
        expect(a).to be_within(err).of(b)
      end

      it "exposes clapack_potrs" do
        pending "potrs requires clapack" unless NMatrix.has_clapack?

        a = NMatrix.new(3, [4, 0, -1,
                            0, 2, 1,
                            0, 0, 1,], dtype: dtype)
        b = NMatrix.new([3, 1], [3, 0, 2], dtype: dtype)

        NMatrix::LAPACK.clapack_potrf(:row, :upper, 3, a, 3)
        NMatrix::LAPACK.clapack_potrs(:row, :upper, 3, 1, a, 3, b, 3)

        x = NMatrix.new([3, 1], [3.5, -5.5, 11], dtype: dtype)

        err = case dtype
                when :float32, :complex64
                  1e-5
                when :float64, :complex128
                  1e-14
        end

        expect(b).to be_within(err).of(x)
      end
    end
  end

  [:float32, :float64, :complex64, :complex128].each do |dtype|
    context dtype do
      it "calculates the singular value decomposition with lapack_gesvd" do
        # example from Wikipedia
        m = 4
        n = 5
        mn_min = [m, n].min
        a = NMatrix.new([m, n], [1, 0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0], dtype: dtype)
        s = NMatrix.new([mn_min], 0, dtype: a.abs_dtype) # s is always real and always returned as float/double, never as complex
        u = NMatrix.new([m, m], 0, dtype: dtype)
        vt = NMatrix.new([n, n], 0, dtype: dtype)

        # This is a pure LAPACK function so it expects column-major functions
        # So we need to transpose the input as well as the output
        a = a.transpose
        NMatrix::LAPACK.lapack_gesvd(:a, :a, m, n, a, m, s, u, m, vt, n, 500)
        u = u.transpose
        vt = vt.transpose

        s_true = NMatrix.new([mn_min], [4, 3, Math.sqrt(5), 0], dtype: a.abs_dtype)
        u_true = NMatrix.new([m, m], [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0], dtype: dtype)
        vt_true = NMatrix.new([n, n], [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, Math.sqrt(0.2), 0, 0, 0, Math.sqrt(0.8), 0, 0, 0, 1, 0, -Math.sqrt(0.8), 0, 0, 0, Math.sqrt(0.2)], dtype: dtype)

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

      it "calculates the singular value decomposition with lapack_gesdd" do
        # example from Wikipedia
        m = 4
        n = 5
        mn_min = [m, n].min
        a = NMatrix.new([m, n], [1, 0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0], dtype: dtype)
        s = NMatrix.new([mn_min], 0, dtype: a.abs_dtype) # s is always real and always returned as float/double, never as complex
        u = NMatrix.new([m, m], 0, dtype: dtype)
        vt = NMatrix.new([n, n], 0, dtype: dtype)

        # This is a pure LAPACK function so it expects column-major functions
        # So we need to transpose the input as well as the output
        a = a.transpose
        NMatrix::LAPACK.lapack_gesdd(:a, m, n, a, m, s, u, m, vt, n, 500)
        u = u.transpose
        vt = vt.transpose

        s_true = NMatrix.new([mn_min], [4, 3, Math.sqrt(5), 0], dtype: a.abs_dtype)
        u_true = NMatrix.new([m, m], [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0], dtype: dtype)
        vt_true = NMatrix.new([n, n], [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, Math.sqrt(0.2), 0, 0, 0, Math.sqrt(0.8), 0, 0, 0, 1, 0, -Math.sqrt(0.8), 0, 0, 0, Math.sqrt(0.2)], dtype: dtype)

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

      it "exposes lapack_geev" do
        n = 3
        a = NMatrix.new([n, n], [-1, 0, 0, 0, 1, -2, 0, 1, -1], dtype: dtype)
        w = NMatrix.new([n], dtype: dtype)
        wi = if a.complex_dtype? # for real dtypes, imaginary parts of eigenvalues are stored in separate vector
          nil
        else
          NMatrix.new([n], dtype: dtype)
        end
        vl = NMatrix.new([n, n], dtype: dtype)
        vr = NMatrix.new([n, n], dtype: dtype)

        # This is a pure LAPACK routine so it expects column-major matrices,
        # so we need to transpose everything.
        a = a.transpose
        NMatrix::LAPACK.lapack_geev(:left, :right, n, a, n, w, wi, vl, n, vr, n, 2 * n)
        vr = vr.transpose
        vl = vl.transpose

        unless a.complex_dtype?
          w += wi * Complex(0, 1)
        end

        w_true = NMatrix.new([n], [Complex(0, 1), -Complex(0, 1), -1], dtype: NMatrix.upcast(dtype, :complex64))
        if a.complex_dtype?
          # For complex types the right/left eigenvectors are stored as columns
          # of vr/vl.
          vr_true = NMatrix.new([n, n], [0, 0, 1,
                                         2 / Math.sqrt(6), 2 / Math.sqrt(6), 0,
                                         Complex(1, -1) / Math.sqrt(6), Complex(1, 1) / Math.sqrt(6), 0,], dtype: dtype)
          vl_true = NMatrix.new([n, n], [0, 0, 1,
                                         Complex(-1, 1) / Math.sqrt(6), Complex(-1, -1) / Math.sqrt(6), 0,
                                         2 / Math.sqrt(6), 2 / Math.sqrt(6), 0,], dtype: dtype)
        else
          # For real types, the real part of the first and second eigenvectors is
          # stored in the first column, the imaginary part of the first (= the
          # negative of the imaginary part of the second) eigenvector is stored
          # in the second column, and the third eigenvector (purely real) is the
          # third column.
          vr_true = NMatrix.new([n, n], [0, 0, 1,
                                         2 / Math.sqrt(6), 0, 0,
                                         1 / Math.sqrt(6), -1 / Math.sqrt(6), 0,], dtype: dtype)
          vl_true = NMatrix.new([n, n], [0, 0, 1,
                                         -1 / Math.sqrt(6), 1 / Math.sqrt(6), 0,
                                         2 / Math.sqrt(6), 0, 0,], dtype: dtype)
        end

        err = case dtype
                when :float32, :complex64
                  1e-6
                when :float64, :complex128
                  1e-15
        end

        expect(w).to be_within(err).of(w_true)
        expect(vr).to be_within(err).of(vr_true)
        expect(vl).to be_within(err).of(vl_true)
      end
    end
  end
end
