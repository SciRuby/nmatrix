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
# == lapacke_spec.rb
#
# Tests for interfaces that are only exposed by nmatrix-lapacke
#

require 'spec_helper'
require "./lib/nmatrix/lapacke"

describe "NMatrix::LAPACK functions implemented with LAPACKE interface" do
  [:float32, :float64, :complex64, :complex128].each do |dtype|
    context dtype do
      it "exposes lapacke_getrf" do
        a = NMatrix.new([3,4], GETRF_EXAMPLE_ARRAY, dtype: dtype)
        ipiv = NMatrix::LAPACK.lapacke_getrf(:row, 3, 4, a, 4)
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

      it "exposes lapacke_getri" do
        a = NMatrix.new(:dense, 3, [1,0,4,1,1,6,-3,0,-10], dtype)
        ipiv = NMatrix::LAPACK::lapacke_getrf(:row, 3, 3, a, 3) # get pivot from getrf, use for getri

        # delta varies for different dtypes
        err = case dtype
                when :float32, :complex64
                  1e-5
                when :float64, :complex128
                  1e-14
              end

        NMatrix::LAPACK::lapacke_getri(:row, 3, a, 3, ipiv)

        b = NMatrix.new(:dense, 3, [-5,0,-2,-4,1,-1,1.5,0,0.5], dtype)
        expect(a).to be_within(err).of(b)
      end

      it "exposes lapacke_getrs with vector solutions" do
        a     = NMatrix.new(3, [-2,4,-3,3,-2,1,0,-4,3], dtype: dtype)
        ipiv  = NMatrix::LAPACK::lapacke_getrf(:row, 3, 3, a, 3)
        b     = NMatrix.new([3,1], [-1, 17, -9], dtype: dtype)

        #be careful! the leading dimenension (lda,ldb) is the number of rows for row-major in LAPACKE. Different from CLAPACK convention!
        NMatrix::LAPACK::lapacke_getrs(:row, false, 3, 1, a, 3, ipiv, b, 1)

        # delta varies for different dtypes
        err = case dtype
                when :float32, :complex64
                  1e-5
                when :float64, :complex128
                  1e-13
              end

        expect(b[0]).to be_within(err).of(5)
        expect(b[1]).to be_within(err).of(-15.0/2)
        expect(b[2]).to be_within(err).of(-13)
      end

      it "exposes lapacke_getrs with matrix solutions" do
        a     = NMatrix.new(3, [-2,4,-3,3,-2,1,0,-4,3], dtype: dtype)
        ipiv  = NMatrix::LAPACK::lapacke_getrf(:row, 3, 3, a, 3)
        b     = NMatrix.new([3,2], [-1, 2, 17, 10, -9, 1], dtype: dtype)

        #be careful! the leading dimenension (lda,ldb) is the number of rows for row-major in LAPACKE. Different from CLAPACK convention!
        NMatrix::LAPACK::lapacke_getrs(:row, false, 3, 2, a, 3, ipiv, b, 2)

        # delta varies for different dtypes
        err = case dtype
                when :float32, :complex64
                  1e-4
                when :float64, :complex128
                  1e-13
              end

        x = NMatrix.new([3,2], [5, -1.5, -7.5, -21.25, -13, -28], dtype: dtype)
        expect(b).to be_within(err).of(x)
      end

      it "exposes lapacke_potrf" do
        # first do upper
        begin
          a = NMatrix.new(:dense, 3, [25,15,-5, 0,18,0, 0,0,11], dtype)
          NMatrix::LAPACK::lapacke_potrf(:row, :upper, 3, a, 3)
          b = NMatrix.new(:dense, 3, [5,3,-1, 0,3,1, 0,0,3], dtype)
          expect(a).to eq(b)
        end

        # then do lower
        a = NMatrix.new(:dense, 3, [25,0,0, 15,18,0,-5,0,11], dtype)
        NMatrix::LAPACK::lapacke_potrf(:row, :lower, 3, a, 3)
        b = NMatrix.new(:dense, 3, [5,0,0, 3,3,0, -1,1,3], dtype)
        expect(a).to eq(b)
      end

      it "exposes lapacke_potri" do
        a = NMatrix.new(3, [4, 0,-1,
                            0, 2, 1,
                            0, 0, 1], dtype: dtype)
        NMatrix::LAPACK::lapacke_potrf(:row, :upper, 3, a, 3)
        NMatrix::LAPACK::lapacke_potri(:row, :upper, 3, a, 3)
        b = NMatrix.new(3, [0.5, -0.5, 1,  0, 1.5, -2,  0, 0, 4], dtype: dtype)
        err = case dtype
                when :float32, :complex64
                  1e-6
                when :float64, :complex128
                  1e-14
              end
        expect(a).to be_within(err).of(b)
      end

      it "exposes lapacke_potrs with vector solution" do
        a = NMatrix.new(3, [4, 0,-1,
                            0, 2, 1,
                            0, 0, 1], dtype: dtype)
        b = NMatrix.new([3,1], [3,0,2], dtype: dtype)

        NMatrix::LAPACK::lapacke_potrf(:row, :upper, 3, a, 3)
        #ldb is different from CLAPACK versions
        NMatrix::LAPACK::lapacke_potrs(:row, :upper, 3, 1, a, 3, b, 1)

        x = NMatrix.new([3,1], [3.5, -5.5, 11], dtype: dtype)

        err = case dtype
                when :float32, :complex64
                  1e-5
                when :float64, :complex128
                  1e-14
              end

        expect(b).to be_within(err).of(x)
      end

      it "exposes lapacke_potrs with matrix solution" do
        a = NMatrix.new(3, [4, 0,-1,
                            0, 2, 1,
                            0, 0, 1], dtype: dtype)
        b = NMatrix.new([3,2], [3,4,
                                0,4,
                                2,0], dtype: dtype)

        NMatrix::LAPACK::lapacke_potrf(:row, :upper, 3, a, 3)
        #ldb is different from CLAPACK versions
        NMatrix::LAPACK::lapacke_potrs(:row, :upper, 3, 2, a, 3, b, 2)

        x = NMatrix.new([3,2], [3.5, 0,
                                -5.5, 4,
                                11, -4], dtype: dtype)

        err = case dtype
                when :float32, :complex64
                  1e-5
                when :float64, :complex128
                  1e-14
              end

        expect(b).to be_within(err).of(x)
      end

      it "calculates the singular value decomposition with lapacke_gesvd" do
        #example from Wikipedia
        m = 4
        n = 5
        mn_min = [m,n].min
        a = NMatrix.new([m,n],[1,0,0,0,2, 0,0,3,0,0, 0,0,0,0,0, 0,4,0,0,0], dtype: dtype)
        s = NMatrix.new([mn_min], 0, dtype: a.abs_dtype) #s is always real and always returned as float/double, never as complex
        u = NMatrix.new([m,m], 0, dtype: dtype)
        vt = NMatrix.new([n,n], 0, dtype: dtype)
        superb = NMatrix.new([mn_min-1], dtype: a.abs_dtype)

        NMatrix::LAPACK.lapacke_gesvd(:row, :a, :a, m, n, a, n, s, u, m, vt, n, superb)

        s_true = NMatrix.new([mn_min], [4,3,Math.sqrt(5),0], dtype: a.abs_dtype)
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

      it "calculates the singular value decomposition with lapacke_gesdd" do
        #example from Wikipedia
        m = 4
        n = 5
        mn_min = [m,n].min
        a = NMatrix.new([m,n],[1,0,0,0,2, 0,0,3,0,0, 0,0,0,0,0, 0,4,0,0,0], dtype: dtype)
        s = NMatrix.new([mn_min], 0, dtype: a.abs_dtype) #s is always real and always returned as float/double, never as complex
        u = NMatrix.new([m,m], 0, dtype: dtype)
        vt = NMatrix.new([n,n], 0, dtype: dtype)

        NMatrix::LAPACK.lapacke_gesdd(:row, :a, m, n, a, n, s, u, m, vt, n)

        s_true = NMatrix.new([mn_min], [4,3,Math.sqrt(5),0], dtype: a.abs_dtype)
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

      it "calculates eigenvalues and eigenvectors using lapacke_geev" do
        n = 3
        a = NMatrix.new([n,n], [-1,0,0, 0,1,-2, 0,1,-1], dtype: dtype)
        w = NMatrix.new([n], dtype: dtype)
        if a.complex_dtype? #for real dtypes, imaginary parts of eigenvalues are stored in separate vector
          wi = nil
        else
          wi = NMatrix.new([n], dtype: dtype)
        end
        vl = NMatrix.new([n,n], dtype: dtype)
        vr = NMatrix.new([n,n], dtype: dtype)

        NMatrix::LAPACK.lapacke_geev(:row, :t, :t, n, a, n, w, wi, vl, n, vr, n)

        if !a.complex_dtype?
          w = w + wi*Complex(0,1)
        end

        w_true = NMatrix.new([n], [Complex(0,1), -Complex(0,1), -1], dtype: NMatrix.upcast(dtype, :complex64))
        if a.complex_dtype?
          #For complex types the right/left eigenvectors are stored as columns
          #of vr/vl.
          vr_true = NMatrix.new([n,n],[0,0,1,
                                       2/Math.sqrt(6),2/Math.sqrt(6),0,
                                       Complex(1,-1)/Math.sqrt(6),Complex(1,1)/Math.sqrt(6),0], dtype: dtype)
          vl_true = NMatrix.new([n,n],[0,0,1,
                                       Complex(-1,1)/Math.sqrt(6),Complex(-1,-1)/Math.sqrt(6),0,
                                       2/Math.sqrt(6),2/Math.sqrt(6),0], dtype: dtype)
        else
          #For real types, the real part of the first and second eigenvectors is
          #stored in the first column, the imaginary part of the first (= the
          #negative of the imaginary part of the second) eigenvector is stored
          #in the second column, and the third eigenvector (purely real) is the
          #third column.
          vr_true = NMatrix.new([n,n],[0,0,1,
                                       2/Math.sqrt(6),0,0,
                                       1/Math.sqrt(6),-1/Math.sqrt(6),0], dtype: dtype)
          vl_true = NMatrix.new([n,n],[0,0,1,
                                       -1/Math.sqrt(6),1/Math.sqrt(6),0,
                                       2/Math.sqrt(6),0,0], dtype: dtype)
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
