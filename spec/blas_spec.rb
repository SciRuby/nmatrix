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
# == blas_spec.rb
#
# Tests for properly exposed BLAS functions.
#

# Can we use require_relative here instead?
require File.join(File.dirname(__FILE__), "spec_helper.rb")

describe NMatrix::BLAS do
  #after :each do
  #  GC.start
  #end

  [:rational32, :rational64, :rational128, :float32, :float64, :complex64, :complex128].each do |dtype|
    context dtype do
      # This is not the same as "exposes cblas trsm", which would be for a version defined in blas.rb (which
      # would greatly simplify the calling of cblas_trsm in terms of arguments, and which would be accessible
      # as NMatrix::BLAS::trsm)
      it "exposes unfriendly cblas_trsm" do
        a     = NMatrix.new(3, [4,-1.quo(2), -3.quo(4), -2, 2, -1.quo(4), -4, -2, -1.quo(2)], dtype: dtype)
        b     = NMatrix.new([3,1], [-1, 17, -9], dtype: dtype)
        NMatrix::BLAS::cblas_trsm(:row, :right, :lower, :transpose, :nonunit, 1, 3, 1.0, a, 3, b, 3)

        # These test results all come from actually running a matrix through BLAS. We use them to ensure that NMatrix's
        # version of these functions (for rationals) give similar results.

        expect(b[0]).to eq(-1.quo(4))
        expect(b[1]).to eq(33.quo(4))
        expect(b[2]).to eq(-13)

        NMatrix::BLAS::cblas_trsm(:row, :right, :upper, :transpose, :unit, 1, 3, 1.0, a, 3, b, 3)

        expect(b[0]).to eq(-15.quo(2))
        expect(b[1]).to eq(5)
        expect(b[2]).to eq(-13)
      end
    end
  end

  [:rational32,:rational64,:rational128].each do |dtype|
    context dtype do
      it "exposes cblas rot"
    end

    context dtype do
      it "exposes cblas rotg"
    end
  end

  [:float32, :float64, :complex64, :complex128, :object].each do |dtype|
    context dtype do

      it "exposes cblas rot" do
        x = NMatrix.new([5,1], [1,2,3,4,5], dtype: dtype)
        y = NMatrix.new([5,1], [-5,-4,-3,-2,-1], dtype: dtype)
        x, y = NMatrix::BLAS::rot(x, y, 1.quo(2), Math.sqrt(3).quo(2), -1)

        expect(x).to be_within(1e-4).of(
                   NMatrix.new([5,1], [-0.3660254037844386, -0.7320508075688772, -1.098076211353316, -1.4641016151377544, -1.8301270189221928], dtype: dtype)
                 )

        expect(y).to be_within(1e-4).of(
                   NMatrix.new([5,1], [-6.830127018922193, -5.464101615137754, -4.098076211353316, -2.732050807568877, -1.3660254037844386], dtype: dtype)
                 )
      end

    end
  end

  [:float32, :float64, :complex64, :complex128, :object].each do |dtype|
    context dtype do

      it "exposes cblas rotg" do
        pending("broken for :object") if dtype == :object

        ab = NMatrix.new([2,1], [6,-8], dtype: dtype)
        c,s = NMatrix::BLAS::rotg(ab)

        if [:float32, :float64].include?(dtype)
          expect(ab[0]).to be_within(1e-6).of(-10)
          expect(ab[1]).to be_within(1e-6).of(-5.quo(3))
          expect(c).to be_within(1e-6).of(-3.quo(5))
        else
          pending "need correct test cases"
          expect(ab[0]).to be_within(1e-6).of(10)
          expect(ab[1]).to be_within(1e-6).of(5.quo(3))
          expect(c).to be_within(1e-6).of(3.quo(5))
        end
        expect(s).to be_within(1e-6).of(4.quo(5))
      end

      # Note: this exposes gemm, not cblas_gemm (which is the unfriendly CBLAS no-error-checking version)
      it "exposes gemm" do
        n = NMatrix.new([4,3], [14.0,9.0,3.0, 2.0,11.0,15.0, 0.0,12.0,17.0, 5.0,2.0,3.0], dtype: dtype)
        m = NMatrix.new([3,2], [12.0,25.0, 9.0,10.0, 8.0,5.0], dtype: dtype)

        #c = NMatrix.new([4,2], dtype)
        r = NMatrix::BLAS.gemm(n, m) #, c)
        #c.should equal(r) # check that both are same memory address

        expect(r).to eq(NMatrix.new([4,2], [273,455,243,235,244,205,102,160], dtype: dtype))
      end


      it "exposes gemv" do
        a = NMatrix.new([4,3], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0], dtype: :float64)
        x = NMatrix.new([3,1], [2.0, 1.0, 0.0], dtype: :float64)

        NMatrix::BLAS.gemv(a, x)
      end

      it "exposes asum" do
        x = NMatrix.new([4,1], [1,2,3,4], dtype: :float64)
        expect(NMatrix::BLAS.asum(x)).to eq(10.0)
      end


      it "exposes nrm2" do
        x = NMatrix.new([4,1], [2,-4,3,5], dtype: :float64)
        expect(NMatrix::BLAS.nrm2(x, 1, 3)).to be_within(1e-10).of(5.385164807134504)
      end

    end
  end
end