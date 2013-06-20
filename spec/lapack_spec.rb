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
# SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
# NMatrix is Copyright (c) 2012, Ruby Science Foundation
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
# == lapack_spec.rb
#
# Tests for properly exposed LAPACK functions.
#

# Can we use require_relative here instead?
require File.join(File.dirname(__FILE__), "spec_helper.rb")

describe NMatrix::LAPACK do
  # where integer math is allowed
  [:byte, :int8, :int16, :int32, :int64, :rational32, :rational64, :rational128, :float32, :float64, :complex64, :complex128].each do |dtype|
    context dtype do
      it "exposes clapack laswp" do
        a = NMatrix.new(:dense, [3,4], [1,2,3,4,5,6,7,8,9,10,11,12], dtype)
        NMatrix::LAPACK::clapack_laswp(3, a, 4, 0, 3, [2,1,3,0], 1)
        b = NMatrix.new(:dense, [3,4], [3,2,4,1,7,6,8,5,11,10,12,9], dtype)
        a.should == b
      end
    end
  end

  # where integer math is not allowed
  [:rational32, :rational64, :rational128, :float32, :float64, :complex64, :complex128].each do |dtype|
    context dtype do
      it "exposes clapack getrf" do
        a = NMatrix.new(:dense, 3, [4,9,2,3,5,7,8,1,6], dtype)
        NMatrix::LAPACK::clapack_getrf(:row, 3, 3, a, 3)

        # delta varies for different dtypes
        err = case dtype
                when :float32, :complex64
                  1e-6
                when :float64, :complex128
                  1e-15
                else
                  1e-64 # FIXME: should be 0, but be_within(0) does not work.
              end

        a[0,0].should == 9 # 8
        a[0,1].should be_within(err).of(2.quo(9)) # 1
        a[0,2].should be_within(err).of(4.quo(9)) # 6
        a[1,0].should == 5 # 1.quo(2)
        a[1,1].should be_within(err).of(53.quo(9)) # 17.quo(2)
        a[1,2].should be_within(err).of(7.quo(53)) # -1
        a[2,0].should == 1 # 3.quo(8)
        a[2,1].should be_within(err).of(52.quo(9))
        a[2,2].should be_within(err).of(360.quo(53))
        # FIXME: these are rounded, == won't work
        #a[2,1].should == 0.544118
        #a[2,2].should == 5.294118
      end

      it "exposes clapack potrf" do
        # first do upper
        begin
          a = NMatrix.new(:dense, 3, [25,15,-5, 0,18,0, 0,0,11], dtype)
          NMatrix::LAPACK::clapack_potrf(:row, :upper, 3, a, 3)
          b = NMatrix.new(:dense, 3, [5,3,-1, 0,3,1, 0,0,3], dtype)
          a.should == b
        rescue NotImplementedError => e
          pending e.to_s
        end

        # then do lower
        a = NMatrix.new(:dense, 3, [25,0,0, 15,18,0,-5,0,11], dtype)
        NMatrix::LAPACK::clapack_potrf(:row, :lower, 3, a, 3)
        b = NMatrix.new(:dense, 3, [5,0,0, 3,3,0, -1,1,3], dtype)
        a.should == b
      end

      # Together, these calls are basically xGESV from LAPACK: http://www.netlib.org/lapack/double/dgesv.f
      it "exposes clapack getrs" do
        a     = NMatrix.new(:dense, 3, [-2,4,-3,3,-2,1,0,-4,3], dtype)
        ipiv  = NMatrix::LAPACK::clapack_getrf(:row, 3, 3, a, 3)
        b     = NVector.new(3, [-1, 17, -9], dtype)

        NMatrix::LAPACK::clapack_getrs(:row, false, 3, 1, a, 3, ipiv, b, 3)

        b[0].should == 5
        b[1].should == -15.quo(2)
        b[2].should == -13
      end

      it "exposes clapack getri" do
        a = NMatrix.new(:dense, 3, [1,0,4,1,1,6,-3,0,-10], dtype)
        ipiv = NMatrix::LAPACK::clapack_getrf(:row, 3, 3, a, 3) # get pivot from getrf, use for getri

        begin
          NMatrix::LAPACK::clapack_getri(:row, 3, a, 3, ipiv)

          b = NMatrix.new(:dense, 3, [-5,0,-2,-4,1,-1,1.5,0,0.5], dtype)
          a.should == b
        rescue NotImplementedError => e
          pending e.to_s
        end
      end


    end
  end
end
