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
# == math_spec.rb
#
# Tests for non-BLAS and non-LAPACK math functions, or for simplified
# versions of unfriendly BLAS and LAPACK functions.
#

# Can we use require_relative here instead?
require File.join(File.dirname(__FILE__), "spec_helper.rb")

describe "math" do

  [:float32, :float64, :complex64, :complex128, :rational32, :rational64, :rational128].each do |dtype|
    context dtype do
      it "should correctly factorize a matrix" do
        m = NMatrix.new(:dense, 3, [4,9,2,3,5,7,8,1,6], dtype)
        a = m.factorize_lu
        a[0,0].should == 8
        a[0,1].should == 1
        a[0,2].should == 6
        a[1,0].should == 0.5
        a[1,1].should == 8.5
        a[1,2].should == -1
        a[2,0].should == 0.375
      end
    end

    context dtype do
      it "should correctly invert a matrix" do
        a = NMatrix.new(:dense, 3, [1,0,4,1,1,6,-3,0,-10], dtype)
        b = NMatrix.new(:dense, 3, [-5,0,-2,-4,1,-1,1.5,0,0.5], dtype)
        a.invert!
        a.should == b
      end
    end
  end

  # TODO: Get it working with ROBJ too
  [:byte,:int8,:int16,:int32,:int64,:float32,:float64,:rational64,:rational128].each do |left_dtype|
    [:byte,:int8,:int16,:int32,:int64,:float32,:float64,:rational64,:rational128].each do |right_dtype|

      # Won't work if they're both 1-byte, due to overflow.
      next if [:byte,:int8].include?(left_dtype) && [:byte,:int8].include?(right_dtype)

      # For now, don't bother testing int-int mult.
      #next if [:int8,:int16,:int32,:int64].include?(left_dtype) && [:int8,:int16,:int32,:int64].include?(right_dtype)
      it "dense handles #{left_dtype.to_s} dot #{right_dtype.to_s} matrix multiplication" do
        #STDERR.puts "dtype=#{dtype.to_s}"
        #STDERR.puts "2"

        nary = if left_dtype.to_s =~ /complex/
                 COMPLEX_MATRIX43A_ARRAY
               elsif left_dtype.to_s =~ /rational/
                 RATIONAL_MATRIX43A_ARRAY
               else
                 MATRIX43A_ARRAY
               end

        mary = if right_dtype.to_s =~ /complex/
                 COMPLEX_MATRIX32A_ARRAY
               elsif right_dtype.to_s =~ /rational/
                 RATIONAL_MATRIX32A_ARRAY
               else
                 MATRIX32A_ARRAY
               end

        n = NMatrix.new([4,3], nary, left_dtype)
        m = NMatrix.new([3,2], mary, right_dtype)

        m.shape[0].should == 3
        m.shape[1].should == 2
        m.dim.should == 2

        n.shape[0].should == 4
        n.shape[1].should == 3
        n.dim.should == 2

        n.shape[1].should == m.shape[0]

        r = n.dot m

        r[0,0].should == 273.0
        r[0,1].should == 455.0
        r[1,0].should == 243.0
        r[1,1].should == 235.0
        r[2,0].should == 244.0
        r[2,1].should == 205.0
        r[3,0].should == 102.0
        r[3,1].should == 160.0

        #r.dtype.should == :float64 unless left_dtype == :float32 && right_dtype == :float32
      end
    end
  end

  [:byte,:int8,:int16,:int32,:int64,:float32,:float64,:rational64,:rational128].each do |left_dtype|
    [:byte,:int8,:int16,:int32,:int64,:float32,:float64,:rational64,:rational128].each do |right_dtype|

      # Won't work if they're both 1-byte, due to overflow.
      next if [:byte,:int8].include?(left_dtype) && [:byte,:int8].include?(right_dtype)

      it "dense handles #{left_dtype.to_s} dot #{right_dtype.to_s} vector multiplication" do
        #STDERR.puts "dtype=#{dtype.to_s}"
        #STDERR.puts "2"
        n = NMatrix.new([4,3], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0], left_dtype)

        m = NVector.new(3, [2.0, 1.0, 0.0], right_dtype)

        m.shape[0].should == 3
        m.shape[1].should == 1

        n.shape[0].should == 4
        n.shape[1].should == 3
        n.dim.should == 2

        n.shape[1].should == m.shape[0]

        r = n.dot m
        # r.class.should == NVector

        r[0,0].should == 4
        r[1,0].should == 13
        r[2,0].should == 22
        r[3,0].should == 31

        #r.dtype.should == :float64 unless left_dtype == :float32 && right_dtype == :float32
      end
    end
  end

end
