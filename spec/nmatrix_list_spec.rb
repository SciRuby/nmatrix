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
# == nmatrix_list_spec.rb
#
# Basic tests for NMatrix's list-of-lists storage type.
#
require "./lib/nmatrix"

describe NMatrix do
  context :list do
    it "should compare with ==" do
      n = NMatrix.new(:list, [3,3,3], :int64)
      m = NMatrix.new(:list, [3,3,3], :int64)
      n.should == m
      n[0,0,0] = 5
      n.should_not == m
      n[0,0,1] = 52
      n[1,2,1] = -4

      m[0,0,0] = 5
      m[0,0,1] = 52
      m[1,2,1] = -4
      n.should == m
    end

    it "should handle missing default value" do
      NMatrix.new(:list, 3, :int8)[0,0].should    == 0
      NMatrix.new(:list, 4, :float64)[0,0].should == 0.0
    end

    it "should allow conversion to a Ruby Hash" do
      n = NMatrix.new(:list, 3, 1, :int64)
      n[0,1] = 50
      h = n.to_h
      h.size.should == 1
      h[0].size.should == 1
      h[0][1].should == 50
      h[0][2].should == 1
      h[1][0].should == 1
    end


    ##TODO: Make this test better. It's not nearly exhaustive enough as is.
    it "should handle recursive removal" do
      n = NMatrix.new(:list, [3,3,3], 0)
      n[0,0,0] = 2
      n[1,1,1] = 1
      n[1,0,0] = 3
      n[0,0,1] = 4

      n[0,0,0].should == 2
      n[1,1,1].should == 1
      n[1,0,0].should == 3
      n[0,0,1].should == 4

      n[1,1,1] = 0
      n[0,0,0].should == 2
      n[1,1,1].should == 0
      n[1,0,0].should == 3
      n[0,0,1].should == 4
    end

    it "should correctly insert a value between the middle and last entries of a three-element list" do
      n = NMatrix.new(:list, 5940, 0, :float64)
      n[0,0] = -7.0710685196786e-01
      n[330,0] = 7.0710685196786e-01
      n[1,0] = -7.0710685196786e-01
      n[2,0] = -7.0710685196786e-01
      n[0,0].should_not == 0
      n[330,0].should_not == 0
      n[2,0].should_not == 0
      n[1,0].should_not == 0
    end
  end
end
