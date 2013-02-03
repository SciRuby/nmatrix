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
# == nvector_spec.rb
#
# Basic tests for NVector.
#

require "./lib/nmatrix"

describe NVector do
  it "initializes" do
    v = NVector.new(5, 0, :float64)
    v.shape[0].should == 5
    v.shape[1].should == 1
  end

  it "permits setting and getting contents" do
    v = NVector.new(5, 0, :float64)
    v[0] = 1.555
    v[0].should == 1.555
  end
  
  it "transpose() changes raw and column stored structure" do
    v = NVector.new 5, :float64
    v = v.transpose
    v.shape[0].should == 1
    v.shape[1].should == 5
    v[0] = 1.555
    v[0].should == 1.555
  end

  it "transpose!() changes destructively its raw and column stored structure" do
    v = NVector.new 5, :float64
    v.transpose!
    v.shape[0].should == 1
    v.shape[1].should == 5
  end

  it "multiply() multiples itself by another NVector" do
    v1 = NVector.new 2, :float64
    v2 = NVector.new 2, :float64
    v1[0] = 1.5
    v1[1] = 2.3
    v2[0] = 1.3
    v2[1] = 2.5
    v1.multiply(v2).should == 12.0
  end

  it "multiply!() multiples destructively itself by another NVector" do
    v1 = NVector.new 2, :float64
    v2 = NVector.new 2, :float64
    v1[0] = 1.5
    v1[1] = 2.3
    v2[0] = 1.3
    v2[1] = 2.5
    v1.multiply!(v2)
    v1.should == 12.0
  end

  it "pretty_print() prints values to standard output with a pretty format" do
    v = NVector.new(5, 0)
    $stdout = StringIO.new
    v.pretty_print
    out = $stdout.string
    $stdout = STDOUT
    out.should == "0  0  0  0  0\n"
  end

  it "inspect() formats the output with inspected, namely human readable format" do
    v = NVector.new(5, 0)
    $stdout = StringIO.new
    p v
    out = $stdout.string
    $stdout = STDOUT
    out.should == "0  0  0  0  0\n"
  end
end
