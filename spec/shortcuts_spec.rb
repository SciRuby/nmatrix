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
# == shortcuts_spec.rb
#
# Basic tests for the shortcuts used in NMatrix and NVector.
#

# Can we use require_relative here instead?
require File.join(File.dirname(__FILE__), "spec_helper.rb")

describe "NMatrix" do
    
  it "zeros() creates a matrix of zeros" do
    m = NMatrix.zeros(3)
    n = NMatrix.new([3, 3], 0)
    
    m.should.eql? n
  end
  
  it "ones() creates a matrix of ones" do
    m = NMatrix.ones(3)
    n = NMatrix.new([3, 3], 1)
    
    m.should.eql? n    
  end
  
  it "eye() creates an identity matrix" do
    m = NMatrix.eye(3)
    identity3 = NMatrix.new([3, 3], [1, 0, 0, 0, 1, 0, 0, 0, 1])
    
    m.should.eql? identity3
  end
  
  it "random() creates a matrix of random numbers" do
    m = NMatrix.random(3)
    m.stype.should == :dense
    m.dtype.should == :float64
  end
  
  it "seq() creates a matrix of integers, sequentially" do
    m = NMatrix.seq(2)
    i = 0
    
    m.each do |elem|
      elem.should == i
      i += 1
    end
  end
  
end

describe "NVector" do
      
  it "zeros() creates a vector of zeros" do
    v = NVector.zeros(4)
    u = NVector.new(4, 0)
    
    v.should.eql? u
  end
  
  it "ones() creates a vector of ones" do
    v = NVector.ones(3)
    u = NVector.new(3, 1)
    
    v.should.eql? u
  end
  
  it "random() creates a vector of random numbers" do
    v = NVector.zeros(4)
    v.dtype.should == :float64
    v.stype.should == :dense
  end
  
  it "seq() creates a vector of integers, sequentially" do
    v = NVector.seq(7)
    i = 0
    
    v.each do |elem|
      elem.should == i
      i += 1
    end
  end
  
  it "linspace() creates a vector with n values equally spaced between a and b" do
    v = NVector.linspace(0, 2, 5)
    i = 0
    
    v.each do |elem|
      elem.should == i * 0.5
      i += 1
    end
  end
  
end

describe "Inline constructor" do
  
  it "creates a NMatrix with the given values" do
    m = NMatrix.new([2, 2], [1, 4, 6, 7])
    n = N[[1, 4], [6, 7]]
    
    m.should.eql? n 
  end
end