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
# Specs for the shortcuts used in NMatrix and in NVector.
#

# Can we use require_relative here instead?
require File.join(File.dirname(__FILE__), "spec_helper.rb")

describe NMatrix do
    
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
    m = NMatrix.random(2)
    
    m.stype.should == :dense
    m.dtype.should == :float64
  end
  
  it "seq() creates a matrix of integers, sequentially" do
    m = NMatrix.seq(2) # 2x2 matrix.
    value = 0
    
    2.times do |i|
      2.times do |j|
        m[i, j].should == value
        value += 1
      end
    end
  end
  
  it "seq() only accepts an integer or a 2-element array as dimension" do
    expect { NMatrix.seq([1, 2, 3]) }.to raise_error
    expect { NMatrix.seq("not an array or integer") }.to raise_error
  end
  
  it "column() returns a NMatrix" do
    m = NMatrix.random(3)
    
    m.column(2).is_a?(NMatrix).should be_true
  end
  
  it "column() accepts a second parameter (only :copy or :reference)" do
    m = NMatrix.random(3)
    
    expect { m.column(1, :copy) }.to_not raise_error
    expect { m.column(1, :reference) }.to_not raise_error
    
    expect { m.column(1, :derp) }.to raise_error
  end
end

describe "NVector" do
      
  it "zeros() creates a vector of zeros" do
    v = NVector.zeros(4)
    
    4.times do |i|
      v[i].should == 0
    end
  end
  
  it "ones() creates a vector of ones" do
    v = NVector.ones(3)
    
    3.times do |i|
      v[i].should == 1
    end
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
  
  it "seq() only accepts integers as dimension" do
    expect { NVector.seq(3) }.to_not raise_error

    expect { NVector.seq([1, 3]) }.to raise_error
    expect { NVector.seq(:wtf) }.to raise_error
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
