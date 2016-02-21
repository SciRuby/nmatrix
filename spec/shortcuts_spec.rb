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
# == shortcuts_spec.rb
#
# Specs for the shortcuts used in NMatrix and in NVector.
#

require 'spec_helper'
require 'pry'

describe NMatrix do
  it "zeros() creates a matrix of zeros" do
    m = NMatrix.zeros(3)
    n = NMatrix.new([3, 3], 0)

    expect(m).to eq n
  end

  it "ones() creates a matrix of ones" do
    m = NMatrix.ones(3)
    n = NMatrix.new([3, 3], 1)

    expect(m).to eq n
  end

  it "eye() creates an identity matrix" do
    m = NMatrix.eye(3)
    identity3 = NMatrix.new([3, 3], [1, 0, 0, 0, 1, 0, 0, 0, 1])

    expect(m).to eq identity3
  end

  it "hilbert() creates an hilbert matrix" do
    m = NMatrix.hilbert(8)
    expect(m[4, 0]).to be_within(0.000001).of(0.2)
    expect(m[4, 1]).to be_within(0.000001).of(0.16666666666666666)
    expect(m[4, 2]).to be_within(0.000001).of(0.14285714285714285)
    expect(m[4, 3]).to be_within(0.000001).of(0.125)

    m = NMatrix.hilbert(3)
    hilbert3 = NMatrix.new([3, 3], [1.0, 0.5, 0.3333333333333333,\
     0.5, 0.3333333333333333, 0.25, 0.3333333333333333, 0.25, 0.2])
    expect(m).to eq hilbert3
    0.upto(2) do |i|
      0.upto(2) do |j|
        expect(m[i, j]).to be_within(0.000001).of(hilbert3[i,j])
      end
    end
  end

  it "inv_hilbert() creates an inverse hilbert matrix" do
    m = NMatrix.inv_hilbert(6)
    inv_hilbert6 = [3360.0,  -88200.0,   564480.0, -1411200.0]
    expect(m[2,0]).to be_within(0.000001).of(inv_hilbert6[0])
    expect(m[2,1]).to be_within(0.000001).of(inv_hilbert6[1])
    expect(m[2,2]).to be_within(0.000001).of(inv_hilbert6[2])
    expect(m[2,3]).to be_within(0.000001).of(inv_hilbert6[3])

    m = NMatrix.inv_hilbert(3)
    inv_hilbert3 = NMatrix.new([3, 3], [  9.0,  -36.0,   30.0, -36.0,  192.0, -180.0, 30.0, -180.0,  180.0] )
    0.upto(2) do |i|
      0.upto(2) do |j|
        expect(m[i, j]).to be_within(0.000001).of(inv_hilbert3[i,j])
      end
    end
  end

  it "diag() creates a matrix with pre-supplied diagonal" do
    arr = [1,2,3,4]
    m = NMatrix.diag(arr)
    expect(m.is_a?(NMatrix)).to be_true
  end

  it "diagonals() contains the seeded values on the diagonal" do
    arr = [1,2,3,4]
    m = NMatrix.diagonals(arr)
    expect(m[0,0]).to eq(arr[0])
    expect(m[1,1]).to eq(arr[1])
    expect(m[2,2]).to eq(arr[2])
    expect(m[3,3]).to eq(arr[3])
  end

  ALL_DTYPES.each do |dtype|
    [:dense, :yale, :list].each do |stype|
      context "#block_diagonal #{dtype} #{stype}" do
        it "block_diagonal() creates a block-diagonal NMatrix" do
          a = NMatrix.new([2,2], [1,2,
                                  3,4])
          b = NMatrix.new([1,1], [123.0])
          c = NMatrix.new([3,3], [1,2,3,
                                  1,2,3,
                                  1,2,3])
          d = Array[ [1,1,1], [2,2,2], [3,3,3] ] 
          e = 12
          m = NMatrix.block_diagonal(a, b, c, d, e, dtype: dtype, stype: stype)
          expect(m).to eq(NMatrix.new([10,10], [1, 2,   0, 0, 0, 0, 0, 0, 0,  0,
                                                3, 4,   0, 0, 0, 0, 0, 0, 0,  0,
                                                0, 0, 123, 0, 0, 0, 0, 0, 0,  0,
                                                0, 0,   0, 1, 2, 3, 0, 0, 0,  0,
                                                0, 0,   0, 1, 2, 3, 0, 0, 0,  0,
                                                0, 0,   0, 1, 2, 3, 0, 0, 0,  0,
                                                0, 0,   0, 0, 0, 0, 1, 1, 1,  0,
                                                0, 0,   0, 0, 0, 0, 2, 2, 2,  0,
                                                0, 0,   0, 0, 0, 0, 3, 3, 3,  0,
                                                0, 0,   0, 0, 0, 0, 0, 0, 0, 12], dtype: dtype, stype: stype))
        end
      end
    end
  end

  context "::random" do
    it "creates a matrix of random numbers" do
      m = NMatrix.random(2)

      expect(m.stype).to eq(:dense)
      expect(m.dtype).to eq(:float64)
    end

    it "creates a complex matrix of random numbers" do
      m = NMatrix.random(2, :dtype => :complex128)
    end

    it "correctly accepts :scale parameter" do
      m = NMatrix.random([2,2], dtype: :byte, scale: 255)
      m.each do |v|
        expect(v).to be >= 0
        expect(v).to be < 255
      end
    end

    it "only accepts an integer or an array as dimension" do
      m = NMatrix.random([2, 2])

      expect(m.stype).to eq(:dense)
      expect(m.dtype).to eq(:float64)

      expect { NMatrix.random(2.0) }.to raise_error
      expect { NMatrix.random("not an array or integer") }.to raise_error
    end
  end
  
  context "::linspace" do
    it "creates a row vector when given only one shape parameter" do
      v = NMatrix.linspace(1, 10, 4)
      #Expect a row vector only
      expect(v.shape.length).to eq(1)
      
      ans = [1.0,4.0,7.0,10.0]

      expect(v[0]).to be_within(0.000001).of(ans[0])
      expect(v[1]).to be_within(0.000001).of(ans[1])
      expect(v[2]).to be_within(0.000001).of(ans[2])
      expect(v[3]).to be_within(0.000001).of(ans[3])
    end
    
    it "creates a matrix of input shape with each entry linearly spaced in row major order" do
      v = NMatrix.linspace(1, Math::PI, [2,2])
      expect(v.dtype).to eq(:float64)

      ans = [1.0, 1.7138642072677612, 2.4277284145355225, 3.1415927410125732]
      
      expect(v[0,0]).to be_within(0.000001).of(ans[0])
      expect(v[0,1]).to be_within(0.000001).of(ans[1])
      expect(v[1,0]).to be_within(0.000001).of(ans[2])
      expect(v[1,1]).to be_within(0.000001).of(ans[3])
    end
  end

  context "::logspace" do
    it "creates a logarithmically spaced vector" do
      v = NMatrix.logspace(1, 2, 10)
      
      expect(v.shape.length).to eq(1)
      
      #Unit test taken from Matlab R2015b output of logspace(1,2,10)
      ans = [10.0000, 12.9155, 16.6810, 21.5443, 27.8256, 35.9381, 46.4159, 59.9484, 77.4264, 100.0000]
      
      expect(v[0].round(4)).to be_within(0.000001).of(ans[0])
      expect(v[1].round(4)).to be_within(0.000001).of(ans[1])
      expect(v[2].round(4)).to be_within(0.000001).of(ans[2])
      expect(v[3].round(4)).to be_within(0.000001).of(ans[3])
      expect(v[4].round(4)).to be_within(0.000001).of(ans[4])
      expect(v[5].round(4)).to be_within(0.000001).of(ans[5])
      expect(v[6].round(4)).to be_within(0.000001).of(ans[6])
      expect(v[7].round(4)).to be_within(0.000001).of(ans[7])
      expect(v[8].round(4)).to be_within(0.000001).of(ans[8])
      expect(v[9].round(4)).to be_within(0.000001).of(ans[9])
    end
    
    it "creates a logarithmically spaced vector bounded by Math::PI if :pi is pre-supplied" do
      v = NMatrix.logspace(1, :pi, 7)
                 
      #Unit test taken from Matlab R2015b output of logspace(1,pi,10)
      ans = [10.0000, 8.2450, 6.7980, 5.6050, 4.6213, 3.8103, 3.1416]
      
      expect(v[0].round(4)).to be_within(0.000001).of(ans[0])
      expect(v[1].round(4)).to be_within(0.000001).of(ans[1])
      expect(v[2].round(4)).to be_within(0.000001).of(ans[2])
      expect(v[3].round(4)).to be_within(0.000001).of(ans[3])
      expect(v[4].round(4)).to be_within(0.000001).of(ans[4])
      expect(v[5].round(4)).to be_within(0.000001).of(ans[5])
      expect(v[6].round(4)).to be_within(0.000001).of(ans[6])
    end    

    it "creates a matrix of input shape with each entry logarithmically spaced in row major order" do
      v = NMatrix.logspace(1, 2, [3,2])
      
      ans = [10.0, 15.8489, 25.1189, 39.8107, 63.0957, 100.0]
      
      expect(v[0,0].round(4)).to be_within(0.000001).of(ans[0])
      expect(v[0,1].round(4)).to be_within(0.000001).of(ans[1])
      expect(v[1,0].round(4)).to be_within(0.000001).of(ans[2])
      expect(v[1,1].round(4)).to be_within(0.000001).of(ans[3])
      expect(v[2,0].round(4)).to be_within(0.000001).of(ans[4])
      expect(v[2,1].round(4)).to be_within(0.000001).of(ans[5])
    end
  end

  it "seq() creates a matrix of integers, sequentially" do
    m = NMatrix.seq(2) # 2x2 matrix.
    value = 0

    2.times do |i|
      2.times do |j|
        expect(m[i,j]).to eq(value)
        value += 1
      end
    end
  end

  it "indgen() creates a matrix of integers as well as seq()" do
    m = NMatrix.indgen(2) # 2x2 matrix.
    value = 0

    2.times do |i|
      2.times do |j|
        expect(m[i, j]).to eq(value)
        value += 1
      end
    end
  end

  it "findgen creates a matrix of floats, sequentially" do
    m = NMatrix.findgen(2) # 2x2 matrix.
    value = 0

    2.times do |i|
      2.times do |j|
        expect(m[i, j]/10).to be_within(Float::EPSILON).of(value.to_f/10)
        value += 1
      end
    end
  end

  it "bindgen() creates a matrix of bytes" do
    m = NMatrix.bindgen(2) # 2x2 matrix.
    value = 0

    2.times do |i|
      2.times do |j|
        expect(m[i, j]).to eq(value)
        value += 1
      end
    end
  end

  it "cindgen() creates a matrix of complexes" do
    m = NMatrix.cindgen(2) # 2x2 matrix.
    value = 0

    2.times do |i|
      2.times do |j|
        expect(m[i, j].real).to be_within(Float::EPSILON).of(value)
        expect(m[i, j].imag).to be_within(Float::EPSILON).of(0.0)
        value += 1
      end
    end
  end

  it "column() returns a NMatrix" do
    m = NMatrix.random(3)

    expect(m.column(2).is_a?(NMatrix)).to be_true
  end

  it "row() returns a NMatrix" do
    m = NMatrix.random(3)

    expect(m.row(2).is_a?(NMatrix)).to be_true
  end

  it "diagonals() creates an NMatrix" do
    arr = [1,2,3,4]
    m = NMatrix.diagonals(arr)
    expect(m.is_a?(NMatrix)).to be_true
  end

  it "diagonals() contains the seeded values on the diagonal" do
    arr = [1,2,3,4]
    m = NMatrix.diagonals(arr)
    expect(m[0,0]).to eq(arr[0])
    expect(m[1,1]).to eq(arr[1])
    expect(m[2,2]).to eq(arr[2])
    expect(m[3,3]).to eq(arr[3])
  end

  context "_like constructors" do
    before :each do
      @nm_1d = NMatrix[5.0,0.0,1.0,2.0,3.0]
      @nm_2d = NMatrix[[0.0,1.0],[2.0,3.0]]
    end

    it "should create an nmatrix of ones with dimensions and type the same as its argument" do
      expect(NMatrix.ones_like(@nm_1d)).to eq NMatrix[1.0, 1.0, 1.0, 1.0, 1.0]
      expect(NMatrix.ones_like(@nm_2d)).to eq NMatrix[[1.0, 1.0], [1.0, 1.0]]
    end

    it "should create an nmatrix of zeros with dimensions and type the same as its argument" do
      expect(NMatrix.zeros_like(@nm_1d)).to eq NMatrix[0.0, 0.0, 0.0, 0.0, 0.0]
      expect(NMatrix.zeros_like(@nm_2d)).to eq NMatrix[[0.0, 0.0], [0.0, 0.0]]
    end
  end

end

describe "NVector" do

  it "zeros() creates a vector of zeros" do
    v = NVector.zeros(4)

    4.times do |i|
      expect(v[i]).to eq(0)
    end
  end

  it "ones() creates a vector of ones" do
    v = NVector.ones(3)

    3.times do |i|
      expect(v[i]).to eq(1)
    end
  end

  it "random() creates a vector of random numbers" do
    v = NVector.random(4)
    expect(v.dtype).to eq(:float64)
    expect(v.stype).to eq(:dense)
  end

  it "seq() creates a vector of integers, sequentially" do
    v = NVector.seq(7)
    expect(v).to eq(NMatrix.new([7,1], [0, 1, 2, 3, 4, 5, 6]))
  end

  it "seq() only accepts integers as dimension" do
    expect { NVector.seq(3) }.to_not raise_error

    expect { NVector.seq([1, 3]) }.to raise_error
    expect { NVector.seq(:wtf) }.to raise_error
  end

  it "indgen() creates a vector of integers as well as seq()" do
    v = NVector.indgen(7)
    expect(v).to eq(NMatrix.new([7,1], [0, 1, 2, 3, 4, 5, 6]))
  end

  it "findgen creates a vector of floats, sequentially" do
    v = NVector.findgen(2)
    expect(v).to eq(NMatrix.new([2,1], [0.0, 1.0]))
  end

  it "bindgen() creates a vector of bytes, sequentially" do
    v = NVector.bindgen(4)
    expect(v).to eq(NMatrix.new([4,1], [0, 1, 2, 3], dtype: :byte))
  end

  it "cindgen() creates a vector of complexes, sequentially" do
    v = NVector.cindgen(2)
    expect(v).to eq(NMatrix.new([2,1], [Complex(0.0, 0.0), Complex(1.0, 0.0)], dtype: :complex64))
  end

  it "linspace() creates a vector with n values equally spaced between a and b" do
    v = NVector.linspace(0, 2, 5)
    expect(v).to eq(NMatrix.new([5,1], [0.0, 0.5, 1.0, 1.5, 2.0]))
  end

  it "logspace() creates a vector with n values logarithmically spaced between decades 10^a and 10^b" do
    v = NVector.logspace(0, 3, 4)
    expect(v).to eq(NMatrix.new([4,1], [1.0, 10.0, 100.0, 1000.0]))
  end
end

describe "Inline constructor" do

  it "creates a NMatrix with the given values" do
    m = NMatrix.new([2, 2], [1, 4, 6, 7])
    n = NMatrix[[1, 4], [6, 7]]

    expect(m).to eq n
  end
end
