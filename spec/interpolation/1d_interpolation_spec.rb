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
# == io_spec.rb
# 
# Tests for 1D interpolation. Only linear interpolation impelemented as of now.

require './lib/nmatrix'

describe NMatrix::Interpolation::OneDimensional do
  before :each do
    @x   = NMatrix.seq [10]
    @y   = @x.exp
    @nd  = NMatrix.new([10,3]).each_column { |c| c[0..9] = @y }
  end

  it "tests if OneDimensional accepts Array inputs" do
    f = NMatrix::Interpolation::OneDimensional.new([0,1,2,3,4,5,6,7,8,9], [1.0, 
      2.718281828459045, 7.38905609893065, 20.085536923187668, 54.598150033144236, 
      148.4131591025766, 403.4287934927351, 1096.6331584284585, 2980.9579870417283, 
      8103.083927575384], {kind: :linear, sorted: true})

    expect(f.interpolate(2.5)).to eq 13.737
  end

  it "tests for linear interpolation for 1-dimensional y values" do
    f = NMatrix::Interpolation::OneDimensional.new(@x, @y, {kind: :linear, 
      precision: 3})

    expect(f.interpolate(2.5))              .to eq 13.737

    expect(f.interpolate([2.5,6.7,0.3,8.6])).to eq [13.737, 888.672, 
      1.515, 6054.234]

    expect(f.interpolate(NMatrix.new([4,1], [2.5,6.7,0.3,8.6]))).to eq NMatrix.new(
      [4], [13.737, 888.672, 1.515, 6054.234])
  end

  it "tests linear interpolation for N-dimensional y values" do

    f = NMatrix::Interpolation::OneDimensional.new(@x,@nd, {kind: :linear, 
      sorted: true, precision: 3})

    expect(f.interpolate(2.5))              .to eq [13.737,13.737,13.737]
    
    expect(f.interpolate([2.5,6.7,0.3,8.6])).to eq \
      [ [13.737  , 13.737  , 13.737 ], 
        [888.672 , 888.672 , 888.672],
        [1.515   , 1.515   , 1.515  ],
        [6054.234, 6054.234, 6054.234 ]
      ]

    expect(f.interpolate(NMatrix.new([4,1], [2.5,6.7,0.3,8.6]))).to eq \
      NMatrix.new([4,3], 
      [ 13.737  , 13.737  , 13.737 , 
        888.672 , 888.672 , 888.672,
        1.515   , 1.515   , 1.515  ,
        6054.234, 6054.234, 6054.234 
      ]) 
    end

  it "tests linear interpolation for N-dimensional y on another axis" do
   f = NMatrix::Interpolation::OneDimensional.new(@x, @nd, {kind: :linear, axis: 1, 
    sorted: true, precision: 3})
   
   expect(f.interpolate(3.5))              .to eq 37.342

   expect(f.interpolate([2.5,6.7,0.3,8.6])).to eq [13.737, 888.672, 
    1.515, 6054.234]

   expect(f.interpolate(NMatrix.new([4,1], [2.5,6.7,0.3,8.6]))).to eq NMatrix.new(
    [4], [13.737, 888.672, 1.515, 6054.234])
  end
end
