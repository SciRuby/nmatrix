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

require 'spec_helper'
require './lib/nmatrix'

describe NMatrix::Interpolation::OneDimensional , :focus => true do
  # TODO : Write doc stating that when user puts in NMatrix for the interpolation
  #   he receives an NMatrix as well

  it "tests for linear interpolation for 1-dimensional y values" do
    x = create_vector :dense
    y = x.exp

    f = NMatrix::Interpolation::OneDimensional.new(x, y, {kind: :linear})

    expect(f.interpolate(2.5))              .to eq 13.737

    expect(f.interpolate([2.5,6.7,0.3,8.6])).to eq [13.7373, 888.6718, 
      1.5155, 6054.2336]

    expect(f.interpolate(NMatrix.new([4,1], [2.5,6.7,0.3,8.6]))).to eq
      NMatrix.new([4,1], [13.7373, 888.6718, 1.5155, 6054.2336])
  end

  it "tests linear interpolation for N-dimensional y values" do
    x = create_vector :dense
    y = NMatrix.new [10,3]

    3.times { |col| y[0..9,col] = x.exp }

    f = NMatrix::Interpolation::OneDimensional.new(x,y, {kind: :linear})

    expect(f.interpolate(2.5))              .to eq [13.737,13.737,13.737]
    
    expect(f.interpolate([2.5,6.7,0.3,8.6])).to eq 
      [ [13.737  ,13.737  , 13.737 ],
        [888.671 ,888.671 ,888.671 ],
        [1.515   ,1.515   ,1.515   ],
        [6054.233,6054.233,6054.233] ]

    expect(f.interpolate(NMatrix.new([4,1], [2.5,6.7,0.3,8.6]))).to eq
      NMatrix.new([4,3], 
        [13.737  ,13.737  , 13.737,
         888.671 ,888.671 ,888.671,
         1.515   ,1.515   ,1.515  ,
         6054.233,6054.233,6054.233 
        ])

   f = NMatrix::Interpolation::OneDimensional.new(x, y, {kind: :linear, axis: 1})
   
   expect(f.interpolate(3.5))              .to eq 37.342
   expect(f.interpolate([2.5,6.7,0.3,8.6])).to eq [13.7373, 888.6718, 
      1.5155, 6054.2336]

   expect(f.interpolate(NMatrix.new([4,1], [2.5,6.7,0.3,8.6]))).to eq
    NMatrix.new([4,1], [13.7373, 888.6718, 1.5155, 6054.2336])
  end
end
