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
# == homogeneous_spec.rb
#
# Specs for the homogeneous transformation matrix methods.
#

require 'spec_helper'
require "./lib/nmatrix/homogeneous.rb"

require 'pry'

describe 'NMatrix' do
  context "#x_rotation" do
    it "should generate a matrix representing a rotation about the x axis" do
      x = NMatrix.x_rotation(Math::PI.quo(6))
      expect(x).to be_within(1e-8).of(NMatrix.new([4,4], [1.0, 0.0, 0.0, 0.0,
                                                      0.0, Math.cos(Math::PI.quo(6)), -0.5, 0.0,
                                                      0.0, 0.5, Math.cos(Math::PI.quo(6)), 0.0,
                                                      0.0, 0.0, 0.0, 1.0] ))
    end
  end


  context "#y_rotation" do
    it "should generate a matrix representing a rotation about the y axis" do
      y = NMatrix.y_rotation(Math::PI.quo(6))
      expect(y).to be_within(1e-8).of(NMatrix.new([4,4], [Math.cos(Math::PI.quo(6)), 0.0, 0.5, 0.0,
                                                      0.0, 1.0, 0.0, 0.0,
                                                     -0.5, 0.0, Math.cos(Math::PI.quo(6)), 0.0,
                                                      0.0, 0.0, 0.0, 1.0] ))
    end
  end

  context "#z_rotation" do
    it "should generate a matrix representing a rotation about the z axis" do
      z = NMatrix.z_rotation(Math::PI.quo(6))
      expect(z).to be_within(1e-8).of(NMatrix.new([4,4], [Math.cos(Math::PI.quo(6)), -0.5, 0.0, 0.0,
                                                      0.5, Math.cos(Math::PI.quo(6)), 0.0, 0.0,
                                                      0.0, 0.0, 1.0, 0.0,
                                                      0.0, 0.0, 0.0, 1.0] ))
    end
  end

  context "#translation" do
    it "should generate a translation matrix from an Array" do
      t = NMatrix.translation([4,5,6])
      expect(t).to be_within(1e-8).of(NMatrix.new([4,4], [1, 0, 0, 4,
                                                      0, 1, 0, 5,
                                                      0, 0, 1, 6,
                                                      0, 0, 0, 1] ))
    end

    it "should generate a translation matrix from x, y, and z values" do
      t = NMatrix.translation(4,5,6)
      expect(t).to be_within(1e-8).of(NMatrix.new([4,4], [1, 0, 0, 4,
                                                      0, 1, 0, 5,
                                                      0, 0, 1, 6,
                                                      0, 0, 0, 1] ))
    end

    it "should generate a translation matrix from an NMatrix with correctly inferred dtype" do
      t = NMatrix.translation(NMatrix.new([3,1], [4,5,6], dtype: :float64) )
      expect(t).to be_within(1e-8).of(NMatrix.new([4,4], [1, 0, 0, 4,
                                                      0, 1, 0, 5,
                                                      0, 0, 1, 6,
                                                      0, 0, 0, 1] ))
      expect(t.dtype).to be(:float64)
    end
  end
end