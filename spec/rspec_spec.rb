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
# == rspec_spec.rb
#
# A spec for testing monkey patches to RSpec for NMatrix.
#
require_relative 'spec_helper'

describe "RSpec" do
  it "should permit #be_within to be used on a dense NMatrix" do
    expect(NMatrix.new([4,1], 1.0, dtype: :complex128, stype: :dense) / 10000.0).to be_within(0.00000001).of(NMatrix.new([4,1], 0.0001, dtype: :float64, stype: :dense))
    expect(NMatrix.new([4,1], 1.0, dtype: :complex128, stype: :dense) / 10000.0).not_to be_within(0.00000001).of(NMatrix.new([4,1], 1.0, dtype: :float64, stype: :dense))
  end
end
