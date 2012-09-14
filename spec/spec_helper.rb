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
# == spec_helper.rb
#
# Common data for testing. 
require "./lib/nmatrix"

MATRIX43A_ARRAY = [14.0, 9.0, 3.0, 2.0, 11.0, 15.0, 0.0, 12.0, 17.0, 5.0, 2.0, 3.0]
MATRIX32A_ARRAY = [12.0, 25.0, 9.0, 10.0, 8.0, 5.0]

COMPLEX_MATRIX43A_ARRAY = MATRIX43A_ARRAY.zip(MATRIX43A_ARRAY.reverse).collect { |ary| Complex(ary[0], ary[1]) }
COMPLEX_MATRIX32A_ARRAY = MATRIX32A_ARRAY.zip(MATRIX32A_ARRAY.reverse).collect { |ary| Complex(ary[0], -ary[1]) }

RATIONAL_MATRIX43A_ARRAY = MATRIX43A_ARRAY.collect { |x| x.to_r }
RATIONAL_MATRIX32A_ARRAY = MATRIX32A_ARRAY.collect { |x| x.to_r }

def create_matrix(stype)
  m = NMatrix.new(stype, [3,3], 0, :int32)

  m[0,0] = 0
  m[0,1] = 1
  m[0,2] = 2
  m[1,0] = 3
  m[1,1] = 4
  m[1,2] = 5
  m[2,0] = 6
  m[2,1] = 7
  m[2,2] = 8

  m
end
