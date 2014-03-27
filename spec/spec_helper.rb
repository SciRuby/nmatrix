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
# == spec_helper.rb
#
# Common data and helper functions for testing.

require "rspec/longrun"
#require "narray/narray"

require "./lib/nmatrix"
require "./lib/nmatrix/rspec"

MATRIX43A_ARRAY = [14.0, 9.0, 3.0, 2.0, 11.0, 15.0, 0.0, 12.0, 17.0, 5.0, 2.0, 3.0]
MATRIX32A_ARRAY = [12.0, 25.0, 9.0, 10.0, 8.0, 5.0]

COMPLEX_MATRIX43A_ARRAY = MATRIX43A_ARRAY.zip(MATRIX43A_ARRAY.reverse).collect { |ary| Complex(ary[0], ary[1]) }
COMPLEX_MATRIX32A_ARRAY = MATRIX32A_ARRAY.zip(MATRIX32A_ARRAY.reverse).collect { |ary| Complex(ary[0], -ary[1]) }

RATIONAL_MATRIX43A_ARRAY = MATRIX43A_ARRAY.collect { |x| x.to_r }
RATIONAL_MATRIX32A_ARRAY = MATRIX32A_ARRAY.collect { |x| x.to_r }

def create_matrix(stype) #:nodoc:
  m = NMatrix.new([3,3], 0, dtype: :int32, stype: stype, default: 0)

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

def create_rectangular_matrix(stype) #:nodoc:
  m = NMatrix.new([5,6], 0, dtype: :int32, stype: stype, default: 0)

  m[0,0] = 1
  m[0,1] = 2
  m[0,2] = 3
  m[0,3] = 4
  m[0,4] = 5
  m[0,5] = 0

  m[1,0] = 6
  m[1,1] = 7
  m[1,2] = 8
  m[1,3] = 9
  m[1,4] = 0
  m[1,5] = 10

  m[2,0] = 11
  m[2,1] = 12
  m[2,2] = 13
  m[2,3] = 0
  m[2,4] = 14
  m[2,5] = 15

  # skip row 3 -- all 0
  m[3,0] = m[3,1] = m[3,2] = m[3,3] = m[3,4] = m[3,5] = 0

  m[4,0] = 16
  m[4,1] = 0
  m[4,2] = 17
  m[4,3] = 18
  m[4,4] = 19
  m[4,5] = 20

  m
end

def create_vector(stype) #:nodoc:
  m = stype == :yale ? NVector.new(stype, 10, :int32) : NVector.new(stype, 10, 0, :int32)

  m[0] = 1
  m[1] = 2
  m[2] = 3
  m[3] = 4
  m[4] = 5
  m[5] = 6
  m[6] = 7
  m[7] = 8
  m[8] = 9
  m[9] = 10

  m
end

# Stupid but independent comparison for slice_spec
def nm_eql(n, m) #:nodoc:
  if n.shape != m.shape
    false
  else # NMatrix
    n.shape[0].times do |i|
      n.shape[1].times do |j|
        if n[i,j] != m[i,j]
          puts "n[#{i},#{j}] != m[#{i},#{j}] (#{n[i,j]} != #{m[i,j]})"
          return false
        end
      end
    end
  end
  true
end

