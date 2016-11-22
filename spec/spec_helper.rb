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

require "../lib/nmatrix"
require "../lib/nmatrix/rspec"

ALL_DTYPES = [:byte,:int8,:int16,:int32,:int64, :float32,:float64, :object,
  :complex64, :complex128]
  
NON_INTEGER_DTYPES = [:float32, :float64, :complex64, :complex128,
  :object]

FLOAT_DTYPES = [:float32, :float64]
  
MATRIX43A_ARRAY = [14.0, 9.0, 3.0, 2.0, 11.0, 15.0, 0.0, 12.0, 17.0, 5.0, 2.0, 3.0]
MATRIX32A_ARRAY = [12.0, 25.0, 9.0, 10.0, 8.0, 5.0]

COMPLEX_MATRIX43A_ARRAY = MATRIX43A_ARRAY.zip(MATRIX43A_ARRAY.reverse).collect { |ary| Complex(ary[0], ary[1]) }
COMPLEX_MATRIX32A_ARRAY = MATRIX32A_ARRAY.zip(MATRIX32A_ARRAY.reverse).collect { |ary| Complex(ary[0], -ary[1]) }

#3x4 matrix used for testing various getrf and LU decomposition functions
GETRF_EXAMPLE_ARRAY = [-1,0,10,4,9,2,3,5,7,8,1,6]
GETRF_SOLUTION_ARRAY = [9.0, 2.0, 3.0, 5.0, 7.0/9, 58.0/9, -4.0/3, 19.0/9, -1.0/9, 1.0/29, 301.0/29, 130.0/29]

TAU_SOLUTION_ARRAY = [1.8571428571428572,1.9938461538461538, 0.0]

GEQRF_SOLUTION_ARRAY =[                -14.0,                -21.0, 14.000000000000002,
                         0.23076923076923078,  -175.00000000000003,  70.00000000000001,
                        -0.15384615384615385, 0.055555555555555546,              -35.0]

R_SOLUTION_ARRAY   = [-159.2388143638353, -41.00131005172065, -56.75123892439876,  -90.75048729628048, 
                                     0.0, 25.137473501580676,  2.073591725046292,   9.790607357775713, 
                                     0.0,                0.0, -20.83259700334131, -17.592414929551445]

Q_SOLUTION_ARRAY_1 = [-0.632455532033676, -0.5209522876558295, -0.3984263084135902,  -0.41214704991068,
                    -0.42783756578748666, -0.20837937347171134, 0.876505919951498, 0.07259770177184455,
                    -0.48364246567281094, 0.8265854747306287,-0.015758658987033422, -0.2873988222474053,
                    -0.42783756578748666,  0.044081783789183565, -0.26971376257215296, 0.8615487797670971]

Q_SOLUTION_ARRAY_2 = [-0.8571428571428572,   0.3942857142857143,  0.33142857142857146, 
                      -0.4285714285714286,  -0.9028571428571428, -0.03428571428571425, 
                       0.28571428571428575, -0.1714285714285714,   0.9428571428571428]

Q_SOLUTION_ARRAY_3 = [-0.7724247413634004, -0.026670393594597247, -0.6345460653374136, 
                      -0.5777485870360393,  -0.38541856437557026,  0.7194853024298236,
                      -0.26375478973384403,   0.9223563413020934, 0.28229805268947933]

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

def integer_dtype? dtype
  [:byte,:int8,:int16,:int32,:int64].include?(dtype)
end

# If a focus: true option is supplied to any test, running `rake spec focus=true`
# will run only the focused tests and nothing else.
if ENV["focus"] == "true"
  RSpec.configure do |c|
    c.filter_run :focus => true
  end
end

