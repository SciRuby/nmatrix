#--
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
# == nmatrix.rb
#
# This file defines the NVector class.
#++

# This is a specific type of NMatrix in which only one dimension is not 1.
# Although it is stored as a dim-2, n x 1, matrix, it acts as a dim-1 vector
# of size n. If the @orientation flag is set to :row, it is stored as 1 x n
# instead of n x 1.
class NVector < NMatrix
  #
  # call-seq:
  #     new(length) -> NVector
  #     new(length, values) -> NVector
  #     new(length, values, dtype) -> NVector
  #
  # Creates a new NVector.
  #
  # * *Arguments* :
  #   - +length+ -> Size of the vector.
  #   - +values+ -> (optional) Initial values of the vector. Default is 0.
  #   - +dtype+ -> (optional) Default is a guess from the +values+.
  # * *Returns* :
  #   -
  #
  def initialize(length, *args)
    super(:dense, [length, 1], *args)
    orientation
  end

  #
  # call-seq:
  #     orientation -> Symbol
  #
  # Orientation defaults to column (e.g., [3,1] is a column of length 3). It
  # may also be row, e.g., for [1,5].
  #
  def orientation
    @orientation ||= :column
  end

  #
  # call-seq:
  #     transpose -> NVector
  #
  # Returns a transposed copy of the vector.
  #
  # * *Returns* :
  #   - NVector containing the transposed vector.
  #
  def transpose
    t = super()
    t.flip!
  end

  #
  # call-seq:
  #     transpose! -> NVector
  #
  # Transpose the vector in-place.
  #
  # * *Returns* :
  #   - NVector containing the transposed vector.
  #
  def transpose!
    super()
    self.flip!
  end

  #
  # call-seq:
  #     multiply(m) ->
  #
  # ...
  #
  # * *Arguments* :
  #   - ++ ->
  # * *Returns* :
  #   -
  #
  def multiply(m)
    t = super(m)
    t.flip!
  end

  #
  # call-seq:
  #     multiply!(m) ->
  #
  # ...
  #
  # * *Arguments* :
  #   - ++ ->
  # * *Returns* :
  #   -
  #
  def multiply!(m)
    super().flip!
  end

  #
  # call-seq:
  #     vector[index] -> element
  #     vector[range] -> NVector
  #
  # Retrieves an element or return a slice.
  #
  # Examples:
  #
  #   u = NVector.new(3, [10, 20, 30])
  #   u[0]              # => 10
  #   u[0] + u[1]       # => 30
  #   u[0 .. 1].shape   # => [2, 1]
  #
  def [](i)
    case @orientation
    when :column;	super(i, 0)
    when :row;	  super(0, i)
    end
  end

  #
  # call-seq:
  #     vector[index] = obj -> obj
  #
  # Stores +value+ at position +index+.
  #
  def []=(i, val)
    case @orientation
    when :column;	super(i, 0, val)
    when :row;	  super(0, i, val)
    end
  end

  #
  # call-seq:
  #     dim -> 1
  #
  # Returns the dimension of a vector, which is 1.
  #
  def dim; 1; end

  # shorthand for the dominant shape component
  def size
    shape[0] > 1 ? shape[0] : shape[1]
  end

  # TODO: Make this actually pretty.
  def pretty_print #:nodoc:
    dim = @orientation == :row ? 1 : 0

    puts (0...shape[dim]).inject(Array.new) { |a, i| a << self[i] }.join('  ')
  end

  protected
  def inspect_helper #:nodoc:
    super() << "orientation:#{self.orientation}"
  end

  #
  # call-seq:
  #     flip_orientation! -> NVector
  #
  # Flip the orientation of the vector.
  #
  # * *Returns* :
  #   - NVector with orientation changed.
  #
  def flip_orientation!
    returning(self) { @orientation = @orientation == :row ? :column : :row }
  end
  alias :flip! :flip_orientation!
end
