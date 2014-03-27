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
# == nmatrix.rb
#
# This file defines the NVector class.
#++

# This is a specific type of NMatrix in which only one dimension is not 1.
# Although it is stored as a dim-2, n x 1, matrix, it acts as a dim-1 vector
# of size n. If the @orientation flag is set to :column, it is stored as n x 1
# instead of 1 x n.
class NVector < NMatrix
  #
  # call-seq:
  #     new(shape) -> NVector
  #     new(stype, shape) -> NVector
  #     new(shape, init) -> NVector
  #     new(:dense, shape, init) -> NVector
  #     new(:list, shape, init) -> NVector
  #     new(shape, init, dtype) -> NVector
  #     new(stype, shape, init, dtype) -> NVector
  #     new(stype, shape, dtype) -> NVector
  #
  # Creates a new NVector. See also NMatrix#initialize for a more detailed explanation of
  # the arguments.
  #
  # * *Arguments* :
  #   - +stype+ -> (optional) Storage type of the vector (:list, :dense, :yale). Defaults to :dense.
  #   - +shape+ -> Shape of the vector. Accepts [n,1], [1,n], or n, where n is a Fixnum.
  #   - +init+ -> (optional) Yale: capacity; List: default value (0); Dense: initial value or values (uninitialized by default).
  #   - +dtype+ -> (optional if +init+ provided) Data type stored in the vector. For :dense and :list, can be inferred from +init+.
  # * *Returns* :
  #   -
  #
  def initialize(*args)
    stype = args[0].is_a?(Symbol) ? args.shift : :dense
    shape = args[0].is_a?(Array) ? args.shift  : [1,args.shift]

    if shape.size != 2 || !shape.include?(1) || shape == [1,1]
      raise(ArgumentError, "shape must be a Fixnum or an Array of positive Fixnums where exactly one value is 1")
    end

    warn "NVector is deprecated"

    super(stype, shape, *args)
  end


  #
  # call-seq:
  #     orientation -> Symbol
  #
  # Orientation defaults to row (e.g., [1,3] is a row of length 3). It
  # may also be column, e.g., for [5,1].
  #
  def orientation
    shape[0] == 1 ? :row : :column
  end

  # Override NMatrix#each_row and #each_column
  def each_column(get_by=:reference, &block) #:nodoc:
    shape[0] == 1 ? self.each(&block) : (yield self)
  end
  def each_row(get_by=:reference, &block) #:nodoc:
    shape[0] == 1 ? (yield self) : self.each(&block)
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
    shape[0] == 1 ? super(0, i) : super(i, 0)
  end

  #
  # call-seq:
  #     vector[index] = obj -> obj
  #
  # Stores +value+ at position +index+.
  #
  def []=(i, val)
    shape[0] == 1 ? super(0, i, val) : super(i, 0, val)
  end

  #
  # call-seq:
  #     dim -> 1
  #
  # Returns the dimension of a vector, which is 1.
  #
  def dim; 1; end

  #
  # call-seq:
  #     size -> Fixnum
  #
  # Shorthand for the dominant shape component
  def size
    shape[0] > 1 ? shape[0] : shape[1]
  end

  #
  # call-seq:
  #     max -> Numeric
  #
  # Return the maximum element.
  def max
    max_so_far = self[0]
    self.each do |x|
      max_so_far = x if x > max_so_far
    end
    max_so_far
  end

  #
  # call-seq:
  #     min -> Numeric
  #
  # Return the minimum element.
  def min
    min_so_far = self[0]
    self.each do |x|
      min_so_far = x if x < min_so_far
    end
    min_so_far
  end


  # TODO: Make this actually pretty.
  def pretty_print(q = nil) #:nodoc:
    dimen = shape[0] == 1 ? 1 : 0

    arr = (0...shape[dimen]).inject(Array.new){ |a, i| a << self[i] }

    if q.nil?
      puts "[" + arr.join("\n") + "]"
    else
      q.group(1, "", "\n") do
        q.seplist(arr, lambda { q.text "  " }, :each)  { |v| q.text v.to_s }
      end
    end
  end

  def inspect #:nodoc:
    original_inspect = super()
    original_inspect = original_inspect[0...original_inspect.size-1]
    original_inspect += " orientation:#{self.orientation}>"
  end

end
