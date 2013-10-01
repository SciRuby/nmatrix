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
# SciRuby is Copyright (c) 2010 - 2013, Ruby Science Foundation
# NMatrix is Copyright (c) 2013, Ruby Science Foundation
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
# == enumerate.rb
#
# Enumeration methods for NMatrix
#++

class NMatrix
  include Enumerable

  ##
  # call-seq:
  #   each -> Enumerator
  #
  # Enumerate through the matrix. @see Enumerable#each
  #
  # For dense, this actually calls a specialized each iterator (in C). For yale and list, it relies upon
  # #each_with_indices (which is about as fast as reasonably possible for C code).
  def each &bl
    if self.stype == :dense
      self.__dense_each__(&bl)
    elsif block_given?
      self.each_with_indices(&bl)
    else # Handle case where no block is given
      Enumerator.new do |yielder|
        self.each_with_indices do |params|
          yielder.yield params
        end
      end
    end
  end

  #
  # call-seq:
  #     flat_map -> Enumerator
  #     flat_map { |elem| block } -> Array
  #
  # Maps using Enumerator (returns an Array or an Enumerator)
  alias_method :flat_map, :map

  ##
  # call-seq:
  #   map -> Enumerator
  #   map { |elem| block } -> NMatrix
  #
  # Returns an NMatrix if a block is given. For an Array, use #flat_map
  #
  # Note that #map will always return an :object matrix, because it has no way of knowing
  # how to handle operations on the different dtypes.
  #
  def map(&bl)
    return enum_for(:map) unless block_given?
    cp = self.cast(dtype: :object)
    cp.map! &bl
    cp
  end

  ##
  # call-seq:
  #   map! -> Enumerator
  #   map! { |elem| block } -> NMatrix
  #
  # Maps in place.
  # @see #map
  #
  def map!
    return enum_for(:map!) unless block_given?
    self.each_stored_with_indices do |e, *i|
      self[*i] = (yield e)
    end
    self
  end


  #
  # call-seq:
  #     each_rank() -> NMatrix
  #     each_rank() { |rank| block } -> NMatrix
  #     each_rank(dimen) -> Enumerator
  #     each_rank(dimen) { |rank| block } -> NMatrix
  #
  # Generic for @each_row, @each_col
  #
  # Iterate through each rank by reference.
  #
  # @param [Fixnum] dimen the rank being iterated over.
  #
  def each_rank(dimen=0, get_by=:reference)
    return enum_for(:each_rank, dimen, get_by) unless block_given?
    (0...self.shape[dimen]).each do |idx|
      yield self.rank(dimen, idx, get_by)
    end
    self
  end
  alias :each_along_dim :each_rank

  #
  # call-seq:
  #     each_row { |row| block } -> NMatrix
  #
  # Iterate through each row, referencing it as an NMatrix slice.
  def each_row(get_by=:reference)
    return enum_for(:each_row, get_by) unless block_given?
    (0...self.shape[0]).each do |i|
      yield self.row(i, get_by)
    end
    self
  end

  #
  # call-seq:
  #     each_column { |column| block } -> NMatrix
  #
  # Iterate through each column, referencing it as an NMatrix slice.
  def each_column(get_by=:reference)
    return enum_for(:each_column, get_by) unless block_given?
    (0...self.shape[1]).each do |j|
      yield self.column(j, get_by)
    end
    self
  end

  #
  # call-seq:
  #     each_layer -> { |column| block } -> ...
  #
  # Iterate through each layer, referencing it as an NMatrix slice.
  #
  # Note: If you have a 3-dimensional matrix, the first dimension contains rows,
  # the second contains columns, and the third contains layers.
  def each_layer(get_by=:reference)
    return enum_for(:each_layer, get_by) unless block_given?
    (0...self.shape[2]).each do |k|
      yield self.layer(k, get_by)
    end
    self
  end


  #
  # call-seq:
  #     each_stored_with_index -> Enumerator
  #
  # Allow iteration across a vector NMatrix's stored values. See also @each_stored_with_indices
  #
  def each_stored_with_index(&block)
    raise(NotImplementedError, "only works for dim 2 vectors") unless self.dim <= 2
    return enum_for(:each_stored_with_index) unless block_given?

    self.each_stored_with_indices do |v, i, j|
      if shape[0] == 1
        yield(v,j)
      elsif shape[1] == 1
        yield(v,i)
      else
        method_missing(:each_stored_with_index, &block)
      end
    end
    self
  end


  ##
  # call-seq:
  #   inject_rank() -> Enumerator
  #   inject_rank(dimen) -> Enumerator
  #   inject_rank(dimen, initial) -> Enumerator
  #   inject_rank(dimen, initial, dtype) -> Enumerator
  #   inject_rank() { |elem| block } -> NMatrix
  #   inject_rank(dimen) { |elem| block } -> NMatrix
  #   inject_rank(dimen, initial) { |elem| block } -> NMatrix
  #   inject_rank(dimen, initial, dtype) { |elem| block } -> NMatrix
  #
  # Reduces an NMatrix using a supplied block over a specified dimension.
  # The block should behave the same way as for Enumerable#reduce.
  #
  # @param [Integer] dimen the dimension being reduced
  # @param [Numeric] initial the initial value for the reduction
  #  (i.e. the usual parameter to Enumerable#reduce).  Supply nil or do not
  #  supply this argument to have it follow the usual Enumerable#reduce
  #  behavior of using the first element as the initial value.
  # @param [Symbol] dtype if non-nil/false, forces the accumulated result to have this dtype
  # @return [NMatrix] an NMatrix with the same number of dimensions as the
  #  input, but with the input dimension now having size 1.  Each element
  #  is the result of the reduction at that position along the specified
  #  dimension.
  #
  def inject_rank(dimen=0, initial=nil, dtype=nil)

    raise(RangeError, "requested dimension (#{dimen}) does not exist (shape: #{shape})") if dimen > self.dim

    return enum_for(:inject_rank, dimen, initial, dtype) unless block_given?

    new_shape = shape
    new_shape[dimen] = 1

    first_as_acc = false

    if initial then
      acc = NMatrix.new(new_shape, initial, :dtype => dtype || self.dtype)
    else
      each_rank(dimen) do |sub_mat|
        acc = (sub_mat.is_a?(NMatrix) and !dtype.nil? and dtype != self.dtype) ? sub_mat.cast(self.stype, dtype) : sub_mat
        break
      end
      first_as_acc = true
    end

    each_rank(dimen) do |sub_mat|
      if first_as_acc
        first_as_acc = false
        next
      end
      acc = yield(acc, sub_mat)
    end

    acc
  end

  alias :reduce_along_dim :inject_rank
  alias :inject_along_dim :inject_rank

end