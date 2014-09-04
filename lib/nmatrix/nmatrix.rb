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
# This file contains those core functionalities which can be
# implemented efficiently (or much more easily) in Ruby (e.g.,
# inspect, pretty_print, element-wise operations).
#++

require_relative './lapack.rb'
require_relative './yale_functions.rb'
require_relative './monkeys'


class NMatrix
  module IO
    extend AutoloadPatch

    module Matlab
      extend AutoloadPatch

      class << self
        def load_mat file_path
          NMatrix::IO::Matlab::Mat5Reader.new(File.open(file_path, "rb+")).to_ruby
        end
        alias :load :load_mat
      end
    end
  end

  class << self
    #
    # call-seq:
    #     load_matlab_file(path) -> Mat5Reader
    #
    # * *Arguments* :
    #   - +file_path+ -> The path to a version 5 .mat file.
    # * *Returns* :
    #   - A Mat5Reader object.
    #
    def load_matlab_file(file_path)
      NMatrix::IO::Mat5Reader.new(File.open(file_path, 'rb')).to_ruby
    end

    #
    # call-seq:
    #     load_pcd_file(path) -> PointCloudReader::MetaReader
    #
    # * *Arguments* :
    #   - +file_path+ -> The path to a PCL PCD file.
    # * *Returns* :
    #   - A PointCloudReader::MetaReader object with the matrix stored in its +matrix+ property
    #
    def load_pcd_file(file_path)
      NMatrix::IO::PointCloudReader::MetaReader.new(file_path)
    end

    #
    # Calculate the size of an NMatrix of a given shape.
    def size(shape)
      shape = [shape,shape] unless shape.is_a?(Array)
      (0...shape.size).inject(1) { |x,i| x * shape[i] }
    end
  end

  # TODO: Make this actually pretty.
  def pretty_print(q) #:nodoc:
    if self.shape.size > 1 and self.shape[1] > 100
      self.inspect.pretty_print(q)
    elsif self.dim > 3 || self.dim == 1
      self.to_a.pretty_print(q)
    else
      # iterate through the whole matrix and find the longest number
      longest = Array.new(self.shape[1], 0)
      self.each_column.with_index do |col, j|
        col.each do |elem|
          elem_len   = elem.inspect.size
          longest[j] = elem_len if longest[j] < elem_len
        end
      end

      if self.dim == 3
        q.group(0, "\n{ layers:", "}") do
          self.each_layer.with_index do |layer,k|
            q.group(0, "\n  [\n", "  ]\n") do
              layer.each_row.with_index do |row,i|
                q.group(0, "    [", "]\n") do
                  q.seplist(self[i,0...self.shape[1],k].to_flat_array, lambda { q.text ", "}, :each_with_index) { |v,j| q.text v.inspect.rjust(longest[j]) }
                end
              end
            end
          end
        end
      else # dim 2
        q.group(0, "\n[\n", "]") do
          self.each_row.with_index do |row,i|
            q.group(1, "  [", "]") do
              q.seplist(self.dim > 2 ? row.to_a[0] : row.to_a, lambda { q.text ", " }, :each_with_index) { |v,j| q.text v.inspect.rjust(longest[j]) }
            end
            q.breakable
          end
        end
      end
    end
  end
  #alias :pp :pretty_print


  #
  # call-seq:
  #     cast(stype, dtype, default) -> NMatrix
  #     cast(stype, dtype) -> NMatrix
  #     cast(stype) -> NMatrix
  #     cast(options) -> NMatrix
  #
  # This is a user-friendly helper for calling #cast_full. The easiest way to call this function is using an
  # options hash, e.g.,
  #
  #     n.cast(:stype => :yale, :dtype => :int64, :default => false)
  #
  # For list and yale, :default sets the "default value" or "init" of the matrix. List allows a bit more freedom
  # since non-zeros are permitted. For yale, unpredictable behavior may result if the value is not false, nil, or
  # some version of 0. Dense discards :default.
  #
  # dtype and stype are inferred from the matrix upon which #cast is called -- so you only really need to provide
  # one. You can actually call this function with no arguments, in which case it functions like #clone.
  #
  # If your dtype is :object and you are converting from :dense to a sparse type, it is recommended that you
  # provide a :default, as 0 may behave differently from its Float, Rational, or Complex equivalent. If no option
  # is given, Fixnum 0 will be used.
  def cast(*params)
    if (params.size > 0 && params[0].is_a?(Hash))
      opts = {
          :stype => self.stype,
          :dtype => self.dtype,
          :default => self.stype == :dense ? 0 : self.default_value
      }.merge(params[0])

      self.cast_full(opts[:stype], opts[:dtype], opts[:default])
    else
      params << self.stype if params.size == 0
      params << self.dtype if params.size == 1
      #HACK: the default value can cause an exception if dtype is not complex
      #and default_value is. (The ruby C code apparently won't convert these.)
      #Perhaps this should be fixed in the C code (in rubyval_to_cval).
      default_value = maybe_get_noncomplex_default_value(params[1])
      params << (self.stype == :dense ? 0 : default_value) if params.size == 2
      self.cast_full(*params)
    end

  end


  #
  # call-seq:
  #     rows -> Integer
  #
  # This shortcut use #shape to return the number of rows (the first dimension)
  # of the matrix.
  #
  def rows
    shape[0]
  end

  #
  # call-seq:
  #     cols -> Integer
  #
  # This shortcut use #shape to return the number of columns (the second
  # dimension) of the matrix.
  #
  def cols
    shape[1]
  end

  #
  # call-seq:
  #     to_hash -> Hash
  #
  # Create a Ruby Hash from an NMatrix.
  #
  def to_hash
    if stype == :yale
      h = {}
      each_stored_with_indices do |val,i,j|
        next if val == 0 # Don't bother storing the diagonal zero values -- only non-zeros.
        if h.has_key?(i)
          h[i][j] = val
        else
          h[i] = {j => val}
        end
      end
      h
    else # dense and list should use a C internal function.
      # FIXME: Write a C internal to_h function.
      m = stype == :dense ? self.cast(:list, self.dtype) : self
      m.__list_to_hash__
    end
  end
  alias :to_h :to_hash


  def inspect #:nodoc:
    original_inspect = super()
    original_inspect = original_inspect[0...original_inspect.size-1]
    original_inspect + " " + inspect_helper.join(" ") + ">"
  end

  def __yale_ary__to_s(sym) #:nodoc:
    ary = self.send("__yale_#{sym.to_s}__".to_sym)

    '[' + ary.collect { |a| a ? a : 'nil'}.join(',') + ']'
  end


  ##
  # call-seq:
  #   integer_dtype?() -> Boolean
  #
  # Checks if dtype is an integer type
  #
  def integer_dtype?
    [:byte, :int8, :int16, :int32, :int64].include?(self.dtype)
  end


  #
  # call-seq:
  #     to_f -> Float
  #
  # Converts an nmatrix with a single element (but any number of dimensions)
  #  to a float.
  #
  # Raises an IndexError if the matrix does not have just a single element.
  #
  def to_f
    raise IndexError, 'to_f only valid for matrices with a single element' unless shape.all? { |e| e == 1 }
    self[*Array.new(shape.size, 0)]
  end

  #
  # call-seq:
  #     to_flat_array -> Array
  #     to_flat_a -> Array
  #
  # Converts an NMatrix to a one-dimensional Ruby Array.
  #
  def to_flat_array
    ary = Array.new(self.size)
    self.each.with_index { |v,i| ary[i] = v }
    ary
  end
  alias :to_flat_a :to_flat_array

  #
  # call-seq:
  #     size -> Fixnum
  #
  # Returns the total size of the NMatrix based on its shape.
  #
  def size
    NMatrix.size(self.shape)
  end


  def to_s #:nodoc:
    self.to_flat_array.to_s
  end

  #
  # call-seq:
  #     nvector? -> true or false
  #
  # Shortcut function for determining whether the effective dimension is less than the dimension.
  # Useful when we take slices of n-dimensional matrices where n > 2.
  #
  def nvector?
    self.effective_dim < self.dim
  end

  #
  # call-seq:
  #     vector? -> true or false
  #
  # Shortcut function for determining whether the effective dimension is 1. See also #nvector?
  #
  def vector?
    self.effective_dim == 1
  end


  #
  # call-seq:
  #     to_a -> Array
  #
  # Converts an NMatrix to an array of arrays, or an NMatrix of effective dimension 1 to an array.
  #
  # Does not yet work for dimensions > 2
  def to_a(dimen=nil)
    if self.dim == 2

      return self.to_flat_a if self.shape[0] == 1

      ary = []
      begin
        self.each_row do |row|
          ary << row.to_flat_a
        end
      #rescue NotImplementedError # Oops. Try copying instead
      #  self.each_row(:copy) do |row|
      #    ary << row.to_a.flatten
      #  end
      end
      ary
    else
      to_a_rec(0)
    end
  end


  #
  # call-seq:
  #     rank(dimension, row_or_column_number) -> NMatrix
  #     rank(dimension, row_or_column_number, :reference) -> NMatrix reference slice
  #
  # Returns the rank (e.g., row, column, or layer) specified, using slicing by copy as default.
  #
  # See @row (dimension = 0), @column (dimension = 1)
  
  def rank(shape_idx, rank_idx, meth = :copy)
    if self.dim != 1 and shape_idx > (self.dim-1)
      raise(RangeError, "#rank call was out of bounds")
    end

    params = Array.new(self.dim)
    params.each.with_index do |v,d|
      params[d] = d == shape_idx ? rank_idx : 0...self.shape[d]
    end

    meth == :reference ? self[*params] : self.slice(*params)
  end

  #
  # call-seq:
  #     column(column_number) -> NMatrix
  #     column(column_number, get_by) -> NMatrix
  #
  # Returns the column specified. Uses slicing by copy as default.
  #
  # * *Arguments* :
  #   - +column_number+ -> Integer.
  #   - +get_by+ -> Type of slicing to use, +:copy+ or +:reference+.
  # * *Returns* :
  #   - A NMatrix representing the requested column as a column vector.
  #
  # Examples:
  #
  #   m = NMatrix.new(2, [1, 4, 9, 14], :int32) # =>  1   4
  #                                                   9  14
  #
  #   m.column(1) # =>   4
  #                     14
  #
  def column(column_number, get_by = :copy)
    rank(1, column_number, get_by)
  end

  alias :col :column

  #
  # call-seq:
  #     row(row_number) -> NMatrix
  #     row(row_number, get_by) -> NMatrix
  #
  # * *Arguments* :
  #   - +row_number+ -> Integer.
  #   - +get_by+ -> Type of slicing to use, +:copy+ or +:reference+.
  # * *Returns* :
  #   - An NMatrix representing the requested row as a row vector.
  #
  def row(row_number, get_by = :copy)
    rank(0, row_number, get_by)
  end


  #
  # call-seq:
  #     reshape(new_shape) -> NMatrix
  #
  # Clone a matrix, changing the shape in the process. Note that this function does not do a resize; the product of
  # the new and old shapes' components must be equal.
  #
  # * *Arguments* :
  #   - +new_shape+ -> Array of positive Fixnums.
  # * *Returns* :
  #   - A copy with a different shape.
  #
  def reshape new_shape,*shapes
    if new_shape.is_a?Fixnum
      newer_shape =  [new_shape]+shapes
    else  # new_shape is an Array
      newer_shape = new_shape
    end
    t = reshape_clone_structure(newer_shape)
    left_params  = [:*]*newer_shape.size
    puts(left_params)
    right_params = [:*]*self.shape.size
    t[*left_params] = self[*right_params]
    t
  end


  #
  # call-seq:
  #     reshape!(new_shape) -> NMatrix
  #     reshape! new_shape  -> NMatrix
  #
  # Reshapes the matrix (in-place) to the desired shape. Note that this function does not do a resize; the product of
  # the new and old shapes' components must be equal.
  #
  # * *Arguments* :
  #   - +new_shape+ -> Array of positive Fixnums.
  #
  def reshape! new_shape,*shapes
    if self.is_ref?
      raise(ArgumentError, "This operation cannot be performed on reference slices")
    else
      if new_shape.is_a?Fixnum
        shape =  [new_shape]+shapes
      else  # new_shape is an Array
        shape = new_shape
      end
      self.reshape_bang(shape)
    end
  end

  #
  # call-seq:
  #     transpose -> NMatrix
  #     transpose(permutation) -> NMatrix
  #
  # Clone a matrix, transposing it in the process. If the matrix is two-dimensional, the permutation is taken to be [1,0]
  # automatically (switch dimension 0 with dimension 1). If the matrix is n-dimensional, you must provide a permutation
  # of +0...n+.
  #
  # * *Arguments* :
  #   - +permutation+ -> Optional Array giving a permutation.
  # * *Returns* :
  #   - A copy of the matrix, but transposed.
  #
  def transpose(permute = nil)
    if self.dim == 1
      return self.clone
    elsif self.dim == 2
      new_shape = [self.shape[1], self.shape[0]]
    elsif permute.nil?
      raise(ArgumentError, "need permutation array of size #{self.dim}")
    elsif permute.sort.uniq != (0...self.dim).to_a
      raise(ArgumentError, "invalid permutation array")
    else
      # Figure out the new shape based on the permutation given as an argument.
      new_shape = permute.map { |p| self.shape[p] }
    end

    if self.dim > 2 # FIXME: For dense, several of these are basically equivalent to reshape.

      # Make the new data structure.
      t = self.reshape_clone_structure(new_shape)

      self.each_stored_with_indices do |v,*indices|
        p_indices = permute.map { |p| indices[p] }
        t[*p_indices] = v
      end
      t
    elsif self.list? # TODO: Need a C list transposition algorithm.
      # Make the new data structure.
      t = self.reshape_clone_structure(new_shape)

      self.each_column.with_index do |col,j|
        t[j,:*] = col.to_flat_array
      end
      t
    else
      # Call C versions of Yale and List transpose, which do their own copies
      self.clone_transpose
    end
  end


  #
  # call-seq:
  #     matrix1.concat(*m2) -> NMatrix
  #     matrix1.concat(*m2, rank) -> NMatrix
  #     matrix1.hconcat(*m2) -> NMatrix
  #     matrix1.vconcat(*m2) -> NMatrix
  #     matrix1.dconcat(*m3) -> NMatrix
  #
  # Joins two matrices together into a new larger matrix. Attempts to determine which direction to concatenate
  # on by looking for the first common element of the matrix +shape+ in reverse. In other words, concatenating two
  # columns together without supplying +rank+ will glue them into an n x 2 matrix.
  #
  # You can also use hconcat, vconcat, and dconcat for the first three ranks. concat performs an hconcat when no
  # rank argument is provided.
  #
  # The two matrices must have the same +dim+.
  #
  # * *Arguments* :
  #   - +matrices+ -> one or more matrices
  #   - +rank+ -> Fixnum (for rank); alternatively, may use :row, :column, or :layer for 0, 1, 2, respectively
  #
  def concat *matrices
    rank = nil
    rank = matrices.pop unless matrices.last.is_a?(NMatrix)

    # Find the first matching dimension and concatenate along that (unless rank is specified)
    if rank.nil?
      rank = self.dim-1
      self.shape.reverse_each.with_index do |s,i|
        matrices.each do |m|
          if m.shape[i] != s
            rank -= 1
            break
          end
        end
      end
    elsif rank.is_a?(Symbol) # Convert to numeric
      rank = {:row => 0, :column => 1, :col => 1, :lay => 2, :layer => 2}[rank]
    end

    # Need to figure out the new shape.
    new_shape = self.shape.dup
    new_shape[rank] = matrices.inject(self.shape[rank]) { |total,m| total + m.shape[rank] }

    # Now figure out the options for constructing the concatenated matrix.
    opts = {stype: self.stype, default: self.default_value, dtype: self.dtype}
    if self.yale?
      # We can generally predict the new capacity for Yale. Subtract out the number of rows
      # for each matrix being concatenated, and then add in the number of rows for the new
      # shape. That takes care of the diagonal. The rest of the capacity is represented by
      # the non-diagonal non-default values.
      new_cap = matrices.inject(self.capacity - self.shape[0]) do |total,m|
        total + m.capacity - m.shape[0]
      end - self.shape[0] + new_shape[0]
      opts = {capacity: new_cap}.merge(opts)
    end

    # Do the actual construction.
    n = NMatrix.new(new_shape, opts)

    # Figure out where to start and stop the concatenation. We'll use NMatrices instead of
    # Arrays because then we can do elementwise addition.
    ranges = self.shape.map.with_index { |s,i| 0...self.shape[i] }

    matrices.unshift(self)
    matrices.each do |m|
      n[*ranges] = m

      # move over by the requisite amount
      ranges[rank]  = (ranges[rank].first + m.shape[rank])...(ranges[rank].last + m.shape[rank])
    end

    n
  end

  def hconcat *matrices
    concat(*matrices, :column)
  end

  def vconcat *matrices
    concat(*matrices, :row)
  end

  def dconcat *matrices
    concat(*matrices, :layer)
  end


  #
  # call-seq:
  #     upper_triangle -> NMatrix
  #     upper_triangle(k) -> NMatrix
  #     triu -> NMatrix
  #     triu(k) -> NMatrix
  #
  # Returns the upper triangular portion of a matrix. This is analogous to the +triu+ method
  # in MATLAB.
  #
  # * *Arguments* :
  #   - +k+ -> Positive integer. How many extra diagonals to include in the upper triangular portion.
  #
  def upper_triangle(k = 0)
    raise(NotImplementedError, "only implemented for 2D matrices") if self.shape.size > 2

    t = self.clone_structure
    (0...self.shape[0]).each do |i|
      if i - k < 0
        t[i, :*] = self[i, :*]
      else
        t[i, 0...(i-k)]             = 0
        t[i, (i-k)...self.shape[1]] = self[i, (i-k)...self.shape[1]]
      end
    end
    t
  end
  alias :triu :upper_triangle


  #
  # call-seq:
  #     upper_triangle! -> NMatrix
  #     upper_triangle!(k) -> NMatrix
  #     triu! -> NMatrix
  #     triu!(k) -> NMatrix
  #
  # Deletes the lower triangular portion of the matrix (in-place) so only the upper portion remains.
  #
  # * *Arguments* :
  #   - +k+ -> Integer. How many extra diagonals to include in the deletion.
  #
  def upper_triangle!(k = 0)
    raise(NotImplementedError, "only implemented for 2D matrices") if self.shape.size > 2

    (0...self.shape[0]).each do |i|
      if i - k >= 0
        self[i, 0...(i-k)] = 0
      end
    end
    self
  end
  alias :triu! :upper_triangle!


  #
  # call-seq:
  #     lower_triangle -> NMatrix
  #     lower_triangle(k) -> NMatrix
  #     tril -> NMatrix
  #     tril(k) -> NMatrix
  #
  # Returns the lower triangular portion of a matrix. This is analogous to the +tril+ method
  # in MATLAB.
  #
  # * *Arguments* :
  #   - +k+ -> Integer. How many extra diagonals to include in the lower triangular portion.
  #
  def lower_triangle(k = 0)
    raise(NotImplementedError, "only implemented for 2D matrices") if self.shape.size > 2

    t = self.clone_structure
    (0...self.shape[0]).each do |i|
      if i + k >= shape[0]
        t[i, :*] = self[i, :*]
      else
        t[i, (i+k+1)...self.shape[1]] = 0
        t[i, 0..(i+k)] = self[i, 0..(i+k)]
      end
    end
    t
  end
  alias :tril :lower_triangle


  #
  # call-seq:
  #     lower_triangle! -> NMatrix
  #     lower_triangle!(k) -> NMatrix
  #     tril! -> NMatrix
  #     tril!(k) -> NMatrix
  #
  # Deletes the upper triangular portion of the matrix (in-place) so only the lower portion remains.
  #
  # * *Arguments* :
  #   - +k+ -> Integer. How many extra diagonals to include in the deletion.
  #
  def lower_triangle!(k = 0)
    raise(NotImplementedError, "only implemented for 2D matrices") if self.shape.size > 2

    (0...self.shape[0]).each do |i|
      if i + k < shape[0]
        self[i, (i+k+1)...self.shape[1]] = 0
      end
    end
    self
  end
  alias :tril! :lower_triangle!


  #
  # call-seq:
  #     layer(layer_number) -> NMatrix
  #     row(layer_number, get_by) -> NMatrix
  #
  # * *Arguments* :
  #   - +layer_number+ -> Integer.
  #   - +get_by+ -> Type of slicing to use, +:copy+ or +:reference+.
  # * *Returns* :
  #   - A NMatrix representing the requested layer as a layer vector.
  #
  def layer(layer_number, get_by = :copy)
    rank(2, layer_number, get_by)
  end



  #
  # call-seq:
  #     shuffle! -> ...
  #     shuffle!(random: rng) -> ...
  #
  # Re-arranges the contents of an NVector.
  #
  # TODO: Write more efficient version for Yale, list.
  # TODO: Generalize for more dimensions.
  def shuffle!(*args)
    method_missing(:shuffle!, *args) if self.effective_dim > 1
    ary = self.to_flat_a
    ary.shuffle!(*args)
    ary.each.with_index { |v,idx| self[idx] = v }
    self
  end


  #
  # call-seq:
  #     shuffle -> ...
  #     shuffle(rng) -> ...
  #
  # Re-arranges the contents of an NVector.
  #
  # TODO: Write more efficient version for Yale, list.
  # TODO: Generalize for more dimensions.
  def shuffle(*args)
    method_missing(:shuffle!, *args) if self.effective_dim > 1
    t = self.clone
    t.shuffle!(*args)
  end


  #
  # call-seq:
  #     sorted_indices -> Array
  #
  # Returns an array of the indices ordered by value sorted.
  #
  def sorted_indices
    return method_missing(:sorted_indices) unless vector?
    ary = self.to_flat_array
    ary.each_index.sort_by { |i| ary[i] }  # from: http://stackoverflow.com/a/17841159/170300
  end


  #
  # call-seq:
  #     binned_sorted_indices -> Array
  #
  # Returns an array of arrays of indices ordered by value sorted. Functions basically like +sorted_indices+, but
  # groups indices together for those values that are the same.
  #
  def binned_sorted_indices
    return method_missing(:sorted_indices) unless vector?
    ary = self.to_flat_array
    ary2 = []
    last_bin = ary.each_index.sort_by { |i| [ary[i]] }.inject([]) do |result, element|
      if result.empty? || ary[result[-1]] == ary[element]
        result << element
      else
        ary2 << result
        [element]
      end
    end
    ary2 << last_bin unless last_bin.empty?
    ary2
  end


  def method_missing name, *args, &block #:nodoc:
    if name.to_s =~ /^__list_elementwise_.*__$/
      raise NotImplementedError, "requested undefined list matrix element-wise operation"
    elsif name.to_s =~ /^__yale_scalar_.*__$/
      raise NotImplementedError, "requested undefined yale scalar element-wise operation"
    else
      super(name, *args, &block)
    end
  end


  def respond_to?(method) #:nodoc:
    if [:shuffle, :shuffle!, :each_with_index, :sorted_indices, :binned_sorted_indices, :nrm2, :asum].include?(method.intern) # vector-only methods
      return vector?
    elsif [:each_layer, :layer].include?(method.intern) # 3-or-more dimensions only
      return dim > 2
    else
      super(method)
    end
  end


  #
  # call-seq:
  #     inject -> symbol
  #
  # This overrides the inject function to use map_stored for yale matrices
  #
  def inject(sym)
    return super(sym) unless self.yale?
    return self.map_stored.inject(sym)
  end


  #
  # call-seq:
  #     clone_structure -> NMatrix
  #
  # This function is like clone, but it only copies the structure and the default value.
  # None of the other values are copied. It takes an optional capacity argument. This is
  # mostly only useful for dense, where you may not want to initialize; for other types,
  # you should probably use +zeros_like+.
  #
  def clone_structure(capacity = nil)
    opts = {stype: self.stype, default: self.default_value, dtype: self.dtype}
    opts = {capacity: capacity}.merge(opts) if self.yale?
    NMatrix.new(self.shape, opts)
  end

  # This is how you write an individual element-wise operation function:
  #def __list_elementwise_add__ rhs
  #  self.__list_map_merged_stored__(rhs){ |l,r| l+r }.cast(self.stype, NMatrix.upcast(self.dtype, rhs.dtype))
  #end
protected

  def inspect_helper #:nodoc:
    ary = []
    ary << "shape:[#{shape.join(',')}]" << "dtype:#{dtype}" << "stype:#{stype}"

    if stype == :yale
      ary <<	"capacity:#{capacity}"

      # These are enabled by the DEBUG_YALE compiler flag in extconf.rb.
      if respond_to?(:__yale_a__)
        ary << "ija:#{__yale_ary__to_s(:ija)}" << "ia:#{__yale_ary__to_s(:ia)}" <<
          "ja:#{__yale_ary__to_s(:ja)}" << "a:#{__yale_ary__to_s(:a)}" << "d:#{__yale_ary__to_s(:d)}" <<
          "lu:#{__yale_ary__to_s(:lu)}" << "yale_size:#{__yale_size__}"
      end

    end

    ary
  end


  # Clone the structure as needed for a reshape
  def reshape_clone_structure(new_shape) #:nodoc:
    raise(ArgumentError, "reshape cannot resize; size of new and old matrices must match") unless self.size == new_shape.inject(1) { |p,i| p *= i }

    opts = {stype: self.stype, default: self.default_value, dtype: self.dtype}
    if self.yale?
      # We can generally predict the change in capacity for Yale.
      opts = {capacity: self.capacity - self.shape[0] + new_shape[0]}.merge(opts)
    end
    NMatrix.new(new_shape, opts)
  end


  # Helper for converting a matrix into an array of arrays recursively
  def to_a_rec(dimen = 0) #:nodoc:
    return self.flat_map { |v| v } if dimen == self.dim-1

    ary = []
    self.each_rank(dimen) do |sect|
      ary << sect.to_a_rec(dimen+1)
    end
    ary
  end


  # NMatrix constructor helper for sparse matrices. Uses multi-slice-setting to initialize a matrix
  # with a given array of initial values.
  def __sparse_initial_set__(ary) #:nodoc:
    self[0...self.shape[0],0...self.shape[1]] = ary
  end


  # Function assumes the dimensions and such have already been tested.
  #
  # Called from inside NMatrix: nm_eqeq
  #
  # There are probably more efficient ways to do this, but currently it's unclear how.
  # We could use +each_row+, but for list matrices, it's still going to need to make a
  # reference to each of those rows, and that is going to require a seek.
  #
  # It might be more efficient to convert one sparse matrix type to the other with a
  # cast and then run the comparison. For now, let's assume that people aren't going
  # to be doing this very often, and we can optimize as needed.
  def dense_eql_sparse? m #:nodoc:
    m.each_with_indices do |v,*indices|
      return false if self[*indices] != v
    end

    return true
  end
  alias :sparse_eql_sparse? :dense_eql_sparse?


  #
  # See the note in #cast about why this is necessary.
  # If this is a non-dense matrix with a complex dtype and to_dtype is
  # non-complex, then this will convert the default value to noncomplex.
  # Returns 0 if dense.  Returns existing default_value if there isn't a
  # mismatch.
  #
  def maybe_get_noncomplex_default_value(to_dtype) #:nodoc:
    default_value = 0
    unless self.stype == :dense then
      if self.dtype.to_s.start_with?('complex') and not to_dtype.to_s.start_with?('complex') then
        default_value = self.default_value.real
      else
        default_value = self.default_value
      end
    end
    default_value
  end

end

require_relative './shortcuts.rb'
require_relative './math.rb'
require_relative './enumerate.rb'
require_relative './interpolation/one_dimensional.rb'