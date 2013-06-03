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
# == nmatrix.rb
#
# This file adds a few additional pieces of functionality (e.g., inspect,
# pretty_print).
#++

require_relative './shortcuts.rb'
require_relative './lapack.rb'

class NMatrix
  # Read and write extensions for NMatrix. These are only loaded when needed.
  #
  module IO
    module Matlab
      class << self
        def load_mat file_path
          NMatrix::IO::Matlab::Mat5Reader.new(File.open(file_path, "rb+")).to_ruby
        end
      end

      # FIXME: Remove autoloads
      autoload :MatReader, 'nmatrix/io/mat_reader'
      autoload :Mat5Reader, 'nmatrix/io/mat5_reader'
    end

    # FIXME: Remove autoloads
    autoload :Market, 'nmatrix/io/market'
  end

  # TODO: Make this actually pretty.
  def pretty_print(q = nil) #:nodoc:
    if dim != 2 || (dim == 2 && shape[1] > 10) # FIXME: Come up with a better way of restricting the display
      inspect
    else

      arr = (0...shape[0]).map do |i|
        ary = []
        (0...shape[1]).each do |j|
          o = begin
                self[i, j]
              rescue ArgumentError
                nil
              end
          ary << (o.nil? ? 'nil' : o)
        end
        ary.inspect
      end

      if q.nil?
        puts "[" + arr.join("\n") + "]"
      else
        q.group(1, "", "\n") do
          q.seplist(arr, lambda { q.text "  " }, :each)  { |v| q.text v.to_s }
        end
      end

    end
  end
  alias :pp :pretty_print

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
    else # dense and list should use a C internal functions.
      to_hash_c
    end
  end
  alias :to_h :to_hash

  #
  # call-seq:
  #     invert! -> NMatrix
  #
  # Use LAPACK to calculate the inverse of the matrix (in-place). Only works on
  # dense matrices.
  #
  # Note: If you don't have LAPACK, e.g., on a Mac, this may not work yet.
  #
  def invert!
    # Get the pivot array; factor the matrix
    pivot = self.getrf!

    # Now calculate the inverse using the pivot array
    NMatrix::LAPACK::clapack_getri(:row, self.shape[0], self, self.shape[0], pivot)

    self
  end

  #
  # call-seq:
  #     invert -> NMatrix
  #
  # Make a copy of the matrix, then invert it (requires LAPACK).
  #
  # * *Returns* :
  #   - A dense NMatrix.
  #
  def invert
    self.cast(:dense, self.dtype).invert!
  end
  alias :inverse :invert

  #
  # call-seq:
  #     getrf! -> NMatrix
  #
  # LU factorization of a general M-by-N matrix +A+ using partial pivoting with
  # row interchanges. Only works in dense matrices.
  #
  # * *Returns* :
  #   - The IPIV vector. The L and U matrices are stored in A.
  # * *Raises* :
  #   - +StorageTypeError+ -> ATLAS functions only work on dense matrices.
  #
  def getrf!
    raise(StorageTypeError, "ATLAS functions only work on dense matrices") unless self.stype == :dense
    NMatrix::LAPACK::clapack_getrf(:row, self.shape[0], self.shape[1], self, self.shape[0])
  end

  #
  # call-seq:
  #     det -> determinant
  #
  # Calculate the determinant by way of LU decomposition. This is accomplished
  # using clapack_getrf, and then by summing the diagonal elements. There is a
  # risk of underflow/overflow.
  #
  # There are probably also more efficient ways to calculate the determinant.
  # This method requires making a copy of the matrix, since clapack_getrf
  # modifies its input.
  #
  # For smaller matrices, you may be able to use +#det_exact+.
  #
  # This function is guaranteed to return the same type of data in the matrix
  # upon which it is called.
  # In other words, if you call it on a rational matrix, you'll get a rational
  # number back.
  #
  # Integer matrices are converted to rational matrices for the purposes of
  # performing the calculation, as xGETRF can't work on integer matrices.
  #
  # * *Returns* :
  #   - The determinant of the matrix. It's the same type as the matrix's dtype.
  # * *Raises* :
  #   - +NotImplementedError+ -> Must be used in 2D matrices.
  #
  def det
    raise(NotImplementedError, "determinant can be calculated only for 2D matrices") unless self.dim == 2

    # Cast to a dtype for which getrf is implemented
    new_dtype = [:byte,:int8,:int16,:int32,:int64].include?(self.dtype) ? :rational128 : self.dtype
    copy = self.cast(:dense, new_dtype)

    # Need to know the number of permutations. We'll add up the diagonals of
    # the factorized matrix.
    pivot = copy.getrf!

    prod = pivot.size % 2 == 1 ? -1 : 1 # odd permutations => negative
    [shape[0],shape[1]].min.times do |i|
      prod *= copy[i,i]
    end

    # Convert back to an integer if necessary
    new_dtype != self.dtype ? prod.to_i : prod
  end

  #
  # call-seq:
  #     complex_conjugate -> NMatrix
  #     complex_conjugate(new_stype) -> NMatrix
  #
  # Get the complex conjugate of this matrix. See also complex_conjugate! for
  # an in-place operation (provided the dtype is already +:complex64+ or
  # +:complex128+).
  #
  # Doesn't work on list matrices, but you can optionally pass in the stype you
  # want to cast to if you're dealing with a list matrix.
  #
  # * *Arguments* :
  #   - +new_stype+ -> stype for the new matrix.
  # * *Returns* :
  #   - If the original NMatrix isn't complex, the result is a +:complex128+ NMatrix. Otherwise, it's the original dtype.
  #
  def complex_conjugate(new_stype = self.stype)
    self.cast(new_stype, NMatrix::upcast(dtype, :complex64)).complex_conjugate!
  end

  #
  # call-seq:
  #     conjugate_transpose -> NMatrix
  #
  # Calculate the conjugate transpose of a matrix. If your dtype is already
  # complex, this should only require one copy (for the transpose).
  #
  # * *Returns* :
  #   - The conjugate transpose of the matrix as a copy.
  #
  def conjugate_transpose
    self.transpose.complex_conjugate!
  end

  #
  # call-seq:
  #     hermitian? -> Boolean
  #
  # A hermitian matrix is a complex square matrix that is equal to its
  # conjugate transpose. (http://en.wikipedia.org/wiki/Hermitian_matrix)
  #
  # * *Returns* :
  #   - True if +self+ is a hermitian matrix, nil otherwise.
  #
  def hermitian?
    return false if self.dim != 2 or self.shape[0] != self.shape[1]

    if [:complex64, :complex128].include?(self.dtype)
      # TODO: Write much faster Hermitian test in C
      self.eql?(self.conjugate_transpose)
    else
      symmetric?
    end
  end

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
  # Successively yields submatrices at each coordinate along a specified
  # dimension.  Each submatrix will have the same number of dimensions as
  # the matrix being iterated, but with the specified dimension's size 
  # equal to 1.
  #
  # @param [Integer] dim the dimension being iterated over.
  #
  def each_along_dim(dim=0) 
    return enum_for(:each_along_dim, dim) unless block_given?
    dims = shape
    shape.each_index { |i| dims[i] = 0...(shape[i]) unless i == dim }
    0.upto(shape[dim]-1) do |i|
      dims[dim] = i
      yield self[*dims]
    end
  end

  ##
  # Reduces an NMatrix using a supplied block over a specified dimension.
  # The block should behave the same way as for Enumerable#reduce.
  #
  # @param [Integer] dim the dimension being reduced
  # @param [Numeric] initial the initial value for the reduction 
  #  (i.e. the usual parameter to Enumerable#reduce).  Supply nil or do not
  #  supply this argument to have it follow the usual Enumerable#reduce 
  #  behavior of using the first element as the initial value.
  # @return [NMatrix] an NMatrix with the same number of dimensions as the
  #  input, but with the input dimension now having size 1.  Each element 
  #  is the result of the reduction at that position along the specified
  #  dimension.
  #
  def reduce_along_dim(dim=0, initial=nil)

    if dim > shape.size then
      raise ArgumentError, "Requested dimension does not exist.  Requested: #{dim}, shape: #{shape}"
    end

    return enum_for(:reduce_along_dim, dim, initial) unless block_given?

    new_shape = shape
    new_shape[dim] = 1

    first_as_acc = false

    if initial then
      acc = NMatrix.new(new_shape, initial)
    else
      each_along_dim(dim) do |sub_mat|
        acc = sub_mat
        break
      end
      first_as_acc = true
    end

    each_along_dim(dim) do |sub_mat|
      if first_as_acc then
        first_as_acc = false
        next
      end
      acc = (yield acc, sub_mat)
    end

    acc

  end

  alias_method :inject_along_dim, :reduce_along_dim

  ##
  # Calculates the mean along the specified dimension.
  #
  # @see #reduce_along_dim
  #
  def mean(dim=0)
    reduce_along_dim(dim, 0.0) do |mean, sub_mat|
      mean + sub_mat/shape[dim]
    end
  end

  ##
  # Calculates the sum along the specified dimension.
  #
  # @see #reduce_along_dim
  def sum(dim=0)
    reduce_along_dim(dim, 0.0) do |sum, sub_mat|
      sum + sub_mat
    end
  end


  ##
  # Calculates the minimum along the specified dimension.
  #
  # @see #reduce_along_dim
  #
  def min(dim=0)
    reduce_along_dim(dim, Float::MAX) do |min, sub_mat|
      min * (min <= sub_mat) + ((min)*0.0 + (min > sub_mat)) * sub_mat
    end
  end

  ##
  # Calculates the maximum along the specified dimension.
  #
  # @see #reduce_along_dim
  #
  def max(dim=0)
    reduce_along_dim(dim, -1.0*Float::MAX) do |max, sub_mat|
      max * (max >= sub_mat) + ((max)*0.0 + (max < sub_mat)) * sub_mat
    end
  end


  ##
  # Calculates the sample variance along the specified dimension.
  #
  # @see #reduce_along_dim
  #
  def variance(dim=0)
    m = mean(dim)
    reduce_along_dim(dim, 0.0) do |var, sub_mat|
      var + (m - sub_mat)*(m - sub_mat)/(shape[dim]-1)
    end
  end

  ##
  # Calculates the sample standard deviation along the specified dimension.
  #
  # @see #reduce_along_dim
  #
  def std(dim=0)
    variance(dim).map! { |e| Math.sqrt(e) }
  end

  ##
  # Converts an nmatrix with a single element (but any number of dimensions)
  #  to a float.
  #
  # Raises an IndexError if the matrix does not have just a single element.
  #
  def to_f
    raise IndexError, 'to_f only valid for matrices with a single element' unless shape.all? { |e| e == 1 }
    self[*Array.new(shape.size, 0)]
  end

  ##
  # @see Enumerable#map
  #
  def map(&bl)
    return enum_for(:map) unless block_given?
    cp = self.dup
    cp.map! &bl
    cp
  end

  ##
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


  class << self
    #
    # call-seq:
    #     load_file(path) -> Mat5Reader
    #
    # * *Arguments* :
    #   - +path+ -> The path to a version 5 .mat file.
    # * *Returns* :
    #   - A Mat5Reader object.
    #
    def load_file(file_path)
      NMatrix::IO::Mat5Reader.new(File.open(file_path, 'rb')).to_ruby
    end

    ##
    # Creates a new matrix of ones with the same dtype and shape as the
    # provided matrix.
    #
    # @param [NMatrix] nm the nmatrix whose dtype and shape will be used
    # @return [NMatrix] a new nmatrix filled with ones.
    #
    def ones_like(nm)
      NMatrix.ones(nm.shape, nm.dtype)
    end

    ##
    # Creates a new matrix of zeros with the same stype, dtype, and shape
    # as the provided matrix.
    #
    # @param [NMatrix] nm the nmatrix whose stype, dtype, and shape will be used
    # @return [NMatrix] a new nmatrix filled with zeros.
    #
    def zeros_like(nm)
      NMatrix.zeros(nm.stype, nm.shape, nm.dtype)
    end
  end

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
end
