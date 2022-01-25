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
# == shortcuts.rb
#
# These are shortcuts for NMatrix and NVector creation, contributed by Daniel
# Carrera (dcarrera@hush.com) and Carlos Agarie (carlos.agarie@gmail.com).
#
# TODO Make all the shortcuts available through modules, allowing someone
# to include them to make "MATLAB-like" scripts.
#
# There are some questions to be answered before this can be done, tho.
#++

class NMatrix

  # Methods for generating magic matrix.
  module MagicHelpers
    class << self
      def odd_magic(nm, shape)
        row = shape - 1
        col = shape / 2       
        nm[row,col] = 1
        (2..shape * shape).each do |index|
          if nm[(row + 1) % shape,(col + 1) % shape] == 0
            row = (row + 1) % shape
            col = (col + 1) % shape
          else
            row = (row - 1 + shape) % shape
          end
            nm[row,col] = index
        end
      end
    
      def doubly_even_magic(nm, shape)
        mini_square_num = shape / 4
        count = 1     
        inv_count = shape * shape
        shape.times do |row|
          shape.times do |col|
            if col >= mini_square_num and col < shape - mini_square_num
              if row >= mini_square_num and row < shape - mini_square_num
      	        nm[row,col] = count
              else 
                nm[row,col] = inv_count
              end
            elsif row < mini_square_num or row >= shape - mini_square_num
              nm[row,col] = count
            else
              nm[row,col] = inv_count
            end
            count += 1
            inv_count -= 1  
          end
        end
      end
    
      def singly_even_magic(nm, shape)
        half_shape = shape / 2
        complementary_pair = (shape - 2) / 4
        swap_col = NMatrix.new([shape])
        index = 0 
        mini_magic = NMatrix.new([half_shape,half_shape], 0, dtype: nm.dtype)
        odd_magic mini_magic, half_shape
        half_shape.times do |row|
          half_shape.times do |col|
            nm[row,col] = mini_magic[row,col]  	
            nm[row + half_shape,col + half_shape] = mini_magic[row,col] + half_shape * half_shape  
            nm[row,col + half_shape] = mini_magic[row,col] + 2 * half_shape * half_shape      
            nm[row + half_shape,col] = mini_magic[row,col] + 3 * half_shape * half_shape       
          end  
        end
  
        (1..complementary_pair).each do |complementary_entry|
          swap_col[index] = complementary_entry
          index += 1
        end
      
        (shape - complementary_pair + 2..shape).each do |center|
          swap_col[index] = center
          index += 1
        end 
      
        (1..half_shape).each do |row|
          (1..index).each do |col|
            temp = nm[row - 1,swap_col[col - 1] - 1]
            nm[row - 1,swap_col[col - 1] - 1] = nm[row + half_shape - 1,swap_col[col - 1] - 1]
            nm[row + half_shape - 1,swap_col[col - 1] - 1] = temp
          end
        end

        temp = nm[complementary_pair,0] 
        nm[complementary_pair,0] = nm[complementary_pair + half_shape,0] 
        nm[complementary_pair + half_shape,0] = temp

        temp = nm[complementary_pair + half_shape,complementary_pair]
        nm[complementary_pair + half_shape,complementary_pair] = nm[complementary_pair,complementary_pair] 
        nm[complementary_pair,complementary_pair] = temp
      end  
    end 
  end 
    
  # call-seq:
  #     m.dense? -> true or false
  #
  # Determine if +m+ is a dense matrix.
  def dense?; return stype == :dense; end

  # call-seq:
  #     m.yale? -> true or false
  #
  # Determine if +m+ is a Yale matrix.
  def yale?;  return stype == :yale; end

  # call-seq:
  #     m.list? -> true or false
  #
  # Determine if +m+ is a list-of-lists matrix.
  def list?;  return stype == :list; end

  class << self
    # call-seq:
    #     NMatrix[Numeric, ..., Numeric, dtype: Symbol] -> NMatrix
    #     NMatrix[Array, dtype: Symbol] -> NMatrix
    #
    # The default value for +dtype+ is guessed from the first parameter. For example:
    #   NMatrix[1.0, 2.0].dtype # => :float64
    #
    # But this is just a *guess*. If the other values can't be converted to
    # this dtype, a +TypeError+ will be raised.
    #
    # You can use the +N+ constant in this way:
    #   N = NMatrix
    #   N[1, 2, 3]
    #
    # NMatrix needs to have a succinct way to create a matrix by specifying the
    # components directly. This is very useful for using it as an advanced
    # calculator, it is useful for learning how to use, for testing language
    # features and for developing algorithms.
    #
    # The NMatrix::[] method provides a way to create a matrix in a way that is compact and
    # natural. The components are specified using Ruby array syntax. Optionally,
    # one can specify a dtype as the last parameter (default is :float64).
    #
    # Examples:
    #
    #   a = N[ 1,2,3,4 ]          =>  1  2  3  4
    #
    #   a = N[ 1,2,3,4, :int32 ]  =>  1  2  3  4
    #
    #   a = N[ [1,2,3], [3,4,5] ] =>  1.0  2.0  3.0
    #                                 3.0  4.0  5.0
    #
    #   a = N[ 3,6,9 ].transpose => 3
    #                               6
    #                               9
    #
    # SYNTAX COMPARISON:
    #
    #   MATLAB:  a = [ [1 2 3] ; [4 5 6] ]   or  [ 1 2 3 ; 4 5 6 ]
    #   IDL:   a = [ [1,2,3] , [4,5,6] ]
    #   NumPy:  a = array( [1,2,3], [4,5,6] )
    #
    #   SciRuby:      a = NMatrix[ [1,2,3], [4,5,6] ]
    #   Ruby array:   a =  [ [1,2,3], [4,5,6] ]
    def [](*params)
      options = params.last.is_a?(Hash) ? params.pop : {}

      # First find the dimensions of the array.
      i = 0
      shape = []
      row = params
      while row.is_a?(Array)
        shape[i] = row.length
        row = row[0]
        i += 1
      end

      # A row vector should be stored as 1xN, not N
      #shape.unshift(1) if shape.size == 1

      # Then flatten the array.
      NMatrix.new(shape, params.flatten, **options)
    end

    #
    # call-seq:
    #    zeros(shape) -> NMatrix
    #    zeros(shape, dtype: dtype) -> NMatrix
    #    zeros(shape, dtype: dtype, stype: stype) -> NMatrix
    #
    # Creates a new matrix of zeros with the dimensions supplied as
    # parameters.
    #
    # * *Arguments* :
    #   - +shape+ -> Array (or integer for square matrix) specifying the dimensions.
    #   - +dtype+ -> (optional) Default is +:float64+
    #   - +stype+ -> (optional) Default is +:dense+.
    # * *Returns* :
    #   - NMatrix filled with zeros.
    #
    # Examples:
    #
    #   NMatrix.zeros(2) # =>  0.0   0.0
    #                          0.0   0.0
    #
    #   NMatrix.zeros([2, 3], dtype: :int32) # =>  0  0  0
    #                                              0  0  0
    #
    #   NMatrix.zeros([1, 5], dtype: :int32) # =>  0  0  0  0  0
    #
    def zeros(shape, opts = {})
      NMatrix.new(shape, 0, {:dtype => :float64}.merge(opts))
    end
    alias :zeroes :zeros

    #
    # call-seq:
    #     ones(shape) -> NMatrix
    #     ones(shape, dtype: dtype, stype: stype) -> NMatrix
    #
    # Creates a matrix filled with ones.
    #
    # * *Arguments* :
    #   - +shape+ -> Array (or integer for square matrix) specifying the shape.
    #   - +opts+ -> (optional) Hash of options from NMatrix#initialize
    # * *Returns* :
    #   - NMatrix filled with ones.
    #
    # Examples:
    #
    #   NMatrix.ones([1, 3]) # =>  1.0   1.0   1.0
    #
    #   NMatrix.ones([2, 3], dtype: :int32) # =>  1  1  1
    #                                             1  1  1
    #
    def ones(shape, opts={})
      NMatrix.new(shape, 1, {:dtype => :float64, :default => 1}.merge(opts))
    end

    # call-seq:
    #   ones_like(nm) -> NMatrix
    #
    # Creates a new matrix of ones with the same dtype and shape as the
    # provided matrix.
    #
    # @param [NMatrix] nm the nmatrix whose dtype and shape will be used
    # @return [NMatrix] a new nmatrix filled with ones.
    #
    def ones_like(nm)
      NMatrix.ones(nm.shape, dtype: nm.dtype, stype: nm.stype, capacity: nm.capacity, default: 1)
    end

    # call-seq:
    #   zeros_like(nm) -> NMatrix
    #
    # Creates a new matrix of zeros with the same stype, dtype, and shape
    # as the provided matrix.
    #
    # @param [NMatrix] nm the nmatrix whose stype, dtype, and shape will be used
    # @return [NMatrix] a new nmatrix filled with zeros.
    #
    def zeros_like(nm)
      NMatrix.zeros(nm.shape, dtype: nm.dtype, stype: nm.stype, capacity: nm.capacity, default: 0)
    end

    #
    # call-seq:
    #     eye(shape) -> NMatrix
    #     eye(shape, dtype: dtype) -> NMatrix
    #     eye(shape, stype: stype, dtype: dtype) -> NMatrix
    #
    # Creates an identity matrix (square matrix rank 2).
    #
    # * *Arguments* :
    #   - +size+ -> Array (or integer for square matrix) specifying the dimensions.
    #   - +dtype+ -> (optional) Default is +:float64+
    #   - +stype+ -> (optional) Default is +:dense+.
    # * *Returns* :
    #   - An identity matrix.
    #
    # Examples:
    #
    #    NMatrix.eye(3) # =>   1.0   0.0   0.0
    #                          0.0   1.0   0.0
    #                          0.0   0.0   1.0
    #
    #    NMatrix.eye(3, dtype: :int32) # =>   1   0   0
    #                                         0   1   0
    #                                         0   0   1
    #
    #    NMatrix.eye(2, dtype: :int32, stype: :yale) # =>   1   0
    #                                                       0   1
    #
    def eye(shape, opts={})
      # Fill the diagonal with 1's.
      m = NMatrix.zeros(shape, {:dtype => :float64}.merge(opts))
      (0...m.shape[0]).each do |i|
        m[i, i] = 1
      end

      m
    end
    alias :identity :eye

    #
    # call-seq:
    #     hilbert(shape) -> NMatrix
    #     hilbert(shape, dtype: dtype) -> NMatrix
    #     hilbert(shape, stype: stype, dtype: dtype) -> NMatrix
    #
    # Creates an hilbert matrix (square matrix).
    #
    # * *Arguments* :
    #   - +size+ -> integer ( for square matrix) specifying the dimensions.
    #   - +dtype+ -> (optional) Default is +:float64+
    #   - +stype+ -> (optional) Default is +:dense+.
    # * *Returns* :
    #   - A hilbert matrix.
    #
    # Examples:
    #
    #    NMatrix.hilbert(3) # =>  1.0     0.5      0.3333333333333333
    #            0.5                         0.3333333333333333    0.25
    #            0.3333333333333333          0.25                  0.2
    #
    def hilbert(shape, opts={})
      m = NMatrix.new([shape,shape], {:dtype => :float64}.merge(opts))
      0.upto(shape - 1) do |i|
        0.upto(i) do |j|
          m[i,j] = 1.0 / (j + i + 1)
          m[j,i] = m[i,j] if i != j
        end
      end
      m
    end

    #
    # call-seq:
    #     inv_hilbert(shape) -> NMatrix
    #     inv_hilbert(shape, dtype: dtype) -> NMatrix
    #     inv_hilbert(shape, stype: stype, dtype: dtype) -> NMatrix
    #
    # Creates an inverse hilbert matrix (square matrix rank 2).
    #
    # * *Arguments* :
    #   - +size+ -> Array (or integer for square matrix) specifying the dimensions.
    #   - +dtype+ -> (optional) Default is +:float64+
    #   - +stype+ -> (optional) Default is +:dense+.
    # * *Returns* :
    #   - A hilbert matrix.
    #
    # Examples:
    #    NMatrix.inv_hilbert(3) # =>   9.0,  -36.0,   30.0
    #                          -36.0,  192.0, -180.0
    #                          30.0, -180.0,  180.0
    #
    #
    def inv_hilbert(shape, opts={})
      opts = {:dtype => :float64}.merge(opts)
      m = NMatrix.new([shape,shape],opts)
      combination = NMatrix.new([2*shape,2*shape],opts)
      #combinations refers to the combination of n things taken k at a time
      0.upto(2*shape-1) do |i|
        0.upto(i) do |j|
          if j != 0 and j != i
            combination[i,j] = combination[i-1,j] + combination[i-1,j-1]
          else
            combination[i,j] = 1
          end
        end
      end

      0.upto(shape-1) do |i|
        0.upto(i) do |j|
          m[i,j] = combination[shape + j,shape - i - 1] * ((i+j)+1) * \
          combination[shape + i,shape - j - 1] * (-1) ** ((i+j)) * \
          combination[(i+j),i] * combination[(i+j),i]
          m[j,i] = m[i,j] if i != j
        end
      end
      m
    end

    #
    # call-seq:
    #     diagonals(array) -> NMatrix
    #     diagonals(array, dtype: dtype, stype: stype) -> NMatrix
    #
    # Creates a matrix filled with specified diagonals.
    #
    # * *Arguments* :
    #   - +entries+ -> Array containing input values for diagonal matrix
    #   - +options+ -> (optional) Hash with options for NMatrix#initialize
    # * *Returns* :
    #   - NMatrix filled with specified diagonal values.
    #
    # Examples:
    #
    #   NMatrix.diagonal([1.0,2,3,4]) # => 1.0 0.0 0.0 0.0
    #                                      0.0 2.0 0.0 0.0
    #                                      0.0 0.0 3.0 0.0
    #                                      0.0 0.0 0.0 4.0
    #
    #   NMatrix.diagonal([1,2,3,4], dtype: :int32) # => 1 0 0 0
    #                                                   0 2 0 0
    #                                                   0 0 3 0
    #                                                   0 0 0 4
    #
    #
    def diagonal(entries, opts={})
      m = NMatrix.zeros(entries.size,
                        {:dtype => guess_dtype(entries[0]), :capacity => entries.size + 1}.merge(opts)
                       )
      entries.each_with_index do |n, i|
        m[i,i] = n
      end
      m
    end
    alias :diag :diagonal
    alias :diagonals :diagonal

    # Generate a block-diagonal NMatrix from the supplied 2D square matrices.
    #
    # * *Arguments*
    #   - +*params+ -> An array that collects all arguments passed to the method. The method
    #                  can receive any number of arguments. Optionally, the last entry of +params+ is 
    #                  a hash of options from NMatrix#initialize. All other entries of +params+ are 
    #                  the blocks of the desired block-diagonal matrix. Each such matrix block can be 
    #                  supplied as a square 2D NMatrix object, or alternatively as an array of arrays 
    #                  (with dimensions corresponding to a square matrix), or alternatively as a number.
    # * *Returns*
    #   - NMatrix of block-diagonal form filled with specified matrices
    #     as the blocks along the diagonal.
    #
    # * *Example* 
    #
    #  a = NMatrix.new([2,2], [1,2,3,4])
    #  b = NMatrix.new([1,1], [123], dtype: :float64)
    #  c = Array.new(2) { [[10,10], [10,10]] }
    #  d = Array[[1,2,3], [4,5,6], [7,8,9]]
    #  m = NMatrix.block_diagonal(a, b, *c, d, 10.0, 11, dtype: :int64, stype: :yale)
    #        => 
    #        [
    #          [1, 2,   0,  0,  0,  0,  0, 0, 0, 0,  0,  0]
    #          [3, 4,   0,  0,  0,  0,  0, 0, 0, 0,  0,  0]
    #          [0, 0, 123,  0,  0,  0,  0, 0, 0, 0,  0,  0]
    #          [0, 0,   0, 10, 10,  0,  0, 0, 0, 0,  0,  0]
    #          [0, 0,   0, 10, 10,  0,  0, 0, 0, 0,  0,  0]
    #          [0, 0,   0,  0,  0, 10, 10, 0, 0, 0,  0,  0]
    #          [0, 0,   0,  0,  0, 10, 10, 0, 0, 0,  0,  0]
    #          [0, 0,   0,  0,  0,  0,  0, 1, 2, 3,  0,  0]
    #          [0, 0,   0,  0,  0,  0,  0, 4, 5, 6,  0,  0]
    #          [0, 0,   0,  0,  0,  0,  0, 7, 8, 9,  0,  0]
    #          [0, 0,   0,  0,  0,  0,  0, 0, 0, 0, 10,  0]
    #          [0, 0,   0,  0,  0,  0,  0, 0, 0, 0,  0, 11]
    #        ]
    #
    def block_diagonal(*params)
      options = params.last.is_a?(Hash) ? params.pop : {}

      params.each_index do |i|
        params[i] = params[i].to_nm if params[i].is_a?(Array) # Convert Array to NMatrix
        params[i] = NMatrix.new([1,1], [params[i]]) if params[i].is_a?(Numeric) # Convert number to NMatrix
      end

      block_sizes = [] #holds the size of each matrix block
      params.each do |b|
        unless b.is_a?(NMatrix)
          raise(ArgumentError, "Only NMatrix or appropriate Array objects or single numbers allowed")
        end
        raise(ArgumentError, "Only 2D matrices or 2D arrays allowed") unless b.shape.size == 2
        raise(ArgumentError, "Only square-shaped blocks allowed") unless b.shape[0] == b.shape[1]
        block_sizes << b.shape[0]
      end

      block_diag_mat = NMatrix.zeros(block_sizes.inject(0,:+), options)
      (0...params.length).each do |n|
        # First determine the size and position of the n'th block in the block-diagonal matrix
        block_size = block_sizes[n]
        block_pos = block_sizes[0...n].inject(0,:+)
        # populate the n'th block in the block-diagonal matrix
        (0...block_size).each do |i|
          (0...block_size).each do |j|
            block_diag_mat[block_pos+i,block_pos+j] = params[n][i,j]
          end
        end
      end

      return block_diag_mat
    end
    alias :block_diag :block_diagonal

    #
    # call-seq:
    #     random(shape) -> NMatrix
    #
    # Creates a +:dense+ NMatrix with random numbers between 0 and 1 generated
    # by +Random::rand+. The parameter is the dimension of the matrix.
    #
    # If you use an integer dtype, make sure to specify :scale as a parameter, or you'll
    # only get a matrix of 0s.
    #
    # * *Arguments* :
    #   - +shape+ -> Array (or integer for square matrix) specifying the dimensions.
    # * *Returns* :
    #   - NMatrix filled with random values.
    #
    # Examples:
    #
    #   NMatrix.random([2, 2]) # => 0.4859439730644226   0.1783195585012436
    #                               0.23193766176700592  0.4503345191478729
    #
    #   NMatrix.random([2, 2], :dtype => :byte, :scale => 255) # => [ [252, 108] [44, 12] ]
    #
    def random(shape, opts={})
      scale = opts.delete(:scale) || 1.0

      if opts[:seed].nil?
        rng = Random.new
      else
        rng = Random.new(opts[:seed])
      end
      

      random_values = []


      # Construct the values of the final matrix based on the dimension.
      if opts[:dtype] == :complex64 || opts[:dtype] == :complex128
        NMatrix.size(shape).times { |i| random_values << Complex(rng.rand(scale), rng.rand(scale)) }
      else
        NMatrix.size(shape).times { |i| random_values << rng.rand(scale) }
      end

      NMatrix.new(shape, random_values, {:dtype => :float64, :stype => :dense}.merge(opts))
    end
    alias :rand :random
    
    # 
    #  call-seq:
    #    magic(shape) -> NMatrix
    #    magic(shape, dtype: dtype) -> NMatrix
    #
    #  The parameter is the dimension of the matrix.
    #
    #  Creates a +:dense+ NMatrix with the following properties:
    #    - An arrangement of the numbers from 1 to n^2 (n-squared) in the matrix, with each number occurring exactly once.
    #    - The sum of the entries of any row, any column, or any main diagonal is the same.
    #    - This sum must be n(n^2+1)/2.
    #   
    #  See: http://www.mathworks.com/help/matlab/ref/magic.html
    # 
    #  * *Arguments* :
    #   - +shape+ -> Array (or integer for square matrix) specifying the dimensions.
    #   - +dtype+ -> (optional) Default is +:float64+
    #  * *Returns* :
    #   - NMatrix with the above given properties.
    #
    #  Examples:
    #       
    #    NMatrix.magic(3) # => [  [4.0, 9.0, 2.0]   [3.0, 5.0, 7.0]   [8.0, 1.0, 6.0] ]
    #    
    #    NMatrix.magic(4, dtype :int32) # => [  [ 1, 15, 14,  4]
    #                                           [12,  6,  7,  9]
    #                                           [ 8, 10, 11,  5]
    #                                           [13,  3,  2, 16] ]
    #                             
    #    NMatrix.magic(6,dtype: :int64) # => [  [31,  9,  2, 22, 27, 20]
    #                                           [ 3, 32,  7, 21, 23, 25]
    #                                           [35,  1,  6, 26, 19, 24]
    #                                           [ 4, 36, 29, 13, 18, 11]
    #                                           [30,  5, 34, 12, 14, 16]
    #                                           [ 8, 28, 33, 17, 10, 15] ]
    #
    def magic(shape, opts={})
      raise(ArgumentError, "shape of two is not allowed") if shape == 2
      nm = NMatrix.new([shape,shape], 0, {:dtype => :float64}.merge(opts))
      if shape % 2 != 0
        MagicHelpers.odd_magic nm, shape
      elsif shape % 4 == 0
        MagicHelpers.doubly_even_magic nm, shape
      else   
        MagicHelpers.singly_even_magic nm, shape
      end
      nm
    end
    
    #
    # call-seq:
    #     linspace(base, limit) -> 1x100 NMatrix
    #     linspace(base, limit, *shape) -> NMatrix
    #
    # Returns an NMatrix with +[shape[0] x shape[1] x .. x shape[dim-1]]+ values of dtype +:float64+ equally spaced from
    # +base+ to +limit+, inclusive.
    #
    # See: http://www.mathworks.com/help/matlab/ref/linspace.html
    #
    # * *Arguments* :
    #   - +base+ -> The first value in the sequence.
    #   - +limit+ -> The last value in the sequence.
    #   - +shape+ -> Desired output shape. Default returns a 1x100 row vector.
    # * *Returns* :
    #   - NMatrix with +:float64+ values.
    #
    # Examples :-
    #
    #   NMatrix.linspace(1,Math::PI, 6)
    #     =>[1.0,
    #        1.4283185005187988,
    #        1.8566370010375977,
    #        2.2849555015563965,
    #        2.7132740020751953,
    #        3.1415927410125732
    #       ]
    #
    #   NMatrix.linspace(1,10, [3,2])
    #     =>[
    #         [              1.0, 2.799999952316284]
    #         [4.599999904632568, 6.400000095367432]
    #         [8.199999809265137,              10.0]
    #       ]
    #
    def linspace(base, limit, shape = [100])
      
      # Convert shape to array format 
      shape = [shape] if shape.is_a? Integer 
      
      #Calculate number of elements 
      count = shape.inject(:*)
            
      # Linear spacing between elements calculated in step
      #   step = limit - base / (count - 1)
      #   [Result Sequence] = [0->N sequence] * step + [Base]
      step = (limit - base) * (1.0 / (count - 1))
      result = NMatrix.seq(shape, {:dtype => :float64}) * step
      result += NMatrix.new(shape, base)
      result
    end

    # call-seq:
    #     logspace(base, limit) -> 1x50 NMatrix with exponent_base = 10 
    #     logspace(base, limit, shape , exponent_base:) -> NMatrix
    #     logspace(base, :pi, n) -> 1xn NMatrix with interval [10 ^ base, Math::PI]
    #
    # Returns an NMatrix with +[shape[0] x shape[1] x .. x shape[dim-1]]+ values of dtype +:float64+ logarithmically spaced from
    # +exponent_base ^ base+ to +exponent_base ^ limit+, inclusive.
    #
    # See: http://www.mathworks.com/help/matlab/ref/logspace.html
    #
    # * *Arguments* :
    #   - +base+ -> exponent_base ** base is the first value in the sequence
    #   - +limit+ -> exponent_base ** limit is the last value in the sequence.
    #   - +shape+ -> Desired output shape. Default returns a 1x50 row vector.
    # * *Returns* :
    #   - NMatrix with +:float64+ values.
    #
    # Examples :-
    #
    #   NMatrix.logspace(1,:pi,7)
    #     =>[
    #         10.0000, 
    #         8.2450, 
    #         6.7980, 
    #         5.6050, 
    #         4.6213, 
    #         3.8103, 
    #         3.1416
    #       ]
    #
    #   NMatrix.logspace(1,2,[3,2])
    #     =>[
    #         [10.0, 15.8489]
    #         [25.1189, 39.8107]
    #         [63.0957, 100.0]
    #       ]
    #
    def logspace(base, limit, shape = [50], exponent_base: 10)

      #Calculate limit for [10 ^ base ... Math::PI] if limit = :pi
      limit = Math.log(Math::PI, exponent_base = 10) if limit == :pi 
      shape = [shape] if shape.is_a? Integer

      #[base...limit]  -> [exponent_base ** base ... exponent_base ** limit]
      result = NMatrix.linspace(base, limit, shape)
      result.map {|element| exponent_base ** element}
    end

    #
    # call-seq:
    #     linspace(base, limit) -> 1x100 NMatrix
    #     linspace(base, limit, *shape) -> NMatrix
    #
    # Returns an NMatrix with +[shape[0] x shape[1] x .. x shape[dim-1]]+ values of dtype +:float64+ equally spaced from
    # +base+ to +limit+, inclusive.
    #
    # See: http://www.mathworks.com/help/matlab/ref/linspace.html
    #
    # * *Arguments* :
    #   - +base+ -> The first value in the sequence.
    #   - +limit+ -> The last value in the sequence.
    #   - +shape+ -> Desired output shape. Default returns a 1x100 row vector.
    # * *Returns* :
    #   - NMatrix with +:float64+ values.
    #
    # Examples :-
    #
    #   NMatrix.linspace(1,Math::PI, 6)
    #     =>[1.0,
    #        1.4283185005187988,
    #        1.8566370010375977,
    #        2.2849555015563965,
    #        2.7132740020751953,
    #        3.1415927410125732
    #       ]
    #
    #   NMatrix.linspace(1,10, [3,2])
    #     =>[
    #         [              1.0, 2.799999952316284]
    #         [4.599999904632568, 6.400000095367432]
    #         [8.199999809265137,              10.0]
    #       ]
    #
    def linspace(base, limit, shape = [100])
      
      # Convert shape to array format 
      shape = [shape] if shape.is_a? Integer 
      
      #Calculate number of elements 
      count = shape.inject(:*)
            
      # Linear spacing between elements calculated in step
      #   step = limit - base / (count - 1)
      #   [Result Sequence] = [0->N sequence] * step + [Base]
      step = (limit - base) * (1.0 / (count - 1))
      result = NMatrix.seq(shape, {:dtype => :float64}) * step
      result += NMatrix.new(shape, base)
      result
    end

    # call-seq:
    #     logspace(base, limit) -> 1x50 NMatrix with exponent_base = 10 
    #     logspace(base, limit, shape , exponent_base:) -> NMatrix
    #     logspace(base, :pi, n) -> 1xn NMatrix with interval [10 ^ base, Math::PI]
    #
    # Returns an NMatrix with +[shape[0] x shape[1] x .. x shape[dim-1]]+ values of dtype +:float64+ logarithmically spaced from
    # +exponent_base ^ base+ to +exponent_base ^ limit+, inclusive.
    #
    # See: http://www.mathworks.com/help/matlab/ref/logspace.html
    #
    # * *Arguments* :
    #   - +base+ -> exponent_base ** base is the first value in the sequence
    #   - +limit+ -> exponent_base ** limit is the last value in the sequence.
    #   - +shape+ -> Desired output shape. Default returns a 1x50 row vector.
    # * *Returns* :
    #   - NMatrix with +:float64+ values.
    #
    # Examples :-
    #
    #   NMatrix.logspace(1,:pi,7)
    #     =>[
    #         10.0000, 
    #         8.2450, 
    #         6.7980, 
    #         5.6050, 
    #         4.6213, 
    #         3.8103, 
    #         3.1416
    #       ]
    #
    #   NMatrix.logspace(1,2,[3,2])
    #     =>[
    #         [10.0, 15.8489]
    #         [25.1189, 39.8107]
    #         [63.0957, 100.0]
    #       ]
    #
    def logspace(base, limit, shape = [50], exponent_base: 10)

      #Calculate limit for [10 ^ base ... Math::PI] if limit = :pi
      limit = Math.log(Math::PI, exponent_base = 10) if limit == :pi 
      shape = [shape] if shape.is_a? Integer

      #[base...limit]  -> [exponent_base ** base ... exponent_base ** limit]
      result = NMatrix.linspace(base, limit, shape)
      result.map {|element| exponent_base ** element}
    end

    #
    # call-seq:
    #     seq(shape) -> NMatrix
    #     seq(shape, options) -> NMatrix
    #     bindgen(shape) -> NMatrix of :byte
    #     indgen(shape) -> NMatrix of :int64
    #     findgen(shape) -> NMatrix of :float32
    #     dindgen(shape) -> NMatrix of :float64
    #     cindgen(shape) -> NMatrix of :complex64
    #     zindgen(shape) -> NMatrix of :complex128
    #     rbindgen(shape) -> NMatrix of :object
    #
    # Creates a matrix filled with a sequence of integers starting at zero.
    #
    # * *Arguments* :
    #   - +shape+ -> Array (or integer for square matrix) specifying the dimensions.
    #   - +options+ -> (optional) Options permissible for NMatrix#initialize
    # * *Returns* :
    #   - NMatrix filled with values 0 through +size+.
    #
    # Examples:
    #
    #   NMatrix.seq(2) # =>   0   1
    #                 2   3
    #
    #   NMatrix.seq([3, 3], dtype: :float32) # =>  0.0  1.0  2.0
    #                                       3.0  4.0  5.0
    #                                       6.0  7.0  8.0
    #
    def seq(shape, options={})

      # Construct the values of the final matrix based on the dimension.
      values = (0 ... NMatrix.size(shape)).to_a

      # It'll produce :int32, except if a dtype is provided.
      NMatrix.new(shape, values, {:stype => :dense}.merge(options))
    end

    {:bindgen => :byte, :indgen => :int64, :findgen => :float32, :dindgen => :float64,
     :cindgen => :complex64, :zindgen => :complex128,
     :rbindgen => :object}.each_pair do |meth, dtype|
      define_method(meth) { |shape| NMatrix.seq(shape, :dtype => dtype) }
    end
  end
end

module NVector #:nodoc:

  class << self
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
    def new(*args)
      stype = args[0].is_a?(Symbol) ? args.shift : :dense
      shape = args[0].is_a?(Array) ? args.shift  : [1,args.shift]

      if shape.size != 2 || !shape.include?(1) || shape == [1,1]
        raise(ArgumentError, "shape must be a Fixnum or an Array of positive Fixnums where exactly one value is 1")
      end

      warn "NVector is deprecated and not guaranteed to work any longer"

      NMatrix.new(stype, shape, *args)
    end

    #
    # call-seq:
    #     zeros(size) -> NMatrix
    #     zeros(size, dtype) -> NMatrix
    #
    # Creates a new vector of zeros with the dimensions supplied as
    # parameters.
    #
    # * *Arguments* :
    #   - +size+ -> Array (or integer for square matrix) specifying the dimensions.
    #   - +dtype+ -> (optional) Default is +:float64+.
    # * *Returns* :
    #   - NVector filled with zeros.
    #
    # Examples:
    #
    #   NVector.zeros(2) # =>  0.0
    #                          0.0
    #
    #   NVector.zeros(3, :int32) # =>  0
    #                                  0
    #                                  0
    #
    def zeros(size, dtype = :float64)
      NMatrix.new([size,1], 0, dtype: dtype)
    end
    alias :zeroes :zeros

    #
    # call-seq:
    #     ones(size) -> NVector
    #     ones(size, dtype) -> NVector
    #
    # Creates a vector of ones with the dimensions supplied as
    # parameters.
    #
    # * *Arguments* :
    #   - +size+ -> Array (or integer for square matrix) specifying the dimensions.
    #   - +dtype+ -> (optional) Default is +:float64+.
    # * *Returns* :
    #   - NVector filled with ones.
    #
    # Examples:
    #
    #   NVector.ones(2) # =>  1.0
    #                         1.0
    #
    #   NVector.ones(3, :int32) # =>  1
    #                                 1
    #                                 1
    #
    def ones(size, dtype = :float64)
      NMatrix.new([size,1], 1, dtype: dtype)
    end

    #
    # call-seq:
    #     random(size) -> NVector
    #
    # Creates a vector with random numbers between 0 and 1 generated by
    # +Random::rand+ with the dimensions supplied as parameters.
    #
    # * *Arguments* :
    #   - +size+ -> Array (or integer for square matrix) specifying the dimensions.
    #   - +opts+ -> (optional) NMatrix#initialize options
    # * *Returns* :
    #   - NVector filled with random numbers generated by the +Random+ class.
    #
    # Examples:
    #
    #   NVector.rand(2) # =>  0.4859439730644226
    #                         0.1783195585012436
    #
    def random(size, opts = {})
      rng = Random.new

      random_values = []
      size.times { |i| random_values << rng.rand }

      NMatrix.new([size,1], random_values, opts)
    end

    #
    # call-seq:
    #     seq(n) -> NVector
    #     seq(n, dtype) -> NVector
    #
    # Creates a vector with a sequence of +n+ integers starting at zero. You
    # can choose other types based on the dtype parameter.
    #
    # * *Arguments* :
    #   - +n+ -> Number of integers in the sequence.
    #   - +dtype+ -> (optional) Default is +:int64+.
    # * *Returns* :
    #   - NVector filled with +n+ integers.
    #
    # Examples:
    #
    #   NVector.seq(2) # =>  0
    #                        1
    #
    #   NVector.seq(3, :float32) # =>  0.0
    #                                  1.0
    #                                  2.0
    #
    def seq(size, dtype = :int64)
      values = (0 ... size).to_a

      NMatrix.new([size,1], values, dtype: dtype)
    end

    #
    # call-seq:
    #     indgen(n) -> NVector
    #
    # Returns an integer NVector. Equivalent to <tt>seq(n, :int32)</tt>.
    #
    # * *Arguments* :
    #   - +n+ -> Size of the sequence.
    # * *Returns* :
    #   - NVector filled with +n+ integers of dtype +:int32+.
    #
    def indgen(n)
      NVector.seq(n, :int32)
    end

    #
    # call-seq:
    #     findgen(n) -> NVector
    #
    # Returns a float NVector. Equivalent to <tt>seq(n, :float32)</tt>.
    #
    # * *Arguments* :
    #   - +n+ -> Size of the sequence.
    # * *Returns* :
    #   - NVector filled with +n+ integers of dtype +:float32+.
    #
    def findgen(n)
      NVector.seq(n, :float32)
    end

    #
    # call-seq:
    #     bindgen(n) -> NVector
    #
    # Returns a byte NVector. Equivalent to <tt>seq(n, :byte)</tt>.
    #
    # * *Arguments* :
    #   - +n+ -> Size of the sequence.
    # * *Returns* :
    #   - NVector filled with +n+ integers of dtype +:byte+.
    #
    def bindgen(n)
      NVector.seq(n, :byte)
    end

    #
    # call-seq:
    #     cindgen(n) -> NVector
    #
    # Returns a complex NVector. Equivalent to <tt>seq(n, :complex64)</tt>.
    #
    # * *Arguments* :
    #   - +n+ -> Size of the sequence.
    # * *Returns* :
    #   - NVector filled with +n+ integers of dtype +:complex64+.
    #
    def cindgen(n)
      NVector.seq(n, :complex64)
    end

    #
    # call-seq:
    #     linspace(a, b) -> NVector
    #     linspace(a, b, n) -> NVector
    #
    # Returns a NVector with +n+ values of dtype +:float64+ equally spaced from
    # +a+ to +b+, inclusive.
    #
    # See: http://www.mathworks.com/help/matlab/ref/linspace.html
    #
    # * *Arguments* :
    #   - +a+ -> The first value in the sequence.
    #   - +b+ -> The last value in the sequence.
    #   - +n+ -> The number of elements. Default is 100.
    # * *Returns* :
    #   - NVector with +n+ +:float64+ values.
    #
    # Example:
    #   x = NVector.linspace(0, Math::PI, 1000)
    #   x.pretty_print
    #     [0.0
    #     0.0031447373909807737
    #     0.006289474781961547
    #     ...
    #     3.135303178807831
    #     3.138447916198812
    #     3.141592653589793]
    #   => nil
    #
    def linspace(a, b, n = 100)
      # Formula: seq(n) * step + a

      # step = ((b - a) / (n - 1))
      step = (b - a) * (1.0 / (n - 1))

      # dtype = :float64 is used to prevent integer coercion.
      result = NVector.seq(n, :float64) * NMatrix.new([n,1], step, dtype: :float64)
      result += NMatrix.new([n,1], a, dtype: :float64)
      result
    end

    #
    # call-seq:
    #     logspace(a, b) -> NVector
    #     logspace(a, b, n) -> NVector
    #
    # Returns a NVector with +n+ values of dtype +:float64+ logarithmically
    # spaced from +10^a+ to +10^b+, inclusive.
    #
    # See: http://www.mathworks.com/help/matlab/ref/logspace.html
    #
    # * *Arguments* :
    #   - +a+ -> The first value in the sequence.
    #   - +b+ -> The last value in the sequence.
    #   - +n+ -> The number of elements. Default is 100.
    # * *Returns* :
    #   - NVector with +n+ +:float64+ values.
    #
    # Example:
    #   x = NVector.logspace(0, Math::PI, 10)
    #   x.pretty_print
    #     [1.0
    #     2.2339109164570266
    #     4.990357982665873
    #     11.148015174505757
    #     24.903672795156997
    #     55.632586516975095
    #     124.27824233101062
    #     277.6265222213364
    #     620.1929186882427
    #     1385.4557313670107]
    #  => nil
    #
    def logspace(a, b, n = 100)
      # Formula: 10^a, 10^(a + step), ..., 10^b, where step = ((b-a) / (n-1)).

      result = NVector.linspace(a, b, n)
      result.each_stored_with_index { |element, i| result[i] = 10 ** element }
      result
    end
  end
end


# This constant is intended as a simple constructor for NMatrix meant for
# experimenting.
#
# Examples:
#
#   a = N[ 1,2,3,4 ]          =>  1  2  3  4
#
#   a = N[ 1,2,3,4, :int32 ]  =>  1  2  3  4
#
#   a = N[ [1,2,3], [3,4,5] ] =>  1  2  3
#                                 3  4  5
#
#   a = N[ 3,6,9 ].transpose => 3
#                               6
#                               9
N = NMatrix
