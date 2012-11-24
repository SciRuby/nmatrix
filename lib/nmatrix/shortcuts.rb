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
# == shortcuts.rb
#
# These are shortcuts for NMatrix and NVector creation, contributed by Daniel
# Carrera (dcarrera@hush.com) and Carlos Agarie (carlos@onox.com.br).

class NMatrix

  class << self
    # zeros() or zeroes()
    #
    #     Creates a new matrix of zeros with the dimensions supplied as
    #     parameters. Optional parameters include:
    #
    #     * A storage type as the first parameter (default is :dense).
    #     * A dtype as the last parameter (default is :float64).
    #
    # Examples:
    #
    #   zeros(2) # =>  0.0   0.0   
    #                  0.0   0.0
    #
    #   zeros([2, 3], :int32) # =>  0  0  0
    #                            0  0  0
    #
    #   zeros(:list, [1, 5], :int32) # =>  0  0  0  0  0
    #
    
    def zeros(*params)
      dtype = params.last.is_a?(Symbol) ? params.pop : :float64
      stype = params.first.is_a?(Symbol) ? params.shift : :dense
      dim = params.first

      NMatrix.new(stype, dim, 0, dtype)
    end

    alias :zeroes :zeros

    # ones()
    #
    #     Creates a :dense matrix of ones with the dimensions supplied
    #     as parameters. Optionaly, one can specify a dtype as the last
    #     parameter (default is :float64).
    #
    # Examples:
    #
    #   ones([1, 3]) # =>  1.0   1.0   1.0
    #
    #   ones([2, 3], :int32) # =>  1  1  1
    #                              1  1  1
    #
    
    def ones(*params)
      dtype = params.last.is_a?(Symbol) ? params.pop : :float64
      dim = params.first

      NMatrix.new(dim, 1, dtype)
    end

    # identity() or eye()
    #
    #     Creates an identity matrix (square matrix rank 2) of the size
    #     supplied as a parameter. Optional parameters include:
    #
    #     * A storage type as the first parameter (default is :dense).
    #     * A dtype as the last parameter (default is :float64).
    #
    # Examples:
    #
    #    eye(3) # =>   1.0   0.0   0.0
    #                  0.0   1.0   0.0
    #                  0.0   0.0   1.0
    #
    #    eye(3, :int32) # =>   1   0   0
    #                          0   1   0
    #                          0   0   1
    #
    #    eye(:yale, 2, :int32) # =>   1   0
    #                                 0   1
    #

    def eye(*params)
      dtype = params.last.is_a?(Symbol) ? params.pop : :float64
      stype = params.first.is_a?(Symbol) ? params.shift : :dense

      dim = params.first
      
      # Fill the diagonal with 1's.
      m = NMatrix.zeros(stype, dim, dtype)
      (0 .. (dim - 1)).each do |i| 
        m[i, i] = 1
      end
      
      m
    end

    alias :identity :eye
    
    # random()
    #
    #     Creates a :dense NMatrix with random numbers between 0 and 1 generated
    #     by Random::rand. The parameter is the dimension of the matrix.
    #
    # Examples:
    #
    #   rand([2, 2]) # => 0.4859439730644226   0.1783195585012436
    #                     0.23193766176700592  0.4503345191478729
    #

    def random(*params)
      dim = params.first
      rng = Random.new

      # Must provide the dimension as an Integer for a square matrix or as an
      # array, e.g. [2, 4, 7].
      unless dim.is_a?(Integer) || dim.is_a?(Array)
        raise ArgumentError, "random() accepts only integers or arrays as \
dimension."
      end

      random_values = []
      
      # Construct the values of the final matrix based on the dimension.
      if dim.is_a?(Integer)
        (dim * dim - 1).times { |i| random_values << rng.rand }
      else
        # Dimensions given by an array. Get the product of the array elements
        # and generate this number of random values.
        dim.reduce(1, :*).times { |i| random_values << rng.rand }
      end

      NMatrix.new(:dense, dim, random_values, :float64)
    end

    # seq()
    #
    #     Creates a :dense NMatrix with a sequence of integers starting at
    #     zero until the matrix is filled. The parameters to the method
    #     are the dimensions of the matrix. Optionaly, one can specify a
    #     dtype as the last parameter (default is :float64).
    #
    # Examples:
    #
    #   seq(2) # =>   0   1
    #                 2   3
    #
    #   seq([3, 3], :float32) # =>  0.0  1.0  2.0
    #                               3.0  4.0  5.0
    #                               6.0  7.0  8.0
    #
    
    def seq(*params)
      dtype = params.last.is_a?(Symbol) ? params.pop : nil
      dim = params.first
      
      # Must provide the dimension as an Integer for a square matrix or as an
      # 2 element array, e.g. [2,4].
      unless dim.is_a?(Integer) || (dim.is_a?(Array) && dim.size < 3)
        raise ArgumentError, "seq() accepts only integers or 2-element arrays \
as dimension."
      end
      
      # Construct the values of the final matrix based on the dimension.
      if dim.is_a?(Integer)
        values = (0 .. (dim * dim - 1)).to_a
      else
        # Dimensions given by a 2 element array.
        values = (0 .. (dim.first * dim.last - 1)).to_a
      end
      
      # It'll produce :int32, except if a dtype is provided.
      NMatrix.new(:dense, dim, values, dtype)
    end

    #########################################
    # FUNCTIONS FOR MATLAB AND IDL REFUGEES #
    #########################################

    #
    # These are functions that replicate existing functionality, but
    # would probably be appreciated by MATLAB or IDL users.
    #

    # indgen() , findgen() , bindgen() , cindgen()
    #
    #      These IDL functions are similar to seq() but less flexible.
    #      They produce one-dimensional vectors:
    #
    #      indgen    -- Integer vector   --   seq(n, :int32)
    #      findgen   -- Float vector     --   seq(n, :float32)
    #      bindgen   -- Byte vector      --   seq(n, :byte)
    #      cindgen   -- Complex vector   --   seq(n, :complex64)
    #

    def indgen(n)
      NMatrix.seq(n, :int32)
    end

    def findgen(n)
      NMatrix.seq(n, :float32)
    end

    def bindgen(n)
      NMatrix.seq(n, :byte)
    end

    def cindgen(n)
      NMatrix.seq(n, :complex64)
    end
  
  end

  #
  # These shortcuts are to be called directly from a NMatrix object, i.e.:
  #
  # >> m = NMatrix.random(3)
  # >> m.column(2)
  #

  # column()
  #
  #     Returns the column specified. The second parameter defaults to
  #     :copy, which returns a copy of the selected column, but it can be
  #     specified as :reference, which will return a reference to it.
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
    unless [:copy, :reference].include?(get_by)
      raise ArgumentError, "column() 2nd parameter must be :copy or :reference"
    end
    
    if get_by == :copy
      self.slice(0 ... self.shape[0], column_number)
    else # by reference
      self[0 ... self.shape[0], column_number]
    end
  end
end

class NVector < NMatrix
  
  class << self
    # zeros() or zeroes()
    #
    #     Creates a new matrix of zeros with the dimensions supplied as
    #     parameters. Optional parameters include:
    #
    #     * A storage type as the first parameter (default is :dense).
    #     * A dtype as the last parameter (default is :float64).
    #
    # Examples:
    #
    #   zeros(2) # =>  0.0   0.0   
    #
    #   zeros(3, :int32) # =>  0  0  0
    #
    
    def zeros(*params)
      dtype = params.last.is_a?(Symbol) ? params.pop : :float64
      dim = params.first

      NVector.new(dim, 0, dtype)
    end

    alias :zeroes :zeros

    # ones()
    #
    #     Creates a :dense matrix of ones with the dimensions supplied
    #     as parameters. Optionaly, one can specify a dtype as the last
    #     parameter (default is :float64).
    #
    # Examples:
    #
    #   ones(3) # =>  1.0   1.0   1.0
    #
    #   ones(2, :int32) # =>  1  1
    #
    
    def ones(*params)
      dtype = params.last.is_a?(Symbol) ? params.pop : :float64
      dim = params.first

      NVector.new(dim, 1, dtype)
    end
    
    # random()
    #
    #     Creates a :dense NMatrix with random numbers between 0 and 1 generated
    #     by Random::rand. The parameter is the dimension of the matrix.
    #
    # Examples:
    #
    #   rand(2) # => 0.4859439730644226   0.1783195585012436
    #

    def random(*params)
      rng = Random.new
      dim = params.first

      random_values = []
      dim.times { |i| random_values << rng.rand }
      
      NVector.new(dim, random_values, :float64)
    end

    # seq()
    #
    #     Creates a :dense NMatrix with a sequence of integers starting at
    #     zero until the matrix is filled. The parameters to the method
    #     are the dimensions of the matrix. Optionaly, one can specify a
    #     dtype as the last parameter (default is :float64).
    #
    # Examples:
    #
    #   seq(2) # =>   0   1
    #
    #   seq(3, :float32) # =>  0.0  1.0  2.0
    #
    
    def seq(*params)
      dtype = params.last.is_a?(Symbol) ? params.pop : nil
      dim = params.first
      
      unless dim.is_a?(Integer)
        raise ArgumentError, "NVector::seq() only accepts integers as \
dimension."
      end
            
      values = (0 .. (dim - 1)).to_a
      
      NVector.new(dim, values, dtype)
    end

    #########################################
    # FUNCTIONS FOR MATLAB AND IDL REFUGEES #
    #########################################

    #
    # These are functions that replicate existing functionality, but
    # would probably be appreciated by MATLAB or IDL users.
    #

    # indgen() , findgen() , bindgen() , cindgen()
    #
    #      These IDL functions are similar to seq() but less flexible.
    #      They produce one-dimensional vectors:
    #
    #      indgen    -- Integer vector   --   seq(n, :int32)
    #      findgen   -- Float vector     --   seq(n, :float32)
    #      bindgen   -- Byte vector      --   seq(n, :byte)
    #      cindgen   -- Complex vector   --   seq(n, :complex64)
    #

    def indgen(n)
      NVector.seq(n, :int32)
    end

    def findgen(n)
      NVector.seq(n, :float32)
    end

    def bindgen(n)
      NVector.seq(n, :byte)
    end

    def cindgen(n)
      NVector.seq(n, :complex64)
    end
   
    # linspace()
    #
    #      This MATLAB function somewhat resembles seq(), but it differs
    #      enough that is likely to be legitimately useful to non-MATLAB
    #      refugees. This function takes three parameter, the last one
    #      being an integer.
    #
    #      linspace( a, b, n )
    #
    #      This returns a vector with n values equally spaced from a to b,
    #      inclusive.
    #      
    #      Following the MATLAB implementation, if n isn't provided it's
    #      assumed to be 100.
    #
    # Ex:  x = linspace(0, pi, 1000)
    #      y = sin(x)
    #

    def linspace(a, b, n = 100)
      # See: http://www.mathworks.com/help/matlab/ref/linspace.html
      # Formula:  seq(n) * step + a
      
      # step = ((b - a) / (n - 1))
      step = (b - a) * (1.0 / (n - 1))
      
      # dtype = :float64 is used to prevent integer coercion.
      result = NVector.seq(n, :float64) * NVector.new(n, step, :float64)
      result += NVector.new(n, a, :float64)
      result
    end
  end
end

# NMatrix needs to have a succinct way to create a matrix by specifying
# the components directly. This is very usefeul for using NMatrix as an
# advanced calculator, it is useful for learning NMatrix and it is also
# useful for testing language features or developing algorithms.
#
# The N[] function provides a way to create a matrix in a way that is
# very short and very natural, simply by specifying the components in
# the traditional Ruby array syntax.  Optionally, one can specify a
# dtype as the last parameter (default is :float64).
#
# a = N[ 1,2,3,4 ]          =>  1.0  2.0  3.0  4.0
#
# a = N[ 1,2,3,4, :int32 ]  =>  1  2  3  4
#
# a = N[ [1,2,3], [3,4,5] ] =>  1.0  2.0  3.0
#                               3.0  4.0  5.0
#
#
# SYNTAX COMPARISON:
#
#     MATLAB:		a = [ [1 2 3] ; [4 5 6] ]   or  [ 1 2 3 ; 4 5 6 ]
#     IDL:			a = [ [1,2,3] , [4,5,6] ]
#     NumPy:		a = array( [1,2,3], [4,5,6] )
#
#     SciRuby:      a = N[ [1,2,3], [4,5,6] ]
#     Ruby array:   a =  [ [1,2,3], [4,5,6] ]
#

class N
  class << self
    def [](*params)
      dtype = params.last.is_a?(Symbol) ? params.pop : nil

      # First find the dimensions of the array.
      i = 0
      dim = []
      foo = params
      while foo.is_a?(Array)
        dim[i] = foo.length
        foo    = foo[0]
        i += 1
      end

      # Then flatten the array.
      NMatrix.new(dim, params.flatten, dtype)
    end
  end
end

# TODO Make all the shortcuts available through modules, allowing someone
# to include them to make "MATLAB-like" scripts.
#
# There are some questions to be answered before this can be done, tho.