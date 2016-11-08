require 'java'
require_relative '../../../ext/nmatrix_java/vendor/commons-math3-3.6.1.jar'
require_relative '../../../ext/nmatrix_java/target/nmatrix.jar'

java_import 'org.apache.commons.math3.linear.ArrayRealVector'
java_import 'org.apache.commons.math3.linear.RealMatrix'
java_import 'org.apache.commons.math3.linear.MatrixUtils'
java_import 'org.apache.commons.math3.linear.DecompositionSolver'
java_import 'org.apache.commons.math3.linear.LUDecomposition'
java_import 'org.apache.commons.math3.linear.QRDecomposition'
java_import 'org.apache.commons.math3.linear.CholeskyDecomposition'
java_import 'MatrixGenerator'
java_import 'ArrayGenerator'
java_import 'MathHelper'
java_import 'ArrayComparator'

class NMatrix
  include_package 'org.apache.commons.math3.analysis.function'
  attr_accessor :shape, :dim, :dtype, :stype, :s

  def initialize(*args)
    if args[-1] == :copy
      @shape = [2,2]
      @s = [0,0,0,0]
      @dim = shape.is_a?(Array) ? shape.length : 2
    else
      if (args.length <= 3)
        @shape = args[0]
        if args[1].is_a?(Array)
          elements = args[1]
          if args.length > 2
            hash = args[2]
            # puts hash
            @dtype = hash[:dtype]
            @stype = hash[:stype]
          else
            @dtype = :double
            @stype = :dense
          end
        else
          # elements = Java::double[shape[0]*shape[1]].new{ Java::Double.NaN }
          if args.length > 1
            if args[1].is_a?(Symbol)
              hash = args[1]
              @dtype = hash[:dtype]
              @stype = hash[:stype]
              elements = Array.new(shape*shape) unless shape.is_a? Array
            else
              elements = Array.new(shape*shape) unless shape.is_a? Array
            end
          end
        end
      else

        offset = 0
        if (!args[0].is_a?(Symbol) && !args[0].is_a?(String))
          @stype = :dense
        else
          offset = 1
          @stype = :dense
          @dtype = args[-1]
        end

        @shape = args[offset]
        elements = args[offset+1]

      end


      @shape = [shape,shape] unless shape.is_a?(Array)
      # @dtype = interpret_dtype(argc-1-offset, argv+offset+1, stype);
      # @dtype = args[:dtype] if args[:dtype]
      @dtype_sym = nil
      @stype_sym = nil
      @default_val_num = nil
      @capacity_num = nil
      @size = (0...@shape.size).inject(1) { |x,i| x * @shape[i] }

      j=0

      if (elements.is_a?(ArrayRealVector))
        @s = elements
      # elsif elements.java_class.to_s == "[D"
      #   @s = ArrayRealVector.new(elements)
      else
        storage = Array.new(size)
        elements = [elements,elements] unless elements.is_a?(Array)
        if size > elements.length
          (0...size).each do |i|
            j=0 unless j!=elements.length
            storage[i] = elements[j]
            j+=1
          end
        else
          storage = elements
        end
        if @dtype == :object
          @s = storage
        else
          @s = ArrayRealVector.new(storage.to_java Java::double)
        end
      end

      @dim = @shape.is_a?(Array) ? @shape.length : 2

    end
  end

  # Needs to be properly implemented
  def self.min_dtype(alpha)
    :int8
  end

  def self.upcast(alpha, beta)
    false
  end

  def clone
    result = create_dummy_nmatrix
    # ArrayRealVector#clone is disable, hence use copy
    # that returns a deep copy of the object.
    result.s = @s.copy
    return result
  end

  def entries
    return @s.toArray.to_a
  end

  def twoDMat
    return MatrixUtils.createRealMatrix MatrixGenerator.getMatrixDouble(self.s.toArray, @shape[0], @shape[1])
  end

  def dtype
    return @dtype
  end

  #FIXME
  def self.guess_dtype arg
    :float32
  end

  def stype
    @stype = :dense
  end

  def cast_full *args
    return self
  end

  def default_value
    return nil
  end

  def __list_default_value__
    #not implemented currently
  end

  def __yale_default_value__
    #not implemented currently
  end

  def [] *args
    return xslice(args)
  end

  def slice(*args)
    return xslice(args)
  end

  def []=(*args)
    to_return = nil
    if args.length > @dim+1
      raise(ArgumentError, "wrong number of arguments (#{args.length} for #{effective_dim(dim+1)})" )
    else
      slice = get_slice(@dim, args, @shape)
      dense_storage_set(slice, args[-1])
      to_return = args[-1]
    end
    return to_return
  end

  def is_ref?

  end

  # def dim
  #   shape.is_a?(Array) ? shape.length : 2
  # end

  alias :dimensions :dim

  def effective_dim(s)
    d = 0
    (0...@dim).each do |i|
      d+=1 unless @shape[i] == 1
    end
    return d
  end

  alias :effective_dimensions :effective_dim



  protected

  def create_dummy_nmatrix
    nmatrix = NMatrix.new(:copy)
    nmatrix.shape = self.shape
    nmatrix.dim = self.dim
    nmatrix.dtype = self.dtype
    nmatrix.stype = self.stype
    return nmatrix
  end

  def __list_to_hash__

  end

  public

  def shape
    @shape
  end

   def supershape s
    if (s[:src] == @s)
      return shape
       # easy case (not a slice)
    else
      @s = s[:src]
    end

    new_shape = Array.new(dim)
    (0...dim).each do |index|
      new_shape[index] = shape[index]
    end

    return new_shape
  end

  def offset
    # ArrayRealVector takes care of the offset value when indexing a Vector.
    # Hence, return 0.
    0
  end

  def det_exact
    # if (:stype != :dense)
    #   raise Exception.new("can only calculate exact determinant for dense matrices")
    #   return nil
    # end

    if (@dim != 2 || @shape[0] != @shape[1])
      raise(ShapeError, "matrices must be square to have a determinant defined")
      return nil
    end
    to_return = nil
    if (dtype == :object)
      # to_return = *reinterpret_cast<VALUE*>(result);
    else
      to_return = LUDecomposition.new(self.twoDMat).getDeterminant()
    end

    return to_return.round(3)
  end

  alias :det :det_exact

  def complex_conjugate!

  end


  protected

  def count_max_elements
    return size
  end

  def reshape_bang arg
    if(@stype == :dense)
      shape_ary = arg
      size = count_max_elements
      new_size = 1
      shape = interpret_shape(shape_ary, dim)

      (0...dim).each do |index|
        new_size *= shape[index]
      end

      if (size == new_size)
        self.shape = shape
        self.dim = dim
        return self
      else
         raise(ArgumentError, "reshape cannot resize; size of new and old matrices must match")
      end
    else
      raise(NotImplementedError, "reshape in place only for dense stype")
    end
  end

  def interpret_shape(shape_ary, dim)
    shape = []

    if shape_ary.is_a?(Array)
      dim = shape_ary.length

      (0...dim).each do |index|
        shape[index] = shape_ary[index].to_i
      end

    elsif shape_ary.is_a?(FIXNUM)
      dim = 2
      shape = Array.new(dim)
      shape[0] = shape_ary.to_i
      shape[1] = shape_ary.to_i
    else
      raise(ArgumentError, "Expected an array of numbers or a single Fixnum for matrix shape")
    end

    return shape
  end


  public

  def each_with_indices
    nmatrix = create_dummy_nmatrix
    stride = get_stride(self)
    offset = 0
    #Create indices and initialize them to zero
    coords = Array.new(dim){ 0 }

    shape_copy =  Array.new(dim)
    (0...size).each do |k|
      dense_storage_coords(nmatrix, k, coords, stride, offset)
      slice_index = dense_storage_pos(coords,stride)
      ary = Array.new
      if (@dtype == :object)
        ary << self.s[slice_index]
      else
        ary << self.s.toArray.to_a[slice_index]
      end
      (0...dim).each do |p|
        ary << coords[p]
      end

      # yield the array which now consists of the value and the indices
      yield(ary)
    end if block_given?

    return nmatrix
  end


  def each_stored_with_indices
    nmatrix = create_dummy_nmatrix
    stride = get_stride(self)
    offset = 0
    #Create indices and initialize them to zero
    coords = Array.new(dim){ 0 }

    shape_copy =  Array.new(dim)

    (0...size).each do |k|
      dense_storage_coords(nmatrix, k, coords, stride, offset)
      slice_index = dense_storage_pos(coords,stride)
      ary = Array.new
      if (@dtype == :object)
        ary << self.s[slice_index]
      else
        ary << self.s.toArray.to_a[slice_index]
      end
      (0...dim).each do |p|
        ary << coords[p]
      end
      # yield the array which now consists of the value and the indices
      yield(ary)
    end if block_given?

    return nmatrix
  end

  def map_stored

  end

  def each_ordered_stored_with_indices

  end


  protected

  def __dense_each__
    nmatrix = create_dummy_nmatrix
    stride = get_stride(self)
    offset = 0
    #Create indices and initialize them to zero
    coords = Array.new(dim){ 0 }

    shape_copy =  Array.new(dim)
    (0...size).each do |k|
      if (@dtype == :object)
        dense_storage_coords(nmatrix, k, coords, stride, offset)
        slice_index = dense_storage_pos(coords,stride)
        yield self.s[slice_index]
      else
        dense_storage_coords(nmatrix, k, coords, stride, offset)
        slice_index = dense_storage_pos(coords,stride)
        yield self.s.toArray.to_a[slice_index]
      end
    end if block_given?
    if (@dtype == :object)
      return @s.to_enum
    else
      return @s.toArray().to_a.to_enum
    end
  end

  def __dense_map__
    nmatrix = create_dummy_nmatrix
    stride = get_stride(self)
    offset = 0
    coords = Array.new(dim){ 0 }
    shape_copy =  Array.new(dim)

    s= Java::double[size].new
    (0...size).each do |k|
      dense_storage_coords(nmatrix, k, coords, stride, offset)
      slice_index = dense_storage_pos(coords,stride)

      y = yield @s.getEntry(slice_index)
      @s.setEntry(slice_index, y)
    end
    nmatrix.s = ArrayRealVector.new s

    return nmatrix
  end

  def __dense_map_pair__

  end

  def __list_map_merged_stored__

  end

  def __list_map_stored__

  end

  def __yale_map_merged_stored__

  end

  def __yale_map_stored__

  end

  def __yale_stored_diagonal_each_with_indices__

  end

  def __yale_stored_nondiagonal_each_with_indices__

  end


  public

  def ==(otherNmatrix)
    result = false
    if (otherNmatrix.is_a?(NMatrix))
      #check dimension
      if (@dim != otherNmatrix.dim)
        raise(ShapeError, "cannot compare matrices with different dimension")
      end
      #check shape
      (0...dim).each do |i|
        if (@shape[i] != otherNmatrix.shape[i])
          raise(ShapeError, "cannot compare matrices with different shapes");
        end
      end

      #check the entries
      if dtype == :object
        result = @s == otherNmatrix.s
      else
        result = ArrayComparator.equals(@s.toArray, otherNmatrix.s.toArray)
      end
    end
    result
  end

  def =~ (other)
    lha = @s.toArray.to_a
    rha = other.s.toArray.to_a
    resultArray = Array.new(lha.length)
    if (other.is_a?(NMatrix))
      #check dimension
      if (@dim != other.dim)
        raise(ShapeError, "cannot compare matrices with different dimension")
        return nil
      end
      #check shape
      (0...dim).each do |i|
        if (@shape[i] != other.shape[i])
          raise(ShapeError, "cannot compare matrices with different shapes");
          return nil
        end
      end
      #check the entries
      (0...lha.length).each do |i|
        resultArray[i] = lha[i] == rha[i] ? true : false
      end
      result = NMatrix.new(:copy)
      result.shape = @shape
      result.dtype = :object
      result.s = resultArray
    end
    result
  end

  def !~ (other)
    lha = @s.toArray.to_a
    rha = other.s.toArray.to_a
    resultArray = Array.new(lha.length)
    if (other.is_a?(NMatrix))
      #check dimension
      if (@dim != other.dim)
        raise(ShapeError, "cannot compare matrices with different dimension")
        return nil
      end
      #check shape
      (0...dim).each do |i|
        if (@shape[i] != other.shape[i])
          raise(ShapeError, "cannot compare matrices with different shapes");
          return nil
        end
      end
      #check the entries
      (0...lha.length).each do |i|
        resultArray[i] = lha[i] != rha[i] ? true : false
      end
      result = NMatrix.new(:copy)
      result.shape = @shape
      result.dtype = :object
      result.s = resultArray
    end
    result
  end

  def <= (other)
    lha = @s.toArray.to_a
    rha = other.s.toArray.to_a
    resultArray = Array.new(lha.length)
    if (other.is_a?(NMatrix))
      #check dimension
      if (@dim != other.dim)
        raise(ShapeError, "cannot compare matrices with different dimension")
        return nil
      end
      #check shape
      (0...dim).each do |i|
        if (@shape[i] != other.shape[i])
          raise(ShapeError, "cannot compare matrices with different shapes");
          return nil
        end
      end
      #check the entries
      (0...lha.length).each do |i|
        resultArray[i] = lha[i] <= rha[i] ? true : false
      end
      result = NMatrix.new(:copy)
      result.shape = @shape
      result.dtype = :object
      result.s = resultArray
    end
    result
  end

  def >= (other)
    lha = @s.toArray.to_a
    rha = other.s.toArray.to_a
    resultArray = Array.new(lha.length)
    if (other.is_a?(NMatrix))
      #check dimension
      if (@dim != other.dim)
        raise(ShapeError, "cannot compare matrices with different dimension")
        return nil
      end
      #check shape
      (0...dim).each do |i|
        if (@shape[i] != other.shape[i])
          raise(ShapeError, "cannot compare matrices with different shapes");
          return nil
        end
      end
      #check the entries
      (0...lha.length).each do |i|
        resultArray[i] = lha[i] >= rha[i] ? true : false
      end
      result = NMatrix.new(:copy)
      result.shape = @shape
      result.dtype = :object
      result.s = resultArray
    end
    result
  end

  def < (other)
    lha = @s.toArray.to_a
    rha = other.s.toArray.to_a
    resultArray = Array.new(lha.length)
    if (other.is_a?(NMatrix))
      #check dimension
      if (@dim != other.dim)
        raise(ShapeError, "cannot compare matrices with different dimension")
        return nil
      end
      #check shape
      (0...dim).each do |i|
        if (@shape[i] != other.shape[i])
          raise(ShapeError, "cannot compare matrices with different shapes");
          return nil
        end
      end
      #check the entries
      (0...lha.length).each do |i|
        resultArray[i] = lha[i] < rha[i] ? true : false
      end
      result = NMatrix.new(:copy)
      result.shape = @shape
      result.dtype = :object
      result.s = resultArray
    end
    result
  end

  def > (other)
    lha = @s.toArray.to_a
    rha = other.s.toArray.to_a
    resultArray = Array.new(lha.length)
    if (other.is_a?(NMatrix))
      #check dimension
      if (@dim != other.dim)
        raise(ShapeError, "cannot compare matrices with different dimension")
        return nil
      end
      #check shape
      (0...dim).each do |i|
        if (@shape[i] != other.shape[i])
          raise(ShapeError, "cannot compare matrices with different shapes");
          return nil
        end
      end
      #check the entries
      (0...lha.length).each do |i|
        resultArray[i] = lha[i] > rha[i] ? true : false
      end
      result = NMatrix.new(:copy)
      result.shape = @shape
      result.dtype = :object
      result.s = resultArray
    end
    result
  end

  # /////////////////////////
  # // Matrix Math Methods //
  # /////////////////////////

  def dot(other)
    result = nil
    if (other.is_a?(NMatrix))
      #check dimension
      if (@shape.length!=2 || other.shape.length!=2)
        raise(NotImplementedError, "please convert array to nx1 or 1xn NMatrix first")
        return nil
      end
      #check shape
      if (@shape[1] != other.shape[0])
        raise(ArgumentError, "incompatible dimensions")
        return nil
      end

      # if(@stype != other.stype)
      #   raise(NotImplementedError, "matrices must have same stype")
      # end

      result = create_dummy_nmatrix
      result.shape = [@shape[0],other.shape[1]]
      twoDMat = self.twoDMat.multiply(other.twoDMat)
      result.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(twoDMat.getData, @shape[0],other.shape[1]))
    else
      raise(ArgumentError, "cannot have dot product with a scalar");
    end
    return result;
  end

  def symmetric?
    return is_symmetric(false)
  end

  def is_symmetric(hermitian)
    is_symmetric = false

    if (@shape[0] == @shape[1] and @dim == 2)
      if @stype == :dense
        if (hermitian)
          #Currently, we are not dealing with complex matrices.
          eps = 0
          is_symmetric = MatrixUtils.isSymmetric(self.twoDMat, eps)
        else
          eps = 0
          is_symmetric = MatrixUtils.isSymmetric(self.twoDMat, eps)
        end

      else
        #TODO: Implement, at the very least, yale_is_symmetric. Model it after yale/transp.template.c.
        raise Exception.new("symmetric? and hermitian? only implemented for dense currently")
      end
    end
    return is_symmetric ? true : false
  end

  def hermitian?
    return is_symmetric(true)
  end

  def capacity

  end

  # // protected methods

  protected

  def __inverse__
    # if (:stype != :dense)
    #   raise Exception.new("needs exact determinant implementation for this matrix stype")
    #   return nil
    # end

    if (@dim != 2 || @shape[0] != @shape[1])
      raise Exception.new("matrices must be square to have an inverse defined")
      return nil
    end
    to_return = nil
    if (dtype == :RUBYOBJ)
      # to_return = *reinterpret_cast<VALUE*>(result);
    else
      to_return = create_dummy_nmatrix
      twoDMat = MatrixUtils.inverse(self.twoDMat)
      to_return.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(twoDMat.getData, @shape[0], @shape[1]))
    end

    return to_return
  end

  def __inverse__!
    # if (:stype != :dense)
    #   raise Exception.new("needs exact determinant implementation for this matrix stype")
    #   return nil
    # end

    if (@dim != 2 || @shape[0] != @shape[1])
      raise Exception.new("matrices must be square to have an inverse defined")
      return nil
    end
    to_return = nil
    if (dtype == :RUBYOBJ)
      # to_return = *reinterpret_cast<VALUE*>(result);
    else
      twoDMat = MatrixUtils.inverse(self.twoDMat)
      @s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(twoDMat.getData, @shape[0], @shape[1]))
    end

    return self
  end

  def __inverse_exact__
    # if (:stype != :dense)
    #   raise Exception.new("needs exact determinant implementation for this matrix stype")
    #   return nil
    # end

    if (@dim != 2 || @shape[0] != @shape[1])
      raise Exception.new("matrices must be square to have an inverse defined")
      return nil
    end
    to_return = nil
    if (dtype == :RUBYOBJ)
      # to_return = *reinterpret_cast<VALUE*>(result);
    else
      to_return = create_dummy_nmatrix
      twoDMat = MatrixUtils.inverse(self.twoDMat)
      to_return.s = ArrayRealVector.new(ArrayGenerator.getArrayDouble(twoDMat.getData, @shape[0], @shape[1]))
    end

    return to_return

  end

  private

  # // private methods

  def __hessenberg__(param)
    raise(NotImplementedError, "Hessenberg Transformer not implemented for NMatrix-JRuby")
  end
end

# load jruby implementation of operators.
require_relative './slice.rb'
require_relative './math.rb'
require_relative './operators.rb'
require_relative './decomposition.rb'
require_relative './error.rb'
require_relative './enumerable.rb'