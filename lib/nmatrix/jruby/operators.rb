class NMatrix

  # A dummy matrix is a matrix without the elements atrribute.
  # NMatrix#create_dummy_matrix prevents creating copies as @s is set explicitly.
  def +(other)
    result = create_dummy_nmatrix
    if (other.is_a?(NMatrix))
      #check dimension
      raise(ShapeError, "Cannot add matrices with different dimension") if (@dim != other.dim)
      #check shape
      (0...dim).each do |i|
        raise(ShapeError, "Cannot add matrices with different shapes") if (@shape[i] != other.shape[i])
      end
      result.s = @s.copy.add(other.s)
    else
      result.s = @s.copy.mapAddToSelf(other)
    end
    result
  end

  def -(other)
    result = create_dummy_nmatrix
    if (other.is_a?(NMatrix))
      #check dimension
      raise(ShapeError, "Cannot subtract matrices with different dimension") if (@dim != other.dim)
      #check shape
      (0...dim).each do |i|
        raise(ShapeError, "Cannot subtract matrices with different shapes") if (@shape[i] != other.shape[i])
      end
      result.s = @s.copy.subtract(other.s)
    else
      result.s = @s.copy.mapSubtractToSelf(other)
    end
    result
  end

  def *(other)
    result = create_dummy_nmatrix
    if (other.is_a?(NMatrix))
      #check dimension
      raise(ShapeError, "Cannot multiply matrices with different dimension") if (@dim != other.dim)
      #check shape
      (0...dim).each do |i|
        raise(ShapeError, "Cannot multiply matrices with different shapes") if (@shape[i] != other.shape[i])
      end
      result.s = @s.copy.ebeMultiply(other.s)
    else
      result.s = @s.copy.mapMultiplyToSelf(other)
    end
    result
  end

  def /(other)
    result = create_dummy_nmatrix
    if (other.is_a?(NMatrix))
      #check dimension
      raise(ShapeError, "Cannot divide matrices with different dimension") if (@dim != other.dim)
      #check shape
      (0...dim).each do |i|
        raise(ShapeError, "Cannot divide matrices with different shapes") if (@shape[i] != other.shape[i])
      end
      result.s = @s.copy.ebeDivide(other.s)
    else
      result.s = @s.copy.mapDivideToSelf(other)
    end
    result
  end

  def ** val
    result = NMatrix.new(:copy)
    result.shape = @shape
    result.dim = @dim
    result.s = @s.copy.mapToSelf(Power.new(val))
    result
  end

  def %(other)
    raise Exception.new("modulus not supported in NMatrix-jruby")
  end

  def atan2(other, scalar=false)
    result = create_dummy_nmatrix
    if scalar
      result.s = ArrayRealVector.new MathHelper.atan2Scalar(other, @s.toArray)
    else
      if other.is_a? NMatrix
        result.s = ArrayRealVector.new MathHelper.atan2(other.s.toArray, @s.toArray)
      else
        result.s = ArrayRealVector.new MathHelper.atan2Scalar2(other, @s.toArray)
      end
    end
    result
  end

  def ldexp(other, scalar=false)
    result = create_dummy_nmatrix
    if scalar
      result.s = ArrayRealVector.new MathHelper.ldexpScalar(other, @s.toArray)
    else
      if other.is_a? NMatrix
        result.s = ArrayRealVector.new MathHelper.ldexp(other.s.toArray, @s.toArray)
      else
        result.s = ArrayRealVector.new MathHelper.ldexpScalar2(other, @s.toArray)
      end
    end
    result
  end

  def hypot(other, scalar=false)
    result = create_dummy_nmatrix
    if scalar
      result.s = ArrayRealVector.new MathHelper.hypotScalar(other, @s.toArray)
    else
      if other.is_a? NMatrix
        result.s = ArrayRealVector.new MathHelper.hypot(other.s.toArray, @s.toArray)
      else
        result.s = ArrayRealVector.new MathHelper.hypotScalar(other, @s.toArray)
      end
    end
    result
  end

  def sin
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Sin.new())
    result
  end

  def cos
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Cos.new())
    result
  end

  def tan
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Tan.new())
    result
  end

  def asin
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Asin.new())
    result
  end

  def acos
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Acos.new())
    result
  end

  def atan
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Atan.new())
    result
  end

  def sinh
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Sinh.new())
    result
  end

  def cosh
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Cosh.new())
    result
  end

  def tanh
    result = NMatrix.new(:copy)
    result.shape = @shape
    result.dim = @dim
    result.s = @s.copy.mapToSelf(Tanh.new())
    result
  end

  def asinh
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Asinh.new())
    result
  end

  def acosh
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Acosh.new())
    result
  end

  def atanh
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Atanh.new())
    result
  end

  def exp
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Exp.new())
    result
  end

  def log(val = :natural)
    result = create_dummy_nmatrix
    if val == :natural
      result.s = @s.copy.mapToSelf(Log.new())
    else
      result.s = ArrayRealVector.new MathHelper.log(val, @s.toArray)
    end
    result
  end

  def log2
    self.log(2)
  end

  def log10
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Log10.new())
    result
  end

  def sqrt
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Sqrt.new())
    result
  end

  def erf
    result = create_dummy_nmatrix
    result.s = ArrayRealVector.new MathHelper.erf(@s.toArray)
    result
  end

  def erfc
    result = create_dummy_nmatrix
    result.s = ArrayRealVector.new MathHelper.erfc(@s.toArray)
    result
  end

  def cbrt
    result = create_dummy_nmatrix
    result.s = @s.copy.mapToSelf(Cbrt.new())
    result
  end

  def gamma
    result = create_dummy_nmatrix
    result.s = ArrayRealVector.new MathHelper.gamma(@s.toArray)
    result
  end

  def -@
    result = create_dummy_nmatrix
    result.s = @s.copy.mapMultiplyToSelf(-1)
    result
  end

  def floor
    result = create_dummy_nmatrix
    # Need to be changed later
    result.dtype = :int16
    result.s = @s.copy.mapToSelf(Floor.new())
    result
  end

  def ceil
    result = create_dummy_nmatrix
    # Need to be changed later
    result.dtype = :int16
    result.s = @s.copy.mapToSelf(Ceil.new())
    result
  end

  def round
    result = create_dummy_nmatrix
    # Need to be changed later
    result.dtype = :int16
    result.s = ArrayRealVector.new MathHelper.round(@s.toArray)
    result
  end

end