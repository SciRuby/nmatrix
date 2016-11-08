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
# == monkeys.rb
#
# Ruby core extensions for NMatrix.
#++

#######################
# Classes and Modules #
#######################

class Array
  # Convert a Ruby Array to an NMatrix.
  #
  # You must provide a shape for the matrix as the first argument.
  #
  # == Arguments:
  # <tt>shape</tt> :: Array describing matrix dimensions (or Fixnum for square).
  #   If not provided, will be intuited through #shape.
  # <tt>dtype</tt> :: Override data type (e.g., to store a Float as :float32
  #   instead of :float64) -- optional.
  # <tt>stype</tt> :: Optional storage type (defaults to :dense)
  def to_nm(shape = nil, dtype = nil, stype = :dense)
    elements = self.dup

    guess_dtype = ->(type) {
      case type
      when Fixnum   then :int64
      when Float    then :float64
      when Complex  then :complex128
      end
    }

    guess_shape = lambda { |shapey; shape|
      # Get the size of the current dimension
      shape = [shapey.size]
      shape << shapey.map {|s|
        if s.respond_to?(:size) && s.respond_to?(:map)
          guess_shape.call(s)
        else
          nil
        end
      }
      if shape.last.any? {|s| (s != shape.last.first) || s.nil?}
        shape.pop
      end
      if (shape.first != shape.last) && shape.last.all? {|s| s == shape.last.first}
        shape[-1] = shape.last.first
      end
      shape.flatten
    }

    unless shape
      shape = guess_shape.call(elements)
      elements.flatten!(shape.size - 1)
      if elements.flatten != elements
        dtype = :object
      else
        dtype ||= guess_dtype[elements[0]]
      end
    end

    dtype ||= guess_dtype[self[0]]

    matrix = NMatrix.new(:dense, shape, elements, dtype)

    if stype != :dense then matrix.cast(stype, dtype) else matrix end
  end
end

class Object #:nodoc:
  def returning(value)
    yield(value)
    value
  end
end


module Math #:nodoc:
  class << self
    NMatrix::NMMath::METHODS_ARITY_2.each do |meth|
      define_method "nm_#{meth}" do |arg0, arg1|
        if arg0.is_a? NMatrix then
          arg0.send(meth, arg1)
        elsif arg1.is_a? NMatrix then
          arg1.send(meth, arg0, true)
        else
          self.send("old_#{meth}".to_sym, arg0, arg1)
        end
      end
      alias_method "old_#{meth}".to_sym, meth
      alias_method meth, "nm_#{meth}".to_sym
    end
  end
end

class String
  def underscore
    self.gsub(/::/, '/').
    gsub(/([A-Z]+)([A-Z][a-z])/,'\1_\2').
    gsub(/([a-z\d])([A-Z])/,'\1_\2').
    tr("-", "_").
    downcase
  end
end

# Since `autoload` will most likely be deprecated (due to multi-threading concerns),
# we'll use `const_missing`. See: https://www.ruby-forum.com/topic/3036681 for more info.
module AutoloadPatch #:nodoc
  def const_missing(name)
    file = name.to_s.underscore
    require "nmatrix/io/#{file}"
    klass = const_get(name)
    return klass if klass
  end
end
