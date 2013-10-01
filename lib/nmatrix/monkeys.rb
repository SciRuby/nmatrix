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
# == monkeys.rb
#
# Ruby core extensions for NMatrix.
#++

require 'nmatrix/math'

#######################
# Classes and Modules #
#######################

class Array
  # Convert a Ruby Array to an NMatrix.
  #
  # You must provide a shape for the matrix as the first argument.
  #
  # == Arguments:
  # <tt>shape</tt> :: Array describing matrix dimensions (or Fixnum for square) -- REQUIRED!
  # <tt>dtype</tt> :: Override data type (e.g., to store a Float as :float32 instead of :float64) -- optional.
  # <tt>stype</tt> :: Optional storage type (defaults to :dense)
  def to_nm(shape, dtype = nil, stype = :dense)
    dtype ||=
      case self[0]
      when Fixnum		then :int64
      when Float		then :float64
      when Rational	then :rational128
      when Complex	then :complex128
      end

    matrix = NMatrix.new(:dense, shape, self, dtype)

    if stype != :dense then matrix.cast(stype, dtype) else matrix end
  end
end

class Object #:nodoc:
  def returning(value)
    yield(value)
    value
  end
end


module Math
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

