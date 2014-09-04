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
# == interpolation/one_dimensional.rb
#
# One dimensional interpolation methods for NMatrix.
# 
#++

require_relative 'base.rb'

class NMatrix
  module Interpolation

    # Implements one dimensional interpolation routines.
    # 
    # ==== Usage
    # 
    # x = NMatrix.seq [10]
    # y = x.exp
    # 
    # f = NMatrix::Interpolation::OneDimensional.new x, y, {kind: :linear, sorted: true}
    # i = f.interpolate 2.5
    # 
    # puts "Interpolated value for 2.5 is #{i}"
    class OneDimensional < Interpolation::Base

      # Constructor for all One Dimensional interpolation operations.
      # 
      # The function values to be supplied to this class are of the form y = f(x). 
      # 
      # Henceforth, y will be referred to as ordinate and x as absicca. If absicca 
      # and ordinate arrays are not of the same length, then the effective size used 
      # for interpolation will be MIN(x.size, y.size).
      # 
      # ==== Arguments
      # 
      # * +x+ -    The collection of absiccas. Must be a 1 D NMatrix or ruby Array.
      # 
      # * +y+ -    The collection of ordinates corresponding to the absicca. 'y' can 
      #            either be a 1D NMatrix or Array OR a 2D NMatrix. In case y contains 
      #            multiple columns, the interpolation is carried out on each column,
      #            unless specified.
      # 
      # * +opts+ - Various options for carrying out the interpolation, to be specified as
      #            a hash.
      # 
      # ==== Options
      # 
      # * +:kind+ - The kind of interpolation that the user wants to perform. Should be 
      #             specified as a symbol. This option is compulsory. Only linear
      #             interpolation supported as of now. 
      # 
      # * +:sorted+ - Set this option as *true* if the absicca collection is supplied in
      #             the arguments in a sorted manner. If not supplied, it will be assumed
      #             that absiccas are not sorted and they will sorted be sorted anyway.
      # 
      # * +:axis+ - In case of a multidimensional ordinate matrix, specify the column over
      #             which interpolation must be performed. *axis* starts indexing from 0 
      #             and should be lower than the number of columns in the ordinate matrix.
      # 
      # * +:precision+ - Specifies the precision of the interpolated values returned. Defaults
      #             to 3.
      # 
      def initialize x, y, opts={}
        super(x,y,opts)
      end

      # Performs the actual interpolation on the value passed as an argument. Kind of 
      # interpolation performed is determined according to what is specified in the 
      # constructor.
      # 
      # ==== Arguments
      # 
      # * +interpolant+ - The value for which the interpolation is to be performed. Can
      #                   either be a Numeric, Array of Numerics or NMatrix. If multidimensional
      #                   NMatrix is supplied then will flatten it and interpolate over
      #                   all its values. Will return answer in the form of an NMatrix if
      #                   *interpolant* is supplied as an NMatrix.
      def interpolate interpolant
        result = case @opts[:kind]
        when :linear
          for_each (interpolant) { |x| linear_interpolation(x)  }
        else
          raise(ArgumentError, "1 D interpolation of kind #{@opts[:kind]} \
           not supported")
        end

        return result.to_nm if interpolant.is_a?(NMatrix)

        result
      end

     private

      def for_each interpolant
        result = []

        if interpolant.kind_of? Numeric
          return yield interpolant
        else
          interpolant.each { |x| result << yield(x) }
        end

        result
      end

      def linear_interpolation interpolant
        # TODO : Make this more efficient by using hunt/locate from NR
        index  = locate(interpolant)
        same   = @x[index] == @x[index+1]
        result = []              

        if (@y.respond_to?(:vector?) and @y.vector?) or
          @y.instance_of?(Array)

          return @y[index] if same
          return _lin_interpolator @y, index, interpolant
        elsif @opts[:axis]

          return @y.column(@opts[:axis])[index] if same
          return _lin_interpolator @y.column(@opts[:axis]), index, interpolant
        else
          @y.each_column do |c|
            result << (same ? c[index] : _lin_interpolator(c, index, interpolant))
          end
        end

        result
      end

      def _lin_interpolator y, index, interpolant
        (y[index] + 
        ((interpolant - @x[index]) / (@x[index + 1] - @x[index])) * 
         (y[index + 1] - y[index])).round(@opts[:precision])
      end
    end
  end
end