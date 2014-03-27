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
# == rspec.rb
#
# Monkey patches for RSpec improving its ability to work well with
# NMatrix (particularly #be_within).
#

require 'rspec'

# Amend RSpec to allow #be_within for matrices.
module RSpec::Matchers::BuiltIn #:nodoc:
  class BeWithin

    def of(expected)
      @expected = expected
      @unit     = ''
      if expected.is_a?(NMatrix)
        @tolerance = if @delta.is_a?(NMatrix)
                       @delta.abs
                     elsif @delta.is_a?(Array)
                       NMatrix.new(:dense, expected.shape, @delta, :object).abs.cast(:dtype => expected.abs_dtype)
                     else
                       (NMatrix.ones_like(expected) * @delta).abs
                     end
      else
        @tolerance = @delta
      end

      self
    end

    def percent_of(expected)
      @expected  = expected
      @unit      = '%'
      @tolerance = @expected.abs * @delta / 100.0 # <- only change is to reverse abs and @delta
      self
    end

    def matches?(actual)
      @actual = actual
      raise needs_expected     unless defined? @expected
      raise needs_subtractable unless @actual.respond_to? :-
      res = (@actual - @expected).abs <= @tolerance

      #if res.is_a?(NMatrix)
      #  require 'pry'
      #  binding.pry
      #end

      res.is_a?(NMatrix) ? !res.any? { |x| !x } : res
    end

  end
end