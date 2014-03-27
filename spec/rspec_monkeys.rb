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
# == rspec_monkeys.rb
#
# A set of monkey patches for RSpec allowing checks of NMatrix types
#

module RSpec::Matchers::BuiltIn
  class BeWithin

    def of(expected)
      @expected = expected
      @unit     = ''
      if expected.is_a?(NMatrix)
        @tolerance = if @delta.is_a?(NMatrix)
                       @delta.clone
                     elsif @delta.is_a?(Array)
                       NMatrix.new(:dense, expected.shape, @delta, expected.dtype)
                     else
                       NMatrix.ones_like(expected) * @delta
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
  end
end