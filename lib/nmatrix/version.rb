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
# SciRuby is Copyright (c) 2010 - 2016, Ruby Science Foundation
# NMatrix is Copyright (c) 2012 - 2016, John Woods and the Ruby Science Foundation
#
# Please see LICENSE.txt for additional copyright notices.
#
# == Contributing
#
# By contributing source code to SciRuby, you agree to be bound by
# our Contributor Agreement:
#
# * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
#++

class NMatrix
  # Note that the format of the VERSION string is needed for NMatrix
  # native IO. If you change the format, please make sure that native
  # IO can still understand NMatrix::VERSION.
  module VERSION #:nodoc:
    MAJOR = 0
    MINOR = 2
    TINY = 5
    #PRE = "a"

    STRING = [MAJOR, MINOR, TINY].compact.join(".")
    #STRING = [MAJOR, MINOR, TINY, PRE].compact.join(".")
  end
end

