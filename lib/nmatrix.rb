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
# == nmatrix.rb
#
# This file loads the C extension for NMatrix and all the ruby files.
#

# For some reason nmatrix.so ends up in a different place during gem build.
if File.exist?("lib/nmatrix/nmatrix.so") #|| File.exist?("lib/nmatrix/nmatrix.bundle")
  # Development
  require "nmatrix/nmatrix.so"
else
  # Gem
  require "nmatrix.so"
end

require 'nmatrix/nmatrix.rb'
require 'nmatrix/version.rb'
#require 'nmatrix/nvector.rb'
require 'nmatrix/blas.rb'
require 'nmatrix/monkeys'
require "nmatrix/shortcuts.rb"
