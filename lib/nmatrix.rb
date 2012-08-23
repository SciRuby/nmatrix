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
# == nmatrix.rb
#
# This file loads the C extension for NMatrix and adds an autoload for the
# NMatrix and NVector classes.

############
# Requires #
############

# NMatrix

# For some reason nmatrix.so ends up in a different place during gem build.
if File.exist? 'lib/nmatrix/nmatrix.so'
	# Development
	require 'nmatrix/nmatrix.so'
else
	# Gem
	require 'nmatrix.so'
end

# Monkey patches.
require 'nmatrix/monkeys'

#############
# Autoloads #
#############

autoload :NMatrix, 'nmatrix/nmatrix'
autoload :NVector, 'nmatrix/nvector'

require "nmatrix/helpers"