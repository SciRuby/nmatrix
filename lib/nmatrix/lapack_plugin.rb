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
# == lapack_plugin.rb
#
# This file `require`s either nmatrix-atlas or nmatrix-lapacke depending on which
# is available.
#
# The idea is that if a developer wants to use a LAPACK feature which is provided
# by both of these gems (e.g. NMatrix#potrf! or NMatrix::LAPACK.geev),
# but doesn't care which one is installed, they can
# just `require 'nmatrix/lapack_plugin'` rather than having to choose between
# `require 'nmatrix/lapacke'` or `require 'nmatrix/lapacke'` 
#++

begin
  require 'nmatrix/atlas'
rescue LoadError
  begin
    require 'nmatrix/lapacke'
  rescue LoadError
    raise(LoadError,"Either nmatrix-atlas or nmatrix-lapacke must be installed")
  end
end
