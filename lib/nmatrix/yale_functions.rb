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
# == yale_functions.rb
#
# This file contains some shortcut functions for the specialty
# Yale matrix extensions (mostly for debugging and experimental
# purposes, but sometimes applicable when you need to speed up
# your code a lot).
#++

module NMatrix::YaleFunctions
  # call-seq:
  #     yale_nd_row_size(i) -> Fixnum
  #
  # Returns the size of a given non-diagonal row.
  def yale_nd_row_size i
    yale_ija(i+1) - yale_ija(i)
  end

  # call-seq:
  #     yale_ja_at(i) -> Array
  #
  # Returns the non-diagonal column indices which are stored in a given row.
  def yale_ja_at i
    yale_nd_row(i, :keys)
  end
  alias :yale_nd_row_as_array :yale_ja_at

  # call-seq:
  #     yale_ja_set_at(i) -> Set
  #
  # Returns the non-diagonal column indices which are stored in a given row, as a Set.
  def yale_ja_set_at i
    require 'set'
    yale_nd_row(i, :keys).to_set
  end
  alias :yale_nd_row_as_set :yale_ja_set_at

  # call-seq:
  #     yale_ja_sorted_set_at -> SortedSet
  #
  # Returns the non-diagonal column indices which are stored in a given row, as a Set.
  def yale_ja_sorted_set_at i
    require 'set'
    SortedSet.new(yale_nd_row(i, :keys))
  end
  alias :yale_nd_row_as_sorted_set :yale_ja_sorted_set_at

  # call-seq:
  #     yale_nd_row_as_hash(i) -> Hash
  #
  # Returns the non-diagonal column indices and entries stored in a given row.
  def yale_nd_row_as_hash i
    yale_nd_row(i, :hash)
  end

  # call-seq:
  #     yale_ja_d_keys_at(i) -> Array
  #
  # Returns the diagonal and non-digonal column indices stored in a given row.
  def yale_ja_d_keys_at i
    ary = yale_nd_row(i, :keys)
    return ary if i >= self.shape[1] || self[i,i] == self.default_value
    ary << i
  end
  alias :yale_row_as_array :yale_ja_d_keys_at

  # call-seq:
  #     yale_ja_d_keys_set_at(i) -> Set
  #
  # Returns the diagonal and non-diagonal column indices stored in a given row.
  def yale_ja_d_keys_set_at i
    require 'set'
    yale_ja_d_keys_at(i).to_set
  end
  alias :yale_row_as_set :yale_ja_d_keys_set_at

  # call-seq:
  #     yale_ja_d_keys_sorted_set_at(i) -> SortedSet
  #
  # Returns the diagonal and non-diagonal column indices stored in a given row.
  def yale_ja_d_keys_sorted_set_at i
    require 'set'
    SortedSet.new(yale_row_as_array(i))
  end
  alias :yale_row_as_sorted_set :yale_ja_d_keys_sorted_set_at

  # call-seq:
  #     yale_row_as_hash(i) -> Hash
  #
  # Returns the diagonal and non-diagonal column indices and entries stored in a given row.
  def yale_row_as_hash i
    h = yale_nd_row(i, :hash)
    return h if i >= self.shape[1] || self[i,i] == self.default_value
    h[i] = self[i,i]
  end
end