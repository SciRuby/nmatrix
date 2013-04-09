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
# == extconf.rb
#
# Configuration file for the NMatrix C/C++ extension.

require 'mkmf'

# Configure our compilers and compiler options.
$CFLAGS							= '-I. -fPIC -Wall -pipe -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector -mtune=native -O3'
CONFIG['CXXFLAGS']	= '-std=c++11'

if clang_path = find_executable('clang') and clang_pp_path = find_executable('clang++')
  CONFIG['CC']	= clang_path
  CONFIG['CXX']	= clang_pp_path
end

# Necessary header files.
find_header('ruby/config.h')

# List the source files that need to be compiled.
$srcs = [
         #	'nmatrix.cpp',
         'ruby_constants.cpp',

         'data/data.cpp',
         'math.cpp',
         'storage/storage.cpp',
         'storage/dense.cpp'
        ]

$objs = $srcs.map { |f| f.sub('.cpp', '.o') }

# Create the actual Makefile.
create_makefile('NMatrix')
