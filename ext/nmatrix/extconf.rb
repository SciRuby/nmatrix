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
# == extconf.rb
#
# This file checks for ATLAS and other necessary headers, and
# generates a Makefile for compiling NMatrix.

require File.expand_path("../../../lib/nmatrix/mkmf", __FILE__)

$INSTALLFILES = [
  ['nmatrix.h'       , '$(archdir)'], 
  ['nmatrix.hpp'     , '$(archdir)'],
  ['nmatrix_config.h', '$(archdir)'], 
  ['nm_memory.h'     , '$(archdir)'],
  ['ruby_constants.h', '$(archdir)']
]

if /cygwin|mingw/ =~ RUBY_PLATFORM
  $INSTALLFILES << ['libnmatrix.a', '$(archdir)']
end

$DEBUG = true
$CFLAGS = ["-Wall -Werror=return-type",$CFLAGS].join(" ")
$CXXFLAGS = ["-Wall -Werror=return-type",$CXXFLAGS].join(" ")
$CPPFLAGS = ["-Wall -Werror=return-type",$CPPFLAGS].join(" ")

# When adding objects here, make sure their directories are included in CLEANOBJS down at the bottom of extconf.rb.
basenames = %w{nmatrix ruby_constants data/data util/io math util/sl_list storage/common storage/storage storage/dense/dense storage/yale/yale storage/list/list}
$objs = basenames.map { |b| "#{b}.o"   }
$srcs = basenames.map { |b| "#{b}.cpp" }

#$libs += " -lprofiler "

create_conf_h("nmatrix_config.h")
create_makefile("nmatrix")

Dir.mkdir("data") unless Dir.exists?("data")
Dir.mkdir("util") unless Dir.exists?("util")
Dir.mkdir("storage") unless Dir.exists?("storage")
Dir.chdir("storage") do
  Dir.mkdir("yale")  unless Dir.exists?("yale")
  Dir.mkdir("list")  unless Dir.exists?("list")
  Dir.mkdir("dense") unless Dir.exists?("dense")
end

# to clean up object files in subdirectories:
open('Makefile', 'a') do |f|
  clean_objs_paths = %w{data storage storage/dense storage/yale storage/list util}.map { |d| "#{d}/*.#{CONFIG["OBJEXT"]}" }
  f.write("CLEANOBJS := $(CLEANOBJS) #{clean_objs_paths.join(' ')}")
end
