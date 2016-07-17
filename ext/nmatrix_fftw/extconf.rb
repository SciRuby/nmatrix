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
# == nmatrix_fftw/extconf.rb
#
# This file checks FFTW3 and other necessary headers/shared objects.

require 'nmatrix/mkmf'

fftw_libdir = RbConfig::CONFIG['libdir']
fftw_incdir = RbConfig::CONFIG['includedir']
fftw_srcdir = RbConfig::CONFIG['srcdir']

$CFLAGS = ["-Wall -Werror=return-type -I$(srcdir)/../nmatrix -I$(srcdir)/lapacke/include",$CFLAGS].join(" ")
$CXXFLAGS = ["-Wall -Werror=return-type -I$(srcdir)/../nmatrix -I$(srcdir)/lapacke/include -std=c++11",$CXXFLAGS].join(" ")
$CPPFLAGS = ["-Wall -Werror=return-type -I$(srcdir)/../nmatrix -I$(srcdir)/lapacke/include -std=c++11",$CPPFLAGS].join(" ")

flags = " --include=#{fftw_incdir} --libdir=#{fftw_libdir}"

if have_library("fftw3")
  $CFLAGS += [" -lfftw3 -lm #{$CFLAGS} #{$flags}"].join(" ")
  dir_config('nmatrix_fftw', fftw_incdir, fftw_libdir)
  dir_config('nmatrix_fftw')
end

create_conf_h("nmatrix_fftw_config.h")
create_makefile("nmatrix_fftw")

# to clean up object files in subdirectories:
open('Makefile', 'a') do |f|
  clean_objs_paths = %w{ }.map { |d| "#{d}/*.#{CONFIG["OBJEXT"]}" }
  f.write("CLEANOBJS := $(CLEANOBJS) #{clean_objs_paths.join(' ')}")
end
