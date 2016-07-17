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

require "nmatrix/mkmf"

#$INSTALLFILES = [['nmatrix.h', '$(archdir)'], ['nmatrix.hpp', '$(archdir)'], ['nmatrix_config.h', '$(archdir)'], ['nm_memory.h', '$(archdir)']]
if /cygwin|mingw/ =~ RUBY_PLATFORM
  #$INSTALLFILES << ['libnmatrix.a', '$(archdir)']
end

$DEBUG = true
#not the right way to add this include directory
$CFLAGS = ["-Wall -Werror=return-type -I$(srcdir)/../nmatrix -I$(srcdir)/lapacke/include",$CFLAGS].join(" ")
$CXXFLAGS = ["-Wall -Werror=return-type -I$(srcdir)/../nmatrix -I$(srcdir)/lapacke/include",$CXXFLAGS].join(" ")
$CPPFLAGS = ["-Wall -Werror=return-type -I$(srcdir)/../nmatrix -I$(srcdir)/lapacke/include",$CPPFLAGS].join(" ")

# When adding objects here, make sure their directories are included in CLEANOBJS down at the bottom of extconf.rb.
# Why not just autogenerate this list from all .c/.cpp files in directory?
basenames = %w{nmatrix_lapacke math_lapacke lapacke}
$objs = basenames.map { |b| "#{b}.o"   }
$srcs = basenames.map { |b| "#{b}.cpp" }

# For some reason, if we try to look for /usr/lib64/atlas on a Mac OS X Mavericks system, and the directory does not
# exist, it will give a linker error -- even if the lib dir is already correctly included with -L. So we need to check
# that Dir.exists?(d) for each.
ldefaults = {lapack: ["/usr/local/lib"].delete_if { |d| !Dir.exists?(d) } }

# It is not clear how this variable should be defined, or if it is necessary at all. 
# See issue https://github.com/SciRuby/nmatrix/issues/403
idefaults = {lapack: [] }

unless have_library("lapack")
  dir_config("lapack", idefaults[:lapack], ldefaults[:lapack])
end

# Order matters here: ATLAS has to go after LAPACK: http://mail.scipy.org/pipermail/scipy-user/2007-January/010717.html
$libs += " -llapack "
#To use the Intel MKL, comment out the line above, and also comment out the bit above with have_library and dir_config for lapack.
#Then add something like the line below (for exactly what linker flags to use see https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor ):
#$libs += " -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential "

create_conf_h("nmatrix_lapacke_config.h")
create_makefile("nmatrix_lapacke")

# to clean up object files in subdirectories:
open('Makefile', 'a') do |f|
  clean_objs_paths = %w{ }.map { |d| "#{d}/*.#{CONFIG["OBJEXT"]}" }
  f.write("CLEANOBJS := $(CLEANOBJS) #{clean_objs_paths.join(' ')}")
end
