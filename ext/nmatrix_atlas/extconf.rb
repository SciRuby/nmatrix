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
$CFLAGS = ["-Wall -Werror=return-type -I$(srcdir)/../nmatrix",$CFLAGS].join(" ")
$CXXFLAGS = ["-Wall -Werror=return-type -I$(srcdir)/../nmatrix",$CXXFLAGS].join(" ")
$CPPFLAGS = ["-Wall -Werror=return-type -I$(srcdir)/../nmatrix",$CPPFLAGS].join(" ")

# When adding objects here, make sure their directories are included in CLEANOBJS down at the bottom of extconf.rb.
# Why not just autogenerate this list from all .c/.cpp files in directory?
basenames = %w{nmatrix_atlas math_atlas}
$objs = basenames.map { |b| "#{b}.o"   }
$srcs = basenames.map { |b| "#{b}.cpp" }

# The next line allows the user to supply --with-atlas-dir=/usr/local/atlas,
# --with-atlas-lib or --with-atlas-include and tell the compiler where to look
# for ATLAS. The same for all the others
#
#dir_config("clapack", ["/usr/local/atlas/include"], [])
#
#

# Is g++ having trouble finding your header files?
# Try this:
#   export C_INCLUDE_PATH=/usr/local/atlas/include
#   export CPLUS_INCLUDE_PATH=/usr/local/atlas/include
# (substituting in the path of your cblas.h and clapack.h for the path I used). -- JW 8/27/12

idefaults = {lapack: ["/usr/include/atlas"],
             cblas: ["/usr/local/atlas/include", "/usr/include/atlas"],
             atlas: ["/usr/local/atlas/include", "/usr/include/atlas"]}

# For some reason, if we try to look for /usr/lib64/atlas on a Mac OS X Mavericks system, and the directory does not
# exist, it will give a linker error -- even if the lib dir is already correctly included with -L. So we need to check
# that Dir.exists?(d) for each.
ldefaults = {lapack: ["/usr/local/lib", "/usr/local/atlas/lib", "/usr/lib64/atlas"].delete_if { |d| !Dir.exists?(d) },
             cblas: ["/usr/local/lib", "/usr/local/atlas/lib", "/usr/lib64/atlas"].delete_if { |d| !Dir.exists?(d) },
             atlas: ["/usr/local/lib", "/usr/local/atlas/lib", "/usr/lib", "/usr/lib64/atlas"].delete_if { |d| !Dir.exists?(d) }}

if have_library("clapack") # Usually only applies for Mac OS X
  $libs += " -lclapack "
end

unless have_library("lapack")
  dir_config("lapack", idefaults[:lapack], ldefaults[:lapack])
end

unless have_library("cblas")
  dir_config("cblas", idefaults[:cblas], ldefaults[:cblas])
end

unless have_library("atlas")
  dir_config("atlas", idefaults[:atlas], ldefaults[:atlas])
end

# If BLAS and LAPACK headers are in an atlas directory, prefer those. Otherwise,
# we try our luck with the default location.
if have_header("atlas/cblas.h")
  have_header("atlas/clapack.h")
else
  have_header("cblas.h")
  have_header("clapack.h")
end


# Although have_func is supposed to take a list as its second argument, I find that it simply
# applies a :to_s to the second arg and doesn't actually check each one. We may want to put
# have_func calls inside an :each block which checks atlas/clapack.h, cblas.h, clapack.h, and
# lastly lapack.h. On Ubuntu, it only works if I use atlas/clapack.h. --@mohawkjohn 8/20/14
have_func("clapack_dgetrf", "atlas/clapack.h")
have_func("clapack_dgetri", "atlas/clapack.h")
have_func("dgesvd_", "clapack.h") # This may not do anything. dgesvd_ seems to be in LAPACK, not CLAPACK.

have_func("cblas_dgemm", "cblas.h")

#have_func("rb_scan_args", "ruby.h")

#find_library("lapack", "clapack_dgetrf")
#find_library("cblas", "cblas_dgemm")
#find_library("atlas", "ATL_dgemmNN")
# Order matters here: ATLAS has to go after LAPACK: http://mail.scipy.org/pipermail/scipy-user/2007-January/010717.html
$libs += " -llapack -lcblas -latlas "
#$libs += " -lprofiler "

create_conf_h("nmatrix_atlas_config.h")
create_makefile("nmatrix_atlas")

# to clean up object files in subdirectories:
open('Makefile', 'a') do |f|
  clean_objs_paths = %w{ }.map { |d| "#{d}/*.#{CONFIG["OBJEXT"]}" }
  f.write("CLEANOBJS := $(CLEANOBJS) #{clean_objs_paths.join(' ')}")
end
