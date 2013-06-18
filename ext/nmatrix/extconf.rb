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
# SciRuby is Copyright (c) 2010 - 2013, Ruby Science Foundation
# NMatrix is Copyright (c) 2013, Ruby Science Foundation
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

require "mkmf"


# Function derived from NArray's extconf.rb.
def have_type(type, header=nil) #:nodoc:
  printf "checking for %s... ", type
  STDOUT.flush

  src = <<"SRC"
#include <ruby.h>
SRC


  src << <<"SRC" unless header.nil?
#include <#{header}>
SRC

  r = try_link(src + <<"SRC")
  int main() { return 0; }
  int t() { #{type} a; return 0; }
SRC

  unless r
    print "no\n"
    return false
  end

  $defs.push(format("-DHAVE_%s", type.upcase))

  print "yes\n"

  return true
end

# Function derived from NArray's extconf.rb.
def create_conf_h(file) #:nodoc:
  print "creating #{file}\n"
  File.open(file, 'w') do |hfile|
    header_guard = file.upcase.sub(/\s|\./, '_')

    hfile.puts "#ifndef #{header_guard}"
    hfile.puts "#define #{header_guard}"
    hfile.puts

    # FIXME: Find a better way to do this:
    if RUBY_VERSION >= '2.0'
      hfile.puts "#define RUBY_2 1"
    end

    for line in $defs
      line =~ /^-D(.*)/
      hfile.printf "#define %s 1\n", $1
    end

    hfile.puts
    hfile.puts "#endif"
  end
end

if RUBY_VERSION < '1.9'
  raise(NotImplementedError, "Sorry, you need at least Ruby 1.9!")
else
  $INSTALLFILES = [['nmatrix.h', '$(archdir)'], ['nmatrix.hpp', '$(archdir)'], ['nmatrix_config.h', '$(archdir)']]
  if /cygwin|mingw/ =~ RUBY_PLATFORM
    $INSTALLFILES << ['libnmatrix.a', '$(archdir)']
  end
end

if /cygwin|mingw/ =~ RUBY_PLATFORM
  CONFIG["DLDFLAGS"] << " --output-lib libnmatrix.a"
end

$DEBUG = true
$CFLAGS = ["-Wall ",$CFLAGS].join(" ")

$srcs = [
         'nmatrix.cpp',
         'ruby_constants.cpp',

         'data/data.cpp',
         'util/math.cpp',
         'util/sl_list.cpp',
         'util/io.cpp',
         'storage/common.cpp',
         'storage/storage.cpp',
         'storage/dense.cpp',
         'storage/yale.cpp',
         'storage/list.cpp'
        ]
# add smmp in to get generic transp; remove smmp2 to eliminate funcptr transp

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


unless have_library("lapack") and have_header("clapack.h")
  dir_config("lapack", ["/usr/include/atlas"], ["/usr/local/lib", "/usr/local/atlas/lib"])
end

unless have_library("cblas") and have_header("cblas.h")
  dir_config("cblas", ["/usr/local/atlas/include", "/usr/include/atlas"], ["/usr/local/lib", "/usr/local/atlas/lib"])
end

unless have_library("atlas")
  dir_config("atlas", ["/usr/local/atlas/include", "/usr/include/atlas"], ["/usr/local/atlas/lib", "/usr/local/lib", "/usr/lib"])
end

#find_library("lapack", "clapack_dgetrf")
have_header("clapack.h")
have_header("cblas.h")

have_func("clapack_dgetrf", "clapack.h")


#find_library("cblas", "cblas_dgemm")
#find_library("atlas", "ATL_dgemmNN")

have_func("cblas_dgemm", "cblas.h")

# Order matters here: ATLAS has to go after LAPACK: http://mail.scipy.org/pipermail/scipy-user/2007-January/010717.html
$libs += " -llapack -lcblas -latlas "

$objs = %w{nmatrix ruby_constants data/data util/io util/math util/sl_list storage/common storage/storage storage/dense storage/yale storage/list}.map { |i| i + ".o" }

#CONFIG['CXX'] = 'clang++'
CONFIG['CXX'] = 'g++'

def find_newer_gplusplus #:nodoc:
  [8,7,6,5,4,3].each do |minor|
    result = `which g++-4.#{minor}`
    next if result.empty?
    CONFIG['CXX'] = "g++-4.#{minor}"
    return CONFIG['CXX']
  end
  false
end

def gplusplus_version #:nodoc:
  `LANG="en_US" #{CONFIG['CXX']} -v 2>&1`.lines.to_a.last.match(/gcc\sversion\s(\d\.\d.\d)/).captures.first
end


if CONFIG['CXX'] == 'clang++'
  $CPP_STANDARD = 'c++11'

else
  version = gplusplus_version
  if version < '4.3.0' && CONFIG['CXX'] == 'g++'  # see if we can find a newer G++, unless it's been overridden by user
    if !find_newer_gplusplus
      raise("You need a version of g++ which supports -std=c++0x or -std=c++11. If you're on a Mac and using Homebrew, we recommend using mac-brew-gcc.sh to install a more recent g++.")
    end
    version = gplusplus_version
  end

  if version < '4.7.0'
    $CPP_STANDARD = 'c++0x'
  else
    $CPP_STANDARD = 'c++11'
  end
end

# For release, these next two should both be changed to -O3.
$CFLAGS += " -O0 -g "
$CPPFLAGS += " -O0 -g -std=#{$CPP_STANDARD} " #-fmax-errors=10 -save-temps

CONFIG['warnflags'].gsub!('-Wshorten-64-to-32', '') # doesn't work except in Mac-patched gcc (4.2)
CONFIG['warnflags'].gsub!('-Wdeclaration-after-statement', '')
CONFIG['warnflags'].gsub!('-Wimplicit-function-declaration', '')

create_conf_h("nmatrix_config.h")
create_makefile("nmatrix")

Dir.mkdir("data") unless Dir.exists?("data")
Dir.mkdir("util") unless Dir.exists?("util")
Dir.mkdir("storage") unless Dir.exists?("storage")

# to clean up object files in subdirectories:
open('Makefile', 'a') do |f|
  f.write <<EOS
CLEANOBJS := $(CLEANOBJS) data/*.#{CONFIG["OBJEXT"]} storage/*.#{CONFIG["OBJEXT"]} util/*.#{CONFIG["OBJEXT"]}
EOS
end
