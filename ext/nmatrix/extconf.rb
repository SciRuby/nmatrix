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
    hfile.puts "#define RUBY_2 1" if RUBY_VERSION >= '2.0'

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
  $INSTALLFILES = [['nmatrix.h', '$(archdir)'], ['nmatrix.hpp', '$(archdir)'], ['nmatrix_config.h', '$(archdir)'], ['nm_memory.h', '$(archdir)']]
  if /cygwin|mingw/ =~ RUBY_PLATFORM
    $INSTALLFILES << ['libnmatrix.a', '$(archdir)']
  end
end

if /cygwin|mingw/ =~ RUBY_PLATFORM
  CONFIG["DLDFLAGS"] << " --output-lib libnmatrix.a"
end

$DEBUG = true
$CFLAGS = ["-Wall -Werror=return-type",$CFLAGS].join(" ")
$CXXFLAGS = ["-Wall -Werror=return-type",$CXXFLAGS].join(" ")
$CPPFLAGS = ["-Wall -Werror=return-type",$CPPFLAGS].join(" ")

# When adding objects here, make sure their directories are included in CLEANOBJS down at the bottom of extconf.rb.
basenames = %w{nmatrix ruby_constants data/data util/io math util/sl_list storage/common storage/storage storage/dense/dense storage/yale/yale storage/list/list}
$objs = basenames.map { |b| "#{b}.o"   }
$srcs = basenames.map { |b| "#{b}.cpp" }

#CONFIG['CXX'] = 'clang++'
CONFIG['CXX'] = 'g++'

def find_newer_gplusplus #:nodoc:
  print "checking for apparent GNU g++ binary with C++0x/C++11 support... "
  [9,8,7,6,5,4,3].each do |minor|
    ver = "4.#{minor}"
    gpp = "g++-#{ver}"
    result = `which #{gpp}`
    next if result.empty?
    CONFIG['CXX'] = gpp
    puts ver
    return CONFIG['CXX']
  end
  false
end

def gplusplus_version
  cxxvar = proc { |n| `#{CONFIG['CXX']} -E -dM - </dev/null | grep #{n}`.chomp.split(' ')[2] }
  major = cxxvar.call('__GNUC__')
  minor = cxxvar.call('__GNUC_MINOR__')
  patch = cxxvar.call('__GNUC_PATCHLEVEL__')

  raise("unable to determine g++ version (match to get version was nil)") if major.nil? || minor.nil? || patch.nil?

  "#{major}.#{minor}.#{patch}"
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
  puts "using C++ standard... #{$CPP_STANDARD}"
  puts "g++ reports version... " + `#{CONFIG['CXX']} --version|head -n 1|cut -f 3 -d " "`
end

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

have_library('lapack_atlas', 'clapack_dgetrf')

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


# For release, these next two should both be changed to -O3.
$CFLAGS += " -O3" #" -O0 -g "
#$CFLAGS += " -static -O0 -g "
$CPPFLAGS += " -O3 -std=#{$CPP_STANDARD}" #" -O0 -g -std=#{$CPP_STANDARD} " #-fmax-errors=10 -save-temps
#$CPPFLAGS += " -static -O0 -g -std=#{$CPP_STANDARD} "

CONFIG['warnflags'].gsub!('-Wshorten-64-to-32', '') # doesn't work except in Mac-patched gcc (4.2)
CONFIG['warnflags'].gsub!('-Wdeclaration-after-statement', '')
CONFIG['warnflags'].gsub!('-Wimplicit-function-declaration', '')

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
