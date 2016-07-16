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
  # Workaround for Windows having a different name for null device
  null_device = (system 'echo "Hi" </dev/null') ? '/dev/null' : 'nul'

  cxxvar = proc { |n| `#{CONFIG['CXX']} -E -dM - <#{null_device} | grep #{n}`.chomp.split(' ')[2] }
  major = cxxvar.call('__GNUC__')
  minor = cxxvar.call('__GNUC_MINOR__')
  patch = cxxvar.call('__GNUC_PATCHLEVEL__')

  raise("unable to determine g++ version (match to get version was nil)") if major.nil? || minor.nil? || patch.nil?

  "#{major}.#{minor}.#{patch}"
end


if CONFIG['CXX'] == 'clang++'
  $CXX_STANDARD = 'c++11'

else
  version = gplusplus_version
  if version < '4.3.0' && CONFIG['CXX'] == 'g++'  # see if we can find a newer G++, unless it's been overridden by user
    if !find_newer_gplusplus
      raise("You need a version of g++ which supports -std=c++0x or -std=c++11. If you're on a Mac and using Homebrew, we recommend using mac-brew-gcc.sh to install a more recent g++.")
    end
    version = gplusplus_version
  end

  if version < '4.7.0'
    $CXX_STANDARD = 'c++0x'
  else
    $CXX_STANDARD = 'c++11'
  end
  puts "using C++ standard... #{$CXX_STANDARD}"
  puts "g++ reports version... " + `#{CONFIG['CXX']} --version|head -n 1|cut -f 3 -d " "`
end

#$libs += " -lprofiler "

# For release, these next two should both be changed to -O3.
$CFLAGS += " -O3 "
#$CFLAGS += " -static -O0 -g "
$CXXFLAGS += " -O3 -std=#{$CXX_STANDARD} " #-fmax-errors=10 -save-temps
#$CXXFLAGS += " -static -O0 -g -std=#{$CXX_STANDARD} "

CONFIG['warnflags'].gsub!('-Wshorten-64-to-32', '') # doesn't work except in Mac-patched gcc (4.2)
CONFIG['warnflags'].gsub!('-Wdeclaration-after-statement', '')
CONFIG['warnflags'].gsub!('-Wimplicit-function-declaration', '')

have_func("rb_array_const_ptr", "ruby.h")
have_macro("FIX_CONST_VALUE_PTR", "ruby.h")
have_macro("RARRAY_CONST_PTR", "ruby.h")
have_macro("RARRAY_AREF", "ruby.h")

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
