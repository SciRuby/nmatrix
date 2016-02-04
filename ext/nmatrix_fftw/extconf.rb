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

require 'mkmf'

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

fftw_libdir = RbConfig::CONFIG['libdir']
fftw_incdir = RbConfig::CONFIG['includedir']
fftw_srcdir = RbConfig::CONFIG['srcdir']

$CFLAGS = ["-Wall -Werror=return-type -I$(srcdir)/../nmatrix -I$(srcdir)/lapacke/include",$CFLAGS].join(" ")
$CXXFLAGS = ["-Wall -Werror=return-type -I$(srcdir)/../nmatrix -I$(srcdir)/lapacke/include -std=c++11",$CXXFLAGS].join(" ")
$CPPFLAGS = ["-Wall -Werror=return-type -I$(srcdir)/../nmatrix -I$(srcdir)/lapacke/include -std=c++11",$CPPFLAGS].join(" ")

CONFIG['CXX'] = 'g++'

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
