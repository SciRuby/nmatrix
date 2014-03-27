#!/bin/bash

brew-path() { brew info $1 | head -n3 | tail -n1 | cut -d' ' -f1; }

# Try using the following for GCC 4.9:
#
#   brew install gmp4 mpfr2 libmpc08 isl011 cloog018
#
#

brew install gcc49 --enable-fortran
# Source for this is: http://stackoverflow.com/questions/19535422/os-x-10-9-gcc-links-to-clang


# You may wish to re-install your Ruby if you're using rbenv. To do
# so, make sure you've installed openssl, readline, and libyaml.
#
# The commands for this are:
#
#    CC=gcc-4.8 RUBY_CONFIGURE_OPTS="--with-openssl-dir=`brew --prefix openssl` --with-readline-dir=`brew --prefix readline` --with-gcc=gcc-4.8 --enable-shared" rbenv install --keep 2.0.0-p247
#
#