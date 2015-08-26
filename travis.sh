#!/bin/bash

set -ev #fail at the first command that returns non-zero exit value

if [ "$TRAVIS_OS_NAME" = "osx" ]
then
  export INSTALL_COMMAND="brew install -y"
else
  export INSTALL_COMMAND="sudo apt-get install -y"
fi

if [ "$1" = "before_install" ]
then
  gem install bundler -v '~> 1.6'
  if [ "$TRAVIS_OS_NAME" = "osx" ]
  then
    brew update
  else
    sudo apt-get update -qq
  fi

  if [ -n "$USE_ATLAS" ]
  then
      $INSTALL_COMMAND libatlas-base-dev
  fi

  # travis-ci runs on Ubuntu 12.04, where the openblas package doesn't
  # provide a liblapack.so, so we test using the blas from openblas
  # and the reference lapack implementation. Not entirely sure if
  # this will work.
  if [ -n "$USE_OPENBLAS" ]
  then
    $INSTALL_COMMAND libopenblas-dev
    # Since we install libopenblas first, liblapack won't try to install
    # libblas (the reference BLAS implementation).
    $INSTALL_COMMAND liblapack-dev
  fi

  if [ -n "$USE_REF" ]
  then
    $INSTALL_COMMAND liblapack-dev
  fi
fi

if [ "$1" = "script" ]
then
  if [ -n "$USE_ATLAS" ]
  then
    # Need to put these commands on separate lines (rather than use &&)
    # so that bash set -e will work.
    bundle exec rake compile nmatrix_plugins=all
    bundle exec rake spec nmatrix_plugins=all
  fi

  if [ -n "$USE_OPENBLAS" ]
  then
    bundle exec rake compile nmatrix_plugins=lapacke
    bundle exec rake spec nmatrix_plugins=lapacke
  fi

  if [ -n "$USE_REF" ]
  then
    bundle exec rake compile nmatrix_plugins=lapacke
    bundle exec rake spec nmatrix_plugins=lapacke
  fi

  if [ -n "$NO_EXTERNAL_LIB" ]
  then
    bundle exec rake compile
    bundle exec rake spec
  fi
fi
