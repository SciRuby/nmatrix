#!/bin/bash

set -ev #fail at the first command that returns non-zero exit value

# Use rbenv on OSX iff ruby_version is given
if [ -n "$ruby_version" -a "$TRAVIS_OS_NAME" = "osx" ]; then
  export PATH="$HOME/.rbenv/bin:$PATH"
  if [ -x $HOME/.rbenv/bin/rbenv ]; then
    eval "$(rbenv init -)"
  fi
  export RBENV_VERSION=$ruby_version
  unset GEM_PATH GEM_HOME
fi

if [ "$1" = "install" ]
then
  bundle install --jobs=3 --retry=3 --path=${BUNDLE_PATH:-vendor/bundle}
fi

if [ "$1" = "before_install" ]
then
  case "$TRAVIS_OS_NAME" in
    linux)
      sudo apt-get update -qq
      ;;
    osx)
      brew update >/dev/null
      ;;
  esac

  # Installing ruby by using rbenv on OSX iff ruby_version is given
  if [ -n "$ruby_version" -a "$TRAVIS_OS_NAME" = "osx" ]; then
    git clone https://github.com/rbenv/rbenv.git ~/.rbenv
    git clone https://github.com/rbenv/ruby-build.git ~/.rbenv/plugins/ruby-build

    eval "$(rbenv init -)"

    # Install ruby
    (
      brew install bison openssl readline libyaml
      brew link --force openssl
      RBENV_VERSION=system
      MAKEOPTS='-j 4'
      CONFIGURE_OPTS="--disable-install-doc --with-out-ext=tk,tk/tkutil --with-opt-dir=/usr/local"
      rbenv install $ruby_version
    )

    gem pristine --all
    gem update --no-document --system
    gem update --no-document
  fi

  gem install --no-document bundler -v '~> 1.6'

  if [ -n "$USE_ATLAS" ]
  then
    case "$TRAVIS_OS_NAME" in
      linux)
        sudo apt-get install -y libatlas-base-dev
        ;;
      osx)
        echo "FIXME: ATLAS on OSX environment is not supported, currently" >2
        exit 1
        ;;
    esac
  fi

  # travis-ci runs on Ubuntu 12.04, where the openblas package doesn't
  # provide a liblapack.so, so we test using the blas from openblas
  # and the reference lapack implementation. Not entirely sure if
  # this will work.
  if [ -n "$USE_OPENBLAS" ]
  then
    case "$TRAVIS_OS_NAME" in
      linux)
        sudo apt-get install -y libopenblas-dev
        # Since we install libopenblas first, liblapack won't try to install
        # libblas (the reference BLAS implementation).
        sudo apt-get install -y liblapack-dev
        ;;
      osx)
        brew install homebrew/science/openblas
        ;;
    esac
  fi

  if [ -n "$USE_REF" ]
  then
    case "$TRAVIS_OS_NAME" in
      linux)
        sudo apt-get install -y liblapack-dev
        ;;
      osx)
        brew install homebrew/dupes/lapack
        ;;
    esac
  fi
fi

if [ "$1" = "script" ]
then
  nmatrix_plugins_opt=''

  if [ -n "$USE_ATLAS" ]
  then
    # Need to put these commands on separate lines (rather than use &&)
    # so that bash set -e will work.
    nmatrix_plugins_opt='nmatrix_plugins=atlas'
  fi

  if [ -n "$USE_OPENBLAS" ]
  then
    nmatrix_plugins_opt='nmatrix_plugins=lapacke'
  fi

  if [ -n "$USE_REF" ]
  then
    nmatrix_plugins_opt='nmatrix_plugins=lapacke'
  fi

  if [ -n "$NO_EXTERNAL_LIB" ]
  then
    nmatrix_plugins_opt=''
  fi

  bundle exec rake travis:env

  bundle exec rake compile $nmatrix_plugins_opt || {
    echo === Contents of mkmf.log ===
    cat tmp/*/nmatrix/*/mkmf.log
    exit 1
  }
  bundle exec rake spec $nmatrix_plugins_opt
fi
