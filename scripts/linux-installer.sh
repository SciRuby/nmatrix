#!/bin/bash
#
#     Linux-Installer.sh
#
#     An install script for Fedora and OpenSUSE Systems.and Debian.
#     Installs to home directory
#
#     Copyright (c) 2014
#
#      Author:
#              Magdalen Berns
#     Contact:
#              <m.berns@thismagpie.com>

bundle_installer()
    if [ -f  Gemfile ] ; then
        cd nmatrix
        bundle install
        bundle exec rake compile -- --with-lapacklib
        bundle exec rake spec
    fi

if [ -x /usr/bin/zypper ] ; then
    echo ""
    echo "Going to install the following packages:"
    echo "gcc gcc-c++ curl and cpupower."
    echo "Login as root, now..."
    sudo zypper in -y gcc gcc-c++ curl cpupower
    sudo cpupower frequency-set -g performance
    echo
    echo "Going to install the following packages:"
    echo "libatlas3 libatlas3-devel blas-devel"
    sudo zypper in -y libatlas3-devel blas-devel ruby-devel lapack-devel
elif [ -f /usr/bin/yum ] ; then
    echo ""
    echo "Login as root, now..."
    sudo yum install gcc gcc-c++ curl cpupower
    sudo cpupower frequency-set -g performance
    echo
    echo "Going to install the following packages:"
    echo "libatlas3 libatlas3-devel blas-devel"
    sudo yum install libatlas3-devel blas-devel ruby-devel lapack-devel
elif [ -f /usr/bin/apt-get ] ; then
    echo ""
    source ~/.bash_profile
    echo "Login as root, now..."
    sudo apt-get install gcc gcc-c++ curl cpupower
    sudo cpupower frequency-set -g performance
    echo
    echo "Going to install the following packages:"
    echo "libatlas3 and libatlas3-dev"
    sudo apt-get install libatlas3-dev ruby-dev lapack-dev
fi

echo "Installation for $s complete!"
echo "Setting PATHS"

if [ -d /usr/include/atlas ] ; then
  echo "Setting CPLUS_INCLUDE_PATH in ~/.bashrc"
  echo "export CPLUS_INCLUDE_PATH=/usr/include/atlas" >> ~/.bashrc
  echo "Setting C_INCLUDE_PATH in ~/.bashrc"
  echo "export C_INCLUDE_PATH=/usr/include/atlas" >> ~/.bashrc
  echo "source ~/.bashrc"
fi

if [ -d ../../nmatrix ]; then
  echo "NMatrix found"
  cd ..
  if [ -f  Gemfile ] ; then
    bundle_installer
  fi
else
  echo "This script was not run from the NMatrix scripts directory."
  while true; do
      read -p "Do you want to install NMatrix here and now? y or n and press [ENTER]:" ans
      case $ans in
        [Yy]* ) if [ ! -d $PWD/nmatrix ] ; then
                      echo ""
                      echo "Cloning NMatrix in $PWD"
                      git clone https://github.com/sciruby/nmatrix
                      bundle_installer
                else
                      echo "Directory Already exists"
                      cd nmatrix
                      bundle_installer
                fi
                break;;
        [Nn]* ) echo ""
                echo "Sorry it is not possible to install NMatrix unless you are willing to clone it first!"; exit;;
      esac
  done
fi
