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

# Check to see whether the installer is zypper.
# If so, then install necessary packages.
if [ -x /usr/bin/zypper ] ; then
    source ~/.profile
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
    source ~/.bashrc
    echo ""
    echo "Login as root, now..."
    sudo yum install gcc gcc-c++ curl cpupower
    sudo cpupower frequency-set -g performance
    echo
    echo "Going to install the following packages:"
    echo "libatlas3 libatlas3-devel blas-devel"
    sudo yum install libatlas3-devel blas-devel ruby-devel lapack-devel
fi

# Check to see whether the installer is apt-get
# If so, then install necessary packages.
if [ -f /usr/bin/apt-get ] ; then
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
    echo "Setting CPLUS_INCLUDE_PATH..."
    echo "export CPLUS_INCLUDE_PATH=/usr/include/atlas" >> ~/.profile
    echo "Setting C_INCLUDE_PATH..."
    echo "export C_INCLUDE_PATH=/usr/include/atlas" >> ~/.profile
    echo "source ~/.profile" >> ~/.bashrc
fi

git clone https://github.com/SciRuby/nmatrix.git

# If gem has a Gemfile then it will be installed
if [ -f  nmatrix/Gemfile ] ; then
    cd nmatrix
    bundle install
    bundle exec rake compile -- --with-lapacklib
    bundle exec rake spec
fi
