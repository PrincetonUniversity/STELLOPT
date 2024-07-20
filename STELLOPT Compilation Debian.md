STELLOPT Compilation on Ubuntu
==============================

This page details how to compile the STELLOPT family of codes on
[Ubuntu](@http://www.ubuntu.com/). In order to do so you will need to
install [GNU based compilers](@http://gcc.gnu.org/) on your Linux
machine. In principle this set of steps should work for any Linux based
distribution which uses [Debian](@https://www.debian.org/) based package
management.

Installation Steps
-----

1\. Verify that your version of Ubuntu is up to date.
2\. Use apt-get to install your packages. 

    sudo apt-get install git
    sudo apt-get install gcc
    sudo apt-get install gfortran
    sudo apt-get install openmpi-bin
    sudo apt-get install libopenmpi-dev
    sudo apt-get install libopenblas-dev
    sudo apt-get install libscalapack-openmpi-dev
    sudo apt-get install libnetcdf-dev
    sudo apt-get install libnetcdff-dev
    sudo apt-get install libhdf5-openmpi-dev
    sudo apt-get install hdf5-tools
    sudo apt-get install python3 pip

3\. Setup your environment variables

    export MACHINE="ubuntu"
    export STELLOPT_PATH=<path to repo directory>


General Notes
-------------
Some of the package names may change over time.

