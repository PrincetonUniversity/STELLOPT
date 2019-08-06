STELLOPT Compilation on Ubuntu
==============================

This page details how to compile the STELLOPT family of codes on
[Ubuntu](@http://www.ubuntu.com/). In order to do so you will need to
install [GNU based compilers](@http://gcc.gnu.org/) on your Linux
machine. In principle this set of steps should work for any Linux based
distribution which uses [Debian](@https://www.debian.org/) based package
management.

------------------------------------------------------------------------

Installation Steps

1.  Verify that your version of Ubuntu is up to date.
2.  Use MacPorts to install your packages. For this example we\'re using
    gfortran 4.8 but you may substitute any other version.

<!-- -->

    sudo apt-get install gfortran-4.8
    sudo apt-get install openmpi-bin
    sudo apt-get install libopenmpi-dev
    sudo apt-get install gfortran
    sudo apt-get install g++
    sudo apt-get install libnetcdf-dev
    sudo apt-get install libhdf5-serial-dev
    sudo apt-get install pgplot5
    sudo apt-get install libblas-dev
    sudo apt-get install liblapack-dev
    sudo apt-get install libncarg-dev

1.  Now create a directory in which you wish to compile STELLOPT and
    place the STELLOPT ZIP file in that directory. Then unzip it.
2.  Unzip stellinstall
3.  Now you should be able to run the setup script and compile the
    STELLOPT family of codes.
