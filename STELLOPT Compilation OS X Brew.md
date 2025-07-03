STELLOPT Compilation on OS X (Brew)
============================

This page details how to compile the STELLOPT family of codes on OS X.
In order to do so you will need to install
[GNU based compilers](http://gcc.gnu.org/) on your Apple machine and a
package manager (such as Homebrew).

------------------------------------------------------------------------

Installation Steps

1\. Download and install [Homebrew](https://brew.sh/).

3\. Using [Homebrew](https://brew.sh/), install the following
packages taking note of the caveates.

    <===gccXX version===>
    #  gcc10 : Working
    #  gcc11 : Issue with hdf5 23.02.2024
    #  gcc12 : Issue with hdf5 23.02.2024
    #  gcc13 : Working as of 14.11.2024
    #  gcc14 : Libraries still missing (14.11.2024)
    brew install git
    brew install gcc
    brew install openmpi
    brew install netcdf
    brew install netcdf-fortran
    brew install fftw
    brew install openblas
    brew install scalapack
    brew install python3
    # All other python modules should be handled durring 
    # install using pip and setuptools

4\. Installing HDF5 with high level language support.

    brew edit hdf5

Edit the `args` portion of the file to include `-DHDF5_BUILD_HL:BOOL=ON`, then save and close the `hdf5.rb` file.

    brew reinstall --build-from-source hdf5

5\. (optional) The following are optional for compiling other codes.

    ############ TRAVIS ############
    brew install netcdf-cxx
    ############ SFINCS ############
    brew install superlu
    brew install petsc

6\. Now pull stellopt with the command 

    git clone git@github.com:PrincetonUniversity/STELLOPT.git


7\. Set the environement variable (note the first one STELLOPT_PATH should be put in your .zshenv file, assuming you use the default ZSH)

    export STELLOPT_PATH=<path to your repository>
    export MACHINE=osx_brew
    build_all

This will build all codes starting with LIBSTELL.  If you see errors durring compilation of LIBSTELL please run `build_all >& stellopt_build.log` and send the full stellopt_build.log file to a developer for debugging.
