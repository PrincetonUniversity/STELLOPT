STELLOPT Compilation on OS X
============================

This page details how to compile the STELLOPT family of codes on OS X.
In order to do so you will need to install
[GNU based compilers](http://gcc.gnu.org/) on your Apple machine and a
package manager (such as MacPorts).

------------------------------------------------------------------------

Installation Steps

1\. Download and install [X Code](https://developer.apple.com/xcode/) for
your version of OS X.

2\. Download and install [MacPorts](https://www.macports.org/). Please
note we\'re assuming you\'ll be using the default location for
installation of the various ports. If not you\'ll need to change the
path to the variable MACPORTS in the setup script (and rezip).

3\. Using [MacPorts](https://www.macports.org/), install the following
ports taking note of the caveates. Places where you see gccXX use whatever
version of GCC you built for. Also, there is an option to use OpenBLAS
or the built-in Accellerate Framework. Note that in going from OS10.15 
to OS26.0 some warning appear.

    <===gccXX version===>
    #  gcc10 : Working
    #  gcc11 : Issue with hdf5 23.02.2024
    #  gcc12 : Issue with hdf5 23.02.2024
    #  gcc13 : Working as of 14.11.2024
    #  gcc14 : Working as of 18.09.2025
    #  gcc15 : No OpenMPI as of 18.09.2025
    export MACOSX_DEPLOYMENT_TARGET=26.0 # Needed at the 10.15->26.0 phase
    sudo port install git
    sudo port install gccXX        (gccXX should be gcc10 or another gcc variant)
    sudo port install openmpi-gccXX +fortran
    # The next two are needed for hdf5 and netcdf to build properly
    sudo port install llvm-18
    sudo port install clang-18 +analyzer +libstdcxx
    sudo port select --set mpi openmpi-gccXX-fortran
    sudo port select --set gcc mp-gccXX
    sudo port select --set clang mp-clang-17
    sudo port select --set llvm mp-llvm-18
    # Because of dependencies if HDF5 can't be built with +openmpi neither can netcdf
    sudo port install hdf5 +fortran +hl +gccXX (+openmpi issue with libc++)
    sudo port install netcdf (+openmpi issue with libc++)
    sudo port install netcdf-fortran +gccXX (+openmpi issue with libc++)
    sudo port install fftw-3 +openmpi +gccXX
    # This section for using OpenBLAS
    sudo port install OpenBLAS +lapack +native +gccXX
    sudo port install scalapack +openmpi +openblas +gccXX
    # This section for using the Accelerate Framework
    sudo port install scalapack +openmpi +accelerate +gccXX
    <===pythonYYY version===>
    sudo port install pythonYYY
    sudo port select --set python pythonYYY
    sudo port select --set python3 pythonYYY
    sudo port install pyYYY-pip
    sudo port select --set pip pipYYY
    sudo port select --set pip3 pipYYY
    # All other python modules should be handled durring 
    # install using pip and setuptools

4\. (optional) The following are optional for compiling other codes.

    ############ TRAVIS ############
    sudo port install netcdf-cxx4 (+openmpi issue with libc++)
    ############ COILOPT++ ############
    sudo port install silo +gccXX -hdf5
    sudo port install gsl +gccXX +optimize
    ############ SFINCS ############
    sudo port install mumps +accelerate +openmpi +gccXX
    #sudo port install superlu +accelerate +gccXX
    sudo port install superlu_dist +accelerate +openmpi +gccXX
    sudo port install petsc +accelerate +openmpi +mumps +superlu_dist +fftw +gccXX (+hdf5 require parallel)
    ############ MANGO ##############
    sudo port install gsl

5\. Now pull stellopt with the command 

    git clone git@github.com:PrincetonUniversity/STELLOPT.git


6\. Set the environement variable (note the first one STELLOPT_PATH should be put in your .zshenv file, assuming you use the default ZSH)

    export STELLOPT_PATH=<path to your repository>
    export MACHINE=macports
    build_all

This will build all codes starting with LIBSTELL.  If you see errors durring compilation of LIBSTELL please run `build_all >& stellopt_build.log` and send the full stellopt_build.log file to a developer for debugging.
