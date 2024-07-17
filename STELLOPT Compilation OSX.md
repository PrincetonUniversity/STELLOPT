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

3\. Use MacPorts to install your packages. Note you should probably used
gcc7 or gcc8 and replace the gccX\'s with the consistent version. I\'ve
also chosen OpenMPI but the sources could also be built with MPICH.
Please note that as of GCC10 the flag `-fallow-argument-mismatch` must be
added at compile time.

    <===gccX version===>
    #  gcc10 : Working
    #  gcc11 : Issue with hdf5 02.23.2024
    #  gcc12 : Issue with hdf5 02.23.2024
    #  gcc13 : Issue with openmpi-gcc12 02.23.2024
    sudo port install gcc10                   (gccX should be gcc7 or gcc8 or another gcc variant)
    sudo port install gcc10-libcxx (for example is gcc10 above then gcc10-libcxx)
    sudo port install openmpi-gcc10 +fortran
    sudo port select --set mpi openmpi-gcc10-fortran
    sudo port select --set gcc mp-gcc10
    sudo port install hdf5 +fortran +gcc10 +hl +openmpi
    sudo port install pnetcdf +gcc10 +openmpi -mpich
    sudo port install netcdf +gcc10
    sudo port install netcdf-fortran +gcc10
    sudo port install fftw-3 +gcc10 +openmpi
    sudo port install OpenBLAS +gcc10 +lapack +native
    sudo port install scalapack +gcc10 +openmpi +openblas
    sudo port install git
    # The following are needed for TRAVIS
    sudo port install netcdf-cxx4 +gcc10
    # The following are needed for COILOPT++
    sudo port install silo +gcc10 -hdf5
    sudo port install gsl +gccX +optimize
    # The following is for Python support (please take note of your python version)
    sudo port install pythonYY +gcc10
    sudo port select --set python pythonYY
    sudo port select --set python3 pythonYY
    sudo port install pyYY-pip
    sudo port select --set pip pipYY
    sudo port select --set pip3 pipYY
    # It is suggested to use pip3 to install the rest of your python packages
    pip3 install scipy h5py pyqt4 matplotlib
    # Other codes may need PETSC (SFINCS) and SLEPC (GENE) macports are broken but you can build yourself
    # PETSC should needs to be built from source
    # ./configure --download-sowing --with-openblas=1 --with-openblas-dir=/opt/local/
    # make PETSC_DIR=<path_to_git_dir>/petsc PETSC_ARCH=arch-darwin-c-debug all
    # export PETSC_DIR=<path_to_git_dir>/petsc
    # export PETSC_ARCH=arch-darwin-c-debug

4\. Now pull stellopt with the command 

    git clone git@github.com:PrincetonUniversity/STELLOPT.git


5\. Set the environement variable

    export MACHINE=macports
    build_all

This will build all codes starting with LIBSTELL.  If you see errors durring compilation of LIBSTELL please run `build_all >& stellopt_build.log` and send the full stellopt_build.log file to a developer for debugging.
