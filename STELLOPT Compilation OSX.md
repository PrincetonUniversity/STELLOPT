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

    <===OSX Version===>
    10.14 Mojave Something appears to be broken in scalapack build (GCC7/GCC8/GCC9)
    <===gccX version===>
    #  gcc7  : Appears to be working
    #  gcc8  : Problems with PetSc so no GENE/SFINCS
    #  gcc9  : Problems with PetSc so no GENE/SFINCS
    #  gcc10 : Not tested 
    sudo port install gccX                   (gccX should be gcc7 or gcc8 or another gcc variant)
    sudo port install openmpi-gccX +fortran
    sudo port select --set mpi openmpi-gccX-fortran
    # Soft link for mpifort is broken for GCC8 as of Feb 20, 2019 ignore next line otherwise
    sudo ln -sf /opt/local/bin/mpifort-openmpi-gcc8 /opt/local/bin/mpifort
    sudo port install hdf5 +fortran +gccX +hl +openmpi
    sudo port install netcdf-fortran +gccX +openmpi
    sudo port install fftw-3 +gccX +openmpi
    sudo port install pgplot +gccX
    sudo port install OpenBLAS +gccX +lapack +native
    sudo port install scalapack +gccX +openmpi +openblas
    # The following are needed for TRAVIS
    sudo port install netcdf-cxx4 +gccX +openmpi
    # The following are needed for COILOPT++
    sudo port install silo +gccX -hdf5
    sudo port install gsl +gccX +optimize
    # The following is for Python support (please take note of your python version)
    sudo port install python38 +gccX
    sudo port select --set python python27
    sudo port select --set python3 python37
    sudo port install py38-pyqt4 +gccX     (Note there is a weird issue with dbus requireing the port to be force installed)
    sudo port install py38-matplotlib +dvipng +pyside +qt4 +gccX
    sudo port install py38-h5py +gccX +openmpi
    # The following are needed for GENE
    #sudo port install petsc +gccX +openmpi +openblas -accelerate +metis +mumps +parmetis +suitesparse +superlu_dist +complex
    #sudo port install slepc +gccX +openmpi +openblas +arpack -accelerate
    # NCARG has in issue with ESMF, but you don't need it for STELLOPT just XGTOVMI
    # I placed it here since it builds upon GSL and FFTW-3
    #sudo port install ncarg +gccX +openmpi

4\. Now the NCARG package broken as of this writing but there\'s a
workaround. Download the tarball from the
[NCAR website](http://www.ncarg.ucar.edu/) the place the tarball in the
appropriate location
([see ticket \#43615](https://trac.macports.org/ticket/42541)). Then
rerun the port install command for ncarg.

5\. Now pull stellopt with the command 

    git clone git@github.com:PrincetonUniversity/STELLOPT.git


6\. Set the environement variable

    export MACHINE=macports
    build_all

This will build all codes starting with LIBSTELL.  If you see errors durring compilation of LIBSTELL please run `build_all >& stellopt_build.log` and send the full stellopt_build.log file to a developer for debugging.
