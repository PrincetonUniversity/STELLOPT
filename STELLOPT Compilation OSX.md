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
    #  gcc7 : Appears to be working
    #  gcc8 : Problems with PetSc so no GENE/SFINCS
    #  gcc9 : Unable to confirm
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
    # The following is for Python support
    sudo port install python37 +gccX
    sudo port select --set python python27
    sudo port select --set python3 python37
    sudo port install py37-pyqt4 +gccX     (Note there is a weird issue with dbus requireing the port to be force installed)
    sudo port install py37-matplotlib +dvipng +pyside +qt4 +gccX
    sudo port install py37-h5py +gcc7 +openmpi
    # The following are needed for GENE
    sudo port install petsc +gccX +openmpi +openblas -accelerate +metis +mumps +parmetis +suitesparse +superlu_dist +complex
    sudo port install slepc +gccX +openmpi +openblas +arpack -accelerate
    # NCARG has in issue with ESMF, but you don't need it for STELLOPT just XGTOVMI
    # I placed it here since it builds upon GSL and FFTW-3
    #sudo port install ncarg +gccX +openmpi

    <===gcc8 version===>

    Tested on: 
    MacBook Pro (Retina, 15-inch, Mid 2015), 2.5 GHz Intel Core i7, 16 GB 1600 MHz DDR3
    AMD Radeon R9 M370X 2048 MG Intel Iris Pro 1536 MB
    macOS Mojave, Version 10.14.2
    ** SEE NOTE ON MOJAVE

    MacPorts installation version: MacPorts-2.5.4-10.14-Mojave 17-35-43-409.pkg

    sudo port install gcc8
    sudo port install openmpi
    sudo port install hdf5 +fortran +gcc8 +hl +openmpi
      (also asks for openmpi-gcc8. say yes.)
     sudo port select --set mpi openmpi-gcc8-fortran
    sudo port install netcdf +gcc8 +openmpi 
    sudo port install netcdf-fortran +gcc8 +openmpi
    sudo port install pgplot 
       (C bindings broken in gcc8)
    sudo port install OpenBLAS +gcc8 
    sudo port install scalapack +gcc8 +openmpi +openblas

    # NOTE!!!!
    # COILOPT++, Python support, and GENE support packaged not installed/tested. 
    # To Do:   Upate for gcc8

    # NOTE2:  the soft-link for mpifort is broken here.  Redirect it to mpifort-openmpi-gcc8 using
    # sudo rm /opt/local/bin/mpifort
    # sudo ln -sf /opt/local/bin/mpifort-openmpi-gcc8 /opt/local/bin/mpifort
    # 
    # NOTE3: if you run into ld x86-64 architecture linking errors, remove all object / mod files (*.o/*.mod) and recompile.

    # MOJAVE SUPPORT: 
    OSX MOJAVE does not support universal build architectures (64-bit only). Change the line in /opt/local/macports.conf to x86_64 or bust.  Also to avoid x86_64 linking errors, you must remove all precompiled *.o files before building clean. 

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
