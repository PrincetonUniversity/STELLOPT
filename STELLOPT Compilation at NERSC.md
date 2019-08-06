STELLOPT Compilation at NERSC
=============================

This page details how to compile the STELLOPT family of codes on
machines at [NERSC](@http://www.nersc.gov/). Note that for the BEAMS3D
code you will need to compile ADAS separately.

------------------------------------------------------------------------

**Edison**

1.  Load the necessary module files

<!-- -->

    module load cray-hdf5-parallel
    module load cray-netcdf-hdf5parallel
    module load fftw/3.3.4.0
    module load cray-petsc-complex
    module load slepc-complex
    module load cray-tpsl
    module load gsl
    module load silo

1.  Now create a directory in which you wish to compile STELLOPT and
    place the STELLOPT ZIP file in that directory. Then unzip it.
2.  Unzip stellinstall
3.  Now you should be able to run the setup script and compile the
    STELLOPT family of codes.

** Hopper **

1.  Load the necessary modules files

<!-- -->

    module swap PrgEnv-pgi PrgEnv-intel
    module load cray-hdf5-parallel
    module load cray-netcdf-hdf5parallel
    module load ncar
    module load fftw/3.3.4.0
    module load cray-petsc-complex
    module load slepc-complex
    module load cray-tpsl
    module load gsl
    module load silo

1.  Now create a directory in which you wish to compile STELLOPT and
    place the STELLOPT ZIP file in that directory. Then unzip it.
2.  Unzip stellinstall
3.  Now you should be able to run the setup script and compile the
    STELLOPT family of codes.

**Cori**

1.  Load the necessary module files (Intel is the default compiler)

<!-- -->

    module load cray-hdf5-parallel
    module load cray-netcdf-hdf5parallel
    module load fftw/3.3.4.5
    module load cray-petsc-complex
    module load slepc-complex
    module load cray-tpsl
    module load gsl
    module load silo

1.  Now create a directory in which you wish to compile STELLOPT and
    place the STELLOPT ZIP file in that directory. Then unzip it.
2.  Unzip stellinstall
3.  Now you should be able to run the setup script and compile the
    STELLOPT family of codes. \-\--
    -   Compiling NTCC on Cori** gmake FORTRAN\_VARIANT=Intel
        INTEL\_CC=Y USEFC=Y \_64=Y GCC\_VERSION=492
        TERMCAP=/usr/lib64/libtermcap.so.2.0.8**
