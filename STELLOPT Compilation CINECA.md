STELLOPT Compilation at CINECA
====================================================

This page details how to compile the STELLOPT family of codes
at [CINECA](@http://www.hpc.cineca.it/). In order
to do so you will need an account on their system. These build
instructions are for the Intel based compilers

Marconi
-----

    module load git
    module load intel/pe-xe-2018--binary
    module load intelmpi/2018--binary
    module load mkl/2018--binary
    module load zlib/1.2.8--gnu--6.1.0
    module load szip/2.1--gnu--6.1.0
    module load hdf5/1.10.4--intelmpi--2018--binary
    module load netcdf/4.6.1--intelmpi--2018--binary
    module load netcdff/4.4.4--intelmpi--2018--binary
    module load netcdf-cxx4/4.3.0--intelmpi--2018--binary
    module load fftw/3.3.7--intelmpi--2018--binary 
    module load petsc/3.13.3--intelmpi--2018--binary

General Notes
-------------

For the make.inc script to detect Marconi you will need
to set `MACHINE='marconi'`.
