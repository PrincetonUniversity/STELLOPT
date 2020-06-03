STELLOPT Compilation at the MPCDF (RZG-MPG Machines)
====================================================

This page details how to compile the STELLOPT family of codes on Hydra
at [Rechenzentrum Garching](@http://www.rzg.mpg.de/) (RZG-MPE). In order
to do so you will need an account on their system. These build
instructions are for the Intel based compilers

Draco
-----

    module load git
    module load intel
    module load mkl
    module load impi
    module load netcdf-mpi
    module load hdf5-mpi
    module load fftw
    module load petsc-cplx
    module load slepc-cplx
    module load nag_flib/intel/mk24

Cobra
-----

    module load git
    module load intel
    module load mkl
    module load impi
    module load hdf5-mpi/1.10.5
    module load netcdf-mpi/4.7.0
    module load fftw-mpi
    module load petsc-complex
    module load slepc-complex

General Notes
-------------

The Cobra and Draco clusters are automatically detected
by the make.inc script.
