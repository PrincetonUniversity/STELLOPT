STELLOPT Compilation at the MPCDF (RZG-MPG Machines)
====================================================

This page details how to compile the STELLOPT family of codes
at [Max-Plank Computational Data Facility (formerly Rechenzentrum Garching)](@http://www.rzg.mpg.de/) (MPCDF). In order to do so you will need an account on their system. These build
instructions are for the Intel based compilers

Draco
-----

    module load git
    module load intel
    module load mkl
    module load impi
    module load netcdf-mpi
    module load hdf5-mpi
    module load fftw-mpi

Cobra
-----

    module load git
    module load intel
    module load mkl
    module load impi
    module load hdf5-mpi/1.10.5
    module load netcdf-mpi/4.7.0
    module load fftw-mpi

Raven
-----

    module load git
    module load intel/19.1.2
    module load mkl
    module load impi/2019.8
    module load netcdf-mpi
    module load hdf5-mpi
    module load fftw-mpi
    module load anaconda/3/2020.02

IPP-HGW Theory (clus47)
-----

    module load intel
    module load mkl
    module load netcdf
    module load hdf5
    module load openmpi
    module load petsc
    module load slepc
    module load nag_flib/intel/mk25

General Notes
-------------

The Cobra and Draco clusters are automatically detected
by the make.inc script. The IPP-HGW theory cluster
requires the user to set the environment variable 
MACHINE equal to 'theoryhgw'.

For the benchmark tests you'll need python3
packages.  It's best to install a copy of miniconda
locally on your account and setup your packages
from there (Draco & Cobra). [miniconda downloads](https://docs.conda.io/en/latest/miniconda.html)
