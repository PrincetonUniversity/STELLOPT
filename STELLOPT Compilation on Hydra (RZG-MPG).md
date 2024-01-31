STELLOPT Compilation at the MPCDF (RZG-MPG Machines)
====================================================

This page details how to compile the STELLOPT family of codes
at [Max-Plank Computational Data Facility (formerly Rechenzentrum Garching)](@http://www.rzg.mpg.de/) (MPCDF). In order to do so you will need an account on their system. These build
instructions are for the Intel based compilers


Cobra
-----

    module load git
    module load intel/2023.1.0.x
    module load mkl/2023.1
    module load impi/2021.9
    module load hdf5-mpi/1.14.1
    module load netcdf-mpi/4.8.1
    module load fftw-mpi

Raven
-----

    module load git
    module load intel/21.4.0
    module load mkl/2021.2
    module load impi/2021.4
    module load netcdf-mpi/4.8.1
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
