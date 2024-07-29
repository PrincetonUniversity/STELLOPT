STELLOPT Compilation at the MPCDF
====================================================

This page details how to compile the STELLOPT family of codes
at [Max-Plank Computational Data Facility (formerly Rechenzentrum Garching)](@https://www.mpcdf.mpg.de/) (MPCDF). In order to do so you will need an account on their system. These build
instructions are for the Intel based compilers


Cobra (Defunct)
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
    module load anaconda/3/2021.11

Viper
-----

    module load git/2.43
    module load intel/2024.0
    module load mkl/2024.0
    module load impi/2021.11
    module load netcdf-mpi/4.9.2
    module load hdf5-mpi/1.14.1
    module load fftw-mpi/3.3.10
    module load anaconda/3/2021.11

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

The Raven and Viper clusters are automatically detected
by the make.inc script. The IPP-HGW theory cluster
requires the user to set the environment variable 
MACHINE equal to 'theoryhgw'.

There is an issue with Anaconda where it links against
it's own libraries which are not the same as the system
libraries. This doesn't affect the Fortran parts of the 
code. The C-Types interface we use in Python to access
routines in libstell.so is affected. The knock-on effect
is that the benchmarks will fail with error messages
about libssh.so. The fix is for the user to add
export LD_PRELOAD=/usr/lib64/libcrypto.so.1.1:/usr/lib64/libssl.so.1.1
to your .bashrc (or equivalent). This is an
Anaconda/MPCDF issue, not realated to STELLOPT directly.

