STELLOPT Compilation with Anaconda
==============================

This page details how to compile the STELLOPT family of codes using
[Anaconda](@https://anaconda.org/). In order to do this you'll need
to install an Anaconda variant.

-----

Installation Steps

1\. Create your environment (here named `stellopt_env`) and install its dependencies

    conda create -n stellopt_env gfortran gxx hdf5[build="mpi*"] hdf5tools netcdf-fortran netcdf4 openblas openmpi-mpifort python scalapack

2\. Activate your environment

    conda activate stellopt_env

3\. Setup your environment variables

    export MACHINE="conda"
    export STELLOPT_PATH=<path to repo directory>

Note: there is some question as to whether scalapack is mpi-enabled.
