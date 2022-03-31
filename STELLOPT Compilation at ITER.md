STELLOPT Compilation on the ITER Cluster
====================================================

This page details how to compile the STELLOPT family of codes on the ITER cluster.
Since this cluster is only for ITER users we assume you know how to use and access it.

ITER Cluster
-----

    export MACHINE="iter"
    export STELLOPT_PATH="$HOME/src/STELLOPT"
    module load IMAS/3.35.0-4.10.0-2020b
    module load netCDF-Fortran/4.5.3-iimpi-2020b
    module load iWrap/0.4.1-intel-2020b
    module load FFTW/3.3.8-iimpi-2020b
    module load XMLlib/3.3.1-intel-2020b

Building IMAS Actors (iWrap)
-------------

To build the IMAS Actors first make sure everything compiles using the ./build_all script.
Then in the IMAS_INTERFACE directory issue the command

    make actor_vmec

This will build the VMEC actor
