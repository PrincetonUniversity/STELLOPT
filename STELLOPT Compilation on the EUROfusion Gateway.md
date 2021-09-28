STELLOPT Compilation on the EUROfusion Gateway
====================================================

This page details how to compile the STELLOPT family of codes
at [EUROfusion Gateway Machine](@https://wiki.eufus.eu/doku.php?id=start). In order to do so you will need an account on their system.  As a first step please note that to use the batch system codes must be on the `/pfs/home/$USER` directory.

Gateway
-----

    export MACHINE="efgateway"
    export STELLOPT_PATH="/pfs/home/$USER/src/STELLOPT"
    module load imasenv
    module load hdf5/1.8.17--intelmpi--2017--binary

General Notes
-------------

Note again the $HOME directory is `/afs/eufus.eu/user/` which is avialable under SLURM.  So you should build all codes in `pfs/home/$USER`.
