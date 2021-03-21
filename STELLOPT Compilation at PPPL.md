STELLOPT Compilation & Usage at the PPPL clusters
========================================

![](images/PPPL-LOGO-FNLWH-GRADIENT_300px_WEB.jpg)

This page details how to compile/use the STELLOPT family of codes at the
[PPPL](@http://www.pppl.gov/) clusters. The proper modules need to be
loaded then compilation can begin.

PPPL has several clusters, "portal", "stellar", and "traverse".
For more information, please with [PPPL Research Computing](https://pppl-intranet.princeton.edu/departments/computing-and-information-technology/research-computing).

------------------------------------------------------------------------

# Compilation

On the PPPL clusters, there are third-party modules maintained by Dr. Caoxiang Zhu.
Here are the information to load the existing modules.

## Portal

Intel compiler with OpenMPI
```batch
module load mod_stellopt
module load stellopt/intel 
```
GCC compiler with OpenMPI
```batch
module load mod_stellopt
module load stellopt/gcc 
```

You can set the env variable `$STELLOPT_PATH` to the folder containing STELLOPT sources.
By default, it will use the makefiles `SHARE/make_pppl_intel.inc` and `SHARE/make_pppl_gcc.inc`, respectively.
You can copy and customize the makefile based on your needs.

**PPPL portal has issues with MPI-shared memory. It causes errors in `BEAMS3D` (see [#109](https://github.com/PrincetonUniversity/STELLOPT/issues/109)) and in free-boundary VMEC (occasionally).**

## Stellar

Intel compiler with Intel MPI
```batch
module use /home/caoxiang/module
module load stellopt/intel 
```
GCC compiler (to be updated)

Stellar is using Intel CPUs, so it is recommended to use the Intel compiler.
The corresponding makefile is `SHARE/make_stellar.inc`.

## Traverse

Traverse is a GPU-based cluster, so it is not recommended to use STELLOPT on it.
If you would like to use STELLLOPT on it, please refer to the general compilation page.

# Usage

If you want to use the compiled executables directly, you can load the following modules directly.
The `/bin/` directory will be added to your path variable as well so codes can be called without specifying the full path. Such as:

	mpirun -np 64 xstelloptv2 input.test

There are more avaliable versions.
You can use `module avail stellopt` to list all the possible versions. 

## Portal

The `develop` branch compiled by Intel compiler with OpenMPI
```batch
module load mod_stellopt
module load stellopt/develop_intel
```
The `develop` branch compiled by GCC compiler with OpenMPI (default)
```batch
module load mod_stellopt
module load stellopt/develop_gcc 
```

## Stellar

The `develop` branch compiled by Intel compiler with Intel MPI (default)
```batch
module use /home/caoxiang/module
module load stellopt/develop_intel
```
GCC compiler (to be updated)

