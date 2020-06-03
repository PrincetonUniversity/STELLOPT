STELLOPT Compilation at the PPPL cluster
========================================

![](images/PPPL-LOGO-FNLWH-GRADIENT_300px_WEB.jpg)

This page details how to compile the STELLOPT family of codes at the
[PPPL](@http://www.pppl.gov/) cluster. The proper modules need to be
loaded then compilation can begin. Please note that if you require
additional module to be loaded you should do this before loading the
compiler. This will prevent the \$PATH variable from searching the wrong
directories.

------------------------------------------------------------------------

### GENERAL Instructions

On the PPPL cluster, STELLOPT is now maintained as an installed module.
Use the following command to load stellopt

    module load stellopt

To load the most current version of the code. This will load all the
necessary modules. The /bin/ directory will be added to your path
variable as well so codes can be called without specifying the full
path. Such as:

	mpirun -np 64 xstelloptv2 input.test

There are more avaliable versions. To list all the available versions,
you can do the following

```
module avail stellopt
```

To check the details of each version, you can use `module what-is
stellopt` or `module show stellopt`.


------------------------------------------------------------------------

### GNU

Load the appropriate module files for the [GNU](https://gcc.gnu.org/)
compiler.

```
     module load gcc/8.1.0
     module load szip
     module load openmpi
     module load gsl
     module load hdf
     module load scalapack
     module load blacs
     module load fftw
     module load hdf5-parallel
     module load curl
     module load netcdf-c
     module load netcdf-fortran
     module load netcdf-cxx4
     module load matlab
     module load python/3.6.4
```

------------------------------------------------------------------------

### Intel

Load the appropriate module files for the [Intel](https://software.intel.com/en-us/fortran-compilers)
compiler.

```
   module load intel
   module load openmpi
   module load szip
   module load hdf
   module load hdf5-parallel
   module load netcdf-c
   module load netcdf-fortran
   module load lapack
   module load petsc_complex
   module load fftw
   module load slepc_complex
   module load scalapack
   module load blacs
   module load python/3.6.4
```

The PPPL cluster should be automatically detected otherwise
please set `MACHINE=pppl_gcc` to properly compile the code.
