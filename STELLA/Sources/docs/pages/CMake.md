---
title: CMake in stella
subtitle: Some notes on using and developing CMake for stella
---

Stella now has (experimental) support for building with CMake.

# Using CMake

## CMake options

All CMake options can be set on the command line with the syntax:
`-D<variable>=value`. Boolean flags can use `on/yes/true/1` or
`off/no/false/0` to turn them on/off -- these options are
case-insensitive.

### Optimisation vs debugging

CMake has a built-in option for setting optimisation or debugging
flags: `CMAKE_BUILD_TYPE`. If this is not set, by default stella will
use `RelWithDebInfo`, which is the equivalent of `-O2 -g`: moderate
optimisation with debug symbols. This doesn't turn on any runtime
checking, it only keeps names of functions, variables, etc. for error
messages and backtraces in the event of a crash.

Set `-DCMAKE_BUILD_TYPE=Release` to turn on full optimisations.

Use `-DCMAKE_BUILD_TYPE=Debug` to turn _off_ optimisations and turn on
various compile- and run-time checks (depending on the compiler).

## Dependencies

Stella has several optional dependencies, the location of which can be
specified by using the `<package>_ROOT` variables. The list of
dependencies and their location variables are as follows:

- MPI: This is the odd-one-out, in that the best way to control which
  MPI implementation is found is via `MPIEXEC_EXECUTABLE`. This is
  already automatically set to the current `mpirun/mpiexec` in your
  `PATH`, so this shouldn't need to be set.
- LAPACK: `LAPACK_ROOT`. On Cray systems that use the Cray Programming
  Environment, this is automatically handled by the Cray Compiler
  Environment and so is not user-controllable.
- FFTW: `FFTW_ROOT`. Stella also searches for the `fftw-wisdom`
  executable in your `PATH` as a first guess
- NetCDF: `netCDF_ROOT` for the C library, and `netCDFFortran_ROOT`
  for the Fortran API. Stella searches for `nc-config` and `nf-config`
  in your `PATH` and uses those to query the netCDF configuration
  
All of these dependencies, with the exception of MPI, can be turned on
or off with the `STELLA_ENABLE_<name>` variable.

# Developing the stella CMake build system

One important consideration when developing stella is that any new
files _must_ be listed in the `STELLA_SOURCES_*` variables: either
`STELLA_SOURCES_f90` _or_ `STELLA_SOURCES_fpp` as appropriate. If you
add a new file and do not add it to exactly one of these variables,
you will get a build error.
