#######################################################################
#     This is a custom file compiling libstell on windows under
#     msys/mingw.  Follow the instructions under README_mscompile.md
#     to setup the environment
#
#
# Note:  -Wl,-z-defs causes the compiler to fail immediately if there are undefined symbols in a library
#        -Wl,--allow-multiple-definition will miss double definition errors in libraries
#        -Wl,--enable-auto-import - allow omittion of decorations in code for dllimport/dllexport
#        -Wl,--enable-runtime-pseudo-relocs -  allow relinking for structs/classes/etc in object files
# -Wl,--export-all-symbols
# -Wl,--whole-archive
# -static-libstdc++
# -static-libgcc
# -static-libgfortran
#
# NOTE:
#  "-Wl,--export-dynamic" is not supported for PX+ type targets (windows+mingw)
#
#
#  gcc9 has a few differences in system calls. We put in some kluges for this
#  by making a WIN64 macro.
#
#######################################################################
#          If you want to explicitly override user given flags
#######################################################################
  FC = gfortran
  CC = gcc

HARD_FFLAGS = \
              -Wl,--export-all-symbols

HARD_CFLAGS = \

HARD_LINK_FLAGS =

#######################################################################
#            Find the source paths of the STELLOPT distribution.
#######################################################################
  CURRENT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
  STELLOPT_DIR:=$(abspath $(CURRENT_DIR)/../)

#######################################################################
#            Define User Specific Output Paths.
#######################################################################
  # Define the output directory.  This should be an absolute path.
  STELLOPT_HOME = $(STELLOPT_DIR)/bin

#######################################################################
#            Define Basic Utilities
#######################################################################
  # Temporaryily copy STELLOPT_HOME to MYHOME since MYHOME is currently
  # used in all the makefiles.
  MYHOME = $(STELLOPT_HOME)

  SHELL = /bin/sh
  PWD1 = `pwd`
  PRECOMP:= cpp -traditional-cpp -E -P -C -DMACOSX
  COMPILE = $(FC) $(HARD_FFLAGS)
  COMPILE_FREE = $(FC) -ffree-form -ffree-line-length-none \
                     -ffixed-line-length-none $(HARD_FFLAGS)
  LINK    = $(FC) $(FLAGS) $(HARD_LINK_FLAGS) -o
  LINK_AR = ar -ruvs
  LINK_C  = $(CC) -shared \
            -Wl,-no_compact_unwind \
            $(HARD_CFLAGS)

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -O2 -fexternal-blas -Wmissing-include-dirs
  FLAGS_D = -g -O0 -fcheck=all

#######################################################################
#            MPI Options
#
# MPI RUN OPTIONS:
#  MPI_RUN_OPTS = -np 4
#  MPI_RUN_OPTS_SM = -np 16
#  MPI_RUN_OPTS_MD = -np 64
#  MPI_RUN_OPTS_LG = -np 256
#
# Let's not pollute the namespace if unnecessary
#######################################################################
  LMPI    = T

ifeq ($(LMPI),T)
  MPI_COMPILE = mpifort $(HARD_FFLAGS)
  MPI_COMPILE_FREE = mpifort -ffree-form -fno-range-check \
                     -ffree-line-length-none -ffixed-line-length-none \
                     $(HARD_FFLAGS)
  MPI_COMPILE_C = mpicc
  MPI_LINK = mpicc -shared -pthread \
             -Wl,-no_compact_unwind \
             $(HARD_CFLAGS)

  MPI_RUN = mpiexec
  MPI_RUN_OPTS = -np 4

  MPI_LIB = -L$(PORTSDIR)/lib/mpich-gcc9 -lmpi -lmpifort
endif

#######################################################################
#            Define LIBS:
#            dependency order matters (highest left, lowest right)
#
#   SCALAPACK depends on FFTW3 and can depend on MPI
#     Make sure that you set those flags appropriately depending on your build
#
#     NOTE: libm is necessary for basic math utilities
#
#######################################################################

  LSCALAPACK = T
  LOPENBLAS = T

ifeq ($(LSCALAPACK),T)
  ifeq ($(LOPENBLAS),T)
    LIBS    = -L$(PORTSDIR)/lib -lscalapack -lopenblas
  else
    LIBS    = -L$(PORTSDIR)/lib -lscalapack -llapack -lblas
  endif
else
  LIBS    = -L$(PORTSDIR)/lib -llapack -lblas
endif

#######################################################################
#            FFTW3 Options
#######################################################################
  LFFTW3 = T

ifeq ($(LFFTW3),T)
  FFTWHOME = $(PORTSDIR)
  FFTW3_INC = -I$(FFTWHOME)/include
  FFTW3_LIB = -L$(FFTWHOME)/lib -lfftw3_threads -lfftw3

  ifeq ($(LMPI),T)
    FFTW3_LIB := -L$(FFTWHOME)/lib -lfftw3_mpi $(FFTW3_LIB)
  endif

endif

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F

ifeq ($(LNAG),T)
  NAG_ROOT =
  NAG_LIB =
#  NAG_LIB = -L$(NAG_ROOT)/lib -lnag_nag
endif

#######################################################################
#            TRAVIS Options
#
#    TRAVIS depends on the MCONF library and the standard C++ library
#    MCONF depends on the stdc++ library and libc
#######################################################################
  LTRAVIS= F

ifeq ($(LTRAVIS),T)
  TRAVIS_DIR = $(TRAVIS_PATH)
  LIB_TRAVIS = libtravis64_sopt.a
  LIB_MCONF  = libmconf64.a
  TRAVIS_LIB = -L$(TRAVIS_DIR)/lib/$(LIB_TRAVIS) \
               -L$(TRAVIS_DIR)/mag_conf/lib/$(LIB_MCONF) \
               -L$(PORTSDIR)/lib -lstdc++ \
               -L$(PORTSDIR)/lib/gcc9 -lgfortran \
               -L/usr/lib -pthread -lm -lc
endif

#######################################################################
#            SFINCS Options
#
#   Note:  PETSC has a bunch of dependencies
#
#   Also: The included variables file defines the PETSC make variables
#    PETSC_LIB, FC_LINKER_FLAGS, and FFLAGS, which are used below in
#    SFINCS_LIB. Headsup for namespace conflicts: overlapping FFLAGS?
#
#   PETSC_DIR should be set as user env variable
#
#######################################################################

  LSFINCS = F

ifeq ($(LSFINCS),T)
  SFINCS_DIR = $(SFINCS_PATH)
  SFINCS_INC = -I$(SFINCS_DIR)
  LIB_SFINCS = libsfincs.a

  include ${PETSC_DIR}/lib/petsc/conf/variables

  SFINCS_LIB = -L$(SFINCS_DIR)/$(LIB_SFINCS) \
              -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lX11
endif

#######################################################################
#            NETCDF Options
#
# NOTE: NETCDF depends on HDF5 (and possibly your pthread library)
#    netcdf-fortran depends on netcdf
#    netcdf-c depends on HDF5 and its libraries
#     (but netcdf-c "knows" where to find HDF5 libraries and libz, libcurl)
#
# Optionally you can use to automatically set the dependencies / order
#    this slows down compilation on my computer
#  NETCDF_INC = $(shell nc-config --fflags)  # fortran flags
#  NETCDF_LIB = $(shell nc-config --flibs)   # fortran libs
#######################################################################
  LNETCDF = T

ifeq ($(LNETCDF),T)
  NETCDF_C_HOME = $(PORTSDIR)
  NETCDF_F_HOME = $(PORTSDIR)
  NETCDF_C = $(NETCDF_C_HOME)
  NETCDF_F = $(NETCDF_F_HOME)
  NETCDF_INC = $(shell nc-config --fflags)
  NETCDF_LIB = $(shell nc-config --flibs)
#               -L$(PORTSDIR)/gcc9 -lgfortran \
#               -L/usr/lib -pthread -lm -lc
endif

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T

ifeq ($(LHDF5),T)
  HDF5_HOME = $(PORTSDIR)
  HDF5_INC = -I$(HDF5_HOME)/include

  HDF5_LIB = -L$(HDF5_HOME)/lib -lhdf5_hl_fortran -lhdf5_fortran -lhdf5_hl \
             -lhdf5 -lz -lcurl \
             -L/usr/lib -pthread -lm -lc
endif

#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = F

ifeq ($(LNTCC),T)
  NTCC_INC =
  NTCC_LIB =
endif

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT = F

ifeq ($(LPGPLOT),T)
  PGPLOT_DIR = $(PORTSDIR)
  PGPLOT_INC = -I$(PORTSDIR)/include
  PGPLOT_LIB = -L$(PGPLOT_DIR) -lpgplot -L/usr/X11R6/lib -lX11
endif

#######################################################################
#             SILO Options
#######################################################################
  LSILO = F

ifeq ($(LSILO),T)
  SILOHOME = $(PORTSDIR)
  SILO_INC = -I$(SILOHOME)/include
  SILO_LIB = -L$(SILOHOME)/lib -lsilo
endif

#######################################################################
#            DKES/NEO Options
#######################################################################
  LDKES = T
  LNEO  = T

#######################################################################
#            GENE Options
#######################################################################
  LGENE = F

ifeq ($(LGENE),T)
  GENE_INC = -I$(GENE_PATH)
  GENE_DIR = $(GENE_PATH)
  LIB_GENE = libgene.a
  GENE_LIB = $(GENE_DIR)/$(LIB_GENE) \
             -L/u/slazerso/src/GENE17_2016/external/pppl_cluster/futils/src -lfutils \
             -L$(FFTWHOME)/lib -lfftw3 \
             -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc \
             -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lX11
endif

#######################################################################
#            COILOPT++ Options
#######################################################################
  LCOILOPT = F

ifeq ($(LCOILOPT),T)
  COILOPT_INC = -I$(COILOPT_PATH)
  COILOPTPP_DIR = $(COILOPT_PATH)
  LIB_COILOPTPP = libcoilopt++.a

  GSLHOME = $(PORTSDIR)
  COILOPT_LIB = $(COILOPT_PATH)/$(LIB_COILOPTPP) \
                -L$(GSLHOME)/lib -lgsl -lgslcblas \
                -L$(PORTSDIR)/lib -lstdc++ \
                -L$(PORTSDIR)/lib/mpich-gcc9 -lmpi_cxx
endif

#######################################################################
#            TERPSICHORE Options
#######################################################################
  LTERPSICHORE= F

ifeq ($(LTERPSICHORE),T)
  TERPSICHORE_INC = -I$(TERPSICHORE_PATH)
  TERPSICHORE_DIR = $(TERPSICHORE_PATH)
  LIB_TERPSICHORE = libterpsichore.a
  TERPSICHORE_LIB = $(TERPSICHORE_DIR)/$(LIB_TERPSICHORE)
endif

#######################################################################
#            REGCOIL Options
#######################################################################
  LREGCOIL= F

ifeq ($(LREGCOIL),T)
  REGCOIL_DIR = $(REGCOIL_PATH)
  REGCOIL_INC = -I$(REGCOIL_DIR)
  LIB_REGCOIL = libregcoil.a
  REGCOIL_LIB = $(REGCOIL_DIR)/$(LIB_REGCOIL) -fopenmp
endif

#######################################################################
#            LIBSTELL Shared Options
#######################################################################
# Using c-compiler for link step:  need to include gfortran library
#    Recursion through libraries is necessary to resolve undefined symbols
#    Order matters! Moving left to right:
#            search LIBS  and   add items to undefined symbol / function table
#            re-search LIBS because there are recursive dependencies
#            check libc resolve some and add more undefined symbols
#            check gfortran lib and resolve more undefined symbols
#            ... add some math/compression undefined symbols
#            check zlib and resolve ... add more from gfortran / libc
#            resolve the new ones from libc
#            check libm for math symbols
#            resolve the new ones from libc
#            finally resolve the c-compiler undefined symbols for gfortran
#
#   The recursive tool in ld is not available on Darwin. We repeat
#   the shared library in make_all.inc after appending code libraries
#
LIB_SHARE = $(LIBS) \
            -L$(PORTSDIR)/lib -lz
ifeq ($(LMPI),T)
  LIB_SHARE += -L$(PORTSDIR)/lib/mpich-gcc9 -lmpi -lmpifort
endif
LIB_SHARE += \
            -L$(PORTSDIR)/lib/gcc9 -lgfortran -lstdc++ \
            -L/usr/lib/ -lm -lpthread -lc

######################################################################
######################################################################
