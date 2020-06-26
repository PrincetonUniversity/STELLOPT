# -*- mode: makefile -*-

#######################################################################
#            Define User Specific Output Paths
#######################################################################
  # Set a default directory if one has not already been defined.
  STELLOPT_HOME ?= $(STELLOPT_PATH)/bin


#######################################################################
#            Define Basic Utilities
#######################################################################
  # Temporaryily copy STELLOPT_HOME to MYHOME since MYHOME is currently
  # used in all the makefiles.
  MYHOME = $(STELLOPT_HOME)
  
  SHELL = /bin/sh
  PWD1 = `pwd`
  PRECOMP:= cpp -traditional-cpp -E -P -C -DMACOSX
  COMPILE = gfortran
  COMPILE_FREE = gfortran -ffree-form -ffree-line-length-none -ffixed-line-length-none
  LINK    = mpif90 $(FLAGS) -o
  LINK_AR = ar -ruvs
  LINK_C  = gcc -shared -Wl,-no_compact_unwind
#  LINK_C  = gcc -shared -Wl,--no-undefined -Wl,-export-dynamic 

#######################################################################
#            Define Compiler Flags
#######################################################################
#  FLAGS_R = -O2 -g -fexternal-blas -fbacktrace -fcheck=all 
  FLAGS_R = -O0 -g -Wmissing-include-dirs -fbacktrace -fcheck=all 
#  FLAGS_R = -O2 -fexternal-blas -Wmissing-include-dirs
#  FLAGS_D = -g -O0 -fcheck=all
  FLAGS_R = $(FLAGS_D)
#  LIBS    = -L/opt/local/lib -lopenblas -lscalapack  
  LIBS    = -L/opt/local/lib -lopenblas 
# Because scalapack is broken in macport for MPI4.0
#  LIBS += /Users/lazerson/Sims_Work/scalapack-2.0.2/libscalapack.a
  LIBS    = -L/opt/local/lib -framework Accelerate -lscalapack   -g

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpif90
  MPI_COMPILE_FREE = mpif90 -ffree-form \
                     -ffree-line-length-none -ffixed-line-length-none
  MPI_COMPILE_C = mpicc 
  MPI_LINK = mpif90
  MPI_LINK = mpif90 -shared  -Wl,-no_compact_unwind
  MPI_RUN = mpiexec
  MPI_RUN_OPTS = -np 2
#  MPI_RUN_OPTS_SM = -np 16
#  MPI_RUN_OPTS_MD = -np 64
#  MPI_RUN_OPTS_LG = -np 256

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB = -L$(NAG_ROOT)/lib -lnag_nag

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = -I/opt/local/include
  NETCDF_LIB = -L/opt/local/lib -lnetcdf -lnetcdff

#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = F
  NTCC_INC = -I/Users/lazerson/Sims_Work/pub/mod
  NTCC_LIB = -L/Users/lazerson/Sims_Work/pub/lib \
             -lpreact -lsigsub -laladdinsub -lpreact -lpspline \
             -lportlib -lezcdf

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  HDF5_INC = -I/opt/local/include
  HDF5_LIB = -L/opt/local/lib -lhdf5 -lhdf5hl_fortran -lhdf5_hl \
             -lhdf5_fortran

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT = F
  PGPLOT_INC = -I/opt/local/include
  PGPLOT_LIB = -L/opt/local/lib -lpgplot -lX11

#######################################################################
#             SILO Options
#######################################################################
  LSILO = F
  SILO_INC = -I/opt/local/include
  SILO_LIB = -L/opt/local/lib -lsilo

#######################################################################
#            DKES/NEO Options
#######################################################################
  LDKES = T
  LNEO  = T

#######################################################################
#            FFTW3 Options
#######################################################################
  LFFTW3 = F
  FFTW3_INC = 
  FFTW3_LIB = 

#######################################################################
#            GENE Options
#######################################################################
  LGENE = F
  GENE_INC = -I$(GENE_PATH)
  GENE_DIR = $(GENE_PATH)
  LIB_GENE = libgene.a
  GENE_LIB = $(GENE_DIR)/$(LIB_GENE) \
             -L/u/slazerso/src/GENE17_2016/external/pppl_cluster/futils/src -lfutils \
             -L$(FFTWHOME)/lib -lfftw3 \
             -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc \
             -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lX11

#######################################################################
#            COILOPT++ Options
#######################################################################
  LCOILOPT = F
  COILOPT_INC = -I$(COILOPT_PATH)
  COILOPTPP_DIR = $(COILOPT_PATH)
  LIB_COILOPTPP = libcoilopt++.a
  COILOPT_LIB = $(COILOPT_PATH)/$(LIB_COILOPTPP) \
                -L$(GSLHOME)/lib -lgsl -lgslcblas -lstdc++ -lmpi_cxx

#######################################################################
#            TERPSICHORE Options
#######################################################################
  LTERPSICHORE= F
  TERPSICHORE_INC = -I$(TERPSICHORE_PATH)
  TERPSICHORE_DIR = $(TERPSICHORE_PATH)
  LIB_TERPSICHORE = libterpsichore.a
  TERPSICHORE_LIB = $(TERPSICHORE_DIR)/$(LIB_TERPSICHORE)

#######################################################################
#            TRAVIS Options
#######################################################################
  LTRAVIS= F
  TRAVIS_DIR = $(TRAVIS_PATH)
  LIB_TRAVIS = libtravis64_sopt.a
  LIB_MCONF  = libmconf64.a
  TRAVIS_LIB = $(TRAVIS_DIR)/lib/$(LIB_TRAVIS) \
               $(TRAVIS_DIR)/magconf/lib/$(LIB_MCONF) -lstdc++

#######################################################################
#            SFINCS Options
#######################################################################

LSFINCS = F
#LSFINCS = T
SFINCS_DIR = /Users/mattland/sfincs/fortran/version3
SFINCS_INC = -I$(SFINCS_DIR)
LIB_SFINCS = libsfincs.a

PETSC_DIR=/opt/local/lib/petsc
include ${PETSC_DIR}/lib/petsc/conf/variables
# The above line defines the PETSC make variables PETSC_LIB, FC_LINKER_FLAGS, and FFLAGS, which are used below in SFINCS_LIB.

SFINCS_LIB = $(SFINCS_DIR)/$(LIB_SFINCS) \
	$(PETSC_LIB) $(FC_LINKER_FLAGS) $(FFLAGS)

#######################################################################
#            MANGO Options
#######################################################################

#LMANGO = F
LMANGO = T

MANGO_DIR = /Users/mattland/mango

# If LMANGO=T, any flags specified in the next 2 lines will be added to MANGO_F_COMPILE_FLAGS and MANGO_F_LINK_FLAGS respectively
# as defined in $(MANGO_DIR)/lib/mangoMakeVariables
MANGO_INC=
MANGO_LIB=

#######################################################################
#            REGCOIL Options
#######################################################################
  LREGCOIL= F
  REGCOIL_DIR = $(REGCOIL_PATH)
  REGCOIL_INC = -I$(REGCOIL_DIR) 
  LIB_REGCOIL = libregcoil.a
  REGCOIL_LIB = $(REGCOIL_DIR)/$(LIB_REGCOIL) -fopenmp

#######################################################################
#            LIBSTELL Shared Options
#######################################################################

LIB_SHARE = $(LIBS)
