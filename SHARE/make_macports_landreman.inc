# -*- mode: makefile -*-

#######################################################################
#            Define User Specific Output Paths
#######################################################################
  # Set a default directory if one has not already been defined.
#  $STELLOPT_HOME ?= $(HOME)/stellopt/bin
  STELLOPT_HOME = $(HOME)/stellopt/bin


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
  LINK    = ld $(FLAGS) -o
  LINK_AR = ar -ruvs
  LINK_C  = gcc -shared -Wl,-no_compact_unwind

#######################################################################
#            Define Compiler Flags
#######################################################################
#  FLAGS_R = -O2 -g -fexternal-blas -fbacktrace -fcheck=all 
#  FLAGS_R = -O2 -fexternal-blas -Wmissing-include-dirs
#  FLAGS_R = -O2 -Wmissing-include-dirs
  FLAGS_R = -g -O0  -Wmissing-include-dirs
  FLAGS_D = -g -O0 -fcheck=all 
#  LIBS    = -L/opt/local/lib -lopenblas -lscalapack  
  LIBS    = -L/opt/local/lib -framework Accelerate -lscalapack   -g

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = /opt/local/bin/mpif90 
  MPI_COMPILE_FREE = /opt/local/bin/mpif90 -ffree-form \
                     -ffree-line-length-none -ffixed-line-length-none
  MPI_COMPILE_C = gcc -shared -Wl,-no_compact_unwind
  MPI_LINK = /opt/local/bin/mpif90
  MPI_RUN = mpiexec
  MPI_RUN_OPTS = -np 4
#  MPI_RUN_OPTS_SM = -np 16
#  MPI_RUN_OPTS_MD = -np 64
#  MPI_RUN_OPTS_LG = -np 256

  # For GCC6
  #MPI_COMPILE = mpif90
  #MPI_COMPILE_FREE = mpif90 -ffree-form \
  #                   -ffree-line-length-none -ffixed-line-length-none
  #MPI_COMPILE_C = mpicc
  #MPI_LINK = mpif90

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

  # FOR GCC6
  #NETCDF_LIB = -L$(NETCDF_HOME)/lib -lnetcdff

#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = F
  NTCC_INC = -I/Users/slazerso/Sims_PPPL/STELLOPT/NTCC/MAC/mod
  NTCC_LIB = -L/Users/slazerso/Sims_PPPL/STELLOPT/NTCC/MAC/lib \
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
  PGPLOT_INC = -I$(PGPLOT_DIR)
  PGPLOT_LIB = -L$(PGPLOT_DIR) -lpgplot -L/usr/X11R6/lib -lX11

#######################################################################
#             SILO Options
#######################################################################
  LSILO = F
  SILO_INC = -I$(SILOHOME)/include
  SILO_LIB = -L$(SILOHOME)/lib -lsilo

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
               $(TRAVIS_DIR)/mag_conf/lib/$(LIB_MCONF) -lstdc++

#######################################################################
#            SFINCS Options
#######################################################################

#LSFINCS = F
LSFINCS = T
SFINCS_DIR = /Users/mattland/sfincs/fortran/version3
SFINCS_INC = -I$(SFINCS_DIR)
LIB_SFINCS = libsfincs.a

PETSC_DIR=/opt/local/lib/petsc
include ${PETSC_DIR}/lib/petsc/conf/variables
# The above line defines the PETSC make variables PETSC_LIB, FC_LINKER_FLAGS, and FFLAGS, which are used below in SFINCS_LIB.

SFINCS_LIB = $(SFINCS_DIR)/$(LIB_SFINCS) \
	$(PETSC_LIB) $(FC_LINKER_FLAGS) $(FFLAGS)

#######################################################################                                                                                       
#            REGCOIL Options                                                                                                                                  
#######################################################################                                                                                       
  LREGCOIL= F
#  REGCOIL_DIR = $(REGCOIL_PATH)
  REGCOIL_DIR = /Users/mattland/regcoil
  REGCOIL_INC = -I$(REGCOIL_DIR)
  LIB_REGCOIL = libregcoil.a
  REGCOIL_LIB = $(REGCOIL_DIR)/$(LIB_REGCOIL) -fopenmp

#######################################################################
#            LIBSTELL Shared Options
#######################################################################
LIB_SHARE = /opt/local/lib/libnetcdf.dylib \
         /opt/local/lib/libnetcdff.dylib \
         /opt/local/lib/libnetcdf_c++4.dylib \
         /opt/local/lib/libhdf5.dylib \
         /opt/local/lib/libhdf5hl_fortran.dylib\
         /opt/local/lib/libhdf5_hl.dylib \
         /opt/local/lib/libhdf5_fortran.dylib\
         /opt/local/lib/gcc6/libgfortran.dylib \
         /opt/local/lib/gcc6/libstdc++.dylib \
         /opt/local/lib/libopenblas.dylib \
         /opt/local/lib/openmpi-gcc6/libmpi.dylib \
         /opt/local/lib/openmpi-gcc6/libmpi_mpifh.12.dylib \
         /opt/local/lib/libz.dylib \
         /usr/lib/libc.dylib \
         /usr/lib/libm.dylib \
         /usr/lib/libpthread.dylib


