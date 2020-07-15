#######################################################################
#            Define User Specific Output Paths
#######################################################################
  # Set a default directory if one has not already been defined.
  STELLOPT_HOME ?= $(HOME)/bin


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
  LINK_C  = gcc -shared -Wl,--no-undefined -Wl,-export-dynamic 

#######################################################################
#            Libraries path
#######################################################################
  SCALAPACK_HOME ?= /usr/local/opt/scalapack
  FFTW_HOME ?= /usr/local/opt/fftw
  HDF5_HOME ?= /usr/local/opt/hdf5
  NETCDF_HOME ?= /usr/local/opt/netcdf
  # netcdf needs manually installed netcdf-frotran; still have some issue
  

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -O2 -g -fbacktrace -fcheck=all 
  FLAGS_D = -O0 -g -fcheck=all -fbacktrace -Wextra \
            -Wtarget-lifetime -fbounds-check -ffpe-trap=zero -finit-real=snan
  LIBS    = -L${SCALAPACK_HOME}/lib -lscalapack

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpifort
  MPI_COMPILE_FREE = mpif90 -ffree-form \
                     -ffree-line-length-none -ffixed-line-length-none
  MPI_COMPILE_C = mpicc 
  MPI_LINK = mpifort
  MPI_LINK = mpifort -shared  -Wl,-no_compact_unwind
  MPI_RUN = mpirun
  MPI_RUN_OPTS = -np 2

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB = -L$(NAG_ROOT)/lib -lnag_nag

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = -I${NETCDF_HOME}/include
  NETCDF_LIB = -L${NETCDF_HOME}/lib -lnetcdf -lnetcdff

#######################################################################
#            FFTW3 Options
#######################################################################
  LFFTW3 = T
  FFTW3_INC = -I${FFTW_HOME}/include
  FFTW3_LIB = -L${FFTW_HOME}/lib -lfftw3

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  HDF5_INC = -I${HDF5_HOME}/include
  HDF5_LIB = -L${HDF5_HOME}/lib -lhdf5hl_fortran -lhdf5_hl \
             -lhdf5_fortran -lhdf5 -lz -ldl -lm

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT = T
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
#            REGCOIL Options
#######################################################################
  LREGCOIL= F
  REGCOIL_DIR = $(REGCOIL_PATH)
  REGCOIL_INC = -I$(REGCOIL_DIR) 
  LIB_REGCOIL = libregcoil.a
  REGCOIL_LIB = $(REGCOIL_DIR)/$(LIB_REGCOIL) -fopenmp

#######################################################################
#            SFINCS Options
#######################################################################

  LSFINCS = F
  SFINCS_DIR = $(SFINCS_PATH)
  SFINCS_INC = -I$(SFINCS_DIR)
  LIB_SFINCS = libsfincs.a
  SFINCS_LIB = $(SFINCS_DIR)/$(LIB_SFINCS) \
             -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lX11

#######################################################################
#            LIBSTELL Shared Options
#######################################################################

LIB_SHARE = -lc -lgfortran -lstdc++ -lmpi -lmpi_mpifh -lz -lc -lm -lpthread $(LIBS) -lc
