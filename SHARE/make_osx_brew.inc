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
  # temp add some checking flags to avoid BEAMS3d issues
  FLAGS_R = -O2 -g -fexternal-blas -march=native
  FLAGS_D = -O0 -g -fexternal-blas -fbacktrace -fcheck=all,no-array-temps -fbounds-check
  LIBS    = -L${SCALAPACK_HOME}/lib -lscalapack -framework Accelerate
  FLAGS_R += -fallow-argument-mismatch
  FLAGS_D += -fallow-argument-mismatch

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpif90
  MPI_COMPILE_FREE = mpif90 -ffree-form \
                     -ffree-line-length-none -ffixed-line-length-none 
  MPI_COMPILE_C = mpicc 
  MPI_LINK = mpif90
  MPI_LINK = mpif90 -Wl,-no_compact_unwind
  MPI_RUN = mpiexec
  MPI_RUN_OPTS = --use-hwthread-cpus

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB = -L$(NAG_ROOT)/lib -lnag_nag

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = $(shell nf-config --fflags)
  NETCDF_LIB = -L${NETCDF_HOME}/lib -lnetcdf -lnetcdff
  NETCDF_LIB = $(shell nf-config --flibs)

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
  HDF5_LIB = -L${HDF5_HOME}/lib -lhdf5_hl_fortran -lhdf5_hl \
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
ifneq ("$(wildcard $(GENE_PATH))","")
  LGENE = T
  GENE_INC = -I$(GENE_PATH)
  GENE_DIR = $(GENE_PATH)
  LIB_GENE = libgene.a
  GENE_LIB = 
else
  LGENE = F
endif

#######################################################################
#            COILOPT++ Options
#######################################################################
ifneq ("$(wildcard $(COILOPT_PATH))","")
  LCOILOPT = T
  COILOPT_INC = -I$(COILOPT_PATH)
  COILOPTPP_DIR = $(COILOPT_PATH)
  LIB_COILOPTPP = libcoilopt++.a
  COILOPT_LIB = $(COILOPT_PATH)/$(LIB_COILOPTPP) \
                -L$(GSLHOME)/lib/x86_64-linux-gnu -lgsl -lgslcblas -lstdc++ -lmpi_cxx
else
  LCOILOPT = F
endif

#######################################################################
#            TERPSICHORE Options
#######################################################################
ifneq ("$(wildcard $(TERPSICHORE_PATH))","")
  LTERPSICHORE= T
  TERPSICHORE_INC = -I$(TERPSICHORE_PATH)
  TERPSICHORE_DIR = $(TERPSICHORE_PATH)
  LIB_TERPSICHORE = libterpsichore.a
  TERPSICHORE_LIB = $(TERPSICHORE_DIR)/$(LIB_TERPSICHORE)
else
  LTERPSICHORE = F
endif

#######################################################################
#            TRAVIS Options
#######################################################################
ifneq ("$(wildcard $(TRAVIS_PATH))","")
  LTRAVIS = T
  TRAVIS_DIR = $(TRAVIS_PATH)
  LIB_TRAVIS = libtravis64_sopt.a
  LIB_MCONF  = libmconf64.a
  LIB_FADDEEVA = libfaddeeva64.a
  TRAVIS_LIB = $(TRAVIS_DIR)lib/$(LIB_TRAVIS) \
               $(TRAVIS_DIR)faddeeva_package/lib/$(LIB_FADDEEVA) \
               $(TRAVIS_DIR)magconf/lib/$(LIB_MCONF) -lc++
else
  LTRAVIS = F
endif

#######################################################################
#            REGCOIL Options
#######################################################################
ifneq ("$(wildcard $(REGCOIL_PATH))","")
  LREGCOIL= T
  REGCOIL_DIR = $(REGCOIL_PATH)
  REGCOIL_INC = -I$(REGCOIL_DIR) 
  LIB_REGCOIL = libregcoil.a
  REGCOIL_LIB = $(REGCOIL_DIR)/$(LIB_REGCOIL) -fopenmp
else
  LREGCOIL = F
endif

#######################################################################
#            SFINCS Options
#######################################################################
ifneq ("$(wildcard $(SFINCS_PATH))","")
  LSFINCS = T
  SFINCS_DIR = $(SFINCS_PATH)
  SFINCS_INC = -I$(SFINCS_DIR)
  LIB_SFINCS = libsfincs.a
  SFINCS_LIB = $(SFINCS_DIR)/$(LIB_SFINCS) \
             -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lX11
else
  LSFINCS = F
endif

#######################################################################
#            Available Energy Options
#######################################################################
ifneq ("$(wildcard $(AEOPT_PATH))","")
  LAEOPT = T
  AEOPT_DIR = $(AEOPT_PATH)
  AEOPT_INC = -I$(AEOPT_DIR) 
  LIB_AEOPT = libtrapAE.a
  AEOPT_LIB = $(AEOPT_PATH)/$(LIB_AEOPT)
else
  LAEOPT = F
endif

#######################################################################
#            MANGO Options
#######################################################################
ifneq ("$(wildcard $(MANGO_PATH))","")
  LMANGO = T
  MANGO_DIR = $(MANGO_PATH)
# Any flags specified on the next 2 lines will be added to the
# MANGO_F_COMPILE_FLAGS and MANGO_F_LINK_FLAGS lines as defined in
# $(MANGO_DIR)/lib/mangoMakeVariables
  MANGO_INC = 
  MANGO_LIB = 
else
  LMANGO = F
endif

#######################################################################
#            LIBSTELL Shared Options
#######################################################################

LIB_SHARE = -lc -lgfortran -lstdc++ -lmpi -lmpi_mpifh -lz -lc -lm -lpthread $(LIBS) -lc
