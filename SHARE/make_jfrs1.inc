#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = $(HOME)/bin
  PRECOMP:= cpp -P -traditional -DLINUX
  COMPILE = ftn
  COMPILE_FREE = ftn
  LINK    = ld $(FLAGS) -o
  LINK_AR = ar -ruv
  LINK_C  = $(CC) -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -O3 -dynamic -fpic -mcmodel=medium -ffast-math \
            -fexternal-blas -fall-intrinsics -fbacktrace
  FLAGS_D = -g -O0 -dynamic -fpic -mcmodel=medium -ffast-math \
            -fexternal-blas -fall-intrinsics -fbacktrace

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = ftn
  MPI_COMPILE_FREE = ftn -ffree-form -ffree-line-length-none \
                        -ffixed-line-length-none
  MPI_COMPILE_C = cc
  MPI_LINK = ftn -fbacktrace
  MPI_RUN = srun
  MPI_RUN_OPTS = -n 40 -c 1 --cpu_bind=cores
  MPI_RUN_OPTS_SM = -n 40 -c 1 --cpu_bind=cores
#  MPI_RUN_OPTS_MD = -n 544 -c 4 --cpu_bind=cores
  MPI_RUN_OPTS_LG = -n 40 -c 1 --cpu_bind=cores

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB =

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC =
  NETCDF_LIB = 

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  HDF5_INC = -I$(HDF5_DIR)/include
  HDF5_LIB = -L$(HDF5_DIR)/lib -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT = F
  PGPLOT_INC =
  PGPLOT_LIB = 

#######################################################################
#             SILO Options
#######################################################################
  LSILO = F
  SILO_INC = 
  SILO_LIB = 

#######################################################################
#            FFTW3 Options
#######################################################################
  LFFTW3 = F
  FFTW3_INC = 
  FFTW3_LIB = 

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
             -L$(FFTW_DIR)/lib -lfftw3 -lfftw3f \
             -L$(SLEPC_DIR)/lib -lslepc 

#######################################################################
#            COILOPT++ Options
#######################################################################
  LCOILOPT = F
  COILOPT_INC = -I$(COILOPT_PATH)
  COILOPTPP_DIR = $(COILOPT_PATH)
  LIB_COILOPTPP = libcoilopt++.a
  COILOPT_LIB = $(COILOPT_PATH)/$(LIB_COILOPTPP) 

#######################################################################
#            TERPSICHORE Options
#######################################################################
  LTERPSICHORE = F
  TERPSICHORE_INC = -I$(TERPSICHORE_PATH)
  TERPSICHORE_DIR = $(TERPSICHORE_PATH)
  LIB_TERPSICHORE = libterpsichore.a
  TERPSICHORE_LIB = $(TERPSICHORE_DIR)/$(LIB_TERPSICHORE)

#######################################################################
#            TRAVIS Options
#######################################################################
  LTRAVIS = F
  TRAVIS_DIR = $(TRAVIS_PATH)
  LIB_TRAVIS = libtravis64_sopt.a
  LIB_MCONF  = libmconf64.a
  TRAVIS_LIB = $(TRAVIS_DIR)/lib/$(LIB_TRAVIS) \
               $(TRAVIS_DIR)/mag_conf/lib/$(LIB_MCONF) -lstdc++
#  TRAVIS_LIB = $(TRAVIS_DIR)/lib/$(LIB_TRAVIS) \
#               $(TRAVIS_DIR)/magconf/lib/$(LIB_MCONF) -lstdc++

#######################################################################
#            SFINCS Options
#######################################################################

  LSFINCS = F
  SFINCS_DIR = $(SFINCS_PATH)
  SFINCS_INC = -I$(SFINCS_DIR)
  LIB_SFINCS = libsfincs.a
  SFINCS_LIB = $(SFINCS_DIR)/$(LIB_SFINCS)

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

