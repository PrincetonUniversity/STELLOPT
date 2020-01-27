#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = $(HOME)/bin_cori
  PRECOMP:= fpp -C -P -DLINUX
  COMPILE = ftn
  COMPILE_FREE = ftn
  LINK    = ftn -Bstatic -o
  LINK_AR = ar -ruv
  LINK_C  = $(CC) -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -fp-model strict -assume noold_unit_star -cxxlib 
  FLAGS_D = -g -dynamic -traceback -fp-model strict -assume noold_unit_star -cxxlib
  LIBS    = 
  FLAGS_R += -g -traceback -dynamic

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = ftn
  MPI_COMPILE_FREE = ftn
  MPI_COMPILE_C = cc
  MPI_LINK = ftn
  MPI_RUN = srun
  MPI_RUN_OPTS = -n 68 -c 4 --cpu_bind=cores
  MPI_RUN_OPTS_SM = -n 68 -c 4 --cpu_bind=cores
#  MPI_RUN_OPTS_MD = -n 544 -c 4 --cpu_bind=cores
  MPI_RUN_OPTS_LG = -n 4352 -c 4 --cpu_bind=cores

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
  HDF5_INC = 
  HDF5_LIB = 

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
  SILO_INC = -I$(SILO_DIR)/include
  SILO_LIB = -L$(SILO_DIR)/lib -lsilo

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
  LCOILOPT = T
  COILOPT_INC = -I$(COILOPT_PATH)
  COILOPTPP_DIR = $(COILOPT_PATH)
  LIB_COILOPTPP = libcoilopt++.a
  COILOPT_LIB = $(COILOPT_PATH)/$(LIB_COILOPTPP) 

#######################################################################
#            TERPSICHORE Options
#######################################################################
  LTERPSICHORE = T
  TERPSICHORE_INC = -I$(TERPSICHORE_PATH)
  TERPSICHORE_DIR = $(TERPSICHORE_PATH)
  LIB_TERPSICHORE = libterpsichore.a
  TERPSICHORE_LIB = $(TERPSICHORE_DIR)/$(LIB_TERPSICHORE)

#######################################################################
#            TRAVIS Options
#######################################################################
  LTRAVIS = T
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

