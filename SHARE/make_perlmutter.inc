#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = /global/cfs/cdirs/m4505/sw/STELLOPT/bin
  PRECOMP:= /lib/cpp -traditional -DLINUX
  COMPILE = $(FC)
  COMPILE_FREE = $(FC)
  LINK    = $(FC) -Bstatic -o
  LINK_AR = ar -ruv
  LINK_C  = $(CC) -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -g -fbacktrace -dynamic -fallow-argument-mismatch -ffree-line-length-0
  FLAGS_D = -g -O0 -fexternal-blas -fcheck=all -fbacktrace -Wextra \
   -Wtarget-lifetime -fbounds-check -ffpe-trap=zero -finit-real=snan \
   -ffree-line-length-0 -dynamic -fallow-argument-mismatch 
  LIBS    = 

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = $(MPIFC)
  MPI_COMPILE_FREE = $(MPIFC)
  MPI_COMPILE_C = $(MPICC)
  MPI_LINK = $(FC)
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
  #NETCDF_DIR = ## Will be set by module 
  NETCDF_INC = -I$(NETCDF_DIR)/include
  NETCDF_LIB = -L$(NETCDF_DIR)/lib -lnetcdff

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  PRECOMP+= -DHDF5_PAR
  #HDF5_DIR = ## Will be set by module 
  HDF5_INC = -I$(HDF5_DIR)/include 
  HDF5_LIB = -L$(HDF5_DIR)/lib -lhdf5_hl -lhdf5_fortran -lhdf5 \
             -lhdf5hl_fortran -lpthread -lz -lm 


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
#            Available Energy Options
#######################################################################
  LAEOPT= F
  AEOPT_DIR = $(AEOPT_PATH)
  AEOPT_INC = -I$(AEOPT_DIR) 
  LIB_AEOPT = libtrapAE.a
  AEOPT_LIB = $(AEOPT_PATH)/$(LIB_AEOPT)

#######################################################################
#            LIBSTELL Shared Options
#######################################################################
LIB_SHARE = $(LIBS)

