# -*- mode: makefile -*-

#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = $(HOME)/bin
  PRECOMP:= fpp -C -P -DLINUX
  COMPILE = ifort
  COMPILE_FREE = ifort
  LINK    = icc -Bstatic -o
  LINK_AR = ar -ruv
  LINK_C  = icc -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -I$(MKL_HOME)/include -O2 -fp-model strict -ip \
            -assume noold_unit_star -g -traceback
  FLAGS_D = -I$(MKL_HOME)/include -fp-model strict -ip \
            -assume noold_unit_star -g -traceback
  FLAGS_D+= -check all -check bounds -check noarg_temp_created \
            -check uninit -ftrapuv -debug all
  LIBS    = -Wl,-rpath,$(MKL_HOME)/lib/intel64 \
            -L$(MKL_HOME)/lib/intel64 -lmkl_scalapack_lp64 \
            -lmkl_intel_lp64 -lmkl_core -lmkl_sequential \
            -lmkl_blacs_intelmpi_lp64  -lpthread -lm -lpthread

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpiifort
  MPI_COMPILE_FREE = mpiifort
  MPI_COMPILE_C = mpiifort 
  MPI_LINK = mpiifort
  MPI_RUN = srun
  MPI_RUN_OPTS = --nodes=1 --ntasks-per-node=32 --time=0:30:00 -p express
  MPI_RUN_OPTS_SM   = --nodes=1 --ntasks-per-node=32 --time=0:30:00 -p express
  MPI_RUN_OPTS_MD   = --nodes=8 --ntasks-per-node=32 --time=0:30:00 -p express
  MPI_RUN_OPTS_LG   = --nodes=16 --ntasks-per-node=32 --time=0:30:00 -p express
#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB = $(NAGFLIB) $(NAGMKLLIB)

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = -I$(NETCDF_HOME)/include 
  NETCDF_LIB = -Wl,-rpath,$(NETCDF_HOME)/lib -L$(NETCDF_HOME)/lib \
               -lnetcdf -lnetcdff 

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  PRECOMP+= -DHDF5_PAR
  HDF5_INC = -I$(HDF5_HOME)/include 
  HDF5_LIB = -Wl,-rpath,$(HDF5_HOME)/lib -L$(HDF5_HOME)/lib \
             -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran 

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
  FFTW3_INC = -I$(FFTW_HOME)/include
  FFTW3_LIB = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f 

#######################################################################
#            DKES/NEO Options
#######################################################################
  LDKES = T
  LNEO  = T

#######################################################################
#            GENE Options
#######################################################################
  LGENE = F
  GENE_DIR = $(GENE_PATH)/bin/obj_draco
  FUTILS_DIR = $(GENE_PATH)/external/draco/futils-gene-mod/src
  LIB_GENE = libgene.a
  LIB_FUTILS = libfutils.a
  GENE_INC = -I$(GENE_DIR) \
             -I$(FUTILS_DIR)
  GENE_LIB = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f -Xlinker \
             -L$(SLEPC_HOME)/lib -lslepc \
             -L$(PETSC_HOME)/lib -lpetsc \
             $(GENE_DIR)/$(LIB_GENE) \
             $(FUTILS_DIR)/$(LIB_FUTILS)

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
  LSFINCS = F
  include ${PETSC_DIR}/lib/petsc/conf/variables
  SFINCS_DIR = $(SFINCS_PATH)/fortran/version3
  SFINCS_INC = -I$(SFINCS_DIR) -I${PETSC_DIR}/include
  LIB_SFINCS = libsfincs.a

  SFINCS_LIB = $(SFINCS_DIR)/$(LIB_SFINCS) \
        $(PETSC_LIB) $(FC_LINKER_FLAGS) $(FFLAGS)

#######################################################################
#            MANGO Options
#######################################################################

LMANGO = F

MANGO_DIR = /u/mlan/mango

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
  REGCOIL_LIB = $(REGCOIL_DIR)/$(LIB_REGCOIL) -qopenmp

#######################################################################
#            LIBSTELL Shared Options
#######################################################################
  LIB_SHARE = $(LIBS) $(NETCDF_LIB) $(HDF5_LIB) \
             -lifport -lifcoremt -limf -lm -lc \
             -lpthread -lirc -lunwind -lirc_s

