#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = /pfs/home/$(USER)/bin
  PRECOMP:= fpp -C -P -DLINUX
  COMPILE = ifort
  COMPILE_FREE = ifort -free
  LINK    = ifort -Bstatic -o
  LINK_AR = ar -ruv
  LINK_C  = icc -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  DIAG_DISABLE = -diag-disable=15009 -check noarg_temp_created
  FLAGS_R = -O2 -fp-model strict -ip -assume noold_unit_star \
		-axMIC-AVX512 $(DIAG_DISABLE) \
		$(shell pkg-config mkl-ifort --cflags)
  FLAGS_D = -g -traceback -O0 -fp-model strict -ip -check all -fp-stack-check \
		-assume noold_unit_star -axMIC-AVX512 $(DIAG_DISABLE) \
		$(shell pkg-config mkl-ifort --cflags)
  # This is a workaround for broken pkg-config
  MKLROOT ?= /cineca/prod/opt/compilers/intel/pe-xe-2017/binary/mkl
  LIBS =  ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a \
           -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a \
                             ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a \
                             ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group \
           -lpthread -lm -ldl

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpiifort
  MPI_COMPILE_FREE = mpiifort -free
  MPI_COMPILE_C = mpicc
  MPI_LINK = mpiifort
  MPI_RUN  = srun
  MPI_RUN_OPTS = --partition=gwdbg \
		--time=0:30:00
  MPI_RUN_OPTS_SM   = $(MPI_RUN_OPTS) --nodes=1 --ntasks-per-node=24
  MPI_RUN_OPTS_MD   = $(MPI_RUN_OPTS) --nodes=2 --ntasks-per-node=48
  MPI_RUN_OPTS_LG   = $(MPI_RUN_OPTS) --nodes=2 --ntasks-per-node=48

#######################################################################
#            IMAS Options
#######################################################################
  LIMAS = T
  IMAS_INC = $(shell pkg-config imas-ifort --cflags) \
             $(shell pkg-config xmllib --cflags)
  IMAS_LIB = $(shell pkg-config imas-ifort --libs) \
             $(shell pkg-config xmllib --libs)

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB = $(NAGFLIB) $(NAGMKLLIB)

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = $(shell pkg-config netcdf-fortran-ifort --cflags)
  NETCDF_LIB = $(shell pkg-config netcdf-fortran-ifort --libs) 

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  PRECOMP+= -DHDF5_PAR
#  HDF5_INC = -I$(HDF5_HOME)/include 
#  HDF5_LIB = -Wl,-rpath,$(HDF5_HOME)/lib -L$(HDF5_HOME)/lib \
#             -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran
  HDF5_INC = $(shell pkg-config hdf5 --cflags)
  HDF5_LIB = $(shell pkg-config hdf5 --libs)

#######################################################################
#            FFTW3 Options
#######################################################################
  LFFTW3 = T
  FFTW3_INC = $(shell pkg-config fftw3-ifort --cflags)
  FFTW3_LIB = $(shell pkg-config fftw3-ifort --libs)

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
  LIB_SHARE = $(LIBS) 

