#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = $(HOME)/bin
  PRECOMP:= fpp -C -P -DLINUX
  COMPILE = ifort
  COMPILE_FREE = ifort -free
  LINK    = ifort -Bstatic -o
  LINK_AR = ar -ruv
  LINK_C  = icc -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  DIAG_DISABLE = -diag-disable=15009,5425 -check noarg_temp_created
  FLAGS_R = -O2 -fp-model strict -ip -assume noold_unit_star \
		$(DIAG_DISABLE) \
		$(shell pkg-config mkl-dynamic-ilp64-seq --cflags)
  FLAGS_D = -g -traceback -O0 -fp-model strict -ip -check all -fp-stack-check \
		-assume noold_unit_star$(DIAG_DISABLE) \
		$(shell pkg-config mkl-dynamic-ilp64-seq --cflags)
  LIBS =  ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a \
    ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a \
    -L${MKLROOT}/lib/intel64 \
    -lmkl_scalapack_ilp64 -lmkl_intel_ilp64 -lmkl_sequential \
    -lmkl_core -lmkl_blacs_intelmpi_ilp64 -lpthread -lm -ldl

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = h5pfc
  MPI_COMPILE_FREE = h5pfc -free
  MPI_COMPILE_C = h5pcc
  MPI_LINK = h5pfc
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
  NETCDF_INC = $(shell pkg-config netcdf-fortran --cflags)
  NETCDF_LIB = $(shell pkg-config netcdf-fortran --libs) 

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  PRECOMP+= -DHDF5_PAR
#  HDF5_INC = -I$(HDF5_HOME)/include 
#  HDF5_LIB = -Wl,-rpath,$(HDF5_HOME)/lib -L$(HDF5_HOME)/lib \
#             -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran
#  HDF5_INC = $(shell pkg-config hdf5 --cflags)
#  HDF5_LIB = $(shell pkg-config hdf5 --libs)

#######################################################################
#            FFTW3 Options
#######################################################################
  LFFTW3 = T
  FFTW3_INC = $(shell pkg-config fftw3 --cflags)
  FFTW3_LIB = $(shell pkg-config fftw3 --libs)

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

