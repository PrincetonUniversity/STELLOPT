#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = $(HOME)/bin
  PRECOMP:= fpp -C -P -DLINUX
  COMPILE = ifort
  COMPILE_FREE = ifort
  LINK    = ftn -Bstatic -o
  LINK_AR = ar -ruv
  LINK_C  = icc -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  DIAG_DISABLE = -diag-disable=15009
  FLAGS_R = -O2 -fp-model strict -ip -assume noold_unit_star \
		-axMIC-AVX512 $(DIAG_DISABLE) \
		-I$(MKL_LIB)
  FLAGS_D = -g -traceback -O0 -fp-model strict -ip \
		-assume noold_unit_star -axMIC-AVX512 $(DIAG_DISABLE) \
		-I$(MKL_INC)
  LIBS    = -L$(MKL_LIB)  -lmkl_intel_lp64 -lmkl_core -lmkl_sequential \
		-lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpiifort
  MPI_COMPILE_FREE = mpiifort
  MPI_COMPILE_C = mpiifort
  MPI_LINK = mpiifort
  MPI_RUN  = srun
  MPI_RUN_OPTS = --partition=skl_fua_dbg \
		--time=0:30:00
  MPI_RUN_OPTS_SM   = $(MPI_RUN_OPTS) --nodes=1 --ntasks-per-node=24
  MPI_RUN_OPTS_MD   = $(MPI_RUN_OPTS) --nodes=2 --ntasks-per-node=48
  MPI_RUN_OPTS_LG   = $(MPI_RUN_OPTS) --nodes=4 --ntasks-per-node=48

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB = $(NAGFLIB) $(NAGMKLLIB)

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = $(shell nc-config --fflags)
  NETCDF_LIB = $(shell nc-config --flibs)

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  HDF5_INC = -I$(HDF5_HOME)/include 
  HDF5_LIB = -Wl,-rpath,$(HDF5_HOME)/lib -L$(HDF5_HOME)/lib \
             -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran 

#######################################################################
#            FFTW3 Options
#######################################################################
  LFFTW3 = T
  FFTW3_INC = -I$(FFTW_HOME)/include
  FFTW3_LIB = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f 

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
  LGENE = F
  GENE_DIR = $(GENE_PATH)/bin/obj_cobra
  FUTILS_DIR = $(GENE_PATH)/external/cobra/futils-gene-mod/src
  LIB_GENE = libgene.a
  LIB_FUTILS = libfutils.a
  GENE_INC = -I$(GENE_DIR) \
             -I$(FUTILS_DIR)
  GENE_LIB = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f \
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
  LTERPSICHORE= T
  TERPSICHORE_INC = -I$(TERPSICHORE_PATH)
  TERPSICHORE_DIR = $(TERPSICHORE_PATH)
  LIB_TERPSICHORE = libterpsichore.a
  TERPSICHORE_LIB = $(TERPSICHORE_DIR)/$(LIB_TERPSICHORE)

#######################################################################
#            TRAVIS Options
#######################################################################
  LTRAVIS= T
  TRAVIS_DIR = $(TRAVIS_PATH)
  LIB_TRAVIS = libtravis64_sopt.a
  LIB_MCONF  = libmconf64.a
  TRAVIS_LIB = $(TRAVIS_DIR)/lib/$(LIB_TRAVIS) \
               $(TRAVIS_DIR)/magconf/lib/$(LIB_MCONF) -lstdc++

#######################################################################
#            SFINCS Options
#######################################################################

  LSFINCS = F
  SFINCS_DIR = $(SFINCS_PATH)
  SFINCS_INC = -I$(SFINCS_DIR) -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
  LIB_SFINCS = libsfincs.a

  SFINCS_LIB = $(SFINCS_DIR)/$(LIB_SFINCS) \
        $(PETSC_LIB) $(FC_LINKER_FLAGS) $(FFLAGS)

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
  LIB_SHARE = $(LIBS) $(NETCDF_LIB) $(HDF5_LIB) 

