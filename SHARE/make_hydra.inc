#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = /hydra/u/sul/bin
  PRECOMP:= fpp -C -P -DLINUX
  COMPILE = ifort
  COMPILE_FREE = ifort
  LINK    = ftn -Bstatic -o
  LINK_AR = ar -ruv
  LINK_C  = $(CC) -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -I$(MKL_HOME)/include -O2 -fp-model strict -xAVX -ip \
            -assume noold_unit_star 
  FLAGS_D = -g -I$(MKL_HOME)/include
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
  MPI_COMPILE_C = 
  MPI_LINK = mpiifort

#######################################################################
#            NAG Options
#######################################################################
  LNAG = T
  NAG_LIB = $(NAGFLIB) $(NAGMKLLIB)

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = -I$(NETCDF_HOME)/include 
  NETCDF_LIB = -Wl,-rpath,$(NETCDF_HOME)/lib -L$(NETCDF_HOME)/lib \
               -lnetcdf -lnetcdff 

#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = F
  NTCC_INC = -I$(NTCCHOME)/mod
  NTCC_LIB = -L$(NTCCHOME)/lib -laladdinsub -ladpak -lcppsub \
             -lcomput -lportlib 

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
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
  GENE_LIB = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f -Xlinker \
             -rpath=$(FFTW_HOME)/lib \
             -Wl,-rpath,$(SLEPC_HOME)/intel14.0-mpi.ibm1308-double-cplx/lib \
             -L$(SLEPC_HOME)/intel14.0-mpi.ibm1308-double-cplx/lib -lslepc \
             -Wl,-rpath,$(PETSC_HOME)/lib -L$(PETSC_HOME)/lib -lpetsc \
             -Wl,-rpath,$(MKL_HOME)/lib/intel64 -L$(MKL_HOME)/lib/intel64 \
             -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core \
             -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm 

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
LIB_SHARE = 

