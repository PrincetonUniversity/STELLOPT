#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = $(HOME)/bin_stellopt_g3_gnu
  PRECOMP:= cpp -traditional -DLINUX
  COMPILE = gfortran
  COMPILE_FREE = gfortran -ffree-form -ffree-line-length-none \
                          -ffixed-line-length-none
  LINK    = gfortran -$(FLAGS) -o
  LINK_AR = ar -ruv
  LINK_C  = gcc -shared -Wl,-no_compact_unwind

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -I${MKLROOT}/include -O2 -fexternal-blas
  FLAGS_D = -I${MKLROOT}/include -O0 -fexternal-blas -g -fbacktrace \
            -fcheck=all -fbounds-check -Wall
  LIBS    =  ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a \
            -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a \
             ${MKLROOT}/lib/intel64/libmkl_sequential.a \
             ${MKLROOT}/lib/intel64/libmkl_core.a \
             ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
            -Wl,--end-group -lpthread -lm -ldl
#  FLAGS_R += -fdefault-integer-8  -m64
#  FLAGS_D += -fdefault-integer-8  -m64

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpifort
  MPI_COMPILE_FREE = mpifort -ffree-form -ffree-line-length-none \
                            -ffixed-line-length-none
  MPI_COMPILE_C = mpigcc
  MPI_LINK = mpigcc
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
  NETCDF_INC = $(shell nc-config --fflags)
  NETCDF_LIB = $(shell nc-config --flibs)

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
  FFTW3_LIB = -Wl,-rpath,$(FFTW_HOME)/lib -L$(FFTW_HOME)/lib \
              -lfftw3 -lfftw3f 

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
  GENE_LIB = $(GENE_DIR)/$(LIB_GENE) \
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
  SFINCS_DIR = $(SFINCS_PATH)/fortran/version3
  SFINCS_INC = -I$(SFINCS_DIR)
  LIB_SFINCS = libsfincs.a
  SFINCS_LIB = $(SFINCS_DIR)/$(LIB_SFINCS) 

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
#            VMECTOOLS Options
#######################################################################
  LVMECTOOLS= T
  VMECTOOLS_DIR = $(VMECTOOLS_PATH)
  VMECTOOLS_INC = -I$(VMECTOOLS_DIR) 
  LIB_VMECTOOLS = libvmectools.a
  VMECTOOLS_LIB = $(VMECTOOLS_DIR)/lib/$(LIB_VMECTOOLS)

#######################################################################
#            LIBSTELL Shared Options
#######################################################################
  LIB_SHARE = $(LIBS) $(NETCDF_LIB) $(HDF5_LIB) 

