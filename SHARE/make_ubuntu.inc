#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  ARCH := $(shell uname -m)
  PWD1 = `pwd`
  MYHOME = $(HOME)/bin
  PRECOMP:= /lib/cpp -traditional -DLINUX
  COMPILE = $(FC)
  COMPILE_FREE = $(FC) -ffree-form -ffree-line-length-none -ffixed-line-length-none
  LINK    = ld $(FLAGS) -o
  LINK_AR = ar -ruv
  LINK_C  = $(CC) -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -O2 -march=native -fexternal-blas -fbacktrace -fcheck=all,no-array-temps
  FLAGS_D = -g -O0 -march=native -fexternal-blas -fbacktrace -fcheck=all
  LIBS    = -L/usr/lib -lblas -llapack -lscalapack-openmpi
  FLAGS_R += -Wno-argument-mismatch
  FLAGS_D += -Wno-argument-mismatch

  # Need to allow argument mismatch for gfortran version 10 and above
  GFORTRAN_VERSION := $(shell $(FC) -dumpversion)
  ifeq "$(shell [ $(GFORTRAN_VERSION) -gt 9 ] && echo true)" "true"
      FLAGS_R += -fallow-argument-mismatch
      FLAGS_D += -fallow-argument-mismatch
  endif

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpif90.openmpi
  MPI_COMPILE_FREE = mpif90 -ffree-form \
                     -ffree-line-length-none -ffixed-line-length-none
  MPI_COMPILE_C = mpicc.openmpi
  MPI_LINK = mpif90.openmpi
  MPI_RUN = mpiexec
  MPI_RUN_OPTS = --use-hwthread-cpus
#  MPI_RUN_OPTS_SM = -np 16
#  MPI_RUN_OPTS_MD = -np 64
#  MPI_RUN_OPTS_LG = -np 256

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
  NETCDF_LIB = $(shell nf-config --flibs)

#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = F
  NTCC_INC = -I/usr/include
  NTCC_LIB = -L/usr/lib -laladdinsub -lr8slatec -ladpak\
             -lcppsub -lcomput -lpspline -lportlib -lezcdf -lmds_sub \
             -lmdstransp -lvaxonly

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  HDF5_INC = -I/usr/include/hdf5/openmpi
  HDF5_LIB = -L/usr/lib/$(ARCH)-linux-gnu/hdf5/openmpi \
             -lhdf5 -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT = F
  PGPLOT_INC = -I$(PGPLOT_DIR)
  PGPLOT_LIB = -L$(PGPLOT_DIR) -lpgplot -L/usr/lib/$(ARCH)-linux-gnu -lX11

#######################################################################
#             SILO Options
#######################################################################
  LSILO = F
  SILO_INC = -I$(SILOHOME)/include
  SILO_LIB = -L$(SILOHOME)/lib/$(ARCH)-linux-gnu -lsiloh5

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
             -L/u/slazerso/src/GENE17_2016/external/pppl_cluster/futils/src -lfutils \
             -L$(FFTWHOME)/lib -lfftw3 \
             -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc \
             -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lX11

#######################################################################
#            COILOPT++ Options
#######################################################################
  LCOILOPT = F
  COILOPT_INC = -I$(COILOPT_PATH)
  COILOPTPP_DIR = $(COILOPT_PATH)
  LIB_COILOPTPP = libcoilopt++.a
  COILOPT_LIB = $(COILOPT_PATH)/$(LIB_COILOPTPP) \
                -L$(GSLHOME)/lib/$(ARCH)-linux-gnu -lgsl -lgslcblas -lstdc++ -lmpi_cxx

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
#            Available Energy Options
#######################################################################
  LAEOPT= F
  AEOPT_DIR = $(AEOPT_PATH)
  AEOPT_INC = -I$(AEOPT_DIR)
  LIB_AEOPT = libtrapAE.a
  AEOPT_LIB = $(AEOPT_PATH)/$(LIB_AEOPT)

#######################################################################
#            System-specific path and binary definitions
#######################################################################
FC = gfortran
BLASHOME = /usr/lib/$(ARCH)-linux-gnu
BLACS_HOME = /usr/lib/$(ARCH)-linux-gnu
SCALAPACK_HOME = /usr/lib/$(ARCH)-linux-gnu
MPIHOME = /usr/lib/$(ARCH)-linux-gnu
NETCDF_HOME = /usr/lib/$(ARCH)-linux-gnu
HDF5_HOME = /usr/lib/$(ARCH)-linux-gnu
PGPLOT_DIR = /usr/lib/$(ARCH)-linux-gnu
SILOHOME = /usr/lib/$(ARCH)-linux-gnu

#######################################################################
#            LIBSTELL Shared Options
#######################################################################
#LIB_SHARE = $(BLASHOME)/lib/$(ARCH)-linux-gnu/libblas.so \
          $(SCALAPACK_HOME)/lib/libscalapack-openmpi.so \
          $(BLACS_HOME)/lib/libblacs-openmpi.so  $(BLACS_HOME)/lib/libblacsCinit-openmpi.so $(BLACS_HOME)/lib/libblacsF77init-openmpi.so \
          $(HDF5_HOME)/lib/$(ARCH)-linux-gnu/libhdf5_hl.so $(HDF5_HOME)/lib/$(ARCH)-linux-gnu/libhdf5_fortran.so $(HDF5_HOME)/lib/$(ARCH)-linux-gnu/libhdf5hl_fortran.so $(HDF5_HOME)/lib/$(ARCH)-linux-gnu/libhdf5.so \
          $(HDF5_HOME)/lib/$(ARCH)-linux-gnu/libmpi_usempif08.so $(HDF5_HOME)/lib/$(ARCH)-linux-gnu/libhdf5_openmpi_fortran.so \
          $(NETCDF_HOME)/lib/$(ARCH)-linux-gnu/libnetcdf.so $(NETCDF_HOME)/lib/$(ARCH)-linux-gnu/libnetcdff.so $(NETCDF_HOME)/lib/$(ARCH)-linux-gnu/libnetcdf_c++.so \
          $(SILOHOME)/lib/$(ARCH)-linux-gnu/libsiloh5.so \
          $(GSLHOME)/lib/$(ARCH)-linux-gnu/libgsl.so \
          $(GCC6_HOME)/libgfortran.so $(GCC6_HOME)/libstdc++.so \
          $(MPIHOME)/lib/$(ARCH)-linux-gnu/libmpi.so $(MPIHOME)/lib/$(ARCH)-linux-gnu/libmpi_mpifh.so \
          /usr/lib/$(ARCH)-linux-gnu/libm.so /usr/lib/liblapack.so /usr/lib/$(ARCH)-linux-gnu/libdl.so
#LIB_SHARE = $(BLASHOME)/libblas.so \
          $(SCALAPACK_HOME)/libscalapack-openmpi.so.1 \
          $(BLACS_HOME)/libblacs-openmpi.so.1 \
          $(BLACS_HOME)/libblacsCinit-openmpi.so.1 \
          $(BLACS_HOME)/libblacsF77init-openmpi.so.1 \
          $(HDF5_HOME)/$(ARCH)-linux-gnu/libhdf5_openmpi_fortran.so \
          $(NETCDF_HOME)/libnetcdf.so \
          $(NETCDF_HOME)/libnetcdff.so \
          $(NETCDF_HOME)/libnetcdf_c++.so \
          $(SILOHOME)/libsiloh5.so \
          $(GSLHOME)/libgsl.so \
          $(GCC49_HOME)/libgfortran.so \
          $(GCC49_HOME)/libstdc++.so \
          $(MPIHOME)/lib/libmpi.so \
          $(MPIHOME)/lib/libmpif77.so \
          /usr/lib/$(ARCH)-linux-gnu/libm.so \
          /usr/lib/liblapack.so \
          /usr/lib/$(ARCH)-linux-gnu/libdl.so

LIB_SHARE = $(BLASHOME)/libblas.so
