#######################################################################
#            Define User Specific Output Paths
#######################################################################
# Set a default directory if one has not already been defined.
STELLOPT_HOME ?= $(HOME)/bin

#######################################################################
#            Define Debian Jessie Specifec locations
#######################################################################
SYSTEM_HOME = /usr
GCC4_HOME = $(SYSTEM_HOME)
MPI_HOME = $(SYSTEM_HOME)
BLAS_HOME = $(SYSTEM_HOME)
BLACS_HOME = $(SYSTEM_HOME)
SCALAPACK_HOME = $(SYSTEM_HOME)
NAG_HOME =
NETCDF_HOME = /usr
NTCC_HOME = 
HDF5_HOME = /usr/lib/x86_64-linux-gnu/hdf5/serial/
PGPLOT_HOME =
SILO_HOME =
GLS_HOME =

GENE_HOME =

#######################################################################
#            Define Basic Utilities
#######################################################################
  # Temporaryily copy STELLOPT_HOME to MYHOME since MYHOME is currently
  # used in all the makefiles.
  MYHOME = $(STELLOPT_HOME)
  
  SHELL = /bin/sh
  PWD1 = `pwd`
  PRECOMP:= /usr/bin/cpp -traditional -DLINUX
  COMPILE_C = gcc
  COMPILE = gfortran
  COMPILE_FREE = $(COMPILE) -ffree-form -ffree-line-length-none -ffixed-line-length-none
  LINK    = ld $(FLAGS) -o
  LINK_AR = ar -ruv
  LINK_C  = $(COMPILE_C) -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
#  FLAGS_R = -O2 -g -fexternal-blas -fbacktrace -fcheck=all 
  FLAGS_R = -g -O2 -fexternal-blas -fbacktrace
  FLAGS_D = -g -O0 -fexternal-blas -fcheck=all 
  LIBS    = -L$(BLAS_HOME)/lib -lopenblas

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = $(MPI_HOME)/bin/mpif90
  MPI_COMPILE_FREE = $(MPI_COMPILE) -ffree-form \
                     -ffree-line-length-none -ffixed-line-length-none
  MPI_COMPILE_C = $(MPI_HOME)/bin/mpicc
  MPI_LINK = $(MPI_HOME)/bin/mpif90
  MPI_RUN = mpiexec
  MPI_RUN_OPTS = -np 4
#  MPI_RUN_OPTS_SM = -np 16
#  MPI_RUN_OPTS_MD = -np 64
#  MPI_RUN_OPTS_LG = -np 256

#######################################################################
#            GCC Options
#######################################################################
  GCC4_LIB = $(GCC4_HOME)/lib/x86_64-linux-gnu/

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB = -L$(NAG_HOME)/lib -lnag_nag

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = -I$(NETCDF_HOME)/include
  NETCDF_LIB = -L$(NETCDF_HOME)/lib -lnetcdf -lnetcdff

#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = F
  NTCC_INC = -I$(NTCC_HOME)/mod
  NTCC_LIB = -L$(NTCC_HOME)/lib -laladdinsub -lr8slatec -ladpak\
             -lcppsub -lcomput -lpspline -lportlib -lezcdf -lmds_sub \
             -lmdstransp -lvaxonly

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  HDF5_INC = -I$(HDF5_HOME)/include
  HDF5_LIB = -L$(HDF5_HOME)/lib -lhdf5_hl -lhdf5_fortran -lhdf5 \
             -lhdf5hl_fortran -lpthread -lz -lm

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT = F
  PGPLOT_INC = -I$(PGPLOT_HOME)
  PGPLOT_LIB = -L$(PGPLOT_HOME) -lpgplot -L/usr/X11R6/lib -lX11

#######################################################################
#             SILO Options
#######################################################################
  LSILO = F
  SILO_INC = -I$(SILO_HOME)/include
  SILO_LIB = -L$(SILO_HOME)/lib -lsilo

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
  GENE_INC = -I$(GENE_HOME)
  GENE_DIR = $(GENE_HOME)
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
                -L$(GSLHOME)/lib -lgsl -lgslcblas -lstdc++ -lmpi_cxx

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
#            LIBSTELL Shared Options
#######################################################################
LIB_SHARE = $(BLAS_HOME)/lib/libopenblas.so \
          $(SCALAPACK_HOME)/lib/libscalapack.a \
          $(BLACS_HOME)/lib/libmpiblacs.a  $(BLACS_HOME)/lib/libmpiblacsCinit.a $(BLACS_HOME)/lib/libmpiblacsF77init.a \
          $(NTCC_HOME)/lib/libaladdinsub.a $(NTCC_HOME)/lib/libr8slatec.a $(NTCC_HOME)/lib/libadpak.a $(NTCC_HOME)/lib/libcppsub.a $(NTCC_HOME)/lib/libcomput.a \
          $(NTCC_HOME)/lib/libpspline.a $(NTCC_HOME)/lib/libportlib.a $(NTCC_HOME)/lib/libezcdf.a $(NTCC_HOME)/lib/libmds_sub.a $(NTCC_HOME)/lib/libmdstransp.a $(NTCC_HOME)/lib/libvaxonly.a \
          $(HDF5_HOME)/lib/libhdf5_hl.so $(HDF5_HOME)/lib/libhdf5_fortran.so $(HDF5_HOME)/lib/libhdf5hl_fortran.so $(HDF5_HOME)/lib/libhdf5.so \
          $(NAG_HOME)/lib/libnag_nag.so  \
          $(NETCDF_HOME)/lib/libnetcdf.so $(NETCDF_HOME)/lib/libnetcdff.so $(NETCDF_HOME)/lib/libnetcdf_c++.so \
          $(SILO_HOME)/lib/libsilo.a \
          $(GSL_HOME)/lib/libgsl.so \
          $(GCC4_LIBlibgfortran.so $(GCC4_LIB/libstdc++.so \
          $(MPI_HOME)/lib/libmpi.so $(MPI_HOME)/lib/libmpi_f77.so $(MPI_HOME)/lib/libmpi_f90.so
          
