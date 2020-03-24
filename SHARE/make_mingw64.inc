#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = /home/weir/bin
  FC = gfortran
  CC = gcc
  PRECOMP:= /mingw64/bin/cpp -traditional -DWIN
  COMPILE = $(FC)
  COMPILE_FREE = $(FC) -ffree-form -ffree-line-length-none -ffixed-line-length-none
  LINK    = ld $(FLAGS) -o
  LINK_AR = ar -ruv
  LINK_C  = $(CC) -shared -Wl,-z-defs -Wl,--export-all-symbols -Wl,--enable-auto-import, -Wl,--whole-archive
#  LINK_C  = $(CC) -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -O2 -fexternal-blas
  FLAGS_D = -g -O0 -fexternal-blas -fbacktrace -fcheck=all 
  LIBS    = -L/mingw64/lib -L/mingw64/bin/libblas.dll -L/mingw64/bin/liblapack.dll
#  LIBS    = -L/mingw64/lib -L/mingw64/bin/libopenblas.dll

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = F
  MPI_COMPILE = mpif90
  MPI_COMPILE_FREE = mpif90 -ffree-form \
                     -ffree-line-length-none -ffixed-line-length-none
  MPI_COMPILE_C = mpicc
  MPI_LINK = mpif90
  MPI_RUN = mpiexec
  MPI_RUN_OPTS = -np 4
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
  NETCDF_INC = -I/mingw64/include
  NETCDF_LIB = -L/mingw64/lib -L/mingw64/bin/libnetcdf-15.dll -L/mingw64/bin/libnetcdff-7.dll

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
  H5DIR = /mingw64
  HDF5_INC = -I$(H5DIR)/include
  HDF5_LIB = -L$(H5DIR)/lib -L/mingw64/bin/libhdf5-0.dll \
             -L/mingw64/bin/libhdf5_hl_fortran-0.dll \
             -L/mingw64/bin/libhdf5_hl-0.dll \
             -L/mingw64/bin/libhdf5_fortran-0.dll

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT = F
  PGPLOT_INC = -I$(PGPLOT_DIR)
  PGPLOT_LIB = -L$(PGPLOT_DIR) -lpgplot -L/usr/lib -lX11

#######################################################################
#             SILO Options
#######################################################################
  LSILO = F
  SILOHOME = /home/weir
  SILO_INC = -I$(SILOHOME)/include
  SILO_LIB = -L$(SILOHOME)/lib -lsiloh5

#######################################################################
#            FFTW3 Options
#######################################################################
  LFFTW3 = T
  FFTWHOME = /mingw64
  FFTW3_INC = -L$(FFTWHOME)/include
  FFTW3_LIB = -L$(FFTWHOME)/bin

#######################################################################
#            DKES/NEO Options
#######################################################################
  LDKES = F
  LNEO  = F

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
  GSLHOME = \mingw64
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
#            LIBSTELL Shared Options
#######################################################################
#LIB_SHARE = $(BLASHOME)/lib/libblas.so \
#          $(SCALAPACK_HOME)/lib/libscalapack-openmpi.so \
#          $(BLACS_HOME)/lib/libblacs-openmpi.so  $(BLACS_HOME)/lib/libblacsCinit-openmpi.so $(BLACS_HOME)/lib/libblacsF77init-openmpi.so \
#          $(HDF5_HOME)/lib/x86_64-linux-gnu/libhdf5_hl.so $(HDF5_HOME)/lib/x86_64-linux-gnu/libhdf5_fortran.so $(HDF5_HOME)/lib/x86_64-linux-gnu/libhdf5hl_fortran.so $(HDF5_HOME)/lib/x86_64-linux-gnu/libhdf5.so \
#          $(HDF5_HOME)/lib/x86_64-linux-gnu/libmpi_usempif08.so $(HDF5_HOME)/lib/x86_64-linux-gnu/libhdf5_openmpi_fortran.so \
#          $(NETCDF_HOME)/lib/x86_64-linux-gnu/libnetcdf.so $(NETCDF_HOME)/lib/x86_64-linux-gnu/libnetcdff.so $(NETCDF_HOME)/lib/x86_64-linux-gnu/libnetcdf_c++.so \
#          $(SILOHOME)/lib/x86_64-linux-gnu/libsiloh5.so \
#          $(GSLHOME)/lib/x86_64-linux-gnu/libgsl.so \
#          $(GCC6_HOME)/libgfortran.so $(GCC6_HOME)/libstdc++.so \
#          $(MPIHOME)/lib/x86_64-linux-gnu/libmpi.so $(MPIHOME)/lib/x86_64-linux-gnu/libmpi_mpifh.so \
#          /usr/lib/x86_64-linux-gnu/libm.so /usr/lib/liblapack.so /usr/lib/x86_64-linux-gnu/libdl.so
LIB_SHARE = /mingw64/bin/libblas.dll \
          /mingw64/bin/libscalapack.dll \
          /mingw64/bin/libhdf5_fortran-0.dll \
          /mingw64/bin/libnetcdf-15.dll \
          /mingw64/bin/libnetcdff-7.dll \
          /mingw64/bin/libnetcdf_c++-1.dll \
          /mingw64/bin/libgsl-25.dll \
          /mingw64/bin/libgfortran-5.dll \
          /mingw64/bin/libstdc++-6.dll \
          /mingw64/bin/liblapack.dll


