#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = /home/weir/src/stellopt/bin
  FC = gfortran
  CC = gcc
  LK = ld
#  LK = $(FC)
  FLAGS = 
  FLAGS += -Wl,--enable-auto-import 
  LKFLAGS = -no-undefined --enable-runtime-pseudo-reloc
  PRECOMP:= cpp -traditional -DWIN64
  COMPILE = $(FC)
  COMPILE_FREE = $(FC) -ffree-form -ffree-line-length-none -ffixed-line-length-none
  LINK    = $(LK) $(LKFLAGS) $(FLAGS) -o
  LINK_AR = ar -ruv
  LINK_C  = $(CC) -shared -Wl,-z-defs
  LINK_C  += $(FLAGS) -static-libgcc -static-libstdc++

# Note:  -Wl,--z-defs causes the compiler to fail immediately if there are undefined symbols in a library
#        -Wl,--allow-multiple-definition will miss double definition errors in libraries
#        -Wl,--enable-auto-import - allow omittion of decorations in code for dllimport/dllexport
#        -Wl,--enable-runtime-pseudo-relocs -  allow relinking for structs/classes/etc in object files 
# -Wl,--export-all-symbols 
# -Wl,--whole-archive

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -O2 -fexternal-blas 
  FLAGS_D = -g -O0 -fexternal-blas -fbacktrace -fcheck=all 
  LIBS    = -L/mingw64/lib -lscalapack -lopenblas

##  LIBS    = -L/mingw64/bin/libblas.dll \
##            -L/mingw64/bin/liblapack.dll 
#  LIBS    = -L/mingw64/bin/libopenblas.dll 
##  LIBS    += -L/mingw64/bin -lmsmpi

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = F
  MPI_COMPILE = mpifort
  MPI_COMPILE_FREE = mpifort -ffree-form -fno-range-check \
                     -ffree-line-length-none -ffixed-line-length-none
  MPI_COMPILE_C = mpicc
  MPI_LINK = mpifort
  MPI_RUN = mpiexec
  MPI_RUN_OPTS = -np 4
#  MPI_RUN_OPTS_SM = -np 16
#  MPI_RUN_OPTS_MD = -np 64
#  MPI_RUN_OPTS_LG = -np 256

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB = 

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NCDIR = /mingw64
  NETCDF_INC = -I$(NCDIR)/include
  NETCDF_LIB = -L$(NCDIR)/lib -lnetcdff -lnetcdf

##  NETCDF_INC = $(shell nc-config --fflags)
##  NETCDF_LIB = $(shell nc-config --flibs)
#  NETCDF_INC = -I$(NCDIR)/include
#  NETCDF_LIB = \
#             -L$(NCDIR)/bin/libnetcdff-7.dll \
#             -L$(NCDIR)/bin/libnetcdf-15.dll

#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = F
  NTCC_INC = 
  NTCC_LIB = 

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  H5DIR = /mingw64
  HDF5_INC = -I$(H5DIR)/include \
             -I$(H5DIR)/include/static
#             -I$(H5DIR)/include/shared
  HDF5_LIB = -L$(H5DIR)/lib -lhdf5 -lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran

#  HDF5_LIB = \
#             -L$(H5DIR)/bin/libhdf5-0.dll \
#             -L$(H5DIR)/bin/libhdf5_hl_fortran-0.dll \
#             -L$(H5DIR)/bin/libhdf5_hl-0.dll \
#             -L$(H5DIR)/bin/libhdf5_fortran-0.dll

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
  SILOHOME = 
  SILO_INC = 
  SILO_LIB = 

#######################################################################
#            FFTW3 Options
#######################################################################
  LFFTW3 = T
  FFTWHOME = /mingw64
  FFTW3_INC = -I$(FFTWHOME)/include
  FFTW3_LIB = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f 

#  FFTW3_LIB = \
#              -L$(FFTWHOME)/bin/libfftw3-3.dll
#              -L$(FFTWHOME)/bin/libfftw3f-3.dll

#######################################################################
#            DKES/NEO Options
#######################################################################
  LDKES = F
  LNEO  = F

#######################################################################
#            GENE Options
#######################################################################
  LGENE = F
  GENE_INC = 
  GENE_DIR = 
  LIB_GENE = 
  GENE_LIB = 

#######################################################################
#            COILOPT++ Options
#######################################################################
  LCOILOPT = F
  GSLHOME = 
  COILOPT_INC = 
  COILOPTPP_DIR = 
  LIB_COILOPTPP = 
  COILOPT_LIB = 

#######################################################################
#            TERPSICHORE Options
#######################################################################
  LTERPSICHORE= F
  TERPSICHORE_INC = 
  TERPSICHORE_DIR = 
  LIB_TERPSICHORE = 
  TERPSICHORE_LIB = 

#######################################################################
#            TRAVIS Options
#######################################################################
  LTRAVIS= F
  TRAVIS_DIR = /home/weir/src/TRAVIS_develop/trunk/TRAVIS-src
  LIB_TRAVIS = libtravis64_sopt.a
  LIB_MCONF  = libmconf64.a
  TRAVIS_LIB = -L/$(TRAVIS_DIR)/lib -ltravis64_sopt \
               -L/$(TRAVIS_DIR)/mag_cnf/lib -lmconf64 

#  TRAVIS_LIB = -L/mingw64/bin/libstdc++-6.dll \
#               -L$(TRAVIS_DIR)/mag_conf/lib/$(LIB_MCONF) \
#               -L$(TRAVIS_DIR)/lib/$(LIB_TRAVIS)

#######################################################################
#            REGCOIL Options
#######################################################################
  LREGCOIL= F
  REGCOIL_DIR = 
  REGCOIL_INC = 
  LIB_REGCOIL = 
  REGCOIL_LIB = 

#######################################################################
#            LIBSTELL Shared Options
#######################################################################
#          /mingw64/bin/libblas.dll \
#          /mingw64/bin/liblapack.dll \

LIB_SHARE = \
          /mingw64/bin/libopenblas.dll \
          /mingw64/bin/libscalapack.dll \
          $(H5DIR)/bin/libhdf5-0.dll \
          $(H5DIR)/bin/libhdf5_hl_fortran-0.dll \
          $(H5DIR)/bin/libhdf5_hl-0.dll \
          $(H5DIR)/bin/libhdf5_fortran-0.dll \
          $(NCDIR)/bin/libnetcdf-15.dll \
          $(NCDIR)/bin/libnetcdff-7.dll

#          $(NCDIR)/bin/libnetcdf_c++4-1.dll \
#          $(TRAVIS_DIR)/mag_conf/lib/$(LIB_MCONF) \
#          $(TRAVIS_DIR)/lib/$(LIB_TRAVIS)

#          /mingw64/bin/libgfortran-5.dll \
#          /mingw64/bin/libstdc++-6.dll \
#          $(FFTWHOME)/bin/libfftw3-3.dll \
#          /mingw64/bin/libstdc++-6.dll \
#          /mingw64/bin/libgsl-25.dll \
#          /mingw64/x86_64-w64-mingw32/lib/libmsvcrt.a \
#          /mingw64/x86_64-w64-mingw32/lib/libmsvcrt-os.a 