#######################################################################
#     This is a custom file compiling libstell on windows under 
#     msys/mingw.  Follow the instructions under README_mscompile.md 
#     to setup the environment
#
#
# Note:  -Wl,-z-defs causes the compiler to fail immediately if there are undefined symbols in a library
#        -Wl,--allow-multiple-definition will miss double definition errors in libraries
#        -Wl,--enable-auto-import - allow omittion of decorations in code for dllimport/dllexport
#        -Wl,--enable-runtime-pseudo-relocs -  allow relinking for structs/classes/etc in object files 
# -Wl,--export-all-symbols 
# -Wl,--whole-archive
# -static-libstdc++
# -static-libgcc 
# -static-libgfortran
#
# NOTE:
#  "-Wl,--export-dynamic" is not supported for PX+ type targets (windows+mingw)
#
#######################################################################
#          If you want to explicitly override user given flags
#######################################################################
#  FC = gfortran
#  CC = gcc
#  FLAGS = -Wl,--enable-auto-import -Wl,--export-all-symbols -Wl,--enable-runtime-pseudo-reloc
#
#######################################################################
#            Define Basic Utilities
#######################################################################
  HARD_FFLAGS =  -Wl,--no-undefined \
                 -Wl,--enable-auto-import \
                 -Wl,--export-all-symbols \
                 -Wl,--enable-runtime-pseudo-reloc

  HARD_CFLAGS = -Wl,-z-defs 
#                    -Wl,--no-undefined \
#                    -Wl,--enable-auto-import \
#                    -Wl,--export-all-symbols \
#                    -Wl,--enable-runtime-pseudo-reloc

  HARD_LINK_FLAGS = -Wl,--no-undefined \
                    -Wl,-no_compact_unwind \
                    -Wl,--enable-auto-import \
                    -Wl,--export-all-symbols \
                    -Wl,--enable-runtime-pseudo-reloc

  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = /home/weir/src/stellopt/bin
  PRECOMP:= cpp -traditional -DWIN64
  COMPILE = $(FC) $(HARD_FFLAGS)
  COMPILE_FREE = $(FC) -ffree-form -ffree-line-length-none \
                       -ffixed-line-length-none $(HARD_FFLAGS)
  LINK    = ld $(FLAGS) $(HARD_LINK_FLAGS) -o
  LINK_AR = ar -ruv
  LINK_C  = $(CC) -shared $(HARD_CFLAGS)


#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -O2 -fexternal-blas 
  FLAGS_D = -g -O0 -fexternal-blas -fbacktrace -fcheck=all 

#######################################################################
#            Define LIBS: 
#            dependency order matters (highest left, lowest right) 
#
#   SCALAPACK depends on FFTW3 and can depend on MPI
#     Make sure that you set those flags appropriately depending on your build
#
#     NOTE: mingw auto-imports libm but includes an empty libm for compatibility
#
#######################################################################

  LSCALAPACK = T
  LOPENBLAS = T

ifeq ($(LSCALAPACK),T)
  ifeq ($(LOPENBLAS),T)
    LIBS    = -L/mingw64/lib -lscalapack -lopenblas
  else
    LIBS    = -L/mingw64/lib -lscalapack -llapack -lblas
  endif
else
  LIBS    = -L/mingw64/lib -llapack -lblas
endif

#  LIBS    = -L/mingw64/lib -lscalapack -lopenblas -lmpi -lwinpthread -lgfortran 
#
#######################################################################
#            FFTW3 Options
#######################################################################
  LFFTW3 = T

ifeq ($(LFFTW3),T)
  FFTWHOME = /mingw64
  FFTW3_INC = -I$(FFTWHOME)/include
  FFTW3_LIB = -L$(FFTWHOME)/lib -lfftw3f -lfftw3
endif

#######################################################################
#            MPI Options
#
# MPI RUN OPTIONS:
#  MPI_RUN_OPTS = -np 4
#  MPI_RUN_OPTS_SM = -np 16
#  MPI_RUN_OPTS_MD = -np 64
#  MPI_RUN_OPTS_LG = -np 256
#
# Let's not pollute the namespace if unnecessary
#######################################################################
  LMPI    = T

ifeq ($(LMPI),T)
  MPI_COMPILE = mpif90 $(HARD_FFLAGS)
  MPI_COMPILE_FREE = mpif90 -ffree-form -fno-range-check \
                     -ffree-line-length-none -ffixed-line-length-none \
                     $(HARD_FFLAGS) 
  MPI_COMPILE_C = mpicc
  MPI_LINK = mpicc -shared -pthread $(HARD_LINK_FLAGS)

  MPI_RUN = mpiexec
  MPI_RUN_OPTS = -np 4

  MPI_LIB = /mingw64/lib -lmpi -lwinpthread
endif

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F

ifeq ($(LNAG),T)
  NAG_ROOT =
  NAG_LIB =
#  NAG_LIB = -L$(NAG_ROOT)/lib -lnag_nag
endif

#######################################################################
#            NETCDF Options
#
# NOTE: NETCDF depends on HDF5 (and possibly your pthread library)
#    netcdf-fortran depends on netcdf
#    netcdf-c depends on HDF5 and its libraries 
#     (but netcdf-c "knows" where to find HDF5 libraries and libz, libcurl)
#
# Optionally you can use to automatically set the dependencies / order
#    this slows down compilation on my computer
#  NETCDF_INC = $(shell nc-config --fflags)  # fortran flags
#  NETCDF_LIB = $(shell nc-config --flibs)   # fortran libs
#######################################################################
  LNETCDF = T

ifeq ($(LNETCDF),T)
  NETCDF_C_HOME = /mingw64
  NETCDF_F_HOME = /mingw64
  NETCDF_C = $(NETCDF_C_HOME)
  NETCDF_F = $(NETCDF_F_HOME)
  NETCDF_INC = -I$(NETCDF_F)/include -I$(NETCDF_C)/include 
  NETCDF_LIB = -L$(NETCDF_F)/lib -lnetcdff -L$(NETCDF_C)/lib -lnetcdf -lnetcdf -lm
endif

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T

ifeq ($(LHDF5),T)
  HDF5_HOME = /mingw64
  HDF5_INC = -I$(HDF5_HOME)/include -I$(HDF5_HOME)/include/static
#  HDF5_INC = -I$(HDF5_HOME)/include -I$(HDF5_HOME)/include/shared

  HDF5_LIB = -L$(HDF5_HOME)/lib -lhdf5_hl_fortran -lhdf5_fortran -lhdf5_hl \
             -lhdf5 -lz -lm -lcurl
endif


#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = F

ifeq ($(LNTCC),T)
  NTCC_INC = 
  NTCC_LIB =

#  NTCC_INC = -I$(NTCC_HOME)/mod -I$(PSPLINE_HOME)/mod
#  NTCC_LIB = -L$(NTCC_HOME)/lib -laladdinsub -lr8slatec -ladpak\
#             -lcppsub -lcomput -lportlib -lezcdf -lmds_sub \
#             -lmdstransp -lvaxonly -L$(PSPLINE_HOME)/lib -lpspline
#
endif

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT = F

ifeq ($(LPGPLOT),T)
  PGPLOT_DIR = 
  PGPLOT_INC = 
  PGPLOT_LIB =
#
#  PGPLOT_LIB = -L$(PGPLOT_DIR) -lpgplot -L/usr/X11R6/lib -lX11
#
endif

#######################################################################
#             SILO Options
#######################################################################
  LSILO = F

ifeq ($(LSILO),T)
  SILOHOME = 
  SILO_INC = 
  SILO_LIB =
#
#  SILO_LIB = -L$(SILOHOME)/lib -lsilo
#
endif

#######################################################################
#            DKES/NEO Options
#######################################################################
  LDKES = F
  LNEO  = T

#######################################################################
#            GENE Options
#######################################################################
  LGENE = F

ifeq ($(LGENE),T)
  GENE_INC = 
  GENE_DIR = 
  LIB_GENE = libgene.a
  GENE_LIB =

#  GENE_LIB = $(GENE_DIR)/$(LIB_GENE) \
#             -L$(FFTW_HOME)/lib -lfftw3 \
#             -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc \
#             -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lX11
#
endif

#######################################################################
#            COILOPT++ Options
#
#
#######################################################################
  LCOILOPT = F

ifeq ($(LCOILOPT),T)
  COILOPT_INC = 
  COILOPTPP_DIR = 
  LIB_COILOPTPP = libcoilopt++.a
  COILOPT_LIB =

#  GSLHOME = /mingw64
#  COILOPT_LIB = $(COILOPT_PATH)/$(LIB_COILOPTPP) \
#                -L$(GSLHOME)/lib -lgsl -lgslcblas -lstdc++
#
endif

#######################################################################
#            TERPSICHORE Options
#######################################################################
  LTERPSICHORE = F

ifeq ($(LTERPSICHORE),T)
  TERPSICHORE_INC = 
  TERPSICHORE_DIR = 
  LIB_TERPSICHORE = libterpsichore.a
  TERPSICHORE_LIB =
#
#  TERPSICHORE_LIB = $(TERPSICHORE_DIR)/$(LIB_TERPSICHORE)
#
endif

#######################################################################
#            TRAVIS Options
#
#    TRAVIS depends on the MCONF library and the standard C++ library
#    MCONF depends on the stdc++ library and libc (msvcrt)
#######################################################################
  LTRAVIS= F

ifeq ($(LTRAVIS),T)
  TRAVIS_DIR = /home/weir/src/TRAVIS_develop/trunk/TRAVIS-src
  LIB_TRAVIS = libtravis64_sopt.a
  LIB_MCONF  = libmconf64.a
  TRAVIS_LIB =
#
#  TRAVIS_LIB = -L/$(TRAVIS_DIR)/lib/$(LIB_TRAVIS) \
#               -L/$(TRAVIS_DIR)/mag_conf/lib/$(LIB_MCONF) \
#               -lstdc++ -lmsvcrt
#
endif

#######################################################################
#            REGCOIL Options
#######################################################################
  LREGCOIL= F

ifeq ($(LREGCOIL),T)
  REGCOIL_DIR = 
  REGCOIL_INC = 
  LIB_REGCOIL = libregcoil.a
  REGCOIL_LIB = 
#
#  REGCOIL_LIB = $(REGCOIL_DIR)/$(LIB_REGCOIL) -fopenmp
#
endif

#######################################################################
#            SFINCS Options
#
#   Note:  PETSC has a bunch of dependencies
#
#######################################################################

  LSFINCS = F

ifeq ($(LSFINCS),T)
  SFINCS_DIR = $(SFINCS_PATH)
  SFINCS_INC = -I$(SFINCS_DIR)
  LIB_SFINCS = libsfincs.a
  SFINCS_LIB =
#
#  SFINCS_LIB = $(SFINCS_DIR)/$(LIB_SFINCS) \
#             -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lX11
#
endif

#######################################################################
#            LIBSTELL Shared Options
#######################################################################
# Using c-compiler for link step:  need to include gfortran library
#    Recursion through libraries is necessary to resolve undefined symbols
#    Order matters! Moving left to right:
#            search LIBS  and   add items to undefined symbol / function table
#            re-search LIBS because there are recursive dependencies
#            check msvcrt (libc in Windows) resolve some and add more undefined symbols
#            check gfortran lib and resolve more undefined symbols ... add some math/compression undefined symbols
#            check zlib and resolve ... add more from gfortran / libc
#            resolve the new ones from msvcrt
#            check libm for math symbols
#            resolve the new ones from msvcrt
#            finally resolve the c-compiler undefined symbols for gfortran
#                   
#   Trying out a recursive tool in ld
#
LIB_SHARE = $(LIBS) -lz -lcurl -lmsvcrt -lm -lgfortran

#LIB_SHARE = $(LIBS) $(LIBS) -lz -lcurl -lmsvcrt -lm -lgfortran -lmsvcrt

######################################################################
######################################################################
#
#LIB_SHARE = -lmsvcrt -lgfortran -lz -lmsvcrt -lm $(LIBS) -lmsvcrt 
#
#LIB_SHARE = -lmsvcrt -lgfortran -lstdc++ -lmpi -lz -lmsvcrt -lm -lwinpthread $(LIBS) -lmsvcrt 
#
#LIB_SHARE = -lgfortran -lstdc++ $(LIBS)
#
#  LIB_SHARE = \
#            -L/mingw64/bin -lmpi \
#            -L$(LIBS) \
#            -L$(FFTW3_LIB) \
#            -L$(NTCC_LIB) \
#            -L$(HDF5_LIB) \
#            -L$(NETCDF_LIB) \
#            -L$(NAG_LIB) \
#            -L$(PGPLOT_LIB) \
#            -L$(SILO_LIB) \
#            -L$(GENE_LIB) \
#            -L$(COILOPT_LIB) \
#            -L$(TERPSICHORE_LIB) \
#            -L$(TRAVIS_LIB) \
#            -L$(REGCOIL_LIB)
#
#  LIB_SHARE = \
#          -L/mingw64/bin/libopenblas.dll \
#          -L/mingw64/bin/libscalapack.dll \
#          -L$(H5DIR)/bin/libhdf5-0.dll \
#          -L$(H5DIR)/bin/libhdf5_hl_fortran-0.dll \
#          -L$(H5DIR)/bin/libhdf5_hl-0.dll \
#          -L$(H5DIR)/bin/libhdf5_fortran-0.dll \
#          -L$(NCDIR)/bin/libnetcdf-15.dll \
#          -L$(NCDIR)/bin/libnetcdff-7.dll
#
#          $(NCDIR)/bin/libnetcdf_c++4-1.dll \
#          $(TRAVIS_DIR)/mag_conf/lib/$(LIB_MCONF) \
#          $(TRAVIS_DIR)/lib/$(LIB_TRAVIS)
#
#          -L/mingw64/bin/libgfortran-5.dll \
#          /mingw64/bin/libstdc++-6.dll \
#          $(FFTWHOME)/bin/libfftw3-3.dll \
#          /mingw64/bin/libstdc++-6.dll \
#          /mingw64/bin/libgsl-25.dll \
#          /mingw64/x86_64-w64-mingw32/lib/libmsvcrt.a \
#          /mingw64/x86_64-w64-mingw32/lib/libmsvcrt-os.a
#
#          /mingw64/bin/libblas.dll \
#          /mingw64/bin/liblapack.dll \
#            -L/mingw64/bin -lgfortran \
#            -L/mingw64/bin/libgfortran-5.dll \
#            -L/mingw64/bin -lgfortran \
# 