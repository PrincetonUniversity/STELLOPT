#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = $(STELLOPT_PATH)/bin
  PRECOMP:= fpp -C -P -DLINUX
  COMPILE = ifort
  COMPILE_FREE = ifort
  LINK    = ld $(FLAGS) -o
  LINK_AR = ar -ruv
  LINK_C  = $(CC) -shared -Wl,--no-undefined -Wl,-export-dynamic 

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -mcmodel=large -O2 -m64 -unroll0 -fno-alias -ip -traceback \
	-init=zero,arrays -I${MKLROOT}/include -mkl
  FLAGS_D = -O0 -g -traceback -check all -check bounds \
	-check noarg_temp_created -check uninit \
	-ftrapuv -init=zero,arrays -debug all -I${MKLROOT}/include -mkl 

  LIBS     =  ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a \
              ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -lpthread -lm -ldl \
	      -L$(SCALAPACK_HOME)/lib -lscalapack 

# BLACS requires full links to library files
  LIBS    += ${BLACS_HOME}/lib/blacs_MPI-LINUX-0.a
  LIBS    += ${BLACS_HOME}/lib/blacsF77init_MPI-LINUX-0.a
  LIBS    += ${BLACS_HOME}/lib/blacsCinit_MPI-LINUX-0.a

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpifort
  MPI_COMPILE_FREE = mpifort
  MPI_COMPILE_C = mpicc
  MPI_LINK = mpifort -shared -fpic
  MPI_RUN = srun
  MPI_RUN_OPTS = -n 64 -p general -t 5:00:00 -J STELLOPT_BENCH

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB = -L$(NAG_ROOT)/lib -lnag_nag

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = -I$(NETCDF_FORTRAN_HOME)/include -I$(NETCDF_C_HOME)/include 
  NETCDF_LIB = -L$(NETCDF_C_HOME)/lib -lnetcdf -L$(NETCDF_FORTRAN_HOME)/lib -lnetcdff

#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = T
  NTCC_INC = -I$(NTCC_HOME)/mod -I$(PSPLINE_HOME)/mod
  NTCC_LIB = -L$(NTCC_HOME)/lib -laladdinsub -lr8slatec -ladpak\
             -lcppsub -lcomput -lportlib -lezcdf -lmds_sub \
             -lmdstransp -lvaxonly -L$(PSPLINE_HOME)/lib -lpspline

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  HDF5_INC =  -I$(HDF5_HOME)/include
  HDF5_LIB = -L$(HDF5_HOME)/lib -lhdf5_hl -lhdf5_fortran -lhdf5 \
             -lhdf5hl_fortran -lpthread -lz -lm

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT = F
  PGPLOT_INC = 
  PGPLOT_LIB = -L$(PGPLOT_DIR) -lpgplot -L/usr/X11R6/lib -lX11

#######################################################################
#             SILO Options
#######################################################################
  LSILO = F
  SILO_INC = -I$(SILOHOME)/include
  SILO_LIB = -L$(SILOHOME)/lib -lsilo

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
             -L$(FFTW_HOME)/lib -lfftw3 \
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
                -L$(GSLHOME)/lib -lgsl -lgslcblas -lstdc++

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
  LTRAVIS= F
  TRAVIS_DIR = $(TRAVIS_PATH)
  LIB_TRAVIS = libtravis64_sopt.a
  LIB_MCONF  = libmconf64.a
  TRAVIS_LIB = $(TRAVIS_DIR)/lib/$(LIB_TRAVIS) \
               $(TRAVIS_DIR)/mag_conf/lib/$(LIB_MCONF) -lstdc++

#######################################################################
#            REGCOIL Options
#######################################################################
  LREGCOIL= T
  REGCOIL_DIR = $(REGCOIL_PATH)
  REGCOIL_INC = -I$(REGCOIL_DIR) 
  LIB_REGCOIL = libregcoil.a
  REGCOIL_LIB = $(REGCOIL_DIR)/$(LIB_REGCOIL) -fopenmp

#######################################################################
#            SFINCS Options
#######################################################################

  LSFINCS = F
  SFINCS_DIR = $(SFINCS_PATH)
  SFINCS_INC = -I$(SFINCS_DIR)
  LIB_SFINCS = libsfincs.a
  SFINCS_LIB = $(SFINCS_DIR)/$(LIB_SFINCS) \
             -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lX11

#######################################################################
#            LIBSTELL Shared Options
#######################################################################
LIB_SHARE = -lc -lgfortran -lstdc++ -lmpi -lmpi_mpifh -lz -lc -lm -lpthread $(LIBS) -lc
