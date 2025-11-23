#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = $(STELLOPT_PATH)/bin
  PRECOMP:= /lib/cpp -traditional -DLINUX
  COMPILE = $(FC)
  COMPILE_FREE = $(FC) -ffree-form -ffree-line-length-none -ffixed-line-length-none
  LINK    = ld $(FLAGS) -o
  LINK_AR = ar -ruv
  LINK_C  = $(CC) -shared -Wl,--no-undefined -Wl,-export-dynamic 

#######################################################################
#            Define Compiler Flags
#######################################################################
  # temp add some checking flags to avoid BEAMS3d issues
  FLAGS_R = -O2 -fexternal-blas -g -fbacktrace 
  FLAGS_D = -g -O0 -fexternal-blas -fcheck=all -fbacktrace -Wextra \
   -Wtarget-lifetime -fbounds-check -ffpe-trap=zero -finit-real=snan

# Crutch for broken module vars
  NETCDF_F = $(NETCDF_FORTRAN_HOME)
  NETCDF_C = $(NETCDF_C_HOME)
  # $(HDF5_HOME) is used.

  LIBS    = -lopenblas -L$(SCALAPACK_HOME)/lib -lscalapack 

# BLACS requires full links to library files
  LIBS    += ${BLACS_HOME}/lib/blacs_MPI-LINUX-0.a
  LIBS    += ${BLACS_HOME}/lib/blacsF77init_MPI-LINUX-0.a
  LIBS    += ${BLACS_HOME}/lib/blacsCinit_MPI-LINUX-0.a

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpifort
  MPI_COMPILE_FREE = mpifort -ffree-form \
                     -ffree-line-length-none -ffixed-line-length-none
  MPI_COMPILE_C = mpicc
  MPI_LINK = mpicc -shared -pthread -Wl,--no-undefined -Wl,-export-dynamic
  MPI_RUN = srun
  MPI_RUN_OPTS = -N 1 -n 32 -p general -t 30:00 -q debug --mem 8G -J STELLOPT_BENCH
  MPI_RUN_OPTS_SM = -N 1 -n 32 -p general -t 30:00 -q debug --mem 8G -J STELLOPT_BENCH
  MPI_RUN_OPTS_MD = -N 2 -n 64 -p general -t 30:00 --mem 8G -J STELLOPT_BENCH
  MPI_RUN_OPTS_LG = -N 4 -n 128 -p general -t 30:00 --mem 8G -J STELLOPT_BENCH

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB = -L$(NAG_ROOT)/lib -lnag_nag

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = -I$(NETCDF_F)/include -I$(NETCDF_C)/include 
  NETCDF_LIB = -L$(NETCDF_C)/lib -lnetcdf -L$(NETCDF_F)/lib -lnetcdff

#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = F
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
LIB_SHARE = -lc -lgfortran -lstdc++ -lmpi -lmpi_mpifh -lz -lc -lm -lpthread $(LIBS) -lc
