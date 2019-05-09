#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = /bin
  PRECOMP:= /usr/bin/cpp -traditional -DLINUX
  COMPILE = $(FC)
  COMPILE_FREE = $(FC) -ffree-form -ffree-line-length-none -ffixed-line-length-none
  LINK    = ld $(FLAGS) -o
  LINK_AR = ar -ruv
  LINK_C  = $(CC) -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = --O3 -fexternal-blas
  FLAGS_D = -g -O0 -fexternal-blas -fbacktrace -fcheck=all 
  LIBS    = -lscalapack -lopenblas

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpifort
  MPI_COMPILE_FREE = mpifort -ffree-form \
                     -ffree-line-length-none -ffixed-line-length-none
  MPI_COMPILE_C = mpicc
  MPI_LINK = mpifort

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB = -L$(NAG_ROOT)/lib -lnag_nag

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = 
  NETCDF_LIB = -lnetcdf -lnetcdff

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
  HDF5_INC =
  HDF5_LIB = -lhdf5 -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT = F
  PGPLOT_INC = -I$(PGPLOT_DIR)
  PGPLOT_LIB = -L$(PGPLOT_DIR) -lpgplot -L/usr/lib/x86_64-linux-gnu -lX11

#######################################################################
#             SILO Options
#######################################################################
  LSILO = F
  SILO_INC = -I$(SILOHOME)/include
  SILO_LIB = -L$(SILOHOME)/lib/x86_64-linux-gnu -lsiloh5

#######################################################################
#            DKES/NEO Options
#######################################################################
  LDKES = F
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
                -L$(GSLHOME)/lib/x86_64-linux-gnu -lgsl -lgslcblas -lstdc++ -lmpi_cxx

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

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#   DO NOT EDIT BELOW HERE DO NOT EDIT BELOW HERE

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

# Some Initial stuff
# HOME_BIN = $(MYHOME)
HOME_BIN = ../../
MOD1_PATH = -I$(MYHOME)/libstell_dir
MOD_PATH = 
CLEAN_ARCH = 
LIB_DIR = ../../LIBSTELL
LIB = libstell.a

# Handle Which Code is being CALLED
ifeq ($(MAKECMDGOALS),xstelloptv2)
  PRECOMP += -DBNORM -DBEAMS3D_OPT -DJINV_OPT -DTXPORT_OPT -DDIAGNO_OPT -DSKS -DSKS2
  VMEC_DIR ?= ../../VMEC2000
  LIB_VMEC= libvmec.a
  MOD_PATH+= -I$(VMEC_DIR)/$(LOCTYPE)
  BEAMS3D_DIR = ../../BEAMS3D
  LIB_BEAMS3D = libbeams3d.a
  MOD_PATH+= -I$(BEAMS3D_DIR)/$(LOCTYPE)
  BOOTSJ_DIR = ../../BOOTSJ
  LIB_BOOTSJ = libbootsj.a
  MOD_PATH+= -I$(BOOTSJ_DIR)/$(LOCTYPE)
  BNORM_DIR = ../../BNORM
  LIB_BNORM = libbnorm.a
  MOD_PATH+= -I$(BNORM_DIR)/$(LOCTYPE)
  BOOZ_DIR = ../../BOOZ_XFORM
  LIB_BOOZ= libbooz.a
  MOD_PATH+= -I$(BOOZ_DIR)/$(LOCTYPE)
  COBRA_DIR = ../../COBRAVMEC
  LIB_COBRA = libcobra.a
  MOD_PATH+= -I$(COBRA_DIR)/$(LOCTYPE)
  DIAGNO_DIR = ../../DIAGNO
  LIB_DIAGNO= libdiagno.a
  MOD_PATH+= -I$(DIAGNO_DIR)/$(LOCTYPE)
  JINV_DIR = ../../J_INVARIANT
  LIB_JINV = libjinv.a
  MOD_PATH+= -I$(JINV_DIR)/$(LOCTYPE)
  MGRID_DIR= ../../MAKEGRID
  LIB_MGRID= libmakegrid.a
  MOD_PATH+= -I$(MGRID_DIR)/$(LOCTYPE)
  LIB_LINK=  $(VMEC_DIR)/$(LOCTYPE)/$(LIB_VMEC)  \
	$(BEAMS3D_DIR)/$(LOCTYPE)/$(LIB_BEAMS3D) \
	$(BOOTSJ_DIR)/$(LOCTYPE)/$(LIB_BOOTSJ) \
	$(BNORM_DIR)/$(LOCTYPE)/$(LIB_BNORM) \
	$(BOOZ_DIR)/$(LOCTYPE)/$(LIB_BOOZ) \
	$(COBRA_DIR)/$(LOCTYPE)/$(LIB_COBRA) \
	$(DIAGNO_DIR)/$(LOCTYPE)/$(LIB_DIAGNO) \
	$(JINV_DIR)/$(LOCTYPE)/$(LIB_JINV) \
	$(MGRID_DIR)/$(LOCTYPE)/$(LIB_MGRID) 
  ifeq ($(LDKES),T)
    DKES_DIR = ../../DKES
    LIB_DKES = libdkes.a
    MOD_PATH += -I$(DKES_DIR)/$(LOCTYPE)
    PRECOMP += -DDKES_OPT
    LIB_LINK += $(DKES_DIR)/$(LOCTYPE)/$(LIB_DKES)
  else
    DKES_DIR = 
    LIB_DKES = 
  endif
  ifeq ($(LNEO),T)
    NEO_DIR = ../../NEO
    LIB_NEO = libneo.a
    MOD_PATH += -I$(NEO_DIR)/$(LOCTYPE)
    PRECOMP += -DNEO_OPT
    LIB_LINK += $(NEO_DIR)/$(LOCTYPE)/$(LIB_NEO)
  else
    NEO_DIR = 
    LIB_NEO = 
  endif
else
  LGENE = F
  LCOILOPT = F
  LTERPSICHORE = F
  LTRAVIS = F
endif
ifeq ($(MAKECMDGOALS),xanimec)
  PRECOMP += -D_ANIMEC
endif
ifeq ($(MAKECMDGOALS),xparvmec)
  PRECOMP += -DSKS  -DSKS2
endif
ifeq ($(MAKECMDGOALS),xsiesta)
  PRECOMP += -DSKS  -DSKS2
endif
ifeq ($(MAKECMDGOALS),xvmec2000)
  PRECOMP += -DSKS  -DSKS2
endif
ifeq ($(MAKECMDGOALS),xgtovmi)
  ifeq ($(LPGPLOT),T)
    MOD1_PATH += $(PGPLOT_INC)
    LIB_LINK += $(PGPLOT_LIB)
  endif
endif
ifeq ($(MAKECMDGOALS),xdes_plot)
  ifeq ($(LPGPLOT),T)
    MOD1_PATH += $(PGPLOT_INC)
    LIB_LINK += $(PGPLOT_LIB)
  endif
endif

# Other Codes
ifeq ($(LGENE),T)
  MOD1_PATH += $(GENE_INC)
  PRECOMP += -DGENE
  MOD_PATH += -I$(GENE_DIR)
else
  GENE_LIB = 
  GENE_DIR = 
  LIB_GENE = 
endif
ifeq ($(LCOILOPT),T)
  MOD1_PATH += $(COILOPT_INC)
  PRECOMP += -DCOILOPTPP
  MOD_PATH += -I$(COILOPTPP_DIR)
else
  COILOPT_LIB = 
  COILOPTPP_DIR = 
  LIB_COILOPTPP = 
endif
ifeq ($(LTERPSICHORE),T)
  MOD1_PATH += $(TERPSICHORE_INC)
  PRECOMP += -DTERPSICHORE
  MOD_PATH += -I$(TERPSICHORE_DIR)
else
  TERPSICHORE_LIB = 
  TERPSICHORE_DIR = 
  LIB_TERPSICHORE = 
endif
ifeq ($(LTRAVIS),T)
  PRECOMP += -DTRAVIS
else
  TRAVIS_LIB = 
  TRAVIS_DIR = 
  LIB_TRAVIS = 
endif
ifeq ($(LREGCOIL),T)
  MOD1_PATH += $(REGCOIL_INC)
  PRECOMP += -DREGCOIL
  MOD_PATH += -I$(REGCOIL_DIR)
else
  REGCOIL_LIB = 
  REGCOIL_DIR = 
  LIB_REGCOIL = 
endif

# Handle other codes
LIB_LINK += $(GENE_LIB) $(COILOPT_LIB) $(TERPSICHORE_LIB) $(TRAVIS_LIB) $(REGCOIL_LIB)

# Add libraries to LIB_LINK
# LIB_LINK+= $(MYHOME)/libstell.a $(LIBS)
# This version avoids referencing the bin directory for the symbolic link to lib.
# Which we can't reference since we now have the awk file in the main directory.
LIB_LINK+= $(LIB_DIR)/$(LOCTYPE)/$(LIB) $(LIBS)

# 3rd Party Libs
ifeq ($(LMPI),T)
  PRECOMP += -DMPI_OPT
  VMEC_DIR = ../../VMEC2000
  COMPILE = $(MPI_COMPILE)
  COMPILE_FREE = $(MPI_COMPILE_FREE)
  LINK    = $(MPI_COMPILE) $(FLAGS) -o 
  LINK_C  = $(MPI_COMPILE_C) -shared -Wl,-z,defs
endif
ifeq ($(LNETCDF),T)
  MOD1_PATH += $(NETCDF_INC)
  PRECOMP += -DNETCDF
  LIB_LINK += $(NETCDF_LIB)
  MOD_PATH += $(NETCDF_INC)
endif
ifeq ($(LHDF5),T)
  MOD1_PATH += $(HDF5_INC)
  PRECOMP += -DLHDF5 -DH5_USE_16_API
  LIB_LINK += $(HDF5_LIB)
  MOD_PATH += $(HDF5_INC)
endif
ifeq ($(LNAG),T)
  PRECOMP += -DNAG
  LIB_LINK += $(NAG_LIB)
endif
ifeq ($(LNTCC),T)
  MOD1_PATH += $(NTCC_INC)
  PRECOMP += -DLADAS -DNTCC
  LIB_LINK += $(NTCC_LIB)
  MOD_PATH += $(NTCC_INC)
endif
ifeq ($(LSILO),T)
  MOD1_PATH += $(SILO_INC)
  PRECOMP += -DSILO_AVAI
  LIB_LINK += $(SILO_LIB)
endif

# This is hear because of MPI above resetting LINK
ifeq ($(MAKECMDGOALS),libstell.a)
  LINK = $(LINK_AR)
  COMPILE += -fPIC
  COMPILE_FREE += -fPIC
  INC_PATH = 
endif
ifeq ($(MAKECMDGOALS),libstell.so)
  LIB_LINK = $(LIB_SHARE)
  COMPILE += -fPIC
  COMPILE_FREE += -fPIC
  INC_PATH = 
endif

#This must come here
MOD_PATH+= -I
