#######################################################################
#   This is universal to all build systems, do not edit unless your
#   adding a NEW code or library to the build.
#   New means NEW not new version or new system.
#######################################################################

# Home bin specifies where the build scripts are located,
# for example awk_cdir.awk.
HOME_BIN = ../../
MOD1_PATH = -I$(MYHOME)/libstell_dir
MOD_PATH = 
CLEAN_ARCH = 
LIB_DIR = ../../LIBSTELL
LIB = libstell.a

# Handle third party code dependencies
ifeq ($(LSFINCS),T)
  include ${PETSC_DIR}/lib/petsc/conf/variables
  SFINCS_INC += -I${PETSC_DIR}/include
  SFINCS_LIB += $(PETSC_LIB) $(FC_LINKER_FLAGS) $(FFLAGS)
endif
ifeq ($(LGENE),T)
  GENE_LIB +=-L$(FFTW_HOME)/lib -lfftw3 -lfftw3f -Xlinker \
             -L$(SLEPC_HOME)/lib -lslepc \
             -L$(PETSC_HOME)/lib -lpetsc
endif

LIMAS?=F

# Handle Which Code is being CALLED
ifeq ($(MAKECMDGOALS),xstelloptv2)
  PRECOMP += -DBNORM -DBEAMS3D_OPT -DJINV_OPT -DTXPORT_OPT -DDIAGNO_OPT -DSKS  -DSKS2
  VMEC_DIR = ../../VMEC2000
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
  FIELDLINES_DIR = ../../FIELDLINES
  LIB_FIELDLINES = libfieldlines.a
  MOD_PATH+= -I$(FIELDLINES_DIR)/$(LOCTYPE)
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
	$(FIELDLINES_DIR)/$(LOCTYPE)/$(LIB_FIELDLINES) \
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
else ifeq ($(MAKECMDGOALS),xthrift)
  PRECOMP += -DSKS  -DSKS2
  VMEC_DIR = ../../VMEC2000
  LIB_VMEC= libvmec.a
  MOD_PATH+= -I$(VMEC_DIR)/$(LOCTYPE)
  BOOTSJ_DIR = ../../BOOTSJ
  LIB_BOOTSJ = libbootsj.a
  MOD_PATH+= -I$(BOOTSJ_DIR)/$(LOCTYPE)
  BOOZ_DIR = ../../BOOZ_XFORM
  LIB_BOOZ= libbooz.a
  MOD_PATH+= -I$(BOOZ_DIR)/$(LOCTYPE)
  DIAGNO_DIR = ../../DIAGNO
  LIB_DIAGNO= libdiagno.a
  MOD_PATH+= -I$(DIAGNO_DIR)/$(LOCTYPE)
  LIB_LINK=  $(VMEC_DIR)/$(LOCTYPE)/$(LIB_VMEC)  \
  $(BOOTSJ_DIR)/$(LOCTYPE)/$(LIB_BOOTSJ) \
  $(BOOZ_DIR)/$(LOCTYPE)/$(LIB_BOOZ) \
  $(DIAGNO_DIR)/$(LOCTYPE)/$(LIB_DIAGNO)
  LGENE = F
  LCOILOPT = F
  LSFINCS = F
  LREGCOIL= F
  LMANGO = F
  LTERPSICHORE = F
else ifeq ($(MAKECMDGOALS),help)
else ifeq ($(MAKECMDGOALS),test_stellopt)
else ifeq ($(MAKECMDGOALS),test_stellopt_coilopt)
else ifeq ($(MAKECMDGOALS),test_stellopt_gene)
else ifeq ($(MAKECMDGOALS),test_stellopt_regcoil)
else ifeq ($(MAKECMDGOALS),test_stellopt_sfincs)
else ifeq ($(MAKECMDGOALS),test_stellopt_terpsichore)
else ifeq ($(MAKECMDGOALS),test_stellopt_travis)
else
  LGENE = F
  LCOILOPT = F
  LSFINCS = F
  LREGCOIL= F
  LMANGO = F
  LTERPSICHORE = F
  LTRAVIS = F
endif
ifeq ($(MAKECMDGOALS),xanimec)
  PRECOMP += -D_ANIMEC
endif
ifeq ($(MAKECMDGOALS),xvmec2000)
  PRECOMP += -DSKS  -DSKS2
endif
ifeq ($(MAKECMDGOALS),xsiesta)
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
ifeq ($(LCOILOPT),T)
  MOD1_PATH += $(COILOPT_INC)
  PRECOMP += -DCOILOPTPP
  MOD_PATH += -I$(COILOPTPP_DIR)
else
  COILOPT_LIB = 
  COILOPTPP_DIR = 
  LIB_COILOPTPP = 
endif
ifeq ($(LGENE),T)
  MOD1_PATH += $(GENE_INC)
  PRECOMP += -DGENE
  MOD_PATH += -I$(GENE_DIR)
else
  GENE_LIB = 
  GENE_DIR = 
  LIB_GENE = 
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
ifeq ($(LSFINCS),T)
  MOD1_PATH += $(SFINCS_INC)
  PRECOMP += -DSFINCS
  MOD_PATH += $(SFINCS_INC)
else
  SFINCS_LIB =
  SFINCS_DIR =
  LIB_SFINCS =
endif
ifeq ($(LMANGO),T)
  # This next line will define the variables MANGO_C_COMPILE_FLAGS, MANGO_F_COMPILE_FLAGS, MANGO_C_LINK_FLAGS, and MANGO_F_LINK_FLAGS
  include $(MANGO_DIR)/lib/mangoMakeVariables
  # Add any flags defined by the individual stellopt make_*.inc files:
  MANGO_INC += $(MANGO_F_COMPILE_FLAGS)
  MANGO_LIB += $(MANGO_F_LINK_FLAGS)
  MOD1_PATH += $(MANGO_INC)
  PRECOMP += -DMANGO
  MOD_PATH += $(MANGO_INC)
else
  MANGO_LIB =
  MANGO_DIR =
  LIB_MANGO =
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
ifeq ($(LAEOPT),T)
  MOD1_PATH += $(AEOPT_INC)
  PRECOMP += -DAEOPT
  MOD_PATH += -I$(AEOPT_DIR)
else
  AEOPT_LIB = 
  AEOPT_DIR = 
  LIB_AEOPT = 
endif

# Handle other codes
LIB_LINK += $(COILOPT_LIB) $(GENE_LIB) $(REGCOIL_LIB) $(SFINCS_LIB) $(MANGO_LIB) $(TERPSICHORE_LIB) $(TRAVIS_LIB) $(AEOPT_LIB)

# Add libraries to LIB_LINK
# LIB_LINK+= $(MYHOME)/libstell.a $(LIBS)
# This version avoids referencing the bin directory for the symbolic link to lib.
# Which we can't reference since we now have the awk file in the main directory.
LIB_LINK+= $(LIB_DIR)/$(LOCTYPE)/$(LIB) $(LIBS)

# 3rd Party Libs
ifeq ($(LMPI),T)
  PRECOMP += -DMPI_OPT
  COMPILE = $(MPI_COMPILE)
  COMPILE_FREE = $(MPI_COMPILE_FREE)
  LINK    = $(MPI_LINK) $(FLAGS) -o 
  LINK_C  = $(MPI_LINK) -shared
endif
ifeq ($(LNETCDF),T)
  MOD1_PATH += $(NETCDF_INC)
  PRECOMP += -DNETCDF
  LIB_LINK += $(NETCDF_LIB)
  MOD_PATH += $(NETCDF_INC)
  LIB_SHARE += $(NETCDF_LIB)
endif
ifeq ($(LFFTW3),T)
  MOD1_PATH += $(FFTW3_INC)
  PRECOMP += -DFFTW3
  LIB_LINK += $(FFTW3_LIB)
  MOD_PATH += $(FFTW3_INC)
  LIB_SHARE += $(FFTW3_LIB)
endif
ifeq ($(LHDF5),T)
  MOD1_PATH += $(HDF5_INC)
  PRECOMP += -DLHDF5
  LIB_LINK += $(HDF5_LIB)
  MOD_PATH += $(HDF5_INC)
  LIB_SHARE += $(HDF5_LIB)
endif
ifeq ($(LNAG),T)
  PRECOMP += -DNAG
  LIB_LINK += $(NAG_LIB)
  LIB_SHARE += $(NAG_LIB)
endif
ifeq ($(LSILO),T)
  MOD1_PATH += $(SILO_INC)
  PRECOMP += -DSILO_AVAI
  LIB_LINK += $(SILO_LIB)
  LIB_SHARE += $(SILO_LIB)
endif
ifeq ($(LIMAS),T)
  MOD1_PATH += $(IMAS_INC)
  PRECOMP += -DIMAS
  LIB_LINK += $(IMAS_LIB)
  LIB_SHARE += $(IMAS_LIB)
endif

# This is hear because of MPI above resetting LINK
ifeq ($(MAKECMDGOALS),libstell.a)
  LINK = $(LINK_AR)
  COMPILE += -fPIC
  COMPILE_FREE += -fPIC
  INC_PATH = 
endif
ifeq ($(MAKECMDGOALS),libstell.so)
  LINK = $(LINK_C) -shared
  LIB_LINK = $(LIB_SHARE)
  COMPILE += -fPIC
  COMPILE_FREE += -fPIC
  INC_PATH = 
endif

# Pass Git information to code
#GIT_BRANCH  = $(shell git branch --show-current) #Only works on version >2.2
GIT_REPO    = $(shell git ls-remote --get-url)
GIT_BRANCH  = $(shell git rev-parse --abbrev-ref HEAD)
GIT_VERSION = $(shell git describe --abbrev=4 --dirty --always --tags)
GIT_HASH    = $(shell git show -s --format=%H)
TIME        = $(shell date +"%d.%m.%Y %H:%M:%S")
PRECOMP += -DGIT_VERSION_EXT="'$(GIT_VERSION)'"
PRECOMP += -DGIT_HASH_EXT="'$(GIT_HASH)'"
PRECOMP += -DGIT_BRANCH_EXT="'$(GIT_BRANCH)'"
PRECOMP += -DGIT_REPO_EXT="'$(GIT_REPO)'"
PRECOMP += -DBUILT_ON_EXT="'$(TIME)'"

#This must come here
MOD_PATH+= -I
