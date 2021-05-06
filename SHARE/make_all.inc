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
SHARED_EXT ?= .so
STATIC_EXT ?= .a
LIB = libstell$(STATIC_EXT)

#######################################################################

# For debugging print out the make-command goals and the machine name
#$(info $$MAKECMDGOALS is [${MAKECMDGOALS}])
#$(info $$MACHINE is [${MACHINE}])
#
# These are old switches that worked but were confusing
#ifeq ($(strip $(MACHINE)),mingw64)
#ifneq (,$(findstring MINGW,$(shell uname -s)))
#
# Debugging stuff
#tstresult = $(findstring MINGW,$(shell uname -s))
#$(info $$tstresult is [${tstresult}])
#
# Can we auto-detect which type of library to build and link against?
#    Two steps:
#       machine type (UNIX/POSIX versus WIN)
#       static versus shared
#
#
# If the machine is using MINGW assume that it is WINDOWS
#ifeq (MINGW,$(findstring MINGW,$(shell uname -s)))
#  LIB := libstell.lib
#  ifeq ($(findstring libstell.dll,$(MAKECMDGOALS)),libstell.dll)
#    LIB := libstell.dll
#  endif
#else
#  LIB := libstell.a
#  ifeq ($(findstring shared,$(MAKECMDGOALS)),shared)
#    LIB := libstell.so
#  endif
#endif
##$(info $$LIB is [${LIB}])

#######################################################################

# Handle Which Code is being CALLED
ifeq ($(MAKECMDGOALS),xstelloptv2)
  PRECOMP += -DBNORM -DBEAMS3D_OPT -DJINV_OPT -DTXPORT_OPT -DDIAGNO_OPT -DSKS  -DSKS2
  VMEC_DIR = ../../VMEC2000
  LIB_VMEC= libvmec$(STATIC_EX)
  MOD_PATH+= -I$(VMEC_DIR)/$(LOCTYPE)
  BEAMS3D_DIR = ../../BEAMS3D
  LIB_BEAMS3D = libbeams3d$(STATIC_EX)
  MOD_PATH+= -I$(BEAMS3D_DIR)/$(LOCTYPE)
  BOOTSJ_DIR = ../../BOOTSJ
  LIB_BOOTSJ = libbootsj$(STATIC_EX)
  MOD_PATH+= -I$(BOOTSJ_DIR)/$(LOCTYPE)
  BNORM_DIR = ../../BNORM
  LIB_BNORM = libbnorm$(STATIC_EX)
  MOD_PATH+= -I$(BNORM_DIR)/$(LOCTYPE)
  BOOZ_DIR = ../../BOOZ_XFORM
  LIB_BOOZ= libbooz$(STATIC_EX)
  MOD_PATH+= -I$(BOOZ_DIR)/$(LOCTYPE)
  COBRA_DIR = ../../COBRAVMEC
  LIB_COBRA = libcobra$(STATIC_EX)
  MOD_PATH+= -I$(COBRA_DIR)/$(LOCTYPE)
  DIAGNO_DIR = ../../DIAGNO
  LIB_DIAGNO = libdiagno$(STATIC_EX)
  MOD_PATH+= -I$(DIAGNO_DIR)/$(LOCTYPE)
  FIELDLINES_DIR = ../../FIELDLINES
  LIB_FIELDLINES = libfieldlines$(STATIC_EX)
  MOD_PATH+= -I$(FIELDLINES_DIR)/$(LOCTYPE)
  JINV_DIR = ../../J_INVARIANT
  LIB_JINV = libjinv$(STATIC_EX)
  MOD_PATH+= -I$(JINV_DIR)/$(LOCTYPE)
  MGRID_DIR= ../../MAKEGRID
  LIB_MGRID = libmakegrid$(STATIC_EX)
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
    LIB_DKES = libdkes$(STATIC_EX)
    MOD_PATH += -I$(DKES_DIR)/$(LOCTYPE)
    PRECOMP += -DDKES_OPT
    LIB_LINK += $(DKES_DIR)/$(LOCTYPE)/$(LIB_DKES)
  else
    DKES_DIR =
    LIB_DKES =
  endif
  ifeq ($(LNEO),T)
    NEO_DIR = ../../NEO
    LIB_NEO = libneo$(STATIC_EX)
    MOD_PATH += -I$(NEO_DIR)/$(LOCTYPE)
    PRECOMP += -DNEO_OPT
    LIB_LINK += $(NEO_DIR)/$(LOCTYPE)/$(LIB_NEO)
  else
    NEO_DIR =
    LIB_NEO =
  endif
else ifeq ($(MAKECMDGOALS),help)
else ifeq ($(MAKECMDGOALS),test_stellopt)
else ifeq ($(MAKECMDGOALS),test_stellopt_coilopt)
else ifeq ($(MAKECMDGOALS),test_stellopt_gene)
else ifeq ($(MAKECMDGOALS),test_stellopt_regcoil)
else ifeq ($(MAKECMDGOALS),test_stellopt_sfincs)
else ifeq ($(MAKECMDGOALS),test_stellopt_terpsichore)
else ifeq ($(MAKECMDGOALS),test_stellopt_travis)
else
  LIB_LINK =
  LGENE = F
  LCOILOPT = F
  LSFINCS = F
  LREGCOIL= F
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

# Handle other codes
LIB_LINK += $(COILOPT_LIB) $(GENE_LIB) $(REGCOIL_LIB) $(SFINCS_LIB) $(TERPSICHORE_LIB) $(TRAVIS_LIB)

# Add libraries to LIB_LINK
# LIB_LINK+= $(MYHOME)/libstell$(STATIC_EX) $(LIBS)
# This version avoids referencing the bin directory for the symbolic link to lib.
# Which we can't reference since we now have the awk file in the main directory.
LIB_LINK+= $(LIB_DIR)/$(LOCTYPE)/$(LIB) $(LIBS)

# 3rd Party Libs
ifeq ($(LMPI),T)
  PRECOMP += -DMPI_OPT
  COMPILE = $(MPI_COMPILE)
  COMPILE_FREE = $(MPI_COMPILE_FREE)
  LIB_LINK += $(MPI_LIB)
  LINK    = $(MPI_COMPILE) $(FLAGS) -o
  LINK_C  = $(MPI_LINK)
  LIB_SHARE += $(MPI_LIB)
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
  PRECOMP += -DLHDF5 -DH5_USE_16_API
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

# This is hear because of MPI above resetting LINK
#ifeq ($(MAKECMDGOALS), $(filter $(MAKECMDGOALS),libstell.a libstell.lib))
ifeq ($(MAKECMDGOALS), libstell$(STATIC_EXT))
  LINK = $(LINK_AR)
  COMPILE += -fPIC
  COMPILE_FREE += -fPIC
  INC_PATH =
endif
#ifeq ($(MAKECMDGOALS), $(filter $(MAKECMDGOALS),libstell.so libstell.dll))
ifeq ($(MAKECMDGOALS), libstell$(SHARED_EXT))
  LINK = $(LINK_C)
  LIB_LINK = $(LIB_SHARE)
  COMPILE += -fPIC
  COMPILE_FREE += -fPIC
  INC_PATH =
endif

#######################################################################
#     ld has a cool feature to deal with recursive dependencies
#
#     not included in the default ld linker on osx with xcode
#
#######################################################################

ifeq ($(shell uname -s),Darwin)
  LIB_LINK := $(LIB_LINK) $(LIB_LINK)
else
  LIB_LINK := -Wl,--start-group $(LIB_LINK) -Wl,--end-group
endif

#######################################################################
#   If compiling a DLL we should also make an import library
#
#######################################################################

ifeq ($(MAKECMDGOALS),libstell.dll)
  LIB_LINK += -Wl,--output-def,libstell.def,--out-implib,libstell.dll.a
endif

#######################################################################

#This must come here
MOD_PATH+= -I

#######################################################################
#######################################################################
