#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
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
#  FLAGS_R = -O2 -g -fexternal-blas -fbacktrace -fcheck=all 
  FLAGS_R = -g -O2 -fexternal-blas -fbacktrace
  FLAGS_D = -g -O0 -fexternal-blas -fcheck=all 
  LIBS    = -L$(BLASHOME)/lib -lopenblas \
            -L$(SCALAPACK_HOME)/lib -lscalapack \
            -L$(BLACS_HOME)/lib -lmpiblacs -lmpiblacsCinit \
            -lmpiblacsF77init -lmpiblacs 

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = $(MPIHOME)/bin/mpif90
  MPI_COMPILE_FREE = $(MPIHOME)/bin/mpif90 -ffree-form \
                     -ffree-line-length-none -ffixed-line-length-none
  MPI_COMPILE_C = $(MPIHOME)/bin/mpicc
  MPI_LINK = $(MPIHOME)/bin/mpif90

  # For GCC6
  #MPI_COMPILE = mpif90
  #MPI_COMPILE_FREE = mpif90 -ffree-form \
  #                   -ffree-line-length-none -ffixed-line-length-none
  #MPI_COMPILE_C = mpicc
  #MPI_LINK = mpif90

#######################################################################
#            NAG Options
#######################################################################
  LNAG = T
  NAG_LIB = -L$(NAG_ROOT)/lib -lnag_nag

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = -I$(NETCDF_HOME)/include
  NETCDF_LIB = -L$(NETCDF_HOME)/lib -lnetcdf -lnetcdff

  # FOR GCC6
  #NETCDF_LIB = -L$(NETCDF_HOME)/lib -lnetcdff

#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = T
  NTCC_INC = -I$(NTCCHOME)/mod
  NTCC_LIB = -L$(NTCCHOME)/lib -laladdinsub -lr8slatec -ladpak\
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
  LPGPLOT = T
  PGPLOT_INC = -I$(PGPLOT_DIR)
  PGPLOT_LIB = -L$(PGPLOT_DIR) -lpgplot -L/usr/X11R6/lib -lX11

#######################################################################
#             SILO Options
#######################################################################
  LSILO = T
  SILO_INC = -I$(SILOHOME)/include
  SILO_LIB = -L$(SILOHOME)/lib -lsilo

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
#            LIBSTELL Shared Options
#######################################################################
LIB_SHARE = $(BLASHOME)/lib/libopenblas.so \
          $(SCALAPACK_HOME)/lib/libscalapack.a \
          $(BLACS_HOME)/lib/libmpiblacs.a  $(BLACS_HOME)/lib/libmpiblacsCinit.a $(BLACS_HOME)/lib/libmpiblacsF77init.a \
          $(NTCCHOME)/lib/libaladdinsub.a $(NTCCHOME)/lib/libr8slatec.a $(NTCCHOME)/lib/libadpak.a $(NTCCHOME)/lib/libcppsub.a $(NTCCHOME)/lib/libcomput.a \
          $(NTCCHOME)/lib/libpspline.a $(NTCCHOME)/lib/libportlib.a $(NTCCHOME)/lib/libezcdf.a $(NTCCHOME)/lib/libmds_sub.a $(NTCCHOME)/lib/libmdstransp.a $(NTCCHOME)/lib/libvaxonly.a \
          $(HDF5_HOME)/lib/libhdf5_hl.so $(HDF5_HOME)/lib/libhdf5_fortran.so $(HDF5_HOME)/lib/libhdf5hl_fortran.so $(HDF5_HOME)/lib/libhdf5.so \
          $(NAG_ROOT)/lib/libnag_nag.so  \
          $(NETCDFHOME)/lib/libnetcdf.so $(NETCDFHOME)/lib/libnetcdff.so $(NETCDFHOME)/lib/libnetcdf_c++.so \
          $(SILOHOME)/lib/libsilo.a \
          $(GSLHOME)/lib/libgsl.so \
          $(GCC4_HOME)/lib64/libgfortran.so $(GCC4_HOME)/lib64/libstdc++.so \
          $(MPIHOME)/lib/libmpi.so $(MPIHOME)/lib/libmpi_f77.so $(MPIHOME)/lib/libmpi_f90.so 
#          /usr/lib/libc.so /usr/lib/libm.so /usr/lib/libpthread.so 

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
  PRECOMP += -DBNORM -DDKES_OPT -DBEAMS3D_OPT -DJINV_OPT -DNEO_OPT -DTXPORT_OPT -DDIAGNO_OPT  -DDKES_OPT  -DNEO_OPT  -DSKS  -DSKS2
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
  DKES_DIR= ../../DKES
  LIB_DKES= libdkes.a
  MOD_PATH+= -I$(DKES_DIR)/$(LOCTYPE)
  JINV_DIR = ../../J_INVARIANT
  LIB_JINV = libjinv.a
  MOD_PATH+= -I$(JINV_DIR)/$(LOCTYPE)
  MGRID_DIR= ../../MAKEGRID
  LIB_MGRID= libmakegrid.a
  MOD_PATH+= -I$(MGRID_DIR)/$(LOCTYPE)
  NEO_DIR = ../../NEO
  LIB_NEO = libneo.a
  MOD_PATH+= -I$(NEO_DIR)/$(LOCTYPE)
  LIB_LINK=  $(VMEC_DIR)/$(LOCTYPE)/$(LIB_VMEC)  \
	$(BEAMS3D_DIR)/$(LOCTYPE)/$(LIB_BEAMS3D) \
	$(BOOTSJ_DIR)/$(LOCTYPE)/$(LIB_BOOTSJ) \
	$(BNORM_DIR)/$(LOCTYPE)/$(LIB_BNORM) \
	$(BOOZ_DIR)/$(LOCTYPE)/$(LIB_BOOZ) \
	$(COBRA_DIR)/$(LOCTYPE)/$(LIB_COBRA) \
	$(DIAGNO_DIR)/$(LOCTYPE)/$(LIB_DIAGNO) \
	$(DKES_DIR)/$(LOCTYPE)/$(LIB_DKES) \
	$(JINV_DIR)/$(LOCTYPE)/$(LIB_JINV) \
	$(MGRID_DIR)/$(LOCTYPE)/$(LIB_MGRID) \
	$(NEO_DIR)/$(LOCTYPE)/$(LIB_NEO)
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

# Handle other codes
LIB_LINK += $(GENE_LIB) $(COILOPT_LIB) $(TERPSICHORE_LIB) $(TRAVIS_LIB)

# Add libraries to LIB_LINK
# LIB_LINK+= $(MYHOME)/libstell.a $(LIBS)
# This version avoids referencing the bin directory for the symbolic link to lib.
# Which we can't reference since we now have the awk file in the main directory.
LIB_LINK+= $(LIB_DIR)/$(LOCTYPE)/$(LIB) $(LIBS)

# 3rd Party Libs
ifeq ($(LMPI),T)
  PRECOMP += -DMPI_OPT
  VMEC_DIR = ../../PARVMEC
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
