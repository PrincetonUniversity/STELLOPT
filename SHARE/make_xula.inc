#######################################################################
#            Define Basic Utilities
#######################################################################

  SHELL = /bin/bash
  PWD1 = `pwd`
  #compiladores=intel   # cambiar para compilar con gnu
  compiladores=intel
  #MYHOME = $(HOME)/STELLOPT/CIEMAT
  MYHOME = $(HOME)/STELLOPT/CIEMAT/xintel
  
#  MYHOME = $(HOME)/STELLOPT/CIEMAT/x$(compiladores)
#  MYHOME = $(HOME)/STELLOPT/CIEMAT/xTemp

  PRECOMP:= /lib/cpp -P -traditional -DLINUX
  COMPILE = $(FC)
  ifeq ($(compiladores),intel)
  	COMPILE_FREE = $(FC) -stand f95 -ip  
  	OPENBLAS_HOME = /mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/openblas-0.3.6-bwopn5fpda3kcxclpyto2sqybautaofg/
  endif 
  ifeq ($(compiladores),gnu)
  	COMPILE_FREE = $(FC) -stand f95 -ip  -ffree-form -ffree-line-length-none -ffixed-line-length-none 
 	OPENBLAS_HOME = /mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/gcc-8.2.0/openblas-0.3.3-5qt5gvtw2ctxm7m6jzmmepr6lnxfjbxz/
  endif
  mpiExe=`which mpif90`
  mpiBin=`dirname $mpiExe`
  MPI_HOME=`dirname $mpiBin`
#
  LINK    = ld -lc $(FLAGS) -o
  LINK_AR = ar -ruv
#  LINK_C  = $(CC) -shared -Wl,-z-defs
  LINK_C  = $(CC) -shared -Wl,--no-undefined -Wl,-export-dynamic 



#######################################################################
#            Define Compiler Flags
# Crutch for broken module vars
  NETCDF_F = $(NETCDF_HOME)
  NETCDF_C = $(NETCDF_HOME)


#######################################################################
#            Define Compiler Flags
#######################################################################
#  FLAGS_R = -O2 -g -fexternal-blas -fbacktrace -fcheck=all 
  #FLAGS_R = -g -O2 -fexternal-blas -fbacktrace -I/usr/include 
  FLAGS_R = -g -O2   -I/usr/include \
            -I$(OPENBLAS_HOME)/include/ -I$(MPI_HOME)/include 

  #FLAGS_D = -g -O0 -fexternal-blas -fcheck=all -I/usr/include
  FLAGS_D = -g -O0  -fcheck=all -I/usr/include -Wextra \
   -Wtarget-lifetime -fbounds-check -ffpe-trap=zero -finit-real=snan\
             -I$(OPENBLAS_HOME)/include/ \
        	-I$(MPI_HOME)/include 

  LIBS    = -L$(OPENBLAS_HOME)/lib \
		-lopenblas -lm -L$(SCALAPACK_HOME)/lib -lscalapack -lc \
             -L$(MPI_HOME)/lib 

  ifeq ($(compiladores),gnu)
	LIBS +=  -lgfortran 
  endif

#-L/lib64/openmpi/lib -lmpiblacs /lib64/libopenblas.so.0 \
 #           /lib64/openmpi/lib/libscalapack.so.2

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpif90
  ifeq ($(compiladores),gnu)
  	MPI_COMPILE_FREE = mpif90 -ffree-form \
                     -ffree-line-length-none -ffixed-line-length-none
  endif
  ifeq ($(compiladores),intel)
	MPI_COMPILE_FREE = mpif90
  endif
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
  NAG_LIB = 

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = -I$(NETCDF_HOME)/include
  NETCDF_LIB = -L$(MPI_HOME)/lib -L$(NETCDF_HOME)/lib -lnetcdff

#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = F
  NTCC_INC = -I$(NTCCHOME)/mod
  NTCC_LIB = -L$(NTCCHOME)/lib -laladdinsub -lr8slatec -ladpak\
             -lcppsub -lcomput -lpspline -lportlib -lezcdf -lmds_sub \
             -lmdstransp -lvaxonly

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  HDF5_INC = -I$(HDF5_HOME)/include -DH5_USE_16_API
  HDF5_LIB = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -ldl -lm -Wl,-rpath,-z -L$(PGPLOT_DIR) -lpgplot -L$(PNG_DIR)/lib -I$(PNG_DIR)/include -lpng
# -Wl,-rpath -Wl,-z,relro
#/usr/lib64/openmpi/lib \
#              /usr/lib64/openmpi/lib/libhdf5hl_fortran.so \
#              /usr/lib64/openmpi/lib/libhdf5_hl.so \
#              /usr/lib64/openmpi/lib/libhdf5_fortran.so \
#              /usr/lib64/openmpi/lib/libhdf5.so \
#              -Wl,-z,relro -lhdf5_fortran -lhdf5 -lz -lsz -lz -ldl -lm -Wl,-rpath -Wl,/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/gcc-8.2.0/openmpi-3.1.3-roaqnl35m7f57g7tluunnp6bdxomkfrg/lib

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT=F
  PGPLOT_INC = -I$(PGPLOT_DIR)/ -I$(PNG_DIR)/include
  PGPLOT_LIB = -L$(PGPLOT_DIR)/ -lpgplot -L$(PNG_DIR)/lib -lpng
#######################################################################
#             SILO Options
#######################################################################
  LSILO = T
  SILO_INC = 
  SILO_LIB = 

#######################################################################
#            DKES/NEO Options
#######################################################################
  LDKES = T
  LNEO  = T

#######################################################################
#            KNOSOS Options
#######################################################################
  LKNOSOS = T
  KNOSOS_PATH = ../../KNOSOS
  KNOSOS_INC = -I$(KNOSOS_PATH)
  KNOSOS_DIR = $(KNOSOS_PATH)
  LIB_KNOSOS = libknosos.a
  KNOSOS_LIB = -lfftw3 

#######################################################################
#            STELLA Options
#######################################################################
  LSTELLA = T
  STELLA_PATH= ../../STELLA
  STELLA_INC = -I$(STELLA_PATH)
  STELLA_DIR = $(STELLA_PATH)
  LIB_STELLA = libstella.a
#  STELLA_LIB = -lfftw3 #$(STELLA_DIR)/$(LIB_STELLA) \
#       -L${PETSC_DIR}/lib -lpetsc \
#       -L/${MKL_HOME}/lib/intel64 -mkl \
#       -lfftw3 -Wl,-rpath,${MKL_HOME}/lib/intel64,-rpath,${PETSC_DIR}/lib

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
  COILOPT_PATH= $(HOME)/STELLOPT/CIEMAT/COILOPT
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
  REGCOIL_DIR = $(HOME)/STELLOPT/CIEMAT/regcoil-master
#$(REGCOIL_PATH)
  REGCOIL_INC = -I$(REGCOIL_DIR) 
  LIB_REGCOIL = libregcoil.a
  REGCOIL_LIB = $(REGCOIL_DIR)/$(LIB_REGCOIL) -fopenmp

#######################################################################
#            LIBSTELL Shared Options
######################################################################
LIB_SHARE = -lz -lc -lm -lpthread $(LIBS) -lc
  
#-lc -lgfortran -lstdc++ -lmpi -lmpi_mpifh -lz -lc -lm -lpthread $(LIBS) -lc


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
  ifeq ($(LKNOSOS),T)
    KNOSOS_DIR = ../../KNOSOS
    LIB_KNOSOS = libknosos.a
    MOD_PATH+= -I$(KNOSOS_DIR)/$(LOCTYPE)
    PRECOMP += -DKNOSOS_OPT
    LIB_LINK += $(KNOSOS_DIR)/$(LOCTYPE)/$(LIB_KNOSOS)
  else
    KNOSOS_DIR =
    LIB_KNOSOS =
  endif
  ifeq ($(LSTELLA),T)
    STELLA_DIR = ../../STELLA
    LIB_STELLA = libstella.a
    MOD_PATH += -I$(STELLA_DIR)/$(LOCTYPE)
    PRECOMP += -DSTELLA_OPT
    LIB_LINK += $(STELLA_DIR)/$(LOCTYPE)/$(LIB_STELLA)
  else
    STELLA_DIR =
    LIB_STELLA =
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
  LGENE = F
  LCOILOPT = F
  LTERPSICHORE = F
  LREGCOIL= F
  LTRAVIS = F
  LSFINCS = F
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
ifeq ($(LSFINCS),T)
  MOD1_PATH += $(SFINCS_INC)
  PRECOMP += -DSFINCS
  MOD_PATH += -I$(SFINCS_DIR)
else
  SFINCS_LIB = 
  SFINCS_DIR = 
  LIB_SFINCS = 
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
LIB_LINK += $(KNOSOS_LIB) $(STELLA_LIB) $(GENE_LIB) $(COILOPT_LIB) $(TERPSICHORE_LIB) $(TRAVIS_LIB) $(REGCOIL_LIB)

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
  LINK    = $(MPI_COMPILE) $(FLAGS) -o 
  LINK_C  = $(MPI_LINK) 
endif
ifeq ($(LNETCDF),T)
  MOD1_PATH += $(NETCDF_INC)
  PRECOMP += -DNETCDF
  LIB_LINK += $(NETCDF_LIB)
  MOD_PATH += $(NETCDF_INC)
  LIB_SHARE += $(NETCDF_LIB)
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
ifeq ($(LNTCC),T)
  MOD1_PATH += $(NTCC_INC)
  PRECOMP += -DLADAS -DNTCC
  LIB_LINK += $(NTCC_LIB)
  MOD_PATH += $(NTCC_INC)
  LIB_SHARE += $(NTCC_LIB)
endif
ifeq ($(LSILO),T)
  MOD1_PATH += $(SILO_INC)
  PRECOMP += -DSILO_AVAI
  LIB_LINK += $(SILO_LIB)
  LIB_SHARE += $(SILO_LIB)
endif

# This is hear because of MPI above resetting LINK
ifeq ($(MAKECMDGOALS),libstell.a)
  LINK = $(LINK_AR)
  COMPILE += -fPIC
  COMPILE_FREE += -fPIC
  INC_PATH = 
endif
ifeq ($(MAKECMDGOALS),libstell.so)
  LINK = $(LINK_C)
  LIB_LINK = $(LIB_SHARE)
  COMPILE += -fPIC
  COMPILE_FREE += -fPIC
  INC_PATH = 
endif

#This must come here
MOD_PATH+= -I
