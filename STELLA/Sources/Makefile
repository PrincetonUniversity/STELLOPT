####################################################################
#
#  Makefile for stella gyrokinetic turbulence code 
#
#  (requires GNU's gmake)
#
GK_PROJECT ?= stella
#
#  Makefile written by Bill Dorland and Ryusuke Numata
#
#  LAST UPDATE: 10/12/09
#
# * Available Compilers (tested on limited hosts)
#   (must be Fortran 95 Standard compliant)
#
# Intel ifort
# GNU's gfortran and g95
# IBM XL Fortran xlf90
# PathScale Compiler Suite pathf90
# The Portland Group pgf90
# NAGWare f95 (v5.1)
# Lahey/Fujitsu Fortran lf95
# 
# * Frequently Tested Hosts, Systems
#
# Standard Linux
# Standard Mac OS X with MacPorts
# Franklin at NERSC and Jaguar at NCCS (Cray XT4 with PGI)
# Bassi at NERSC (IBM Power 5 with IBM XL Fortran)
# Ranger (Sun Constellation Linux Cluster with Intel)
#
# * Switches:
#
# Here is a list of switches with simple explanation.
# In the brackets, values accepted are shown,
# where "undefined" means blank.
# Switches with (bin) are checked if they are defined or not
# What values they have do not matter.
# Be careful that DEBUG=off means DEBUG=on.
#
# turns on debug mode (bin)
DEBUG ?=
# turns on test mode (bin)
TEST ?=
# optimization (on,aggressive,undefined)
OPT ?= on
# prevents linking with shared libraries (bin)
STATIC ?=
# promotes precisions of real and complex (bin)
DBLE ?= on
# turns on distributed memory parallelization using MPI (bin)
USE_MPI ?= on
# which FFT library to use (fftw,fftw3,mkl_fftw,undefined) 
USE_FFT ?= fftw3
# uses netcdf library (bin)
USE_NETCDF ?= on
# uses parallel netcdf library
USE_PARALLEL_NETCDF ?=
# uses hdf5 library (bin)
USE_HDF5 ?= 
# Use local random number generator (mt,undefined)
# see also README
USE_LOCAL_RAN ?=
# Use local special functions (bin)
USE_LOCAL_SPFUNC ?=
# Use nag libraray (spfunc,undefined)
USE_NAGLIB ?=
# link to sfincs library at compilation
USE_SFINCS ?=
# Use LAPACK, needed for test particle collisions
USE_LAPACK ?= on
# does the compiler support the iso_c_binding features of Fortran? 
# (needed for local parallel LU decomposition) 
HAS_ISO_C_BINDING ?= on
#
# * Targets:
#
#  depend: generate dependency
#  test_make: print values of variables
#  clean: clean up
#  distclean: does "make clean" + removes platform links & executables
#  tar: pack
#
############################################################### DEFAULT VALUE
#
# These variables can be set in platform-dependent Makefile.
#

MAKE		= make
CPP		= cpp
#CPPFLAGS	= -C -P -traditional
CPPFLAGS	= -P -traditional
export FC	= f90
#MPIFC		?= mpif90-mpich-gcc48
export MPIFC	?= mpif90
#export MPIFC	?= mpifort
H5FC		?= h5fc
H5FC_par	?= h5pfc
F90FLAGS	= 
F90OPTFLAGS	=
CC		= cc
#MPICC		?= mpicc-mpich-gcc48
MPICC		?= mpicc
H5CC		?= h5cc
H5CC_par	?= h5pcc
CFLAGS		=
COPTFLAGS 	=
LD 		= $(FC)
LDFLAGS 	= $(F90FLAGS)
ARCH 		= ar
ARCHFLAGS 	= cr
RANLIB		= ranlib
AWK 		= awk
PERL		= perl
FORD       ?= ford

MPI_INC	?=
MPI_LIB ?=
FFT_INC ?=
FFT_LIB ?=
NETCDF_INC ?=
NETCDF_LIB ?=
HDF5_INC ?=
HDF5_LIB ?=
NAG_LIB ?=
NAG_PREC ?= dble
SFINCS_LIB ?=
SFINCS_INC ?=
PETSC_LIB ?=
PETSC_INC ?=
LIBSTELL_LIB ?=

# Record the top level path. Note we don't just use $(PWD) as this
# resolves to the directory from which make was invoked. The approach
# taken here ensures that GK_HEAD_DIR is the location of this
# Makefile. Note the realpath call removes the trailing slash so
# later we need to add a slash if we want to address subdirectories
GK_THIS_MAKEFILE := $(abspath $(lastword $(MAKEFILE_LIST)))
GK_HEAD_DIR := $(realpath $(dir $(GK_THIS_MAKEFILE)))

######################################################### PLATFORM DEPENDENCE

# compile mode switches (DEBUG, TEST, OPT, STATIC, DBLE)
# must be set before loading Makefile.$(GK_SYSTEM) because they may affect
# compiler options.
# However, Makefile.local may override some options set in Makefile.$(GK_SYSTEM),
# thus it is included before and after Makefile.$(GK_SYSTEM)
sinclude Makefile.local

# include system-dependent make variables
ifndef GK_SYSTEM
	ifdef SYSTEM
$(warning SYSTEM environment variable is obsolete)
$(warning use GK_SYSTEM instead)
	GK_SYSTEM = $(SYSTEM)
	else
$(error GK_SYSTEM is not set)
	endif
endif
include Makefile.$(GK_SYSTEM)

sinclude Makefile.local

#############################################################################

export F90FLAGS
export NETCDF_INC
export NETCDF_LIB
export UTILS=$(GK_HEAD_DIR)/utils
export GEO=$(GK_HEAD_DIR)/geo
export VMEC=$(GEO)/vmec_interface
export LIBSTELL=$(VMEC)/mini_libstell
LIBSTELL_LIB=$(LIBSTELL)/mini_libstell.a

ifeq ($(MAKECMDGOALS),depend)
# must invoke full functionality when make depend
	MAKE += USE_HDF5=on USE_FFT=fftw3 USE_NETCDF=on USE_MPI=on \
		USE_LOCAL_BESSEL=on USE_LOCAL_RAN=mt
endif

ifdef USE_HDF5
	ifndef USE_MPI
$(error Currently, USE_HDF5 works with USE_MPI)
	endif
endif

ifndef USE_FFT
$(warning USE_FFT is off)
$(warning Be sure that nonlinear run makes no sense)
endif

ifdef HAS_ISO_C_BINDING
	CPPFLAGS += -DISO_C_BINDING
endif

ifdef USE_MPI
	FC = $(MPIFC)
	CC = $(MPICC)
	CPPFLAGS += -DMPI
endif
ifeq ($(USE_FFT),fftw)
	CPPFLAGS += -DFFT=_FFTW_
	ifeq ($(FFT_LIB),)
		FFT_LIB = -lfftw -lrfftw
	endif
endif

ifeq ($(USE_FFT),fftw3)
	CPPFLAGS += -DFFT=_FFTW3_
	ifeq ($(FFT_LIB),)
		FFT_LIB = -lfftw3
	endif
endif

ifeq ($(USE_FFT),mkl_fftw)
	CPPFLAGS += -DFFT=_FFTW_
endif

ifdef USE_NETCDF
	ifeq ($(NETCDF_LIB),)
		NETCDF_LIB = -lnetcdf
        endif
	CPPFLAGS += -DNETCDF
endif
ifdef USE_LAPACK
        LAPACK_LIB ?= -llapack
	CPPFLAGS += -DLAPACK
endif
ifdef USE_HDF5
	ifdef USE_MPI
		FC = $(H5FC_par)
		CC = $(H5CC_par)
		ifdef USE_PARALLEL_NETCDF
			CPPFLAGS += -DNETCDF_PARALLEL
		endif

	else
		FC = $(H5FC)
		CC = $(H5CC)
	endif
	CPPFLAGS += -DHDF
endif
ifeq ($(USE_LOCAL_RAN),mt)
	CPPFLAGS += -DRANDOM=_RANMT_
endif
ifdef USE_LOCAL_SPFUNC
	CPPFLAGS += -DSPFUNC=_SPLOCAL_
else
	ifeq ($(findstring spfunc,$(USE_NAGLIB)),spfunc)
		CPPFLAGS += -DSPFUNC=_SPNAG_
	endif
endif
ifdef USE_NAGLIB
	ifeq ($(NAG_PREC),dble)
		ifndef DBLE
$(warning Precision mismatch with NAG libarray)	
		endif
		CPPFLAGS += -DNAG_PREC=_NAGDBLE_
	endif
	ifeq ($(NAG_PREC),sngl)
		ifdef DBLE
$(warning Precision mismatch with NAG libarray)	
		endif
		CPPFLAGS += -DNAG_PREC=_NAGSNGL_
	endif
endif
ifdef USE_SFINCS
	CPPFLAGS += -DUSE_SFINCS
endif

LIBS	+= $(DEFAULT_LIB) $(MPI_LIB) $(FFT_LIB) $(NETCDF_LIB) $(HDF5_LIB) \
		$(NAG_LIB) $(SFINCS_LIB) $(PETSC_LIB) $(LIBSTELL_LIB) \
		$(LAPACK_LIB)
F90FLAGS+= $(F90OPTFLAGS)
INC_FLAGS= $(DEFAULT_INC) $(MPI_INC) $(FFT_INC) $(NETCDF_INC) $(HDF5_INC) \
	   $(SFINCS_INC) $(PETSC_INC) $(LAPACK_INC)
CFLAGS += $(COPTFLAGS)

DATE=$(shell date +%y%m%d)
TARDIR=$(GK_PROJECT)_$(DATE)
TOPDIR=$(CURDIR)
ifeq ($(notdir $(CURDIR)), $(UTILS))
	TOPDIR=$(subst /$(UTILS),,$(CURDIR))
endif
ifeq ($(notdir $(CURDIR)), $(GEO))
	TOPDIR=$(subst /$(GEO),,$(CURDIR))
endif
ifneq ($(TOPDIR),$(CURDIR))
	SUBDIR=true
endif

VPATH = $(UTILS):$(GEO):$(VMEC)
# this just removes non-existing directory from VPATH
VPATH_tmp := $(foreach tmpvp,$(subst :, ,$(VPATH)),$(shell [ -d $(tmpvp) ] && echo $(tmpvp)))
VPATH = .:$(shell echo $(VPATH_tmp) | sed "s/ /:/g")
#
ifdef SUBDIR
	VPATH +=:..
endif
DEPEND=Makefile.depend
DEPEND_CMD=$(PERL) fortdep

# most common include and library directories
DEFAULT_INC_LIST = . $(UTILS) $(LIBSTELL) $(VMEC) $(GEO)
DEFAULT_LIB_LIST =
DEFAULT_INC=$(foreach tmpinc,$(DEFAULT_INC_LIST),$(shell [ -d $(tmpinc) ] && echo -I$(tmpinc)))
DEFAULT_LIB=$(foreach tmplib,$(DEFAULT_LIB_LIST),$(shell [ -d $(tmplib) ] && echo -L$(tmplib)))

# list of intermediate f90 files generated by preprocessor
F90FROMFPP = $(patsubst %.fpp,%.f90,$(notdir $(wildcard *.fpp */*.fpp)))

####################################################################### RULES

#.SUFFIXES:
.SUFFIXES: .fpp .f90 .c .o

.f90.o:
	$(FC) $(F90FLAGS) $(INC_FLAGS) -c $<
.fpp.f90:
	$(CPP) $(CPPFLAGS) $< $@
.F90.o:
	$(FC) $(F90FLAGS) $(CPPFLAGS) $(INC_FLAGS) -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<

##################################################################### TARGETS

# .DEFAULT_GOAL works for GNU make 3.81 (or higher)
# For 3.80 or less, see all target
.DEFAULT_GOAL := $(GK_PROJECT)_all
ifeq ($(notdir $(CURDIR)),utils)
	.DEFAULT_GOAL := utils_all
endif
ifeq ($(notdir $(CURDIR)),geo)
	.DEFAULT_GOAL := geo_all
endif

.PHONY: all $(GK_PROJECT)_all

all: $(.DEFAULT_GOAL)

include $(DEPEND)

sinclude Makefile.target_$(GK_PROJECT)

# Include unit test makefile, empty target so Make doesn't attempt to
# build the file
tests/unit/Makefile:
include tests/unit/Makefile

tests/integrated/Makefile:
include tests/integrated/Makefile

check: check-unit check-integrated

############################################################### SPECIAL RULES

# comment this out to keep intermediate .f90 files
#.PRECIOUS: $(F90FROMFPP)

.INTERMEDIATE: $(GK_PROJECT)_transforms.f90 $(GK_PROJECT)_io.f90 $(GK_PROJECT)_save.f90 \
		mp.f90 fft_work.f90 response_matrix.f90 multibox.f90 sources.f90 \
		fields.f90 mp_lu_decomposition.f90

############################################################# MORE DIRECTIVES

.PHONY: depend clean distclean tar test_make

depend:
	@$(DEPEND_CMD) -m "$(MAKE)" -1 -o -v=0 $(VPATH)

clean:
	-rm -f *.o *.mod *.g90 *.h core */core *~
	-rm -f $(GEO)/*.o $(GEO)/*~
	-rm -f Makefiles/*~
	-rm -f $(UTILS)/*.o $(UTILS)/*~
	$(MAKE) -C $(VMEC) clean

cleanlib:
	-rm -f *.a

distclean: unlink clean cleanlib

%.o: %.mod

tar:
	@[ ! -d $(TARDIR) ] || echo "ERROR: directory $(TARDIR) exists. Stop."
	@[ -d $(TARDIR) ] || $(MAKE) tar_exec

test_make:
	@echo GK_SYSTEM is $(GK_SYSTEM)
	@echo .DEFAULT_GOAL is $(.DEFAULT_GOAL)
	@echo VPATH is $(VPATH)
	@echo CURDIR is $(CURDIR)
	@echo TOPDIR is $(TOPDIR)
	@echo NETCDF_DIR is $(NETCDF_DIR)
	@echo FFTW_DIR is $(FFTW_DIR)
	@echo
	@echo Compile mode:
	@echo  DEBUG is $(DEBUG)
	@echo  TEST is $(TEST)
	@echo  OPT is $(OPT)
	@echo  STATIC is $(STATIC)
	@echo  DBLE is $(DBLE)
	@echo
	@echo Functions:
	@echo  USE_MPI is $(USE_MPI)
	@echo  USE_FFT is $(USE_FFT)
	@echo  USE_NETCDF is $(USE_NETCDF)
	@echo  USE_HDF5 is $(USE_HDF5)
	@echo  USE_LOCAL_RAN is $(USE_LOCAL_RAN)
	@echo  USE_LOCAL_SPFUNC is $(USE_LOCAL_SPFUNC)
	@echo  USE_NAGLIB is $(USE_NAGLIB)
	@echo  USE_LAPACK is $(USE_LAPACK)
	@echo  DEFAULT_LIB is $(DEFAULT_LIB)
	@echo  MPI_LIB is $(MPI_LIB)
	@echo
	@echo FC is $(FC)
	@echo F90FLAGS is $(F90FLAGS)
	@echo F90OPTFLAGS is $(F90OPTFLAGS)
	@echo CC is $(CC)
	@echo CFLAGS is $(CFLAGS)
	@echo COPTFLAGS is $(COPTFLAGS)
	@echo LD is $(LD)
	@echo LDFLAGS is $(LDFLAGS)
	@echo CPP is $(CPP)
	@echo CPPFLAGS is $(CPPFLAGS)
	@echo LIBS is $(LIBS)
	$(MAKE) -C $(VMEC) test_make

unlink:
	-rm -f $(F90FROMFPP)

revision:
	@LANG=C svn info | awk '{if($$1=="Revision:") printf("%20d",$$2) }' > Revision

TAGS:	*.f90 *.fpp */*.f90 */*.fpp
	etags $^

############################################################# Documentation

ifneq ("$(wildcard $(shell which $(FORD) 2>/dev/null))","")
check_ford_install:
	@echo "Using ford at $(shell which $(FORD))"
else
check_ford_install:
	@echo "Ford command $(FORD) not in path -- is it installed?\\n\\tConsider installing with 'pip install --user ford' and add ${HOME}/.local/bin to PATH" ; which $(FORD)
endif

GIT_VERSION := $(shell git describe --tags --long --match "v*" --first-parent HEAD)

doc: docs/stella_docs.md check_ford_install
	$(FORD) $(INC_FLAGS) -r $(GIT_VERSION) $<

cleandoc:
	@echo "FORD docs"
	-rm -rf docs/html
