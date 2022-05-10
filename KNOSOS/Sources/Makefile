#################################################################################
############################# KNOSOS Makefile ####################################
##################################################################################


###### Detect system #############################################################

HOSTNAME := $(shell uname -n|cut -f1 -d.)
NODENAME := $(shell uname -n|rev|cut -f2 -dl|rev)
SHOSTNAME := $(shell uname -n|cut -f1 -du)

ifeq ($(SHOSTNAME), r000)
CASE=marconi
endif
ifeq ($(HOSTNAME), xula01)
CASE=xula
endif
ifeq ($(HOSTNAME), xula02)
CASE=xula
endif
ifeq ($(HOSTNAME), euler01)
CASE=euler
endif
ifeq ($(HOSTNAME), euler02)
CASE=euler
endif
ifeq ($(NODENAME), eu)
CASE=euler
endif
ifeq ($(HOSTNAME), dirac)
CASE=dirac
endif
ifeq ($(HOSTNAME), fusc3)
CASE=fusc3
endif
ifeq ($(HOSTNAME), macfus33)
CASE=desktop
EXT=Darwin
endif
ifeq ($(HOSTNAME), MBdJLVG)
CASE=desktop
EXT=Darwin
endif
ifeq ($(HOSTNAME), xUdJLVBG)
CASE=desktop
EXT=Linux
endif
ifeq ($(HOSTNAME), task3d-l)
CASE=task3d-l
endif
ifeq ($(HOSTNAME), ncd-ge0)
CASE=ncd-ge0
endif
ifeq ($(HOSTNAME), egcalc-vpn)
CASE=egcalc-vpn
endif
ifeq ($(HOSTNAME), draco01)
CASE=draco
endif
ifeq ($(HOSTNAME), draco02)
CASE=draco
endif


###### if running at marconi #####################################

ifeq ($(CASE), marconi)
COMPILER =mpiifort
#MPI_DIR  =
#PETSC_DIR=
#FFT_DIR  =
#MKL_DIR  =
#COMP_LIB =
NETCDF_DIR=${NETCDFF_HOME}
FFLAGS   =-fpp -DIPPorNIFS -DMPIandPETSc -O3 -xHost -shared-intel -warn all -mcmodel large --with-debugging=0
IFLAGS   =-I${NETCDF_DIR}/include/ -I${NETCDF_INCLUDE} -I$(PETSC_HOME)/include -I${LAPACK_HOME}/include -I${MKL_HOME}/include -I${PETSC_DIR}/lib/petsc/conf/variables -I${PETSC_DIR}/lib/petsc/conf/rules -I${PETSC_DIR}/include -I${PETSC_DIR}/include/petsc -qopenmp -I${PETSC_DIR}/include
LFLAGS   = -L${NETCDF_DIR}/lib/ -lnetcdff -L${NETCDF_HOME}/lib/ -lnetcdf -L$(PETSC_HOME)/lib -lpetsc -L${LAPACK_HOME}/lib -llapack -L${MKL_HOME}/lib/intel64 -lmkl_intel_lp64 
endif

###### if running at draco #####################################

ifeq ($(CASE), draco)
COMPILER =mpiifort
#MPI_DIR  =
#PETSC_DIR=
#FFT_DIR  =
#MKL_DIR  =
#COMP_LIB =
FFLAGS   =-fpp -DIPPorNIFS -DMPIandPETSc -O3 -xHost -shared-intel -warn all -mcmodel large --with-debugging=0
IFLAGS   =-I${PETSC_DIR}/lib/petsc/conf/variables -I${PETSC_DIR}/lib/petsc/conf/rules -I${PETSC_DIR}/include -I${PETSC_DIR}/include/petsc -qopenmp -I${PETSC_DIR}/include
LFLAGS   =-L${PETSC_DIR}/lib -lpetsc -L/${MKL_HOME}/lib/intel64 -mkl -Wl,-rpath,${MKL_HOME}/lib/intel64,-rpath,${PETSC_DIR}/lib
endif

###### if running at TASK3D-L #####################################

ifeq ($(CASE), task3d-l)
COMPILER =mpiifort
#MPI_DIR  =
PETSC_DIR=/usr/local/petsc-3.11.2
#FFT_DIR  =
#MKL_DIR  =
#COMP_LIB =
FFLAGS   =-fpp -DMPIandPETSc -DIPPorNIFS -O3 -xHost -shared-intel -warn all -mcmodel large --with-debugging=0 
IFLAGS   =-I${PETSC_DIR}/include -I${PETSC_DIR}/include/petsc 
LFLAGS   =-mkl -lpetsc -L${PETSC_DIR}/lib 
endif

###### if running at mimizuku #####################################
ifeq ($(CASE), ncd-ge0)
COMPILER =mpif90
#MPI_DIR  =
PETSC_DIR=/usr/local/petsc
#FFT_DIR  =
#MKL_DIR  =
#COMP_LIB =
NETCDF=/usr/local/netcdf4
HDF5=/usr/local/hdf5
FFLAGS   =-fpp -DMPIandPETSc -DIPPorNIFS -O3 -xHost -shared-intel -warn all -mcmodel large --with-debugging=0
IFLAGS   =-I${PETSC_DIR}/include -I${PETSC_DIR}/include/petsc -I$(NETCDF)/include 
LFLAGS   =-mkl=sequential -lpetsc -L${PETSC_DIR}/lib -L$(NETCDF)/lib -lnetcdff -lnetcdf -L$(HDF5)/lib -lhdf5 -lhdf5_hl
endif

###### if running at egcalc-vpn #####################################

ifeq ($(CASE), egcalc-vpn)
COMPILER =mpif90
MPI_DIR  =/home/LHD/satake/local/lib
PETSC_DIR=/home/LHD/satake/local/petsc-3.11.2
FFT_DIR  =
MKL_DIR  =
COMP_LIB =
FFLAGS   =-fpp -DMPIandPETSc -DIPPorNIFS -O3 -xHost -shared-intel -warn all -mcmodel large -align array64byte --with-debugging=0
#FFLAGS   =-fpp -O3 -xHost -shared-intel -warn all -mcmodel large --with-debugging=0
IFLAGS   =-I${PETSC_DIR}/include -I${PETSC_DIR}/include/petsc
LFLAGS   =-mkl -lpetsc -L${PETSC_DIR}/lib
endif

###### if running at Xula cluster at CIEMAT #####################################

ifeq ($(CASE), xula)
COMPILER =mpif90
#MPI_DIR  =
#PETSC_DIR=
#FFT_DIR  =
#MKL_DIR  =
#COMP_LIB =
NETCDF_DIR=/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-17.0.8/netcdf-fortran-4.4.4-4me7upiib2q3vpkd2qzvtuaxpruy7ndo
#NECDF_DIR=${NETCDF}
#FFLAGS   =-fpp -DIPPorNIFS -DMPIandPETSc -check all -O3 -xHost -shared-intel -warn all -mcmodel large --with-debugging=0 --with-scalar-type=real --with-precision=double
FFLAGS   =-fpp -DIPPorNIFS -DMPIandPETSc -O3 -xHost -shared-intel -warn all -mcmodel large --with-debugging=0 --with-scalar-type=real --with-precision=double
IFLAGS   =-I${PETSC_DIR}/lib/petsc/conf/variables -I${PETSC_DIR}/lib/petsc/conf/rules -I${PETSC_DIR}/include -I${PETSC_DIR}/include/petsc -qopenmp -I$(NETCDF_DIR)/include 
LFLAGS   =-L${PETSC_DIR}/lib -lpetsc -L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf -L/${MKL_HOME}/lib/intel64 -mkl -Wl,-rpath,${MKL_HOME}/lib/intel64,-rpath,${PETSC_DIR}/lib
endif


###### if running at euler cluster at CIEMAT #####################################

ifeq ($(CASE), euler)
COMPILER =mpif90 -traceback
MPI_DIR  =/opt/ofed_1.3.1/mpi/intel/openmpi-1.2.6/lib64/
PETSC_DIR=/home/localsoft/petsc-3.7.4/real
FFT_DIR  =/home/localsoft/intel-12/composerxe-2011.3.174/mkl/include/fft_w
MKL_DIR  =/home/roldan/gcc/520euler/lib64/
NETCDF=/soft/netcdf-4.0.1-intel/include
NETCDF2=/home/localsoft/netcdf-4.0.1-intel/lib
#COMP_LIB =
FFLAGS   =-cpp -DMPIandPETSc=.TRUE. -O3 -shared-intel -warn all -mcmodel large --with-debugging=0 
IFLAGS   =-I${PETSC_DIR}/include -I${PETSC_DIR}/include/petsc -I${FFT_DIR} -I${NETCDF} 
LFLAGS   =-mkl -L${MKL_DIR} -lpetsc -L${PETSC_DIR}/lib -L${MPI_DIR} -lnetcdf -L${NETCDF2} 
endif

###### if running at dirac computer at CIEMAT #####################################

ifeq ($(CASE), dirac)
COMPILER =ifort
#MPI_DIR  =
#PETSC_DIR=
#FFT_DIR  =
#MKL_DIR  =
#COMP_LIB =
FFLAGS   =-cpp -O3 -shared-intel -warn all
#IFLAGS   =
LFLAGS   =-mkl
endif

###### if running at fusc computer at CIEMAT #####################################

ifeq ($(CASE), fusc3)
COMPILER =ifort
#MPI_DIR  =
#PETSC_DIR=
#FFT_DIR  =
#MKL_DIR  =
COMP_LIB =LIBS/libnag_i4.a
FFLAGS   =-cpp -DNAG=.TRUE. -O3 -shared-intel -warn all
#IFLAGS   =
LFLAGS   =-mkl
endif

###### if running at Velasco's eMAC, Macbook or xUbuntu at CIEMAT ################

ifeq ($(CASE), desktop)
COMPILER =gfortran -fbounds-check
#MPI_DIR  =
#PETSC_DIR=
#FFT_DIR  =
#MKL_DIR  =
COMP_LIB =LIBS/libfftw3_${EXT}.a
FFLAGS   =-cpp --warn-all 
#IFLAGS   =
LFLAGS   =-llapack -lblas
endif

###### Targets, pre-requisites and commands ######################################

knosos.x: \
	global.o knosos_stellopt.o \
	amb_and_qn.o ana_ntv.o banana.o coefficients.o configuration.o database.o  \
	errors.o fast_ion_confinement.o input.o knosos.o lagrange.o low_collisionality.o \
	low_collisionality_new.o mpi_and_end.o \
	plasmas.o plateau.o pfirsch_schlueter.o times.o trace_impurities.o transport.o

	${COMPILER} -o knosos.x global.o knosos_stellopt.o \
	amb_and_qn.o ana_ntv.o banana.o coefficients.o configuration.o database.o \
	errors.o fast_ion_confinement.o input.o knosos.o lagrange.o low_collisionality.o \
	low_collisionality_new.o mpi_and_end.o \
	plasmas.o plateau.o pfirsch_schlueter.o times.o trace_impurities.o transport.o \
	${COMP_LIB} ${LFLAGS} ${IFLAGS} ${FFLAGS} 

global.o: global.f90
	${COMPILER} ${FFLAGS} -c global.f90

knosos_stellopt.o: knosos_stellopt.f90
	${COMPILER} ${FFLAGS} -c knosos_stellopt.f90

amb_and_qn.o: amb_and_qn.f90
	${COMPILER} ${FFLAGS} -c amb_and_qn.f90

ana_ntv.o: ana_ntv.f90
	${COMPILER} ${FFLAGS} -c ana_ntv.f90

banana.o: banana.f90
	${COMPILER} ${FFLAGS} ${LFLAGS} -c banana.f90

coefficients.o: coefficients.f90
	${COMPILER} ${FFLAGS} -c coefficients.f90

configuration.o: configuration.f90
	${COMPILER} ${FFLAGS} ${IFLAGS} -c configuration.f90

database.o: database.f90
	${COMPILER} ${FFLAGS} -c database.f90

errors.o: errors.f90
	${COMPILER} ${FFLAGS} -c errors.f90

input.o: input.f90
	${COMPILER} ${FFLAGS} -c input.f90

knosos.o: knosos.f90
	${COMPILER} ${FFLAGS} ${IFLAGS} -c knosos.f90

lagrange.o: lagrange.f90
	${COMPILER} ${FFLAGS} -c lagrange.f90

low_collisionality.o: low_collisionality.f90
	${COMPILER} ${FFLAGS} ${IFLAGS} -c low_collisionality.f90

low_collisionality_new.o: low_collisionality_new.f90
	${COMPILER} ${FFLAGS} ${IFLAGS} -c low_collisionality_new.f90

fast_ion_confinement.o: fast_ion_confinement.f90
	${COMPILER} ${FFLAGS} ${IFLAGS} -c fast_ion_confinement.f90

mpi_and_end.o: mpi_and_end.f90
	${COMPILER} ${FFLAGS} -c mpi_and_end.f90

plasmas.o: plasmas.f90
	${COMPILER} ${FFLAGS} -c plasmas.f90

plateau.o: plateau.f90
	${COMPILER} ${FFLAGS} -c plateau.f90

pfirsch_schlueter.o: pfirsch_schlueter.f90
	${COMPILER} ${FFLAGS} ${IFLAGS} -c pfirsch_schlueter.f90

times.o: times.f90
	${COMPILER} ${FFLAGS} -c times.f90

trace_impurities.o: trace_impurities.f90
	${COMPILER} ${FFLAGS} -c trace_impurities.f90

transport.o: transport.f90
	${COMPILER} ${FFLAGS} -c transport.f90

clean:
	rm knosos.x *genmod* *.o global.mod knosos_stellopt.mod
