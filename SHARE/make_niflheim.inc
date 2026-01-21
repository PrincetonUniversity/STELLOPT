CURRENT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
include $(CURRENT_DIR)/make_redhat.inc
LHDF5 = T
  HDF5_DIR = /home/modules/software/HDF5/1.14.6-gompi-2025b
  HDF5_INC = -I${HDF5_DIR}/include
  HDF5_LIB = -L${HDF5_DIR}/lib \
              -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5_hl_fortran
LIBS    = -L/home/modules/software/OpenBLAS/0.3.30-GCC-14.3.0/lib -lopenblas \
            -L/home/modules/software/ScaLAPACK/2.2.2-gompi-2025b-fb/lib -lscalapack
