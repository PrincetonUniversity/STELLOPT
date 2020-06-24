FROM ubuntu:latest
MAINTAINER Caoxiang Zhu <czhu@pppl.gov> & STELLOPT developers

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -q update && apt-get -y install \
    gfortran g++ libopenmpi-dev openmpi-bin \
    libnetcdf-dev libnetcdff-dev libhdf5-openmpi-dev hdf5-tools \
    libblas-dev liblapack-dev libscalapack-openmpi-dev \
    python3 python3-pip \
    make curl git-all

# Install required python packages
RUN python3 -m pip install --upgrade pip
RUN pip install numpy scipy matplotlib notebook h5py xarray coilpy

# Set commands
CMD ["/bin/bash"]
