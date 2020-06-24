FROM ubuntu:latest
MAINTAINER Caoxiang Zhu <czhu@pppl.gov> & STELLOPT developers

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -q update && apt-get -y install \
    gfortran g++ libopenmpi-dev openmpi-bin \
    libnetcdf-dev libnetcdff-dev libhdf5-openmpi-dev hdf5-tools \
    libblas-dev liblapack-dev libscalapack-openmpi-dev \
    python3 python3-numpy python3-h5py make curl \
    git-all

WORKDIR /home/STELLOPT
COPY . /home/STELLOPT

# Set commands
CMD ["/bin/bash"]
