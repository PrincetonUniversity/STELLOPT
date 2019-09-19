# Derived from docs.docker.com/get-started.part2

# Use an official Ubuntu runtime as a parent image
FROM ubuntu:latest

# Set the working directory to STELLOPT
WORKDIR /STELLOPT

# Copy the current directory contents into the container at STELLOPT
COPY . /STELLOPT

# Install all the stuff we need and build STELLOPT
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    openmpi-bin libopenmpi-dev \
    libnetcdf-dev libnetcdff-dev \
    liblapack-dev \
    libblacs-mpi-dev \
    libscalapack-mpi-dev \
    libhdf5-openmpi-dev \
    pgplot5 \
    libncarg-dev \
    libgsl-dev \
    git \
    vim-tiny && \
    cp --remove-destination SHARE/make_docker.inc make.inc && \
    ./build_all 2>&1 | tee build.out

# Run bash
CMD ["/bin/bash"]

