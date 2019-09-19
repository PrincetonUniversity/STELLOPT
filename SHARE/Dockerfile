# Derived from docs.docker.com/get-started.part2

# Use an official Ubuntu runtime as a parent image
FROM ubuntu:latest

# Set the working directory to STELLOPT
WORKDIR /home/stellopt/STELLOPT

# Copy the current directory contents into the container at STELLOPT
COPY . /home/stellopt/STELLOPT

# Install all the stuff we need, build STELLOPT, and setup non-root user
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
    ./build_all 2>&1 | tee build.out && \
    useradd stellopt && chown -R stellopt:stellopt /home/stellopt

USER stellopt

# Run bash
CMD ["/bin/bash"]

