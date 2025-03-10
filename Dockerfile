# Build with docker buildx build --tag libstell .
# Debian image of install
FROM lazerson/stellopt-compile:latest
LABEL version="1.0"
LABEL description="Docker build of LIBSTELL and intitial file copy."

# Set Environment variables
ARG CODE="LIBSTELL"
ARG MYHOME="/home/STELLOPT/bin"
ENV MACHINE="debian"
ENV STELLOPT_PATH=/home/STELLOPT
RUN echo Building ${CODE} for ${MACHINE} in Docker

# Set the working directory
WORKDIR $STELLOPT_PATH

# Copy the entire local directory to the container
COPY . $STELLOPT_PATH

# Compile STELLOPT
RUN cd $STELLOPT_PATH  && ./build_all -j1 $CODE 2>&1 | tee log.build

# Copy all built executables onto global path
RUN cp -RP ${STELLOPT_PATH}/bin/* /usr/local/bin/

# If you run this 
CMD ["/bin/bash"]
