FROM ubuntu:latest
MAINTAINER Samuel Lazerson <lazersos@gmail.com> & STELLOPT developers

# Install gfortran
RUN apt-get update
RUN apt-get install -y gfortran
RUN apt-get install -y make 
RUN apt-get install -y openmpi
RUN apt-get install -y git

# Set the working directory
WORKDIR /home/STELLOPT

# Copy the current directory contents into the container at /usr/src/app
COPY . /home/STELLOPT

# Compile STELLOPT
ENV MACHINE="docker"
ENV STELLOPT_PATH=/home/STELLOPT
RUN echo $STELLOPT_PATH
RUN cd $STELLOPT_PATH  && ./build_all -j4 2>&1 | tee log.build 
RUN chmod -R 777 ${STELLOPT_PATH}/BENCHMARKS
RUN cp -RP ${STELLOPT_PATH}/bin/* /usr/local/bin/

# add user
#RUN apt-get -y install sudo
#RUN useradd visitor && echo "visitor:visitor" | chpasswd && adduser visitor sudo
#WORKDIR /home/visitor
#USER visitor

# Set commands
#CMD ["/bin/bash"]
