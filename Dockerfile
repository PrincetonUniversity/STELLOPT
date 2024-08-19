# Debian image of install
FROM lazerson/stellopt-compile:latest
LABEL maintainer="Samuel Lazerson <lazersos@gmail.com> & STELLOPT developers"
LABEL version="1.0"
LABEL description="Debian image for building STELLOPT on docker."

# Set the working directory
WORKDIR /home/STELLOPT

# Copy the current directory contents into the container at /usr/src/app
COPY . /home/STELLOPT

# Set Environment variables
ARG MACHINE="debian"
ARG MYHOME="/home/STELLOPT/bin"
ENV STELLOPT_PATH=/home/STELLOPT
RUN echo $STELLOPT_PATH

# Compile STELLOPT
RUN cd $STELLOPT_PATH  && ./build_all -j1 2>&1 | tee log.build 
RUN chmod -R 777 ${STELLOPT_PATH}/BENCHMARKS
RUN cp -RP ${STELLOPT_PATH}/bin/* /usr/local/bin/

# add user
#RUN apt-get -y install sudo
#RUN useradd visitor && echo "visitor:visitor" | chpasswd && adduser visitor sudo
#WORKDIR /home/visitor
#USER visitor

# Set commands
#CMD ["/bin/bash"]
