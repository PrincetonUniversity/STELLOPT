STELLOPT Compilation with Docker
==============================

This page details how to compile the STELLOPT family of codes using
[Docker](@https://www.docker.com/). In order to do this you'll need
to install the Docker Desktop program.

Base Image
-----

A base Debian image to build upon is provided [stellopt-compile](@https://hub.docker.com/repository/docker/lazerson/stellopt-compile/general)

This image contains the basic libraries to compile all the codes in STELLOPT.

LIBSTELL Image
-----

In the Git repository you'll fine a Dockerfile which allows one to 
build LIBSTELL and copy in most of the files. You can build this 
image with the command

    docker buildx build --file Dockerfile --tag libstell

It is important you set the tag as libstell.

STELLOPT Image
-----
In the Git repository you'll find a Dockerfile_STELLOPT which allows
you to build the entire STELLOPT suite using

    docker buildx build --file Dockerfile_STELLOPT --tag stellopt

Updating a single code
-----

Sometimes you may want to rebuild a single code after modification.
To do this we provide the Dockerfile_single_code file. The file looks like

    # Your stellopt image which should match your tag
    FROM stellopt:latest
    LABEL description="STELLOPT build for a specific code"
    # Set what code you wish to compile
    ARG CODE="DKES"
    # Set Environment variables
    ARG MYHOME="/home/STELLOPT/bin"
    RUN echo Building ${CODE} for ${MACHINE}
    # Set the working directory
    WORKDIR $STELLOPT_PATH
    # Copy the local code modifications to STELLOPT
    COPY ${CODE} ${STELLOPT_PATH}
    # Compile STELLOPT
    RUN cd $STELLOPT_PATH  && ./build_all -j1 $CODE 2>&1 | tee log.build 
    # Copy all built executables onto global path
    RUN cp -RP ${STELLOPT_PATH}/bin/* /usr/local/bin/
    # Set commands
    CMD ["/bin/bash"]

The line CODE="DKES" can be changed to whichever code you wish to
build. It can be built with (use proper tag name convention)

    docker buildx build --file Dockerfile_single_code --tag dkes  

Running a code
-----

You can enter the image to run codes using

    docker run -it stellopt:latest bash

Note on multiple code builds
-----

If you wish to build multiple codes, this is possible using the CODE argument.
However, keep in mind that the COPY command needs to be modified. For example to 
compile DKES and PENTA.  The file would need to look like:

    # Your stellopt image which should match your tag
    FROM stellopt:latest
    LABEL description="STELLOPT build for a specific code"
    # Set what code you wish to compile
    ARG CODE="DKES PENTA"
    # Set Environment variables
    ARG MYHOME="/home/STELLOPT/bin"
    RUN echo Building ${CODE} for ${MACHINE}
    # Set the working directory
    WORKDIR $STELLOPT_PATH
    # Copy the local code modifications to STELLOPT
    COPY DKES ${STELLOPT_PATH}
    COPY PENTA ${STELLOPT_PATH}
    # Compile STELLOPT
    RUN cd $STELLOPT_PATH  && ./build_all -j1 $CODE 2>&1 | tee log.build 
    # Copy all built executables onto global path
    RUN cp -RP ${STELLOPT_PATH}/bin/* /usr/local/bin/
    # Set commands
    CMD ["/bin/bash"]
