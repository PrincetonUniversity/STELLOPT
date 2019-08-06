STELLOPT Compilation on Centos
==============================

This page details how to compile the STELLOPT family of codes on
[CentOS](@https://www.centos.org/). In order to do so you will need to
install [GNU based compilers](@http://gcc.gnu.org/) on your Linux
machine. In principle this set of steps should work for any Linux based
distribution which uses the
[Redhat Package Manager (RPM)](@http://rpm.org/) based package
management.

------------------------------------------------------------------------

Installation Steps

1.  Verify that your version of CentOS is up to date.
2.  Use [yum](@http://yum.baseurl.org/) to install your packages. For
    this example we\'re using gfortran 4.4 but you may substitute any
    other version.

<!-- -->

    yum install gcc*
    yum install openmpi*
    yum install lapack*

1.  Now we need to install external packages. Pleas visit the
    [Fedora page](@http://fedoraproject.org/wiki/EPEL) for the latest
    stable version of the EPEL.

<!-- -->

    wget http://mirrors.mit.edu/epel/6/i386/epel-release-6-8.noarch.rpm
    rpm -Uvh epel-release-6-8.noarch.rpm
    yum -assumeyes install netcdf*
    yum install ncl*
    yum install hdf5*

1.  Before compiling you\'ll need to load openmpi

<!-- -->

    module load openmpi-x86_64

1.  Now create a directory in which you wish to compile STELLOPT and
    place the STELLOPT ZIP file in that directory. Then unzip it.
2.  You\'ll need to install the
    [PSPLINE library](@http://w3.pppl.gov/ntcc/PSPLINE/) so create a new
    directory call PSPLINE and download the PSPLINE tarball into it.
    Then run the following commands (note use gmake, make causes
    problems)

<!-- -->

    mkdir PSPLINE
    cd PSPLINE
    cp <PATH_TO_PSLINE>/pspline.tar.gz .
    tar xvf pspline.tar.gz
    gmake FORTRAN_VARIANT=GCC LIBDIR=/usr/lib64 NETCDF_DIR=/usr/

1.  Unzip stellinstall
2.  Now you should be able to run the setup script and compile the
    STELLOPT family of codes.
