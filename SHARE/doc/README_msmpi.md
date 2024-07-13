

You have to manually compile the fortran module to get the necessary *.mod files for your projects!

Installed MSMPI to MSYS in an MSYS2 terminal
$ pacman -S mingw-w64-x86_64-msmpi 

Here is a copy of my mingw64/lib/pkgconfig/msmpi.pc 
	prefix=/mingw64
	libdir=${prefix}/lib
	includedir=${prefix}/include
	Name: msmpi
	URL: https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi
	Version: 10.1.1
	Description: Microsoft MPI SDK (mingw-w64)
	Cflags: -I${includedir}
	Libs: -L${libdir} -lmsmpi

NOTE: My installed version of MSMPI (10.1.2) does not match that installed to MSYS

Here is 
==============================================
# NEWER VERSIONS OF MSMPI

1) Install MSMPI
Newest versions:
- https://www.microsoft.com/en-us/download/details.aspx?id=100593

install MSMPI from Microsoft.com ... this is for version 10.1.1
- https://download.microsoft.com/download/2/9/e/29efe9b1-16d7-4912-a229-6734b0c4e235/msmpisdk.msi
- https://download.microsoft.com/download/2/9/e/29efe9b1-16d7-4912-a229-6734b0c4e235/msmpisetup.exe

2) Make sure that mpiexec.exe is available its source is added to your windows path

3) Open a PowerShell or Windows Terminal  (not mingw or msys2!)

cd (your msys root)/mingw64/include   [for me this was C:/Programs/msys64/mingw64/include]

gfortran mpi.f90 -c -fno-range-check

ar cr libmpi.a mpi.o

	C:\Programs\msys64\mingw64\include>ls *mpi*
	H5FDmpi.h   hcompi.h  mpi.f90  mpi.mod  mpi_base.mod       mpi_sizeofs.mod  mpifptr.h
	H5FDmpio.h  libmpi.a  mpi.h    mpi.o    mpi_constants.mod  mpif.h

rm mpi.o 

mv libmpi.a ../lib/

	C:\Programs\msys64\mingw64\include>ls ../lib/*mpi*
	../lib/libmpi.a  ../lib/libmsmpi.a

==============================================
# THESE INSTRUCTIONS WORK WITH OLD VERSIONS OF MSMPI

cd ~/src/msmpi     # (or any worthy directory)
cp /mingw64/lib libmsmpi.a .
cp /mingw64/include/mpifptr.h .
cp /mingw64/include/mpif.h .
cp /mingw64/include/mpi.f90 .

gfortran -c -D_WIN64 -D INT_PTR_KIND\(\)=8 -fno-range-check mpi.f90

ls
example.F90  mpi.f90  mpi.o         mpi_constants.mod  mpif.h
libmsmpi.a   mpi.mod  mpi_base.mod  mpi_sizeofs.mod    mpifptr.h

gfortran -fno-range-check -ffree-form -o example.exe -D_WIN64 -D INT_PTR_KIND\(\)=8 example.F90 libmsmpi.a

ls
example.exe  libmsmpi.a  mpi.mod  mpi_base.mod       mpi_sizeofs.mod  mpifptr.h
example.F90  mpi.f90     mpi.o    mpi_constants.mod  mpif.h

./example.exe
STOP FAILED

mpiexec -n 2 example.exe
 PASSED

cp *.mod /mingw64/include/

=======================================

where example.F90 is:

program example

#if defined USE_MPI_MODULE
    use mpi
    implicit none
#else
    implicit none
#include "mpif.h"
#endif

    integer :: ierr, rank, num_ranks
    logical :: test_ok

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_ranks, ierr)
    call MPI_FINALIZE(ierr)
  
    test_ok = (num_ranks == 2 .and. (rank ==0 .or. rank == 1))

    if (test_ok) then
      if (rank == 0) print *, 'PASSED'
    else
      stop "FAILED"
    endif

end program