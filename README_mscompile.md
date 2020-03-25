This is an experimental branch for compiling LIBSTELL under msys-mingw 
on windows.

=========================
Current issue:

gfortran -fPIC -O2 -fexternal-blas -I/home/weir/bin/libstell_dir -I/mingw64/include -L/mingw64/include -I/mingw64/include -I/mingw64/include -L/mingw64/include -I/mingw64/include -I.. -c ../Sources/Modules/vmec_input.f
vmec_input.c:262:47-68:

Warning: Array reference at (1) out of bounds (0 < 1) in loop beginning at (2)
/mingw64/bin/cpp -traditional -DWIN -DNETCDF -DFFTW3 -DLHDF5 -DH5_USE_16_API ../Sources/Modules/mgrid_mod.f
gfortran -fPIC -O2 -fexternal-blas -I/home/weir/bin/libstell_dir -I/mingw64/include -L/mingw64/include -I/mingw64/include -I/mingw64/include -L/mingw64/include -I/mingw64/include -I.. -c ../Sources/Modules/mgrid_mod.f
mgrid_mod.c:783:19:

Error: Symbol 'temp_rank' at (1) has no IMPLICIT type
mgrid_mod.c:385:40:

Error: Symbol 'mpi_integer' at (1) has no IMPLICIT type
make[2]: *** [makelibstell:38: mgrid_mod.o] Error 1
make[2]: Leaving directory '/home/weir/src/stellopt/LIBSTELL/Release'
make[1]: *** [makefile:20: release] Error 2
make[1]: Leaving directory '/home/weir/src/stellopt/LIBSTELL'
make: *** [makefile:41: clean_release] Error 2


- This error appears to be pop up because the MPI environment doesn't declare
the data type for the "comm" variable in the new mpi_inc 



=========================
Environment setup:
=========================


Install MSYS with mingw64 and GCC9:

http://sourceforge.net/projects/msys2/files/Base/

Open an MSYS2 shell:
pacman -Syuu   # update the package list, will take a while
(you'll have to exit the installer)

close the shell and reopen the MSYS2 shell:
pacman -Syuu   # Run it a second time

close the shell and reopen the MSYS2 shell:
pacman -S base-devel git mingw-w64-x86_64-toolchain   # will take a while

# Make sure some other stuff is installed too
pacman -S wget

close the shell and open an MSYS2-MINGW64 shell:
gcc -v    # test gcc

Add /msys64/mingw64 and /msys64 to your user's environment PATH variable

============

Install the background modules that are necessary / generally helpful

pacman -S mingw-w64-x86_64-fltk mingw-w64-x86_64-f2c
pacman -S mingw-w64-x86_64-lapack mingw-w64-x86_64-openblas
pacman -S mingw-w64-x86_64-gsl mingw-w64-x86_64-fgsl mingw-w64-x86_64-fftw
pacman -S mingw-w64-x86_64-msmpi  #  SEE NOTE ON MPI BELOW
pacman -S mingw-w64-x86_64-hdf4 mingw-w64-x86_64-hdf5

============

SCALAPACK is available from a developer that is pushing to the main MSYS-REPO

cd ~
mkdir src
cd src
git clone git@github.com:okhlybov/MINGW-packages MINGW-branched
git checkout scalapack
cd MINGW-branched/mingw-w64-scalapack
makepkg-mingw -sCLf
pacman -U mingw-w64-x86_64-scalapack-2.0.2-1-any.pkg.tar.zst

=====================

I had to write a few custom builds to make netcdf-c, netcdf-cxx, and 
netcdf-fortran to compile together.

These are located under the STELLOPT/SHARE/pkgbuilds
(WARNING: I set the prefix for installation to /mingw64)

I didn't bother to write new check-sums for the files, so go to each 
directory, and type this:

cd SHARE/pkgbuilds/mingw-w64-netcdf
makepkg-mingw -sCLf --skipchecksums # download everything, configure and compile it, then create a package
pacman -U mingw-w64-x86_64-netcdf-4.7.3-1-any.pkg.tar.zst  # installs

cd SHARE/pkgbuilds/mingw-w64-netcdf-cxx
makepkg-mingw -sCLf --skipchecksums
pacman -U mingw-w64-x86_64-netcdf-cxx-4.3.0-1-any.pkg.tar.zst

cd SHARE/pkgbuilds/mingw-w64-netcdf-fortran
makepkg-mingw -sCLf --skipchecksums
pacman -U mingw-w64-x86_64-netcdf-fortran-4.5.2-1-any.pkg.tar.zst

NOTE: the netcdf-fortran installation might throw an error because HDF5 includes a netcdf header in the same spot
- I used the netcdf-fortran one and overwrote the HDF5 one 
(etc. In principal, the netcdf-fortran one contains the stuff from netcdf-c already and is the more complete of the two)


====================================================================================

Now you can go and try to run the build_msys script on this branch of STELLOPT
-- change the user to your local user name: replace mine "weir" with yours
-- You also have to copy the "awk_cdir.awk" file to the folder above this one (home dir)
-- the other changes I have made in this branch should be with you already...

====================================================================================
====================================================================================
NOTE ON MPI: Get ready to test openmpi and microsoft-mpi

- GCC-9 comes with its own implementation of openmp
- you can test that installation with this shamelessly stolen test code (works for me)

cd ~
mkdir omp_hello
cd omp_hello

wget https://computing.llnl.gov/tutorials/openMP/samples/C/omp_hello.c
gcc -fopenmp omp_hello.c -o omp_hello     # generate the executable file omp_hello.exe

run it:
./omp_hello       # By default, gcc creates 1 thread for each core. My PC has 2 physical cores (4 logical cores under hyperthreading)

# The order of the five output lines will be random
Hello World from thread = 1
Hello World from thread = 2
Hello World from thread = 3
Hello World from thread = 0
Number of threads = 4

export OMP_NUM_THREADS=8  # explicitly set 8 threads by the OMP_NUM_THREADS environment variable
./omp_hello
Hello World from thread = 3
Hello World from thread = 0
Number of threads = 8
Hello World from thread = 6
Hello World from thread = 5
Hello World from thread = 4
Hello World from thread = 1
Hello World from thread = 2
Hello World from thread = 7

======

You can either use the MSYS version of MSMPI, or you can install it yourself.  I have tested both, and there
was no difference in the IMPLICIT type error above. The tests work for both.
----- the MSYS version comes with headers for fortran mpi 

Instructions for installing MS-MPI yourself:

1) Download ms-mpi from microsoft (version 7): https://www.microsoft.com/en-us/download/details.aspx?id=52981

2) install both the SDK and MPI ( msmpisdk.msi and MSMpiSetup.exe )

3) Close and reopen your MSYS2-MINGW64 shell:

printenv | grep "WIN\|MSMPI"
.. should list the environment variables for ms-mpi

Copy the headers and libraries into your msys2 installation:

mkdir ~/msmpi                     # create a temporary folder under your home directory
cd ~/msmpi                        # enter the folder
cp "$MSMPI_LIB64/msmpi.lib" .     # copy msmpi.lib to ~/msmpi/; the import library, which is a placeholder for dll
cp "$WINDIR/system32/msmpi.dll" . # copy msmpi.dll to ~/msmpi/; the runtime library
gendef msmpi.dll                  # generate msmpi.def. For 32-bit, use: gendef -a msmpi.dll, which specifies the stdcall format
dlltool -d msmpi.def -D msmpi.dll -l libmsmpi.a   # generate the (static) library file libmsmpi.a
cp libmsmpi.a /mingw64/lib        # copy this library file to where g++ looks for them;
                                    # try "g++ --print-search-dirs"
cp "$MSMPI_INC/mpi.h" .           # copy the header file mpi.h to ~/msmpi/

open mpi.h under ~/msmpi in an editor, look for “typedef __int64 MPI_Aint”. Just above it, add a new line with “#include <stdint.h>” (without the quotes), which we need for the definition __int64.

cd ~
mkdir mpi_hello
cd mpi_hello

4) Test your msmpi installation

wget https://raw.githubusercontent.com/wesleykendall/mpitutorial/gh-pages/tutorials/mpi-hello-world/code/mpi_hello_world.c
gcc mpi_hello_world.c -lmsmpi -o mpi_hello	# -lmsmpi: links with the msmpi library, the file libmsmpi.a that we generated above

./mpi_hello                       # run it with 1 process
export PATH="$MSMPI_BIN":$PATH    # add MSMPI_BIN (where mpiexec.exe is) to PATH
mpiexec -n 4 mpi_hello.exe        # run it with 4 processes

======

If you run into problems, you might have to manually install some items if you have a weird system:

Install CMAKE

open a new msys2 shell:

pacman -S  mingw-w64-x86_64-cmake                     # install CMake

open a new mingw64 shell and test CMAKE:

cd ~
git clone https://github.com/bast/cmake-example.git   # download a demo code
cd cmake-example
git checkout 7931bf4                                # revert to a 2016 version; (courtesy of Filip Sund) the 2018 version will let
                                                    # CMake download the 2018 googletest, which has a bug in "internal\gtest-port.h."
mkdir build                                         # you cannot build a project in the source-code folder, so we create a subfolder.
cd build                                            # ensure this folder is empty, no previous CMake* files or folder.
cmake .. -G "MSYS Makefiles" -DCMAKE_INSTALL_PREFIX=$MINGW_PREFIX     # ".." locates CMakeLists.txt, "MSYS Makefiles" is needed MSYS2,                                                                       # -D... tells "make install" where to install files.
ls                                                  # CMake creates many files including Makefile for MSYS2's make command
make                                                # create three exe files
ls *.exe
./hello.x                                           # you should get Hello World
./main.x
./unit_tests

====================================================================
Compile blas from source as a test that your system works:

wget http://www.netlib.org/blas/blas-3.6.0.tgz
tar xf blas-3.6.0.tgz
cd BLAS-3.6.0
gfortran.exe -c *.f       # compile each .f file and produce a .o fils
                            # (you can also add the Optimization switch -O3).
ar rv libblas.a *.o       # combine all the .o files into a library file.
cp libblas.a /mingw64/lib # copy the library file to where g++ looks for them;
                            # try "g++ --print-search-dirs"

Compile a dll to make sure you can do that too:
gfortran.exe -shared -o libblas.dll -Wl,--out-implib=libblas.dll.a -Wl,--export-all-symbols -Wl,--enable-auto-import -Wl,--whole-archive *.o -Wl,--no-whole-archive -lgfortran

cp libblas.dll /mingw64/lib/

rm *.o

Overwrite the blas file that you installed with the faster version: Openblas 
--- make sure to close your MSYS2 MINGW64 shell, then open an MSYS2 shell:
pacman -S mingw-w64-x86_64-openblas

close your msys2 shell

=======

Install LAPACK AND BLAS

Open a new msys2 shell:

git clone https://github.com/msys2/MINGW-packages.git # clone the scripts
cd MINGW-packages/mingw-w64-lapack                    # locate the LAPACK script

- edit the PKGBUILD file to use version 3.8 (3.9 has an error):
change "pkgver=3.9.0" to "pkgver=3.8.0"

Now the hash won't match anymore, so install it without checking the keys:
makepkg-mingw --skipchecksums        # build BLAS and LAPACK
pacman -U mingw-w64-x86_64-lapack*   # install BLAS and LAPACK

close your msys2 shell and open a mingw64 shell:
cd ~
wget https://www.math.ucla.edu/~wotaoyin/software/lapack_test.cpp
g++ lapack_test.cpp -llapack -o lapack_test     # build
./lapack_test                                   # run

close your msys2 shell and your mingw64 shell


If you choose to do this you  should reinstall openblas + scalapack to make 
sure you have optimized libraries.

===============================================================
On-going projects:


SILO    - UNSUCCESSFUL SO FAR: stupid symlinks in tarballs break my PKGBUILD solution.
... working around that problem caused another of course: 
        - SILO tries to identify the computer by searching a bunch of lists. Doesn't catch our environment.

# From scratch (Work-in-progress):
wget https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2-bsd.tar.gz

edit the "src/Makefile.am" before configuration:
- replace the lines:   
     libsilo_la_LDFLAGS = -avoid-version 
with:
     libsilo_la_LDFLAGS = -no-undefined -avoid-version 

mkdir build
cd build

# Without a prefix, the default is the build/bin etc. directories
# if you want to strip the unused variables, that is where you should do it pre-install
../configure --enable-shared --enable-static --enable-fortran
#../configure --prefix=/mingw64 --enable-shared

make

ERROR:
make[4]: Leaving directory '/home/weir/src/silo-4.10.2-bsd/build/src'
make[3]: Leaving directory '/home/weir/src/silo-4.10.2-bsd/build/src'
make[2]: Leaving directory '/home/weir/src/silo-4.10.2-bsd/build/src'
Making all in tools
make[2]: Entering directory '/home/weir/src/silo-4.10.2-bsd/build/tools'
Making all in .
make[3]: Entering directory '/home/weir/src/silo-4.10.2-bsd/build/tools'
make[3]: Nothing to be done for 'all-am'.
make[3]: Leaving directory '/home/weir/src/silo-4.10.2-bsd/build/tools'
Making all in silock
make[3]: Entering directory '/home/weir/src/silo-4.10.2-bsd/build/tools/silock'
gcc -DHAVE_CONFIG_H -I. -I../../../tools/silock -I../..  -I../../src/silo -I../../../src/silo -I/mingw64/include   -fPIC -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -Wdeclaration-after-statement -MT silock.o -MD -MP -MF .deps/silock.Tpo -c -o silock.o ../../../tools/silock/silock.c
../../../tools/silock/silock.c: In function 'handleInvalidValue':
../../../tools/silock/silock.c:116:7: error: unknown type name 'fpclass_t'; did you mean 'fpclassify'?
  116 |    {  fpclass_t theClass = fpclass(value);
      |       ^~~~~~~~~
      |       fpclassify
make[3]: *** [Makefile:416: silock.o] Error 1
make[3]: Leaving directory '/home/weir/src/silo-4.10.2-bsd/build/tools/silock'
make[2]: *** [Makefile:387: all-recursive] Error 1
make[2]: Leaving directory '/home/weir/src/silo-4.10.2-bsd/build/tools'
make[1]: *** [Makefile:437: all-recursive] Error 1
make[1]: Leaving directory '/home/weir/src/silo-4.10.2-bsd/build'
make: *** [Makefile:366: all] Error 2


make check


make install


PETSC    - UNSUCCESSFUL SO FAR, still working on SILO 

--- You have to install python2 for this package ---

wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.12.4.tar.gz


======

openmpi issues:

wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.3.tar.gz

from the openmpi-4.0.3 README:
 The symbols that now no longer appear by default in Open MPI's mpi.h
  are:

  - MPI_Address (replaced by MPI_Get_address)
  - MPI_Errhandler_create (replaced by MPI_Comm_create_errhandler)
  - MPI_Errhandler_get (replaced by MPI_Comm_get_errhandler)
  - MPI_Errhandler_set (replaced by MPI_Comm_set_errhandler)
  - MPI_Type_extent (replaced by MPI_Type_get_extent)
  - MPI_Type_hindexed (replaced by MPI_Type_create_hindexed)
  - MPI_Type_hvector (replaced by MPI_Type_create_hvector)
  - MPI_Type_lb (replaced by MPI_Type_get_extent)
  - MPI_Type_struct (replaced by MPI_Type_create_struct)
  - MPI_Type_ub (replaced by MPI_Type_get_extent)
  - MPI_LB (replaced by MPI_Type_create_resized)
  - MPI_UB (replaced by MPI_Type_create_resized)
  - MPI_COMBINER_HINDEXED_INTEGER
  - MPI_COMBINER_HVECTOR_INTEGER
  - MPI_COMBINER_STRUCT_INTEGER
  - MPI_Handler_function (replaced by MPI_Comm_errhandler_function)

  Although these symbols are no longer prototyped in mpi.h, they
  are still present in the MPI library in Open MPI v4.0.1 and later
  releases of the v4.0.x release stream. This enables legacy MPI
  applications to link and run successfully with
  Open MPI v4.0.x, even though they will fail to compile.

  All that being said, if you are unable to immediately update your
  application to stop using these legacy MPI-1 symbols, you can
  re-enable them in mpi.h by configuring Open MPI with the
  --enable-mpi1-compatibility flag.

  NOTE: Open MPI v4.0.0 had an error where these symbols were not
        included in the library if configured without --enable-mpi1-compatibility
        (see https://github.com/open-mpi/ompi/issues/6114).
        This is fixed in v4.0.1, where --enable-mpi1-compatibility
        flag only controls what declarations are present in the MPI header.


./configure CC=gcc CXX=g++ FC=gfortran --prefix=<directory>  --enable-static --enable-shared --enable-mpi-fortran=all
--with-slurm --enable-mpi-cxx PKG_CONFIG=/mingw64/bin/pkg-config

-finline-functions -fno-strict-aliasing added to CFLAGS
-finline-functions aded to CPPFLAGS

INTEGER*16 MPI datatype unsupported

pacman -S mingw-w64-x86_64-perl

pacman -S mingw-w64-cross-binutils




