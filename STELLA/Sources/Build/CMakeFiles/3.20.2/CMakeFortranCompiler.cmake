set(CMAKE_Fortran_COMPILER "/usr/bin/f95")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "4.8.5")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/usr/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/usr/bin/gcc-ranlib")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/usr/lib/gcc/x86_64-redhat-linux/4.8.5/finclude;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/netlib-lapack-3.8.0-fcypwwwujl67z2sdytr3dpok3auhlgtu/include;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/libxml2-2.9.8-rhrn72yfyzknktaybz3xbvq5ud5xe7nt/include/libxml2;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/openmpi-4.0.3-pknvyjqy5moogwfifjfnmducss4tk6fr/include;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/zlib-1.2.11-kdh3v4z6c36m67oovrg55ba53abj3wev/include;/mnt/lustre/localsoft/mkl/compilers_and_libraries_2020.0.166/linux/mkl/include;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/openssl-1.0.2o-xwk65x3bzwz5rpv2ojs6hcptrqgcfbto/include;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/fftw-3.3.8-ymonpurwnnwrjsw3zw4kh5d7brp2sbn6/include;/mnt/lustre/localsoft/netcdf-v4.6.1-intel19/include;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/openmpi-3.1.3-xnwkw7bcmpsrdymg77auoyudr2pr3mu5/include;/mnt/lustre/localsoft/intel/intel19/Fortran/compilers_and_libraries_2019.5.281/linux/mkl/include;/mnt/lustre/localsoft/intel/intel19/C/compilers_and_libraries_2019.5.281/linux/ipp/include;/mnt/lustre/localsoft/intel/intel19/C/compilers_and_libraries_2019.5.281/linux/mkl/include;/mnt/lustre/localsoft/intel/intel19/C/compilers_and_libraries_2019.5.281/linux/pstl/include;/mnt/lustre/localsoft/intel/intel19/C/compilers_and_libraries_2019.5.281/linux/tbb/include;/mnt/lustre/localsoft/intel/intel19/C/compilers_and_libraries_2019.5.281/linux/daal/include;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/gcc-8.2.0/freetype-2.9.1-3cmv5kwusp6ze4bweuopbyhdybogjka2/include/freetype2;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libxml2-2.9.8-bqpperpne5khnokjv4jgy4c27xmefqpp/include/libxml2;/usr/lib/gcc/x86_64-redhat-linux/4.8.5/include;/usr/local/include;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;gcc_s;gcc;quadmath;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/netlib-lapack-3.8.0-fcypwwwujl67z2sdytr3dpok3auhlgtu/lib64;/usr/lib/gcc/x86_64-redhat-linux/4.8.5;/usr/lib64;/lib64;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/openmpi-4.0.3-pknvyjqy5moogwfifjfnmducss4tk6fr/lib;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/zlib-1.2.11-kdh3v4z6c36m67oovrg55ba53abj3wev/lib;/mnt/lustre/localsoft/mkl/compilers_and_libraries_2020.0.166/linux/compiler/lib/intel64_lin;/mnt/lustre/localsoft/mkl/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64_lin;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/openssl-1.0.2o-xwk65x3bzwz5rpv2ojs6hcptrqgcfbto/lib;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/fftw-3.3.8-ymonpurwnnwrjsw3zw4kh5d7brp2sbn6/lib;/mnt/lustre/localsoft/netcdf-v4.6.1-intel19/lib;/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/openmpi-3.1.3-xnwkw7bcmpsrdymg77auoyudr2pr3mu5/lib;/mnt/lustre/localsoft/intel/intel19/Fortran/compilers_and_libraries_2019.5.281/linux/mpi/intel64/libfabric/lib;/mnt/lustre/localsoft/intel/intel19/Fortran/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64_lin;/mnt/lustre/localsoft/intel/intel19/Fortran/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin;/mnt/lustre/localsoft/intel/intel19/C/compilers_and_libraries_2019.5.281/linux/mpi/intel64/libfabric/lib;/mnt/lustre/localsoft/intel/intel19/C/compilers_and_libraries_2019.5.281/linux/ipp/lib/intel64;/mnt/lustre/localsoft/intel/intel19/C/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64_lin;/mnt/lustre/localsoft/intel/intel19/C/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin;/mnt/lustre/localsoft/intel/intel19/C/compilers_and_libraries_2019.5.281/linux/tbb/lib/intel64/gcc4.7;/mnt/lustre/localsoft/intel/intel19/C/compilers_and_libraries_2019.5.281/linux/daal/lib/intel64_lin;/mnt/lustre/localsoft/intel/intel19/C/compilers_and_libraries_2019.5.281/linux/tbb/lib/intel64_lin/gcc4.4;/usr/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
