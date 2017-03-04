!*******************************************************************************
!>  @file mpi_inc.f
!>  @brief Contains module @ref mpi_inc.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Umbrella module avoid multiple inlcudes of the mpif.h header.
!*******************************************************************************

      MODULE mpi_inc
      USE mpi_params, ONLY: MPI_COMM_PARVMEC
      IMPLICIT NONE
      !INTEGER, PARAMETER :: MPI_COMM_PARVMEC=452
#if defined(MPI_OPT)
      INCLUDE 'mpif.h'
#endif

      END MODULE
