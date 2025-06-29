!-----------------------------------------------------------------------
!     Program:       FIELDLINES
!     Authors:       S. Lazerson
!     Date:          02/21/2012
!     Description:   The FIELDLINES code is a versatile fieldlines
!                    tracing package for MHD Equilibria.  In general,
!                    the code follows fieldlines on an R,phi,Z grid
!                    similar to an mgrid file.  It's main function is
!                    to produce Poincare plots but it can also be used
!                    to visualize the 3D trajectories of fieldlines.
!     References:
!-----------------------------------------------------------------------
      PROGRAM FIELDLINES
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE fieldlines_runtime
      USE fieldlines_interface_mod
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          numargs      Number of input arguments
!          i            Index
!          arg_len      Length of input strings
!          arg1         Input file
!          args         Input arguments
!-----------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------
!     Begin Program
!-----------------------------------------------------------------------
      
      myworkid = master
      ierr_mpi = MPI_SUCCESS

      ! Setup MPI
      CALL fieldlines_init_mpi

      ! Nullify pointers
      CALL fieldlines_init_pointers

      ! Setup HDF5
      CALL fieldlines_init_hdf5

      ! Handle the command line
      CALL fieldlines_init_commandline

      ! Output the header information
      CALL fieldlines_output_header

      ! Initialize the Calculation
      CALL fieldlines_init

      ! Handle different run type
      IF (lemc3 .or. lbfield_only .or. lafield_only) nruntype=runtype_norun
      SELECT CASE(nruntype)
         CASE(runtype_old)
            CALL fieldlines_follow
         CASE(runtype_full)
            IF (lverb) WRITE(6,'(A)') '===========ROUGH GRID=========='
            IF (lverb) WRITE(6,'(A,F8.5,A,F8.5,A)') '   EDGE_INNER[R,Z]   = [',&
                     r_start(1),',',z_start(1),']'
            IF (lverb) WRITE(6,'(A,F8.5,A,F8.5,A)') '   EDGE_OUTER[R,Z]   = [',&
                     MAXVAL(r_start,MASK = r_start > 0),',',MAXVAL(z_start,MASK = r_start > 0),']'
#if defined(MPI_OPT)
            CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init:rough_follow',ierr_mpi)
#endif
            CALL fieldlines_follow  ! This call on field grid
#if defined(MPI_OPT)
            CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init:follow',ierr_mpi)
#endif
            CALL fieldlines_init_subgrid
            CALL fieldlines_follow  ! This call on subgrid grid
            !CALL fieldlines_periodic_orbits  ! This call on subgrid grid
            !IF (lverb) CALL fieldlines_calc_surface_fit(25)
         CASE(runtype_backflow)
            IF (lverb) WRITE(6,'(A)') '===========Wall Hits=========='
            CALL fieldlines_follow
            CALL fieldlines_init_backflow
         CASE(runtype_gridgen)
            IF (lverb) WRITE(6,'(A)') '===========Grid Generation=========='
            CALL fieldlines_gridgen
         CASE(runtype_norun)
      END SELECT

      ! Output Date
      CALL fieldlines_write

      ! Clean up
      CALL fieldlines_cleanup
     
!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM FIELDLINES
