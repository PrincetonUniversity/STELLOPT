!-----------------------------------------------------------------------
!     Program:       BEAMS3D
!     Authors:       M. McMillan S. Lazerson
!     Date:          06/20/2012
!     Description:   The BEAMS3D code performs Monte-Carlo particle
!                    simulations on an R-phi-Z cylindrical grid.
!     References:
!-----------------------------------------------------------------------
PROGRAM BEAMS3D
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE beams3d_runtime
    USE beams3d_interface_mod

    !-----------------------------------------------------------------------
    !     Local Variables
    !-----------------------------------------------------------------------
    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !     Begin Program
    !-----------------------------------------------------------------------

    ! Setup MPI
    CALL beams3d_init_mpi

    ! Setup HDF5
    CALL beams3d_init_hdf5

    ! Initialize constanst
    CALL beams3d_init_constants

    ! Handle the command line
    CALL beams3d_init_commandline

    ! Output the header information
    CALL beams3d_output_header

    ! Initialize the Calculation
    CALL beams3d_init

    ! Follow Fieldlines
    CALL beams3d_follow_gc

    ! Write Ouput
    CALL beams3d_write('TRAJECTORY_PARTIAL')
    IF (lascot) THEN
        IF (lascotfl) THEN
            CALL beams3d_write_ascoth5('FIELDLINES')
        ELSE
            CALL beams3d_write_ascoth5('MARKER')
        END IF
    END IF
    IF (lascot4) CALL beams3d_write_ascoth4('MARKER')

    ! Write diagnostics stuff
    CALL beams3d_diagnostics

    !Write Fidasim Distribution function
     IF (lfidasim) THEN
        IF (lverb) THEN
            WRITE(6, '(A)') '----- WRITING FIDASIM DISTRIBUTION -----'
         END IF
        ! IF (ldepo) THEN
        !     CALL beams3d_write_fidasim('DISTRIBUTION_GC_MC') !not implemented yet
        ! ELSE
            CALL beams3d_write_fidasim('DISTRIBUTION_GC_F') !Should stay here as it alters dist5d_prof
        ! END IF
     END IF

    ! Clean up
    CALL beams3d_cleanup

    !-----------------------------------------------------------------------
    !     End Program
    !-----------------------------------------------------------------------
END PROGRAM BEAMS3D
