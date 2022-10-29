!-----------------------------------------------------------------------
!     Program:       THRIFT
!     Authors:       L. van Ham, S. Lazerson
!     Date:          11/XX/2022
!     Description:   The THRIFT code performs current profile evolution
!                    calculations for three dimensional equilbria.
!     References:
!-----------------------------------------------------------------------
PROGRAM THRIFT
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE thrift_runtime
    USE thrift_interface_mod

    !-----------------------------------------------------------------------
    !     Local Variables
    !-----------------------------------------------------------------------
    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !     Begin Program
    !-----------------------------------------------------------------------

    ! Setup MPI
    CALL thrift_init_mpi

    ! Setup HDF5
    CALL thrift_init_hdf5

    ! Initialize constanst
    CALL thrift_init_constants

    ! Handle the command line
    CALL thrift_init_commandline

    ! Output the header information
    CALL thrift_output_header

    ! Initialize the Calculation
    !CALL thrift_init

    ! ----work here

    ! Clean up
    CALL thrift_cleanup

    !-----------------------------------------------------------------------
    !     End Program
    !-----------------------------------------------------------------------
END PROGRAM THRIFT
