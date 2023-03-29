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
    WRITE(6,*) 'setup mpi'
    CALL thrift_init_mpi
    

    ! Setup HDF5
    WRITE(6,*) 'setup hdf5'
    CALL thrift_init_hdf5

    ! Initialize constants
    WRITE(6,*) 'setup constants'
    CALL thrift_init_constants

    ! Handle the command line
    WRITE(6,*) 'commandline'
    CALL thrift_init_commandline

    ! Output the header information
    WRITE(6,*) 'header'
    CALL thrift_output_header

    ! Initialize the Calculation
    WRITE(6,*) 'init'
    CALL thrift_init

    ! ----work here
    WRITE(6,*) 'evolve'
    CALL thrift_evolve

    ! ----Write output
    CALL thrift_write

    ! Clean up
    CALL thrift_cleanup

    !-----------------------------------------------------------------------
    !     End Program
    !-----------------------------------------------------------------------
END PROGRAM THRIFT
