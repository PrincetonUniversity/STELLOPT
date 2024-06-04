!-----------------------------------------------------------------------
!     Program:       WOUT2INDATA
!     Authors:       J. Schilling
!     Date:          06/04/2024
!     Description:   The WOUT2INDATA code reads a VMEC2000 INDATA
!                    namelist and a VMEC2000 WOUT file from
!                    the current LIBSTELL library and writes out
!                    a new INDATA namelist with the plasma boundary
!                    from the WOUT file.
!     References:
!-----------------------------------------------------------------------
      PROGRAM WOUT2INDATA
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE vmec_input
      USE read_wout_mod
      USE safe_open_mod
!-----------------------------------------------------------------------
!     Local Variables
!          numargs      Number of input arguments
!          i            Index
!          arg_len      Length of input strings
!          arg1         Input file
!          args         Input arguments
!-----------------------------------------------------------------------
      IMPLICIT NONE
      integer                                      :: numargs,i,ier,&
                                                      istat, js, mn, &
                                                      fid
      integer, parameter                           :: arg_len = 256
      character*(arg_len)                          :: wout_file
      character*(arg_len)                          :: input_ext
      character*(arg_len),allocatable,dimension(:) :: args
!-----------------------------------------------------------------------
!     Begin Program
!-----------------------------------------------------------------------
      ! First Handle the input arguments


      CALL GETCARG(1, wout_file, numargs)
      wout_file = TRIM(wout_file)

      ! Read the INDATA namelist.

      ! Read the wout file.
      CALL read_wout_file(wout_file,ier)
      IF (ier /=0) THEN
         WRITE(*,*) 'Error Reading File: ',TRIM(wout_file)
         STOP
      END IF

      ! Replace the boundary geometry in the INDATA namelist
      ! with the plasma boundary coefficients from the wout file.



      ! Write the new INDATA namelist.



!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM WOUT2INDATA
