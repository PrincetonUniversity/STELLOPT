      SUBROUTINE interactive_input()
!  interactive_input
!    For task MGRID, input parameters when prompted. This is for
!    backwards compatibility with the earlier version of makegrid.

      USE write_mgrid, only: mgrid_ext, mgrid_mode, lstell_sym,                & 
     &   rmin, rmax, zmin, zmax, kp, ir, jz,                                   & 
     &   mgrid_file

      USE makegrid_global, only: task

      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: numargs, i
      CHARACTER(LEN=100) :: arg1
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      
      WRITE (6, 220, advance='no')                                             & 
     &  ' Enter extension of "coils" file     : '
      READ (*, *) mgrid_ext

      WRITE (6, '(a,/,a)', advance='no')                                       &  
     &  ' Scale (S) bfield to unit current/turn OR',                           & 
     &  ' use raw (R) currents from coils file: '
      READ (*, *) mgrid_file
      IF (mgrid_file(1:1) == 'R' .or. mgrid_file(1:1) == 'r')
     &   mgrid_mode = 'R'

      WRITE (6, 220, advance='no')                                             & 
     &  ' Assume stellarator symmetry (Y/N)?  : '
      READ (*, *) mgrid_file
      IF (mgrid_file(1:1) == 'Y' .or. mgrid_file(1:1) == 'y')
     &   lstell_sym = .true.

      WRITE (6, 220, advance='no')                                             & 
     &  ' Enter rmin (min radial grid dimension)  : '
      READ (*, *) rmin

      WRITE (6, 220, advance='no')                                             &  
     &  ' Enter rmax (max radial grid dimension)  : '
      READ (*, *) rmax

      WRITE (6, 220, advance='no')                                             & 
     &  ' Enter zmin (min vertical grid dimension): '
      READ (*, *) zmin

      WRITE (6, 220, advance='no')                                             &  
     &  ' Enter zmax (max vertical grid dimension): '
      READ (*, *) zmax

      WRITE (6, 220, advance='no')                                             & 
     &  ' Enter number of toroidal planes/period  : '
      READ (*, *) kp

      WRITE (6, 220, advance='no')                                             & 
     &  ' Enter number of r (radial) mesh points  : '
      READ (*, *) ir

      WRITE (6, 220, advance='no')                                             & 
     &  ' Enter number of z mesh points  : '
      READ (*, *) jz
      task='MGRID'
      
 220  FORMAT(a)
      
      END SUBROUTINE interactive_input