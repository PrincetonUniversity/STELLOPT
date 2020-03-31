!*******************************************************************************
!  File packafile.f
!      Module for use in program v3post
!      Packs the responses into a single binary file as  NETCDF files
!         for the coil responses and the plasma responses are read.
! packed file format: (intgerer*4 and real(rprec) logicals and strings)
!	       n_diagn_c n_field_cg
!	       cl_response(1:n_diagn_c,1:n_field_cg)
!	       ir jz kp kp_store, 
!	       rmin,rmax,zmin,zmax
!	       n_field_periods
!	       lstell_sym
!	        REPEAT n_diagn_c TIMES:
!	        diag# shortname
!	        component#, array(ir,jz,kp_store)
! THIS PRODUCES A FILES SIZE OF ABOUT 7.9 Mb / DIAGNOSTIC
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!    This module uses the following modules:
!       stel_kinds
!       stel_constants
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  Initial Coding - Ed Lazarus 03.31.2005
!-------------------------------------------------------------------------------
!   USAGE : Called within read_response if pack=T in V3POST namelist
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!   COMMENTS I A response to the failure to successfully read prfun files when
!		many parallel proccesses are active
!-------------------------------------------------------------------------------
!
!
!*******************************************************************************
!  MODULE read_response
!    
! SECTION I.	 VARIABLE DECLARATIONS
! SECTION II.	 INTERFACE BLOCKS
! SECTION III.	 COIL PACKING
! SECTION IV.	 PLASMA PACKING
!*******************************************************************************

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

      MODULE packafile

! work_arrays avoids premature size definition of arrays
      USE stel_kinds
      USE stel_constants
      
       
!  Other variables
      INTEGER :: npack = 0, nhere=1
      INTEGER :: one1 = 1, two2 = 2, three3 = 3
      LOGICAL :: pack_in_read = .false.
      ChARACTER(len=80) :: packafilename
!*******************************************************************************
! SECTION II. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------

      CONTAINS

!*******************************************************************************
! SECTION III. COIL RESPONSE FUNCTION READ
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE open_packafile (filename)
      USE safe_open_mod
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(len=*) :: filename
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: istat


      CALL safe_open(npack, istat, filename, 
     1	'unknown', 'unformatted')
      IF (istat .ne. 0) STOP 
     1	'safe_open iostat != 0 in open_packafile'
! closed on exit from V3POST

      END SUBROUTINE open_packafile

      SUBROUTINE read_packafile (filename)
      USE safe_open_mod
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(len=*) :: filename
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: istat


      CALL safe_open(npack, istat, filename,
     1  'old', 'unformatted')
      IF (istat .ne. 0) STOP
     1  'safe_open iostat != 0 in read_packafile'
! closed on exit from V3POST

      END SUBROUTINE read_packafile


      END MODULE packafile

