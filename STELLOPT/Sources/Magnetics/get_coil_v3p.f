      SUBROUTINE get_coil_v3p (coilfilename)
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          get_coil_v3p computes external coil contribution        **
!**                                                                  **
      USE stel_kinds
      USE stel_constants
      USE v3post_rfun
!DEC$ IF DEFINED (MPI_OPT)
      USE read_response
!DEC$ ELSE
      USE read_response_nompi
!DEC$ ENDIF
      USE safe_open_mod
      USE read_wout_mod, ONLY: mgrid_mode, extcur, ns
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*) , INTENT(in) :: coilfilename

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER(iprec) :: i, j=100, istat

!  Declaration of derived type, for coil response function
      TYPE(clresfun) crf
!  Local declaration of coil response pointer, sizes
!  JDH 06.28.2003. Note that n_diagn_c is in module v3post_rfun,
!     and n_diagn_c needs to be conveyed to get_plasma_v3p
      REAL(rprec), DIMENSION(:,:), POINTER :: cl_response
      INTEGER(iprec) :: n_field_cg

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!----------------------------------------------------------------------
!-- Get external coil response functions                             --
!-- Read external coil currents from wout file                       --
!----------------------------------------------------------------------
      IF (mgrid_mode == 'N') THEN
         WRITE (6, '(a,/,a)')' Old-style extcur array in wout file',
     1         ' Assuming "Scaled" mode for extcur data'
         mgrid_mode = 'S'
      ENDIF
!DEC$ IF DEFINED (MPI_OPT)
      CALL cdf_crfun_read(TRIM(coilfilename),crf,
     &                    istat, mgrid_mode)
!DEC$ ELSE
      CALL cdf_crfun_read_nompi(TRIM(coilfilename),crf,     
     &    istat, mgrid_mode)
!DEC$ ENDIF

      IF (istat .ne. 0) THEN
         WRITE (6, *)' In get_coil_v3p, there was an error reading',
     1           ' cdf file: ', TRIM(coilfilename)
         STOP
      ENDIF

!  Point cl_response to component from derived type
      cl_response => crf%rdiag_coilg

!  Define local variable (n_field_cg) and module v3post_rfun variable
!    (n_diagn_c) for array sizes
      n_diagn_c = SIZE(cl_response,1)
      n_field_cg = SIZE(cl_response,2)
      IF (n_diagn_c .ne. crf%n_diagn_c) STOP 'get_coil_v3p #1'
      IF (n_field_cg .ne. crf%n_field_cg) STOP 'get_coil_v3p #2'

!----------------------------------------------------------------------
      IF (.not. ALLOCATED(extcur)) 
     1   STOP 'In get_coil_v3p, extcur not allocated'
      IF (SIZE(extcur) .ne. n_field_cg) WRITE(6,*)
     1    'Warning: external coil/wout nextcur mismatch'

!----------------------------------------------------------------------
!-- Allocate flux loop and magnetic probe signals                    --
!----------------------------------------------------------------------
      IF (ALLOCATED(signal_diag)) DEALLOCATE (signal_diag)
      ALLOCATE (signal_diag(n_diagn_c))
!----------------------------------------------------------------------
!-- Compute magnetic responses due to external coils                 --
!----------------------------------------------------------------------
      signal_diag(1:n_diagn_c)%cal = zero
      DO i = 1, MIN(n_field_cg, SIZE(extcur))
        signal_diag(1:n_diagn_c)%cal =
     &            signal_diag(1:n_diagn_c)%cal
     &          + cl_response(1:n_diagn_c,i)*extcur(i)
      ENDDO
      signal_diag(1:n_diagn_c)%cext = signal_diag(1:n_diagn_c)%cal

      RETURN

!----------------------------------------------------------------------
!----------------------------------------------------------------------

      ENTRY no_coil (coilfilename)
! set n_diagn_c
      CALL safe_open(j, istat, TRIM(coilfilename), 'old', 'formatted')
      IF (istat .ne. 0) STOP 'safe_open'
      DO WHILE (istat .eq. 0)
        READ (j, *, iostat=istat, end=10005) n_diagn_c
10005 CONTINUE
      ENDDO
      IF (ALLOCATED(signal_diag)) DEALLOCATE (signal_diag)
      ALLOCATE (signal_diag(n_diagn_c))
      signal_diag(1:n_diagn_c)%cal = zero
      signal_diag(1:n_diagn_c)%cext = zero

      CLOSE (j)
      END SUBROUTINE get_coil_v3p
