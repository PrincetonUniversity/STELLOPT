
SUBROUTINE neo_init(npsi)
! Initialization Routine
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE
  INTEGER, INTENT(out)       :: npsi
  INTEGER                    :: imn
! **********************************************************************
! Read input from data file and allocate necessary arrays
! **********************************************************************
!  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_read'
  CALL neo_read
!  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_read'
! **********************************************************************
  npsi = ns
! **********************************************************************
! Allocate and prepare necessary arrays
! **********************************************************************
!  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_prep'
  CALL neo_prep
!  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_prep'
! **********************************************************************
! Calculation of rt0 and bmref
! **********************************************************************
  rt0=0
  bmref=0
  DO imn=1,mnmax
     IF(ixm(imn).EQ.0 .AND. ixn(imn).EQ.0) THEN
        rt0 = rmnc(1,imn)
        bmref = bmnc(1,imn)

        rt0_g = rt0
        bmref_g = bmref
     ENDIF
  ENDDO
  IF(rt0.EQ.ZERO .OR. bmref.EQ.ZERO) THEN
    WRITE (w_us,*) ' NEO_INIT: Fatal problem setting rt0 or bmref'
    STOP
  ENDIF
!
  nper = nfp
! **********************************************************************
  w_u6_open = 0
! **********************************************************************
  RETURN
END SUBROUTINE neo_init
