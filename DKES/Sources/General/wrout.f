      SUBROUTINE wrout(fz1s, fz1c, fz3s, fz3c, srces)

c  This SUBROUTINE creates the file DKESF containing the magnetic field
c  and the final solution distributions.

C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vnamecl2
      USE safe_open_mod
      USE dkes_input
      USE dkes_realspace, ONLY: mvalue, nvalue, mpnt, borbi1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(mpnt,0:lalpha) :: fz1s, fz1c, fz3s, fz3c
      REAL(rprec), DIMENSION(mpnt,4,2,2) :: srces
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: idkes = 7, iobin = 66
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, n, mn, ll, l, iodkes = idkes
C-----------------------------------------------

      CALL safe_open(iodkes, mn, 'DKESF', 'unknown', 'formatted')
      IF (mn .ne. 0) STOP 'ERROR OPENING DKESF FILE'

      WRITE (iodkes, 100) chip, psip, nzperiod, cmul1, efield1, wtov,
     1   weov, mpolb, ntorb, lfout, mpnt, blabl(ibbi)
  100 FORMAT(/8x,'CHIP',8x,'PSIP',4x,'NZPERIOD',8x,'CMUL',6x,'EFIELD',4x
     1   ,'OMEGAT/V',4x,'OMEGAE/V',7x,'MPOLB',7x,'NTORB',7x,'LFOUT',2x,
     2   'BLOCK SIZE'//1p,2e12.4,i12,4e12.4,4i12,2/4x,'M',4x,'N',14x,a3,
     3   18x,'BTHETA',7x,'BZETA'/)

!     WRITE FIRST ROW

      DO mn = 1, mpnt
         m = mvalue(mn)
         n = nvalue(mn)
         IF (m.ne.0 .or. n.ne.0) CYCLE
         WRITE (iodkes, 200) m, n, borbi1(mn), btheta, bzeta
         EXIT
      END DO
  200 FORMAT(2i5,5x,1pe12.4,12x,2e12.4)
     
      DO mn = 1, mpnt
         m = mvalue(mn)
         n = nvalue(mn)
         IF (m.ne.0 .and. n.eq.0 .and. borbi1(mn).ne.zero) THEN
            WRITE (iodkes, 300) m, n, borbi1(mn)
         END IF
      END DO
  300 FORMAT(2i5,5x,1pe12.4)

      DO n = 1, ntorb
         DO mn = 1, mpnt
           IF (nvalue(mn).eq.n .and. borbi1(mn).ne.zero) THEN
               WRITE (iodkes, 300) mvalue(mn),n,borbi1(mn)
           END IF
         END DO
      END DO

      IF (ipmb .eq. 1) THEN
         fz1c(:,:lfout) = zero
         fz3c(:,:lfout) = zero
      ENDIF
      IF (ipmb .eq. 2) THEN
         fz1s(:,:lfout) = zero
         fz3s(:,:lfout) = zero
      ENDIF

      n = nzperiod
      IF (n .eq. 0) n = 1

!
!     WRITE "CONSTRAINT" PHANTOM DISTRIBUTION (l = -1, remember l is offset by 1 here)
!
      WRITE (iodkes, 400)
  400 FORMAT(/,'  Constraint Phantom Distribution',//,
     &   4x,'M',4x,'N',13x,'F11C',8x,'F11S',8x,'F13C',8x,'F13S'/)
      WRITE (iodkes, 500) (mvalue(mn),nvalue(mn),fz1s(mn,0),fz1c(mn,0),
     1   fz3s(mn,0),fz3c(mn,0),mn=1,mpnt)
  500 FORMAT(2i5,5x,1p,4e12.4)

!
!     WRITE DISTRIBUTION FUNCTION FOURIER COEFFIENTS FOR EACH LEGENDRE PITCH
!     IN NAMELIST, SET idisk:
!        idisk = 0 : DKES computes all l values.
!        idisk = 1 : DKES computes only l = 0,1,2,3,4,5 (default)
!
!     THESE ARE COEFFICIENTS OF ORTHONORMALIZED LEGENDRE FUNCTIONS and FOURIER HARMONICS
!
      OPEN (unit=iobin, file='dkes_bin', status='replace', 
     &      form='unformatted')
      WRITE (iodkes, 600)
  600 FORMAT(/, '  Fourier Harmonics of Distribution Functions',//,
     1   4x,'M',4x,'N',5x,'F01C',8x,'F01S',8x,'F03C',8x,
     1   'F03S',14x,'S01C',8x,'S01S',8x,'S03C',8x,'S03S'/)

      WRITE (iobin) mpnt, lfout
      l = -1
      DO ll = 1, lfout
         l = l + 1
         WRITE (iodkes, 650) l
         WRITE (iobin) l
         WRITE (iobin) (mvalue(mn), nvalue(mn), fz1c(mn,ll), 
     &               fz1s(mn,ll), fz3c(mn,ll), fz3s(mn,ll), mn=1,mpnt)
         IF (ll > 4) WRITE (iodkes, 700) (mvalue(mn),nvalue(mn),
     1      fz1c(mn,ll),fz1s(mn,ll),fz3c(mn,ll),fz3s(mn,ll),mn=1,mpnt)
         IF (ll <= 4)
     1      WRITE (iodkes, 800) (mvalue(mn),nvalue(mn),
     2      fz1c(mn,ll),fz1s(mn,ll),fz3c(mn,ll),fz3s(mn,ll),
     3      srces(mn,ll,1,2),srces(mn,ll,1,1),srces(mn,ll,2,2),
     4      srces(mn,ll,2,1),mn=1,mpnt)
      END DO

  650 FORMAT(/,2x,'PITCH L =',i5,/) 
  700 FORMAT(2i5,1p,4e12.4)
  800 FORMAT(2i5,1p,4e12.4,6x,4e12.4)

      CLOSE(iodkes)
      CLOSE(iobin)

      END SUBROUTINE wrout
