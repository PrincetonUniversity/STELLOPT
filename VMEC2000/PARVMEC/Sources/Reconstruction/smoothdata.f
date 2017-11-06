      SUBROUTINE smoothdata(nwout)
      USE vmec_main
      USE vsvd
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nwout
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: ndata_elems = 11
c-----------------------------------------------
c   L o c a l   V a r i a b l e s
c-----------------------------------------------
      INTEGER :: npts, i, ndata1
      REAL(rprec), DIMENSION(2*ns-1) :: hmid,ymid,y2mid,
     1  wmid,tenmid,requal
      REAL(rprec), DIMENSION(:), POINTER :: datap

c-----------------------------------------------
!
!       spline output data onto equally-spaced r-mesh for plotting
!
      npts = 2*ns - 1
      IF (npts .le. 1) RETURN

      DO i = 1, npts
         IF (i .le. ns) curint(i) = twopi/mu0*buco(ns + 1 - i)
         IF (i .gt. ns) curint(i) = twopi/mu0*buco(i - ns + 1)
      END DO

      DO i = 1, npts
         wmid(i) = 1.0
         tenmid(i) = 0.1
         requal(i) = rmid(1) + ((i - 1)*(rmid(npts)-rmid(1)))/
     1      (npts - 1)
      END DO

      DO ndata1 = 1, ndata_elems
         SELECT CASE (ndata1)
         CASE (1)
            datap => datamse
         CASE (2)
            datap => qmid
         CASE (3)
            datap => shear
         CASE (4)
            datap => presmid
         CASE (5)
            datap => alfa
         CASE (6)
            datap => curmid
         CASE (7)
            datap => curint
         CASE (8)
            datap => psimid
         CASE (9)
            datap => ageo
         CASE (10)
            datap => volpsi
         CASE (11)
            datap => phimid
         END SELECT
         CALL setspline (rmid, wmid, datap, hmid, ymid, y2mid, tenmid,
     1      tenmid(1), npts, natur)
         DO i = 1, npts
            CALL splint (rmid, ymid, y2mid, npts, requal(i), datap(i))
         END DO
      END DO
!
!     WRITE out splined data
!
      WRITE (nwout, 703) (datamse(i),requal(i),qmid(i),shear(i),presmid
     1   (i),alfa(i),curmid(i),i=1,npts)
      WRITE (nwout, 703) (rsort(i),ATAN(datastark(isortr(i)))/dcon,ABS(
     1   qmeas(i)),i=1,imse2 - 1)
      WRITE (nwout, 703) (rthom(i),datathom(i),i=1,itse)
  703 FORMAT(5e20.13)

      IF (lmac) THEN
        WRITE (nmac, 705)
  705 FORMAT(//,' FOLLOWING DATA EQUALLY SPACED IN R-MIDPLANE'//,
     1   '        RMID       J-PHI       SHEAR        QMID',
     2   '   MSE-PITCH     PRESMID         PSI        AMID',
     3   '      VOLUME         PHI',/,
     4   '         [M]    [A/M**2]                        ',
     5   '       [Deg]        [Pa]        [Wb]         [M]',
     6   '      [M**3]        [Wb]',/)
        WRITE (nmac, 707) (requal(i),curmid(i),shear(i),qmid(i),
     1    datamse(i),presmid(i),psimid(i),ageo(i),volpsi(i),
     2    phimid(i),i=1,npts)
        WRITE (nmac, 709) phimid(2*ns-1), psimid(2*ns-1)
      END IF
  707 FORMAT(1p,10e12.3)
  709 FORMAT(/,' phi-edge =',t16,1p,e12.3,t40,'psi-edge =',t56,e12.3)

      END SUBROUTINE smoothdata
