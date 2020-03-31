      SUBROUTINE transpmn(pmns, bsubtmnc, bsubzmnc, pmnc, bsubtmns, 
     1              bsubzmns, xm, xn, gpsi, ipsi, mnmax, js, lasym)
      USE stel_kinds
      USE booz_params, ONLY: ns
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: mnmax, js
      REAL(rprec), DIMENSION(mnmax) :: pmns, bsubtmnc, bsubzmnc
      REAL(rprec), DIMENSION(mnmax) :: pmnc, bsubtmns, bsubzmns
      REAL(rprec), DIMENSION(mnmax) :: xm, xn
      REAL(rprec), DIMENSION(ns)    :: gpsi, ipsi
      LOGICAL, INTENT(in) :: lasym
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn
      REAL(rprec) :: rad_cur
C-----------------------------------------------
!
!     COMPUTE THE PART OF Pmns,c WHICH IS INDEPENDENT OF Lmns,c (the COVARIANT source terms
!     in Eq.10). THE NET P IN REAL SPACE IS GIVEN AS FOLLOWS:
!
!     P(FINAL) = {SUM(m,n)[pmns(m,n)*SIN(arg) + pmnc(m,n)*COS(arg)] - Ipsi*Lambda} / D
!
!     WHERE arg = m*theta - n*zeta, D = gpsi + iota*Ipsi
!
      DO mn = 1,mnmax
        IF (NINT(xm(mn)) .ne. 0) THEN
          pmns(mn) = bsubtmnc(mn)/xm(mn)
          IF (lasym) pmnc(mn) = -bsubtmns(mn)/xm(mn)
        ELSE IF (NINT(xn(mn)) .ne. 0) THEN
          pmns(mn) = -bsubzmnc(mn)/xn(mn)
          IF (lasym) pmnc(mn) =  bsubzmns(mn)/xn(mn)
        ELSE
          pmns(mn) = 0
          IF (lasym) pmnc(mn) = 0
          gpsi(js) = bsubzmnc(mn)
          Ipsi(js) = bsubtmnc(mn)
        ENDIF

!       DIAGNOSTIC: CHECK IF RADIAL CURRENT VANISHES
        CYCLE !COMMENT THIS FOR DIAGNOSTIC DUMP
        IF (xm(mn).ne.0 .or. xn(mn).ne.0) THEN
           rad_cur = xn(mn)*bsubtmnc(mn)+xm(mn)*bsubzmnc(mn)
           rad_cur = rad_cur/(ABS(Ipsi(js))+ABS(gpsi(js)))
           IF (ABS(rad_cur) .gt. 1.E-10_dp) THEN
              PRINT *,'m: ',xm(mn),' n: ',xn(mn),' Radial Current: ',
     1                rad_cur
           END IF
        END IF

      END DO

      END SUBROUTINE transpmn
