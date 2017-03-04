      SUBROUTINE getrz(rmc,rms,zmc,zms,r0c,z0c,rhoc,rhos,m,mrho_in)
      USE Vname1
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER m, mrho_in
      REAL(rprec), DIMENSION(0:mrho-1) :: rhoc, rhos
      REAL(rprec) :: rmc, rms, zmc, zms, r0c, z0c
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mrho1
C-----------------------------------------------

      mrho1 = mrho_in - 1
      rhos(0) = zero            !Enforce constraint on toroidal angle

      IF (m .eq. 0) THEN
         rmc = r0c
         zmc = z0c
         rms = zero
         zms = zero
      ELSE IF (m .lt. mrho1) THEN
         rmc = (t1m(m)*rhoc(m-1) + t2m(m)*rhoc(m+1))
         zms = (t1m(m)*rhoc(m-1) - t2m(m)*rhoc(m+1))
         rms = (t1m(m)*rhos(m-1) + t2m(m)*rhos(m+1))
         zmc =-(t1m(m)*rhos(m-1) - t2m(m)*rhos(m+1))
      ELSE                      !Can change highest m constraints here...
         rmc = t1m(m)*rhoc(m-1) * HB_Parameter
         zms = rmc              * HB_Parameter
         rms = t1m(m)*rhos(m-1) * HB_Parameter
         zmc =-rms              * HB_Parameter
      ENDIF

      END SUBROUTINE getrz
