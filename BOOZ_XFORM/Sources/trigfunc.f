      SUBROUTINE trigfunc (uang, vang, cosm, sinm, cosn, sinn, 
     1                     mpol, ntor, nznt)
      USE stel_kinds
      USE booz_params, ONLY: nfp
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: mpol, ntor, nznt
      REAL(rprec), DIMENSION(nznt), INTENT(in) :: uang, vang
      REAL(rprec), DIMENSION(nznt,0:mpol), INTENT(out) :: cosm, sinm
      REAL(rprec), DIMENSION(nznt,0:ntor), INTENT(out) :: cosn, sinn
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, n
C-----------------------------------------------

      cosm(:,0) = 1
      sinm(:,0) = 0
      cosm(:,1) = COS(uang)
      sinm(:,1) = SIN(uang)

      DO m = 2,mpol
         cosm(:,m) = cosm(:,m-1)*cosm(:,1)
     1             - sinm(:,m-1)*sinm(:,1)
         sinm(:,m) = sinm(:,m-1)*cosm(:,1)
     1             + cosm(:,m-1)*sinm(:,1)
      END DO

      cosn(:,0) = 1
      sinn(:,0) = 0
      IF (ntor .ge. 1) THEN
         cosn(:,1) = COS(vang*nfp)
         sinn(:,1) = SIN(vang*nfp)
      END IF

      DO n = 2,ntor
         cosn(:,n) = cosn(:,n-1)*cosn(:,1)
     1             - sinn(:,n-1)*sinn(:,1)
         sinn(:,n) = sinn(:,n-1)*cosn(:,1)
     1             + cosn(:,n-1)*sinn(:,1)
      END DO

      END SUBROUTINE trigfunc