      SUBROUTINE chisq_kappa(target_in, sigma, ivar, num, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: bigno, lextcur, mpol1_opt, mnmax_opt
      USE boozer_params, ONLY: rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ivar, nopt
      INTEGER :: num
      REAL(kind=rprec), INTENT(in) :: sigma, target_in
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(kind=rprec), PARAMETER :: zero = 0, c1p5 = 1.5_dp,
     1    xtpi=3.141592654
      INTEGER, PARAMETER :: nthsrch = 50
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, jmax, k
      REAL(rprec) rmin, rmax, zmax, znew, thsrch, th0, dth
      REAL(rprec), DIMENSION(0:mpol1_opt) :: rbc, zbs
C-----------------------------------------------
      IF ( ABS(sigma) >= bigno) RETURN

      num = num + 1

      IF (nopt .gt. 0) THEN
         index_array(num) = ivar

         wegt(num) = sigma
         chisq_target(num) = target_in

!    find the n=0 terms
         rbc = 0
         zbs = 0

         DO j=1, mnmax_opt
            IF( NINT(xn_bdy(j)) == 0 ) THEN
               i = NINT(xm_bdy(j))
               rbc(i) = rmnc_bdy(j)
               zbs(i) = zmns_bdy(j)
            ENDIF
         ENDDO

         rmax = SUM(rbc(0:mpol1_opt))
         rmin = SUM(rbc(0:mpol1_opt:2)) - SUM(rbc(1:mpol1_opt:2))
         zmax = 0

         th0 = 0
         dth = xtpi/REAL(nthsrch+1,kind=rprec)
         jmax = 0

         DO k=1, 2
           DO j=1,nthsrch
             thsrch=th0 + j*dth

             znew = SUM(zbs(0:mpol1_opt)*
     1             COS(thsrch*(/ (i, i=0, mpol1_opt) /) ))

             IF( znew > zmax ) THEN
                jmax = j - 1
                zmax = znew
             ENDIF
           ENDDO

           th0 = th0 + jmax*dth
           dth = (2*dth)/(nthsrch+1)
         ENDDO

         chisq_match(num) = 2*zmax/(rmax - rmin)

      ELSE
         IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
      ENDIF

      END SUBROUTINE chisq_kappa
