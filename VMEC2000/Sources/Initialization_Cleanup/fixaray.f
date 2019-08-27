      SUBROUTINE fixaray
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: jmin2, mscale, nscale, 
     &                       mnyq, nnyq, signgs
#ifdef _HBANGLE
      USE angle_constraints, ONLY: init_multipliers
#endif
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(dp), PARAMETER :: two=2, pexp=4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER  :: i, m, j, n, mn, mn1, nmin0, istat1, istat2
      INTEGER  :: mnyq0, nnyq0
      REAL(dp) :: argi, arg, argj, dnorm, tfixon, tfixoff
C-----------------------------------------------
!
!     INDEX OF LOCAL VARIABLES
!
!     mscale   array for norming theta-trig functions (internal use only)
!              so that the discrete SUM[cos(mu)*cos(m'u)] = .5 delta(m,m')
!     nscale   array for norming zeta -trig functions (internal use only)

!
!    COMPUTE TRIGONOMETRIC FUNCTION ARRAYS
!    NOTE: ARRAYS ALLOCATED HERE ARE GLOBAL AND ARE DEALLOCATED IN FILEOUT
!    NOTE: NEED 2 X NYQUIST FOR FAST HESSIAN CALCULATIONS 
!
      mnyq0  = ntheta1/2
      nnyq0  = nzeta/2
 
      mnyq = MAX(0, 2*mnyq0, 2*mpol1)
      nnyq = MAX(0, 2*nnyq0, 2*ntor)

      mnmax_nyq = nnyq/2 + 1 + mnyq*(nnyq + 1)/2

      ALLOCATE(cosmu(ntheta3,0:mnyq),  sinmu(ntheta3,0:mnyq),
     1         cosmum(ntheta3,0:mnyq), sinmum(ntheta3,0:mnyq),
     2         cosmui(ntheta3,0:mnyq), cosmumi(ntheta3,0:mnyq),
     2         cosmui3(ntheta3,0:mnyq),cosmumi3(ntheta3,0:mnyq),
     3         sinmui(ntheta3,0:mnyq), sinmumi(ntheta3,0:mnyq),
     4         cosnv(nzeta,0:nnyq),    sinnv(nzeta,0:nnyq),
     5         cosnvn(nzeta,0:nnyq),   sinnvn(nzeta,0:nnyq),
     6         cos01(nznt), sin01(nznt), stat=istat1 )
      ALLOCATE(xm(mnmax), xn(mnmax), ixm(mnsize), jmin3(0:mnsize-1),
     1         xm_nyq(mnmax_nyq), xn_nyq(mnmax_nyq),
     2         mscale(0:mnyq), nscale(0:nnyq), stat=istat2)

      IF (istat1.ne.0) THEN
         STOP 'allocation error in fixaray: istat1'
      END IF
      IF (istat2.ne.0) THEN
         STOP 'allocation error in fixaray: istat2'
      END IF

      dnorm = one/(nzeta*(ntheta2 - 1))
      IF (lasym) THEN
         dnorm = one/(nzeta*ntheta3)     !Fix, SPH012314
      END IF

      mscale(0) = 1;  nscale(0) = 1                      
!     mscale(0) = osqrt2;  nscale(0) = osqrt2    !versions < 6.9, incorrectly used osqrt2

      mscale(1:mnyq) = mscale(0)/osqrt2
      nscale(1:nnyq) = nscale(0)/osqrt2
      r0scale = mscale(0)*nscale(0)

!
!     GENERALLY, ONLY NEED THIS FROM 1, ntheta2 EXCEPT IN GETBRHO ROUTINE
!
      DO i = 1, ntheta3
         argi = twopi*(i - 1)/ntheta1
         DO m = 0, mnyq
            arg = argi*m
            cosmu(i,m) = COS(arg)*mscale(m)
            sinmu(i,m) = SIN(arg)*mscale(m)
            cosmui(i,m) = dnorm*cosmu(i,m)
            cosmui3(i,m) = cosmui(i,m)          !Use this if integration over FULL 1,ntheta3 interval 
            sinmui(i,m) = dnorm*sinmu(i,m)
            IF (i.EQ.1 .OR. i.EQ.ntheta2) THEN
               cosmui(i,m) = cosmui(i,m)/2
            END IF
            IF (ntheta2 .EQ. ntheta3) THEN
               cosmui3(i,m) = cosmui(i,m)
            END IF
            cosmum(i,m) = cosmu(i,m)*(m)
            sinmum(i,m) = -sinmu(i,m)*(m)
            cosmumi(i,m) = cosmui(i,m)*(m)
            cosmumi3(i,m) = cosmui3(i,m)*m
            sinmumi(i,m) = -sinmui(i,m)*(m)
         END DO
      END DO

      DO j = 1, nzeta
         argj = twopi*(j - 1)/nzeta
         DO n = 0, nnyq
            arg = argj*(n)
            cosnv(j,n) = COS(arg)*nscale(n)
            sinnv(j,n) = SIN(arg)*nscale(n)
            cosnvn(j,n) = cosnv(j,n)*(n*nfp)
            sinnvn(j,n) = -sinnv(j,n)*(n*nfp)
         END DO
      END DO

!
!     R,Z,L / s**(m/2) ARE LINEAR NEAR ORIGIN
!
      mn = 0
      mn1 = 0
      DO m = 0, mpol1
         xmpq(m,1) = m*(m - 1)
         xmpq(m,2) = m**pexp
         xmpq(m,3) = m**(pexp + 1)
         DO n = 0, ntor
            jmin3(mn) = jmin2(m)
            mn = mn + 1
            ixm(mn) = m
         END DO
         nmin0 = -ntor
         IF (m .eq. 0) nmin0 = 0
         DO n = nmin0, ntor
            mn1 = mn1 + 1
            xm(mn1) = m
            xn(mn1) = n*nfp
         END DO
      END DO

      IF (mn1 .ne. mnmax) THEN
         STOP 'mn1 != mnmax'
      END IF

!
!     COMPUTE NYQUIST-SIZED ARRAYS FOR OUTPUT.
!     RESTORE m,n Nyquist TO 1 X ... (USED IN WROUT, JXBFORCE)
!      mnyq = mnyq0;  nnyq = nnyq0
      mnyq = mnyq/2
      nnyq = nnyq/2

      mn1 = 0
      DO m = 0, mnyq
         nmin0 = -nnyq
         IF (m .eq. 0) nmin0 = 0
         DO n = nmin0, nnyq
            mn1 = mn1 + 1
            xm_nyq(mn1) = m
            xn_nyq(mn1) = n*nfp
         END DO
      END DO

      IF (mn1 .ne. mnmax_nyq) THEN
         STOP 'mn1 != mnmax_nyq'
      END IF

      mn = 0
      m = 1
      DO i = 1, ntheta3
         argi = twopi*(i - 1)/ntheta1
         DO j = 1, nzeta
            mn = mn + 1
            cos01(mn) =  m*COS(m*argi)*mscale(m)
            sin01(mn) = -m*SIN(m*argi)*mscale(m)
         END DO
      END DO

      faccon(0) = zero
      faccon(mpol1) = zero
      faccon(1:mpol1 - 1) = -0.25_dp*signgs/xmpq(2:mpol1,1)**2

#ifdef _HBANGLE
      CALL init_multipliers
#endif

      IF (lrecon) THEN
         CALL getgreen
      END IF

      END SUBROUTINE fixaray
