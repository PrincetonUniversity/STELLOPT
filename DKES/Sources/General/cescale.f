      SUBROUTINE cescale(srces)

c  This subroutine performs cmul and efield scaling, and calculates the
c  dominant - diagonal scaling arrays.

C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vnamecl2
      USE dkes_input, ONLY: lalpha, psip
      USE dkes_realspace
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: l0 = 1, l1 = 2, l2 = 3, l3 = 4,
     1   epar = 2, pgrad = 1, plus = 1, minus = 2
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(mpnt,lsource,2,2), INTENT(out) :: srces
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn
      REAL(rprec) :: cmul2, eoc, e2oc, SQRT10, bmod
C-----------------------------------------------
      SQRT10 = SQRT(10._dp)

      IF (psip .ne. zero) THEN
         weov = efield1/psip                               !!EXB drift frequency (normed to 1/v)
      ELSE
         weov = efield1/EPSILON(psip)
      END IF

      cmul2 = cmul1*cmul1
      wcyclo = 9.58e7_dp * b00                          !!b00 in [Tesla]
      vthermi = 9.79e3_dp * SQRT(2._dp) * SQRT(1.e3_dp) !!vi for Ti=1keV, [m/s]

c  Spitzer FUNCTION contribution to conductivity
c  Note that fspitzer = qB/nu/SQRT(bsqav), S3 = qB*jacobian/SQRT(bsqav), so that
c  g33s ~ INT(q**2)/nu * vp.
      g33s = one/(cmul1*rt3o2**2)                       !!1./rt3o2**2 = 1/(2/3) from pitch integral of q**2

      eoc  = -efield1/cmul1
      e2oc = -efield1*eoc

      cols = cmul1*cols0(1:lalpha)
      al1 = al01/cmul1
      al2 = al02/cmul1
      al3 = e2oc*al03
      al4 = al04/cmul1
      bl1 = eoc*bl01
      bl2 = eoc*bl02
      bl3 = eoc*bl03
      bl4 = eoc*bl04
      cl1 = cl01/cmul1
      cl2 = cl02/cmul1
      cl3 = cl03/cmul1
      cl4 = cl04/cmul1

      bmod = SQRT(bsqav)
      s1cs1 = s1cs10/cols(l2)

c     e-field dependent particle conservation matrix elements
      diagle = 0
      DO mn = 1,mpnt
         diagle(mn,mn,1) = efield1*exbgrad(mn)
         diagle(mn,mn,2) =-efield1*exbgrad(mn)
      END DO

c     itype = pgrad  (density gradient sources)
      srces(:,l1,pgrad,minus) = omgl(l2)*auxs1(:,1)/cols(l2)
      srces(:,l2,pgrad,minus) = auxs1(:,2)*eoc/cols0(l2)
      srces(:,l3,pgrad,minus) = omgl(l3)*auxs1(:,3)/cols(l2)


c     itype = epar (parallel E-field sources)
      srces(:,l1,epar,plus) = auxs3p(:,1)*eoc/bmod                       !!l=1 component of s3
      srces(:,l2,epar,plus) = 3*omgl(l2)*auxs3p(:,2)/cmul1/bmod          !!l=2 component of s3
      srces(:,l1,epar,minus) = auxs3m(:)/bmod

      END SUBROUTINE cescale
