      SUBROUTINE vacuum(rmnc, rmns, zmns, zmnc, xm, xn, 
     1                  plascur, rbtor, wint, ns, ivac_skip, ivac, 
     2                  mnmax, ier_flag, lscreen)
      USE vacmod
      USE vparams, ONLY: nthreed, zero, one, mu0
      USE vmec_params, ONLY: norm_term_flag, phiedge_error_flag
      USE vmec_input, ONLY: lrfp         ! JDH Added 2013-11-25, to test for RFP
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ns, ivac_skip, ivac, mnmax, ier_flag
      REAL(rprec) :: plascur, rbtor
      REAL(rprec), DIMENSION(mnmax), INTENT(in) ::
     1   rmnc, rmns, zmns, zmnc, xm, xn
      REAL(rprec), DIMENSION(*), INTENT(in) :: wint
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, n, n1, m, i, info
      REAL(rprec), DIMENSION(:), POINTER :: potcos, potsin
      REAL(rprec), ALLOCATABLE :: bsubu(:), bsubv(:), potu(:), potv(:)
      REAL(rprec), ALLOCATABLE :: amatrix(:)
      REAL(rprec):: dn2, dm2, cosmn, sinmn, huv, hvv,
     1    det, bsupu, bsupv, bsubuvac, fac
C-----------------------------------------------
!
!     THIS ROUTINE COMPUTES .5 * B**2 ON THE VACUUM / PLASMA SURFACE
!     BASED ON THE PROGRAM BY P. MERKEL [J. Comp. Phys. 66, 83 (1986)]
!     AND MODIFIED BY W. I. VAN RIJ AND S. P. HIRSHMAN (1987)

!     THE USER MUST SUPPLY THE FILE << MGRID >> WHICH INCLUDES THE MAGNETIC
!     FIELD DATA TO BE READ BY THE SUBROUTINE BECOIL
!     THE "VACUUM.INC" FILE IS DEFINED IN VMEC.UNIX
!
!
      ier_flag = norm_term_flag

      IF (.not.ALLOCATED(potvac)) STOP 'POTVAC not ALLOCATED in VACCUM'

      ALLOCATE (amatrix(mnpd2*mnpd2), bsubu(nuv2), bsubv(nuv2),
     1    potu(nuv2), potv(nuv2), stat = i)
      IF (i .ne. 0) STOP 'Allocation error in vacuum'

      potsin => potvac(1:mnpd)
      potcos => potvac(1+mnpd:)

      ALLOCATE (bexu(nuv2), bexv(nuv2), bexn(nuv2),
     1     bexni(nuv2), r1b(nuv), rub(nuv2), rvb(nuv2),
     2     z1b(nuv), zub(nuv2), zvb(nuv2), auu(nuv2), auv(nuv2),
     3     avv(nuv2), snr(nuv2), snv(nuv2), snz(nuv2), drv(nuv2),
     4     guu_b(nuv2), guv_b(nuv2), gvv_b(nuv2), rzb2(nuv),
     5     rcosuv(nuv), rsinuv(nuv), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in vacuum'

!
!       INDEX OF LOCAL VARIABLES
!
!       rmnc,rmns,zmns,zmnc:     Surface Fourier coefficients (m,n) of R,Z
!       xm,xn:     m, n values corresponding to rc,zs array
!       bsqvac:    B**2/2 at the vacuum INTERFACE
!       plascur:   net toroidal current
!       rbtor  :   net (effective) poloidal current (loop integrated R*Btor)
!       mnmax:     number of R, Z modes in Fourier series of R,Z
!       ivac_skip: regulates whether full (=0) or incremental (>0)
!                 update of matrix elements is necessary
!
!
!       compute and store mean magnetic fields (due to
!       toroidal plasma current and EXTERNAL tf-coils)
!       note: these are fixed for a constant current iteration
!
!       bfield = rbtor*grad(zeta) + plascur*grad("theta") - grad(potential)
!
!       where "theta" is computed using Biot-Savart law for filaments
!       Here, the potential term is needed to satisfy B dot dS = 0 and has the form:
!
!       potential = SUM potsin*SIN(mu - nv) + potcos*COS(mu - nv)
!

      IF (.not. ALLOCATED(tanu)) CALL precal
      CALL surface (rmnc, rmns, zmns, zmnc, xm, xn, mnmax)
      CALL bextern (plascur, wint, ns)
!
!     Determine scalar magnetic potential POTVAC
!
      CALL scalpot (potvac, amatrix, wint, ns, ivac_skip)
      CALL solver (amatrix, potvac, mnpd2, 1, info)
      IF (info .ne. 0) STOP 'Error in solver in VACUUM'
!
!       compute tangential covariant (sub u,v) and contravariant
!       (super u,v) magnetic field components on the plasma surface
!
      potu(:nuv2) = zero;  potv(:nuv2) = zero

      mn = 0
      DO n = -nf, nf
         dn2 = -(n*nfper)
         n1 = ABS(n)
         DO m = 0, mf
            mn = mn + 1
            dm2 = m
            DO i = 1, nuv2
               cosmn = cosu1(i,m)*cosv1(i,n1) + csign(n)*
     1                 sinu1(i,m)*sinv1(i,n1)
               potu(i) = potu(i) + dm2*potsin(mn)*cosmn
               potv(i) = potv(i) + dn2*potsin(mn)*cosmn
            END DO
            IF (.not.lasym) CYCLE
            DO i = 1, nuv2
               sinmn = sinu1(i,m)*cosv1(i,n1) - csign(n)*
     1                 cosu1(i,m)*sinv1(i,n1)
               potu(i) = potu(i) - dm2*potcos(mn)*sinmn
               potv(i) = potv(i) - dn2*potcos(mn)*sinmn
            END DO
         END DO
      END DO

      DO i = 1, nuv2
         bsubu(i) = potu(i) + bexu(i)                     !Covariant components
         bsubv(i) = potv(i) + bexv(i)
         huv = p5*guv_b(i)*(nfper)
         hvv = gvv_b(i)*(nfper*nfper)
         det = one/(guu_b(i)*hvv-huv*huv)
         bsupu = (hvv*bsubu(i)-huv*bsubv(i))*det          !Contravariant components
         bsupv = ((-huv*bsubu(i))+guu_b(i)*bsubv(i))*det
         bsqvac(i) = p5*(bsubu(i)*bsupu + bsubv(i)*bsupv)       !.5*|Bvac|**2
         brv(i) = rub(i)*bsupu + rvb(i)*bsupv
         bphiv(i) = r1b(i)*bsupv
         bzv(i) = zub(i)*bsupu + zvb(i)*bsupv
      END DO

!
!       PRINT OUT VACUUM PARAMETERS
!
      IF (ivac .eq. 0) THEN
         ivac = ivac + 1
         IF (lscreen) WRITE (*, 200) nfper, mf, nf, nu, nv
         WRITE (nthreed, 200) nfper, mf, nf, nu, nv
  200    FORMAT(/,2x,'In VACUUM, np =',i3,2x,'mf =',i3,2x,'nf =',i3,
     1      ' nu =',i3,2x,'nv = ',i4)
         bsubuvac = SUM(bsubu(:nuv2)*wint(ns:ns*nuv2:ns))*signgs*pi2    !-plasma current/pi2
         bsubvvac = SUM(bsubv(:nuv2)*wint(ns:ns*nuv2:ns))
         fac = 1.e-6_dp/mu0
         IF (lscreen )WRITE (*,1000) bsubuvac*fac,
     1       plascur*fac, bsubvvac, rbtor
         WRITE (nthreed, 1000) bsubuvac*fac, plascur*fac,
     1       bsubvvac, rbtor
 1000    FORMAT(2x,'2*pi * a * -BPOL(vac) = ',1p,e10.2,
     1      ' TOROIDAL CURRENT = ',e10.2,/,2x,'R * BTOR(vac) = ',
     2      e10.2,' R * BTOR(plasma) = ',e10.2)
!  JDH Add test for RFP. 2013-11-25
!         IF (rbtor*bsubvvac .lt. zero) ier_flag = phiedge_error_flag
!         IF (ABS((plascur - bsubuvac)/rbtor) .gt. 1.e-2_dp)
!     1      ier_flag = 10
         IF (rbtor*bsubvvac .lt. zero) THEN
            IF (lrfp) THEN
 1100    FORMAT(2x,'lrfp is TRUE. Ignore phiedge sign problem')
               IF (lscreen) WRITE(*,1100) 
               WRITE(nthreed,1100) 
            ELSE
               ier_flag = phiedge_error_flag
            ENDIF
         ENDIF
         IF (ABS((plascur - bsubuvac)/rbtor) .gt. 1.e-2_dp) THEN
            IF (lrfp) THEN
1200     FORMAT(2x,'lrfp is TRUE. Proceed with convergence')
               IF (lscreen) WRITE(*,1200) 
               WRITE(nthreed,1200) 
            ELSE
               ier_flag = 10
            ENDIF
         ENDIF
!  END JDH Add test for RFP. 2013-11-25
      ENDIF

      IF (ALLOCATED(bexu))
     1    DEALLOCATE (bexu, bexv, bexn, bexni, r1b, rub, rvb, z1b, zub,
     2    zvb, auu, auv, avv, snr, snv, snz, drv, guu_b, guv_b, gvv_b,
     3    rzb2, rcosuv, rsinuv, stat=i)
      IF (i .ne. 0) STOP 'Deallocation error in vacuum'

      DEALLOCATE (amatrix, bsubu, bsubv, potu, potv, stat = i)

      END SUBROUTINE vacuum
