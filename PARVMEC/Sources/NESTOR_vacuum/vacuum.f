#if defined(SKS)
      SUBROUTINE vacuum_par(rmnc, rmns, zmns, zmnc, xm, xn, 
     1                  plascur, rbtor, wint, ns, ivac_skip, ivac, 
     2                  mnmax, ier_flag, lscreen)
      USE vacmod
      USE vparams, ONLY: nthreed, zero, one, mu0
      USE vmec_params, ONLY: norm_term_flag, phiedge_error_flag
      USE vmec_input, ONLY: lrfp         ! JDH Added 2013-11-25, to test for RFP
      USE vmec_main, ONLY: nznt, irst
      USE parallel_include_module
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ns, ivac_skip, ivac, mnmax, ier_flag
      REAL(rprec) :: plascur, rbtor
      REAL(rprec), DIMENSION(mnmax), INTENT(in) ::
     1   rmnc, rmns, zmns, zmnc, xm, xn
      REAL(rprec), DIMENSION(nuv3), INTENT(in) :: wint
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, n, n1, m, i, info, j
      REAL(rprec), DIMENSION(:), POINTER :: potcos, potsin
      REAL(rprec), ALLOCATABLE :: potu(:), potv(:)
      REAL(rprec), ALLOCATABLE :: amatrix(:)
      REAL(rprec):: dn2, dm2, cosmn, sinmn, huv, hvv,
     1              det, bsubuvac, fac, ton, toff

      REAL(rprec) :: tmp1(2), tmp2(2)
      REAL(rprec) :: skston, skstoff
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
      CALL second0(skston)
      IF (.not.ALLOCATED(potvac)) STOP 'POTVAC not ALLOCATED in VACCUM'

      ALLOCATE (amatrix(mnpd2*mnpd2), potu(nuv3), potv(nuv3), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in vacuum'

      potsin => potvac(1:mnpd)
      potcos => potvac(1+mnpd:)

      ALLOCATE (bexu(nuv3), bexv(nuv3), bexn(nuv3),
     1     bexni(nuv3), r1b(nuv), rub(nuv3), rvb(nuv3),
     2     z1b(nuv), zub(nuv3), zvb(nuv3), auu(nuv3), auv(nuv3),
     3     avv(nuv3), snr(nuv3), snv(nuv3), snz(nuv3), drv(nuv3),
     4     guu_b(nuv3), guv_b(nuv3), gvv_b(nuv3), rzb2(nuv),
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
      CALL second0(ton)
      IF (.not. ALLOCATED(tanu)) CALL precal_par (wint)
      CALL surface_par (rmnc, rmns, zmns, zmnc, xm, xn, mnmax)
      CALL second0(toff)
      timer_vac(tsurf) = timer_vac(tsurf) + (toff-ton)

      ton = toff
      CALL bextern_par (plascur, wint,lscreen)
      CALL second0(toff)
      timer_vac(tbext) = timer_vac(tbext)+(toff-ton)
!
!     Determine scalar magnetic potential POTVAC
!
      CALL second0(ton)
      CALL scalpot_par (potvac, amatrix, wint, ivac_skip)
      CALL second0(toff)
      timer_vac(tscal) = timer_vac(tscal) + (toff-ton)

      ton = toff
      CALL solver (amatrix, potvac, mnpd2, 1, info)
      CALL second0(toff)
      timer_vac(tsolver) = timer_vac(tsolver) + (toff-ton)
      solver_time = solver_time + (toff - ton)

      IF (info .ne. 0) STOP 'Error in solver in VACUUM'
!
!       compute tangential covariant (sub u,v) and contravariant
!       (super u,v) magnetic field components on the plasma surface
!
      potu(nuv3min:nuv3max) = 0;  potv(nuv3min:nuv3max) = 0

      mn = 0
      DO n = -nf, nf
         dn2 = -(n*nfper)
         n1 = ABS(n)
         DO m = 0, mf
            mn = mn + 1
            dm2 = m
            j = 0
            DO i = nuv3min, nuv3max
               j = j+1
               cosmn = potsin(mn)*cosmni(j,mn)/(pi2*pi2*wint(i))
               potu(i) = potu(i) + dm2*cosmn
               potv(i) = potv(i) + dn2*cosmn
            END DO
            IF (.not.lasym) CYCLE
            j = 0
            DO i = nuv3min, nuv3max
               j = j+1
               sinmn = potcos(mn)*sinmni(j,mn)/(pi2*pi2*wint(i))
               potu(i) = potu(i) - dm2*sinmn
               potv(i) = potv(i) - dn2*sinmn
            END DO
         END DO
      END DO

      DO i = nuv3min, nuv3max
         bsubu_sur(i) = potu(i) + bexu(i)                               !Covariant components
         bsubv_sur(i) = potv(i) + bexv(i)
         huv = p5*guv_b(i)*(nfper)
         hvv = gvv_b(i)*(nfper*nfper)
         det = one/(guu_b(i)*hvv-huv*huv)
         bsupu_sur(i) = (hvv*bsubu_sur(i)-huv*bsubv_sur(i))*det         !Contravariant components
         bsupv_sur(i) = ((-huv*bsubu_sur(i))+guu_b(i)*bsubv_sur(i))*det
         bsqvac(i) = p5*(bsubu_sur(i)*bsupu_sur(i) 
     1             +     bsubv_sur(i)*bsupv_sur(i))                     !.5*|Bvac|**2
         brv(i) = rub(i)*bsupu_sur(i) + rvb(i)*bsupv_sur(i)
         bphiv(i) = r1b(i)*bsupv_sur(i)
         bzv(i) = zub(i)*bsupu_sur(i) + zvb(i)*bsupv_sur(i)
      END DO

!
!       PRINT OUT VACUUM PARAMETERS
!
      IF (ivac .eq. 0) THEN
         ivac = ivac + 1
         IF (vrank .EQ. 0) THEN
            IF (lscreen) WRITE (*, 200) nfper, mf, nf, nu, nv
            WRITE (nthreed, 200) nfper, mf, nf, nu, nv
         END IF
  200    FORMAT(/,2x,'In VACUUM, np =',i3,2x,'mf =',i3,2x,'nf =',i3,
     1      ' nu =',i3,2x,'nv = ',i4)

         bsubuvac = 0
         bsubvvac = 0
         DO i=nuv3min, nuv3max
           bsubuvac = bsubuvac + bsubu_sur(i)*wint(i)
           bsubvvac = bsubvvac + bsubv_sur(i)*wint(i)
         END DO
         tmp1(1)=bsubuvac; tmp1(2)=bsubvvac
         CALL second0(ton)
         IF (vlactive) THEN
           CALL MPI_Allreduce(tmp1,tmp2,2,MPI_REAL8,MPI_SUM,VAC_COMM,
     1                        MPI_ERR)
         END IF
         CALL second0(toff)
         allreduce_time = allreduce_time + (toff - ton)
         bsubuvac = tmp2(1); bsubvvac = tmp2(2)
         bsubuvac = bsubuvac*signgs*pi2

         fac = 1.e-6_dp/mu0
         IF (vrank .EQ. 0) THEN
         IF (lscreen ) WRITE (*,1000) bsubuvac*fac,
     1       plascur*fac, bsubvvac, rbtor
         WRITE (nthreed, 1000) bsubuvac*fac, plascur*fac,
     1       bsubvvac, rbtor
         END IF
 1000    FORMAT(2x,'2*pi * a * -BPOL(vac) = ',1p,e10.2,
     1      ' TOROIDAL CURRENT = ',e10.2,/,2x,'R * BTOR(vac) = ',
     2      e10.2,' R * BTOR(plasma) = ',e10.2)
!  JDH Add test for RFP. 2013-11-25
!         IF (rbtor*bsubvvac .lt. zero) ier_flag = phiedge_error_flag
!         IF (ABS((plascur - bsubuvac)/rbtor) .gt. 1.e-2_dp)
!     1      ier_flag = 10
         IF (rbtor*bsubvvac .lt. zero) THEN
            IF (lrfp) THEN
               IF (vrank .EQ. 0) THEN
               IF (lscreen) WRITE(*,1100) 
               WRITE(nthreed,1100) 
               END IF
            ELSE
               ier_flag = phiedge_error_flag
            ENDIF
         ENDIF
         IF (ABS((plascur - bsubuvac)/rbtor) .gt. 1.e-2_dp) THEN
            IF (lrfp) THEN
               IF (vrank .EQ. 0) THEN
               IF (lscreen) WRITE(*,1200) 
               WRITE(nthreed,1200) 
               END IF
            ELSE
               ier_flag = 10
            ENDIF
         ENDIF
!  END JDH Add test for RFP. 2013-11-25
      ENDIF
1100  FORMAT('lrfp is TRUE. Ignore phiedge sign problem')
1200  FORMAT('lrfp is TRUE. Proceed with convergence')

      IF (ALLOCATED(bexu))
     1    DEALLOCATE (bexu, bexv, bexn, bexni, r1b, rub, rvb, z1b, zub,
     2    zvb, auu, auv, avv, snr, snv, snz, drv, guu_b, guv_b, gvv_b,
     3    rzb2, rcosuv, rsinuv, stat=i)
      IF (i .ne. 0) STOP 'Deallocation error in vacuum'

      DEALLOCATE (amatrix, potu, potv, stat=i)
      IF (i .ne. 0) STOP 'Deallocation error in vacuum'
      
      CALL second0(ton)
      IF (vlactive) THEN
        CALL MPI_Allgatherv(MPI_IN_PLACE,numjs_vac,MPI_REAL8,bsqvac,
     1          counts_vac,disps_vac,MPI_REAL8,VAC_COMM,MPI_ERR)
      END IF
      CALL second0(toff)
      timer_vac(tallgv) = timer_vac(tallgv) + (toff-ton)

      skstoff = toff
      vacuum_time = vacuum_time + (skstoff - skston)
      END SUBROUTINE vacuum_par
#endif


      SUBROUTINE vacuum(rmnc, rmns, zmns, zmnc, xm, xn, 
     1                  plascur, rbtor, wint, ns, ivac_skip, ivac, 
     2                  mnmax, ier_flag, lscreen)
      USE vacmod
      USE vparams, ONLY: nthreed, zero, one, mu0
      USE vmec_params, ONLY: norm_term_flag, phiedge_error_flag
      USE vmec_input, ONLY: lrfp         ! JDH Added 2013-11-25, to test for RFP
#if defined(SKS)
      USE parallel_include_module
#endif
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
      REAL(rprec), ALLOCATABLE :: potu(:), potv(:)
      REAL(rprec), ALLOCATABLE :: amatrix(:)
      REAL(rprec):: dn2, dm2, cosmn, sinmn, huv, hvv,
     1    det, bsubuvac, fac
#if defined(SKS)
      REAL(rprec) :: skston, skstoff, ton, toff
#endif
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
#if defined(SKS)
      CALL second0(skston)
#endif
      ier_flag = norm_term_flag

      IF (.not.ALLOCATED(potvac)) STOP 'POTVAC not ALLOCATED in VACCUM'

      ALLOCATE (amatrix(mnpd2*mnpd2), potu(nuv3), potv(nuv3), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in vacuum'

      potsin => potvac(1:mnpd)
      potcos => potvac(1+mnpd:)

      ALLOCATE (bexu(nuv3), bexv(nuv3), bexn(nuv3),
     1     bexni(nuv3), r1b(nuv), rub(nuv3), rvb(nuv3),
     2     z1b(nuv), zub(nuv3), zvb(nuv3), auu(nuv3), auv(nuv3),
     3     avv(nuv3), snr(nuv3), snv(nuv3), snz(nuv3), drv(nuv3),
     4     guu_b(nuv3), guv_b(nuv3), gvv_b(nuv3), rzb2(nuv),
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
#if defined(SKS)
      CALL second0(ton)
#endif
      CALL solver (amatrix, potvac, mnpd2, 1, info)
#if defined(SKS)
      CALL second0(toff)
      s_solver_time = s_solver_time + (toff - ton)
#endif

      IF (info .ne. 0) STOP 'Error in solver in VACUUM'
!
!       compute tangential covariant (sub u,v) and contravariant
!       (super u,v) magnetic field components on the plasma surface
!
      potu(:nuv3) = zero;  potv(:nuv3) = zero

      mn = 0
      DO n = -nf, nf
         dn2 = -(n*nfper)
         n1 = ABS(n)
         DO m = 0, mf
            mn = mn + 1
            dm2 = m
            DO i = 1, nuv3
               cosmn = cosu1(i,m)*cosv1(i,n1) + csign(n)*
     1                 sinu1(i,m)*sinv1(i,n1)
               potu(i) = potu(i) + dm2*potsin(mn)*cosmn
               potv(i) = potv(i) + dn2*potsin(mn)*cosmn
            END DO
            IF (.not.lasym) CYCLE
            DO i = 1, nuv3
               sinmn = sinu1(i,m)*cosv1(i,n1) - csign(n)*
     1                 cosu1(i,m)*sinv1(i,n1)
               potu(i) = potu(i) - dm2*potcos(mn)*sinmn
               potv(i) = potv(i) - dn2*potcos(mn)*sinmn
            END DO
         END DO
      END DO

      DO i = 1, nuv3
         bsubu_sur(i) = potu(i) + bexu(i)                 !Covariant components
         bsubv_sur(i) = potv(i) + bexv(i)
         huv = p5*guv_b(i)*(nfper)
         hvv = gvv_b(i)*(nfper*nfper)
         det = one/(guu_b(i)*hvv-huv*huv)
         bsupu_sur(i) = (hvv*bsubu_sur(i)-huv*bsubv_sur(i))*det  !Contravariant components
         bsupv_sur(i) = ((-huv*bsubu_sur(i))+guu_b(i)*bsubv_sur(i))*det
         bsqvac(i) = p5*(bsubu_sur(i)*bsupu_sur(i) 
     1                 + bsubv_sur(i)*bsupv_sur(i))       !.5*|Bvac|**2
         brv(i) = rub(i)*bsupu_sur(i) + rvb(i)*bsupv_sur(i)
         bphiv(i) = r1b(i)*bsupv_sur(i)
         bzv(i) = zub(i)*bsupu_sur(i) + zvb(i)*bsupv_sur(i)
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
         bsubuvac = SUM(bsubu_sur(:nuv3)*wint(ns:ns*nuv3:ns))*signgs*pi2    !-plasma current/pi2
         bsubvvac = SUM(bsubv_sur(:nuv3)*wint(ns:ns*nuv3:ns))

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
1100  FORMAT('lrfp is TRUE. Ignore phiedge sign problem')
               IF (lscreen) WRITE(*,1100) 
               WRITE(nthreed,1100) 
            ELSE
               ier_flag = phiedge_error_flag
            ENDIF
         ENDIF
         IF (ABS((plascur - bsubuvac)/rbtor) .gt. 1.e-2_dp) THEN
            IF (lrfp) THEN
1200  FORMAT('lrfp is TRUE. Proceed with convergence')
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

      DEALLOCATE (amatrix, potu, potv, stat = i)

#if defined(SKS)
      CALL second0(skstoff)
      s_vacuum_time = s_vacuum_time + (skstoff - skston)
#endif
      END SUBROUTINE vacuum
