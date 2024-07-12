      SUBROUTINE mercier(gsqrt, bsq, bdotj, iotas, wint,
     1           r1, rt, rz, zt, zz, bsubu, vp, phips, pres, ns, nznt)
      USE safe_open_mod
      USE vmercier
      USE vmec_input, ONLY: input_extension
      USE vparams, ONLY: one, zero, twopi, nmercier0
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ns, nznt
      REAL(rprec), DIMENSION(ns,nznt), INTENT(in) ::
     1     gsqrt, bsq
      REAL(rprec), DIMENSION(ns,nznt), INTENT(inout) :: bdotj
      REAL(rprec), DIMENSION(ns*nznt), INTENT(in) :: wint, bsubu
      REAL(rprec), DIMENSION(ns,nznt,0:1), INTENT(in) ::
     1     r1, rt, rz, zt, zz
      REAL(rprec), DIMENSION(ns), INTENT(in) ::
     1     iotas, vp, phips, pres
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 0.5_dp, two = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ns1, i, imercier0, nmerc = nmercier0, nrzt
      REAL(rprec) :: sign_jac, hs, sqs, denom
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1     gpp, gsqrt_full, b2
      REAL(rprec), DIMENSION(nznt) :: gtt, rtf, ztf,
     1   rzf, zzf, r1f, jdotb, ob2, b2i
      REAL(rprec), DIMENSION(ns) :: vp_real, phip_real,
     1   shear, vpp, presp, torcur, ip, sj, tpp, tjb, tbb, tjj
      CHARACTER(LEN=120) mercier_file
C-----------------------------------------------
      nrzt = ns*nznt
      mercier_file = 'mercier.'//TRIM(input_extension)
      CALL safe_open (nmerc, imercier0, mercier_file, 'replace',
     1   'formatted')
      IF (imercier0 .ne. 0) RETURN


      ALLOCATE (gpp(ns,nznt), gsqrt_full(ns,nznt), b2(ns,nznt), stat=i)
      IF (i .ne. 0) STOP 'allocation error in Mercier'

!
!     SCALE VP, PHIPS TO REAL UNITS (VOLUME, TOROIDAL FLUX DERIVATIVES)
!     AND PUT GSQRT IN ABS UNITS (SIGNGS MAY BE NEGATIVE)
!     NOTE: VP has (coming into this routine) the sign of the jacobian multiplied out
!          i.e., vp = signgs*<gsqrt>
!     THE SHEAR TERM MUST BE MULTIPLIED BY THE SIGN OF THE JACOBIAN
!     (OR A BETTER SOLUTION IS TO RETAIN THE JACOBIAN SIGN IN ALL TERMS, INCLUDING
!      VP, THAT DEPEND EXPLICITLY ON THE JACOBIAN. WE CHOOSE THIS LATTER METHOD...)
!
!     COMING INTO THIS ROUTINE, THE JACOBIAN(gsqrt) = 1./(grad-s . grad-theta X grad-zeta)
!     WE CONVERT THIS FROM grad-s to grad-phi DEPENDENCE BY DIVIDING gsqrt by PHIP_real
!
!     NOTE: WE ARE USING 0 < s < 1 AS THE FLUX VARIABLE, BEING CAREFUL
!     TO KEEP d(phi)/ds == PHIP_real FACTORS WHERE REQUIRED
!     THE V'' TERM IS d2V/d(PHI)**2, PHI IS REAL TOROIDAL FLUX
!
!     SHEAR = d(iota)/d(phi)   :  FULL MESH
!     VPP   = d(vp)/d(phi)     :  FULL MESH
!     PRESP = d(pres)/d(phi)   :  FULL MESH  (PRES IS REAL PRES*mu0)
!     IP    = d(Itor)/d(phi)   :  FULL MESH
!
!     ON ENTRY, BDOTJ = Jacobian * J*B  ON THE FULL RADIAL GRID
!               BSQ = 0.5*|B**2| + p IS ON THE HALF RADIAL GRID
!

      ns1 = ns - 1
      IF (ns1 .le. 0) RETURN
      hs = one/ns1
      sign_jac = zero
      IF (gsqrt(ns,1) .ne. zero)
     1    sign_jac = ABS(gsqrt(ns,1))/gsqrt(ns,1)

      IF (sign_jac .eq. zero) RETURN
      phip_real = twopi * phips * sign_jac
!
!     NOTE: phip_real should be > 0 to get the correct physical sign of REAL-space gradients
!     For example, grad-p, grad-Ip, etc. However, with phip_real defined this way,
!     Mercier will be correct
!
      vp_real(2:ns) = sign_jac*(twopi*twopi)*vp(2:ns)/phip_real(2:ns)  !!dV/d(PHI) on half mesh

!
!     COMPUTE INTEGRATED TOROIDAL CURRENT
!
      DO i = 2,ns
         torcur(i)=sign_jac*twopi*SUM(bsubu(i:nrzt:ns)*wint(i:nrzt:ns))
      END DO

!
!     COMPUTE SURFACE AVERAGE VARIABLES ON FULL RADIAL MESH
!
      DO i = 2,ns1
        phip_real(i) = p5*(phip_real(i+1) + phip_REAL(i))
        denom     = one/(hs*phip_real(i))
        shear(i)  = (iotas(i+1) - iotas(i))*denom       !!d(iota)/d(PHI)
        vpp(i)    = (vp_real(i+1) - vp_real(i))*denom   !!d(VP)/d(PHI)
        presp(i)  = (pres(i+1) - pres(i))*denom         !!d(p)/d(PHI)
        ip(i)     = (torcur(i+1) - torcur(i))*denom     !!d(Itor)/d(PHI)
      END DO

!
!     COMPUTE GPP == |grad-phi|**2 = PHIP**2*|grad-s|**2           (on full mesh)
!             Gsqrt_FULL = JACOBIAN/PHIP == jacobian based on flux (on full mesh)
!

      DO i = 2, ns1
        gsqrt_full(i,:) = p5*(gsqrt(i,:) + gsqrt(i+1,:))
        bdotj(i,:) = bdotj(i,:)/gsqrt_full(i,:)
        gsqrt_full(i,:) = gsqrt_full(i,:)/phip_real(i)
        sj(i) = hs*(i-1)
        sqs = SQRT(sj(i))
        rtf(:) = rt(i,:,0) + sqs*rt(i,:,1)
        ztf(:) = zt(i,:,0) + sqs*zt(i,:,1)
        gtt(:) = rtf(:)*rtf(:) + ztf(:)*ztf(:)
        rzf(:) = rz(i,:,0) + sqs*rz(i,:,1)
        zzf(:) = zz(i,:,0) + sqs*zz(i,:,1)
        r1f(:) = r1(i,:,0) + sqs*r1(i,:,1)
        gpp(i,:) = gsqrt_full(i,:)**2/(gtt(:)*r1f(:)**2 +
     1             (rtf(:)*zzf(:) - rzf(:)*ztf(:))**2)     !!1/gpp
      END DO

!
!     COMPUTE SURFACE AVERAGES OVER dS/|grad-PHI|**3 => |Jac| du dv / |grad-PHI|**2
!     WHERE Jac = gsqrt/phip_real
!
      DO i = 2,ns
        b2(i,:) = two*(bsq(i,:) - pres(i))
      END DO

      DO i = 2,ns1
        b2i(:) = p5*(b2(i+1,:) + b2(i,:))
        ob2(:) = gsqrt_full(i,:)/b2i(:)
        tpp(i) = SUM(ob2(:)*wint(i:nrzt:ns))                !<1/B**2>
        ob2(:) = b2i(:) * gsqrt_full(i,:) * gpp(i,:)
        tbb(i) = SUM(ob2(:)*wint(i:nrzt:ns))                !<b*b/|grad-phi|**3>
        jdotb(:) = bdotj(i,:) * gpp(i,:) * gsqrt_full(i,:)
        tjb(i) = SUM(jdotb(:)*wint(i:nrzt:ns))              !<j*b/|grad-phi|**3>
        jdotb(:) = jdotb(:) * bdotj(i,:) / b2i(:)
        tjj(i) = SUM(jdotb(:)*wint(i:nrzt:ns))              !<(j*b)2/b**2*|grad-phi|**3>
      END DO

      DEALLOCATE (gpp, gsqrt_full, b2, stat=i)

!
!     REFERENCE: BAUER, BETANCOURT, GARABEDIAN, MHD Equilibrium and Stability of Stellarators
!     We break up the Omega-subs into a positive shear term (Dshear) and a net current term, Dcurr
!     Omega_subw == Dwell and Omega-subd == Dgeod (geodesic curvature, Pfirsch-Schluter term)
!
!     Include (eventually) Suydam for reference (cylindrical limit)
!

      WRITE(nmerc,90)
 90   FORMAT(6x,'S',10x,'PHI',9x,'IOTA',8x,'SHEAR',7x,' VP ',8x,'WELL',
     1       8x,'ITOR',7x,'ITOR''',7x,'PRES',7x,'PRES''',/,120('-'))

      DO i = 2,ns1
         sqs = p5*(vp_real(i) + vp_real(i+1))*sign_jac
         IF (sqs .eq. zero) CYCLE
         WRITE(nmerc,100) sj(i), hs*SUM(phip_real(2:i)),
     1   p5*(iotas(i+1)+iotas(i)), shear(i)/sqs,
     2   sqs, -vpp(i)*sign_jac,
     3   p5*(torcur(i) + torcur(i+1)), ip(i)/sqs,
     4   p5*(pres(i) + pres(i+1)), presp(i)/sqs
      END DO

 100  FORMAT(1p,10e12.4)

      WRITE(nmerc,190)
 190  FORMAT(/,6x,'S',8x,'DMerc',8x,'DShear',7x,'DCurr',7x,'DWell',
     1     7x,'Dgeod',/,100('-'))

      DO i = 2,ns1
         tpp(i) = (twopi*twopi)*tpp(i)
         tjb(i) = (twopi*twopi)*tjb(i)
         tbb(i) = (twopi*twopi)*tbb(i)
         tjj(i) = (twopi*twopi)*tjj(i)
         Dshear(i) = shear(i) * shear(i)/4
         Dcurr(i)  =-shear(i) * (tjb(i) - ip(i) *tbb(i))
         Dwell(i)  = presp(i) * (vpp(i) - presp(i) *tpp(i))*tbb(i)
         Dgeod(i)  = tjb(i) *tjb(i)  - tbb(i) *tjj(i)
         DMerc(i)  = Dshear(i) + Dcurr(i) + Dwell(i) + Dgeod(i)
         WRITE(nmerc,100) sj(i), Dmerc(i), Dshear(i),
     1         Dcurr(i), Dwell(i), Dgeod(i)
      END DO

      CLOSE (nmerc)

      END SUBROUTINE mercier
