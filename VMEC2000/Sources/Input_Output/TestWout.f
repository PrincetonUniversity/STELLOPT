      SUBROUTINE TestWout(rzl_array, br, bz, bsupu, bsupv)
      USE stel_kinds, ONLY: rprec, dp
      USE stel_constants, ONLY: twopi
      USE read_wout_mod, ONLY: read_wout_file, read_wout_deallocate, 
     1                         tosuvspace, currumnc, currvmnc,
     2                         jcuru, jcurv, xm_nyq, xn_nyq, 
     3                         ns_w => ns, mnmax_w => mnmax_nyq,
     4                         xm_w => xm_nyq, xn_w => xn_nyq
      USE vmec_params, ONLY: mscale, nscale, ntmax
      USE vmec_dim, ONLY: ns, ntheta3, ntheta1, nrzt
      USE vmec_input, ONLY: mpol, nzeta, ntor, lasym, nfp
      USE vmec_main, ONLY: hs, lthreed
      USE realspace, ONLY: r1, z1, shalf, sqrts
      USE vforces
      USE vmec_utils, ONLY: getbcyl
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(in)  :: 
     1             rzl_array(ns,0:ntor,0:mpol,2*ntmax)
      REAL(rprec), DIMENSION(ns,nzeta,ntheta3), INTENT(in) 
     1                         :: bsupu, bsupv
      REAL(rprec), DIMENSION(nrzt), INTENT(in) :: br, bz
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, js, k, l, istat, mn, mn0
      REAL(rprec) :: si, ui, vi, r12, z12,
     1               sflx, uflx, rfull, zfull
      REAL(rprec) :: br12, bz12, bphi12, 
     1               br1f, bz1f, bphi1f, br2f, bz2f, bphi2f
      LOGICAL :: ls_mesh = .false.               !=T, full mesh, =F, half mesh
C-----------------------------------------------
!      CALL read_wout_file("tichman", istat)
      CALL read_wout_file("fix", istat)
!      CALL read_wout_file("atf14", istat)

!
!     CHECK ON currumn, currvmn vs jcuru, jcurv
!
      DO mn = 1, mnmax_w
         IF (NINT(xm_w(mn)).eq.0 .and. NINT(xn_w(mn)).eq.0) THEN
            mn0 = mn
            EXIT
         END IF
      END DO   

      WRITE (30, *)' mnmax_nyq = ', mnmax_w,' ns = ', ns
      WRITE (30, *)
      WRITE (30, *)'  JS     JCURU       CURRU       JCURV       CURRV'
      DO js = 2, ns_w
         WRITE (30, 123) js, jcuru(js), currumnc(mn0,js),
     1                       jcurv(js), currvmnc(mn0,js)
      END DO
 123  FORMAT (i5, 1p,4e12.3)   

      RETURN
      IF (istat .ne. 0) RETURN

!      gsqrt => azmn_o    !gsqrt
!      r12 => armn_o; bphi => czmn_o
!      bsupu => crmn_e   
!      bsupv => czmn_e   

      js = 3*ns/4
!      js = 2                   !!Check this: not working as well as it should at js=2, js>=3 is EXCELLENT!!!
      si = hs*(js-1.5_dp)

      IF (ls_mesh) THEN
        WRITE(33,'(3a)')
     1  '   L  sflx-del  uflx-del      br1/2     bphi1/2       bz1/2',
     2  '     brf-wout  bphif-wout    bzf-wout    brf-vmec  bphif-vmec',
     3  '    bzf-vmec'
      END IF

      DO k = 1, nzeta
         vi = (twopi*(k-1))/nzeta
         DO j = 1, ntheta3
            ui = (twopi*(j-1))/ntheta1
            l = js + ns*(k-1 + nzeta*(j-1))

            IF (ls_mesh) THEN
            
               rfull = r1(l,0) + sqrts(l)*r1(l,1)
               zfull = z1(l,0) + sqrts(l)*z1(l,1)           
!
!          TEST WOUT-BASED VERSION
!
               CALL GetBcyl(rfull, vi/nfp, zfull, br1f, bphi1f, bz1f, 
     1                      sflx, uflx, istat)
               IF (istat .ne. 0) PRINT *,' WOUT-BASED ISTAT = ', istat
!
!          TEST VMEC-BASED VERSION 
!
               CALL GetBcyl(rfull, vi/nfp, zfull, br2f, bphi2f, bz2f, 
     1           sflx, uflx, bsupu, bsupv, rzl_array, ns, ntor, mpol, 
     2           ntmax, nzeta, ntheta3, nfp, mscale, nscale, 
     3           lthreed, lasym, istat)
               IF (istat .ne. 0) PRINT *,' VMEC-BASED ISTAT = ', istat
 
               WRITE (33, 1225) l, sflx-si, MOD(uflx-ui,twopi), 
     1            br(l), czmn_o(l), bz(l),
     2            br1f, bphi1f, bz1f, br2f, bphi2f, bz1f
               
 1225   FORMAT(i4, 1p,2e10.2, 9e12.3)

            ELSE

            z12 = (z1(l,0) + z1(l-1,0) + shalf(l)*
     1            (z1(l,1) + z1(l-1,1)))/2
            r12 = armn_o(l)

!
!           TEST WOUT-BASED VERSION
!
            CALL GetBcyl(r12, vi/nfp, z12, br12, bphi12, bz12, 
     1                   sflx, uflx, istat)

            IF (istat .eq. 0) THEN
               WRITE (33, 1224) l, sflx-si, MOD(uflx-ui,twopi), 
     1            (br(l) - br12)/bphi12, (czmn_o(l)-bphi12)/bphi12, 
     2            (bz(l) - bz12)/bphi12
            ELSE 
               WRITE (33, *)' istat = ', istat
            END IF


!
!          TEST VMEC-BASED VERSION (OVERLOADED, SAME FUNCTION CALL, DIFF ARG.)
!
            CALL GetBcyl(r12, vi/nfp, z12, br12, bphi12, bz12, sflx, 
     1           uflx, bsupu, bsupv, rzl_array, ns, ntor, mpol, 
     2           ntmax, nzeta, ntheta3, nfp, mscale, nscale, 
     3           lthreed, lasym, istat)
            IF (istat .eq. 0) THEN
               WRITE (34, 1224) l, sflx-si, MOD(uflx-ui,twopi), 
     1            (br(l) - br12)/bphi12, (czmn_o(l)-bphi12)/bphi12, 
     2            (bz(l) - bz12)/bphi12
            ELSE
               WRITE (34, *)' istat = ', istat
            END IF

            END IF

!            CALL tosuvspace(si,ui,vi,gsqrt=gsqrt1)
!            WRITE(33, 1223) l, gsqrt(l), gsqrt1,
!     1                 ABS(gsqrt(l)-gsqrt1)/ABS(gsqrt1)
         END DO
      END DO
 1223 FORMAT(' l = ',i4,' GSQRT = ',1p,e12.4,' FROM WOUT = ',e12.4,
     1      ' REL.ERROR = ',e10.2)
 1224 FORMAT(' l = ',i4, ' del-SFLX = ', 1p,e10.2, 
     1                   ' del-UFLX = ', e10.2,
     2                   ' del-BR = ', e12.4,' del-BPHI = ', 
     3            e12.4, ' del-BZ = ', e12.4)

      CALL read_wout_deallocate

      END SUBROUTINE TestWout
