
      SUBROUTINE scaleup_boundary(dr,mnmax,xm,xn,rmnc_p,zmns_p,rmnc_ws,zmns_ws,&
                                   rmns_p,zmnc_p,rmns_ws,zmnc_ws)
      USE stel_kinds, ONLY: rprec
      USE fieldlines_grid
      implicit none
!-----------------------------------------------------------------------
!     Input Variables
!-----------------------------------------------------------------------
      INTEGER :: mnmax
      INTEGER, DIMENSION(mnmax) :: xm, xn
      REAL(rprec) :: dr
      REAL(rprec), DIMENSION(mnmax) :: rmnc_p, zmns_p
      REAL(rprec), DIMENSION(mnmax) :: rmnc_ws, zmns_ws
      REAL(rprec), INTENT(in), OPTIONAL :: rmns_p(mnmax), zmnc_p(mnmax)
      REAL(rprec), INTENT(out), OPTIONAL :: rmns_ws(mnmax), zmnc_ws(mnmax)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec) :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: m, n, mn, n1, i, nuvb, iflag, nub, nvb
      integer :: ku, kv, k, nmax, mmax, mpol, ntor, iout
      real(rprec) :: cosnv, sinnv, alub, alvb, twopi, rtemp, ztemp
      real(rprec) :: rtemps,ztempc, dr_hold, dist
      real(rprec) :: cosmu, sinmu, theta, zeta, cosmn1, sinmn1, dnorm,&
         b1, c1, r1, re, ae, error, Map_v, sep_tol
      real(rprec), allocatable, dimension(:) :: rbn, zbn
!-----------------------------------------------
      external Map_v
      
      ! Handle passing various arrays though torlines_backgroud
      mnmax_m = mnmax
      dr_m = dr
      IF (ALLOCATED(rmnc_m)) DEALLOCATE(rmnc_m)
      IF (ALLOCATED(zmns_m)) DEALLOCATE(zmns_m)
      IF (ALLOCATED(ixm_m)) DEALLOCATE(ixm_m)
      IF (ALLOCATED(ixn_m)) DEALLOCATE(ixn_m)
      IF (ALLOCATED(rmns_m)) DEALLOCATE(rmns_m)
      IF (ALLOCATED(zmnc_m)) DEALLOCATE(zmnc_m)
      ALLOCATE(rmnc_m(mnmax_m))
      ALLOCATE(zmns_m(mnmax_m))
      ALLOCATE(ixm_m(mnmax_m))
      ALLOCATE(ixn_m(mnmax_m))
      ALLOCATE(rmns_m(mnmax_m))
      ALLOCATE(zmnc_m(mnmax_m))
      
      
      rmnc_m = 0.0; zmns_m = 0.0; rmns_m = 0.0; zmnc_m = 0.0
      rmnc_m = rmnc_p
      zmns_m = zmns_p
      ixm_m = xm
      ixn_m = -xn/nfp_m  ! Assume VMEC format
      IF (PRESENT(rmns_p)) THEN
         rmns_m = rmns_p
         zmnc_m = zmnc_p
      END IF

         
      
!
!     Attempts to scale up boundary by a uniform distance coil_separation
!     On entry, rmnc_b, zmns_b contain the boundary coefficients
!     On exit,  they contain the scaled up boundary coefficients
!
      mpol = MAXVAL(ixm_m)
      ntor = MAXVAL(ixn_m)
      nub  = 6*mpol
      nvb  = 6*ntor
      IF (nvb <= 1) nvb = 1
      nuvb = nub * nvb
      twopi = 8*atan(one)
      alub = twopi/nub
      alvb = twopi/nvb                    !!Note: ixn has nfp factored out already...
      
     

!      allocate (rmnc_ws(nuvb), zmns_ws(nuvb), rmns_ws(nuvb),
!     1          zmnc_ws(nuvb), ixm_ws(nuvb), ixn_ws(nuvb), stat=iflag)
!      if (iflag .ne. 0)stop 'Allocation error in bnorm scaleup_boundary'

   
      ! Check sign of Jacobian
      u0 = twopi/360
      v0 = 0.0; rb_ws = 0.0; zb_ws = 0.0
      isgn_m = 1
      DO mn = 1, mnmax_m
         IF (ixm_m(mn) == 0) CYCLE
         rb_ws = rb_ws + rmnc_m(mn)*cos(ixm_m(mn)*u0)&
                       + rmns_m(mn)*sin(ixm_m(mn)*u0)
         zb_ws = zb_ws + zmns_m(mn)*sin(ixm_m(mn)*u0)&
                       + zmnc_m(mn)*cos(ixm_m(mn)*u0)
      END DO
      IF (rb_ws*zb_ws < 0) isgn_m = -1
         
      IF (ALLOCATED(rbn)) DEALLOCATE(rbn,zbn)
      allocate (rbn(nuvb), zbn(nuvb), stat=iflag)
      if (iflag .ne. 0)stop 'Allocation error in bnorm scaleup_boundary'

      re = 10*epsilon(re)
      ae = re
      i = 0

      IF (ntor /= 0) THEN
         do ku = 1, nub
           u0 = alub*(ku - 1)
           do kv = 1, nvb
             i = i + 1
             v0 = alvb*(kv - 1)                                            !!Np*(Real toroidal angle)
             b1 = v0 - twopi/6
             c1 = v0 + twopi/6
             r1 = v0
             iflag = 0
   !         Map v0 at plasma boundary to vb_ws at winding surface
             call fzero(Map_v, b1, c1, r1, re, ae, iflag)
             error = abs(Map_v (b1))
!             if (iflag .gt. 2) print *,'  i = ', i,' iflag = ', iflag,&
!             ' v0 = ', v0,' vb_ws = ', vb_ws, ' v1 = ', b1,&
!             ' Error mapping to Winding Surface in BNORM code ', error
             rbn(i) = rb_ws
             zbn(i) = zb_ws
           end do
         end do
      ELSE
         do ku = 1, nub
           u0 = alub*(ku - 1)
           do kv = 1, nvb
             i = i + 1
             v0 = alvb*(kv - 1)    
             CALL normal_vector(u0,v0)
             rbn(i) = rb_ws
             zbn(i) = zb_ws
           end do
         end do
      END IF

!
!     FFT new surface (use NESCOIL convention, mu + nv)
!
   
      rmnc_ws = zero;  zmns_ws = zero
      sep_tol = dr/1000
      
      DO mn = 1, mnmax_m
         m = ixm_m(mn)
         n = ixn_m(mn)
         dnorm = one/nuvb
         if (m.ne.0 .or. n.ne.0) dnorm = 2*dnorm
         i = 0
         rtemp = zero;  ztemp = zero; rtemps = zero; ztempc = zero;
         do ku = 1, nub
            theta= alub*(ku-1)
            cosmu = cos(m*theta)*dnorm
            sinmu = sin(m*theta)*dnorm
            do kv = 1, nvb
               i = i + 1
               zeta = alvb*(kv-1)
               cosnv = cos(n*zeta)
               sinnv = sin(n*zeta)
               cosmn1 = cosmu*cosnv - sinmu*sinnv          !cos(mu+nv) NESCOIL CONVENTION
               sinmn1 = sinmu*cosnv + cosmu*sinnv          !sin(mu+nv)
               !cosmn1 = cosmu*cosnv + sinmu*sinnv          !cos(mu-nv) VMEC CONVENTION
               !sinmn1 = sinmu*cosnv - cosmu*sinnv          !sin(mu-nv)
               rtemp  = rtemp  + rbn(i) * cosmn1
               ztemp  = ztemp  + zbn(i) * sinmn1
               rtemps = rtemps + rbn(i) * sinmn1
               ztempc = ztempc + zbn(i) * cosmn1
            end do
         end do
         if (abs(rtemp).lt.sep_tol .and. abs(ztemp).lt.sep_tol) CYCLE
         rmnc_ws(mn) = rtemp
         zmns_ws(mn) = ztemp
         IF (PRESENT(rmns_ws)) THEN
            rmns_ws(mn) = rtemps
            zmnc_ws(mn) = ztempc
         END IF
      END DO

      deallocate (rbn, zbn)

      end subroutine scaleup_boundary
