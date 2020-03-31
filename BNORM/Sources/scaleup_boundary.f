
      subroutine scaleup_boundary(nub, nvb)
      use neswrite
!      use meshes, only: md, nd
      use normal_info, only: u0, v0, rb_ws, zb_ws, vb_ws
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nub, nvb
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: m, n, mn, n1, i, nuvb, iflag
      integer :: ku, kv, k, nmax, mmax
      real(rprec) :: cosnv, sinnv, alub, alvb, twopi, rtemp, ztemp
      real(rprec) :: rtemps,ztempc
      real(rprec) :: cosmu, sinmu, theta, zeta, cosmn1, sinmn1, dnorm,
     1   b1, c1, r1, re, ae, error, Map_v, sep_tol
      real(rprec), allocatable, dimension(:) :: rbn, zbn
C-----------------------------------------------
      external Map_v
!
!     Attempts to scale up boundary by a uniform distance coil_separation
!     On entry, rmnc_b, zmns_b contain the boundary coefficients
!     On exit,  they contain the scaled up boundary coefficients
!
      nuvb = nub * nvb
      twopi = 8*atan(one)
      alub = twopi/nub
      alvb = twopi/nvb                    !!Note: ixn has nfp factored out already...

      allocate (rmnc_ws(nuvb), zmns_ws(nuvb), rmns_ws(nuvb),
     1          zmnc_ws(nuvb), ixm_ws(nuvb), ixn_ws(nuvb), stat=iflag)
      if (iflag .ne. 0)stop 'Allocation error in bnorm scaleup_boundary'

      allocate (rbn(nuvb), zbn(nuvb), stat=iflag)
      if (iflag .ne. 0)stop 'Allocation error in bnorm scaleup_boundary'

      re = 10*epsilon(re)
      ae = re
      i = 0

      do ku = 1, nub
        u0 = alub*(ku - 1)
        do kv = 1, nvb
          i = i + 1
          v0 = alvb*(kv - 1)                                            !!Np*(Real toroidal angle)
          b1 = v0 - twopi/6
          c1 = v0 + twopi/6
          r1 = v0
!         Map v0 at plasma boundary to vb_ws at winding surface
          call fzero(Map_v, b1, c1, r1, re, ae, iflag)
          error = abs(Map_v (b1))
          if (iflag .gt. 2) print *,'  i = ', i,' iflag = ', iflag,
     1    ' v0 = ', v0,' vb_ws = ', vb_ws, ' v1 = ', b1,
     2    ' Error mapping to Winding Surface in BNORM code ', error

          rbn(i) = rb_ws
          zbn(i) = zb_ws
        end do
      end do

!
!     FFT new surface (use NESCOIL convention, mu + nv)
!
      mmax = min(md, (nub-1)/2)
      nmax = min(nd, abs(nvb-1)/2)
      mnmax_ws = 0
      rmnc_ws = zero;  zmns_ws = zero
      sep_tol = 1.e-3_dp*coil_separation

      mloop: do m = 0, mmax
         nloop: do n = -nmax, nmax
            if (m.eq.0 .and. n.gt.0) cycle nloop
            dnorm = one/nuvb
            if (m.ne.0 .or. n.ne.0) dnorm = 2*dnorm
            i = 0
            rtemp = zero;  ztemp = zero
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
                  rtemp  = rtemp  + rbn(i) * cosmn1
                  ztemp  = ztemp  + zbn(i) * sinmn1
                  rtemps = rtemps + rbn(i) * sinmn1
                  ztempc = ztempc + zbn(i) * cosmn1
               end do
            end do
            if (abs(rtemp).lt.sep_tol .and. abs(ztemp).lt.sep_tol)
     1      cycle nloop
            mnmax_ws = mnmax_ws+1
            rmnc_ws(mnmax_ws) = rtemp
            zmns_ws(mnmax_ws) = ztemp
            rmns_ws(mnmax_ws) = rtemps
            zmnc_ws(mnmax_ws) = ztempc
            ixm_ws(mnmax_ws) = m
            ixn_ws(mnmax_ws) = n
         end do nloop
      end do mloop

      deallocate (rbn, zbn)

      end subroutine scaleup_boundary
