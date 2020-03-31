
! ------------------------------------------------------------
      subroutine min_xerr_phimn(method, eps, nu, nv, nuvh,
     1   nbdim, ndim, bfn, ben, dsur, svdv,svdw,nsvdw )
!................................................................
c     Purpose:
c     Calculate phi(m,n) which minimize xerr (rather than berr)
c................................................................

      use stel_kinds
      USE NumParams
      use OutCtrl, ONLY: w_xerr
      USE LoopCtrl
      use Vvacuum2, ONLY: cl1, iota_edge, phip_edge, ms1, ns1
      USE Vprecal1, ONLY: np
      implicit none
c................................................................
C   D u m m y   A r g u m e n t s
c................................................................
      integer method, numsvd, nu, nv, nuvh, nbdim, ndim, nsvdw
      real(rprec) eps
      real(rprec), dimension(nuvh,nbdim) :: bfn
      real(rprec), dimension(nuvh) :: ben, dsur
      real(rprec), dimension(ndim) :: svdw
      real(rprec), dimension(ndim,ndim) :: svdv
c................................................................
C   L o c a l   V a r i a b l e s
c................................................................
      integer :: mnl, mf, nf, nvh, mnf,
     1   ia, istat, i, m, n, k, kv, nwt, ku, j, nn, ns
      real(rprec) :: epsq, du, hval, u, res,
     1   rmsxerr, absxerr
      real(rprec), allocatable, dimension(:,:,:) :: smunv, cmunv
      real(rprec), allocatable, dimension(:,:) ::
     1    cu, cv, su, sv, xg
      real(rprec), allocatable, dimension(:) ::
     1    lmn, luv, resmn, xh, gval
c................................................................
c
c  Inside-nescoil processor to calculate phi(m,n) current potentials
c  by minimizing displacements by changing Nescoil Green functions to
c  those needed for calculating displacements from phimn.
c
c  Given ben, bfn in (theta,phi), and
c        fourier coeffs lmn of lambda and g from vmec,
c  Calculate the modified Green functions, and
c  solve in fourier space to get phi(m,n)
c...
c  t=theta, p=phi, u = t + lambda(t,p), v = p
c  Straight-line system is (u,v):
c  u   = theta + lambda( theta, phi ),  v = phi
c  t   = theta,  p = phi
c...
c  Inputs:
c  method  - 0:use given svdv to calc xerr, else solve for svdv
c  np  - number of field periods
c  eps  - Resonance broadening parameter
c
c  nu   - number of points used in u (poloidal) grid
c  nv   - number of points used in v (toroidal) grid in 1 period
c  nuvh - number of points in half period + edge = nu*(1+nv/2)
c  nbdim-
c  ben  -
c  bfn  -
c  dsur -
c...
c  Outputs:
c  phimn - Fourier coeffs of phi which minimize xerror
c  numsvd- Number of svd weights to use
c...
c  Local Arrays and variables:
c
c  Following come from VMEC:
c  ml   - number of m values for the lambda coefficients (0:ml)
c  nl   - number of n values for the lambda coefficients (0:nl)
c  mnl  - length of lmn array = 1+nl+ml*(2*nl+1))
c  lmn  - Fourier coefficients of the stream function lambda from VMEC
c  phip_edge - radial derivative of toroidal flux at this surface
c  iota_edge - iota on this plasma surface
c
c  Following are used for fourier transform calculations:
c  smunv, cmunv, su, cu, sv, cv
c
c................................................................
c
c  Allocate local storage for arrays on nuvh spatial halfgrid
c  and for cos/sin factors on the largest fourier grid mf,nf
      mf = nu/2
      nf = nv/2
      mnf = 1 + nf + mf*(2*nf + 1)
      nvh = 1 + nv/2                             !half+centerline
      nuvh = nu*nvh
      istat=0
      if (.not.allocated(smunv)) allocate(
     1   smunv(nuvh,0:mf,-nf:nf), cmunv(nuvh,0:mf,-nf:nf),
     2   cu(nu,0:mf),cv(nvh,0:nf), su(nu,0:mf),sv(nvh,0:nf),
     3   luv(nuvh), resmn(mnf), gval(ndim),
     4   xg(mnf,ndim), xh(mnf), stat=istat)
      if (istat .ne. 0) stop 'allocation error in min__xerr_phimn_nes'
      luv(:) = zero
      resmn(:) = zero
      gval(:) = zero
      xh(:) = zero
      smunv(:nuvh,0:mf,-nf:nf) = zero
      cmunv(:nuvh,0:mf,-nf:nf) = zero
      cu(:nu,0:mf) = zero
      cv(:nvh,0:nf) = zero
      su(:nu,0:mf) = zero
      sv(:nvh,0:nf) = zero
      xg(:mnf,:ndim) = zero

c  Calculate all the resonance factors 1/(m*iota+n*np) now, so
c  they can be used to turn dmn into xmn = dmn * resmn
       if (eps .lt. zero) then   !no resonance modification
          resmn(:mnf) = one
       else                    !standard resonance modification
          epsq = eps*eps
          resmn(1) = 1
          k = 0                !m=n=0 is never used here
          do n = 1, nf
             k = k + 1
            resmn(k) = real(n*np,rprec)/((n*np)**2 + epsq)
          enddo
          do m = 1, mf
             do n = -nf, nf
                k = k + 1
                res = (m*iota_edge + n*np)
                resmn(k) = res / (res**2 + epsq)
             enddo
          enddo
       endif

c............................................................
c  Calculate all the sin/cos factors needed on fourier grid mf,nf
c  for the theta, phi VMEC system, NOT the straight-line one
c  Note: For VMEC, m*u-n*np*v convention means np should be negative
c  These can be used on smaller lmn grids too, so saves time
      call trigfact(nu,nv,nvh,nuvh,mf,nf,cmunv,smunv,cu,su,cv,sv)

c  Following are all in VMEC co-ordinates, not the straight-line ones
c  so do them now before sin/cos(mu) are chaged later

c  Turn incoming lambda coeffs lmn into luv(theta,phi) with sine
c  in the theta, phi VMEC system, not the straight-line one
      mnl = 1+ns1 + ms1*(2*ns1 + 1)
      istat=0
      if (.not.allocated(lmn)) allocate (lmn(mnl), stat=istat)
      if (istat .ne. 0) stop 'allocation error in min_xerr_phimn'
      lmn(:) = zero
      k = 0
      do n = 0, ns1              !First the m = 0 coeffs
         k = k + 1
         lmn(k) = -cl1(0,-n)     !Because bnorm wrote only -n for m=0
      enddo                      !and we need only +n for m=0
      do m = 1, ms1              !Then the m > 0 ones
         do n = -ns1, ns1
            k = k + 1
            lmn(k) = cl1(m,n)    !here pack the cl1 array
         enddo
      enddo

      call fmn_to_uv(nu, nv, nuvh, luv, ms1, ns1, mnl, lmn,
     1     mf, nf, smunv)

c  Multiply ben, bfn by dsur/phip
      ben(:nuvh) = ben(:nuvh)*dsur(:nuvh)/phip_edge
      do i = 1, ndim
         bfn(:nuvh,i) = bfn(:nuvh,i)*(dsur(:nuvh)/phip_edge)
      end do
c............................................................
c  Note: Everything is calculated at (theta,phi) grid in VMEC
c
c  Integrate over u,v by doing xmn = xmn + f[ku,kv], i.e,
c  first calculate dmn from b dot ds
c  dmn = sum_ku_kv[sin(m[u+luv(kuv)]+n*np*v) * (bds[kuv])] *(du*dv/phip)
c  and then calculate xmn from dmn
c  xmn = dmn *(m*iota+n*np)/( (m*iota+n*np)^2 + eps^2 )
c
c  Note: You need to recaculate only sin/cos(m[u+luv(kuv)]), since
c        sin/cos(n*np*v) were already calculated above
c  Note: It is better to put the ku,kv loops outside the m,n loops
c        because at each u,v, u+lam(u,v) can be calculated once and
c        from it cos/sin(m*[u+lam(u,v)]) can be calculated by making
c        only one call to cos/sin(u+lam(u,v)) and then iterating
c  Note: sin( m[theta+lambda] + n*np*phi) and
c        cos( m[theta+lambda] + n*np*phi) will be calculated here.
c        They will replace the smunv, cmunv calculated earlier,
c         i.e., sin/cos( m*theta + n*np*phi)
c
      xg(:mnf,:ndim) = zero            !zero all fourier coeffs of modified
      xh(:mnf) = zero                  !Green functions
      du = pi2/nu                      !Step size in theta
      su(:nu,0) = zero                 !for m=0
      cu(:nu,0) = one
      i = 0                            !i is real space index 1,nuvh
c......
      do kv = 1, nvh                   !See surface_plas.f, this is the same
c
c       Check for centerline (v=1/2 or kv=1+nv/2)
c       use value of v because it tkaes care of both odd and even nv
         if ((kv.eq.1) .or. (2*(kv-1).eq.nv)) then
            nwt = 1                     !centerline gets only one weight
         else
            nwt = 2                   !all others carry twice the weight
         endif
c.......
         do ku = 1, nu
            i = i + 1                 !i goes from 1 to nuvh=nu*(1+nv/2)
            hval = nwt*ben(i)         !do this outside m,n loops to save time
            gval(:ndim) = nwt*bfn(i,:ndim)
c
c         Re-calculate sin/cos(m*u) for this u since it differs from
c          theta by u = theta + lambda: straight coordinate
c         Note: cos/sin(n*v) need not change, use precalculated value
c         Note: This calculation cannot be done outside the
c               v loop since lambda is a function of both u and v
c
c         Find the straight-line u at this theta,phi = theta+lambda:
            u = du*(ku - 1) + luv(i)             !u = theta + lambda
c
c         Make only one call per u,v to cos/sin and then iterate
c         calculate cos/sin(m*[theta+lambda]):
            su(ku,1) = sin(u)                    !sin/cos(theta+lambda)
            cu(ku,1) = cos(u)
            do m = 2, mf                     !Now calculate sin/cos(m*u)
               n = m - 1
               su(ku,m) = su(ku,1)*cu(ku,n) + cu(ku,1)*su(ku,n)
               cu(ku,m) = cu(ku,1)*cu(ku,n) - su(ku,1)*su(ku,n)
            end do
c
c         Now assemble everything at this (u,v) for all (m,n), i.e.,
c         1. Calculate sin(m*[u+luv(u,v)]+n*np*v)
c            at this u,v for all m,n in (mf,nf) fourier space
c         2. Then add the contribution to xmn from this u,v point
c
c         Note: Since luv(u,v) depends on both u and v, this cannot
c               be done outside the ku,kv loops, must be done here
c         Skip the m=0,n=0 case since 1/(m*iota+n) is singular, so
c         we set xmn(m=0,n=0) = 0, i.e., f.s.average of xerr = 0
c         Note: Only cu,cv, su,sv WILL be changed here, but
c               smunv and cmunv will NOT be changed. So from now on,
c               cu,cv, su,sv will be for the straight-line system, and
c               smunv and cmunv will be for theta, phi VMEC system
c         Note: - sign in summation since bmn's are derivs of cos
            j = 1                          !j = m,n space index 1 to mnf
            do n = 1, nf                         !first the m = 0 cases
               j = j + 1           !skips j=1, starts from 2 for m=0,n=1
c           for m=0, m*theta = m*u = 0 in both systems
c           Add to xg and xh integrals from this theta,phi point
               xh(j) = xh(j) + hval*sv(kv,n)
               xg(j,:ndim) = xg(j,:ndim) - gval(:ndim)*sv(kv,n)
            end do
c
            do m = 1, mf                         !Next the m>0 cases
               do n = -nf, nf
                  j = j + 1
c           Save cmunv and smunv in straight-line system (u,v)
c           for use in calculating xtp and dxdt, dxdp later.
                  nn = abs(n)     !Need this since su/cu are over n=0:nf
                  ns = sign(1,n)
                  cmunv(i,m,n)=cu(ku,m)*cv(kv,nn)-su(ku,m)*sv(kv,nn)*ns
                  smunv(i,m,n)=su(ku,m)*cv(kv,nn)+cu(ku,m)*sv(kv,nn)*ns
c            Add to dmn integral from this theta,phi point
                  xh(j) = xh(j) + hval*smunv(i,m,n)
                  xg(j,:ndim) = xg(j,:ndim) - gval(:ndim)*smunv(i,m,n)
               end do
            end do
c
         end do                                  !ku loop
      end do                                     !kv loop
c.......
c      Multiply xg, xh by resonance and overall normalization factor,
c      Note: the 1/(2*pi)^2 is obsent since u,v go over 0 to 1
c                                             !2 since only m>0 are used
      resmn(:mnf) = 2*resmn(:mnf)/(nu*nv)
      xh(:mnf) = xh(:mnf)*resmn(:mnf)
      do i = 1, ndim
         xg(:mnf,i) = resmn(:mnf)*xg(:mnf,i)
      end do

c............................................................
c....  Calculate phi(m,n) by SVD solving xg * phimn = -hmn
c      Note: negative sign in xh has already been absorbed
c      Note: this will be a square matrix unless lamdamn
c            had more coeffs than mf=nu/2 and nf=nv/2, where
c            nu,nv are points on plasma surface 1 period
c            Since lamda is a function defined olny on these
c            points, it should NOT have more fourier coeffs.
c            mnf should be = nuvh, but ndim can be smaller
       write(inesc, '(a)') 'mf, nf, mnf, nuvh, ndim ='
       write(inesc,*)  mf, nf, mnf, nuvh, ndim
       if(method.ne.0) call svd_solve(mnf, ndim, mnf, ndim,
     1                    xg, xh, svdv, svdw, nsvdw)

c............................................................
c      Calculate RMS Xerror scan over all weights kept.
c      svdv(:,i) = phi_i = ith solution (with i svd weights kept,
c      xg * phi_i - xh = error due to ith solution
c      Report its RMS value when summed over all mnf (m,n) modes
c      Use the no longer needed resmn(:mnf) array for temp storage
       write (inesc, '(a)') '---- Xerror table ----'
       write (inesc, '(a)')
     1'   isvd     weight(i),     rms Xerr(i),      abs Xerr(i)'
       do i = nsvdw, 1, -1
          resmn(:mnf)=Matmul(xg(:mnf,:ndim),svdv(:ndim,i))-xh(:mnf)
          rmsxerr=sqrt( sum(resmn(:mnf)*resmn(:mnf)) /mnf )
          absxerr=sum( abs(resmn(:mnf)) ) / mnf
          write(inesc,"(i5,3g25.15)") i, svdw(i), rmsxerr, absxerr
       enddo
       write (inesc, '(a)') '---- end Xerror table ----'

c............................................................
c      Restore ben, bfn by multiplying by phip/dsur
c      for later Berr calculation using phimn solutions found here
       ben(:nuvh) = phip_edge * ben(:nuvh) / dsur(:nuvh)
         do i = 1, ndim
            bfn(:nuvh,i) = phip_edge * bfn(:nuvh,i) / dsur(:nuvh)
         enddo

c............................................................
c      Deallocate everything
      if (iloop .le. 0) deallocate( smunv,cmunv, cu,cv, su,sv,
     1  luv, lmn, resmn, gval, xg, xh, stat=istat)

      end subroutine min_xerr_phimn
