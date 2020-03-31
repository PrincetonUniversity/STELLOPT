
!-------------------------------------------------------------
      subroutine ab_matrix(ab, bfn, dumv, ndim, mdspwa, curwt)
!................................................................
c     Purpose:
c     Calculate ab matrix:
c     use MDSPWA  to control power of dsur,
c         CURWT   to control amount of J minimization
c................................................................
C   M o d u l e s
C-----------------------------------------------
      use Vmeshes
      use NumParams
      USE Vvacuum3
      USE Vprecal1
      USE Vprecal3
      USE Vprecal6, ONLY: comu, simu
      USE Vprecal7, ONLY: conv, sinv
      USE Vsurfac14, ONLY: ben
      USE Vsurfac13, ONLY: dsur1
      USE Vsurface2, ONLY: xu, yu, xv, ru, rv
      USE Vsurface3, ONLY: yv, zu, zv
      USE Vsurface7, ONLY: snx, sny, snz
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ndim, mdspwa
      real(rprec) :: ab(ndim, ndim+1), bfn(nuvh1, ndim),
     1   dumv(nuvh1), curwt
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer::mexp,l,lp,i,j,n,m,ndimp,jv,ju,m1,n1,m2,n2
      integer :: icm(ndim), icn(ndim)
      real(rprec) :: curwt2, cupnp, eux, euy, euz, evx, evy, evz,
     1   cmunv1, f1x, f1y, f1z, ex, ey, ez, cmunv2, f2x, f2y, f2z
C-----------------------------------------------

      ab(:ndim,:ndim+1) = zero

c     Calculate the dsur factor for integral over plasma surface
c     Also include the muliplicity factor(i) (1 or 2) in dumv(i)
      if (mdspwa >= 2) then
         mexp = mdspwa - 2                       !exponent in dsur**mexp
         dumv(:nuvh1) = factor(:nuvh1)*dsur1(:nuvh1)**mexp
      else
         dumv(:nuvh1) = one             !here just sum|berr|^2 minimized
      endif

c...1.First calculate the G and H parts for standard nescoil
c     Calculate ab matrix G part= Tr(G)*G for LSQ case
c     by integrating over nuvh1 on plasma surface
c     This makes the matrix square (ndim,ndim)
      do l = 1, ndim
         do lp = l, ndim
          ab(l,lp) = ab(l,lp) +
     1     sum( bfn(:nuvh1,l) *bfn(:nuvh1,lp) *dumv(:nuvh1) )
         end do
      end do

c     Calculate ab matrix h part (last column(ab)=Tr(G)*bnorm)
c     by integrating over nuvh1
      do l = 1, ndim
         ab(l,ndim+1) = ab(l,ndim+1) -
     1     sum( bfn(:nuvh1,l) *ben(:nuvh1) *dumv(:nuvh1) )
      end do

c...2.Add jminimization parts Fi and Ei to ab matrix
      if (curwt > 0) then
         write (inesc, '(a,1pe16.8)')
     1     'Jmin targeted Using SVD Weight = ',curwt
         curwt2 = 2*curwt

c        Calculate the dsur factor for integral over coil surface
c        Also include the factor(i) (1 or 2) and curwt here into dumv(i)
c        Also fold in the 1/|eu x ev| = 1/dsur into dumv(i)
         if (mdspwa >= 2) then
            mexp = (-1) + .5_dp*(mdspwa - 2)    !divide by 2 to avoid sqrt

            do i = 1, nu
               dumv(i) = curwt*
     1            (snx(i)*snx(i)+sny(i)*sny(i)+snz(i)*snz(i))**mexp
               j = i + nuv/2
               dumv(j) = curwt*
     1            (snx(j)*snx(j)+sny(j)*sny(j)+snz(j)*snz(j))**mexp
            end do

            do i = nu+1, nuv/2
               dumv(i)= CURWT2 *
     1         (snx(i)*snx(i)+sny(i)*sny(i)+snz(i)*snz(i))**mexp
            enddo

         else
            dumv(:nuvh1) = one       !use just sum||^2 , not integral ds
         endif

c        Fill index arrays to get m,n from imn
c        ndim -> im,in map in fourier space
         i = 0
         do n = 1, nf                        !first m=0 cases (nf total)
            i = i + 1
            icm(i) = 0
            icn(i) = n
         end do
         do m = 1, mf             !next m>0 cases  ( mf*(2*nf+1) total )
            do n = -nf, nf
               i = i + 1
               icm(i) = m
               icn(i) = n
            end do
         end do

c        Add F and E to the G and H parts: Do not store F or E matrices
         cupnp = cup/np
         ndimp = ndim + 1
         i = 0
         do jv = 1, 1 + nv/2                     !half v + centerline
            do ju = 1, nu          !all of u. Total = nu*(1+nv/2) = nuvh
               i = i + 1
               eux = xu(i)                       !Tanget vector eu
               euy = yu(i)
               euz = zu(i)
               evx = xv(i)                       !Tanget vector ev
               evy = yv(i)
               evz = zv(i)

               do l = 1, ndim               !go over all ndim m,n values
                  m1 = icm(l)
                  n1 = iabs(icn(l))
                  cmunv1 = pi2*(comu(ju,m1)*conv(jv,n1)
     1                   -      simu(ju,m1)*sinv(jv,n1))
                  f1x = n1*eux - m1*evx
                  f1y = n1*euy - m1*evy
                  f1z = n1*euz - m1*evz

                  ex = cut*evx - cupnp*eux
                  ey = cut*evy - cupnp*euy
                  ez = cut*evz - cupnp*euz

c                 Add the F * E parts (x,y,z) to ab matrix:
                  ab(l,ndimp) = ab(l,ndimp) + dumv(i)*cmunv1*
     1               (f1x*ex + f1y*ey + f1z*ez)

                  do lp = l, ndim         !go over all ndim m',n' values
                     m2 = icm(lp)
                     n2 = iabs(icn(lp))
                     cmunv2 = pi2*
     1                 (comu(ju,m2)*conv(jv,n2)-simu(ju,m2)*sinv(jv,n2))
                     f2x = n2*eux - m2*evx
                     f2y = n2*euy - m2*evy
                     f2z = n2*euz - m2*evz

c                    Add the F * F parts (x,y,z) to ab matrix:
                     ab(l,lp) = ab(l,lp) + dumv(i)*cmunv1*cmunv2*
     1                  (f1x*f2x + f1y*f2y + f1z*f2z)
                  end do
               end do
            end do
         end do

      endif                       !End addition of jminimization section

c...3.Finally, Symmetrize ab matrix G part
c     no need to do center spine, so ndim-1
      do lp = 1, ndim - 1
         ab(1+lp:ndim,lp) = ab(lp,1+lp:ndim)
      end do

      end subroutine ab_matrix
