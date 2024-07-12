
! ----------------------------------------------------------------------
      subroutine fourier(ndim)
! ----------------------------------------------------------------------
c                                                             01.01.89
c     purpose:
c
c
C-----------------------------------------------
      use Vmeshes
      use NumParams
      USE Vfourier2
      USE Vsolver1, ONLY: a, work
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ndim
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nel, i, ni, nim, k, inc, jump, nuf, lot
C-----------------------------------------------
c
c   1. fourier transformation with respect to u
c
      nel = nv*ndim
      do i = nel, 1, -1
         ni = i*(nu + 2)
         nim = i*nu + 1
         a(ni)   = zero
         a(ni-1) = zero
         do k = 1, nu
            a(ni-1-k) = a(nim-k)
         end do
      end do
c
      inc = 1
      jump = nu + 2
      call fft991 (a, work, trigs, ifax, inc, jump, nu, nel, -1)
c
      nuf = 2*(mf + 1)
      do i = 1, nel
         ni = (i - 1)*(nu + 2)
         nim = (i - 1)*nuf
         do k = 1, nuf
            a(nim+k) = a(ni+k)
         end do
      end do
c
c   2.   fourier transformation with respect to v
c
      inc = mf + 1
      jump = inc*nv
      lot = ndim
      do i = 1, mf + 1
         call cfft99(a(2*i-1),work,trigsv,ifaxv,inc,jump,nv,lot,-1)
      end do

      end subroutine fourier
