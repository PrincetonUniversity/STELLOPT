
! *************************************************************************
!  LIMSEP
!
!  This subroutine computes minimum separation distance (signed) between
!  the plasma and a set of limiting surfaces, described as an array of
!  piecewise linear curves at specified toroidal angles.
!
!  This subroutine is modified from SURFSEP, which does a similar calculation
!  to a bounding surface specified by Fourier Coefficients.
!
!  M. Zarnstorff   July 2002
! *************************************************************************

      subroutine limsep(mnmax_pl, r_pl, z_pl, xm_pl, xn_pl, nphi_lim,
     1    phi_lim, nrz_max, r_lim, z_lim, nper, distance, ierr)
      use stel_kinds
!      use cmnf1
      use cmnf2
!      use geom

      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer mnmax_pl, nphi_lim, nrz_max, nper, ierr
      real(rprec) :: distance
      real(rprec), dimension(mnmax_pl) :: r_pl, z_pl, xm_pl, xn_pl
      real(rprec), dimension(nphi_lim) :: phi_lim
      real(rprec), dimension(nrz_max, nphi_lim) :: r_lim, z_lim
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, parameter :: igrid = 64

      integer :: i, ip, it, il, iu, np
      real(rprec) :: d, du, u2, v2, u2_new, d_min
      real(rprec), dimension(igrid) :: r2, z2, dist

      integer :: nlim, it_min, ip_min, il_min
      integer, dimension(1) :: iresult, iresult2

!-----------------------------------------------
!
      ierr = 0                                   ! presume success to start

!
!     set np to igrid
      np = igrid

!
!     compute the total number of modes for each surface
      nmn2 = mnmax_pl
!
!     equate the common block variables with the arguments
!
      allocate (xm2(nmn2), xn2(nmn2), rmn2(nmn2), zmn2(nmn2),
     1          stat=ierr)
      if (ierr .ne. 0) return
!
!
!     equate the common block variables with the arguments
!     Surface 2 - Plasma
!
      xm2(:nmn2) = xm_pl(:nmn2)
      xn2(:nmn2) = xn_pl(:nmn2)
      rmn2(:nmn2) = r_pl(:nmn2)
      zmn2(:nmn2) = z_pl(:nmn2)
!
!     write (6,*) 'done copying data', nmn2,nmn1
!
      distance = 1.0e30_dp

      du = 1.0_dp/np

      d_min = 1.0e30_dp

      it_min = -1
      ip_min = -1
      il_min = -1

!
!  grided search to find nearest segment and plasma location
!
      do it = 1, nphi_lim
         if( phi_lim(it) > 360 .or. phi_lim(it) < -360 ) exit
         if( r_lim(1,it) == 0 .or. r_lim(2,it) == 0) exit

         do ip = 3, nrz_max
            if( r_lim(ip, it) == 0 ) then
                nlim = ip - 1
                exit
            endif
         enddo

         v2 = modulo((nper*phi_lim(it))/360, 1._rprec)
         

!
!   get (r,z) points on plasma surface
!
         do ip = 1, np
            u2 = (ip - 1)*du
            call dmnf2 (u2, v2, r2(ip), z2(ip))
!            write(*,*) '|ip=',ip,'| u2=',u2,'| v2=',v2,'| r2(ip)=', 
!     1                  r2(ip),'| z2(ip)=',z2(ip)
         enddo

         do il = 2, nlim
            do ip = 1, np
               call lim_dist( r2(ip), z2(ip), r_lim(il-1, it),
     1                        z_lim(il-1, it), dist(ip))
            enddo
            iresult = minloc(dist)

            if( minval(dist) < d_min) then
               d_min = minval(dist)
               it_min = it
               iresult = minloc(dist)
               ip_min = iresult(1)
               il_min = il
            endif
         enddo
!         write(*,*) '| il=',il,'| d_min=',d_min,'| it_min=',it_min,
!     1              '| iresult=',iresult,'| ip_min=',ip_min,
!     2              '| il_min=',il_min
      enddo
!      write(*,*) 'Miniminum Distance at: ','phi= ',phi_lim(it_min),
!     1           '   u= ',(ip_min-1)*du,'   r2= ',r2(ip_min),
!     2           '   z2=',z2(ip_min),'   r_lim= ',r_lim(il_min,it_min),
!     3           '   z_lim=',z_lim(il_min,it_min)

!
!  if no valid limiter segments specified, we are done
!
      if( it == -1) then
         distance = 0
         ierr = 1
         return
      endif

!
!     refine the point of minimum distance from the plasma
!     to the limiter segment by binary search
!
!      print *,'min found:',phi_lim(it_min),d_min,il_min
      v2 = modulo((nper*phi_lim(it_min))/360, 1._rprec)

      u2 = (ip_min - 1)*du
      u2_new = u2

      do ip = 1, 10
         du = du/2

         do iu = -1, 1, 2
            call dmnf2( u2+du*iu, v2, r2(1), z2(1))

            call lim_dist( r2(1), z2(1), r_lim(il_min-1, it_min),
     1                        z_lim(il_min-1, it_min), dist(1))

            if( dist(1) < d_min ) then
               u2_new = u2 + du*iu
               d_min = dist(1)
            endif
         enddo

         u2 = u2_new
      enddo

      distance = d_min

      deallocate(xm2, xn2, rmn2, zmn2)

      end subroutine limsep
