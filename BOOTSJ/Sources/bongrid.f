
      subroutine bongrid(irho, ians, ihere)

c  First time through, form the THETA-ZETAH grid
c  Every time through, evaluate B and QSAFETY on the grid for this RHO.
c  Here, m is the poloidal mode number and n (and nh, for n-hat) is the
c  toroidal mode number.

C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer irho, ians, ihere
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j, nh, m, mbuse1, imax, jmax
!      integer, save :: ihere = 0
      integer :: ij_max(2)
      real(rprec), save :: twopi, dth, dzetah, dthm, dzetahm
      real(rprec) :: d, dbmintok, sbmaxtok,
     1 b, del
C-----------------------------------------------

      if (ihere .eq. 0) then

c  allocate all of the variables needed on the theta-zeta grid

         call allocate_angles

c-----------------------------------------------------------------------
c  Form the THETA-ZETAH grid.

         twopi = 8*atan(one)
         dth = twopi/nthetah
         dzetah = twopi/nzetah

         do i = 1, nthetah
            theta(i) = (i - 1)*dth
         end do
         do j = 1, nzetah
            zetah(j) = (j - 1)*dzetah
         end do

!   Form a finer mesh to evaluate fine grained bmax.

         dthm = dth/(nthetahm-1)
         dzetahm = dzetah/(nzetahm-1)

!  load sin and cosine arrays

         do j = 1, nzetah
            do nh = 0, nbuse
               sinnj(nh,j) = sin(zetasign*nh*zetah(j))
               cosnj(nh,j) = cos(zetasign*nh*zetah(j))
            enddo
         enddo
         do i = 1, nthetah
            do m = -mbuse, mbuse
               sinmi(m,i) = sin(m*theta(i))
               cosmi(m,i) = cos(m*theta(i))
            end do
         end do
         ihere = 1
      endif          !end of ihere initial evaluation

c  Now evaluate B on the theta-zetah grid for this RHO using the epsilons just
c  found.  Loop over theta and phihat, summing up all the (m,nh) terms in B.

      do j = 1, nzetah
         do i = 1, nthetah
            b = zero
            do m = -mbuse, mbuse
               do nh = 0, nbuse
                  b = b + amnfit(irho,m,nh)*
     1               (cosmi(m,i)*cosnj(nh,j)-sinmi(m,i)*sinnj(nh,j))
                  IF (lasym_bootsj) b = b + amnfit2(irho,m,nh)*
     1               (sinmi(m,i)*cosnj(nh,j)+cosmi(m,i)*sinnj(nh,j))
               end do
            end do
            bfield(i,j) = abs(b)
         end do
      end do

c   find max of b on global mesh

      ij_max = maxloc(bfield)
      imax = ij_max(1)
      jmax = ij_max(2)

c  use the theta and zeta from this search as the center of a finer search

c  first form the grid

      thetam(1) = theta(imax) - dth/2
      zetahm(1) = zetah(jmax) - dzetah/2

      do i = 2, nthetahm
         thetam(i) = thetam(i-1) + dthm
      enddo
      do j = 2, nzetahm
         zetahm(j) = zetahm(j-1) + dzetahm
      enddo

c  load the sines and cosines on the finer mesh

      do j = 1, nzetahm
         do nh = 0, nbuse
            sinnjm(nh,j) = sin(zetasign*nh*zetahm(j))
            cosnjm(nh,j) = cos(zetasign*nh*zetahm(j))
         enddo
      enddo
      do i = 1, nthetahm
         do m = -mbuse, mbuse
            sinmim(m,i) = sin(m*thetam(i))
            cosmim(m,i) = cos(m*thetam(i))
         end do
      end do

c  evaluate b on the finer mesh

      do j = 1, nzetahm
         do i = 1, nthetahm
            b = zero
            do m = -mbuse, mbuse
               do nh = 0, nbuse
                  b = b + amnfit(irho,m,nh)*
     1              (cosmim(m,i)*cosnjm(nh,j)-sinmim(m,i)*sinnjm(nh,j))
                  IF (lasym_bootsj) b = b + amnfit2(irho,m,nh)*
     1               (sinmim(m,i)*cosnjm(nh,j)+cosmim(m,i)*sinnjm(nh,j))
               end do
            end do
            bfieldm(i,j) = abs(b)
         end do
      end do


c- evaluate bmax1(irho), thetamax(irho), b2avg(irho), and zetahmax(irho)
c  based on finer mesh evaluation

      ij_max = maxloc(bfieldm)
      imax = ij_max(1)
      jmax = ij_max(2)
      thetamax(irho) = thetam(imax)
      zetahmax(irho) = zetahm(jmax)
      bmax1(irho) = bfieldm(imax,jmax)

c  evaluate jacobian.  Leave off flux surface quantites.  Evaluate
c  the sum for later use.  Note that the value is not scaled to
c  Bmax.

      gsqrt_b(:nthetah,:nzetah) = one/bfield(:nthetah,:nzetah)**2
      sum_gsqrt_b = sum(gsqrt_b)

c  find b2avg LAB--boozer paper uses both bzero2 (1/leading coefficient of 1/b**2 expansion
c  and <b**2>.  They are the same (in boozer coordinates).  Both result from
c  flux surface averages. I will use <b2> everywhere

      b2avg(irho) = sum(bfield**2 * gsqrt_b)/sum_gsqrt_b

c  Scale the array BFIELD so that it contains B/Bmax instead of B.

      bfield(:nthetah,:nzetah) = bfield(:nthetah,:nzetah)/bmax1(irho)
      where(bfield .gt. one) bfield = one
      b2obm = bfield**2

c   pressure related derivatives are needed
c   first calculate difference denominators for pressure related derivatives
c   difference array has differences on the full mesh

c   use a parabolic fit near rho = 0 to get derivatives at 1st half mesh

      if(irho .eq. 1) then
         drho = (rhoar(2)**2 - rhoar(1)**2)/(2*rhoar(1))

c   use slope of last two points at outer rho point

      elseif(irho .eq. irup) then
         drho = 0.5_dp*(d_rho(irho)+d_rho(irho-1))

c  all other points

      else
         drho = d_rho(irho) + 0.5_dp*(d_rho(irho+1)+d_rho(irho-1))
      endif

c  evaluate Electron temperature gradients in Kev

      if (irho .ne. 1 .and. irho .ne. irup) temperho1 =
     1   (tempe1(irho+1)-tempe1(irho-1))/drho
      if (irho .eq. 1) temperho1 = (tempe1(irho+1)-tempe1(irho))/drho
      if (irho .eq. irup) temperho1=(tempe1(irho)-tempe1(irho-1))/drho

c evaluate Ion temperature gradients in Kev

      if (irho .ne. 1 .and. irho .ne. irup) tempirho1 =
     1  (tempi1(irho+1)-tempi1(irho-1))/drho
      if (irho .eq. 1)tempirho1 = (tempi1(irho+1)-tempi1(irho))/drho
      if (irho .eq. irup)temperho1=(tempi1(irho)-tempi1(irho-1))/drho

c  evaluate electron density gradients in 10**20 m-3

      if (irho .ne. 1 .and. irho .ne. irup) densrho1 =
     1  (dense(irho+1)-dense(irho-1))/drho
      if (irho .eq. 1) densrho1 = (dense(irho+1)-dense(irho))/drho
      if (irho .eq. irup) densrho1 = (dense(irho)-dense(irho-1))/drho

c  write out the lower order Boozer coefficients

      mbuse1 = mbuse
      if (mbuse > 5) mbuse1 = 5            !mbuse1 > 5 will not fit on page

      write (ians, 400)
  400 format(/' nh ',$)
      do m = -mbuse1, mbuse1
         write (ians, 402) m
      end do
  402 format('   m=',i2,'   ',$)
      write (ians, '(a1)') ' '
c
      do nh = 0, nbuse
         write (ians, 406) nh, (amnfit(irho,m,nh),m=(-mbuse1),mbuse1)
      end do
  406 format(1x,i2,1p13e10.3)

c  Calculate the fraction trapped and fraction passing for the "equivalent"

      call tok_fraction(fttok(irho),irho)

      fptok(irho) = one - fttok(irho)

      end subroutine bongrid
