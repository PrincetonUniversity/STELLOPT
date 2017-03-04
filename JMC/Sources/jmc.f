
      program jmc

      use stel_kinds
      use read_boozer_mod
      implicit none
c-----------------------------------------------
c   L o c a l   P a r a m e t e r s
c-----------------------------------------------
      real (kind=rprec) ::  hs, tflux, dmu0
      integer :: k, i, j, i_surf, numargs, m, n 
      integer :: i_alloc
      integer :: ns_booz_max, mb_retained, nb_retained, m0b, n0b,
     >           mnmax, ns, nfp 
      real (kind=rprec) :: twopi
      integer, dimension(:), allocatable :: ixm, ixn
      real (kind=rprec), dimension(:), allocatable ::
     >                  es, iota, pprime, curr_pol, curr_tor,
     >                  sqrtg00
      real (kind=rprec), dimension(:,:), allocatable ::
     >                  rmnc, zmnc, lmnc, bmnc
      character*(120) :: booz_input_file, arg1

      TWOPI = 8.*atan(1.)
      dmu0 = 2.0*twopi*1.0e-7
c
c     Read data from the boozmn file and allocate storage:
c
      call getcarg(1, arg1, numargs)
      booz_input_file = trim(arg1)
     
      k=0
      call read_boozer_file(trim(booz_input_file),k)
      if (k .ne. 0) stop 'Error reading boozmn file'
c
c     Read control file
c
      open(unit=5, file="mkjmc_param.in",status="unknown",
     >                                   iostat = k)
      if (k .ne. 0) stop "Error when open mkjmc_param.in"
      read(5,*) ns_booz_max
      read(5,*) mb_retained
      read(5,*) nb_retained
c     if( ns_booz_max .ne. ns_b-2 ) stop "ns \= ns_booz_max +2"
c
      m0b = mboz_b-1
      n0b = nboz_b
      mnmax = mnboz_b
      ns = ns_b
      nfp = nfp_b
      tflux = abs( phi_b(ns_b) )

!     truncate modes
!     assume m=0, mboz
!            n=-nboz, nboz
!     discard m> mb_retained
!             iabs(n)> nb_retained
!     perhaps better to retain largest mnmax modes
!             for a reference surfaces

      if (mnboz_b .ne. (mboz_b-1)*(2*nboz_b+1)+nboz_b+1 )
     >      stop "boozer modes inconsistent"

      k = 0
      if (mb_retained .le. 0) mb_retained = m0b
      if (nb_retained .le. 0) nb_retained = n0b 
      do j = 1, mnboz_b
      if( ixm_b(j) .gt. mb_retained ) cycle
      if( iabs(ixn_b(j)/nfp_b) .gt. nb_retained ) cycle
      k = k + 1
      enddo

      mnmax = k

! **********************************************************************
! Allocate storage arrays
! **********************************************************************

      ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
      IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

      ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),                
     >     pprime(ns), sqrtg00(ns), stat = i_alloc)
      IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

      ALLOCATE(rmnc(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),        
     >     bmnc(ns,mnmax),                                             
     >     stat = i_alloc)
      IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
! **********************************************************************


      k = 0
      do j = 1, mnboz_b
      if( ixm_b(j) .gt. mb_retained ) cycle
      if( iabs(ixn_b(j)/nfp_b) .gt. nb_retained ) cycle
      k = k + 1
      ixm(k) = ixm_b(j)
      ixn(k) = ixn_b(j)
      enddo

      hs = 1./float( ns_b -1 )
      do i = 3, ns

       k = 0
       do j = 1, mnboz_b

       if( ixm_b(j) .gt. mb_retained ) cycle
       if( iabs(ixn_b(j)/nfp_b) .gt. nb_retained ) cycle
       k = k + 1

        rmnc(i,k) = rmnc_b(j,i)
        zmnc(i,k) = zmns_b(j,i)
c
c       note: zeta_b - zeta from booz_xform, 

        lmnc(i,k) = -pmns_b(j,i)*nfp_b/twopi
        bmnc(i,k) = bmnc_b(j,i)

           if (ixm_b(j).eq.0 .and. ixn_b(j).eq.0)            
     >         sqrtg00(i) = gmnc_b(j,i)

       enddo

        es(i) = (float(i)-1.5)*hs
        iota(i) = iota_b(i)
c
c       note sqrtg00 are filled with zeros
c
        if( i .lt. ns ) then
        pprime(i) = (pres_b(i+1) - pres_b(i-1))/(2.*hs)
        else
        pprime(ns) = 2.0*pprime(ns-1)-pprime(ns-2)
        endif

c       sqrtg00(i) = 0.0

        curr_pol(i) = bvco_b(i)*twopi/float(nfp_b)
        curr_tor(i) = buco_b(i)*twopi
c
c       new F and I function for currents

c       curr_pol(i) = bvco_b(i)
c       curr_tor(i) = buco_b(i)

      enddo

      call read_boozer_deallocate
c
c     write jmc Bmn file
c
      m0b = 0
      n0b = 0
      do k = 1, mnmax
      if( ixm(k) .gt. m0b ) m0b = ixm(k)
      ixn(k) = -ixn(k)/nfp
      if( iabs(ixn(k)) .gt. n0b ) n0b = iabs(ixn(k))
      enddo

      open(unit=9,file="Bmns."//trim(arg1),status="unknown")
      write(9,*)     ' m0b  n0b  nsurfaces  nper flux'
      write(9,'(4i5,1x,f10.5)') m0b, n0b, ns-2, nfp, tflux

c      /* Output the mapping */

       do i = 3, ns
         write(9,*) '    s          iota        curr_pol     curr_tor ',
     .              '    pprime       sqrt g(0,0)'
         write(9,'(1p6e12.4)') es(i),iota(i)/nfp,curr_pol(i),
     .                         curr_tor(i),pprime(i)*dmu0,sqrtg00(i)
         write(9,*)
     .   '   m     n        r          z    (phib-phi)*nper/twopi',
     .   '    bmn'
       do k = 1, mnmax
       m = ixm(k)
       n = ixn(k)
       if( m .eq. 0 .and. n .ne. 0 ) then
         write(9,'(2i5,1p4e16.8)') m, n, rmnc(i,k)/2., zmnc(i,k)/2.,
     .                             lmnc(i,k)/2., bmnc(i,k)/2.
         write(9,'(2i5,1p4e16.8)') m, -n, rmnc(i,k)/2., -zmnc(i,k)/2.,
     .                             -lmnc(i,k)/2., bmnc(i,k)/2.
         else
         write(9,'(2i5,1p4e16.8)') m, n, rmnc(i,k), zmnc(i,k),
     .                             lmnc(i,k), bmnc(i,k)
       endif
       enddo
       enddo

      close(9)
      stop
      end

        subroutine getcarg(iseq, arg, numargs)
        integer:: iseq, numchars, numargs
        character*(*) :: arg

        if(iseq .eq. 1) numargs = iargc()
        numchars = getarg(iseq, arg)

        return
        end
