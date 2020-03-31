
! ----------------------------------------------------------------------
      subroutine read_nescoil_input
!................................................................
c                                                             01/01/89
c     purpose:
c     Read input file with ALL nescoil settings
c     Write ALL inputs into output file so all settings will be with
c     the answers for later use, and also the output file can be used
c     as the input file to recreate an old run later.
c
c     Some of these dimensions used to be in meshes file in old version.
c     Now they are all dynamically allocated arrays
c................................................................
      USE stel_kinds
      use safe_open_mod
      use SvdCtrl, ONLY: mstrt, mstep, mkeep, mdspw, curwt, trgwt
      use OutCtrl
      use LoopCtrl
      use Vmeshes
!      USE Vprecal1, ONLY: np
      USE Vvacuum1, ONLY: cr, cz, ms, ns
      USE Vvacuum2, ONLY: ms1, ns1, cr1, cz1, cl1,
     1    iota_edge, phip_edge, curpol
      USE Vvacuum3
      USE Vprecal1
      USE Vvacuum4, ONLY: cr2, cz2, ms2, ns2
      USE Vvacuum5, ONLY: cr3, cz3, ms3, ns3
      USE Vvacuum6, ONLY: cf, sf, msf, nsf
      use Vbnorm, ONLY : extension
      use NumParams
      implicit none
c................................................................
C   L o c a l   V a r i a b l e s
c................................................................
      integer, parameter :: nescoil0 = 7
      integer :: m, n, ntotal, k, mr, nr, istat, numargs, iunit
      real(rprec) :: crin, czin, clin, crsin, czcin, clcin
      character*200 :: arg1
      logical :: lasym

      lasym = .false.
      clcin = 0
      czcin = 0
      crsin = 0
c ----------------------------------------------------------------------
c read data from command line, xnescoil nescin.extension
C-----------------------------------------------
      call getcarg(1, arg1, numargs)
      if (numargs .lt.1 .or. arg1 .eq. '-h' .or. arg1 .eq. '/h') then
         print *,
     1   ' Syntax:  xnescoil nescin.extension'
         stop 'Invalid command line'
      endif

      numargs = index(arg1,'nescin.')
      if (numargs .le. 0) then
         extension = trim(arg1)
      else
         extension = arg1(numargs+7:len_trim(arg1))
      end if

      iunit = nescoil0
      call safe_open(iunit, istat, trim(arg1), 'old', 'formatted')
      if (istat .ne. 0) then
         print *,' Type xnescoil -h for proper syntax'
         stop 'Error opening input file in nescoil'
      endif

      call safe_open(inesc, istat, 'nescout.'//extension,
     1   'unknown', 'formatted')
      if (istat .ne. 0) stop 'Error opening nescout file'

c................................................................
c read spatial dimensions
      read (iunit, *, err=99)
      read (iunit, *, err=99)
      read (iunit, *, err=99) nu, nv, nu1, nv1, npol, ntor
      read (iunit, *, err=99)
      write(inesc, 10) '----- Grid Spatial Dimensions -----'
      write(inesc, 10) 'nu, nv, nu1, nv1, npol, ntor'
      write(inesc,"(6i6,l)")  nu, nv, nu1, nv1, npol, ntor, lasym
      
      if (lasym) stop "Non-Stellarator Symmetry not supported!"

c................................................................
c read fourier dimensions
      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99) mf, nf, md, nd
      read (iunit, *, err=99)
      write (inesc, 10) '----- Fourier Dimensions -----'
      write (inesc, 10) 'mf, nf, md, nd'
      write (inesc,"(4i6)")  mf, nf, md, nd

      nmax  = 3*nv*(nu+2)+1
      mnd   = (md + 1)*(2*nd + 1)
      nuv   = nu*nv
      nuv1  = nu1*nv1
      nuvh  = nuv/2 + nu
      nuvh1 = nuv1/2 + nu1

c................................................................
c read plasma parameters from vmec
      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99) np, iota_edge, phip_edge, curpol
      read (iunit, *, err=99)
      write (inesc, 10) '----- Plasma information from VMEC -----'
      write (inesc, 10) 'np, iota_edge, phip_edge, curpol'
      write (inesc,"(i6,3g25.16)")  np, iota_edge, phip_edge, curpol

c................................................................
c read currents and nescoil controls

      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99) cut, cup, ibex
      read (iunit, *, err=99)
      write (inesc, 10) '----- Current Controls -----'
      write (inesc, 10) 'cut, cup, ibex'
      write (inesc,"(2g25.16,i6)")  cut, cup, ibex

c................................................................
c read svd and resonance control switches
      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99) mstrt,mstep,mkeep,mdspw,curwt,trgwt
      read (iunit, *, err=99)
      write (inesc, 10) '----- SVD controls -----'
      write (inesc, 10) 'mstrt, mstep, mkeep, mdspw, curwt, trgwt'
      write (inesc,"(4i6,2g25.16)")  mstrt,mstep,mkeep,mdspw,curwt,trgwt

c................................................................
c read detailed output writing control switches
      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99)
     1    w_psurf, w_csurf, w_bnuv, w_jsurf, w_xerr, w_svd
      read (iunit, *, err=99)
      write (inesc, 10) '----- Output controls -----'
      write (inesc, 10)
     1  'w_psurf, w_csurf, w_bnuv, w_jsurf, w_xerr, w_svd'
      write (inesc,"(6i6)") w_psurf,w_csurf,w_bnuv,w_jsurf,w_xerr,w_svd

c................................................................
c Allocate surface info arrays
      if (.not.allocated(cr))
     1 allocate (cr(0:md,-nd:nd), cz(0:md,-nd:nd), cr1(0:md,-nd:nd),
     1          cz1(0:md,-nd:nd), cr2(0:md,-nd:nd), cz2(0:md,-nd:nd),
     2          cr3(0:md,-nd:nd), cz3(0:md,-nd:nd), cf(0:md,-nd:nd),
     3          sf(0:md,-nd:nd), cl1(0:md,-nd:nd))
      cr1(0:md,(-nd):nd) = zero
      cz1(0:md,(-nd):nd) = zero
      cl1(0:md,(-nd):nd) = zero

c................................................................
c   read  plasma boundary data
      ms1 = 0
      ns1 = 0
      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99) ntotal
      write (inesc, 10) '----- Plasma Surface -----'
      write (inesc, 10) 'Number of fourier modes in table'
      write (inesc, *)  ntotal

      read (iunit, *,err=99)
      read (iunit, *,err=99)
      write (inesc, 10)
     1   '----- Plasma boundary fourier coefficients  -----'
      write (inesc, 10)
     1   '   m    n        R(m,n)     Z(m,n)    Lamda(m,n)'

      do k = 1, ntotal
         read (iunit, *,err=99) mr, nr, crin, czin, clin,
     1                                 crsin,czcin,clcin
         if (mr.lt.0 .or. mr.gt.md .or. nr.lt.(-nd) .or. nr.gt.nd) then
            write (inesc, 10) 'error(s) in input plasma fourier coeffs:'
            write (inesc, '(2(a,i4))') 'Found m = ', mr, ' and n = ', nr
            if (mr.lt.0) then
               write (inesc, '(a,i4,a)') 'Cannot have m = ', mr, ' < 0'
            endif
            if (mr.gt.md) then
               write (inesc, 20) 'm > md = ', md
               write (inesc, 20) 'change md to ', mr
            else
               if (nr.lt.(-nd)) write (inesc, 20) 'n < -nd = ', -nd
               if (nr.gt.nd) write (inesc, 20) 'n > nd = ', nd
               write (inesc, 20) 'change nd to ', abs(nr)
            endif
            write (inesc, 10) 'Correct nescoil input file and rerun'
            stop 'error in nescoil input file plas rzmn'
         endif
         write (inesc,"(2i4,6g20.10)") mr, nr, crin, czin, clin,
     1                                 crsin,czcin,clcin
         cr1(mr,nr) = crin
         cz1(mr,nr) = czin
         cl1(mr,nr) = clin
         ms1 = max(ms1,abs(mr))
         ns1 = max(ns1,abs(nr))
      end do

      read (iunit, *, err=99)

c................................................................
c read coil current surface data
c
      cr(0:md,(-nd):nd) = zero
      cz(0:md,(-nd):nd) = zero
c
      ms = 0
      ns = 0
      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99) ntotal
      write (inesc, 10) '----- Coil Surface -----'
      write (inesc, 10) 'Number of fourier modes in table'
      write (inesc,*)  ntotal

      read (iunit, *,err=99)
      read (iunit, *,err=99)
      write (inesc, 10) '----- Coil surface fourier coefficients -----'
      write (inesc, 10) '    m    n         R(m,n)         Z(m,n)'

      do k = 1, ntotal
         read (iunit, *,err=99) mr, nr, crin, czin, crsin, czcin
         if (mr.lt.0 .or. mr.gt.md .or. nr.lt.(-nd) .or. nr.gt.nd) then
            write (inesc, 10) 'error(s) in input plasma fourier coeffs:'
            write (inesc, '(2(a,i4))') 'Found m = ', mr, ' and n = ', nr
            if (mr.lt.0) then
               write (inesc,'(a,i4,a)') 'Cannot have m = ', mr, ' < 0'
            endif
            if (mr.gt.md) then
               write (inesc, 20) 'm > md = ', md
               write (inesc, 20) 'change md to ', mr
            else
               if (nr.lt.(-nd)) write (inesc, 20) 'n < -nd = ', -nd
               if (nr.gt.nd) write (inesc, 20) 'n > nd = ', nd
               write (inesc, 20) 'change nd to ', abs(nr)
            endif
            write (inesc, 10) 'Correct nescoil input file and rerun'
            stop 'error in nescoil input file coil surface rzmn'
         endif
         write (inesc,"(2i4,4g20.10)") mr, nr, crin, czin, crsin, czcin
         cr(mr,nr) = crin
         cz(mr,nr) = czin
         ms = max(ms,abs(mr))
         ns = max(ns,abs(nr))
      end do

c................................................................
c     Take care of n=0 coeffs
      cr1(0:ms1,0) = .5_dp*cr1(0:ms1,0)
      cz1(0:ms1,0) = .5_dp*cz1(0:ms1,0)
      cr(0:ms,0) = .5_dp*cr(0:ms,0)
      cz(0:ms,0) = .5_dp*cz(0:ms,0)

c................................................................
c     Write global info based on settings into output
      write (inesc, 10)
     1  '----- end inputs, begin outputs. Nescoil Version 1.0 -----'
      if( is_zero(cup) .and. is_zero(cut) ) then
         write (inesc, 10) '----- Solving for Saddle coils -----'
      else
         write (inesc, 10) '----- Solving for Modular coils -----'
      endif

      if( ibex .eq. 0 ) then
         write (inesc, 10) '----- No background coils used -----'
      else
         write (inesc, 10) '----- Background coils used -----'
      endif

 10   format(a)
 20   format(a, i4)

      return

c................................................................
c     error branch, help user
 99   call nescoil_help
      stop 'error while reading input file in  read_nescoil_input'

      end subroutine  read_nescoil_input
