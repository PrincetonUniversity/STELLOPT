
! ----------------------------------------------------------------------
!     ROUTINES FOR WRITING NESCOIL INPUT FILE (Added 9/98: SPH)
! ----------------------------------------------------------------------
      subroutine bn_write_nescoil_input(extension)
      use meshes
      use bnvariables
      use neswrite, only: coil_separation, iota_edge, phip_edge, mnmax,
     1   ixm, ixn, nfp, ntheta, nzeta, rmnc, zmns, rmns, zmnc,
     2  ixm_ws, ixn_ws, mnmax_ws, rmnc_ws, zmns_ws, lasym_bn,
     3  rmns_ws, zmnc_ws
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*) :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: ibex = 0
      integer, parameter :: nescoil0 = 15
      integer, parameter :: nmod_seg = 64, nfilaments = 10
      integer :: nu_plasma, nv_plasma
      real(rprec) :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat, mn, m, n, iunit
      real(rprec) :: cut, cup
C-----------------------------------------------
      iunit = nescoil0
      call safe_open(iunit, istat, 'nescin.' // extension, 'replace',
     1      'formatted')
      if (istat .ne. 0)
     1      stop 'Unable to open nescoil input file in bnorm'

!-----------------------------------------------
!     Write data for NESCOIL. First 3 lines (field periods) are not in
!     original NESCOIL code, but are needed here.
!-----------------------------------------------
      cut = zero         !Make these readable from input (wout file)
      cup = one          !1 for modular, 0 for saddles...

      write (iunit, 10) '------ Spatial dimensions ----'
      write (iunit, 10) 'nu, nv, nu1, nv1, npol, ntor, lasym_bn'
      nu_plasma = nu; nv_plasma = nv
      write (iunit,*) nu, nv, nu_plasma, nv_plasma, nmod_seg,
     1                nfilaments, lasym_bn
      write (iunit,*)

      write (iunit, 10) '------ Fourier Dimensions ----'
      write (iunit, 10) 'mf, nf, md, nd (max in surf and bnorm files)'
      write (iunit, *) mf, nf, md, nd
      write (iunit,*)

      write (iunit, 10) '------ Plasma information from VMEC ----'
      write (iunit, 10) 'np     iota_edge       phip_edge       curpol'
      write (iunit, *) nfp, iota_edge, phip_edge, curpol
      write (iunit,*)

      write (iunit, 10) '------ Current Controls ----'
      write (iunit, 10) 'cut  cup  ibex(=1,use fixed background coils)'
      write (iunit,*) cut, cup, ibex
      write (iunit,*)

      write (iunit, 10) '------ SVD controls -----'
      write (iunit, 10) 'mstrt, mstep, mkeep, mdspw, curwt, trgwt'
      write (iunit,*) 0,     0,      0,     4,    0.0_dp,   0.0_dp
      write (iunit,*)

      write (iunit, 10) '------ Output controls -----'
      write (iunit, 10) 'w_psurf w_csurf w_bnuv w_jsurf w_xerr w_svd'
      write (iunit,*) 0,       0,       0,      0,       0,      0
      write (iunit,*)

      write (iunit, 10) '------ Plasma Surface ---- '
      write (iunit, 10) 'Number of fourier modes in table'
      write (iunit,*) mnmax
      write (iunit, 10) 'Table of fourier coefficients'
      write (iunit, 10) 'm,n,crc,czs,cls,crs,czc,clc'
      do mn = 1, mnmax
         m = ixm(mn)
         n = ixn(mn)        !(nfp divided out already,
                            ! and n->-n for nescoil convention)
         write (iunit,'(x,2i6,1p6e20.12)') m, n,
     1         cr(m,n), cz(m,n), cl(m,n), crs(m,n), czc(m,n), clc(m,n)
      end do

      ntheta = 2*(mb+1) + 6
      ntheta = max(64, ntheta)
      nzeta  = 2*nb + 4
      nzeta  = max(32, nzeta)

      allocate (rmnc(mnmax), zmns(mnmax), rmns(mnmax), zmnc(mnmax),
     1          stat=istat)
      if (istat .ne. 0) stop 'Allocation error in bnorm bn_write'

      rmnc = zero; zmns = zero; rmns = zero; zmnc = zero;

      do mn = 1, mnmax
         m = ixm(mn)
         n = ixn(mn)              !!Nescoil convention: mu + nv
         rmnc(mn) = cr(m,n)
         zmns(mn) = cz(m,n)
      end do
      
      if (lasym_bn) then
         do mn = 1, mnmax
            rmns(mn) = crs(m,n)
            zmnc(mn) = czc(m,n)
         enddo
      endif

      call scaleup_boundary(ntheta, nzeta)
      
      if (.not.lasym_bn) then
         rmns_ws = zero
         zmnc_ws = zero
      endif

      write (iunit,*)
      write (iunit, '(a, 1pe20.12, a)')
     1  '------ Current Surface: Coil-Plasma separation = ',
     2  coil_separation,' -----'
      write (iunit, 10) 'Number of fourier modes in table'
      write (iunit,*) mnmax_ws
      write (iunit, 10) 'Table of fourier coefficients'
      write (iunit, 10) 'm,n,crc2,czs2,crs2,czc2'

      do mn = 1, mnmax_ws
        write (iunit,'(x,2i6,1p4e20.12)') ixm_ws(mn),
     1       ixn_ws(mn), rmnc_ws(mn), zmns_ws(mn),
     2                   rmns_ws(mn), zmnc_ws(mn)
      end do

      close (iunit)
 10   format (a)


      deallocate (rmnc_ws, zmns_ws, rmns_ws, zmnc_ws, ixm_ws, ixn_ws)
      end subroutine bn_write_nescoil_input
