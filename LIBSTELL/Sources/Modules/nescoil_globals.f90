!-----------------------------------------------------------------------
!     Module:        nescoil_globals
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          03/20/2026
!     Description:   This module contains the NESCOIL global variables
!                    needed by the input namelist.
!-----------------------------------------------------------------------
      MODULE nescoil_globals
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE

      ! Moved from num
      INTEGER :: inesc = 6

      ! Moved from Outctrl
      INTEGER :: w_psurf, w_csurf, w_bnuv, w_jsurf, w_xerr, w_svd

      ! Moved from SvdCtrl
      integer :: mstrt, mstep, mkeep, mdspw, noaccuracy
      real(rprec) ::  curwt, trgwt

      ! Moved from Vmeshes
      integer :: nu, nv, nuv, nuvh, nu1, nv1, nuv1, nuvh1, &
                 mf, nf, npol, ntor, md, nd, nmax, mnd

      ! Moved from Vprecal1
      integer :: np

      ! Moved from Vvacuum1
      integer   ms, ns
      real(rprec), dimension(:,:), allocatable :: cr, cz

      ! Moved from Vvacuum2
      integer   ms1, ns1
      real(rprec), dimension(:,:), allocatable :: cr1, cz1, cl1
      real(rprec) :: iota_edge, phip_edge, curpol

      ! Moved from Vvacuum3
      integer   ibex
      real(rprec) ::  cup, cut

      ! Moved from Vvacuum4
      real(rprec), dimension(:,:), allocatable :: cr2, cz2

      ! Moved from Vvacuum5
      real(rprec), dimension(:,:), allocatable :: cr3, cz3

      ! Moved from Vvacuum6
      real(rprec), dimension(:,:), allocatable :: cf, sf


      END MODULE nescoil_globals
