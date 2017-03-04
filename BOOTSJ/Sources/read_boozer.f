
      subroutine read_boozer(extension)
      use read_boozer_mod
      use parambs
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*) :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat
C-----------------------------------------------
!
!     NOTE: read_boozer_deallocate is called later, after bmn_b are assigned
!
      IF (.not. ALLOCATED(idx_b)) THEN
         call read_boozer_file(extension, istat)
         if (istat .ne. 0) then
            print *,'Error reading boozer file in BOOTSJ, istat = ',
     1               istat
            stop
         end if
      END IF

      irdim   = ns_b
      irup    = ns_b - 1                !No. points in radial (half-mesh) profiles
      periods = nfp_b                   !No. field periods
      lasym_bootsj = lasym_b

c  now that we know the number of radial points,  radial quantites
c  can be allocated

      call allocate_radial

c  LAB--change the indexing on aipsi and gpsi to relect half mesh status
      aipsi(1:irup)    = buco_b(2:ns_b)                 !Boozer I
      gpsi (1:irup)    = bvco_b(2:ns_b)                 !Boozer g
      qsafety(1:irup)  = one/(iota_b(2:ns_b) +
     1                        sign(1.0e-14_dp,iota_b(2:ns_b)))
      pres1(1:irup)    = pres_b(2:ns_b)
      betar(1:irup)    = beta_b(2:ns_b)
      idx(1:irup)      = idx_b(2:ns_b)
      flux(2:ns_b)     = phi_b(2:ns_b)
      phip(1:irup)     = phip_b(2:ns_b)


      sign_jacobian = one   !version 6.1 phi_b has the sign of the physical flux
c     phip_b retains the internal vmec convention.
      if( gpsi(irup)*phip_b(ns_b) <= zero) sign_jacobian = -one
      flux(1) = zero

      end subroutine read_boozer
