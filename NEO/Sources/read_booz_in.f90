
subroutine read_booz_in

      use read_boozer_mod, dp1 => dp
      use neo_input
      use neo_units
      use neo_control
      use neo_exchange
      implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real (rprec) ::  hs
      integer :: k, i, j, i_surf, mn0
      integer :: i_alloc
      character*(120) :: booz_input_file
!
!     Read data from the boozmn file and allocate storage:
!
      booz_input_file = extension

      call read_boozer_file (booz_input_file,k)
      if (k .ne. 0) stop 'Error reading boozmn file'
      
      ! If LASYM_B detected then stop code.
      lasym=.false.
      IF (lasym_b) THEN
         lasym=.true.
!         STOP 'ERROR: LASYM_B Detected in BOOZMN File.'
      END IF

      m0b = mboz_b-1
      n0b = nboz_b
      mnmax = mnboz_b

      ns = ns_b
      nfp = nfp_b
      flux = abs( phi_b(ns_b) )
      max_n_mode = max_n_mode * nfp

!     truncate modes
!     assume m=0, mboz
!            n=-nboz, nboz
!     discard m> mb_retained
!             iabs(n)> nb_retained
!     perhaps better to retain largest mnmax modes
!             for a reference surfaces

      if (mnboz_b .ne. (mboz_b-1)*(2*nboz_b+1)+nboz_b+1 )   &
           stop "boozer modes inconsistent"

      k = 0
      if (max_m_mode .le. 0) max_m_mode = m0b
      if (max_n_mode .le. 0) max_n_mode = n0b * nfp
      do j = 1, mnboz_b
        if( ixm_b(j) .gt. max_m_mode ) cycle
        if( iabs(ixn_b(j) ) .gt. max_n_mode ) cycle
        k = k + 1
      enddo

      mnmax = k
      m_max = min0(m0b, max_m_mode)
      m_max = m_max + 1
      n_max = min0(n0b, max_n_mode/nfp)
      n_max = 2*n_max+1
      if ((no_fluxs.gt.0) .and. (no_fluxs.le.ns_b)) ns = no_fluxs

! **********************************************************************
! Allocate storage arrays
! **********************************************************************

      ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
      IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

      ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
      IF(i_alloc /= 0)                                              &
           STOP 'Allocation for integer arrays with pointers failed!'

      ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
      IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

      ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),        &
           pprime(ns), sqrtg00(ns), stat = i_alloc)
      IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

      ALLOCATE(rmnc(ns,mnmax), zmns(ns,mnmax), lmns(ns,mnmax),      &
           bmnc(ns,mnmax), stat = i_alloc)
      IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
      
      IF (lasym) THEN
         ALLOCATE(rmns(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),&
                  bmns(ns,mnmax), stat=i_alloc)
        IF(i_alloc /= 0) STOP &
           'Allocation for fourier arrays (2) failed!'
      END IF
! **********************************************************************

      k = 0;   mn0 = 0
      do j = 1, mnboz_b
         if( ixm_b(j) .gt. max_m_mode ) cycle
         if( iabs(ixn_b(j) ) .gt. max_n_mode ) cycle
         if (ixm_b(j).eq.0 .and. ixn_b(j).eq.0) mn0 = j
         k = k + 1
         ixm(k) = ixm_b(j)
         ixn(k) = ixn_b(j)
      enddo

      if (mn0 .eq. 0) stop 'M=0, N=0 MUST BE INCLUDE IN BOOZER SPECTRUM!'

      hs = ONE/(ns_b -1)
      do i = 1, ns                  !!NEED TO INCLUDE IN LOOP CASE FOR NO_FLUXS < 0 (SPH)

         if (no_fluxs > 0) then
            i_surf = fluxs_arr(i)   !!NOT ALLOCATED YET IF NO_FLUXS<0
         else
            i_surf = i
         end if
         k = 0
         do j = 1, mnboz_b

           if( ixm_b(j) .gt. max_m_mode ) cycle
           if( abs(ixn_b(j)) .gt. max_n_mode ) cycle
           k = k + 1

           rmnc(i,k) = rmnc_b(j,i_surf)
           zmns(i,k) = zmns_b(j,i_surf)
!
!       note: zeta - zeta_b from booz_xform,
!             Note: in BOOZMN file, zeta_b - zeta is stored
!             Also, lmnc is normalized to one field period
!

           lmns(i,k) = -pmns_b(j,i_surf)*nfp_b/twopi
           bmnc(i,k) = bmnc_b(j,i_surf)

           if (ixm_b(j).eq.0 .and. ixn_b(j).eq.0)            &
               sqrtg00(i) = gmnc_b(j,i_surf)
               
           if (lasym) then
              rmns(i,k) = rmns_b(j,i_surf)
              zmnc(i,k) = zmnc_b(j,i_surf)
              lmnc(i,k) = -pmnc_b(j,i_surf)*nfp_b/twopi
              bmns(i,i) = bmns_b(j,i_surf)
              if (ixm_b(j).eq.0 .and. ixn_b(j).eq.0)            &
               sqrtg00(i) = gmnc_b(j,i_surf) + gmns_b(j,i_surf)
           end if

         enddo

         if (rmnc(i,mn0) .eq. ZERO) then
            write (w_us, *)' The surface i = ', i_surf, ' is absent in the BOOZMN file'
            stop
         end if

         es(i) = (i_surf-1.5_dp)*hs
         iota(i) = iota_b(i_surf)
         if (i_surf .lt. ns_b) then
            pprime(i) = (pres_b(i_surf+1) - pres_b(i_surf))/hs
         else
            pprime(i) = 0
         end if
!
!       these are the F and I functions, respectively (currents divided by 2*pi, sign)
!
         curr_pol(i) = bvco_b(i_surf)
         curr_tor(i) = buco_b(i_surf)

      enddo


      call read_boozer_deallocate

END SUBROUTINE READ_BOOZ_IN
