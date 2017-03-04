!-----------------------------------------------------------------------
!     Subroutine:    stelltran_calcboot
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          09/03/2015
!     Description:   Calculate bootstrap current
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_calcboot
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime
      USE stelltran_vars
      USE stelltran_equilutils
      ! BOOZ_XFORM LIBRARIES
      USE boozer_utils, ONLY: load_boozer
      USE booz_params, mboz_xboozer => mboz, nboz_xboozer => nboz,&
                       lscreen_xboozer => lscreen, nfp_xboozer => nfp, &
                       ns_xboozer => ns, lasym_xboozer => lasym_b
      USE read_wout_mod, ONLY: input_extension_wout => input_extension,&
                               ns_vmec => ns, aspect_vmec => aspect,&
                               rmax_vmec => rmax_surf, rmin_vmec=> rmin_surf,&
                               betaxis_vmec => betaxis
      USE safe_open_mod, ONLY: safe_open
      USE booz_persistent
      USE read_boozer_mod
      ! BOOTSJ LIBRARIES
      USE bootsj_input
      use parambs, lscreen_bootsj=>lscreen, phip_bootsj=>phip, &
                   gpsi_bootsj=>gpsi
      use vmec0
      use trig
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      ! BOOZER related
      INTEGER ::  ier, ik, ij, iunit, nu, nv, u, v, mn
      INTEGER,SAVE ::  irun_setup_booz = 0
      INTEGER, ALLOCATABLE :: im(:), in(:)
      REAL(rprec), ALLOCATABLE :: xu(:), xv(:)
      REAL(rprec), ALLOCATABLE :: r_temp(:,:,:), z_temp(:,:,:), p_temp(:,:,:), g_temp(:,:,:), b_temp(:,:,:)
      REAL(rprec), ALLOCATABLE :: ru(:,:,:), rv(:,:,:)
      REAL(rprec), ALLOCATABLE :: zu(:,:,:), zv(:,:,:)
      REAL(rprec), ALLOCATABLE :: pu(:,:,:), pv(:,:,:)
      REAL(rprec), ALLOCATABLE :: fmn_temp(:,:)
      ! Bootstrap related
      integer, parameter :: nfax = 13
      INTEGER ::  i, j, ntheta, nzeta, nlis
      INTEGER ::  dex_ion, dex_zeff
      REAL(rprec) :: t1, t2
      integer, dimension(nfax) :: ifaxu, ifaxv
      integer :: ntrigu, ntrigv, i1
      integer :: ntheta_min, nzeta_min, ir, m, n
      integer :: irho, irho1, ierr, ijbs, ians, ians_plot
      real(rprec), dimension(:), allocatable :: cputimes
      real(rprec) :: time1, timecpu, unit, file, status, err,&
         time2, r, x, al31t, gradbs1, gradbs2,&
         gradbs3, gradbs4,  al31s, aibstot, a
      real(rprec) :: a1, tempe0, tempi0, pres10, pres0
      real(rprec), dimension(:), allocatable :: work
      integer :: ihere = 0
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      real(rprec) , EXTERNAL :: al31
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
	
      ! First transform to boozer coords
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000')
            lscreen_xboozer = lveryfirst_pass
            ! We need to pass mboz and nboz to the boozer routines
            mboz_xboozer = mboz
            nboz_xboozer = nboz
            ! Now read the wout file
            IF (ALLOCATED(lsurf_boz)) DEALLOCATE(lsurf_boz)
            ALLOCATE(lsurf_boz(ns_vmec))
            lsurf_boz = .true.
            iunit = -COUNT(lsurf_boz) ! This tricks the code into not reading the wout or input file
            ier = 0
            CALL read_wout_booz(TRIM(proc_string),iunit,ier)
            !CLOSE(iunit)
            IF (ier .ne. 0) PRINT *,'read_wout_booz_error'
            ! Calculate Boozer Coordinates
            DO ik = 2, ns_vmec
               IF (lsurf_boz(ik)) CALL boozer_coords(ik,irun_setup_booz)
            END DO
            ! Output the data
            CALL write_boozmn(TRIM(proc_string))
            ! Load the boozer data from memory
            mnboz_b = mnboz
            mboz_b  = mboz_xboozer
            nboz_b  = nboz_xboozer
            nfp_b   = nfp_xboozer
            ns_b    = ns_xboozer
            aspect_b = aspect_vmec
            rmax_b   = rmax_vmec
            rmin_b   = rmin_vmec
            betaxis_b = betaxis_vmec
            IF (ALLOCATED(idx_b)) DEALLOCATE(idx_b); ALLOCATE(idx_b(ns_b))
            IF (ALLOCATED(iota_b)) DEALLOCATE(iota_b); ALLOCATE(iota_b(ns_b))
            IF (ALLOCATED(pres_b)) DEALLOCATE(pres_b); ALLOCATE(pres_b(ns_b))
            IF (ALLOCATED(phip_b)) DEALLOCATE(phip_b); ALLOCATE(phip_b(ns_b))
            IF (ALLOCATED(phi_b)) DEALLOCATE(phi_b); ALLOCATE(phi_b(ns_b))
            IF (ALLOCATED(beta_b)) DEALLOCATE(beta_b); ALLOCATE(beta_b(ns_b))
            IF (ALLOCATED(buco_b)) DEALLOCATE(buco_b); ALLOCATE(buco_b(ns_b))
            IF (ALLOCATED(bvco_b)) DEALLOCATE(bvco_b); ALLOCATE(bvco_b(ns_b))
            idx_b    = 1
            idx_b(1) = 0
            iota_b = hiota
            pres_b = pres
            phip_b = phip
            phi_b  = phi
            beta_b = beta_vol
            buco_b = buco
            bvco_b = bvco
            IF (ALLOCATED(ixm_b)) DEALLOCATE(ixm_b); ALLOCATE(ixm_b(mnboz_b))   
            IF (ALLOCATED(ixn_b)) DEALLOCATE(ixn_b); ALLOCATE(ixn_b(mnboz_b))
            ixm_b = NINT(xmb)
            ixn_b = NINT(xnb)
            IF (ALLOCATED(bmnc_b)) DEALLOCATE(bmnc_b); ALLOCATE(bmnc_b(mnboz_b,ns_b))
            IF (ALLOCATED(rmnc_b)) DEALLOCATE(rmnc_b); ALLOCATE(rmnc_b(mnboz_b,ns_b))  
            IF (ALLOCATED(zmns_b)) DEALLOCATE(zmns_b); ALLOCATE(zmns_b(mnboz_b,ns_b))  
            IF (ALLOCATED(pmns_b)) DEALLOCATE(pmns_b); ALLOCATE(pmns_b(mnboz_b,ns_b))  
            IF (ALLOCATED(gmnc_b)) DEALLOCATE(gmnc_b); ALLOCATE(gmnc_b(mnboz_b,ns_b))  
            bmnc_b = 0.0; rmnc_b = 0.0; zmns_b = 0.0; pmns_b = 0.0; gmnc_b = 0.0
            lasym_b = .FALSE.
            IF (lasym_xboozer) THEN
               lasym_b = .TRUE.
               IF (ALLOCATED(bmns_b)) DEALLOCATE(bmns_b); ALLOCATE(bmns_b(mnboz_b,ns_b))
               IF (ALLOCATED(rmns_b)) DEALLOCATE(rmns_b); ALLOCATE(rmns_b(mnboz_b,ns_b))  
               IF (ALLOCATED(zmnc_b)) DEALLOCATE(zmnc_b); ALLOCATE(zmnc_b(mnboz_b,ns_b))  
               IF (ALLOCATED(pmnc_b)) DEALLOCATE(pmnc_b); ALLOCATE(pmnc_b(mnboz_b,ns_b))  
               IF (ALLOCATED(gmns_b)) DEALLOCATE(gmns_b); ALLOCATE(gmns_b(mnboz_b,ns_b))  
               bmns_b = 0.0; rmns_b = 0.0; zmnc_b = 0.0; pmnc_b = 0.0; gmns_b = 0.0
            END IF
            ! The internal BOOZER variables are in packed form (meaning that they run
            ! over the number of surfaces calculated not the total number of
            ! equilibrium surfaces.  Thus we need to get them in the read_boozer_mod
            ! form where the radial index runs over the total number of ns_b
            ! surfaces.
            ij = 1
            DO ik = 1, ns_b
              IF (idx_b(ik) == 0) CYCLE
              bmnc_b(:,ik) = bmncb(:,ij)
              rmnc_b(:,ik) = rmncb(:,ij)
              zmns_b(:,ik) = zmnsb(:,ij)
              pmns_b(:,ik) = pmnsb(:,ij)
              gmnc_b(:,ik) = gmncb(:,ij)
              IF (lasym_xboozer) THEN
                 bmns_b(:,ik) = bmnsb(:,ij)
                 rmns_b(:,ik) = rmnsb(:,ij)
                 zmnc_b(:,ik) = zmncb(:,ij)
                 pmnc_b(:,ik) = pmncb(:,ij)
                 gmns_b(:,ik) = gmnsb(:,ij)
              END IF
              ij = ij + 1
            END DO
            ! Free the memory
            CALL free_mem_boozer
            irun_setup_booz = 0
      END SELECT
      
      ! Now calculate the bootstrap current
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000')
            ! Get the data into bootsj
            call second0 (time1)
            timecpu = time1
            lscreen_bootsj = lveryfirst_pass
            CALL read_boozer('')
            ians = 12
            CALL safe_open(ians, ierr, 'answers.'//trim(proc_string), 'replace','formatted')
            ians_plot = 14
            CALL safe_open(ians_plot, ierr, 'answers_plot.'//trim(proc_string),'replace', 'formatted')
            ijbs = 15
            CALL safe_open(ijbs, ierr, 'jBbs.'//trim(proc_string),'replace', 'formatted')
            ! Assume the bootsj namelist has been read in
            if(damp_bs .lt. 0.0) then !in this case no damp_bs was read in
              if(damp .gt. 0.0) then
                 damp_bs = sqrt(damp)
              else
                 damp_bs = 0.001_dp      !in this case no damp was read in
              endif
            endif
            ! Dimensionality Check
            if(nboz_b .lt. nbuse) nbuse = nboz_b
            if(mboz_b .lt. mbuse) mbuse = mboz_b
            nzeta_min = 2*nbuse + 1
            ntheta_min = 2*mbuse + 1
            do i = 0, 6
               nzetah = 4*2**i
               if(nzetah .gt. nzeta_min) exit
               nzetah = 2*2**i * 3
               if(nzetah .gt. nzeta_min) exit
            enddo
            do i = 0, 6
               nthetah = 4*2**i
               if(nthetah .gt. ntheta_min) exit
               nthetah = 2*2**i * 3
               if(nthetah .gt. ntheta_min) exit
            enddo
            ! Convert bmn's to amnfit's
            amnfit = 0.0
            lsurf = .false.
            status = tiny(a1)
            do ir = 1, irup
               do mn = 1,mnboz_b
                  m = ixm_b(mn)
                  n = ixn_b(mn)/nfp_b
                  if (m .gt. mboz_b) stop 'boozmn indexing conflict, m'
                  if (abs(n) .gt. nboz_b) stop 'boozmn indexing conflict, n'
                  if (n.lt.0) then
                     m = -m
                     n = -n
                  end if
                  if (m.eq.0 .and. n.eq.0 .and. bmnc_b(mn,ir).gt.status) lsurf(ir) = .true.
                  amnfit(ir,m,n) = bmnc_b(mn,ir+1)     !!2nd half grid == 1st grid pt. here
                  IF (lasym_b) amnfit2(ir,m,n) = bmns_b(mn,ir+1) 
               end do
            end do
            ! Note we don't deallocate the boozer coordinates
            zeff1 = max(1.0_rprec,zeff1)
            ! Setup Mesh
            psimax = maxval(abs(flux))
            if(flux(irup) .lt. 0.0) psimax = -psimax
            do ir = 1, irup
              rhoar(ir) = 0.5_dp*(flux(ir) + flux(ir+1))/psimax
              d_rho(ir) = (flux(ir+1) - flux(ir))/psimax
            end do
            if (iotasign .lt. 0) then
               qsafety(:irup) = iotasign*qsafety(:irup)
            endif
            aiogar(:irup) = aipsi(:irup)/(gpsi_bootsj(:irup)+1.0e-36_dp)
            call positiv (pres1, irup, 2) !to be sure that arrays are positive
            call positiv (betar, irup, 2)
            ! Now we get TI and TE on VMEC mesh (our way)
            ! Notes: This part of the code wants quantities in 10^20 [m^-3] and
            !        [keV]
            DO ir = 1, irup	   
               CALL eval_prof_spline(prof_length,te(:,1),te(:,2),rhoar(ir),tempe1(ir),ier)
            END DO
            DO ir = 1, irup
               CALL eval_prof_spline(prof_length,ne(:,1),ne(:,2),rhoar(ir),dense(ir),ier)
            END DO
            IF (teti<=0) THEN
               DO ir = 1, irup
                  CALL eval_prof_spline(prof_length,ti(:,1),ti(:,2),rhoar(ir),tempi1(ir),ier)
               END DO
            ELSE
               tempi1 = tempe1 / teti
            END IF
            dense = dense + 1E-36_dp
            tempe1 = tempe1/1000.         ! [eV] to [keV]
            tempi1 = tempi1/1000.         ! [eV] to [keV]
            dense  = dense/(1.0E+20)       ! [m^-3] to 10^20 [m^-3]
            tempe0 = tempe1(1)            !central electron temperature in keV
            tempi0 = tempi1(1)                 !central ion temperature in keV
            pres10 = pres1(1)                    !central total pressure in Pa
            call positiv (tempe1, irup, 2)
            call positiv (tempi1, irup, 2)
            IF (ALLOCATED(work)) DEALLOCATE(work)
            allocate(work(irdim))
            call smooth1 (dense, 1, irup, work, 0.0)
            call positiv (dense, irup, 2)
            i1 = irup - 1
            a = tempe1(irup) + tempi1(irup)/zeff1
            a1 = tempe1(i1) + tempi1(i1)/zeff1
            dense(irup) = dense(i1)*a1*betar(irup)/(a*betar(i1)+1.E-36_dp)
            DO ir = 1, irup
               CALL eval_prof_spline(prof_length,zeff(:,1),zeff(:,2),rhoar(ir),zeff1,ier)
               densi(ir) = dense(ir)/zeff1
            END DO
            CALL eval_prof_spline(prof_length,zeff(:,1),zeff(:,2),rhoar(1),zeff1,ier)
            !dex_zeff = MINLOC(zeff_aux_s(2:),DIM=1)
            !IF (dex_zeff > 4) THEN
            !   DO ir = 1, irup
            !      CALL get_equil_zeff(rhoar(ir),TRIM(zeff_type),zeff1,ier)
            !      densi(ir) = dense(ir)/zeff1
            !   END DO
            !   CALL get_equil_zeff(rhoar(1),TRIM(zeff_type),zeff1,ier)
            !ELSE
            !   densi(:irup) = dense(:irup)/zeff1
            !END IF
            dens0 = dense(1)                         !central electron density
            DEALLOCATE(work)
            ntrigu = 3*nthetah/2 + 1
            ntrigv = 2*nzetah
            IF (ALLOCATED(cputimes)) DEALLOCATE(cputimes)
            IF (ALLOCATED(dmn)) DEALLOCATE(dmn)
            IF (ALLOCATED(fmn)) DEALLOCATE(fmn)
            IF (ALLOCATED(rfmn)) DEALLOCATE(rfmn)
            IF (ALLOCATED(alpha1mn)) DEALLOCATE(alpha1mn)
            IF (ALLOCATED(trigsv)) DEALLOCATE(trigsv)
            IF (ALLOCATED(trigsu)) DEALLOCATE(trigsu)
            allocate (cputimes(irup))
            allocate (dmn(-mbuse:mbuse,0:nbuse), fmn(-mbuse:mbuse,0:nbuse),&
                rfmn(-mbuse:mbuse,0:nbuse),alpha1mn(-mbuse:mbuse,0:nbuse),&
                trigsv(ntrigv),trigsu(ntrigu), stat=irho)
            if (irho .ne. 0) stop 'allocation error in bootsj main'
            ! Now we need to calculate indexes
            ! Now loop over surfaces
            idx(1:irup) = 1
            l_boot_all = .TRUE.
            ifaxu = 1; ifaxv = 1
            call fftfax_g (nthetah, ifaxu, trigsu) ! Initialize trigsu
            call cftfax_g (nzetah, ifaxv, trigsv)  ! Initialize trigsv
            IF (lveryfirst_pass) THEN
               WRITE(6,'(A)')             '==========================================='
               WRITE(6,'(A,A,A)')         '=========  B O O T S J (v.',version_(1:4),')  =========='
               WRITE(6,'(2X,2(A,I3.3))')   'M_CALC = ',mbuse,';   M_BOOZER = ',mboz 
               WRITE(6,'(2X,2(A,I3.3))')   'N_CALC = ',nbuse,';   N_BOOZER = ',nboz
               WRITE(6,'(2X,2(A,I3.3))')   'NTHETA = ',nthetah,';   NZETA    = ',nzetah
               WRITE(6,'(2X,A,F10.4)')   'ZEFF =',zeff1
               WRITE(6,'(2X,A,F10.4,A)') 'TE0  =',tempe0,' [keV]'
               WRITE(6,'(2X,A,F10.4,A)') 'TI0  =',tempi0,' [keV]'
               WRITE(6,'(2X,A,F10.4,A)') 'NE0  =',dense(1),'x10^20 [m^-3]'
               WRITE(6,'(2X,A,F10.4,A)') 'NI0  =',densi(1),'x10^20 [m^-3]'
               IF (l_boot_all) WRITE(6,'(2X,A)') '<FULL CURRENT CALCULATTION>'
               WRITE(6,'(A)') '-------------------------------------------'
               WRITE(6,'(A)') '   dex      rho      Te[keV]     Ti[keV]       Ne           Ni        BETA       J_BOOT     TOK_FRAC'
               CALL FLUSH(6)             
            END IF
            DO irho = 1, irup
               IF (idx(irho) .eq. 0) CYCLE
               call second0 (time2)
               timecpu = time2 - time1
               cputimes(irho) = timecpu
               irho1 = irho -1
               r     = SQRT(rhoar(irho) + 1.0E-36_dp)
               CALL bongrid(irho,ians,ihere)
               x     = fttok(irho)/(fptok(irho)+1.0E-36_dp)
               al31t = al31(x,zeff1,alphae,alphai)
               CALL grad(gradbs1,gradbs2,gradbs3,gradbs4,irho)
               bsdenste(irho) = gradbs1*al31t          !due to dens gradient
               bsdensti(irho) = gradbs2*al31t          !due to dens gradient
               bstempte(irho) = gradbs3*al31t          !due to temp gradient
               bstempti(irho) = gradbs4*al31t          !due to temp gradient
               dibst(irho) = bsdenste(irho) + bsdensti(irho) + bstempte(irho) + bstempti(irho) !total Jbst
               IF (l_boot_all) THEN
                  IF (irho .eq. 1) THEN
                     aibst(1) = dibst(1)*d_rho(1)
                  ELSE
                     aibst(irho) = aibst(irho1)+dibst(irho)*d_rho(irho)
                  END IF
               END IF
               CALL denmf (trigsu, trigsv, ifaxu, ifaxv, irho)
               CALL caprsh2(irho)
               IF (h2(irho) == 0.0) CYCLE
               CALL woflam (trigsu, trigsv, ifaxu, ifaxv, irho)
               CALL othersums(irho)
               amain(irho) = 0.5_rprec*(1.0_rprec - aiogar(irho)/qsafety(irho)) + 0.5_rprec*(1.0_rprec + aiogar(irho)/qsafety(irho))*h2(irho)
               gbsnorm(irho) = amain(irho) + other1(irho) +  aiterm1(irho)
               x = ftrapped(irho)/(fpassing(irho)+1.E-36_dp)
               al31s = al31(x,zeff1,alphae,alphai)
               call grad (gradbs1, gradbs2, gradbs3, gradbs4, irho)
               bsdense(irho) = gbsnorm(irho)*gradbs1*al31s
               bsdensi(irho) = gbsnorm(irho)*gradbs2*al31s
               bstempe(irho) = gbsnorm(irho)*gradbs3*al31s
               bstempi(irho) = gbsnorm(irho)*gradbs4*al31s
               dibs(irho) = bsdense(irho) + bsdensi(irho) + bstempe(irho) + bstempi(irho)
               ajBbs(irho) = (2.0e6_dp)*dmu0*dibs(irho)*(pres1(irho)/betar(irho))/psimax 
               IF (l_boot_all) THEN
                  IF (irho .eq. 1) THEN
                     aibs(1) = dibs(1)*d_rho(1)
                  ELSE
                     aibs(irho) = aibs(irho1)+dibs(irho)*d_rho(irho)
                  END IF
               END IF
               bsnorm(irho) = dibs(irho)/(dibst(irho)+1.E-36_dp)
               call second0 (time2)
               timecpu = time2 - time1
               cputimes(irho) = timecpu
               IF (lscreen_bootsj) WRITE(6,'(2X,I3,8(2X,E11.4))') irho,rhoar(irho),tempe1(irho),tempi1(irho),dense(irho),densi(irho),betar(irho),ajBbs(irho),bsnorm(irho)
               CALL FLUSH(6) 
            END DO
            IF (lscreen_bootsj .and. l_boot_all) WRITE(6,'(1X,A,E22.14,A)') 'Total bootstrap current: ',aibs(irup),' [MA]'
            l_boot_all = .FALSE. ! To Prevent output from writing it to the screen again.
            IF (lveryfirst_pass) THEN
               WRITE(6,'(A)')         '==========================================='
            END IF
            dibs = dibs*1E6 ! Get jboot in A.
            CALL output(cputimes,aibstot,ijbs,ians,ians_plot)
            IF (ALLOCATED(cputimes)) DEALLOCATE(cputimes)
            CLOSE(ians)
            CLOSE(ians_plot)
            CLOSE(ijbs)
      END SELECT
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE stelltran_calcboot
