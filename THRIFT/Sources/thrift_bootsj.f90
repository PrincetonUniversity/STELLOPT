!-----------------------------------------------------------------------
!     Subroutine:    thrift_bootsj
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/06/2012
!     Description:   This subroutine calculates the boostrap current
!-----------------------------------------------------------------------
      SUBROUTINE thrift_bootsj(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_input_mod
      USE thrift_vars, nrho_thrift => nrho
      USE thrift_profiles_mod
      USE thrift_equil
      USE thrift_funcs
      USE mpi_params
      USE mpi_inc
      USE EZspline
      USE EZspline_obj
      ! BOOTSJ LIBRARIES
      USE bootsj_input
      use parambs, lscreen_bootsj=>lscreen
      use vmec0
      use read_boozer_mod
      use trig
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      integer, parameter :: nfax = 13
      INTEGER ::  ier, i, j, ik, ntheta, nzeta, nlis
      INTEGER ::  dex_ion, dex_zeff
      INTEGER :: mystart,myend, chunk, numprocs_local
      INTEGER, ALLOCATABLE :: mnum(:)
      integer, dimension(nfax) :: ifaxu, ifaxv
      integer :: ntrigu, ntrigv, i1
      integer :: ntheta_min, nzeta_min, ir, mn, m, n
      integer :: irho, irho1, ierr, iunit, ijbs, ians, ians_plot
      real(rprec), dimension(:), allocatable :: cputimes
      real(rprec) :: time1, timecpu, unit, file, status, err,&
         time2, r, x, al31t, gradbs1, gradbs2, temp, &
         gradbs3, gradbs4,  al31s, aibstot, a, roa, zeff_temp,&
         s_val, rho_val
      real(rprec) :: a1, tempe0, tempi0, pres10, pres0, jdotb
      real(rprec), dimension(:), allocatable :: work
      integer :: ihere = 0
      CHARACTER(LEN=32) :: temp_str
      ! Helpers to get dI/ds
      REAL(rprec) :: vp, dPhidrho
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rho_temp, dIds_temp,j_temp
      TYPE(EZspline1_r8) :: dIds_spl
      INTEGER :: bcs0(2)
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      real(rprec) , EXTERNAL :: al31
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' --------------------  BOOTSTRAP CALCULATION USING BOOTSJ  -------------------'
      IF (lvmec) THEN
            ! Get the data into bootsj
            call second0 (time1)
            timecpu = time1
            lscreen_bootsj = lscreen
            ! Open Files for Output
            IF (myworkid == master) THEN
               CALL read_boozer('')
               ians = 12
               CALL safe_open(ians, ierr, 'answers.'//trim(proc_string), 'replace','formatted')
               ians_plot = 14
               CALL safe_open(ians_plot, ierr, 'answers_plot.'//trim(proc_string),'replace', 'formatted')
               ijbs = 15
               CALL safe_open(ijbs, ierr, 'jBbs.'//trim(proc_string),'replace', 'formatted')
            ELSE
               ians = 12
               WRITE(temp_str,'(I8.8)') myworkid
               CALL safe_open(ians, ierr, 'answers_'//TRIM(temp_str)//'.'//trim(proc_string),'replace', 'formatted')
            END IF
#if defined(MPI_OPT)
            CALL BCAST_BOOTSJ_INPUT(master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_COMM_SIZE( MPI_COMM_MYWORLD, numprocs_local, ierr_mpi )
            CALL MPI_BCAST(mnboz_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(mboz_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(nboz_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(nfp_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(ns_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(irdim,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(irup,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(periods,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(lasym_b,1,MPI_LOGICAL,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(lasym_bootsj,1,MPI_LOGICAL,master,MPI_COMM_MYWORLD,ierr_mpi)
            IF (myworkid /= master) CALL allocate_radial
            CALL MPI_BCAST(flux,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(qsafety,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(aipsi,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(idx,irup,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(gpsi,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(pres1,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(betar,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(phip,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(sign_jacobian,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            IF (.not. ALLOCATED(ixm_b)) ALLOCATE(ixm_b(mnboz_b))
            IF (.not. ALLOCATED(ixn_b)) ALLOCATE(ixn_b(mnboz_b))
            IF (.not. ALLOCATED(bmnc_b)) ALLOCATE(bmnc_b(mnboz_b,ns_b))
            IF (.not. ALLOCATED(bmns_b) .and. lasym_b) ALLOCATE(bmns_b(mnboz_b,ns_b))
            CALL MPI_BCAST(ixm_b,mnboz_b,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(ixn_b,mnboz_b,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            ik = SIZE(bmnc_b)
            CALL MPI_BCAST(bmnc_b,ik,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            IF (lasym_b) CALL MPI_BCAST(bmns_b,ik,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            ! Handle unbroadcast variables
            dibs=0; aibs = 0; dibst = 0; aibst = 0; bsdense = 0; bsdensi = 0; bstempe = 0
            bstempi = 0; bsdenste = 0; bsdensti = 0; bstempte = 0; bstempti = 0
            capr = 0; caps = 0; h2 = 0; ftrapped = 0; fpassing =0; epsttok = 0
            fttok = 0; gbsnorm = 0; aiterm1 = 0; other1 = 0; rhoar = 0; bsnorm = 0
            fptok = 0; amain = 0; bmax1 = 0; thetamax = 0; zetahmax = 0; ajBbs = 0
            d_rho = 0; b2avg = 0 

            ! Divide up work
            IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
            ALLOCATE(mnum(numprocs_local))
            mnum=0
            i = 1
            DO
               IF (SUM(mnum,DIM=1) == irup) EXIT
               IF (i > numprocs_local) i = 1
               mnum(i) = mnum(i) + 1
               i=i+1
            END DO
            mystart = 1
            DO i = 1, myworkid
               mystart = SUM(mnum(1:i))+1
            END DO
            myend = mystart + mnum(myworkid+1) - 1
            IF (myend < mystart) myend = mystart
            IF (mnum(myworkid+1) == 0) mystart = myend + 1
            DEALLOCATE(mnum)

#endif
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
            aiogar(:irup) = aipsi(:irup)/(gpsi(:irup)+1.0e-36_dp)
            call positiv (pres1, irup, 2) !to be sure that arrays are positive
            call positiv (betar, irup, 2)
            ! Now we get TI and TE on VMEC mesh (our way)
            ! Notes: This part of the code wants quantities in 10^20 [m^-3] and
            !        [keV]
            ier = 0
            IF (myworkid == master) THEN
               DO ir = 1, irup
                  roa = sqrt(rhoar(ir))
                  CALL get_prof_te(roa, THRIFT_T(mytimestep), tempe1(ir))
                  CALL get_prof_ne(roa, THRIFT_T(mytimestep), dense(ir))
                  CALL get_prof_ti(roa, THRIFT_T(mytimestep), 1, tempi1(ir))
                  dense(ir) = dense(ir) + 1.E-36_dp
               END DO
            END IF
#if defined(MPI_OPT)
            CALL MPI_BCAST(tempe1,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(dense,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(tempi1,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
#endif
            !WHERE (tempe1 <10)   tempe1 = 10
            !WHERE (tempi1 <10)   tempi1 = 10
            !WHERE (dense  <1E17) tempi1 = 1E17
            tempe1 = tempe1/1000.         ! [eV] to [keV]
            tempi1 = tempi1/1000.         ! [eV] to [keV]
            dense  = dense/(1.0E+20)       ! [m^-3] to 10^20 [m^-3]
            tempe0 = tempe1(1)            !central electron temperature in keV
            tempi0 = tempi1(1)                 !central ion temperature in keV
            pres10 = pres1(1)                    !central total pressure in Pa
            pres0 = 1.6022E4_DP
            call positiv (tempe1, irup, 2)
            call positiv (tempi1, irup, 2)
            IF (ALLOCATED(work)) DEALLOCATE(work)
            allocate(work(irdim))
            call smooth1 (dense, 1, irup, work, 0.0)
            call positiv (dense, irup, 2)
            DEALLOCATE(work)
            ! This was to control the density at the edge but can have
            ! some seriously disasterous concequeneces.
            !i1 = irup - 1
            !a = tempe1(irup) + tempi1(irup)/zeff1
            !a1 = tempe1(i1) + tempi1(i1)/zeff1
            !dense(irup) = dense(i1)*a1*betar(irup)/(a*betar(i1)+1.E-36_dp)
            !IF (myworkid == master) WRITE(6,*) a, a1, betar(i1)
            IF (myworkid == master) THEN
               roa = 0
               CALL get_prof_zeff(roa, THRIFT_T(mytimestep), zeff1)
               DO ir = 1, irup
                  roa = sqrt(rhoar(ir))
                  CALL get_prof_ni(roa, THRIFT_T(mytimestep), 1, densi(ir))
               END DO
               densi  = densi/(1.0E+20)       ! [m^-3] to 10^20 [m^-3]
            END IF
#if defined(MPI_OPT)
            CALL MPI_BCAST(zeff1,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(densi,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
#endif
            dens0 = dense(1)                         !central electron density
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
            cputimes=0; dmn=0; fmn=0; rfmn=0; alpha1mn=0; trigsv=0; trigsu=0
            ! Now we need to calculate indexes
            ! Now loop over surfaces
            ihere = 0
            idx = 0
            DO irho=1,irup
               idx(irho) = 1
            END DO
            l_boot_all = .TRUE.
            do irho=1, irup
               if(idx(irho) .eq. 0) l_boot_all = .false.
            enddo
            call fftfax_g (nthetah, ifaxu, trigsu) ! Initialize trigsu
            call cftfax_g (nzetah, ifaxv, trigsv)  ! Initialize trigsv
            IF (lscreen) THEN
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
            bsnorm =0; capr = 0; caps = 0; ftrapped =0; h2 =0;
            amain =0; aiterm1 = 0; other1 =0;
            cputimes = 0; dibs=0; aibs=0; gbsnorm =0;
            bsdenste = 0; bsdensti = 0; bstempte = 0; bstempti = 0;
            bsdense = 0; bsdensi = 0; bstempe = 0; bstempi = 0;
            dibst = 0; aibst = 0; amain = 0; ajBbs = 0;
            DO irho = mystart, myend
               IF (idx(irho) .eq. 0) CYCLE
#if defined(MPI_OPT)
               time1 = MPI_WTIME()
#else
               call second0 (time1)
#endif
               !timecpu = time2 - time1
               !cputimes(irho) = timecpu
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
#if defined(MPI_OPT)
               time2 = MPI_WTIME()
#else
               call second0 (time2)
#endif
               timecpu = time2 - time1
               cputimes(irho) = time2 - time1
               IF (lscreen_bootsj) WRITE(6,'(2X,I3,8(2X,E11.4))') irho,rhoar(irho),tempe1(irho),tempi1(irho),dense(irho),densi(irho),betar(irho),ajBbs(irho),bsnorm(irho)
               CALL FLUSH(6) 
            END DO
#if defined(MPI_OPT)
            CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
            IF (myworkid == master) THEN
               CALL MPI_REDUCE(MPI_IN_PLACE,dibs,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,ajBbs,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,bsnorm,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,capr,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,caps,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,ftrapped,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,h2,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,amain,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,aiterm1,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,other1,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,gbsnorm,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,thetamax,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,zetahmax,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,bmax1,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,fttok,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,fptok,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,fpassing,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,cputimes,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,aibs,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,bsdense,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,bsdensi,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,bstempe,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,bstempi,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,aibst,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            ELSE
               CALL MPI_REDUCE(dibs,dibs,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(ajBbs,ajBbs,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(bsnorm,bsnorm,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(capr,capr,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(caps,caps,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(ftrapped,ftrapped,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(h2,h2,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(amain,amain,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(aiterm1,aiterm1,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(other1,other1,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(gbsnorm,gbsnorm,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(thetamax,thetamax,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(zetahmax,zetahmax,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(bmax1,bmax1,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(fttok,fttok,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(fptok,fptok,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(fpassing,fpassing,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(cputimes,cputimes,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(aibs,aibs,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(bsdense,bsdense,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(bsdensi,bsdensi,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(bstempe,bstempe,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(bstempi,bstempi,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(aibst,aibst,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               IF (ALLOCATED(cputimes)) DEALLOCATE(cputimes)
               IF (ALLOCATED(bmnc_b)) DEALLOCATE(bmnc_b)
               IF (ALLOCATED(bmns_b)) DEALLOCATE(bmns_b)
               CALL FLUSH(ians)
               CLOSE(UNIT=ians,STATUS='DELETE')
               CALL deallocate_all
               RETURN
            END IF
#endif
            !  Now compute values for BOOTSJ
            ! rhoar : norm toroidal flux (s)
            ! diBs : dI/ds - what VMEC needs
            ! But we need j = dI/ds*Aminor/dVds
            ! But we need j = dI/ds * ds/dA
            !               = dI/ds / (dA/ds)
            !               = dI/ds / (dV/ds / (2*pi*Rmajor)) 
            !               = dI/ds * 2 * pi * Rmajor / dV/ds
            IF (myworkid == master) THEN
               ALLOCATE(rho_temp(irup+2),dIds_temp(irup+2))
               rho_temp(1)        = 0.0
               rho_temp(2:irup+1) = rhoar
               rho_temp(irup+2)   = 1.0
               dIds_temp(2:irup+1)   = dibs*1E6 ! dibs is in MA
               dIds_temp(1)          = 2*dIds_temp(2)-dIds_temp(3)
               dIds_temp(irup+2)     = 2*dIds_temp(irup+1)-dIds_temp(irup)
               bcs0=(/ 0, 0/)
               CALL EZspline_init(dIds_spl,irup+2,bcs0,ier)
               dIds_spl%x1        = sqrt(rho_temp)
               dIds_spl%isHermite = 1
               CALL EZspline_setup(dIds_spl,dIds_temp,ier,EXACT_DIM=.true.)
               DEALLOCATE(rho_temp,dIds_temp)

               ! Calculate J in s space = dI/ds * 1/(pi*a^2)
               ! dIds is in s space but grids dont necessarily agree
               ALLOCATE(j_temp(nsj))
               DO i = 1, nsj
                  s_val = THRIFT_S(i)
                  CALL EZspline_interp(dIds_spl,s_val,temp,ier)
                  j_temp(i) = temp/(pi*eq_Aminor**2) ! for some reason 'pi' is an ambigious reference
               END DO
               CALL EZspline_free(dIds_spl,ier)

               ! Convert to J in rho space
               CALL Js_to_Jrho(j_temp, THRIFT_JBOOT(:,mytimestep))
               DEALLOCATE(j_temp)

            END IF

            !  Output to screen the total bootstrap current
            IF ((myend .lt. irup).and.(lscreen_bootsj)) THEN
               DO irho=myend+1,irup
                  IF (idx(irho) .eq. 0) CYCLE
                  WRITE(6,'(2X,I3,8(2X,E11.4))') irho,rhoar(irho),tempe1(irho),tempi1(irho),dense(irho),densi(irho),betar(irho),ajBbs(irho),bsnorm(irho)
               END DO
            END IF
            ! Recompute AIBS
            IF (l_boot_all) THEN
               aibs(1) = dibs(1)*d_rho(1)
               DO irho = 2,irup
                  aibs(irho) = aibs(irho-1) + dibs(irho)*d_rho(irho)
               END DO
            END IF
            IF (lscreen_bootsj .and. l_boot_all) WRITE(6,'(1X,A,E22.14,A)') 'Total bootstrap current: ',aibs(irup),' [MA]'
            l_boot_all = .FALSE. ! To Prevent output from writing it to the screen again.
            IF (lscreen) THEN
               WRITE(6,'(A)')         '==========================================='
            END IF
            CALL output(cputimes,aibstot,ijbs,ians,ians_plot)
            IF (ALLOCATED(cputimes)) DEALLOCATE(cputimes)
            CLOSE(ians)
            CLOSE(ians_plot)
            CLOSE(ijbs)
      END IF
      IF (lscreen) WRITE(6,'(a)') ' -------------------  BOOTSJ BOOTSTRAP CALCULATION DONE  ---------------------'
      RETURN
  90  format(5e16.8)
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_bootsj
