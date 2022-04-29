!-----------------------------------------------------------------------
!     Module:        beams3d_input_mod
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de) M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          06/20/2012
!     Description:   This module contains the FIELDLINES input namelist and
!                    subroutine which initializes and reads the
!                    FIELDLINES input namelist.
!-----------------------------------------------------------------------
      MODULE beams3d_input_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime
      USE beams3d_lines, ONLY: nparticles, ns_prof1, ns_prof2, ns_prof3, &
                               ns_prof4, ns_prof5, partvmax
      USE beams3d_grid, ONLY: nr, nphi, nz, rmin, rmax, zmin, zmax, &
                              phimin, phimax, vc_adapt_tol, nte, nne, nti,&
                              nzeff, npot, plasma_mass, plasma_Zavg, &
                              plasma_Zmean, therm_factor, &
                              B_kick_min, B_kick_max, freq_kick, E_kick, &
                              rmin_fida, rmax_fida, zmin_fida, zmax_fida, phimin_fida, phimax_fida, &
                              raxis_fida, zaxis_fida, phiaxis_fida, nr_fida, nphi_fida, nz_fida, &
                              nenergy_fida, npitch_fida, energy_fida, pitch_fida
      USE safe_open_mod, ONLY: safe_open
      USE mpi_params
      USE mpi_inc

!-----------------------------------------------------------------------
!     Module Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      ! These are helpers to give the ns1_prof variables user friendly names
      INTEGER :: nrho_dist, ntheta_dist, nzeta_dist, nvpara_dist, nvperp_dist
      REAL(rprec) :: temp
!-----------------------------------------------------------------------
!     Input Namelists
!         &beams3d_input
!            nr             Number of radial gridpoints
!            nphi           Number of toroidal gridpoints
!            nz             Number of vertical gridpoints
!            rmin           Minimum radial extent of grid [m]
!            rmax           Maximum radial extent of grid [m]
!            phimin         Minimum toroidal extent of grid [radians]
!            phimax         Maximum toroidal extent of grid [radians]
!            zmin           Minimum vertical extent of grid [m]
!            zmax           Maximum vertical extent of grid [m]
!            mu             Diffusion coefficient
!            r_start        Radial starting locations for fieldlines [m]
!            phi_start      Toroidal starting locations for fieldlines [radians]
!            phi_end        Toroidal ending locations for fieldlines [radians]
!            z_start        Vertical starting locations for fieldlines [m]
!            npoinc         Number of points (per field period) in which data is saved
!            dphi           Fieldlines following stepsize [radians]
!            follow_tol     Tollerance for fieldline following (LSODE and NAG)
!            vc_adapt_tol   Tollerance for adaptive integration using Virtual casing
!                           (note set to negative value to use non-adaptive integration)
!            int_type       Field line integration method
!                           'NAG','LSODE','RKH68'
!            plasma_mass    Mean plasma mass in [kg]
!            plasma_Zavg    <Z> = sum(n_k*Z_k^2)/sum(n_k*Z_k)
!            plasma_Zmean   [Z] = sum(n_k*Z_k^2*(m_k/plasma_mass))/sum(n_k*Z_k)
!
!            NOTE:  Some grid parameters may be overriden (such as
!                   phimin and phimax) to properly represent a given
!                   field period.
!-----------------------------------------------------------------------
      NAMELIST /beams3d_input/ nr, nphi, nz, rmin, rmax, zmin, zmax, &
                               phimin, phimax, nparticles_start, &
                               r_start_in, phi_start_in, z_start_in, &
                               vll_start_in, npoinc, follow_tol, &
                               t_end_in, mu_start_in, charge_in, &
                               mass_in, Zatom_in, vc_adapt_tol,  &
                               int_type, Adist_beams, Asize_beams, &
                               Div_beams, E_beams, Dex_beams, &
                               mass_beams, charge_beams, Zatom_beams, &
                               r_beams, z_beams, phi_beams, TE_AUX_S, &
                               TE_AUX_F, NE_AUX_S, NE_AUX_F, TI_AUX_S, &
                               TI_AUX_F, POT_AUX_S, POT_AUX_F, &
                               NI_AUX_S, NI_AUX_F, NI_AUX_Z, NI_AUX_M, &
                               ZEFF_AUX_S, ZEFF_AUX_F, P_beams, &
                               ldebug, ne_scale, te_scale, ti_scale, &
                               zeff_scale, plasma_mass, plasma_Zavg, &
                               plasma_Zmean, therm_factor, &
                               fusion_scale, nrho_dist, ntheta_dist, & 
                               nzeta_dist, nvpara_dist, nvperp_dist, &
                               partvmax, lendt_m, te_col_min, &
                               B_kick_min, B_kick_max, freq_kick, E_kick, &
                               rmin_fida, rmax_fida, zmin_fida, &
                               zmax_fida,phimin_fida, phimax_fida, &
                               raxis_fida, zaxis_fida, phiaxis_fida, &
                               nr_fida, nphi_fida, nz_fida, nenergy_fida, &
                               npitch_fida, energy_fida, pitch_fida
      
!-----------------------------------------------------------------------
!     Subroutines
!         read_beams3d_input:   Reads beams3d_input namelist
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE read_beams3d_input(filename, istat)
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(out) :: istat
      LOGICAL :: lexist
      INTEGER :: iunit, local_master, i1
      CHARACTER(LEN=1000) :: line
      ! Initializations
      local_master = 0
      nr     = 101
      nphi   = 360
      nz     = 101
      rmin   =  0.0_rprec
      rmax   =  1.0_rprec
      zmin   = -1.0_rprec
      zmax   =  1.0_rprec
      phimin =  0.0_rprec
      phimax =  pi2
      nparticles_start = 10

      r_start_in   = -1.0
      z_start_in   = -1.0
      phi_start_in = -1.0
      vll_start_in = -1.0
      t_end_in     = -1.0
      mu_start_in  = -1.0
      mass_in      = -1.0
      charge_in    = -1.0
      Zatom_in     = -1.0

      Adist_beams = 1.0_rprec
      Asize_beams = -1.0_rprec
      Div_beams = 1.0_rprec
      Dex_beams = -1
      r_beams = 1.0_rprec
      z_beams = 0.0_rprec
      phi_beams = 0.0_rprec
      E_beams = 0.0_rprec
      mass_beams = 1.0_rprec
      charge_beams = 0.0_rprec
      Zatom_beams = 1.0_rprec
      P_beams = 0.0_rprec
      TE_AUX_S = -1
      TE_AUX_F = -1
      NE_AUX_S = -1
      NE_AUX_F = -1
      TI_AUX_S = -1
      TI_AUX_F = -1
      ZEFF_AUX_S = -1
      ZEFF_AUX_F = -1
      POT_AUX_S = -1
      POT_AUX_F = -1
      NI_AUX_S = -1
      NI_AUX_F = 0
      NI_AUX_Z = 0
      NI_AUX_M = 0
      npoinc = 1
      follow_tol   = 1.0D-7
      vc_adapt_tol = 1.0D-5
      int_type = "LSODE"
      ldebug = .false.
      ne_scale = 1.0
      te_scale = 1.0
      ti_scale = 1.0
      zeff_scale = 1.0
      fusion_scale = 1.0
      plasma_Zmean = 1.0
      plasma_Zavg  = 1.0
      plasma_mass = 1.6726219E-27 ! Assume Hydrogen
      therm_factor = 1.5 ! Factor at which to thermalize particles
      lendt_m = 0.05 ! Max distance a particle travels
      te_col_min = 10 ! Min electron temperature to consider in collisions

      ! Kick model defaults
      B_kick_min = -1.0 ! T
      B_kick_max = 0.0 ! T
      freq_kick = 38.5E6 ! Hz
      E_kick = 100 !V/m

      ! Distribution Function Defaults
      nrho_dist = 64
      ntheta_dist=8
      nzeta_dist=4
      nvpara_dist=32
      nvperp_dist=16
      partvmax = 0 ! Allows user to set value

      !FIDASIM defaults
      rmin_fida = 0.0
      zmin_fida = 0.0
      phimin_fida = 0.0
      rmax_fida = 0.0
      zmax_fida = 0.0
      phimax_fida = 0.0
      nr_fida = 0
      nphi_fida = 0
      nz_fida = 0
      nenergy_fida = 0
      npitch_fida = 0


      ! Read namelist
!      IF (ithread == local_master) THEN
         istat=0
         iunit=12
         INQUIRE(FILE=TRIM(filename),EXIST=lexist)
         IF (.not.lexist) stop 'Could not find input file'
         CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
         IF (istat /= 0) CALL handle_err(NAMELIST_READ_ERR,'beams3d_input in: input.'//TRIM(id_string),istat)
         READ(iunit,NML=beams3d_input,IOSTAT=istat)
         IF (istat /= 0) THEN
            backspace(iunit)
            read(iunit,fmt='(A)') line
            write(6,'(A)') 'Invalid line in namelist: '//TRIM(line)
            CALL handle_err(NAMELIST_READ_ERR,'beams3d_input in: input.'//TRIM(id_string),istat)
         END IF
         CLOSE(iunit)

         ! Update dist function sizes
         ns_prof1=nrho_dist
         ns_prof2=ntheta_dist
         ns_prof3=nzeta_dist
         ns_prof4=nvpara_dist
         ns_prof5=nvperp_dist

         NE_AUX_F = NE_AUX_F*ne_scale
         TE_AUX_F = TE_AUX_F*te_scale
         TI_AUX_F = TI_AUX_F*ti_scale
         ZEFF_AUX_F = ZEFF_AUX_F*zeff_scale
         lbeam = .true.; lkick = .false.
         IF (r_start_in(1) /= -1.0) lbeam = .false.
         IF (lfusion .or. lrestart_particles) lbeam = .false.
         IF (lbbnbi) lbeam = .true.
         IF (lbeam) lcollision = .true.
         IF (B_kick_min >=0 ) lkick = .true.
         nbeams = 0
         DO WHILE ((Asize_beams(nbeams+1) >= 0.0).and.(nbeams<MAXBEAMS))
            nbeams = nbeams + 1
         END DO
         IF (lbbnbi) THEN
            nbeams = 0
            DO WHILE ((Dex_beams(nbeams+1) > 0).and.(nbeams<MAXBEAMS))
               nbeams = nbeams + 1
            END DO
            IF (nbeams == 0)  CALL handle_err(BAD_BEAMDEX_ERR,'beams3d_input in: input.'//TRIM(id_string),nbeams)
         END IF
         IF (lfusion) THEN
            r_start_in = -1
            nbeams = 4
            IF (lfusion_alpha) nbeams = 1
         END IF
         nte = 0
         DO WHILE ((TE_AUX_S(nte+1) >= 0.0).and.(nte<MAXPROFLEN))
            nte = nte + 1
         END DO
         nne = 0
         DO WHILE ((NE_AUX_S(nne+1) >= 0.0).and.(nne<MAXPROFLEN))
            nne = nne + 1
         END DO
         nti = 0
         DO WHILE ((TI_AUX_S(nti+1) >= 0.0).and.(nti<MAXPROFLEN))
            nti = nti + 1
         END DO
         nzeff = 0
         DO WHILE ((ZEFF_AUX_S(nzeff+1) >= 0.0).and.(nzeff<MAXPROFLEN))
            nzeff = nzeff + 1
         END DO
         npot = 0
         DO WHILE ((POT_AUX_S(npot+1) >= 0.0).and.(npot<MAXPROFLEN))
            npot = npot + 1
         END DO

         ! Handle multiple ion species
         IF (ANY(NI_AUX_S >0)) THEN
            nzeff = 0
            DO WHILE ((NI_AUX_S(nzeff+1) >= 0.0).and.(nzeff<MAXPROFLEN))
               nzeff = nzeff + 1
            END DO
            ! Now calc Zeff(1)
            DO i1 = 1, nzeff
               ZEFF_AUX_S(i1) = NI_AUX_S(i1)
               temp = SUM(NI_AUX_F(:,i1)*NI_AUX_Z(:))
               IF (temp > 0) THEN
                  ZEFF_AUX_F(i1) = MAX(SUM(NI_AUX_F(:,i1)*NI_AUX_Z(:)*NI_AUX_Z(:))/temp,1.0)
               ELSE
                  ZEFF_AUX_F(i1) = 1
               END IF
            END DO
            plasma_mass = SUM(NI_AUX_F(:,1)*NI_AUX_M*NI_AUX_M)/(SUM(NI_AUX_F(:,1)*NI_AUX_M))
            plasma_Zavg = SUM(NI_AUX_F(:,1)*NI_AUX_Z*NI_AUX_Z)/(SUM(NI_AUX_F(:,1)*NI_AUX_Z)) ! Note this is just Zeff
            plasma_Zmean = SUM(NI_AUX_F(:,1)*NI_AUX_Z*NI_AUX_Z*NI_AUX_M)/(SUM(NI_AUX_F(:,1)*NI_AUX_Z)*plasma_mass)
         ELSEIF (lfusion) THEN ! Assume 50/50 D T
            nzeff=nne
            NI_AUX_S = NE_AUX_S
            NI_AUX_F(1,:) = 0.5*NE_AUX_F
            NI_AUX_F(2,:) = 0.5*NE_AUX_F
            NI_AUX_M(1) = 3.3435837724E-27;   NI_AUX_Z(1) = 1
            NI_AUX_M(2) = 5.008267217094E-27; NI_AUX_Z(2) = 1 
            ! Now calc Zeff(1)
            DO i1 = 1, nzeff
               ZEFF_AUX_S(i1) = NI_AUX_S(i1)
               temp = SUM(NI_AUX_F(:,i1)*NI_AUX_Z(:))
               IF (temp > 0) THEN
                  ZEFF_AUX_F(i1) = MAX(SUM(NI_AUX_F(:,i1)*NI_AUX_Z(:)*NI_AUX_Z(:))/temp,1.0)
               ELSE
                  ZEFF_AUX_F(i1) = 1
               END IF
            END DO
            plasma_mass = SUM(NI_AUX_F(:,1)*NI_AUX_M*NI_AUX_M)/(SUM(NI_AUX_F(:,1)*NI_AUX_M))
            plasma_Zavg = SUM(NI_AUX_F(:,1)*NI_AUX_Z*NI_AUX_Z)/(SUM(NI_AUX_F(:,1)*NI_AUX_Z)) ! Note this is just Zeff
            plasma_Zmean = SUM(NI_AUX_F(:,1)*NI_AUX_Z*NI_AUX_Z*NI_AUX_M)/(SUM(NI_AUX_F(:,1)*NI_AUX_Z)*plasma_mass)
         ELSEIF (nne > 0) THEN ! Ni=Ne, Z=Zeff
            nzeff = nne
            NI_AUX_S = NE_AUX_S
            NI_AUX_F(1,:) = NE_AUX_F ! NI=NE
            NI_AUX_Z(1) = NINT(plasma_Zavg)
            NI_AUX_M(1) = plasma_mass
            DO i1 = 1, nzeff
               ZEFF_AUX_S(i1) = NI_AUX_S(i1)
               ZEFF_AUX_F(i1) = SUM(NI_AUX_F(:,i1)*NI_AUX_Z(:)*NI_AUX_Z(:))/SUM(NI_AUX_F(:,i1)*NI_AUX_Z(:))
            END DO
         ELSE
            nzeff = 6
            ZEFF_AUX_S(1:6) = (/0.0,0.2,0.4,0.6,0.8,1.0/)
            ZEFF_AUX_F(1:6) = (/1.0,1.0,1.0,1.0,1.0,1.0/)
            NI_AUX_S(1:6)   = (/0.0,0.2,0.4,0.6,0.8,1.0/)
            NI_AUX_F(:,1:6) = 0
            NI_AUX_Z(1) = NINT(plasma_Zavg)
            NI_AUX_M(1) = plasma_mass
         END IF

         IF (lfidasim) THEN
            IF (rmin_fida == 0.0) rmin_fida = rmin
            IF (zmin_fida .eq. 0.0) zmin_fida = zmin
            IF (phimin_fida .eq. 0.0) phimin_fida = phimin
            IF (rmax_fida .eq. 0.0) rmax_fida = rmax
            IF (zmax_fida .eq. 0.0) zmax_fida = zmax
            IF (phimax_fida .eq. 0.0) phimax_fida = phimax
            IF (nr_fida .eq. 0) nr_fida = nr
            IF (nphi_fida .eq. 0) nphi_fida = nphi
            IF (nz_fida .eq. 0) nz_fida = nz
            !nenergy_fida = ns_prof4 !should stay this way!
            !npitch_fida = ns_prof5
            IF (nenergy_fida .eq. 0) nenergy_fida = ns_prof4
            IF (npitch_fida .eq. 0) npitch_fida = ns_prof5
         END IF

         nparticles = 0
         DO WHILE ((r_start_in(nparticles+1) >= 0.0).and.(nparticles<MAXPARTICLES))
            nparticles = nparticles + 1
         END DO
!      END IF

#if defined(HDF5_PAR)
      ! Makes sure that NPARTICLES is divisible by the number of processes
      ! Needed for HDF5 parallel writes.
      IF (lbeam .or. lfusion) THEN
         i1 = nparticles_start/nprocs_beams
         IF (i1*nprocs_beams .ne. nparticles_start) THEN
            nparticles_start = (i1+1)*nprocs_beams
         END IF
      END IF
#endif

      END SUBROUTINE read_beams3d_input

      SUBROUTINE write_beams3d_namelist(iunit_out, istat)
      INTEGER, INTENT(in) :: iunit_out
      INTEGER, INTENT(out) :: istat
      INTEGER :: ik, n
      CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outflt  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outexp  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outcmp  = "(2x,A,1X,'=','(',i3,',',i3,')')"
      CHARACTER(LEN=*), PARAMETER :: outstr  = "(2X,A,1X,'=',1X,'''',A,'''')"
      CHARACTER(LEN=*), PARAMETER :: onevar  = "(2X,A,1X,'=',1X,L1,2(2X,A,1X,'=',1X,ES22.12E3))"
      CHARACTER(LEN=*), PARAMETER :: vecvar  = "(2X,A,'(',I3.3,')',1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: vecvar2  = "(2X,A,'(',I3.3,',',I3.3,')',1X,'=',1X,ES22.12E3)"
      istat = 0
      WRITE(iunit_out,'(A)') '&BEAMS3D_INPUT'
      WRITE(iunit_out,'(A)') '!---------- General Parameters ------------'
      WRITE(iunit_out,outint) 'NR',nr
      WRITE(iunit_out,outint) 'NZ',nz
      WRITE(iunit_out,outint) 'NPHI',nphi
      WRITE(iunit_out,outflt) 'RMIN',rmin
      WRITE(iunit_out,outflt) 'RMAX',rmax
      WRITE(iunit_out,outflt) 'ZMIN',zmin
      WRITE(iunit_out,outflt) 'ZMAX',zmax
      WRITE(iunit_out,outflt) 'PHIMIN',phimin
      WRITE(iunit_out,outflt) 'PHIMAX',phimax
      WRITE(iunit_out,outint) 'NPOINC',npoinc
      WRITE(iunit_out,outstr) 'INT_TYPE',TRIM(int_type)
      WRITE(iunit_out,outflt) 'FOLLOW_TOL',follow_tol
      WRITE(iunit_out,outflt) 'VC_ADAPT_TOL',vc_adapt_tol
      WRITE(iunit_out,outint) 'NPARTICLES_START',nparticles_start
      WRITE(iunit_out,'(A)') '!---------- Plasma Parameters ------------'
      WRITE(iunit_out,outflt) 'PLASMA_MASS',plasma_mass
      WRITE(iunit_out,outflt) 'PLASMA_ZAVG',plasma_zavg
      WRITE(iunit_out,outflt) 'PLASMA_ZMEAN',plasma_zmean
      WRITE(iunit_out,outflt) 'THERM_FACTOR',therm_factor
      WRITE(iunit_out,'(A)') '!---------- Distribution Parameters ------------'
      WRITE(iunit_out,outint) 'NRHO_DIST',ns_prof1
      WRITE(iunit_out,outint) 'NTHETA_DIST',ns_prof2
      WRITE(iunit_out,outint) 'NZETA_DIST',ns_prof3
      WRITE(iunit_out,outint) 'NVPARA_DIST',ns_prof4
      WRITE(iunit_out,outint) 'NVPERP_DIST',ns_prof5
      WRITE(iunit_out,outflt) 'PARTVMAX',partvmax
      IF (B_kick_min>0) THEN
         WRITE(iunit_out,'(A)') '!---------- Kick Model Parameters ------------'
         WRITE(iunit_out,outflt) 'E_KICK',E_kick
         WRITE(iunit_out,outflt) 'FREQ_KICK',freq_kick
         WRITE(iunit_out,outflt) 'B_KICK_MIN',B_kick_min
         WRITE(iunit_out,outflt) 'B_KICK_MAX',B_kick_max
      END IF
      IF (lbeam) THEN
         WRITE(iunit_out,"(A)") '!---------- Profiles ------------'
         WRITE(iunit_out,outflt) 'NE_SCALE',NE_SCALE
         WRITE(iunit_out,outflt) 'TE_SCALE',TE_SCALE
         WRITE(iunit_out,outflt) 'TI_SCALE',TI_SCALE
         WRITE(iunit_out,outflt) 'ZEFF_SCALE',ZEFF_SCALE
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'NE_AUX_S',(ne_aux_s(n), n=1,nne)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'NE_AUX_F',(ne_aux_f(n), n=1,nne)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'TE_AUX_S',(te_aux_s(n), n=1,nte)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'TE_AUX_F',(te_aux_f(n), n=1,nte)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'TI_AUX_S',(ti_aux_s(n), n=1,nti)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'TI_AUX_F',(ti_aux_f(n), n=1,nti)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'ZEFF_AUX_S',(zeff_aux_s(n), n=1,nzeff)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'ZEFF_AUX_F',(zeff_aux_f(n), n=1,nzeff)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'POT_AUX_S',(zeff_aux_s(n), n=1,npot)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'POT_AUX_F',(zeff_aux_f(n), n=1,npot)
         DO n = 1, nbeams
            WRITE(iunit_out,"(A,I2.2)") '!---- BEAM #',n
            IF (dex_beams(n)>0) &
               WRITE(iunit_out,vecvar) 'DEX_BEAMS',n,dex_beams(n)
            WRITE(iunit_out,vecvar) 'T_END_IN',n,t_end_in(n)
            WRITE(iunit_out,vecvar) 'DEX_BEAMS',n,dex_beams(n)
            WRITE(iunit_out,vecvar) 'DIV_BEAMS',n,div_beams(n)
            WRITE(iunit_out,vecvar) 'ADIST_BEAMS',n,adist_beams(n)
            WRITE(iunit_out,vecvar) 'ASIZE_BEAMS',n,asize_beams(n)
            WRITE(iunit_out,vecvar) 'MASS_BEAMS',n,mass_beams(n)
            WRITE(iunit_out,vecvar) 'ZATOM_BEAMS',n,zatom_beams(n)
            WRITE(iunit_out,vecvar) 'CHARGE_BEAMS',n,charge_beams(n)
            WRITE(iunit_out,vecvar) 'E_BEAMS',n,e_beams(n)
            WRITE(iunit_out,vecvar) 'P_BEAMS',n,p_beams(n)
            WRITE(iunit_out,vecvar2) 'R_BEAMS',n,1,r_beams(n,1)
            WRITE(iunit_out,vecvar2) 'PHI_BEAMS',n,1,phi_beams(n,1)
            WRITE(iunit_out,vecvar2) 'Z_BEAMS',n,1,z_beams(n,1)
            WRITE(iunit_out,vecvar2) 'R_BEAMS',n,2,r_beams(n,2)
            WRITE(iunit_out,vecvar2) 'PHI_BEAMS',n,2,phi_beams(n,2)
            WRITE(iunit_out,vecvar2) 'Z_BEAMS',n,2,z_beams(n,2)
         END DO
      ELSE
         n = COUNT(r_start_in > 0)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'R_START_IN',(r_start_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'Z_START_IN',(z_start_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'PHI_START_IN',(phi_start_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'VLL_START_IN',(vll_start_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'MU_START_IN',(mu_start_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'MASS_IN',(mass_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'CHARGE_IN',(charge_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'ZATOM_IN',(zatom_in(ik), ik=1,n)
         n = COUNT(t_end_in > -1)
         WRITE(iunit_out,"(2X,A,1X,'=',I0,'*',ES22.12E3)") 'T_END_IN',n,MAXVAL(t_end_in)
      END IF
      WRITE(iunit_out,'(A)') '/'

      END SUBROUTINE write_beams3d_namelist

      SUBROUTINE BCAST_BEAMS3D_INPUT(local_master,comm,istat)
      USE mpi_inc
      IMPLICIT NONE
      
      INTEGER, INTENT(inout) :: comm
      INTEGER, INTENT(in)    :: local_master
      INTEGER, INTENT(inout) :: istat
      IF (istat .ne. 0) RETURN
#if defined(MPI_OPT)
      CALL MPI_BCAST(lbeam, 1, MPI_LOGICAL, local_master, comm,istat)
      CALL MPI_BCAST(nr,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nphi,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nz,1,MPI_INTEGER, local_master, comm,istat)


      CALL MPI_BCAST(ns_prof1,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(ns_prof2,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(ns_prof3,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(ns_prof4,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(ns_prof5,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(partvmax,1,MPI_REAL8, local_master, comm,istat)

      CALL MPI_BCAST(nbeams,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nparticles_start,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nparticles,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(npoinc,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(rmin,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(rmax,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(zmin,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(zmax,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(phimin,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(phimax,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(vc_adapt_tol,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(plasma_mass,1,MPI_REAL8, local_master, comm,istat)

      CALL MPI_BCAST(nte,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nne,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nti,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(TE_AUX_S,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(TE_AUX_F,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(NE_AUX_S,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(NE_AUX_F,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(TI_AUX_S,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(TI_AUX_F,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(te_scale,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(ne_scale,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(ti_scale,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(zeff_scale,1,MPI_REAL8, local_master, comm,istat)

      CALL MPI_BCAST(t_end_in,MAXPARTICLES,MPI_REAL8, local_master, comm,istat)

      IF (lbeam) THEN
          CALL MPI_BCAST(Adist_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(Dex_beams,MAXBEAMS,MPI_INTEGER, local_master, comm,istat)
          CALL MPI_BCAST(Asize_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(Div_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(E_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(r_beams,MAXBEAMS*2,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(z_beams,MAXBEAMS*2,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(phi_beams,MAXBEAMS*2,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(mass_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(charge_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(Zatom_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
      ELSE
          CALL MPI_BCAST(r_start_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(z_start_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(phi_start_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(mu_start_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(vll_start_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(mass_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(charge_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(Zatom_in,nparticles,MPI_REAL8, local_master, comm,istat)
      END IF

      CALL MPI_BCAST(follow_tol,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(int_type, 256, MPI_CHARACTER, local_master, comm,istat)
#endif
      END SUBROUTINE BCAST_BEAMS3D_INPUT

      END MODULE beams3d_input_mod
