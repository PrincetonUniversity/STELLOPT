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
      USE beams3d_globals
      USE safe_open_mod, ONLY: safe_open
      USE mpi_params
      USE mpi_inc

!-----------------------------------------------------------------------
!     Module Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      ! These are helpers to give the ns1_prof variables user friendly names
      INTEGER :: nrho_dist, ntheta_dist, nzeta_dist, nvpara_dist, nvperp_dist, nphi_dist
      REAL(rprec) :: temp
      ! These are helpers for backwards compatibility all values here
      ! will be ignored elsewhere in the code.
      REAL(rprec) :: plasma_zavg ! 
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
!            Zeff           <Z> = sum(n_k*Z_k^2)/sum(n_k*Z_k)
!            plasma_Zmean   [Z] = sum(n_k*Z_k^2*(plasma_mass/m_k))/sum(n_k*Z_k)
!
!            NOTE:  Some grid parameters may be overriden (such as
!                   phimin and phimax) to properly represent a given
!                   field period.
!-----------------------------------------------------------------------
      NAMELIST /beams3d_input/ nr, nphi, nz, rmin, rmax, zmin, zmax, &
                               phimin, phimax, nparticles_start,&
                               r_start_in, phi_start_in, z_start_in, &
                               vll_start_in, npoinc, follow_tol, &
                               t_end_in, mu_start_in, charge_in, &
                               mass_in, Zatom_in, weight_in, &
                               vc_adapt_tol,  &
                               int_type, Adist_beams, Asize_beams, &
                               Div_beams, E_beams, Dex_beams, &
                               mass_beams, charge_beams, Zatom_beams, &
                               r_beams, z_beams, phi_beams, s_max, TE_AUX_S, &
                               TE_AUX_F, NE_AUX_S, NE_AUX_F, TI_AUX_S, &
                               TI_AUX_F, POT_AUX_S, POT_AUX_F, &
                               NI_AUX_S, NI_AUX_F, NI_AUX_Z, NI_AUX_M, &
                               ZEFF_AUX_S, ZEFF_AUX_F, P_beams, &
                               ldebug, ne_scale, te_scale, ti_scale, &
                               zeff_scale, &
                               plasma_zavg, plasma_mass, plasma_Zmean, &
                               therm_factor, fusion_scale, &
                               nrho_dist, ntheta_dist, & 
                               nzeta_dist, nphi_dist, nvpara_dist, nvperp_dist, &
                               partvmax, rho_max_dist, lendt_m, te_col_min, &
                               B_kick_min, B_kick_max, freq_kick, E_kick,&
                               vr_start_in, vphi_start_in, vz_start_in, &
                               rho_fullorbit, duplicate_factor, &
                               B_kick_min, B_kick_max, freq_kick, E_kick, &
                               rmin_fida, rmax_fida, zmin_fida, &
                               zmax_fida,phimin_fida, phimax_fida, &
                               nr_fida, nphi_fida, nz_fida, nenergy_fida, &
                               npitch_fida, t_fida
      
!-----------------------------------------------------------------------
!     Subroutines
!         init_beams3d_input:   Initializes the namelist
!         read_beams3d_input:   Reads beams3d_input namelist
!-----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE init_beams3d_input
      IMPLICIT NONE
      pi2 = 8.0 * ATAN(1.0)
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

      r_start_in    = -1.0
      z_start_in    = -1.0
      phi_start_in  = -1.0
      vll_start_in  = -1.0
      vr_start_in   =  0.0
      vphi_start_in =  0.0
      vz_start_in   =  0.0
      t_end_in      = -1.0
      mu_start_in   = -1.0
      mass_in       = -1.0
      charge_in     = -1.0
      Zatom_in      = -1.0
      weight_in     =  1.0

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
      s_max = 1.0_rprec
      s_max_te = 0.0_rprec
      s_max_ti = 0.0_rprec
      s_max_ne = 0.0_rprec
      s_max_zeff = 0.0_rprec
      s_max_pot = 0.0_rprec
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
      dexionT = 1
      dexionD = 2
      npoinc = 1
      follow_tol   = 1.0D-9
      vc_adapt_tol = 1.0D-5
      int_type = "LSODE"
      ldebug = .false.
      ne_scale = 1.0
      te_scale = 1.0
      ti_scale = 1.0
      zeff_scale = 1.0
      fusion_scale = 1.0
      plasma_Zmean = 1.0
      plasma_mass = 1.6726219E-27 ! Assume Hydrogen
      therm_factor = 1.5 ! Factor at which to thermalize particles
      lendt_m = 0.05 ! Max distance a particle travels
      te_col_min = 10 ! Min electron temperature to consider in collisions

      ! Kick model defaults
      B_kick_min = -1.0 ! T
      B_kick_max = 0.0 ! T
      freq_kick = 38.5E6 ! Hz
      E_kick = 100 !V/m

      ! Full Oribt model
      rho_fullorbit = 1.0E10 ! Default to off
      duplicate_factor = 1 ! No particle duplication

      ! Distribution Function Defaults
      nrho_dist = 64
      ntheta_dist=8
      nzeta_dist = -1  ! Kept for historical reasons
      nphi_dist  = 4
      nvpara_dist=32
      nvperp_dist=16
      partvmax = 0 ! Allows user to set value
      rho_max_dist = -1

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
      t_fida = 0.0
      RETURN
      END SUBROUTINE init_beams3d_input
      
      SUBROUTINE read_beams3d_input(filename, istat)
         IMPLICIT NONE
         CHARACTER(*), INTENT(in) :: filename
         INTEGER, INTENT(out) :: istat
         LOGICAL :: lexist
         INTEGER :: iunit, local_master, i1, ik
         CHARACTER(LEN=1000) :: line
      ! Initializations
      local_master = 0


      ! Read namelist
         IF (filename /= 'IMAS') THEN
            istat=0
            iunit=12
            INQUIRE(FILE=TRIM(filename),EXIST=lexist)
            IF (.not.lexist) stop 'Could not find input file'
            CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
            IF (istat /= 0) THEN
               WRITE(6,'(A)') 'ERROR opening file: ',TRIM(filename)
               CALL FLUSH(6)
               STOP
            END IF
            READ(iunit,NML=beams3d_input,IOSTAT=istat)
            IF (istat /= 0) THEN
               WRITE(6,'(A)') 'ERROR reading namelist BEAMS3D_INPUT from file: ',TRIM(filename)
               backspace(iunit)
               read(iunit,fmt='(A)') line
               write(6,'(A)') 'Invalid line in namelist: '//TRIM(line)
               CALL FLUSH(6)
               STOP
            END IF
            CLOSE(iunit)
         END IF

         ! Update dist function sizes
         ns_prof1=nrho_dist
         ns_prof2=ntheta_dist
         ns_prof3=MAX(nzeta_dist,nphi_dist)
         ns_prof4=nvpara_dist
         ns_prof5=nvperp_dist

         ! Fix s_max_dist
         IF (rho_max_dist < 0) THEN
            rho_max_dist = DBLE(ns_prof1)/DBLE(ns_prof1-1)
         END IF

         NE_AUX_F = NE_AUX_F*ne_scale
         TE_AUX_F = TE_AUX_F*te_scale
         TI_AUX_F = TI_AUX_F*ti_scale
         ZEFF_AUX_F = ZEFF_AUX_F*zeff_scale
         lbeam = .true.; lkick = .false.; lgcsim = .true.
         
         IF (r_start_in(1) /= -1.0) lbeam = .false.
         IF (lfusion .or. lrestart_particles) lbeam = .false.
         IF (lbbnbi) lbeam = .true.
         IF (lbeam) lcollision = .true.
         IF (B_kick_min >=0 ) lkick = .true.
         nbeams = 0
         DO ik = 1, MAXBEAMS
            IF (Asize_beams(ik) >= 0.0) nbeams = nbeams + 1
         END DO
         IF (lbbnbi) THEN
            nbeams = 0
            DO ik = 1, MAXBEAMS
               IF (Dex_beams(ik) > 0) nbeams = nbeams + 1
            END DO
            IF (nbeams == 0) THEN
               WRITE(6,'(A)') 'BEAMLET beam model requested but nbeams==0'
               WRITE(6,'(A)') '  Check DEX_BEAMS is set in BEAMS3D_INPUT'
               CALL FLUSH(6)
               STOP
            END IF
         END IF
         IF (lfusion) THEN
            r_start_in = -1
            nbeams = 0
            IF (lfusion_alpha) nbeams = nbeams + 1
            IF (lfusion_tritium) nbeams = nbeams + 1
            IF (lfusion_proton) nbeams = nbeams + 1
            IF (lfusion_He3) nbeams = nbeams + 1
         END IF
         nte = 0
         DO ik = 1, MAXPROFLEN
            IF (TE_AUX_S(ik) >= 0.0) nte = nte+1
         END DO
         IF (nte > 0) s_max_te = TE_AUX_S(nte)
         nne = 0
         DO ik = 1, MAXPROFLEN
            IF (NE_AUX_S(ik) >= 0.0) nne = nne+1
         END DO
         IF (nne > 0) s_max_ne = NE_AUX_S(nne)
         nti = 0
         DO ik = 1, MAXPROFLEN
            IF (TI_AUX_S(ik) >= 0.0) nti = nti+1
         END DO
         IF (nti > 0) s_max_ti = NE_AUX_S(nti)
         nzeff = 0
         DO ik = 1, MAXPROFLEN
            IF (ZEFF_AUX_S(ik) >= 0.0) nzeff = nzeff+1
         END DO
         IF (nzeff > 0) s_max_zeff=ZEFF_AUX_S(nzeff)
         npot = 0
         DO ik = 1, MAXPROFLEN
            IF (POT_AUX_S(ik) >= 0.0) npot = npot+1
         END DO
         IF (npot > 0)  s_max_pot = POT_AUX_S(npot)
         ! Handle multiple ion species
         IF (ANY(NI_AUX_S >0)) THEN
            nzeff = 0
            DO ik = 1, MAXPROFLEN
               IF (NI_AUX_S(ik) >= 0.0) nzeff = nzeff+1
            END DO
            IF (nzeff > 0) s_max_zeff=ZEFF_AUX_S(nzeff)
            ! Now calc Zeff(1)
            DO ik = 1, nzeff
               ZEFF_AUX_S(ik) = NI_AUX_S(ik)
               temp = SUM(NI_AUX_F(:,ik)*NI_AUX_Z(:))
               IF (temp > 0) THEN
                  ZEFF_AUX_F(ik) = MAX(SUM(NI_AUX_F(:,ik)*NI_AUX_Z(:)*NI_AUX_Z(:))/temp,1.0)
               ELSE
                  ZEFF_AUX_F(ik) = 1
               END IF
            END DO
            plasma_mass = SUM(NI_AUX_F(:,1)*NI_AUX_M*NI_AUX_M)/(SUM(NI_AUX_F(:,1)*NI_AUX_M))
            plasma_Zmean = SUM(NI_AUX_F(:,1)*NI_AUX_Z*NI_AUX_Z*plasma_mass/NI_AUX_M,DIM=1,MASK=(NI_AUX_M>1E-27))/(SUM(NI_AUX_F(:,1)*NI_AUX_Z))
            ! Set indices for T and D
            DO ik = 1, NION
               IF ((NI_AUX_Z(ik) == 1) .and. (NINT(NI_AUX_M(ik)*6.02214076208E+26) == 3)) dexionT = ik
               IF ((NI_AUX_Z(ik) == 1) .and. (NINT(NI_AUX_M(ik)*6.02214076208E+26) == 2)) dexionD = ik
               IF ((NI_AUX_Z(ik) == 2) .and. (NINT(NI_AUX_M(ik)*6.02214076208E+26) == 3)) dexionHe3 = ik
            END DO
            !WRITE(6,*) ' Tritium index: ',dexionT
            !WRITE(6,*) ' Deuturium index: ',dexionD
            s_max_zeff=ZEFF_AUX_S(nzeff+1)
         ELSEIF (lfusion) THEN ! Assume 50/50 D T
            nzeff=nne
            NI_AUX_S = NE_AUX_S
            NI_AUX_F(1,:) = 0.5*NE_AUX_F
            NI_AUX_F(2,:) = 0.5*NE_AUX_F
            NI_AUX_M(1) = 3.3435837724E-27;   NI_AUX_Z(1) = 1
            NI_AUX_M(2) = 5.008267217094E-27; NI_AUX_Z(2) = 1 
            ! Now calc Zeff(1)
            DO ik = 1, nzeff
               ZEFF_AUX_S(ik) = NI_AUX_S(ik)
               temp = SUM(NI_AUX_F(:,ik)*NI_AUX_Z(:))
               IF (temp > 0) THEN
                  ZEFF_AUX_F(ik) = MAX(SUM(NI_AUX_F(:,ik)*NI_AUX_Z(:)*NI_AUX_Z(:))/temp,1.0)
               ELSE
                  ZEFF_AUX_F(ik) = 1
               END IF
            END DO
            plasma_mass = SUM(NI_AUX_F(:,1)*NI_AUX_M*NI_AUX_M)/(SUM(NI_AUX_F(:,1)*NI_AUX_M))
            plasma_Zmean = SUM(NI_AUX_F(:,1)*NI_AUX_Z*NI_AUX_Z*plasma_mass/NI_AUX_M,DIM=1,MASK=(NI_AUX_M>1E-27))/(SUM(NI_AUX_F(:,1)*NI_AUX_Z))
         ELSEIF (nne > 0) THEN ! Ni=Ne, Z=Zeff
            NI_AUX_Z(1) = 1 ! Assume Hydrogen Plasma
            NI_AUX_M(1) = plasma_mass
            NI_AUX_S = NE_AUX_S
            NI_AUX_F(1,:) = NE_AUX_F
            ! First check if user provided ZEFF
            IF (.not. ANY(ZEFF_AUX_S >0)) THEN
               ! NI=NE
               ! Default ZEFF_AUX_S
               nzeff = 6
               ZEFF_AUX_S(1:6) = (/0.0,0.2,0.4,0.6,0.8,1.0/)
               ZEFF_AUX_F(1:6) = (/1.0,1.0,1.0,1.0,1.0,1.0/)
            END IF
         ELSE
            nzeff = 6
            ZEFF_AUX_S(1:6) = (/0.0,0.2,0.4,0.6,0.8,1.0/)
            ZEFF_AUX_F(1:6) = (/1.0,1.0,1.0,1.0,1.0,1.0/)
            NI_AUX_S(1:6)   = (/0.0,0.2,0.4,0.6,0.8,1.0/)
            NI_AUX_F(:,1:6) = 0
            NI_AUX_Z(1) = 1
            NI_AUX_M(1) = plasma_mass
         END IF
         s_max_zeff=ZEFF_AUX_S(nzeff)

         nparticles = 0
         DO ik = 1, MAXPARTICLES
            IF (r_start_in(ik) >= 0.0) nparticles = nparticles + 1
         END DO

#if !defined(NAG)
      IF (int_type=='NAG') THEN
         int_type = 'LSODE'
         IF (lverb) THEN
            WRITE(6,*) '======================================='
            WRITE(6,*) '  INT_TYPE = NAG in input but BEAMS3D'
            WRITE(6,*) '  is not linked to NAG library.'
            WRITE(6,*) '  Using INT_TYPE = LSODE instead.'
            WRITE(6,*) '======================================='
         END IF
      END IF
#endif

      END SUBROUTINE read_beams3d_input

      SUBROUTINE write_beams3d_namelist(iunit_out, istat)
      INTEGER, INTENT(in) :: iunit_out
      INTEGER, INTENT(out) :: istat
      INTEGER :: ik, n, l
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
      WRITE(iunit_out,'(A)') '!---------- Background Grid Parameters ------------'
      WRITE(iunit_out,outint) 'NR',nr
      WRITE(iunit_out,outint) 'NZ',nz
      WRITE(iunit_out,outint) 'NPHI',nphi
      WRITE(iunit_out,outflt) 'RMIN',rmin
      WRITE(iunit_out,outflt) 'RMAX',rmax
      WRITE(iunit_out,outflt) 'ZMIN',zmin
      WRITE(iunit_out,outflt) 'ZMAX',zmax
      WRITE(iunit_out,outflt) 'PHIMIN',phimin
      WRITE(iunit_out,outflt) 'PHIMAX',phimax
      WRITE(iunit_out,outflt) 'VC_ADAPT_TOL',vc_adapt_tol
      WRITE(iunit_out,'(A)') '!---------- Marker Tracking Parameters ------------'
      WRITE(iunit_out,outstr) 'INT_TYPE',TRIM(int_type)
      WRITE(iunit_out,outflt) 'FOLLOW_TOL',follow_tol
      WRITE(iunit_out,outint) 'NPOINC',npoinc
      WRITE(iunit_out,outint) 'NPARTICLES_START',nparticles_start
      WRITE(iunit_out,outflt) 'LENDT_M',lendt_m
      WRITE(iunit_out,outflt) 'RHO_FULLORBIT',rho_fullorbit
      WRITE(iunit_out,outint) 'DUPLICATE_FACTOR',duplicate_factor
      WRITE(iunit_out,'(A)') '!---------- Distribution Parameters ------------'
      WRITE(iunit_out,outint) 'NRHO_DIST',ns_prof1
      WRITE(iunit_out,outint) 'NTHETA_DIST',ns_prof2
      WRITE(iunit_out,outint) 'NPHI_DIST',ns_prof3
      WRITE(iunit_out,outint) 'NVPARA_DIST',ns_prof4
      WRITE(iunit_out,outint) 'NVPERP_DIST',ns_prof5
      WRITE(iunit_out,outflt) 'PARTVMAX',partvmax
      WRITE(iunit_out,outflt) 'RHO_MAX_DIST',rho_max_dist
      IF (B_kick_min>0) THEN
         WRITE(iunit_out,'(A)') '!---------- Kick Model Parameters ------------'
         WRITE(iunit_out,outflt) 'E_KICK',E_kick
         WRITE(iunit_out,outflt) 'FREQ_KICK',freq_kick
         WRITE(iunit_out,outflt) 'B_KICK_MIN',B_kick_min
         WRITE(iunit_out,outflt) 'B_KICK_MAX',B_kick_max
      END IF
      WRITE(iunit_out,'(A)') '!---------- Plasma Parameters ------------'
      WRITE(iunit_out,outflt) 'PLASMA_MASS',plasma_mass
      WRITE(iunit_out,outflt) 'PLASMA_ZMEAN',plasma_zmean
      WRITE(iunit_out,outflt) 'THERM_FACTOR',therm_factor
      WRITE(iunit_out,"(A)") '!---------- Profiles ---------------------'
      WRITE(iunit_out,outflt) 'NE_SCALE',NE_SCALE
      WRITE(iunit_out,outflt) 'TE_SCALE',TE_SCALE
      WRITE(iunit_out,outflt) 'TI_SCALE',TI_SCALE
      WRITE(iunit_out,outflt) 'ZEFF_SCALE',ZEFF_SCALE
      WRITE(iunit_out,outflt) 'THERM_FACTOR',therm_factor
      ik = COUNT(ne_aux_s >= 0)
      IF (ik > 0) THEN
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'NE_AUX_S',(ne_aux_s(n), n=1,ik)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'NE_AUX_F',(ne_aux_f(n), n=1,ik)
      END IF
      ik = COUNT(te_aux_s >= 0)
      IF (ik > 0) THEN
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'TE_AUX_S',(te_aux_s(n), n=1,ik)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'TE_AUX_F',(te_aux_f(n), n=1,ik)
      END IF
      ik = COUNT(ti_aux_s >= 0)
      IF (ik > 0) THEN
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'TI_AUX_S',(ti_aux_s(n), n=1,ik)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'TI_AUX_F',(ti_aux_f(n), n=1,ik)
      END IF
      ik = COUNT(ni_aux_s >= 0)
      IF (ik > 0) THEN
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'NI_AUX_M',(ni_aux_m(n), n=1,NION)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,I0))") 'NI_AUX_Z',(ni_aux_z(n), n=1,NION)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'NI_AUX_S',(ni_aux_s(n), n=1,ik)
         DO n = 1, NION
            IF (ANY(NI_AUX_F(n,:)>0)) THEN
               WRITE(iunit_out,"(2X,A,'(',I1.1,',:) =',4(1X,ES22.12E3))") 'NI_AUX_F',n,(ni_aux_f(n,l), l=1,ik)
            END IF
         END DO
      END IF
      ik = COUNT(zeff_aux_s >= 0)
      IF (ik > 0) THEN
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'ZEFF_AUX_S',(zeff_aux_s(n), n=1,ik)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'ZEFF_AUX_F',(zeff_aux_f(n), n=1,ik)
      END IF
      ik = COUNT(pot_aux_s >= 0)
      IF (ik > 0) THEN
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'POT_AUX_S',(pot_aux_s(n), n=1,ik)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'POT_AUX_F',(pot_aux_f(n), n=1,ik)
      END IF
      IF (lbeam) THEN
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
         WRITE(iunit_out,"(A,I2.2)") '!---------- Markers ----------------------'
         n = COUNT(r_start_in > 0)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'R_START_IN',(r_start_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'Z_START_IN',(z_start_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'PHI_START_IN',(phi_start_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'VLL_START_IN',(vll_start_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'MU_START_IN',(mu_start_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'MASS_IN',(mass_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'CHARGE_IN',(charge_in(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'ZATOM_IN',(zatom_in(ik), ik=1,n)
         IF (ANY(weight_in /= 1)) WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'WEIGHT_IN',(weight_in(ik), ik=1,n)
         n = COUNT(t_end_in > -1)
         WRITE(iunit_out,"(2X,A,1X,'=',I6,'*',ES19.12E3)") 'T_END_IN',n,MAXVAL(t_end_in)
      END IF
      WRITE(iunit_out,'(A)') '/'

      END SUBROUTINE write_beams3d_namelist

      SUBROUTINE write_beams3d_namelist_byfile(filename)
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER :: iunit, istat
      LOGICAL :: lexists
      
      iunit = 100
      istat = 0
      INQUIRE(FILE=TRIM(filename),exist=lexists)
      IF (lexists) THEN
         OPEN(unit=iunit, file=TRIM(filename), iostat=istat, status="old", position="append")
      ELSE
         OPEN(unit=iunit, file=TRIM(filename), iostat=istat, status="new")
      END IF
      IF (istat .ne. 0) RETURN
      CALL write_beams3d_namelist(iunit,istat)
      CLOSE(iunit)

      RETURN
      END SUBROUTINE write_beams3d_namelist_byfile

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
      CALL MPI_BCAST(lendt_m,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(rho_fullorbit,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(duplicate_factor,1,MPI_INTEGER, local_master, comm,istat)

      CALL MPI_BCAST(nte,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nne,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nti,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(TE_AUX_S,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(TE_AUX_F,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(NE_AUX_S,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(NE_AUX_F,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(TI_AUX_S,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(TI_AUX_F,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(NI_AUX_S,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(NI_AUX_F,MAXPROFLEN*NION,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(te_scale,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(ne_scale,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(ti_scale,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(zeff_scale,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(fusion_scale,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(therm_factor,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(s_max,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(s_max_ne,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(s_max_te,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(s_max_ti,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(s_max_zeff,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(s_max_pot,1,MPI_REAL8, local_master, comm,istat)

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
          CALL MPI_BCAST(weight_in,nparticles,MPI_REAL8, local_master, comm,istat)
      END IF

      CALL MPI_BCAST(follow_tol,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(int_type, 256, MPI_CHARACTER, local_master, comm,istat)

      CALL MPI_BCAST(E_kick,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(freq_kick,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(B_kick_min,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(B_kick_max,1,MPI_REAL8, local_master, comm,istat)
#endif
      END SUBROUTINE BCAST_BEAMS3D_INPUT

      END MODULE beams3d_input_mod
