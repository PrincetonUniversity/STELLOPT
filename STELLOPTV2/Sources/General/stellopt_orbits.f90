!-----------------------------------------------------------------------
!     Subroutine:    stellopt_orbits
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/29/2014
!     Description:   This subroutine calculates the loss of particles
!                    from an equilibrium using the BEAMS3D code.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_orbits(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime, ONLY:  proc_string, bigno, rprec, pi2
      USE equil_utils, ONLY: get_equil_RZ, get_equil_Bflx,&
                             get_equil_ne, get_equil_te, get_equil_ti
      USE equil_vals, ONLY: nfp, rho, orbit_lost_frac
      USE stellopt_targets, ONLY: sigma_orbit, vll_orbit, mu_orbit,&
            nu_orbit, nv_orbit, nsd, np_orbit, mass_orbit, Z_orbit,&
            vperp_orbit, nsd
      USE stellopt_vars, ONLY: ne_type, te_type, ti_type, ne_norm
      USE read_wout_mod, ONLY:  rmax_surf, rmin_surf, zmax_surf
!DEC$ IF DEFINED (BEAMS3D_OPT)
      ! BEAMS3D Libraries
      USE beams3d_runtime, ONLY: nparticles_start, &
            vll_start_in, R_start_in, Z_start_in, PHI_start_in, mu_start_in, &
            mu_start_in, charge_in, mass_in, t_end_in, Zatom_in, &
            TE_AUX_S_BEAMS => TE_AUX_S, TE_AUX_F_BEAMS => TE_AUX_F, &
            NE_AUX_S_BEAMS => NE_AUX_S, NE_AUX_F_BEAMS => NE_AUX_F, &
            TI_AUX_S_BEAMS => TI_AUX_S, TI_AUX_F_BEAMS => TI_AUX_F, &
            BEAMS3D_VERSION
      USE beams3d_grid, ONLY: nte, nne, nti, rmin, rmax, zmin, zmax, &
                              phimin, phimax
      USE beams3d_lines, ONLY: lost_lines
!DEC$ ENDIF
!DEC$ IF DEFINED (MPI_OPT)
      USE mpi_params     
!DEC$ ENDIF
      !USE safe_open_mod
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(inout) :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER :: dex, ik, nparts, u, v, p
      DOUBLE PRECISION :: s_val, u_val, v_val, Ro, Zo, s_min, s_max, tf,&
                     Bs, Bu, Bv, modb
      INTEGER, PARAMETER :: NBEAM_PROF = 50
      CHARACTER(256) :: beam_str
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
!DEC$ IF DEFINED (BEAMS3D_OPT)
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  ORBIT CALCULATION  -------------------------'

      ! Initialize the grid
      s_val = (rmax_surf - rmin_surf)/20 ! 10% scaling of width
      rmin = rmin_surf - s_val
      rmax = rmax_surf + s_val
      zmin = -zmax_surf - s_val
      zmax = zmax_surf + s_val
      phimin = 0
      phimax = pi2/nfp



      ! Initialize the particle starting points
      nparts = 1
      s_min = 1
      s_max = 0
      tf    = MAXVAL(t_end_in)
      DO ik = 1, nsd
         IF (sigma_orbit(ik) .ge. bigno) CYCLE
         DO v = 1, nv_orbit
            DO u = 1, nu_orbit
               s_val = rho(ik)
               s_min = MIN(s_min,s_val)
               s_max = MAX(s_max,s_val)
               u_val = pi2*(u-1)/nu_orbit
               v_val = pi2*(v-1)/(nv_orbit-1)/2
               CALL get_equil_RZ(s_val,u_val,v_val,Ro,Zo, iflag)
               CALL get_equil_Bflx(s_val,u_val,v_val,Bs,Bu,Bv,iflag,modb)
               DO p = 1, np_orbit
                  R_start_in(nparts) = Ro
                  Z_start_in(nparts) = Zo
                  PHI_start_in(nparts) = v_val/nfp
                  vll_start_in(nparts) = vll_orbit(p)
                  Zatom_in(nparts)  = Z_orbit
                  mass_in(nparts)   = mass_orbit
                  charge_in(nparts) = Z_orbit*1.602176565E-19 ! [C] NIST(http://physics.nist.gov/cgi-bin/cuu/Value?e)
                  !charge_in(nparts) = 2*1.602176565E-19 ! He4 [C] NIST(http://physics.nist.gov/cgi-bin/cuu/Value?e)
                  !mass_in(nparts)   = 6.64465675E-27 ! He4 (alpha) [kg] NIST (http://physics.nist.gov/cgi-bin/cuu/Value?mal)
                  !Zatom_in(nparts)  = 2.0 !He4 (alpha)
                  IF (vperp_orbit(p) .ne. 0) THEN
                     mu_start_in(nparts) = mass_orbit*vperp_orbit(p)*vperp_orbit(p)/(2*modb)
                  ELSE
                     mu_start_in(nparts) = mu_orbit(p)
                  END IF
                  t_end_in(nparts) = tf
                  nparts = nparts + 1
               END DO
            END DO
         END DO
      END DO
      nparticles_start = nparts - 1
      IF (lscreen) THEN
         WRITE(6, '(/,a,f5.2)') 'BEAMS3D Version ', BEAMS3D_VERSION
         WRITE(6,'(A)') '----- Particle Initialization -----'
         WRITE(6,'(A,F9.5,A,F9.5,A,I6)') '   S   = [',s_min,',',s_max,'];   NS:   ',COUNT(sigma_orbit .lt. bigno)
         WRITE(6,'(A,F9.5,A,F9.5,A,I4)') '   U   = [',0.0,',',pi2*(nu_orbit-1)/nu_orbit,'];   NU:   ',nu_orbit
         WRITE(6,'(A,F9.5,A,F9.5,A,I4)') '   V   = [',0.0,',',pi2*(nv_orbit-1)/nv_orbit/2,'];   NV:   ',nv_orbit
         WRITE(6,'(A,ES10.2,A,ES10.2,A,I4)') '   V_||= [',MINVAL(vll_start_in(1:np_orbit)),',',MAXVAL(vll_start_in(1:np_orbit)),'];  NP:   ',np_orbit
         WRITE(6,'(A,F9.5,A,ES10.2,A,I4)') '   Mu  = [',MINVAL(mu_start_in(1:np_orbit)),',',MAXVAL(mu_start_in(1:np_orbit)),'];  NP:   ',np_orbit
         CALL FLUSH(6)
      END IF

      ! Initialize temperature profiles
      nne = NBEAM_PROF; nte = NBEAM_PROF; nti = NBEAM_PROF
      DO ik = 1, NBEAM_PROF
         s_val = DBLE(ik-1)/DBLE(NBEAM_PROF-1)
         NE_AUX_S_BEAMS(ik) = s_val
         TE_AUX_S_BEAMS(ik) = s_val
         TI_AUX_S_BEAMS(ik) = s_val
         CALL get_equil_ne(s_val,TRIM(ne_type),v_val,iflag)
         NE_AUX_F_BEAMS(ik) = v_val
         CALL get_equil_te(s_val,TRIM(te_type),v_val,iflag)
         TE_AUX_F_BEAMS(ik) = v_val
         CALL get_equil_ti(s_val,TRIM(ti_type),v_val,iflag)
         TI_AUX_F_BEAMS(ik) = v_val
      END DO

      IF (lscreen) THEN
         WRITE(6,'(A)') '----- Profile Initialization -----'
         WRITE(6,'(A,F7.2,A,F7.2,A,I4)') '   Ne  = [',MINVAL(NE_AUX_F_BEAMS)/1E19,',',MAXVAL(NE_AUX_F_BEAMS)/1E19,'] 10^19 [m^-3];  Nne:   ',nne
         WRITE(6,'(A,F7.2,A,F7.2,A,I4)') '   Te  = [',MINVAL(TE_AUX_F_BEAMS)/1000,',',MAXVAL(TE_AUX_F_BEAMS)/1000,'] [keV];  Nte:   ',nte
         WRITE(6,'(A,F7.2,A,F7.2,A,I4)') '   Ti  = [',MINVAL(TI_AUX_F_BEAMS)/1000,',',MAXVAL(TI_AUX_F_BEAMS)/1000,'] [keV];  Nti:   ',nti
         CALL FLUSH(6)
      END IF



      ! Follow the particles
      beam_str = 'beams3d'
      CALL stellopt_paraexe(beam_str,proc_string,lscreen)
      IF (lscreen) WRITE(6,'(a)') '------------------------  ORBIT CALCULATION (DONE)  ----------------------'
      CALL FLUSH(6)

      ! Now analyze the results
      IF (ALLOCATED(orbit_lost_frac)) DEALLOCATE(orbit_lost_frac)
      ALLOCATE(orbit_lost_frac(nsd))
      p = nv_orbit*nu_orbit*np_orbit
      u = 1
      v = p
      orbit_lost_frac = 0
      IF (lscreen) WRITE(6,'(A)') '    ns     flux     Lost(%)'
      DO ik = 1, nsd
         IF (sigma_orbit(ik) .ge. bigno) CYCLE
         orbit_lost_frac(ik) = DBLE(COUNT(lost_lines(u:v)))/DBLE(p)
         IF (lscreen) WRITE(6,'(3X,I3,3X,F9.5,3X,F5.1)') ik,rho(ik),orbit_lost_frac(ik)*100
         u = v + 1
         v = u + p - 1
      END DO
      IF (ALLOCATED(lost_lines)) DEALLOCATE(lost_lines)
      iflag = 0


!DEC$ ENDIF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_orbits
