!-----------------------------------------------------------------------
!     Module:        beams3d_write_ascoth5
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          08/21/2019
!     Description:   This subroutine outputs simulation data for the
!                    ASCOT5 code.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_write_ascoth5(write_type)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
#ifdef LHDF5
      USE hdf5
      USE ez_hdf5
#endif 
      USE beams3d_lines
      USE beams3d_grid, ONLY: nr, nphi, nz, B_R, B_PHI, B_Z, raxis, &
                                 zaxis, phiaxis, S_ARR, U_ARR, POT_ARR, &
                                 ZEFF_ARR, TE, TI, NE, req_axis, zeq_axis, npot, &
                                 POT_SPL_S, ezspline_interp, phiedge_eq, TE_spl_s, &
                                 NE_spl_s, TI_spl_s, ZEFF_spl_s, nne, nte, nti, nzeff, &
                                 plasma_mass, reff_eq, therm_factor, &
                                 X_BEAMLET, Y_BEAMLET, Z_BEAMLET, &
                                 NX_BEAMLET, NY_BEAMLET, NZ_BEAMLET
      USE beams3d_runtime, ONLY: id_string, npoinc, nbeams, beam, t_end, lverb, &
                                    lvmec, lpies, lspec, lcoil, lmgrid, lbeam, lplasma_only, &
                                    lvessel, lvac, lbeam_simple, handle_err, nparticles_start, &
                                    HDF5_OPEN_ERR,HDF5_WRITE_ERR,&
                                    HDF5_CLOSE_ERR, BEAMS3D_VERSION, weight, e_beams, p_beams,&
                                    charge, Zatom, mass, ldepo, v_neut, &
                                    lcollision, pi, pi2, t_end_in, nprocs_beams, &
                                    div_beams, mass_beams, Zatom_beams, dex_beams, &
                                    qid_str_saved, lascotfl
      USE safe_open_mod, ONLY: safe_open
      USE wall_mod, ONLY: nface,nvertex,face,vertex,ihit_array, wall_free, machine_string
      USE beams3d_write_par
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Input Variables
!          write_type  Type of write to preform
!-----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN):: write_type
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          iunit        File ID
!-----------------------------------------------------------------------
      INTEGER :: ier, iunit, i, d1, d2, d3, k, k1, k2, kmax ,ider, nphi1
      INTEGER(HID_T) :: options_gid, bfield_gid, efield_gid, plasma_gid, &
                        neutral_gid, wall_gid, marker_gid, qid_gid, &
                        nbi_gid, inj_gid, boozer_gid, mhd_gid
      INTEGER, ALLOCATABLE, DIMENSION(:) :: itemp
      REAL :: qid_flt
      DOUBLE PRECISION :: rho_temp, s_temp, rho_max, dbl_temp, gammarel, v_total
      DOUBLE PRECISION, ALLOCATABLE :: rtemp(:,:,:), r1dtemp(:)
      CHARACTER(LEN=10) ::  qid_str
      CHARACTER(LEN=8) :: temp_str8, inj_str8


      !DOUBLE PRECISION, PARAMETER :: e_charge      = 1.60217662E-19 !e_c
      DOUBLE PRECISION, PARAMETER :: c_speed       = 2.99792458E+08 !Speed of light [m/s]
      DOUBLE PRECISION, PARAMETER :: e_charge      = 1.602176565e-19 !e_c
      DOUBLE PRECISION, PARAMETER :: inv_amu       = 6.02214076208E+26 ! 1./AMU [1/kg]
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      qid_str='0003271980'
      SELECT CASE (TRIM(write_type))
         CASE('INIT')
            IF (myworkid == master) THEN
               CALL open_hdf5('ascot5_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.true.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'ascot5_'//TRIM(id_string)//'.h5',ier)

               ! Define rho_max for use later on
               rho_max = SQRT(MAXVAL(MAXVAL(MAXVAL(S_ARR,3),2),1))*1.2

               !--------------------------------------------------------------
               !           OPTIONS
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'options', options_gid, ier)
               CALL RANDOM_NUMBER(qid_flt)
               WRITE(qid_str,'(i10.10)') FLOOR(qid_flt*1.0E9)
               qid_str_saved = qid_str
               CALL write_att_hdf5(options_gid,'active',qid_str,ier)
               CALL h5gcreate_f(options_gid,'opt_'//qid_str, qid_gid, ier)
               CALL DATE_AND_TIME(DATE=temp_str8)
               CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
               CALL write_att_hdf5(qid_gid,'description','Data initialized from BEAMS3D',ier)
               i = 2; IF (lascotfl) i =4
               CALL write_var_hdf5(qid_gid,'SIM_MODE',ier,DBLVAR=DBLE(i))
               i = 0; IF (lascotfl) i =1
               CALL write_var_hdf5(qid_gid,'ENABLE_ADAPTIVE',ier,DBLVAR=DBLE(i))
               CALL write_var_hdf5(qid_gid,'RECORD_MODE',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'FIXEDSTEP_USE_USERDEFINED',ier,DBLVAR=DBLE(1))
               CALL write_var_hdf5(qid_gid,'FIXEDSTEP_USERDEFINED',ier,DBLVAR=DBLE(1.0E-8))
               CALL write_var_hdf5(qid_gid,'FIXEDSTEP_GYRODEFINED',ier,DBLVAR=DBLE(16))
               CALL write_var_hdf5(qid_gid,'ADAPTIVE_TOL_ORBIT',ier,DBLVAR=DBLE(1.0E-8))
               CALL write_var_hdf5(qid_gid,'ADAPTIVE_TOL_CCOL',ier,DBLVAR=DBLE(0.01))
               CALL write_var_hdf5(qid_gid,'ADAPTIVE_MAX_DRHO',ier,DBLVAR=DBLE(1.0))
               CALL write_var_hdf5(qid_gid,'ADAPTIVE_MAX_DPHI',ier,DBLVAR=DBLE(2))
               CALL write_var_hdf5(qid_gid,'ENABLE_ORBIT_FOLLOWING',ier,DBLVAR=DBLE(1))
               CALL write_var_hdf5(qid_gid,'ENABLE_MHD',ier,DBLVAR=DBLE(0))
               IF (lcollision .and. .not. lascotfl) THEN
                  CALL write_var_hdf5(qid_gid,'ENABLE_COULOMB_COLLISIONS',ier,DBLVAR=DBLE(1))
                  CALL write_var_hdf5(qid_gid,'DISABLE_ENERGY_CCOLL',ier,DBLVAR=DBLE(0))
                  CALL write_var_hdf5(qid_gid,'DISABLE_PITCH_CCOLL',ier,DBLVAR=DBLE(0))
                  CALL write_var_hdf5(qid_gid,'DISABLE_GCDIFF_CCOLL',ier,DBLVAR=DBLE(1))
               ELSE
                  CALL write_var_hdf5(qid_gid,'ENABLE_COULOMB_COLLISIONS',ier,DBLVAR=DBLE(0))
                  CALL write_var_hdf5(qid_gid,'DISABLE_ENERGY_CCOLL',ier,DBLVAR=DBLE(1))
                  CALL write_var_hdf5(qid_gid,'DISABLE_PITCH_CCOLL',ier,DBLVAR=DBLE(1))
                  CALL write_var_hdf5(qid_gid,'DISABLE_GCDIFF_CCOLL',ier,DBLVAR=DBLE(1))
               END IF
               CALL write_var_hdf5(qid_gid,'DISABLE_FIRSTORDER_GCTRANS',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_SIMTIMELIM',ier,DBLVAR=DBLE(1))
               CALL write_var_hdf5(qid_gid,'ENDCOND_CPUTIMELIM',ier,DBLVAR=DBLE(1))
               CALL write_var_hdf5(qid_gid,'ENDCOND_RHOLIM',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_ENERGYLIM',ier,DBLVAR=DBLE(1))
               CALL write_var_hdf5(qid_gid,'ENDCOND_WALLHIT',ier,DBLVAR=DBLE(1))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MAXORBS',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MAX_SIMTIME',ier,DBLVAR=DBLE(2*MAXVAL(t_end_in)))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MAX_CPUTIME',ier,DBLVAR=DBLE(14000))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MAX_MILEAGE',ier,DBLVAR=DBLE(14000))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MIN_RHO',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MIN_ENERGY',ier,DBLVAR=DBLE(10))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MIN_THERMAL',ier,DBLVAR=DBLE(therm_factor))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MAX_POLOIDALORBS',ier,DBLVAR=DBLE(100))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MAX_TOROIDALORBS',ier,DBLVAR=DBLE(100))
               CALL write_var_hdf5(qid_gid,'ENABLE_DIST_RHO6D',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENABLE_DIST_6D',ier,DBLVAR=DBLE(0))
               i = 1; IF (lascotfl) i =0
               CALL write_var_hdf5(qid_gid,'ENABLE_DIST_RHO5D',ier,DBLVAR=DBLE(i))
               CALL write_var_hdf5(qid_gid,'ENABLE_DIST_5D',ier,DBLVAR=DBLE(i))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_R',ier,DBLVAR=raxis(1))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_R',ier,DBLVAR=raxis(nr))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_Z',ier,DBLVAR=zaxis(1))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_Z',ier,DBLVAR=zaxis(nz))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_PHI',ier,DBLVAR=phiaxis(1)*180/pi)
               CALL write_var_hdf5(qid_gid,'DIST_MAX_PHI',ier,DBLVAR=DBLE(360))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_RHO',ier,DBLVAR=DBLE(0))
               IF (lplasma_only) THEN
                  CALL write_var_hdf5(qid_gid,'ENDCOND_MAX_RHO',ier,DBLVAR=DBLE(1.0))
                  CALL write_var_hdf5(qid_gid,'DIST_MAX_RHO',ier,DBLVAR=DBLE(1.0))
               ELSE
                  CALL write_var_hdf5(qid_gid,'ENDCOND_MAX_RHO',ier,DBLVAR=DBLE(4.00))
                  CALL write_var_hdf5(qid_gid,'DIST_MAX_RHO',ier,DBLVAR=DBLE(1.5))
               ENDIF
               CALL write_var_hdf5(qid_gid,'DIST_MIN_THETA',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_THETA',ier,DBLVAR=DBLE(360))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_TIME',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_TIME',ier,DBLVAR=DBLE(MAXVAL(t_end_in)))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_CHARGE',ier,DBLVAR=DBLE(-2))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_CHARGE',ier,DBLVAR=DBLE(2))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_R',ier,DBLVAR=DBLE(32))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_Z',ier,DBLVAR=DBLE(32))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_RHO',ier,DBLVAR=DBLE(ns_prof1))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_THETA',ier,DBLVAR=DBLE(ns_prof2))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_PHI',ier,DBLVAR=DBLE(ns_prof3))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_TIME',ier,DBLVAR=DBLE(1))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_CHARGE',ier,DBLVAR=DBLE(1))
               i = 0; IF (lascotfl) i =1
               CALL write_var_hdf5(qid_gid,'ENABLE_ORBITWRITE',ier,DBLVAR=DBLE(i))
               i = 1; IF (lascotfl) i =0
               CALL write_var_hdf5(qid_gid,'ORBITWRITE_MODE',ier,DBLVAR=DBLE(i))
               CALL write_var_hdf5(qid_gid,'ORBITWRITE_NPOINT',ier,DBLVAR=DBLE(npoinc))
               CALL write_var_hdf5(qid_gid,'ORBITWRITE_INTERVAL',ier,DBLVAR=DBLE(MAXVAL(t_end_in))/DBLE(npoinc))
               CALL write_var_hdf5(qid_gid,'ORBITWRITE_TOROIDALANGLES',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ORBITWRITE_POLOIDALANGLES',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENABLE_TRANSCOEF',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'TRANSCOEF_INTERVAL',ier,DBLVAR=DBLE(-1))
               CALL write_var_hdf5(qid_gid,'TRANSCOEF_NAVG',ier,DBLVAR=DBLE(5))
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(options_gid, ier)

               !--------------------------------------------------------------
               !           B-FIELD
               !  NOTE On PSI:
               !     In B_STS B_STS_eval_rho defines psi as
               !         rho[0] = sqrt( (psi - Bdata->psi0) / delta );
               !     So psi is TOROIDAL FLUX.
               !--------------------------------------------------------------
               nphi1 = nphi-1 ! confirmed with Poincare plot
               CALL h5gcreate_f(fid,'bfield', bfield_gid, ier)
               CALL RANDOM_NUMBER(qid_flt)
               WRITE(qid_str,'(i10.10)') FLOOR(qid_flt*1.0E9)
               CALL write_att_hdf5(bfield_gid,'active',qid_str,ier)
               CALL h5gcreate_f(bfield_gid,'B_STS_'//qid_str, qid_gid, ier)
               CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
               CALL write_att_hdf5(qid_gid,'description','Data initialized from BEAMS3D',ier)
               CALL write_var_hdf5(qid_gid,'b_nr',ier,INTVAR=nr)
               CALL write_var_hdf5(qid_gid,'b_nphi',ier,INTVAR=nphi1)
               CALL write_var_hdf5(qid_gid,'b_nz',ier,INTVAR=nz)
               CALL write_var_hdf5(qid_gid,'b_rmin',ier,DBLVAR=raxis(1))
               CALL write_var_hdf5(qid_gid,'b_rmax',ier,DBLVAR=raxis(nr))
               CALL write_var_hdf5(qid_gid,'b_zmin',ier,DBLVAR=zaxis(1))
               CALL write_var_hdf5(qid_gid,'b_zmax',ier,DBLVAR=zaxis(nz))
               CALL write_var_hdf5(qid_gid,'b_phimin',ier,DBLVAR=phiaxis(1)*180/pi)
               CALL write_var_hdf5(qid_gid,'b_phimax',ier,DBLVAR=phiaxis(nphi)*180/pi)
               CALL write_var_hdf5(qid_gid,'psi_nr',ier,INTVAR=nr)
               CALL write_var_hdf5(qid_gid,'psi_nphi',ier,INTVAR=nphi1)
               CALL write_var_hdf5(qid_gid,'psi_nz',ier,INTVAR=nz)
               CALL write_var_hdf5(qid_gid,'psi_rmin',ier,DBLVAR=raxis(1))
               CALL write_var_hdf5(qid_gid,'psi_rmax',ier,DBLVAR=raxis(nr))
               CALL write_var_hdf5(qid_gid,'psi_zmin',ier,DBLVAR=zaxis(1))
               CALL write_var_hdf5(qid_gid,'psi_zmax',ier,DBLVAR=zaxis(nz))
               CALL write_var_hdf5(qid_gid,'psi_phimin',ier,DBLVAR=phiaxis(1)*180/pi)
               CALL write_var_hdf5(qid_gid,'psi_phimax',ier,DBLVAR=phiaxis(nphi)*180/pi)
               CALL write_var_hdf5(qid_gid,'axis_nphi',ier,INTVAR=nphi1)
               CALL write_var_hdf5(qid_gid,'axis_phimin',ier,DBLVAR=phiaxis(1)*180/pi)
               CALL write_var_hdf5(qid_gid,'axis_phimax',ier,DBLVAR=phiaxis(nphi)*180/pi)
               CALL write_var_hdf5(qid_gid,'toroidalPeriods',ier,INTVAR=FLOOR(pi2/phiaxis(nphi)))
               CALL write_var_hdf5(qid_gid,'axisr',nphi,ier,DBLVAR=req_axis(1:nphi1))
               CALL write_var_hdf5(qid_gid,'axisz',nphi,ier,DBLVAR=zeq_axis(1:nphi1))
               !CALL write_var_hdf5(qid_gid,'psi0',ier,DBLVAR=DBLE(-1.0E-3)) ! This is because B_STS isn't that 
               CALL write_var_hdf5(qid_gid,'psi0',ier,DBLVAR=DBLE(0)) ! This is because B_STS isn't that 
               CALL write_var_hdf5(qid_gid,'psi1',ier,DBLVAR=DBLE(phiedge_eq))
               ALLOCATE(rtemp(nr,nphi1,nz))
               rtemp = B_R(1:nr,1:nphi1,1:nz)
               CALL write_var_hdf5(qid_gid,'br',nr,nphi1,nz,ier,DBLVAR=rtemp)
               rtemp = B_PHI(1:nr,1:nphi1,1:nz)
               CALL write_var_hdf5(qid_gid,'bphi',nr,nphi1,nz,ier,DBLVAR=rtemp)
               rtemp = B_Z(1:nr,1:nphi1,1:nz)
               CALL write_var_hdf5(qid_gid,'bz',nr,nphi1,nz,ier,DBLVAR=rtemp)
               rtemp = S_ARR(1:nr,1:nphi1,1:nz)
               DO k = 1, nphi1
                  rtemp(:,i,:) = rtemp(:,i,:) - MINVAL(MINVAL(rtemp(:,i,:),DIM=2),DIM=1)
               END DO
               CALL write_var_hdf5(qid_gid,'psi',nr,nphi1,nz,ier,DBLVAR=rtemp*DBLE(phiedge_eq))
               DEALLOCATE(rtemp)
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(bfield_gid, ier)

               !--------------------------------------------------------------
               !           E-FIELD
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'efield', efield_gid, ier)
               CALL RANDOM_NUMBER(qid_flt)
               WRITE(qid_str,'(i10.10)') FLOOR(qid_flt*1.0E9)
               CALL write_att_hdf5(efield_gid,'active',qid_str,ier)
               CALL h5gcreate_f(efield_gid,'E_1DS_'//qid_str, qid_gid, ier)
               CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
               CALL write_att_hdf5(qid_gid,'description','Data initialized from BEAMS3D',ier)
               CALL write_var_hdf5(qid_gid,'rhomin',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'rhomax',ier,DBLVAR=DBLE(rho_max))
               CALL write_var_hdf5(qid_gid,'reff',ier,DBLVAR=DBLE(reff_eq))
               ! Values must be equidistant in rho.
               IF (npot < 1) THEN ! Because we can run with E=0
                  ALLOCATE(r1dtemp(5))
                  r1dtemp = 0
                  CALL write_var_hdf5(qid_gid,'nrho',ier,INTVAR=5)
                  CALL write_var_hdf5(qid_gid,'dvdrho',5,ier,DBLVAR=r1dtemp)
                  DEALLOCATE(r1dtemp)
               ELSE
                  ALLOCATE(r1dtemp(nr))
                  r1dtemp = 0
                  DO i = 1, nr
                     rho_temp = rho_max*DBLE(i-1)/DBLE(nr-1)
                     s_temp = rho_temp*rho_temp
                     ider = 1
                     IF (s_temp<=1) CALL EZspline_derivative1_r8(POT_spl_s,ider,s_temp,r1dtemp(i),ier)
                     ! df/drho = df/ds * ds/drho = df/ds * 2*rho = df/ds * 2 * SQRT(s)
                     r1dtemp(i) = -r1dtemp(i)*2*rho_temp/reff_eq
                  END DO
                  CALL write_var_hdf5(qid_gid,'nrho',ier,INTVAR=nr)
                  CALL write_var_hdf5(qid_gid,'dvdrho',nr,ier,DBLVAR=r1dtemp)
                  DEALLOCATE(r1dtemp)
               END IF
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(efield_gid, ier)


               !--------------------------------------------------------------
               !           PLASMA
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'plasma', plasma_gid, ier)
               CALL RANDOM_NUMBER(qid_flt)
               WRITE(qid_str,'(i10.10)') FLOOR(qid_flt*1.0E9)
               CALL write_att_hdf5(plasma_gid,'active',qid_str,ier)
               ! 1DS acts screwy in the vacuum region
               CALL h5gcreate_f(plasma_gid,'plasma_1DS_'//qid_str, qid_gid, ier)
               !CALL h5gcreate_f(plasma_gid,'plasma_1D_'//qid_str, qid_gid, ier)
               CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
               CALL write_att_hdf5(qid_gid,'description','Data initialized from BEAMS3D',ier)
               CALL write_var_hdf5(qid_gid,'nion',ier,INTVAR=1)
               CALL write_var_hdf5(qid_gid,'nrho',ier,INTVAR=nr)
               CALL write_var_hdf5(qid_gid,'znum',ier,INTVAR=1)
               CALL write_var_hdf5(qid_gid,'anum',ier,INTVAR=NINT(plasma_mass*inv_amu))
               CALL write_var_hdf5(qid_gid,'charge',ier,INTVAR=1)
               CALL write_var_hdf5(qid_gid,'mass',ier,DBLVAR=plasma_mass*inv_amu)
               ALLOCATE(rtemp(nr,5,1))
               ! Only for 1DS
               CALL write_var_hdf5(qid_gid,'rhomin',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'rhomax',ier,DBLVAR=DBLE(rho_max))
               rtemp = 0
               rtemp(:,5,1) = 1
               DO i = 1, nr
                  rtemp(i,1,1)=rho_max*DBLE(i-1)/DBLE(nr-1)
               END DO
               d1 = COUNT(rtemp(:,1,1) <= 1)+1
               IF (nte > 0)   CALL EZspline_interp( TE_spl_s,   nr, rtemp(:,1,1)**2, rtemp(:,2,1), ier)
               IF (nne > 0)   CALL EZspline_interp( NE_spl_s,   nr, rtemp(:,1,1)**2, rtemp(:,3,1), ier)
               IF (nti > 0)   CALL EZspline_interp( TI_spl_s,   nr, rtemp(:,1,1)**2, rtemp(:,4,1), ier)
               IF (nzeff > 0) CALL EZspline_interp( ZEFF_spl_s, nr, rtemp(:,1,1)**2, rtemp(:,5,1), ier)
               rtemp(d1:,2:4,1) = 0
               rtemp(:,5,1)=rtemp(:,3,1)/rtemp(:,5,1)
               WHERE(rtemp(:,2,1) < 10) rtemp(:,2,1)=10
               WHERE(rtemp(:,4,1) < 10) rtemp(:,4,1)=10
               WHERE(rtemp(:,3,1) < 1.0E18) rtemp(:,3,1)=1.0E18
               WHERE(rtemp(:,5,1) < 1.0E18) rtemp(:,5,1)=1.0E18
               CALL write_var_hdf5( qid_gid, 'rho',          nr, ier, DBLVAR=rtemp(1:nr,1,1))
               CALL write_var_hdf5( qid_gid, 'etemperature', nr, ier, DBLVAR=rtemp(1:nr,2,1))
               CALL write_var_hdf5( qid_gid, 'edensity',     nr, ier, DBLVAR=rtemp(1:nr,3,1))
               CALL write_var_hdf5( qid_gid, 'itemperature', nr, ier, DBLVAR=rtemp(1:nr,4,1))
               CALL write_var_hdf5( qid_gid, 'idensity',     nr, ier, DBLVAR=rtemp(1:nr,5,1))
               ! Need to add idensity
               DEALLOCATE(rtemp)
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(plasma_gid, ier)


               !--------------------------------------------------------------
               !           NEUTRAL
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'neutral', neutral_gid, ier)
               CALL RANDOM_NUMBER(qid_flt)
               WRITE(qid_str,'(i10.10)') FLOOR(qid_flt*1.0E9)
               CALL write_att_hdf5(neutral_gid,'active',qid_str,ier)
               CALL h5gcreate_f(neutral_gid,'N0_3D_'//qid_str, qid_gid, ier)
               CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
               CALL write_att_hdf5(qid_gid,'description','Data initialized from BEAMS3D',ier)
               CALL write_var_hdf5(qid_gid,'nr',ier,INTVAR=nr)
               CALL write_var_hdf5(qid_gid,'nphi',ier,INTVAR=nphi1)
               CALL write_var_hdf5(qid_gid,'nz',ier,INTVAR=nz)
               CALL write_var_hdf5(qid_gid,'rmin',ier,DBLVAR=raxis(1))
               CALL write_var_hdf5(qid_gid,'rmax',ier,DBLVAR=raxis(nr))
               CALL write_var_hdf5(qid_gid,'zmin',ier,DBLVAR=zaxis(1))
               CALL write_var_hdf5(qid_gid,'zmax',ier,DBLVAR=zaxis(nz))
               CALL write_var_hdf5(qid_gid,'phimin',ier,DBLVAR=phiaxis(1)*180/pi)
               CALL write_var_hdf5(qid_gid,'phimax',ier,DBLVAR=phiaxis(nphi)*180/pi)
               CALL write_var_hdf5(qid_gid,'nspecies',ier,INTVAR=1)
               CALL write_var_hdf5(qid_gid,'anum',ier,INTVAR=1)
               CALL write_var_hdf5(qid_gid,'znum',ier,INTVAR=1)
               CALL write_var_hdf5(qid_gid,'maxwellian',ier,INTVAR=1)
               ALLOCATE(rtemp(nr,nphi1,nz))
               rtemp = 0
               CALL write_var_hdf5(qid_gid,'density',nr,nphi1,nz,ier,DBLVAR=rtemp)
               CALL write_var_hdf5(qid_gid,'temperature',nr,nphi1,nz,ier,DBLVAR=rtemp)
               DEALLOCATE(rtemp)
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(neutral_gid, ier)


               !--------------------------------------------------------------
               !           BOOZER (DUMMY)
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'boozer', boozer_gid, ier)
               CALL RANDOM_NUMBER(qid_flt)
               WRITE(qid_str,'(i10.10)') FLOOR(qid_flt*1.0E9)
               CALL write_att_hdf5(boozer_gid,'active',qid_str,ier)
               CALL h5gcreate_f(boozer_gid,'Boozer_'//qid_str, qid_gid, ier)
               CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
               CALL write_att_hdf5(qid_gid,'description','Dummy data initialized from BEAMS3D',ier)
               CALL write_var_hdf5(qid_gid,'nr',ier,INTVAR=nr)
               CALL write_var_hdf5(qid_gid,'nz',ier,INTVAR=nz)
               CALL write_var_hdf5(qid_gid,'rmin',ier,DBLVAR=raxis(1))
               CALL write_var_hdf5(qid_gid,'rmax',ier,DBLVAR=raxis(nr))
               CALL write_var_hdf5(qid_gid,'zmin',ier,DBLVAR=zaxis(1))
               CALL write_var_hdf5(qid_gid,'zmax',ier,DBLVAR=zaxis(nz))
               CALL write_var_hdf5(qid_gid,'npsi',ier,INTVAR=ns_prof1)
               CALL write_var_hdf5(qid_gid,'psimin',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'psimax',ier,DBLVAR=ABS(phiedge_eq))
               CALL write_var_hdf5(qid_gid,'psi0',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'psi1',ier,DBLVAR=ABS(phiedge_eq))
               CALL write_var_hdf5(qid_gid,'ntheta',ier,INTVAR=ns_prof2)
               CALL write_var_hdf5(qid_gid,'nthetag',ier,INTVAR=ns_prof2)
               CALL write_var_hdf5(qid_gid,'r0',ier,DBLVAR=req_axis(1))
               CALL write_var_hdf5(qid_gid,'z0',ier,DBLVAR=zeq_axis(1))
               CALL write_var_hdf5(qid_gid,'nrzs',ier,INTVAR=ns_prof1)
               ALLOCATE(rtemp(ns_prof1,1,1))
               rtemp = 0;
               CALL write_var_hdf5(qid_gid,'rs',ns_prof1,ier,DBLVAR=rtemp(:,1,1))
               CALL write_var_hdf5(qid_gid,'zs',ns_prof1,ier,DBLVAR=rtemp(:,1,1))
               DEALLOCATE(rtemp)
               ALLOCATE(rtemp(nr,nz,1))
               rtemp =0
               CALL write_var_hdf5(qid_gid,'psi_rz',nr,nz,ier,DBLVAR=rtemp(:,:,1))
               DEALLOCATE(rtemp)
               ALLOCATE(rtemp(ns_prof1,ns_prof2,1))
               rtemp =0
               CALL write_var_hdf5(qid_gid,'theta_psithetageom',ns_prof1,ns_prof2,ier,DBLVAR=rtemp(:,:,1))
               CALL write_var_hdf5(qid_gid,'nu_psitheta',ns_prof1,ns_prof2,ier,DBLVAR=rtemp(:,:,1))
               DEALLOCATE(rtemp)
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(boozer_gid, ier)


               !--------------------------------------------------------------
               !           MHD (DUMMY)
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'mhd', mhd_gid, ier)
               CALL RANDOM_NUMBER(qid_flt)
               WRITE(qid_str,'(i10.10)') FLOOR(qid_flt*1.0E9)
               CALL write_att_hdf5(mhd_gid,'active',qid_str,ier)
               CALL h5gcreate_f(mhd_gid,'MHD_STAT_'//qid_str, qid_gid, ier)
               CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
               CALL write_att_hdf5(qid_gid,'description','Dummy data initialized from BEAMS3D',ier)
               CALL write_var_hdf5(qid_gid,'nmode',ier,INTVAR=2)
               CALL write_var_hdf5(qid_gid,'nrho',ier,INTVAR=6)
               CALL write_var_hdf5(qid_gid,'rhomin',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'rhomax',ier,DBLVAR=DBLE(1))
               CALL write_var_hdf5(qid_gid,'nmodes',2,ier,INTVAR=(/1,2/))
               CALL write_var_hdf5(qid_gid,'mmodes',2,ier,INTVAR=(/3,4/))
               CALL write_var_hdf5(qid_gid,'amplitude',2,ier,DBLVAR=DBLE((/0.0,0.0/)))
               CALL write_var_hdf5(qid_gid,'omega',2,ier,DBLVAR=DBLE((/1.0,1.5/)))
               CALL write_var_hdf5(qid_gid,'phase',2,ier,DBLVAR=DBLE((/0.0,0.78525/)))
               ALLOCATE(rtemp(6,2,1))
               rtemp = 1;
               CALL write_var_hdf5(qid_gid,'alpha',6,2,ier,DBLVAR=rtemp(:,:,1))
               CALL write_var_hdf5(qid_gid,'phi',6,2,ier,DBLVAR=rtemp(:,:,1))
               DEALLOCATE(rtemp)
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(mhd_gid, ier)


               !--------------------------------------------------------------
               !           WALL
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'wall', wall_gid, ier)
               CALL RANDOM_NUMBER(qid_flt)
               WRITE(qid_str,'(i10.10)') FLOOR(qid_flt*1.0E9)
               CALL write_att_hdf5(wall_gid,'active',qid_str,ier)
               CALL h5gcreate_f(wall_gid,'wall_3D_'//qid_str, qid_gid, ier)
               CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
               CALL write_att_hdf5(qid_gid,'description','BEAMS3D: '//TRIM(machine_string),ier)
               CALL write_var_hdf5(qid_gid,'nelements',ier,INTVAR=nface)
               ALLOCATE(rtemp(3,nface,3))
               DO i = 1, nface
                  d1 = face(i,1)
                  d2 = face(i,2)
                  d3 = face(i,3)
                  ! X
                  rtemp(1,i,1) = vertex(d1,1)
                  rtemp(2,i,1) = vertex(d2,1)
                  rtemp(3,i,1) = vertex(d3,1)
                  ! Y
                  rtemp(1,i,2) = vertex(d1,2)
                  rtemp(2,i,2) = vertex(d2,2)
                  rtemp(3,i,2) = vertex(d3,2)
                  ! Z
                  rtemp(1,i,3) = vertex(d1,3)
                  rtemp(2,i,3) = vertex(d2,3)
                  rtemp(3,i,3) = vertex(d3,3)
               END DO
               CALL write_var_hdf5( qid_gid, 'x1x2x3', 3, nface, ier, DBLVAR=rtemp(:,:,1))
               CALL write_var_hdf5( qid_gid, 'y1y2y3', 3, nface, ier, DBLVAR=rtemp(:,:,2))
               CALL write_var_hdf5( qid_gid, 'z1z2z3', 3, nface, ier, DBLVAR=rtemp(:,:,3))
               DEALLOCATE(rtemp)
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(wall_gid, ier)

               ! Close file
               CALL close_hdf5(fid,ier)
               IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'ascot5_'//TRIM(id_string)//'.h5',ier)
            END IF
         CASE('BBNBI')
            !--------------------------------------------------------------
            !           BEAM Model
            !           Note BEAMS3D assumes each beam represents an energy
            !--------------------------------------------------------------
            IF (myworkid == master) THEN
               CALL open_hdf5('ascot5_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.false.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'ascot5_'//TRIM(id_string)//'.h5',ier)
               CALL h5gcreate_f(fid,'nbi', nbi_gid, ier)
               CALL RANDOM_NUMBER(qid_flt)
               WRITE(qid_str,'(i10.10)') FLOOR(qid_flt*1.0E9)
               CALL write_att_hdf5(nbi_gid,'active',qid_str,ier)
               CALL h5gcreate_f(nbi_gid,'nbi_'//qid_str, qid_gid, ier)
               CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
               CALL write_att_hdf5(qid_gid,'description','Data initialized from BEAMS3D',ier)
               CALL write_var_hdf5(qid_gid,'ninj',ier,INTVAR=nbeams)
               DO i = 1, nbeams
                  d1 = Dex_beams(i)
                  WRITE(inj_str8,'(i8)') i
                  CALL h5gcreate_f(qid_gid,'inj'//ADJUSTL(inj_str8), inj_gid, ier)
                  CALL write_var_hdf5(inj_gid,'id',ier,INTVAR=i)
                  CALL write_var_hdf5(inj_gid,'power',ier,DBLVAR=P_BEAMS(i))
                  CALL write_var_hdf5(inj_gid,'energy',ier,DBLVAR=E_BEAMS(i))
                  CALL write_var_hdf5(inj_gid,'efrac',3,ier,DBLVAR=DBLE((/1, 0, 0/)))
                  CALL write_var_hdf5(inj_gid,'div_h',ier,DBLVAR=div_beams(i))
                  CALL write_var_hdf5(inj_gid,'div_v',ier,DBLVAR=div_beams(i))
                  dbl_temp = 0
                  CALL write_var_hdf5(inj_gid,'div_halo_frac',ier,DBLVAR=dbl_temp)
                  CALL write_var_hdf5(inj_gid,'div_halo_v',ier,DBLVAR=dbl_temp)
                  CALL write_var_hdf5(inj_gid,'div_halo_h',ier,DBLVAR=dbl_temp)
                  CALL write_var_hdf5(inj_gid,'anum',ier,INTVAR=NINT(mass_beams(i)*inv_amu))
                  CALL write_var_hdf5(inj_gid,'znum',ier,INTVAR=NINT(Zatom_beams(i)))
                  CALL write_var_hdf5(inj_gid,'mass',ier,DBLVAR=mass_beams(i))
                  ! Now Beamlets
                  d2 = SIZE(X_BEAMLET(d1,:))
                  ALLOCATE(r1dtemp(d2))

                  CALL write_var_hdf5(inj_gid,'nbeamlet',ier,INTVAR=d2)
                  r1dtemp = X_BEAMLET(d1,:)
                  CALL write_var_hdf5(inj_gid,'beamletx',d2,ier,DBLVAR=r1dtemp)
                  r1dtemp = Y_BEAMLET(d1,:)
                  CALL write_var_hdf5(inj_gid,'beamlety',d2,ier,DBLVAR=r1dtemp)
                  r1dtemp = Z_BEAMLET(d1,:)
                  CALL write_var_hdf5(inj_gid,'beamletz',d2,ier,DBLVAR=r1dtemp)
                  r1dtemp = NX_BEAMLET(d1,:)
                  CALL write_var_hdf5(inj_gid,'beamletdx',d2,ier,DBLVAR=r1dtemp)
                  r1dtemp = NY_BEAMLET(d1,:)
                  CALL write_var_hdf5(inj_gid,'beamletdy',d2,ier,DBLVAR=r1dtemp)
                  r1dtemp = NZ_BEAMLET(d1,:)
                  CALL write_var_hdf5(inj_gid,'beamletdz',d2,ier,DBLVAR=r1dtemp)
                  
                  DEALLOCATE(r1dtemp)

                  CALL h5gclose_f(inj_gid, ier)
               END DO
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(nbi_gid, ier)

               ! Close File
               CALL close_hdf5(fid,ier)
               IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'ascot5_'//TRIM(id_string)//'.h5',ier)
            END IF

         CASE('MARKER')
            ! For now we assume we will only run this in deposition mode
            ! so that means we look at index (1,i) since
            !    0: starting point
            !    1: Ionization point (or wall)
            !    2: Gyrocenter
            d1 = LBOUND(R_lines,DIM=2)
            d2 = UBOUND(R_lines,DIM=2)
            d3 = 0
            ALLOCATE(itemp(nprocs_beams))
            itemp = 0
            IF (lbeam) d3 = 2
            itemp(myworkid+1) = COUNT(end_state(d1:d2).lt.3) ! Don't count shinethrough particles or port particles
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, itemp, nprocs_beams, MPI_INTEGER, MPI_SUM, MPI_COMM_BEAMS, ierr_mpi)
            k2 = SUM(itemp(1:myworkid+1))
            k1 = k2 - itemp(myworkid+1) + 1
            kmax = SUM(itemp)
            DEALLOCATE(itemp)
            IF (myworkid == master) THEN
               CALL open_hdf5('ascot5_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.false.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'ascot5_'//TRIM(id_string)//'.h5',ier)
               !--------------------------------------------------------------
               !           MARKER
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'marker', marker_gid, ier)
               CALL RANDOM_NUMBER(qid_flt)
               WRITE(qid_str,'(i10.10)') FLOOR(qid_flt*1.0E9)
               CALL write_att_hdf5(marker_gid,'active',qid_str,ier)
               CALL h5gcreate_f(marker_gid,'gc_'//qid_str, qid_gid, ier)
               CALL write_var_hdf5(qid_gid,'n',ier,INTVAR=kmax)
               CALL DATE_AND_TIME(DATE=temp_str8)
               CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
               CALL write_att_hdf5(qid_gid,'description','Data initialized from BEAMS3D',ier)
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(marker_gid, ier)
               !--------------------------------------------------------------
               !           Update options
               !--------------------------------------------------------------
               CALL h5gopen_f(fid,'options', options_gid, ier)
               CALL h5gopen_f(options_gid,'opt_'//qid_str_saved, qid_gid, ier)
               CALL write_var_hdf5(qid_gid,'DIST_MIN_PR',ier,DBLVAR=DBLE(-partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_PR',ier,DBLVAR=DBLE(partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_PPHI',ier,DBLVAR=DBLE(-partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_PPHI',ier,DBLVAR=DBLE(partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_PZ',ier,DBLVAR=DBLE(-partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_PZ',ier,DBLVAR=DBLE(partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_PPA',ier,DBLVAR=DBLE(-partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_PPA',ier,DBLVAR=DBLE(partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_PPE',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_PPE',ier,DBLVAR=DBLE(partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_PR',ier,DBLVAR=DBLE(ns_prof4))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_PPHI',ier,DBLVAR=DBLE(ns_prof4))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_PZ',ier,DBLVAR=DBLE(ns_prof4))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_PPA',ier,DBLVAR=DBLE(ns_prof4))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_PPE',ier,DBLVAR=DBLE(ns_prof5))
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(options_gid, ier)
               ! Close file
               CALL close_hdf5(fid,ier)
               IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'ascot5_'//TRIM(id_string)//'.h5',ier)
            END IF
            CALL MPI_BCAST(qid_str,10,MPI_CHARACTER,master,MPI_COMM_BEAMS,ierr_mpi)
            ALLOCATE(rtemp(k1:k2,13,1))
            k = k1
            DO i = d1, d2
               IF (end_state(i) .ge. 3) CYCLE
               rtemp(k,1,1) = R_lines(d3,i)
               rtemp(k,2,1) = PHI_lines(d3,i)*180/pi
               rtemp(k,3,1) = Z_lines(d3,i)
               dbl_temp     = 2*B_lines(d3,i)*moment_lines(d3,i)/mass(i) ! V_perp^2
               v_total      = SQRT(dbl_temp+vll_lines(d3,i)*vll_lines(d3,i))
               gammarel     = SQRT(1.0 / ( (1.0+v_total/c_speed) * (1.0-v_total/c_speed) ))
               rtemp(k,4,1) = (gammarel-1.0)*mass(i)*c_speed*c_speed/e_charge
               rtemp(k,5,1) = vll_lines(d3,i)/v_total
               !rtemp(k,4,1) = 0.5*mass(i)*(vll_lines(d3,i)*vll_lines(d3,i)+dbl_temp)/e_charge
               !rtemp(i,4,1) = 0.5*mass(i)*(vll_lines(d3,i)*vll_lines(d3,i)+dbl_temp)
               !rtemp(k,5,1) = vll_lines(d3,i)/SQRT(dbl_temp+vll_lines(d3,i)*vll_lines(d3,i)) ! pitch
               rtemp(k,6,1) = 0          ! zeta
               rtemp(k,7,1) = mass(i)*inv_amu ! mass
               rtemp(k,8,1) = Zatom(i)
               rtemp(k,9,1) = NINT(mass(i)*inv_amu) ! Anum
               rtemp(k,10,1) = Zatom(i)
               rtemp(k,11,1) = weight(i) ! weight
               rtemp(k,12,1) = 0.0 ! time
               rtemp(k,13,1) = k
               k=k+1
            END DO
            WHERE(rtemp(:,4,1)==0) rtemp(:,5,1)=0.0 ! avoid NAN pitch angle
            WHERE(rtemp(:,1,1)<raxis(1)) rtemp(:,1,1)=raxis(1) ! avoid NAN pitch angle
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/gc_'//qid_str//'/r',      DBLVAR=rtemp(:,1,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/gc_'//qid_str//'/phi',    DBLVAR=rtemp(:,2,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/gc_'//qid_str//'/z',      DBLVAR=rtemp(:,3,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/gc_'//qid_str//'/energy', DBLVAR=rtemp(:,4,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/gc_'//qid_str//'/pitch',  DBLVAR=rtemp(:,5,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/gc_'//qid_str//'/zeta',   DBLVAR=rtemp(:,6,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/gc_'//qid_str//'/mass',   DBLVAR=rtemp(:,7,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/gc_'//qid_str//'/charge', INTVAR=INT(rtemp(:,8,1)))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/gc_'//qid_str//'/anum',   INTVAR=INT(rtemp(:,9,1)))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/gc_'//qid_str//'/znum',   INTVAR=INT(rtemp(:,10,1)))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/gc_'//qid_str//'/weight', DBLVAR=rtemp(:,11,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/gc_'//qid_str//'/time',   DBLVAR=rtemp(:,12,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/gc_'//qid_str//'/id',     INTVAR=INT(rtemp(:,13,1)))
            DEALLOCATE(rtemp)

         CASE('FIELDLINES')
            ! For now we assume we will only run this in deposition mode
            ! so that means we look at index (1,i) since
            !    0: starting point
            !    1: Ionization point (or wall)
            !    2: Gyrocenter
            d1 = LBOUND(R_lines,DIM=2)
            d2 = UBOUND(R_lines,DIM=2)
            d3 = 0
            ALLOCATE(itemp(nprocs_beams))
            itemp = 0
            IF (lbeam) d3 = 2
            itemp(myworkid+1) = COUNT(end_state(d1:d2).lt.3) ! Don't count shinethrough particles
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, itemp, nprocs_beams, MPI_INTEGER, MPI_SUM, MPI_COMM_BEAMS, ierr_mpi)
            k2 = SUM(itemp(1:myworkid+1))
            k1 = k2 - itemp(myworkid+1) + 1
            kmax = SUM(itemp)
            DEALLOCATE(itemp)
            IF (myworkid == master) THEN
               CALL open_hdf5('ascot5_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.false.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'ascot5_'//TRIM(id_string)//'.h5',ier)
               !--------------------------------------------------------------
               !           MARKER
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'marker', marker_gid, ier)
               CALL RANDOM_NUMBER(qid_flt)
               WRITE(qid_str,'(i10.10)') FLOOR(qid_flt*1.0E9)
               CALL write_att_hdf5(marker_gid,'active',qid_str,ier)
               CALL h5gcreate_f(marker_gid,'fl_'//qid_str, qid_gid, ier)
               CALL write_var_hdf5(qid_gid,'n',ier,INTVAR=kmax)
               CALL DATE_AND_TIME(DATE=temp_str8)
               CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
               CALL write_att_hdf5(qid_gid,'description','Data initialized from BEAMS3D',ier)
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(marker_gid, ier)
               !--------------------------------------------------------------
               !           Update options
               !--------------------------------------------------------------
               CALL h5gopen_f(fid,'options', options_gid, ier)
               CALL h5gopen_f(options_gid,'opt_'//qid_str_saved, qid_gid, ier)
               CALL write_var_hdf5(qid_gid,'DIST_MIN_PR',ier,DBLVAR=DBLE(-partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_PR',ier,DBLVAR=DBLE(partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_PPHI',ier,DBLVAR=DBLE(-partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_PPHI',ier,DBLVAR=DBLE(partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_PZ',ier,DBLVAR=DBLE(-partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_PZ',ier,DBLVAR=DBLE(partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_PPA',ier,DBLVAR=DBLE(-partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_PPA',ier,DBLVAR=DBLE(partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_MIN_PPE',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'DIST_MAX_PPE',ier,DBLVAR=DBLE(partpmax))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_PR',ier,DBLVAR=DBLE(ns_prof4))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_PPHI',ier,DBLVAR=DBLE(ns_prof4))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_PZ',ier,DBLVAR=DBLE(ns_prof4))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_PPA',ier,DBLVAR=DBLE(ns_prof4))
               CALL write_var_hdf5(qid_gid,'DIST_NBIN_PPE',ier,DBLVAR=DBLE(ns_prof5))
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(options_gid, ier)
               ! Close file
               CALL close_hdf5(fid,ier)
               IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'ascot5_'//TRIM(id_string)//'.h5',ier)
            END IF
            CALL MPI_BCAST(qid_str,10,MPI_CHARACTER,master,MPI_COMM_BEAMS,ierr_mpi)
            ALLOCATE(rtemp(k1:k2,7,1))
            k = k1
            DO i = d1, d2
               IF (end_state(i) .ge. 3) CYCLE
               rtemp(k,1,1) = R_lines(d3,i)
               rtemp(k,2,1) = PHI_lines(d3,i)*180/pi
               rtemp(k,3,1) = Z_lines(d3,i)
               dbl_temp     = 2*B_lines(d3,i)*moment_lines(d3,i)/mass(i) ! V_perp^2
               rtemp(k,4,1) = vll_lines(d3,i)/SQRT(dbl_temp+vll_lines(d3,i)*vll_lines(d3,i)) ! pitch
               rtemp(k,5,1) = weight(i) ! weight
               rtemp(k,6,1) = 0.0 ! time
               rtemp(k,7,1) = k
               k=k+1
            END DO
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/fl_'//qid_str//'/r',      DBLVAR=rtemp(:,1,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/fl_'//qid_str//'/phi',    DBLVAR=rtemp(:,2,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/fl_'//qid_str//'/z',      DBLVAR=rtemp(:,3,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/fl_'//qid_str//'/pitch',  DBLVAR=rtemp(:,4,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/fl_'//qid_str//'/weight', DBLVAR=rtemp(:,5,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/fl_'//qid_str//'/time',   DBLVAR=rtemp(:,6,1))
            CALL beams3d_write1d_parhdf5( 1, kmax, k1, k2, '/marker/fl_'//qid_str//'/id',     INTVAR=INT(rtemp(:,7,1)))
            DEALLOCATE(rtemp)
      END SELECT

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_write_ascoth5

