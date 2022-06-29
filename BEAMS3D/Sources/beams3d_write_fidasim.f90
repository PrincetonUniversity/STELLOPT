!-----------------------------------------------------------------------
!     Module:        beams3d_write_fidasim
!     Authors:       D. Kulla (david.kulla@ipp.mpg.de)
!     Date:          04/01/2021
!     Description:   This subroutine outputs simulation data for the
!                    FIDASIM code.
!-----------------------------------------------------------------------
SUBROUTINE beams3d_write_fidasim(write_type)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
   USE stel_kinds, ONLY: rprec
#ifdef LHDF5
      USE hdf5
      USE ez_hdf5
#endif
   USE beams3d_lines, ONLY: ns_prof1, ns_prof2, ns_prof3, &
      ns_prof4, ns_prof5, dist5d_prof, &
      partvmax, dist5D_fida, &
      h2_prof, h3_prof, h4_prof, h5_prof, &
      nsh_prof4,  my_end, r_h, p_h, z_h
   USE beams3d_grid, ONLY: nr, nphi, nz, B_R, B_PHI, B_Z, raxis, &
      zaxis, phiaxis, S_ARR, U_ARR, POT_ARR, &
      ZEFF_ARR, TE, TI, NE, npot, nte, nti, nzeff, &
      X_BEAMLET, Y_BEAMLET, Z_BEAMLET, &
      NX_BEAMLET, NY_BEAMLET, NZ_BEAMLET, &
      POT4D, NE4D, TE4D, TI4D, ZEFF4D, &
      BR4D, BPHI4D, BZ4D, &
      hr, hp, hz, hri, hpi, hzi, S4D, U4D, X4D, Y4D, &
      rmin, rmax, zmin, zmax, phimin, phimax, raxis, zaxis, phiaxis, &
      rmin_fida, rmax_fida, zmin_fida, zmax_fida, phimin_fida, phimax_fida, &
      raxis_fida, zaxis_fida, phiaxis_fida, nr_fida, nphi_fida, nz_fida, &
      nenergy_fida, npitch_fida, energy_fida, pitch_fida, t_fida
   USE beams3d_runtime, ONLY: id_string, npoinc, nbeams, beam, t_end, lverb, &
      lvmec, lpies, lspec, lcoil, lmgrid, lbeam, lplasma_only, &
      lvessel, lvac, lbeam_simple, handle_err, nparticles_start, &
      HDF5_OPEN_ERR,HDF5_WRITE_ERR,&
      HDF5_CLOSE_ERR, BEAMS3D_VERSION, weight, e_beams, p_beams,&
      charge, Zatom, mass, ldepo, v_neut, &
      lcollision, pi, pi2, t_end_in, nprocs_beams, &
      div_beams, mass_beams, Zatom_beams, dex_beams, &
      lascotfl, R_beams, PHI_beams, Z_beams, &
      lsplit, lfidasim2
   USE safe_open_mod, ONLY: safe_open
   !USE wall_mod, ONLY: nface,nvertex,face,vertex,ihit_array, wall_free, machine_string
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
   INTEGER :: ier, iunit, i, j, d1, d2, d3, k, k1, k2, kmax ,ider, nphi1, &
      l, m, n, b, i3, j3, k3
   INTEGER(HID_T) :: options_gid, bfield_gid, efield_gid, plasma_gid, &
      neutral_gid, wall_gid, marker_gid, qid_gid, &
      nbi_gid, inj_gid, boozer_gid, mhd_gid, temp_gid
   INTEGER, ALLOCATABLE, DIMENSION(:) :: itemp
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mask

   REAL*8 :: fvalE(1,3), fval(1), fval2(1), xparam, yparam, zparam
   REAL(rprec) :: jac, v_parr, v_perp, pitch, v, phi_temp, r_temp, z_temp
   REAL(rprec), DIMENSION(4) :: rt,zt,pt
   REAL(rprec), DIMENSION(:,:,:,:,:,:), POINTER :: dist5d_temp
   !REAL*8 :: xparam, yparam, zparam!, hx, hy, hz, hxi, hyi, hzi
   !REAL*8 :: fval(1), fval2(1)

   DOUBLE PRECISION         :: x0, y0, z0, rho_temp, s_temp, dbl_temp, gammarel, v_total, vol
   DOUBLE PRECISION, ALLOCATABLE :: rtemp(:,:,:), rtemp2(:,:,:), rtemp3(:,:,:), rtemp4(:,:,:), r1dtemp(:), r1dtemp2(:), r2dtemp(:,:), r4dtemp(:,:,:,:)

   CHARACTER(LEN=8) :: temp_str8, inj_str8

   INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/), ictE(8)=(/0,1,1,1,0,0,0,0/)
   REAL*8, PARAMETER :: one = 1
   DOUBLE PRECISION, PARAMETER :: c_speed       = 2.99792458E+08 !Speed of light [m/s]
   DOUBLE PRECISION, PARAMETER :: e_charge      = 1.602176565e-19 !e_c
   DOUBLE PRECISION, PARAMETER :: inv_amu       = 6.02214076208E+26 ! 1./AMU [1/kg]
   DOUBLE PRECISION, PARAMETER :: zero          = 0.0D0 ! 0.0
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
   IF (myworkid == master) THEN
      SELECT CASE (TRIM(write_type))
       CASE('INIT')
!           WRITE(327,*) rmin_fida, rmin, zmin_fida, zmin, nr_fida, nz_fida

         nphi1 = nphi-1 ! confirmed with Poincare plot
         CALL open_hdf5('fidasim_'//TRIM(id_string)//'_distribution.h5',fid,ier,LCREATE=.true.)
         CALL h5gopen_f(fid,'/', qid_gid, ier)
         CALL write_att_hdf5(qid_gid,'data_source','Data initialized from BEAMS3D',ier)
         CALL write_var_hdf5(qid_gid,'type',ier,INTVAR=1)
         CALL h5dopen_f(fid, 'type', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'description','Distribution type: 1="Guiding Center Density Function", 2="Guiding Center Monte Carlo", 3="Full Orbit Monte Carlo"',ier)
         CALL h5dclose_f(temp_gid,ier)


         ! FIDASIM GRID
         CALL write_var_hdf5(qid_gid,'nenergy',ier,INTVAR=nenergy_fida)
         CALL h5dopen_f(fid, 'nenergy', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','-',ier)
         CALL write_att_hdf5(temp_gid,'description','Number of energy values',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'npitch',ier,INTVAR=npitch_fida)
         CALL h5dopen_f(fid, 'npitch', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','-',ier)
         CALL write_att_hdf5(temp_gid,'description','Number of pitch values',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'nr',ier,INTVAR=nr_fida)
         CALL h5dopen_f(fid, 'nr', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','-',ier)
         CALL write_att_hdf5(temp_gid,'description','Number of R values',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'nphi',ier,INTVAR=nphi_fida)
         CALL h5dopen_f(fid, 'nphi', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','-',ier)
         CALL write_att_hdf5(temp_gid,'description','Number of Phi values',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'nz',ier,INTVAR=nz_fida)
         CALL h5dopen_f(fid, 'nz', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','-',ier)
         CALL write_att_hdf5(temp_gid,'description','Number of Z values',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'time',ier,DBLVAR=DBLE(t_fida))
         CALL h5dopen_f(fid, 'time', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','s',ier)
         CALL write_att_hdf5(temp_gid,'description','Distribution time',ier)
         CALL h5dclose_f(temp_gid,ier)

         CALL write_var_hdf5(qid_gid,'r',nr_fida, ier,DBLVAR=DBLE(raxis_fida*100)) !convert from m to cm
         CALL h5dopen_f(fid, 'r', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','cm',ier)
         CALL write_att_hdf5(temp_gid,'description','Radius',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'phi',nphi_fida, ier,DBLVAR=DBLE(phiaxis_fida))
         CALL h5dopen_f(fid, 'phi', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','rad',ier)
         CALL write_att_hdf5(temp_gid,'description','Toroidal angle',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'z',nz_fida, ier,DBLVAR=DBLE(zaxis_fida*100)) !convert from m to cm
         CALL h5dopen_f(fid, 'z', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','cm',ier)
         CALL write_att_hdf5(temp_gid,'description','Z',ier)
         CALL h5dclose_f(temp_gid,ier)


         ALLOCATE(r2dtemp(nr_fida,nz_fida))
         r2dtemp = SPREAD(raxis_fida, 2, nz_fida)
         CALL write_var_hdf5(fid,'r2d',nr_fida, nz_fida, ier,DBLVAR=DBLE(r2dtemp*100))
         CALL h5dopen_f(fid, 'r2d', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','cm',ier)
         CALL write_att_hdf5(temp_gid,'description','Radius grid: R(r,z)',ier)
         CALL h5dclose_f(temp_gid,ier)
         r2dtemp = SPREAD(zaxis, 1, nr_fida)
         CALL write_var_hdf5(fid,'z2d',nr_fida, nz_fida, ier,DBLVAR=DBLE(r2dtemp*100))
         CALL h5dopen_f(fid, 'z2d', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','cm',ier)
         CALL write_att_hdf5(temp_gid,'description','Z grid: Z(r,z)',ier)
         CALL h5dclose_f(temp_gid,ier)
         DEALLOCATE(r2dtemp)


         ! Close file
         CALL close_hdf5(fid,ier)
         IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'fidasim_'//TRIM(id_string)//'_distribution.h5',ier)

         CALL open_hdf5('fidasim_'//TRIM(id_string)//'_equilibrium.h5',fid,ier,LCREATE=.true.)
         IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'fidasim_'//TRIM(id_string)//'_equilibrium.h5',ier)


         !--------------------------------------------------------------
         !           EQUILIBRIUM - FIELDS
         !--------------------------------------------------------------
         CALL h5gcreate_f(fid,'fields', qid_gid, ier)
         !CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
         CALL write_att_hdf5(qid_gid,'data_source','Data initialized from BEAMS3D',ier)


         !           GRID
         CALL write_var_hdf5(qid_gid,'nr',ier,INTVAR=nr_fida)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nr',ier)
         CALL h5dopen_f(qid_gid, 'nr', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','-',ier)
         CALL write_att_hdf5(temp_gid,'description','Number of R values',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'nphi',ier,INTVAR=nphi_fida)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nphi',ier)
         CALL h5dopen_f(qid_gid, 'nphi', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','-',ier)
         CALL write_att_hdf5(temp_gid,'description','Number of Phi values',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'nz',ier,INTVAR=nz_fida)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nz',ier)
         CALL h5dopen_f(qid_gid, 'nz', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','-',ier)
         CALL write_att_hdf5(temp_gid,'description','Number of Z values',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'time',ier,DBLVAR=DBLE(0)) !!!!Assumes steady state/dummy
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'time',ier)
         CALL h5dopen_f(qid_gid, 'time', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','s',ier)
         CALL write_att_hdf5(temp_gid,'description','Distribution time',ier)
         CALL h5dclose_f(temp_gid,ier)

         CALL write_var_hdf5(qid_gid,'r',nr_fida, ier,DBLVAR=DBLE(raxis_fida*100)) !convert from m to cm
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'r',ier)
         CALL h5dopen_f(qid_gid, 'r', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','cm',ier)
         CALL write_att_hdf5(temp_gid,'description','Radius',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'phi',nphi_fida, ier,DBLVAR=phiaxis_fida)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'phi',ier)
         CALL h5dopen_f(qid_gid, 'phi', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','rad',ier)
         CALL write_att_hdf5(temp_gid,'description','Toroidal angle',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'z',nz_fida, ier,DBLVAR=DBLE(zaxis_fida*100)) !convert from m to cm
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'z',ier)
         CALL h5dopen_f(qid_gid, 'z', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','cm',ier)
         CALL write_att_hdf5(temp_gid,'description','Z',ier)
         CALL h5dclose_f(temp_gid,ier)


         ALLOCATE(r2dtemp(nr_fida,nz_fida))
         r2dtemp = SPREAD(raxis_fida, 2, nz_fida)
         CALL write_var_hdf5(qid_gid,'r2d',nr_fida, nz_fida, ier,DBLVAR=DBLE(r2dtemp*100))
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'r2d',ier)
         CALL h5dopen_f(qid_gid, 'r2d', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','cm',ier)
         CALL write_att_hdf5(temp_gid,'description','Radius grid: R(r,z)',ier)
         CALL h5dclose_f(temp_gid,ier)
         r2dtemp = SPREAD(zaxis_fida, 1, nr_fida)
         CALL write_var_hdf5(qid_gid,'z2d',nr_fida, nz_fida, ier,DBLVAR=DBLE(r2dtemp*100))
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'z2d',ier)
         CALL h5dopen_f(qid_gid, 'z2d', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','cm',ier)
         CALL write_att_hdf5(temp_gid,'description','Z grid: Z(r,z)',ier)
         CALL h5dclose_f(temp_gid,ier)
         DEALLOCATE(r2dtemp)

         !           MASK
         ALLOCATE(mask(nr_fida,nz_fida))
         mask = 1
         CALL write_var_hdf5(qid_gid,'mask',nr_fida,nz_fida,ier,INTVAR=mask) !PLACEHOLDER, should be "Boolean mask that indicates where the fields are well defined", Dim: [nr,nz]
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'mask',ier)
         CALL h5dopen_f(qid_gid, 'mask', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'description','Boolean mask that indicates where the fields are well defined',ier)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'mask',ier)
         CALL h5dclose_f(temp_gid,ier)
         DEALLOCATE(mask)

         !--------------------------------------------------------------
         !           B-FIELD
         !  NOTE On PSI:
         !     In B_STS B_STS_eval_rho defines psi as
         !         rho[0] = sqrt( (psi - Bdata->psi0) / delta );
         !     So psi is TOROIDAL FLUX.
         !--------------------------------------------------------------


         ALLOCATE(rtemp(nr_fida,nz_fida, nphi_fida))
         ALLOCATE(rtemp2(nr_fida,nz_fida, nphi_fida))
         ALLOCATE(rtemp3(nr_fida,nz_fida, nphi_fida))
         DO l = 1,nr_fida
            DO n = 1,nz_fida
               DO m = 1,nphi_fida            
                  ! Eval Spline
                  i = MIN(MAX(COUNT(raxis < raxis_fida(l)),1),nr-1)
                  j = MIN(MAX(COUNT(phiaxis < phiaxis_fida(m)),1),nphi-1)
                  k = MIN(MAX(COUNT(zaxis < zaxis_fida(n)),1),nz-1)
                  xparam = (raxis_fida(l) - raxis(i)) * hri(i)
                  yparam = (phiaxis_fida(m) - phiaxis(j)) * hpi(j)
                  zparam = (zaxis_fida(n) - zaxis(k)) * hzi(k)
                  CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                                 hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                                 BR4D(1,1,1,1),nr,nphi,nz)
                  rtemp(l,n,m) = fval(1)
                  CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                                 hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                                 BPHI4D(1,1,1,1),nr,nphi,nz)
                  rtemp2(l,n,m) = fval(1)
                  CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                                 hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                                 BZ4D(1,1,1,1),nr,nphi,nz)
                  rtemp3(l,n,m) = fval(1)

               END DO
            END DO
         END DO
         

         CALL write_var_hdf5(qid_gid,'br',nr_fida,nz_fida,nphi_fida,ier,DBLVAR=rtemp)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'br',ier)
         CALL h5dopen_f(qid_gid, 'br', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','T',ier)
         CALL write_att_hdf5(temp_gid,'description','Magnetic field in the r-direction: Br(r,z,phi)',ier)
         CALL h5dclose_f(temp_gid,ier)

         CALL write_var_hdf5(qid_gid,'bt',nr_fida,nz_fida,nphi_fida,ier,DBLVAR=rtemp2)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'bt',ier)
         CALL h5dopen_f(qid_gid, 'bt', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','T',ier)
         CALL write_att_hdf5(temp_gid,'description','Magnetic field in the theta/torodial-direction: Bt(r,z,phi)',ier)
         CALL h5dclose_f(temp_gid,ier)

         CALL write_var_hdf5(qid_gid,'bz',nr_fida,nz_fida,nphi_fida,ier,DBLVAR=rtemp3)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'bz',ier)
         CALL h5dopen_f(qid_gid, 'bz', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','T',ier)
         CALL write_att_hdf5(temp_gid,'description','Magnetic field in the z-direction: Bz(r,z,phi)',ier)
         CALL h5dclose_f(temp_gid,ier)
         DEALLOCATE(rtemp)
         DEALLOCATE(rtemp2)
         DEALLOCATE(rtemp3)

         ! ALLOCATE(rtemp(nr,nz, nphi))
         ! rtemp = reshape(B_R(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/)) !switch phi and z
         ! CALL write_var_hdf5(qid_gid,'br',nr,nz,nphi,ier,DBLVAR=rtemp)
         ! CALL h5dopen_f(qid_gid, 'br', temp_gid, ier)
         ! CALL write_att_hdf5(temp_gid,'units','T',ier)
         ! CALL write_att_hdf5(temp_gid,'description','Magnetic field in the r-direction: Br(r,z,phi)',ier)
         ! CALL h5dclose_f(temp_gid,ier)

         ! rtemp = reshape(B_PHI(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/))
         ! CALL write_var_hdf5(qid_gid,'bt',nr,nz,nphi,ier,DBLVAR=rtemp)
         ! CALL h5dopen_f(qid_gid, 'bt', temp_gid, ier)
         ! CALL write_att_hdf5(temp_gid,'units','T',ier)
         ! CALL write_att_hdf5(temp_gid,'description','Magnetic field in the theta/torodial-direction: Bt(r,z,phi)',ier)
         ! CALL h5dclose_f(temp_gid,ier)

         ! rtemp = reshape(B_Z(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/))
         ! CALL write_var_hdf5(qid_gid,'bz',nr,nz,nphi,ier,DBLVAR=rtemp)
         ! CALL h5dopen_f(qid_gid, 'bz', temp_gid, ier)
         ! CALL write_att_hdf5(temp_gid,'units','T',ier)
         ! CALL write_att_hdf5(temp_gid,'description','Magnetic field in the z-direction: Bz(r,z,phi)',ier)
         ! CALL h5dclose_f(temp_gid,ier)
         ! DEALLOCATE(rtemp)

         !--------------------------------------------------------------
         !           E-FIELD - NEEDS CHECKING
         !--------------------------------------------------------------

         ! Values must be equidistant in rho.
         IF (npot < 1) THEN ! Because we can run with E=0
            ALLOCATE(rtemp(nr_fida,nz_fida,nphi_fida))
            rtemp = 0.0
            CALL write_var_hdf5(qid_gid,'er',nr_fida,nz_fida,nphi_fida, ier,DBLVAR=rtemp)
            CALL write_var_hdf5(qid_gid,'et',nr_fida,nz_fida,nphi_fida, ier,DBLVAR=rtemp)
            CALL write_var_hdf5(qid_gid,'ez',nr_fida,nz_fida,nphi_fida, ier,DBLVAR=rtemp)
            DEALLOCATE(rtemp)
         ELSE
            ALLOCATE(r4dtemp(nr_fida,nphi_fida,nz_fida,3))
            ALLOCATE(r1dtemp(3))
            r1dtemp = 1
            DO l = 1,nr_fida
               DO n = 1,nz_fida
                  DO m = 1,nphi_fida            
                     ! Eval Spline
                     ! Get the gridpoint info (this is possible since all grids are the same)
                     i = MIN(MAX(COUNT(raxis < raxis_fida(l)),1),nr-1)
                     j = MIN(MAX(COUNT(phiaxis < phiaxis_fida(m)),1),nphi-1)
                     k = MIN(MAX(COUNT(zaxis < zaxis_fida(n)),1),nz-1)
                     xparam = (raxis_fida(l) - raxis(i)) * hri(i)
                     yparam = (phiaxis_fida(m) - phiaxis(j)) * hpi(j)
                     zparam = (zaxis_fida(n) - zaxis(k)) * hzi(k)
                     ! Evaluate the Splines
                     CALL R8HERM3FCN(ictE,1,1,fvalE,i,j,k,xparam,yparam,zparam,& !evaluate at grid points
                        hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                        POT4D(1,1,1,1),nr,nphi,nz)
                     r1dtemp(1:3) =-fvalE(1,1:3)
                     r4dtemp(l,n,m,1:3) = -r1dtemp(1:3)
                  END DO
               END DO
            END DO
            CALL write_var_hdf5(qid_gid,'er',nr_fida,nz_fida,nphi_fida, ier,DBLVAR=r4dtemp(:,:,:,1))
            IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'er',ier)
            CALL h5dopen_f(qid_gid, 'er', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','V/m',ier)
            CALL write_att_hdf5(temp_gid,'description','Electric field in the r-direction: Er(r,z,phi)',ier)
            CALL h5dclose_f(temp_gid,ier)

            CALL write_var_hdf5(qid_gid,'et',nr_fida,nz_fida,nphi_fida, ier,DBLVAR=r4dtemp(:,:,:,2))
            IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'et',ier)
            CALL h5dopen_f(qid_gid, 'et', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','V/m',ier)
            CALL write_att_hdf5(temp_gid,'description','Electric field in the toroidal phi-direction: Et(r,z,phi)',ier)
            CALL h5dclose_f(temp_gid,ier)

            CALL write_var_hdf5(qid_gid,'ez',nr_fida,nz_fida,nphi_fida, ier,DBLVAR=r4dtemp(:,:,:,3))
            IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ez',ier)
            CALL h5dopen_f(qid_gid, 'ez', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','V/m',ier)
            CALL write_att_hdf5(temp_gid,'description','Electric field in the z-direction: Ez(r,z,phi)',ier)
            CALL h5dclose_f(temp_gid,ier)
            DEALLOCATE(r1dtemp)
            DEALLOCATE(r4dtemp)
         END IF

         CALL h5gclose_f(qid_gid, ier)

         !--------------------------------------------------------------
         !           PLASMA
         !--------------------------------------------------------------
         CALL h5gcreate_f(fid,'plasma', qid_gid, ier)
         !CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
         CALL write_att_hdf5(qid_gid,'data_source','Data initialized from BEAMS3D ',ier)
         CALL write_att_hdf5(qid_gid,'description','no bulk plasma rotation/flow',ier)

         !           GRID
         CALL write_var_hdf5(qid_gid,'nr',ier,INTVAR=nr_fida)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nr',ier)
         CALL h5dopen_f(qid_gid, 'nr', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','-',ier)
         CALL write_att_hdf5(temp_gid,'description','Number of R values',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'nphi',ier,INTVAR=nphi_fida)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nphi',ier)
         CALL h5dopen_f(qid_gid, 'nphi', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','-',ier)
         CALL write_att_hdf5(temp_gid,'description','Number of Phi values',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'nz',ier,INTVAR=nz_fida)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nz',ier)
         CALL h5dopen_f(qid_gid, 'nz', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','-',ier)
         CALL write_att_hdf5(temp_gid,'description','Number of Z values',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'time',ier,DBLVAR=DBLE(0)) !!!!Assumes steady state/dummy
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'time',ier)
         CALL h5dopen_f(qid_gid, 'time', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','s',ier)
         CALL write_att_hdf5(temp_gid,'description','Distribution time',ier)
         CALL h5dclose_f(temp_gid,ier)

         CALL write_var_hdf5(qid_gid,'r',nr_fida, ier,DBLVAR=DBLE(raxis_fida*100)) !convert from m to cm
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'r',ier)
         CALL h5dopen_f(qid_gid, 'r', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','cm',ier)
         CALL write_att_hdf5(temp_gid,'description','Radius',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'phi',nphi_fida, ier,DBLVAR=phiaxis_fida)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'phi',ier)
         CALL h5dopen_f(qid_gid, 'phi', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','rad',ier)
         CALL write_att_hdf5(temp_gid,'description','Toroidal angle',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'z',nz_fida, ier,DBLVAR=DBLE(zaxis_fida*100)) !convert from m to cm
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'z',ier)
         CALL h5dopen_f(qid_gid, 'z', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','cm',ier)
         CALL write_att_hdf5(temp_gid,'description','Z',ier)
         CALL h5dclose_f(temp_gid,ier)


         ALLOCATE(r2dtemp(nr_fida,nz_fida))
         r2dtemp = SPREAD(raxis_fida, 2, nz_fida)
         CALL write_var_hdf5(qid_gid,'r2d',nr_fida, nz_fida, ier,DBLVAR=DBLE(r2dtemp*100))
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'r2d',ier)
         CALL h5dopen_f(qid_gid, 'r2d', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','cm',ier)
         CALL write_att_hdf5(temp_gid,'description','Radius grid: R(r,z)',ier)
         CALL h5dclose_f(temp_gid,ier)
         r2dtemp = SPREAD(zaxis_fida, 1, nr_fida)
         CALL write_var_hdf5(qid_gid,'z2d',nr_fida, nz_fida, ier,DBLVAR=DBLE(r2dtemp*100))
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'z2d',ier)
         CALL h5dopen_f(qid_gid, 'z2d', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','cm',ier)
         CALL write_att_hdf5(temp_gid,'description','Z grid: Z(r,z)',ier)
         CALL h5dclose_f(temp_gid,ier)
         DEALLOCATE(r2dtemp)

         !           MASK
         ALLOCATE(mask(nr_fida,nz_fida))
         mask = 1
         CALL write_var_hdf5(qid_gid,'mask',nr_fida,nz_fida,ier,INTVAR=mask) !PLACEHOLDER, should be "Boolean mask that indicates where the fields are well defined", Dim: [nr,nz]
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'mask',ier)
         CALL h5dopen_f(qid_gid, 'mask', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'description','Boolean mask that indicates where the fields are well defined',ier)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'mask',ier)
         CALL h5dclose_f(temp_gid,ier)
         DEALLOCATE(mask)

         !           PLASMA ROTATION/FLOW
         ALLOCATE(rtemp(nr_fida,nz_fida,nphi_fida))
         rtemp = 0.0
         CALL write_var_hdf5(qid_gid,'vr',nr_fida,nz_fida, nphi_fida,ier,DBLVAR=rtemp)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'vr',ier)
         CALL h5dopen_f(qid_gid, 'vr', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'description','Bulk plasma flow in the r-direction: Vr(r,z,phi)',ier)
         CALL write_att_hdf5(temp_gid,'units','cm/s',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'vt',nr_fida,nz_fida, nphi_fida,ier,DBLVAR=rtemp)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'vt',ier)
         CALL h5dopen_f(qid_gid, 'vt', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'description','Bulk plasma flow in the toroidal phi-direction: Vphi(r,z,phi)',ier)
         CALL write_att_hdf5(temp_gid,'units','cm/s',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(qid_gid,'vz',nr_fida,nz_fida, nphi_fida,ier,DBLVAR=rtemp)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'vz',ier)
         CALL h5dopen_f(qid_gid, 'vz', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'description','Bulk plasma flow in the z-direction: Vz(r,z,phi)',ier)
         CALL write_att_hdf5(temp_gid,'units','cm/s',ier)
         CALL h5dclose_f(temp_gid,ier)
         DEALLOCATE(rtemp)

         ALLOCATE(rtemp(nr_fida,nz_fida, nphi_fida))
         ALLOCATE(rtemp2(nr_fida,nz_fida, nphi_fida))
         ALLOCATE(rtemp3(nr_fida,nz_fida, nphi_fida))
         ALLOCATE(rtemp4(nr_fida,nz_fida, nphi_fida))

         DO l = 1,nr_fida
            DO n = 1,nz_fida
               DO m = 1,nphi_fida            
                  ! Eval Spline
                  ! Get the gridpoint info (this is possible since all grids are the same)
                  i = MIN(MAX(COUNT(raxis < raxis_fida(l)),1),nr-1)
                  j = MIN(MAX(COUNT(phiaxis < phiaxis_fida(m)),1),nphi-1)
                  k = MIN(MAX(COUNT(zaxis < zaxis_fida(n)),1),nz-1)
                  xparam = (raxis_fida(l) - raxis(i)) * hri(i)
                  yparam = (phiaxis_fida(m) - phiaxis(j)) * hpi(j)
                  zparam = (zaxis_fida(n) - zaxis(k)) * hzi(k)
                  ! Evaluate the Splines
                  CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                                 hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                                 TE4D(1,1,1,1),nr,nphi,nz)
                  rtemp(l,n,m) = max(fval(1),zero)
                  CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                                 hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                                 NE4D(1,1,1,1),nr,nphi,nz)
                  rtemp2(l,n,m) = max(fval(1),zero)
                  CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                                 hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                                 TI4D(1,1,1,1),nr,nphi,nz)
                  rtemp3(l,n,m) = max(fval(1),zero)
                  CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                                 hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                                 ZEFF4D(1,1,1,1),nr,nphi,nz)
                  rtemp4(l,n,m) = max(fval(1),one)
               END DO
            END DO
         END DO

         CALL write_var_hdf5(qid_gid,'te',nr_fida,nz_fida,nphi_fida, ier,DBLVAR=DBLE(rtemp/1000))
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'te',ier)
         CALL h5dopen_f(qid_gid, 'te', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','kev',ier)
         CALL write_att_hdf5(temp_gid,'description','Electron Temperature: Ti(r,z,phi)',ier)
         CALL h5dclose_f(temp_gid,ier)

         CALL write_var_hdf5(qid_gid,'dene',nr_fida,nz_fida,nphi_fida, ier,DBLVAR=DBLE(rtemp2/1000000)) !convert from m^-3 to cm^-3
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'dene',ier)
         CALL h5dopen_f(qid_gid, 'dene', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','cm^-3',ier)
         CALL write_att_hdf5(temp_gid,'description','Electron Number Density: Dene(r,z, phi)',ier)
         CALL h5dclose_f(temp_gid,ier)

         CALL write_var_hdf5(qid_gid,'ti',nr_fida,nz_fida,nphi_fida, ier,DBLVAR=DBLE(rtemp3/1000))
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ti',ier)
         CALL h5dopen_f(qid_gid, 'ti', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','keV',ier)
         CALL write_att_hdf5(temp_gid,'description','Ion Temperature: Ti(r,z,phi)',ier)
         CALL h5dclose_f(temp_gid,ier)

         CALL write_var_hdf5(qid_gid,'zeff',nr_fida,nz_fida,nphi_fida, ier,DBLVAR=DBLE(rtemp4))
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'zeff',ier)
         CALL h5dopen_f(qid_gid, 'zeff', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'units','-',ier)
         CALL write_att_hdf5(temp_gid,'description','Effective Nuclear Charge: Zeff(r,z,phi)',ier)
         CALL h5dclose_f(temp_gid,ier)
         DEALLOCATE(rtemp)
         DEALLOCATE(rtemp2)
         DEALLOCATE(rtemp3)
         DEALLOCATE(rtemp4)

         ! rtemp = reshape(NE(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/))
         ! CALL write_var_hdf5(qid_gid,'dene',nr,nz,nphi, ier,DBLVAR=DBLE(rtemp/1000000)) !convert from m^-3 to cm^-3
         ! CALL h5dopen_f(qid_gid, 'dene', temp_gid, ier)
         ! CALL write_att_hdf5(temp_gid,'units','cm^-3',ier)
         ! CALL write_att_hdf5(temp_gid,'description','Electron Number Density: Dene(r,z, phi)',ier)
         ! CALL h5dclose_f(temp_gid,ier)
         ! rtemp = reshape(TE(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/))
         ! CALL write_var_hdf5(qid_gid,'te',nr,nz,nphi, ier,DBLVAR=DBLE(rtemp/1000))
         ! CALL h5dopen_f(qid_gid, 'te', temp_gid, ier)
         ! CALL write_att_hdf5(temp_gid,'units','kev',ier)
         ! CALL write_att_hdf5(temp_gid,'description','Electron Temperature: Ti(r,z,phi)',ier)
         ! CALL h5dclose_f(temp_gid,ier)
         ! rtemp = reshape(TI(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/))
         ! CALL write_var_hdf5(qid_gid,'ti',nr,nz,nphi, ier,DBLVAR=DBLE(rtemp/1000))
         ! CALL h5dopen_f(qid_gid, 'ti', temp_gid, ier)
         ! CALL write_att_hdf5(temp_gid,'units','keV',ier)
         ! CALL write_att_hdf5(temp_gid,'description','Ion Temperature: Ti(r,z,phi)',ier)
         ! CALL h5dclose_f(temp_gid,ier)
         ! rtemp = reshape(ZEFF_ARR(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/))
         ! CALL write_var_hdf5(qid_gid,'zeff',nr,nz,nphi, ier,DBLVAR=DBLE(rtemp))
         ! CALL h5dopen_f(qid_gid, 'zeff', temp_gid, ier)
         ! CALL write_att_hdf5(temp_gid,'units','-',ier)
         ! CALL write_att_hdf5(temp_gid,'description','Effective Nuclear Charge: Zeff(r,z,phi)',ier)
         ! CALL h5dclose_f(temp_gid,ier)
         ! DEALLOCATE(rtemp)


         CALL h5gclose_f(qid_gid, ier)
         ! Close file
         CALL close_hdf5(fid,ier)
         IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'fidasim_'//TRIM(id_string)//'_equilibrium.h5',ier)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       CASE('DENF')
         ! Distribution needs to be in density [part/m^3] here - called from distnorm (before velocity space normalization)
         !ALLOCATE(r4dtemp(nbeams,nr_fida,nphi_fida,nz_fida))
         ALLOCATE(rtemp(nr_fida,nz_fida,nphi_fida))
         rtemp = 0.0
         IF (lfidasim2) THEN

            DO i=1,nr_fida
               vol = r_h * z_h * p_h * (raxis_fida(i) + r_h / 2.0)
               !WRITE(327,*) i, vol
               !CALL FLUSH(327)
               dist5d_fida(:,i,:,:,:,:) = dist5d_fida(:,i,:,:,:,:) / vol
               DO j = 1, nz_fida
                  DO k=1,nphi_fida
                     rtemp(i,j,k) = SUM(dist5d_fida(:,i,j,k,:,:))
                  END DO
               END DO
            END DO
         ELSE 
         !convert to r z phi
         DO i=1,nr_fida
            DO k = 1, nz_fida
               DO j=1,nphi_fida
                  !convert i,j,k to distribution functionlfidasim indices l,m,n
                  !determine beams3d-grid indices
                  i3 = MIN(MAX(COUNT(raxis < raxis_fida(i)),1),nr-1)
                  j3 = MIN(MAX(COUNT(phiaxis < phiaxis_fida(j)),1),nphi-1)
                  k3 = MIN(MAX(COUNT(zaxis < zaxis_fida(k)),1),nz-1)
                  !setup interpolation
                  xparam = (raxis_fida(i) - raxis(i3)) * hri(i3)
                  yparam = (phiaxis_fida(j) - phiaxis(j3)) * hpi(j3)
                  zparam = (zaxis_fida(k) - zaxis(k3)) * hzi(k3)
                  CALL R8HERM3FCN(ict,1,1,fval,i3,j3,k3,xparam,yparam,zparam,&!maybe switch to x/y interpolation?
                     hr(i3),hri(i3),hp(j3),hpi(j3),hz(k3),hzi(k3),&
                     S4D(1,1,1,1),nr,nphi,nz)
                  y0 = fval(1)
                  CALL R8HERM3FCN(ict,1,1,fval2,i3,j3,k3,xparam,yparam,zparam,&
                     hr(i3),hri(i3),hp(j3),hpi(j3),hz(k3),hzi(k3),&
                     U4D(1,1,1,1),nr,nphi,nz)
                  z0 = fval2(1)

                  IF (z0 < 0) z0 = z0 + pi2
                  x0    = phiaxis_fida(j)
                  IF (x0 < 0) x0 = x0 + pi2
                  ! Calc dist func bins
                  l = MAX(MIN(CEILING(SQRT(y0)*ns_prof1     ), ns_prof1), 1) ! Rho Bin
                  m = MAX(MIN(CEILING( z0*h2_prof           ), ns_prof2), 1) ! U Bin
                  n = MAX(MIN(CEILING( x0*h3_prof           ), ns_prof3), 1) ! V Bin
                   IF (y0 .GT. 1.05) THEN
                      rtemp(i,k,j) = 0 !distribution is 0 outside plasma
                   ELSE
                     rtemp(i,k,j) = SUM(dist5d_prof(:,l,m,n,:,:))!output in r-z-phi
                   END IF
               END DO
            END DO
         END DO
      END IF
         CALL open_hdf5('fidasim_'//TRIM(id_string)//'_distribution.h5',fid,ier,LCREATE=.false.)
         CALL h5gopen_f(fid,'/', qid_gid, ier)
         CALL write_var_hdf5(qid_gid,'denf',nr_fida,nz_fida,nphi_fida,ier,DBLVAR=DBLE(rtemp/1000000.0)) !in cm^3
         CALL h5dopen_f(qid_gid, 'denf', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'description','Fast-ion density (nr_fida,nz_fida,nphi_fida)',ier)
         CALL write_att_hdf5(temp_gid,'units','part/(cm^3)',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL h5gclose_f(qid_gid, ier)
         DEALLOCATE(rtemp)
         CALL close_hdf5(fid,ier)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       CASE('DISTRIBUTION_GC_F')
         ALLOCATE(pitch_fida(npitch_fida))
         ALLOCATE(energy_fida(nenergy_fida))
         ! Do volume normalization
         !CALL beams3d_distnorm !TODO: check if this has already been done
         !jac = MAXVAL(mass_beams);
         FORALL(i = 1:nenergy_fida) energy_fida(i) = (i) / REAL(nenergy_fida+1) * 0.5 * mass_beams(1) * partvmax * partvmax /e_charge / 1000.0 !Potential error when different beam species are used!
         FORALL(i = 1:npitch_fida) pitch_fida(i) = (i) / REAL(npitch_fida+1) * 2.0 - 1.0

         CALL open_hdf5('fidasim_'//TRIM(id_string)//'_distribution.h5',fid,ier,LCREATE=.false.)
         CALL write_var_hdf5(fid,'energy',nenergy_fida,ier,DBLVAR=energy_fida) ! in keV
         CALL h5dopen_f(fid, '/energy', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'description','Energy array',ier)
         CALL write_att_hdf5(temp_gid,'units','keV',ier)
         CALL h5dclose_f(temp_gid,ier)
         CALL write_var_hdf5(fid,'pitch',npitch_fida,ier,DBLVAR=pitch_fida)
         CALL h5dopen_f(fid, '/pitch', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'description','Pitch array',ier)
         CALL write_att_hdf5(temp_gid,'units','-',ier)
         CALL h5dclose_f(temp_gid,ier)

         ! Do phase space change of coordinates
         !Allocate with Radial-like dimensions for clean transfer and to avoid explicitly looping over every element
         
         IF (lfidasim2) THEN
            ALLOCATE(dist5d_temp(nbeams, nenergy_fida, npitch_fida,nr_fida,nz_fida,nphi_fida)) !need temp as velocity bins are in vll/vperp initially
            dist5d_temp = 0
         ELSE
            ALLOCATE(dist5d_fida(nbeams,ns_prof1, ns_prof2, ns_prof3, nenergy_fida, npitch_fida)) !nenergy and npitch are always aligned to distribution
         END IF
         DO b=1,nbeams
            DO d1 = 1, nenergy_fida
               DO d2 = 1, npitch_fida
                  v = SQRT(2 * energy_fida(d1) *1000.0 * e_charge / mass_beams(b))
                  IF (v .gt. partvmax) THEN !) .or. (v .eq. 0.0))
                     IF (.not. lfidasim2) dist5d_fida(b,:,:,:,d1,d2) = 0
                  ELSE
                     pitch = pitch_fida(d2)
                     v_parr = pitch * v
                     v_perp = SQRT(1- pitch * pitch) * v
                     !determine beams3d-grid indices (velocity space)
                     i3 = MAX(MIN(1+nsh_prof4+FLOOR(h4_prof*v_parr), ns_prof4), 1) ! vll
                     j3 = MAX(MIN(CEILING(v_perp*h5_prof         ), ns_prof5), 1) ! Vperp
                     jac = 1 / (mass_beams(b) * SQRT(1-pitch * pitch))
                     !jac = MIN(MAX(jac, 0.0),1.0E)
                     !jac = jac / pi / REAL(1000000) * e_charge / 1000 !convert to TRANSP convention and 1/cm^3/keV
                     jac = v / mass_beams(b) * e_charge / REAL(1000) ! * pi2
                     IF (lfidasim2) THEN
                        dist5d_temp(b,d1,d2,:,:,:) = dist5d_fida(b,:,:,:,i3,j3) * jac
                     ELSE   
                        dist5d_fida(b,:,:,:,d1,d2) = dist5d_prof(b,:,:,:,i3,j3) * jac ! conversion to final grid comes in next steps
                     END IF

                  END IF
               END DO
            END DO
         END DO

         IF (.not. lfidasim2) THEN
         ALLOCATE(dist5d_temp(nbeams,ns_prof1, ns_prof2, ns_prof3, nenergy_fida, npitch_fida))
         dist5d_temp(:,:,:,:,:,:) = dist5d_fida(:,:,:,:,:,:)
         DEALLOCATE(dist5d_fida) 
         !Interpolate rho u v to r phi z distribution function (nearest neighbor at the moment)
         !Now allocate with correct dimensions
         ALLOCATE(dist5d_fida(nbeams,nenergy_fida,npitch_fida,nr_fida,nz_fida,nphi_fida))
         DO b=1,nbeams
            DO i=1,nr_fida
               DO k = 1, nz_fida
                  DO j=1,nphi_fida
                     !convert i,j,k to distribution function indices l,m,n
                     !determine beams3d-grid indices
                     i3 = MIN(MAX(COUNT(raxis < raxis_fida(i)),1),nr-1)
                     j3 = MIN(MAX(COUNT(phiaxis < phiaxis_fida(j)),1),nphi-1)
                     k3 = MIN(MAX(COUNT(zaxis < zaxis_fida(k)),1),nz-1)
                     !setup interpolation
                     xparam = (raxis_fida(i) - raxis(i3)) * hri(i3)
                     yparam = (phiaxis_fida(j) - phiaxis(j3)) * hpi(j3)
                     zparam = (zaxis_fida(k) - zaxis(k3)) * hzi(k3)
                     CALL R8HERM3FCN(ict,1,1,fval,i3,j3,k3,xparam,yparam,zparam,&!maybe switch to x/y interpolation?
                        hr(i3),hri(i3),hp(j3),hpi(j3),hz(k3),hzi(k3),&
                        S4D(1,1,1,1),nr,nphi,nz)
                     y0 = fval(1)
                     CALL R8HERM3FCN(ict,1,1,fval2,i3,j3,k3,xparam,yparam,zparam,&
                        hr(i3),hri(i3),hp(j3),hpi(j3),hz(k3),hzi(k3),&
                        U4D(1,1,1,1),nr,nphi,nz)
                     z0 = fval2(1)

                     IF (z0 < 0) z0 = z0 + pi2
                     IF (x0 < 0) x0 = x0 + pi2
                     IF (y0 < 0) y0 = -y0
                     ! Calc dist func bins
                     x0    = phiaxis_fida(j)
                     l = MAX(MIN(CEILING(SQRT(y0)*ns_prof1     ), ns_prof1), 1) ! Rho Bin
                     m = MAX(MIN(CEILING( z0*h2_prof           ), ns_prof2), 1) ! U Bin
                     n = MAX(MIN(CEILING( x0*h3_prof           ), ns_prof3), 1) ! V Bin
                     IF (y0 .GT. 1.05) THEN !might introduce a small deviation here
                        dist5d_fida(b,:,:,i,k,j) = 0 !distribution is 0 outside plasma
                     ELSE
                        dist5d_fida(b,:,:,i,k,j) = dist5d_temp(b,l,m,n,:,:) !output in r-z-phi
                     END IF

                  END DO
               END DO
            END DO
         END DO
      END IF

         CALL open_hdf5('fidasim_'//TRIM(id_string)//'_distribution.h5',fid,ier,LCREATE=.false.)
         IF (ASSOCIATED(dist5d_fida)) THEN
            IF (lsplit) THEN
               CALL write_var_hdf5(fid,'f',nbeams,nenergy_fida,npitch_fida,nr_fida,nz_fida,nphi_fida,ier,DBLVAR=dist5d_fida)  
            ELSE IF (lfidasim2) THEN
               CALL write_var_hdf5(fid,'f',nenergy_fida,npitch_fida,nr_fida,nz_fida,nphi_fida,ier,DBLVAR=SUM(dist5d_temp, DIM=1))
            ELSE
               CALL write_var_hdf5(fid,'f',nenergy_fida,npitch_fida,nr_fida,nz_fida,nphi_fida,ier,DBLVAR=SUM(dist5d_fida, DIM=1))
            END IF
            CALL h5dopen_f(fid, '/f', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'description','Distribution Function (nenergy_fida,npitch_fida,nr_fida,nz_fida,nphi_fida)',ier)
            CALL write_att_hdf5(temp_gid,'units','part/(cm^3 keV)',ier)
            CALL h5dclose_f(temp_gid,ier)
            IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'dist_fida',ier)
         END IF
         CALL close_hdf5(fid,ier)
         IF (.not. lfidasim2) DEALLOCATE(dist5d_fida)
         DEALLOCATE(pitch_fida)
         DEALLOCATE(energy_fida)
         !DEALLOCATE(dist5d_temp)


         ! INQUIRE(FILE=TRIM(fidasim_input_dat),EXIST=lexist)
         ! IF (.not.lexist) THEN
         !    istat=-1
         !    WRITE(6,*) " ERROR: Could not find file: "//TRIM(fidasim_input_dat)
         !    RETURN
         ! END IF
         ! INQUIRE(FILE=TRIM(fidasim_geometry),EXIST=lexist)
         ! IF (.not.lexist) THEN
         !    istat=-1
         !    WRITE(6,*) " ERROR: Could not find file: "//TRIM(fidasim_geometry)
         !    RETURN
         ! END IF

       CASE('DISTRIBUTION_GC_MC')
       CASE('DISTRIBUTION_FO')

      END SELECT
   END IF !myworkid==master
   RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
END SUBROUTINE beams3d_write_fidasim
