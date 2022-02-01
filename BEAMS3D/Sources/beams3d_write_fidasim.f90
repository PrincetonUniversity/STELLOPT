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
      USE beams3d_lines
      USE beams3d_grid, ONLY: nr, nphi, nz, B_R, B_PHI, B_Z, raxis, &
                                 zaxis, phiaxis, S_ARR, U_ARR, POT_ARR, &
                                 ZEFF_ARR, TE, TI, NE, req_axis, zeq_axis, npot, &
                                 POT_SPL, ezspline_interp, phiedge_eq, TE_spl, &
                                 NE_spl, TI_spl, ZEFF_spl, nne, nte, nti, nzeff, &
                                 plasma_mass, reff_eq, therm_factor, &
                                 X_BEAMLET, Y_BEAMLET, Z_BEAMLET, &
                                 NX_BEAMLET, NY_BEAMLET, NZ_BEAMLET, &
                                 POT4D, NE4D, TE4D, TI4D, ZEFF4D, &
                                 hr, hp, hz, hri, hpi, hzi
      USE beams3d_runtime, ONLY: id_string, npoinc, nbeams, beam, t_end, lverb, &
                                    lvmec, lpies, lspec, lcoil, lmgrid, lbeam, lplasma_only, &
                                    lvessel, lvac, lbeam_simple, handle_err, nparticles_start, &
                                    HDF5_OPEN_ERR,HDF5_WRITE_ERR,&
                                    HDF5_CLOSE_ERR, BEAMS3D_VERSION, weight, e_beams, p_beams,&
                                    charge, Zatom, mass, ldepo, v_neut, &
                                    lcollision, pi, pi2, t_end_in, nprocs_beams, &
                                    div_beams, mass_beams, Zatom_beams, dex_beams, &
                                    qid_str_saved, lascotfl, R_beams, PHI_beams, Z_beams
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
      INTEGER :: ier, iunit, i, j, d1, d2, d3, k, k1, k2, kmax ,ider, nphi1
      INTEGER(HID_T) :: options_gid, bfield_gid, efield_gid, plasma_gid, &
                        neutral_gid, wall_gid, marker_gid, qid_gid, &
                        nbi_gid, inj_gid, boozer_gid, mhd_gid, temp_gid
      INTEGER, ALLOCATABLE, DIMENSION(:) :: itemp
      REAL :: qid_flt
      REAL*8 :: fvalE(1,3)
      DOUBLE PRECISION :: rho_temp, s_temp, dbl_temp, gammarel, v_total!, hx, hy, hz, hxi, hyi, hzi, xparam, yparam, zparam
      DOUBLE PRECISION, ALLOCATABLE :: rtemp(:,:,:), r1dtemp(:), r2dtemp(:,:), Efield4D(:,:,:,:)
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mask
      CHARACTER(LEN=10) ::  qid_str
      CHARACTER(LEN=8) :: temp_str8, inj_str8

      INTEGER, parameter :: ictE(8)=(/0,1,1,1,0,0,0,0/)
      REAL*8, PARAMETER :: one = 1
      !DOUBLE PRECISION, PARAMETER :: e_charge      = 1.60217662E-19 !e_c
      DOUBLE PRECISION, PARAMETER :: c_speed       = 2.99792458E+08 !Speed of light [m/s]
      DOUBLE PRECISION, PARAMETER :: e_charge      = 1.602176565e-19 !e_c
      DOUBLE PRECISION, PARAMETER :: inv_amu       = 6.02214076208E+26 ! 1./AMU [1/kg]
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      SELECT CASE (TRIM(write_type))
         CASE('INIT')
            IF (myworkid == master) THEN
            nphi1 = nphi-1 ! confirmed with Poincare plot
            CALL open_hdf5('fidasim_'//TRIM(id_string)//'_distribution.h5',fid,ier,LCREATE=.true.)
            CALL h5gopen_f(fid,'/', qid_gid, ier)   
            CALL write_att_hdf5(qid_gid,'data_source','Data initialized from BEAMS3D',ier)
            CALL write_var_hdf5(qid_gid,'type',ier,INTVAR=2)
            CALL h5dopen_f(fid, 'type', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'description','Distribution type: 1="Guiding Center Density Function", 2="Guiding Center Monte Carlo", 3="Full Orbit Monte Carlo"',ier)
            CALL h5dclose_f(temp_gid,ier)
            


            !           GRID
            CALL write_var_hdf5(qid_gid,'nr',ier,INTVAR=nr)
            CALL h5dopen_f(fid, 'nr', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','-',ier)
            CALL write_att_hdf5(temp_gid,'description','Number of R values',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'nphi',ier,INTVAR=nphi)
            CALL h5dopen_f(fid, 'nphi', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','-',ier)
            CALL write_att_hdf5(temp_gid,'description','Number of Phi values',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'nz',ier,INTVAR=nz)
            CALL h5dopen_f(fid, 'nz', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','-',ier)
            CALL write_att_hdf5(temp_gid,'description','Number of Z values',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'time',ier,DBLVAR=DBLE(0)) !!!!Assumes steady state/dummy
            CALL h5dopen_f(fid, 'time', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','s',ier)
            CALL write_att_hdf5(temp_gid,'description','Distribution time',ier)
            CALL h5dclose_f(temp_gid,ier)

            CALL write_var_hdf5(qid_gid,'r',nr, ier,DBLVAR=DBLE(raxis*100)) !convert from m to cm
            CALL h5dopen_f(fid, 'r', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm',ier)
            CALL write_att_hdf5(temp_gid,'description','Radius',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'phi',nphi, ier,DBLVAR=phiaxis)
            CALL h5dopen_f(fid, 'phi', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','rad',ier)
            CALL write_att_hdf5(temp_gid,'description','Toroidal angle',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'z',nz, ier,DBLVAR=DBLE(zaxis*100)) !convert from m to cm
            CALL h5dopen_f(fid, 'z', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm',ier)
            CALL write_att_hdf5(temp_gid,'description','Z',ier)
            CALL h5dclose_f(temp_gid,ier)


            ALLOCATE(r2dtemp(nr,nz))
            r2dtemp = SPREAD(raxis, 2, nz)
            CALL write_var_hdf5(fid,'r2d',nr, nz, ier,DBLVAR=DBLE(r2dtemp*100))
            CALL h5dopen_f(fid, 'r2d', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm',ier)
            CALL write_att_hdf5(temp_gid,'description','Radius grid: R(r,z)',ier)
            CALL h5dclose_f(temp_gid,ier)
            r2dtemp = SPREAD(zaxis, 1, nr)
            CALL write_var_hdf5(fid,'z2d',nr, nz, ier,DBLVAR=DBLE(r2dtemp*100))
            CALL h5dopen_f(fid, 'z2d', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm',ier)
            CALL write_att_hdf5(temp_gid,'description','Z grid: Z(r,z)',ier)
            CALL h5dclose_f(temp_gid,ier)
            DEALLOCATE(r2dtemp)


            ! Close file
            CALL close_hdf5(fid,ier)
            IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'ascot5_'//TRIM(id_string)//'_distribution.h5',ier)
               
            CALL open_hdf5('fidasim_'//TRIM(id_string)//'_equilibrium.h5',fid,ier,LCREATE=.true.)
            IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'ascot5_'//TRIM(id_string)//'_equilibrium.h5',ier)


            !--------------------------------------------------------------
            !           EQUILIBRIUM - FIELDS
            !--------------------------------------------------------------
            CALL h5gcreate_f(fid,'fields', qid_gid, ier)
            !CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
            CALL write_att_hdf5(qid_gid,'data_source','Data initialized from BEAMS3D',ier)


           !           GRID
            CALL write_var_hdf5(qid_gid,'nr',ier,INTVAR=nr)
            CALL h5dopen_f(qid_gid, 'nr', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','-',ier)
            CALL write_att_hdf5(temp_gid,'description','Number of R values',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'nphi',ier,INTVAR=nphi)
            CALL h5dopen_f(qid_gid, 'nphi', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','-',ier)
            CALL write_att_hdf5(temp_gid,'description','Number of Phi values',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'nz',ier,INTVAR=nz)
            CALL h5dopen_f(qid_gid, 'nz', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','-',ier)
            CALL write_att_hdf5(temp_gid,'description','Number of Z values',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'time',ier,DBLVAR=DBLE(0)) !!!!Assumes steady state/dummy
            CALL h5dopen_f(qid_gid, 'time', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','s',ier)
            CALL write_att_hdf5(temp_gid,'description','Distribution time',ier)
            CALL h5dclose_f(temp_gid,ier)

            CALL write_var_hdf5(qid_gid,'r',nr, ier,DBLVAR=DBLE(raxis*100)) !convert from m to cm
            CALL h5dopen_f(qid_gid, 'r', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm',ier)
            CALL write_att_hdf5(temp_gid,'description','Radius',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'phi',nphi, ier,DBLVAR=phiaxis)
            CALL h5dopen_f(qid_gid, 'phi', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','rad',ier)
            CALL write_att_hdf5(temp_gid,'description','Toroidal angle',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'z',nz, ier,DBLVAR=DBLE(zaxis*100)) !convert from m to cm
            CALL h5dopen_f(qid_gid, 'z', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm',ier)
            CALL write_att_hdf5(temp_gid,'description','Z',ier)
            CALL h5dclose_f(temp_gid,ier)


            ALLOCATE(r2dtemp(nr,nz))
            r2dtemp = SPREAD(raxis, 2, nz)
            CALL write_var_hdf5(qid_gid,'r2d',nr, nz, ier,DBLVAR=DBLE(r2dtemp*100))
            CALL h5dopen_f(qid_gid, 'r2d', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm',ier)
            CALL write_att_hdf5(temp_gid,'description','Radius grid: R(r,z)',ier)
            CALL h5dclose_f(temp_gid,ier)
            r2dtemp = SPREAD(zaxis, 1, nr)
            CALL write_var_hdf5(qid_gid,'z2d',nr, nz, ier,DBLVAR=DBLE(r2dtemp*100))
            CALL h5dopen_f(qid_gid, 'z2d', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm',ier)
            CALL write_att_hdf5(temp_gid,'description','Z grid: Z(r,z)',ier)
            CALL h5dclose_f(temp_gid,ier)
            DEALLOCATE(r2dtemp)

            !           MASK
            ALLOCATE(mask(nr,nz))
            mask = 1
            CALL write_var_hdf5(qid_gid,'mask',nr,nz,ier,INTVAR=mask) !PLACEHOLDER, should be "Boolean mask that indicates where the fields are well defined", Dim: [nr,nz]
            
            CALL h5dopen_f(qid_gid, 'mask', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'description','Boolean mask that indicates where the fields are well defined',ier)
            CALL h5dclose_f(temp_gid,ier)
            DEALLOCATE(mask)
            
            !--------------------------------------------------------------
            !           B-FIELD
            !  NOTE On PSI:
            !     In B_STS B_STS_eval_rho defines psi as
            !         rho[0] = sqrt( (psi - Bdata->psi0) / delta );
            !     So psi is TOROIDAL FLUX.
            !--------------------------------------------------------------


            ALLOCATE(rtemp(nr,nz, nphi))
            rtemp = reshape(B_R(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/)) !switch phi and z            
            CALL write_var_hdf5(qid_gid,'br',nr,nz,nphi,ier,DBLVAR=rtemp) 
            CALL h5dopen_f(qid_gid, 'br', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','T',ier)
            CALL write_att_hdf5(temp_gid,'description','Magnetic field in the r-direction: Br(r,z,phi)',ier)
            CALL h5dclose_f(temp_gid,ier)

            rtemp = reshape(B_PHI(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/))            
            CALL write_var_hdf5(qid_gid,'bt',nr,nz,nphi,ier,DBLVAR=rtemp)
            CALL h5dopen_f(qid_gid, 'bt', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','T',ier)
            CALL write_att_hdf5(temp_gid,'description','Magnetic field in the theta/torodial-direction: Bt(r,z,phi)',ier)
            CALL h5dclose_f(temp_gid,ier)

            rtemp = reshape(B_Z(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/))            
            CALL write_var_hdf5(qid_gid,'bz',nr,nz,nphi,ier,DBLVAR=rtemp)
            CALL h5dopen_f(qid_gid, 'bz', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','T',ier)
            CALL write_att_hdf5(temp_gid,'description','Magnetic field in the z-direction: Bz(r,z,phi)',ier)
            CALL h5dclose_f(temp_gid,ier)
            DEALLOCATE(rtemp)

                  !--------------------------------------------------------------
               !           E-FIELD - NEEDS CHECKING
               !--------------------------------------------------------------

               ! Values must be equidistant in rho.
            IF (npot < 1) THEN ! Because we can run with E=0
                  ALLOCATE(rtemp(nr,nz,nphi))
                  rtemp = 0.0
                  CALL write_var_hdf5(qid_gid,'er',nr,nz,nphi, ier,DBLVAR=rtemp)
                  CALL write_var_hdf5(qid_gid,'et',nr,nz,nphi, ier,DBLVAR=rtemp)
                  CALL write_var_hdf5(qid_gid,'ez',nr,nz,nphi, ier,DBLVAR=rtemp)
                  DEALLOCATE(rtemp)
            ELSE
                  ALLOCATE(Efield4D(nr,nphi,nz,3))
                  ALLOCATE(r1dtemp(3))
                  r1dtemp = 1
                  DO i = 1, nr-1 !correct?
                        DO j=1, nphi-1
                              DO k=1, nz-1
                                    !hx     = raxis(i+1) - raxis(i)
                                    !hy     = phiaxis(j+1) - phiaxis(j)
                                    !hz     = zaxis(k+1) - zaxis(k)
                                    !hxi    = one / hx
                                    !hyi    = one / hy
                                    !hzi    = one / hz
                                    CALL R8HERM3FCN(ictE,1,1,fvalE,i,j,k,1,1,1,& !evaluate at grid points
                                    hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                                    POT4D(1,1,1,1),nr,nphi,nz)
                                    r1dtemp(1:3) =-fvalE(1,1:3)
                                    Efield4D(i,j,k,1:3) = r1dtemp
                                    !ider=1
                                    !CALL EZspline_gradient3_r8(POT_spl, ider, raxis(i), phiaxis(j), zaxis(k),r1dtemp, ier)
                                    Efield4D(i,j,k,1:3) = -r1dtemp(1:3)
                              END DO
                        END DO
                  END DO
                  !CALL EZspline_gradient3_r8(POT_spl, ider, nr, nphi, nz, raxis, phiaxis, zaxis,Efield4D(1:nr,1:nphi,1:nz,1:3), ier)
                  ALLOCATE(rtemp(nr,nz,nphi))
                  rtemp = reshape(Efield4D(1:nr,1:nphi,1:nz,1), shape(rtemp), order=(/1, 3, 2/)) 
                  CALL write_var_hdf5(qid_gid,'er',nr,nz,nphi, ier,DBLVAR=rtemp)
                  CALL h5dopen_f(qid_gid, 'er', temp_gid, ier)
                  CALL write_att_hdf5(temp_gid,'units','V/m',ier)
                  CALL write_att_hdf5(temp_gid,'description','Electric field in the r-direction: Er(r,z,phi)',ier)
                  CALL h5dclose_f(temp_gid,ier)
                  rtemp = reshape(Efield4D(1:nr,1:nphi,1:nz,2), shape(rtemp), order=(/1, 3, 2/)) 
                  CALL write_var_hdf5(qid_gid,'et',nr,nz,nphi, ier,DBLVAR=rtemp)
                  CALL h5dopen_f(qid_gid, 'et', temp_gid, ier)
                  CALL write_att_hdf5(temp_gid,'units','V/m',ier)
                  CALL write_att_hdf5(temp_gid,'description','Electric field in the toroidal phi-direction: Et(r,z,phi)',ier)
                  CALL h5dclose_f(temp_gid,ier)
                  rtemp = reshape(Efield4D(1:nr,1:nphi,1:nz,3), shape(rtemp), order=(/1, 3, 2/)) 
                  CALL write_var_hdf5(qid_gid,'ez',nr,nz,nphi, ier,DBLVAR=rtemp)
                  CALL h5dopen_f(qid_gid, 'ez', temp_gid, ier)
                  CALL write_att_hdf5(temp_gid,'units','V/m',ier)
                  CALL write_att_hdf5(temp_gid,'description','Electric field in the z-direction: Ez(r,z,phi)',ier)
                  CALL h5dclose_f(temp_gid,ier)
                  DEALLOCATE(r1dtemp)
                  DEALLOCATE(rtemp)
                  DEALLOCATE(Efield4D)
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
            CALL write_var_hdf5(qid_gid,'nr',ier,INTVAR=nr)
            CALL h5dopen_f(qid_gid, 'nr', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','-',ier)
            CALL write_att_hdf5(temp_gid,'description','Number of R values',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'nphi',ier,INTVAR=nphi)
            CALL h5dopen_f(qid_gid, 'nphi', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','-',ier)
            CALL write_att_hdf5(temp_gid,'description','Number of Phi values',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'nz',ier,INTVAR=nz)
            CALL h5dopen_f(qid_gid, 'nz', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','-',ier)
            CALL write_att_hdf5(temp_gid,'description','Number of Z values',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'time',ier,DBLVAR=DBLE(0)) !!!!Assumes steady state/dummy
            CALL h5dopen_f(qid_gid, 'time', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','s',ier)
            CALL write_att_hdf5(temp_gid,'description','Distribution time',ier)
            CALL h5dclose_f(temp_gid,ier)

            CALL write_var_hdf5(qid_gid,'r',nr, ier,DBLVAR=DBLE(raxis*100)) !convert from m to cm
            CALL h5dopen_f(qid_gid, 'r', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm',ier)
            CALL write_att_hdf5(temp_gid,'description','Radius',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'phi',nphi, ier,DBLVAR=phiaxis)
            CALL h5dopen_f(qid_gid, 'phi', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','rad',ier)
            CALL write_att_hdf5(temp_gid,'description','Toroidal angle',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'z',nz, ier,DBLVAR=DBLE(zaxis*100)) !convert from m to cm
            CALL h5dopen_f(qid_gid, 'z', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm',ier)
            CALL write_att_hdf5(temp_gid,'description','Z',ier)
            CALL h5dclose_f(temp_gid,ier)


            ALLOCATE(r2dtemp(nr,nz))
            r2dtemp = SPREAD(raxis, 2, nz)
            CALL write_var_hdf5(qid_gid,'r2d',nr, nz, ier,DBLVAR=DBLE(r2dtemp*100))
            CALL h5dopen_f(qid_gid, 'r2d', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm',ier)
            CALL write_att_hdf5(temp_gid,'description','Radius grid: R(r,z)',ier)
            CALL h5dclose_f(temp_gid,ier)
            r2dtemp = SPREAD(zaxis, 1, nr)
            CALL write_var_hdf5(qid_gid,'z2d',nr, nz, ier,DBLVAR=DBLE(r2dtemp*100))
            CALL h5dopen_f(qid_gid, 'z2d', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm',ier)
            CALL write_att_hdf5(temp_gid,'description','Z grid: Z(r,z)',ier)
            CALL h5dclose_f(temp_gid,ier)
            DEALLOCATE(r2dtemp)

            !           MASK
            ALLOCATE(mask(nr,nz))
            mask = 1
            CALL write_var_hdf5(qid_gid,'mask',nr,nz,ier,INTVAR=mask) !PLACEHOLDER, should be "Boolean mask that indicates where the fields are well defined", Dim: [nr,nz]            
            CALL h5dopen_f(qid_gid, 'mask', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'description','Boolean mask that indicates where the fields are well defined',ier)
            CALL h5dclose_f(temp_gid,ier)
            DEALLOCATE(mask)

            !           PLASMA ROTATION/FLOW
            ALLOCATE(rtemp(nr,nz,nphi))
            rtemp = 0.0
            CALL write_var_hdf5(qid_gid,'vr',nr,nz, nphi,ier,DBLVAR=rtemp)            
            CALL h5dopen_f(qid_gid, 'vr', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'description','Bulk plasma flow in the r-direction: Vr(r,z,phi)',ier)
            CALL write_att_hdf5(temp_gid,'units','cm/s',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'vt',nr,nz, nphi,ier,DBLVAR=rtemp)
            CALL h5dopen_f(qid_gid, 'vt', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'description','Bulk plasma flow in the toroidal phi-direction: Vphi(r,z,phi)',ier)
            CALL write_att_hdf5(temp_gid,'units','cm/s',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'vz',nr,nz, nphi,ier,DBLVAR=rtemp)
            CALL h5dopen_f(qid_gid, 'vz', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'description','Bulk plasma flow in the z-direction: Vz(r,z,phi)',ier)
            CALL write_att_hdf5(temp_gid,'units','cm/s',ier)
            CALL h5dclose_f(temp_gid,ier)

            rtemp = reshape(NE(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/)) 
            CALL write_var_hdf5(qid_gid,'dene',nr,nz,nphi, ier,DBLVAR=DBLE(rtemp/1000000)) !convert from m^-3 to cm^-3
            CALL h5dopen_f(qid_gid, 'dene', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm^-3',ier)
            CALL write_att_hdf5(temp_gid,'description','Electron Number Density: Dene(r,z, phi)',ier)
            CALL h5dclose_f(temp_gid,ier)
            rtemp = reshape(TE(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/)) 
            CALL write_var_hdf5(qid_gid,'te',nr,nz,nphi, ier,DBLVAR=DBLE(rtemp/1000))
            CALL h5dopen_f(qid_gid, 'te', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','kev',ier)
            CALL write_att_hdf5(temp_gid,'description','Electron Temperature: Ti(r,z,phi)',ier)
            CALL h5dclose_f(temp_gid,ier)
            rtemp = reshape(TI(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/)) 
            CALL write_var_hdf5(qid_gid,'ti',nr,nz,nphi, ier,DBLVAR=DBLE(rtemp/1000))
            CALL h5dopen_f(qid_gid, 'ti', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','keV',ier)
            CALL write_att_hdf5(temp_gid,'description','Ion Temperature: Ti(r,z,phi)',ier)
            CALL h5dclose_f(temp_gid,ier)
            rtemp = reshape(ZEFF_ARR(1:nr,1:nphi,1:nz), shape(rtemp), order=(/1, 3, 2/)) 
            CALL write_var_hdf5(qid_gid,'zeff',nr,nz,nphi, ier,DBLVAR=DBLE(rtemp))
            CALL h5dopen_f(qid_gid, 'zeff', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','-',ier)
            CALL write_att_hdf5(temp_gid,'description','Effective Nuclear Charge: Zeff(r,z,phi)',ier)
            CALL h5dclose_f(temp_gid,ier)
            DEALLOCATE(rtemp)


            CALL h5gclose_f(qid_gid, ier)
                        ! Close file
            CALL close_hdf5(fid,ier)
            IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'ascot5_'//TRIM(id_string)//'_equilibrium.h5',ier)
            END IF

      CASE('DISTRIBUTION_GC_F')
            CALL open_hdf5('fidasim_'//TRIM(id_string)//'_distribution.h5',fid,ier,LCREATE=.false.) 
            IF (ASSOCIATED(dist5d_prof)) THEN
                  CALL write_var_hdf5(fid,'dist_prof',nbeams,ns_prof1,ns_prof2,ns_prof3,ns_prof4,ns_prof5,ier,DBLVAR=dist5d_prof,&
                                      ATT='Distribution Function [part/(m^3/s^3)] no physical volume (nbeam,nrho,npol,ntor,nvll,nvperp)',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'dist_prof',ier)
            END IF
      CASE('DISTRIBUTION_GC_MC')
      CASE('DISTRIBUTION_FO')

      END SELECT

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_write_fidasim
