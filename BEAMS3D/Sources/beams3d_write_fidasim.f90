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
USE beams3d_lines, ONLY: ns_prof1, ns_prof2, ns_prof3, ns_prof4, &
                                    ns_prof4, ns_prof5, dist5d_prof, &
                                    partvmax, dist5D_fida, &
                                    h2_prof, h3_prof, h4_prof, h5_prof, my_end
      USE beams3d_grid, ONLY: nr, nphi, nz, B_R, B_PHI, B_Z, raxis, &
                                 zaxis, phiaxis, S_ARR, U_ARR, POT_ARR, &
                                 ZEFF_ARR, TE, TI, NE, req_axis, zeq_axis, npot, &
                                 POT_SPL, ezspline_interp, phiedge_eq, TE_spl, &
                                 NE_spl, TI_spl, ZEFF_spl, nne, nte, nti, nzeff, &
                                 plasma_mass, reff_eq, therm_factor, &
                                 X_BEAMLET, Y_BEAMLET, Z_BEAMLET, &
                                 NX_BEAMLET, NY_BEAMLET, NZ_BEAMLET, &
                                 POT4D, NE4D, TE4D, TI4D, ZEFF4D, &
                                 hr, hp, hz, hri, hpi, hzi, S4D, U4D, &
                                 rmin, rmax, zmin, zmax, phimin, phimax, raxis, zaxis, phiaxis, &
                                 rmin_fida, rmax_fida, zmin_fida, zmax_fida, phimin_fida, phimax_fida, &
                                 raxis_fida, zaxis_fida, phiaxis_fida, nr_fida, nphi_fida, nz_fida
      USE beams3d_runtime, ONLY: id_string, npoinc, nbeams, beam, t_end, lverb, &
                                    lvmec, lpies, lspec, lcoil, lmgrid, lbeam, lplasma_only, &
                                    lvessel, lvac, lbeam_simple, handle_err, nparticles_start, &
                                    HDF5_OPEN_ERR,HDF5_WRITE_ERR,&
                                    HDF5_CLOSE_ERR, BEAMS3D_VERSION, weight, e_beams, p_beams,&
                                    charge, Zatom, mass, ldepo, v_neut, &
                                    lcollision, pi, pi2, t_end_in, nprocs_beams, &
                                    div_beams, mass_beams, Zatom_beams, dex_beams, &
                                    lascotfl, R_beams, PHI_beams, Z_beams
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
      REAL*8 :: fvalE(1,3)
      !REAL*8, POINTER :: hr_fida(:), hp_fida(:), hz_fida(:)
      !REAL*8, POINTER :: hri_fida(:), hpi_fida(:), hzi_fida(:)
      REAL(rprec), DIMENSION(4) :: rt,zt,pt
      REAL*8 :: xparam, yparam, zparam !, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(1)
      INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
      DOUBLE PRECISION         :: x0,y0,z0
      REAL(rprec) :: jac, v_parr, v_perp, pitch
      DOUBLE PRECISION :: rho_temp, s_temp, dbl_temp, gammarel, v_total!, hx, hy, hz, hxi, hyi, hzi, xparam, yparam, zparam
      DOUBLE PRECISION, ALLOCATABLE :: rtemp(:,:,:), r1dtemp(:), r2dtemp(:,:), Efield4D(:,:,:,:)
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mask
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
      IF (myworkid == master) THEN
      SELECT CASE (TRIM(write_type))
         CASE('INIT')
            rmin_fida = rmin
            zmin_fida = zmin
            phimin_fida = phimin
            rmax_fida = rmax
            zmax_fida = zmax
            phimax_fida = phimax/5.
            nr_fida = nr
            nphi_fida = nphi
            nz_fida = nz


            ALLOCATE(raxis_fida(nr_fida))
            ALLOCATE(zaxis_fida(nz_fida))
            ALLOCATE(phiaxis_fida(nphi_fida))
            FORALL(i = 1:nr_fida) raxis_fida(i) = (i-1)*(rmax_fida-rmin_fida)/(nr_fida-1) + rmin_fida
            FORALL(i = 1:nz_fida) zaxis_fida(i) = (i-1)*(zmax_fida-zmin_fida)/(nz_fida-1) + zmin_fida
            FORALL(i = 1:nphi_fida) phiaxis_fida(i) = (i-1)*(phimax_fida-phimin_fida)/(nphi_fida-1) + phimin_fida
            ! Setup grid helpers
            ! Note: All helpers are defined in terms of differences on half grid
            !       so values are indexed from 1 to n-1.  Which we store at n
            !        i = MIN(MAX(COUNT(raxis < r_temp),1),nr-1)
            !        hr(i) = raxis(i+1) - raxis(i)
            !        hri    = one / hr
            !FORALL(i = 1:nr_fida-1) hr_fida(i) = raxis_fida(i+1) - raxis_fida(i)
            !FORALL(i = 1:nz_fida-1) hz_fida(i) = zaxis_fida(i+1) - zaxis_fida(i)
            !FORALL(i = 1:nphi_fida-1) hp_fida(i) = phiaxis_fida(i+1) - phiaxis_fida(i)
            !hri_fida = one / hr_fida
            !hpi_fida = one / hp_fida
            !hzi_fida = one / hz_fida

            nphi1 = nphi-1 ! confirmed with Poincare plot
            CALL open_hdf5('fidasim_'//TRIM(id_string)//'_distribution.h5',fid,ier,LCREATE=.true.)
            CALL h5gopen_f(fid,'/', qid_gid, ier)   
            CALL write_att_hdf5(qid_gid,'data_source','Data initialized from BEAMS3D',ier)
            CALL write_var_hdf5(qid_gid,'type',ier,INTVAR=2)
            CALL h5dopen_f(fid, 'type', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'description','Distribution type: 1="Guiding Center Density Function", 2="Guiding Center Monte Carlo", 3="Full Orbit Monte Carlo"',ier)
            CALL h5dclose_f(temp_gid,ier)
            


            ! FIDASIM GRID
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
            CALL write_var_hdf5(qid_gid,'time',ier,DBLVAR=DBLE(0)) !!!!Assumes steady state/dummy
            CALL h5dopen_f(fid, 'time', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','s',ier)
            CALL write_att_hdf5(temp_gid,'description','Distribution time',ier)
            CALL h5dclose_f(temp_gid,ier)

            CALL write_var_hdf5(qid_gid,'r',nr, ier,DBLVAR=DBLE(raxis_fida*100)) !convert from m to cm
            CALL h5dopen_f(fid, 'r', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','cm',ier)
            CALL write_att_hdf5(temp_gid,'description','Radius',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'phi',nphi, ier,DBLVAR=DBLE(phiaxis_fida)) !Hard-coded toroidal slice for now
            CALL h5dopen_f(fid, 'phi', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'units','rad',ier)
            CALL write_att_hdf5(temp_gid,'description','Toroidal angle',ier)
            CALL h5dclose_f(temp_gid,ier)
            CALL write_var_hdf5(qid_gid,'z',nz, ier,DBLVAR=DBLE(zaxis_fida*100)) !convert from m to cm
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
            

      CASE('DISTRIBUTION_GC_F')
      ! Do volume normalization
        CALL beams3d_distnorm !TODO: check if this has already been done
        !Apply jacobian for transformation from v_parr/v_perp to E,p+
          DO b = 1, nbeams
              DO i = 1, ns_prof4
                    DO j = 1, ns_prof5
                          v_parr = REAL(i) / ns_prof4 * partvmax ! full grid or half grid, which alignment?
                          v_perp = REAL(j) / ns_prof5 * partvmax
                          pitch = MAX(MIN(v_parr / SQRT(v_parr * v_parr + v_perp * v_perp),1.),0.)
                          jac = 1 / (mass_beams(b) * SQRT(1-pitch * pitch))
                          jac = jac * 2 / 1000000 * e_charge !convert to TRANSP convention? and 1/cm^3/eV
                          dist5d_prof(b,:,:,:,i,j) = dist5d_prof(b,:,:,:,i,j) * jac !probably should inicate somewhere that this has been done, or put it in another variable
                    END DO
              END DO
        END DO
  

        ALLOCATE(dist5d_fida(nbeams,nr_fida,nz_fida,nphi_fida,ns_prof4,ns_prof5))
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
                        CALL R8HERM3FCN(ict,1,1,fval,i3,j3,k3,xparam,yparam,zparam,&
                                    hr(i3),hri(i3),hp(j3),hpi(j3),hz(k3),hzi(k3),&
                                    S4D(1,1,1,1),nr,nphi,nz)
                        y0 = fval(1)
                        CALL R8HERM3FCN(ict,1,1,fval,i3,j3,k3,xparam,yparam,zparam,&
                                    hr(i3),hri(i3),hp(j3),hpi(j3),hz(k3),hzi(k3),&
                                    U4D(1,1,1,1),nr,nphi,nz)
                        z0 = fval(1)
                              
                                    ! Calc dist func bins
                        x0    = phiaxis_fida(j)
                        IF (x0 < 0) x0 = x0 + pi2
                        !vperp = SQRT(2*moment*fval(1)/mymass)
                        l = MAX(MIN(CEILING(SQRT(y0)*ns_prof1     ), ns_prof1), 1) ! Rho Bin
                        m = MAX(MIN(CEILING( z0*h2_prof           ), ns_prof2), 1) ! U Bin
                        n = MAX(MIN(CEILING( x0*h3_prof           ), ns_prof3), 1) ! V Bin
                        !d4 = MAX(MIN(1+nsh_prof4+FLOOR(h4_prof*q(4)), ns_prof4), 1) ! vll
                        !d5 = MAX(MIN(CEILING(vperp*h5_prof         ), ns_prof5), 1) ! Vperp
                        
                        dist5d_fida(b,i,k,j,:,:) = dist5d_prof(b,l,m,n,:,:) !output in r-z-phi
                        
  
                        END DO
                  END DO
            END DO
      END DO
            !dist5d_fida = reshape(dist5d_fida, [nbeams,nr_fida,nz_fida,nphi_fida,ns_prof4,ns_prof5], order=(/1, 2, 4, 3, 5, 6/)) 
            CALL open_hdf5('fidasim_'//TRIM(id_string)//'_distribution.h5',fid,ier,LCREATE=.false.) 
            IF (ASSOCIATED(dist5d_prof)) THEN
                  CALL write_var_hdf5(fid,'f',nbeams,nr_fida,nz_fida,nphi_fida,ns_prof4,ns_prof5,ier,DBLVAR=dist5d_fida,&
                                      ATT='Distribution Function [part/(cm^3 eV)] (nr_fida,nz_fida,nphi_fida,ns_prof4,ns_prof5)',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'dist_fida',ier)
            END IF

      CASE('DISTRIBUTION_GC_MC')
      CASE('DISTRIBUTION_FO')

      END SELECT
      END IF !myworkid==master
      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_write_fidasim
