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
                                 NE_spl_s, TI_spl_s, ZEFF_spl_s
      USE beams3d_runtime, ONLY: id_string, npoinc, nbeams, beam, t_end, lverb, lflux, &
                                    lvmec, lpies, lspec, lcoil, lmgrid, lbeam, &
                                    lvessel, lvac, lbeam_simple, handle_err, nparticles_start, &
                                    HDF5_OPEN_ERR,HDF5_WRITE_ERR,&
                                    HDF5_CLOSE_ERR, BEAMS3D_VERSION, weight, e_beams, p_beams,&
                                    charge, Zatom, mass, ldepo, v_neut,lcollision, pi, pi2
      USE safe_open_mod, ONLY: safe_open
      USE wall_mod, ONLY: nface,nvertex,face,vertex,ihit_array, wall_free
      USE beams3d_write_par
      USE mpi_params
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
      INTEGER :: ier, iunit, i, d1, d2, d3
      INTEGER(HID_T) :: options_gid, bfield_gid, efield_gid, plasma_gid, &
                        neutral_gid, wall_gid, marker_gid, qid_gid
      CHARACTER(LEN=10) ::  qid_str
      DOUBLE PRECISION, ALLOCATABLE :: rtemp(:,:,:), r1dtemp(:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      qid_str='0000000000'
      SELECT CASE (TRIM(write_type))
         CASE('INIT')
            IF (myworkid == master) THEN
               CALL open_hdf5('beams3d_ascot5_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.true.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'beams3d_ascot5_'//TRIM(id_string)//'.h5',ier)
               ! Write the ID
               !CALL write_var_hdf5(fid,'/options/',nz,ier,DBLVAR=zaxis,ATT='Vertical Axis [rad]',ATT_NAME='description')

               !--------------------------------------------------------------
               !           OPTIONS
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'options', options_gid, ier)
               CALL write_att_hdf5(options_gid,'active',qid_str,ier)
               CALL h5gcreate_f(options_gid,'opt_'//qid_str, qid_gid, ier)
               CALL write_var_hdf5(qid_gid,'SIM_MODE',ier,DBLVAR=DBLE(1))
               CALL write_var_hdf5(qid_gid,'ENABLE_ADAPTIVE',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'RECORD_MODE',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'FIXEDSTEP_USE_USERDEFINED',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'FIXEDSTEP_USERDEFINED',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'FIXEDSTEP_GYRODEFINED',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ADAPTIVE_TOL_ORBIT',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ADAPTIVE_TOL_CCOL',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ADAPTIVE_MAX_DRHO',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ADAPTIVE_MAX_DPHI',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENABLE_ORBIT_FOLLOWING',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENABLE_COULOMB_COLLISIONS',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'DISABLE_FIRSTORDER_GCTRANS',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'DISABLE_ENERGY_CCOLL',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'DISABLE_PITCH_CCOLL',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'DISABLE_GCDIFF_CCOLL',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_SIMTIMELIM',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_CPUTIMELIM',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_RHOLIM',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_ENERGYLIM',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_WALLHIT',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MAXORBS',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MAX_SIMTIME',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MAX_CPUTIME',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MAX_RHO',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MIN_RHO',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MIN_ENERGY',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MIN_THERMAL',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MAX_POLOIDALORBS',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENDCOND_MAX_TOROIDALORBS',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENABLE_DIST_5D',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENABLE_DIST_6D',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENABLE_DIST_RHO5D',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENABLE_DIST_RHO6D',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'ENABLE_ORBITWRITE',ier,DBLVAR=DBLE(0))
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(options_gid, ier)

               !--------------------------------------------------------------
               !           B-FIELD
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'bfield', bfield_gid, ier)
               CALL write_att_hdf5(bfield_gid,'active',qid_str,ier)
               CALL h5gcreate_f(bfield_gid,'B_STS_'//qid_str, qid_gid, ier)
               CALL write_var_hdf5(qid_gid,'b_nr',ier,INTVAR=nr)
               CALL write_var_hdf5(qid_gid,'b_nphi',ier,INTVAR=nphi)
               CALL write_var_hdf5(qid_gid,'b_nz',ier,INTVAR=nz)
               CALL write_var_hdf5(qid_gid,'b_rmin',ier,DBLVAR=raxis(1))
               CALL write_var_hdf5(qid_gid,'b_rmax',ier,DBLVAR=raxis(nr))
               CALL write_var_hdf5(qid_gid,'b_zmin',ier,DBLVAR=zaxis(1))
               CALL write_var_hdf5(qid_gid,'b_zmax',ier,DBLVAR=zaxis(nz))
               CALL write_var_hdf5(qid_gid,'b_phimin',ier,DBLVAR=phiaxis(1)*180/pi)
               CALL write_var_hdf5(qid_gid,'b_phimax',ier,DBLVAR=phiaxis(nphi)*180/pi)
               CALL write_var_hdf5(qid_gid,'psi_nr',ier,INTVAR=nr)
               CALL write_var_hdf5(qid_gid,'psi_nphi',ier,INTVAR=nphi)
               CALL write_var_hdf5(qid_gid,'psi_nz',ier,INTVAR=nz)
               CALL write_var_hdf5(qid_gid,'psi_rmin',ier,DBLVAR=raxis(1))
               CALL write_var_hdf5(qid_gid,'psi_rmax',ier,DBLVAR=raxis(nr))
               CALL write_var_hdf5(qid_gid,'psi_zmin',ier,DBLVAR=zaxis(1))
               CALL write_var_hdf5(qid_gid,'psi_zmax',ier,DBLVAR=zaxis(nz))
               CALL write_var_hdf5(qid_gid,'psi_phimin',ier,DBLVAR=phiaxis(1)*180/pi)
               CALL write_var_hdf5(qid_gid,'psi_phimax',ier,DBLVAR=phiaxis(nphi)*180/pi)
               CALL write_var_hdf5(qid_gid,'axis_nphi',ier,INTVAR=nphi)
               CALL write_var_hdf5(qid_gid,'axis_phimin',ier,DBLVAR=phiaxis(1)*180/pi)
               CALL write_var_hdf5(qid_gid,'axis_phimax',ier,DBLVAR=phiaxis(nphi)*180/pi)
               CALL write_var_hdf5(qid_gid,'toroidalPeriods',ier,INTVAR=FLOOR(pi2/phiaxis(nphi)))
               CALL write_var_hdf5(qid_gid,'axisr',nphi,ier,DBLVAR=req_axis)
               CALL write_var_hdf5(qid_gid,'axisz',nphi,ier,DBLVAR=zeq_axis)
               CALL write_var_hdf5(qid_gid,'psi0',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'psi1',ier,DBLVAR=DBLE(phiedge_eq))
               ALLOCATE(rtemp(nr,nphi,nz))
               rtemp = RESHAPE(B_R,(/nr,nphi,nz/),ORDER=(/1,2,3/))
               CALL write_var_hdf5(qid_gid,'br',nr,nphi,nz,ier,DBLVAR=rtemp)
               rtemp = RESHAPE(B_PHI,(/nr,nphi,nz/),ORDER=(/1,2,3/))
               CALL write_var_hdf5(qid_gid,'bphi',nr,nphi,nz,ier,DBLVAR=rtemp)
               rtemp = RESHAPE(B_Z,(/nr,nphi,nz/),ORDER=(/1,2,3/))
               CALL write_var_hdf5(qid_gid,'bz',nr,nphi,nz,ier,DBLVAR=rtemp)
               rtemp = RESHAPE(S_ARR,(/nr,nphi,nz/),ORDER=(/1,2,3/))
               CALL write_var_hdf5(qid_gid,'psi',nr,nphi,nz,ier,DBLVAR=rtemp)
               DEALLOCATE(rtemp)
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(bfield_gid, ier)

               !--------------------------------------------------------------
               !           E-FIELD
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'efield', efield_gid, ier)
               CALL write_att_hdf5(efield_gid,'active',qid_str,ier)
               CALL h5gcreate_f(efield_gid,'E_1DS_'//qid_str, qid_gid, ier)
               IF (npot < 1) THEN ! Because we can run with E=0
                  ALLOCATE(r1dtemp(5))
                  r1dtemp = 0
                  CALL write_var_hdf5(qid_gid,'nrho',ier,INTVAR=5)
                  CALL write_var_hdf5(qid_gid,'dvdrho',5,ier,DBLVAR=r1dtemp)
                  DEALLOCATE(r1dtemp)
               ELSE
                  ! This is glitchy since we don't require POT_AUX_S to be equidistant
                  ! should really spline to new grid
                  ALLOCATE(r1dtemp(nr))
                  DO i = 1, nr
                    CALL EZspline_interp(POT_spl_s,DBLE(i-1)/DBLE(nr-1),r1dtemp(i),ier)
                  END DO
                  CALL write_var_hdf5(qid_gid,'nrho',ier,INTVAR=nr)
                  CALL write_var_hdf5(qid_gid,'dvdrho',nr,ier,DBLVAR=r1dtemp)
                  DEALLOCATE(r1dtemp)
               END IF
               CALL write_var_hdf5(qid_gid,'rhomin',ier,DBLVAR=DBLE(0))
               CALL write_var_hdf5(qid_gid,'rhomax',ier,DBLVAR=DBLE(1))
               CALL write_var_hdf5(qid_gid,'reff',ier,DBLVAR=DBLE(phiedge_eq))
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(efield_gid, ier)


               !--------------------------------------------------------------
               !           PLASMA
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'plasma', plasma_gid, ier)
               CALL write_att_hdf5(plasma_gid,'active',qid_str,ier)
               CALL h5gcreate_f(plasma_gid,'plasma_1D_'//qid_str, qid_gid, ier)
               CALL write_var_hdf5(qid_gid,'nion',ier,INTVAR=1)
               CALL write_var_hdf5(qid_gid,'nrho',ier,INTVAR=nr)
               CALL write_var_hdf5(qid_gid,'znum',ier,INTVAR=1)
               CALL write_var_hdf5(qid_gid,'anum',ier,INTVAR=1)
               CALL write_var_hdf5(qid_gid,'charge',ier,INTVAR=1)
               CALL write_var_hdf5(qid_gid,'mass',ier,INTVAR=1)
               ALLOCATE(rtemp(nr,5,1))
               DO i = 1, nr
                 rtemp(i,1,1)=DBLE(i-1)/DBLE(nr-1)
                 CALL EZspline_interp( TE_spl_s,   rtemp(i,1,1), rtemp(i,2,1), ier)
                 CALL EZspline_interp( NE_spl_s,   rtemp(i,1,1), rtemp(i,3,1), ier)
                 CALL EZspline_interp( TI_spl_s,   rtemp(i,1,1), rtemp(i,4,1), ier)
                 CALL EZspline_interp( ZEFF_spl_s, rtemp(i,1,1), rtemp(i,5,1), ier)
                 rtemp(i,5,1)=rtemp(i,3,1)/rtemp(i,5,1)
               END DO
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
               CALL write_att_hdf5(neutral_gid,'active',qid_str,ier)
               CALL h5gcreate_f(neutral_gid,'N0_3D_'//qid_str, qid_gid, ier)
               CALL write_var_hdf5(qid_gid,'nr',ier,INTVAR=nr)
               CALL write_var_hdf5(qid_gid,'nphi',ier,INTVAR=nphi)
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
               ALLOCATE(rtemp(nr,nphi,nz))
               rtemp = RESHAPE(S_ARR,(/nr,nphi,nz/),ORDER=(/1,2,3/))
               rtemp = 0
               CALL write_var_hdf5(qid_gid,'density',nr,nphi,nz,ier,DBLVAR=rtemp)
               CALL write_var_hdf5(qid_gid,'temperature',nr,nphi,nz,ier,DBLVAR=rtemp)
               DEALLOCATE(rtemp)
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(neutral_gid, ier)


               !--------------------------------------------------------------
               !           WALL
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'wall', wall_gid, ier)
               CALL write_att_hdf5(wall_gid,'active',qid_str,ier)
               CALL h5gcreate_f(wall_gid,'wall_3D_'//qid_str, qid_gid, ier)
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
               IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'beams3d_ascot5_'//TRIM(id_string)//'.h5',ier)
            END IF
         CASE('MARKER')
            IF (myworkid == master) THEN
               CALL open_hdf5('beams3d_ascot5_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.false.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'beams3d_ascot5_'//TRIM(id_string)//'.h5',ier)
               !--------------------------------------------------------------
               !           MARKER
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'marker', marker_gid, ier)
               CALL write_att_hdf5(marker_gid,'active',qid_str,ier)
               CALL h5gcreate_f(marker_gid,'gc_'//qid_str, qid_gid, ier)
               CALL write_var_hdf5(qid_gid,'n',ier,INTVAR=nparticles)
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(marker_gid, ier)
               CALL close_hdf5(fid,ier)
               IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'beams3d_ascot5_'//TRIM(id_string)//'.h5',ier)
            END IF
            ! For now we assume we will only run this in deposition mode
            ! so that means we look at index (1,i) since
            !    0: starting point
            !    1: Ionization point (or wall)
            !    2: Gyrocenter
            d1 = LBOUND(R_lines,DIM=2)
            d2 = UBOUND(R_lines,DIM=2)
            ALLOCATE(rtemp(d1:d2,13,1))
            DO i = d1, d2
               rtemp(i,1,1) = R_lines(2,i)
               rtemp(i,2,1) = PHI_lines(2,i)
               rtemp(i,3,1) = Z_lines(2,i)
               rtemp(i,4,1) = 0.5*mass(i)*vll_lines(2,i)*vll_lines(2,i)
               rtemp(i,5,1) = 0.0 ! pitch
               rtemp(i,6,1) = 0.0 ! zeta
               rtemp(i,7,1) = mass(i) ! mass
               rtemp(i,8,1) = charge(i)
               rtemp(i,9,1) = mass(i)/1.6726231000E-27 ! Anum
               rtemp(i,10,1) = Zatom(i)
               rtemp(i,11,1) = 1.0 ! weight
               rtemp(i,12,1) = 0.0 ! time
               rtemp(i,13,1) = i
            END DO
            CALL beams3d_write1d_parhdf5( 1, nparticles, d1, d2, '/marker/gc_'//qid_str//'/r',      DBLVAR=rtemp(:,1,1))
            CALL beams3d_write1d_parhdf5( 1, nparticles, d1, d2, '/marker/gc_'//qid_str//'/phi',    DBLVAR=rtemp(:,2,1))
            CALL beams3d_write1d_parhdf5( 1, nparticles, d1, d2, '/marker/gc_'//qid_str//'/z',      DBLVAR=rtemp(:,3,1))
            CALL beams3d_write1d_parhdf5( 1, nparticles, d1, d2, '/marker/gc_'//qid_str//'/energy', DBLVAR=rtemp(:,4,1))
            CALL beams3d_write1d_parhdf5( 1, nparticles, d1, d2, '/marker/gc_'//qid_str//'/pitch',  DBLVAR=rtemp(:,5,1))
            CALL beams3d_write1d_parhdf5( 1, nparticles, d1, d2, '/marker/gc_'//qid_str//'/zeta',   DBLVAR=rtemp(:,6,1))
            CALL beams3d_write1d_parhdf5( 1, nparticles, d1, d2, '/marker/gc_'//qid_str//'/mass',   DBLVAR=rtemp(:,7,1))
            CALL beams3d_write1d_parhdf5( 1, nparticles, d1, d2, '/marker/gc_'//qid_str//'/charge', DBLVAR=rtemp(:,8,1))
            CALL beams3d_write1d_parhdf5( 1, nparticles, d1, d2, '/marker/gc_'//qid_str//'/anum',   INTVAR=INT(rtemp(:,9,1)))
            CALL beams3d_write1d_parhdf5( 1, nparticles, d1, d2, '/marker/gc_'//qid_str//'/znum',   INTVAR=INT(rtemp(:,10,1)))
            CALL beams3d_write1d_parhdf5( 1, nparticles, d1, d2, '/marker/gc_'//qid_str//'/weight', DBLVAR=rtemp(:,11,1))
            CALL beams3d_write1d_parhdf5( 1, nparticles, d1, d2, '/marker/gc_'//qid_str//'/time',   DBLVAR=rtemp(:,12,1))
            CALL beams3d_write1d_parhdf5( 1, nparticles, d1, d2, '/marker/gc_'//qid_str//'/id',     INTVAR=INT(rtemp(:,13,1)))
            DEALLOCATE(rtemp)
      END SELECT

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_write_ascoth5
