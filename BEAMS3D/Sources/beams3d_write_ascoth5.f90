!-----------------------------------------------------------------------
!     Module:        beams3d_write_ascoth5
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          08/21/2019
!     Description:   This subroutine outputs simulation data for the
!                    ASCOT5 code.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_write_ascoth5
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
                                 ZEFF_ARR, TE, TI, NE
      USE beams3d_runtime, ONLY: id_string, npoinc, nbeams, beam, t_end, lverb, lflux, &
                                    lvmec, lpies, lspec, lcoil, lmgrid, lbeam, &
                                    lvessel, lvac, lbeam_simple, handle_err, nparticles_start, &
                                    HDF5_OPEN_ERR,HDF5_WRITE_ERR,&
                                    HDF5_CLOSE_ERR, BEAMS3D_VERSION, weight, e_beams, p_beams,&
                                    charge, Zatom, mass, ldepo, v_neut,lcollision, pi, pi2
      USE safe_open_mod, ONLY: safe_open
      USE wall_mod, ONLY: nface,nvertex,face,vertex,ihit_array, wall_free
      USE mpi_params
!-----------------------------------------------------------------------
!     Input Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          iunit        File ID
!-----------------------------------------------------------------------
      INTEGER :: ier, iunit
      INTEGER(HID_T) :: options_gid, bfield_gid, efield_gid, plasma_gid, &
                        neutral_gid, wall_gid, marker_gid, qid_gid
      INTEGER(HID_T) :: atype_id
      INTEGER(HID_T) :: temp_sid
      INTEGER(HID_T) :: temp_aid
      INTEGER(HID_T) :: temp_did
      INTEGER(HSIZE_T), DIMENSION(1) :: h5_1d_help = (/1/)
      INTEGER(SIZE_T) :: help_int
      CHARACTER(LEN=10) ::  qid_str
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (myworkid == master) THEN
         qid_str='0000000000'
         CALL open_hdf5('beams3d_ascot5_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.true.)
         IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'beams3d_ascot5_'//TRIM(id_string)//'.h5',ier)
         ! Write the ID
         !CALL write_var_hdf5(fid,'/options/',nz,ier,DBLVAR=zaxis,ATT='Vertical Axis [rad]',ATT_NAME='description')

         ! Create the options Groups
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
         CALL h5gclose_f(qid_gid, ier)
         CALL h5gclose_f(bfield_gid, ier)

         CALL h5gcreate_f(fid,'efield', efield_gid, ier)
         CALL write_att_hdf5(efield_gid,'active',qid_str,ier)
         CALL h5gcreate_f(efield_gid,'opt_'//qid_str, qid_gid, ier)
         CALL h5gclose_f(qid_gid, ier)
         CALL h5gclose_f(efield_gid, ier)

         CALL h5gcreate_f(fid,'plasma', plasma_gid, ier)
         CALL write_att_hdf5(plasma_gid,'active',qid_str,ier)
         CALL h5gcreate_f(plasma_gid,'opt_'//qid_str, qid_gid, ier)
         CALL h5gclose_f(qid_gid, ier)
         CALL h5gclose_f(plasma_gid, ier)

         CALL h5gcreate_f(fid,'neutral', neutral_gid, ier)
         CALL write_att_hdf5(neutral_gid,'active',qid_str,ier)
         CALL h5gcreate_f(neutral_gid,'opt_'//qid_str, qid_gid, ier)
         CALL h5gclose_f(qid_gid, ier)
         CALL h5gclose_f(neutral_gid, ier)

         CALL h5gcreate_f(fid,'wall', wall_gid, ier)
         CALL write_att_hdf5(wall_gid,'active',qid_str,ier)
         CALL h5gcreate_f(wall_gid,'opt_'//qid_str, qid_gid, ier)
         CALL h5gclose_f(qid_gid, ier)
         CALL h5gclose_f(wall_gid, ier)

         CALL h5gcreate_f(fid,'marker', marker_gid, ier)
         CALL write_att_hdf5(marker_gid,'active',qid_str,ier)
         CALL h5gcreate_f(marker_gid,'opt_'//qid_str, qid_gid, ier)
         CALL h5gclose_f(qid_gid, ier)
         CALL h5gclose_f(marker_gid, ier)

         ! Close the HDF5 file
         CALL close_hdf5(fid,ier)
         IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'beams3d_ascot5_'//TRIM(id_string)//'.h5',ier)
         WRITE(6,*) 'DONE DONE DONE'
         STOP
      END IF

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_write_ascoth5
