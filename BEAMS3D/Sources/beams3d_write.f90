!-----------------------------------------------------------------------
!     Module:        beams3d_write
!     Authors:       S. Lazerson (lazerson@pppl.gov) M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          06/21/2012
!     Description:   This subroutine outputs the particle trajectory data to an
!                    HDF5 file or binary file.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_write(write_type)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
#if defined(LHDF5)
      USE ez_hdf5
#endif
      USE beams3d_lines
      USE beams3d_grid, ONLY: nr, nphi, nz, B_R, B_PHI, B_Z, raxis, &
                                 zaxis, phiaxis, S_ARR, U_ARR, POT_ARR, &
                                 ZEFF_ARR, TE, TI, NE, wall_load, wall_shine, &
                                 plasma_mass, plasma_Zmean, &
                                 B_kick_min, B_kick_max, freq_kick, E_kick, NI
      USE beams3d_runtime, ONLY: id_string, npoinc, nbeams, beam, t_end, lverb, &
                                    lvmec, lpies, lspec, lcoil, lmgrid, lbeam, lascot, &
                                    lvessel, lvac, lbeam_simple, handle_err, nparticles_start, &
                                    HDF5_OPEN_ERR,HDF5_WRITE_ERR,&
                                    HDF5_CLOSE_ERR, BEAMS3D_VERSION, weight, e_beams, p_beams,&
                                    charge, Zatom, mass, ldepo, v_neut,lcollision, lfusion, &
                                    leqdsk, eqdsk_string, lhint, lhitonly, lkick, NION
      USE safe_open_mod, ONLY: safe_open
      USE wall_mod, ONLY: nface,nvertex,face,vertex,ihit_array
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
      INTEGER :: ier, iunit
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (myworkid == master) THEN
         SELECT CASE (TRIM(write_type))
            CASE('GRID_INIT')
#if defined(LHDF5)
               CALL open_hdf5('beams3d_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.true.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'beams3d_'//TRIM(id_string)//'.h5',ier)
               CALL write_scalar_hdf5(fid,'VERSION',ier,DBLVAR=BEAMS3D_VERSION,ATT='Version Number',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'VERSION',ier)
               CALL write_scalar_hdf5(fid,'lvmec',ier,BOOVAR=lvmec,ATT='VMEC input',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lvmec',ier)
               CALL write_scalar_hdf5(fid,'lpies',ier,BOOVAR=lpies,ATT='PIES input',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lpies',ier)
               CALL write_scalar_hdf5(fid,'lspec',ier,BOOVAR=lspec,ATT='SPEC input',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lspec',ier)
               CALL write_scalar_hdf5(fid,'leqdsk',ier,BOOVAR=leqdsk,ATT='EQDSK input',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'leqdsk',ier)
               CALL write_scalar_hdf5(fid,'lhint',ier,BOOVAR=leqdsk,ATT='HINT input',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lhint',ier)
               CALL write_scalar_hdf5(fid,'lcoil',ier,BOOVAR=lcoil,ATT='Coil input',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lcoil',ier)
               CALL write_scalar_hdf5(fid,'lmgrid',ier,BOOVAR=lmgrid,ATT='MGRID input',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lmgrid',ier)
               CALL write_scalar_hdf5(fid,'lvessel',ier,BOOVAR=lvessel,ATT='Vessel input',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lvessel',ier)
               CALL write_scalar_hdf5(fid,'lvac',ier,BOOVAR=lvac,ATT='Vacuum calc',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lvac',ier)
               CALL write_scalar_hdf5(fid,'lbeam',ier,BOOVAR=lbeam,ATT='Neutral Beam Calc',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lbeam',ier)
               CALL write_scalar_hdf5(fid,'lbeam_simple',ier,BOOVAR=lbeam_simple,ATT='Simple Beam Energy',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lbeam_simple',ier)
               CALL write_scalar_hdf5(fid,'ldepo',ier,BOOVAR=ldepo,ATT='Only Deposition Flag',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ldepo',ier)
               CALL write_scalar_hdf5(fid,'lkick',ier,BOOVAR=lkick,ATT='Kick Model Flag',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lkick',ier)
               CALL write_scalar_hdf5(fid,'lcollision',ier,BOOVAR=lcollision,ATT='Collisionall Operators Flag',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lcollision',ier)
               CALL write_scalar_hdf5(fid,'lascot',ier,BOOVAR=lascot,ATT='ASCOT5 Output Flag',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lascot',ier)
               CALL write_scalar_hdf5(fid,'lfusion',ier,BOOVAR=lfusion,ATT='Fusion Birth Model Flag',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lfusion',ier)
               CALL write_scalar_hdf5(fid,'lhitonly',ier,BOOVAR=lhitonly,ATT='Flag for only saving wall hits',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lhitonly',ier)
               CALL write_scalar_hdf5(fid,'nr',ier,INTVAR=nr,ATT='Number of Radial Gridpoints',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nr',ier)
               CALL write_scalar_hdf5(fid,'nphi',ier,INTVAR=nphi,ATT='Number of Toroidal Gridpoints',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nphi',ier)
               CALL write_scalar_hdf5(fid,'nz',ier,INTVAR=nz,ATT='Number of Vertical Gridpoints',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nz',ier)
               CALL write_scalar_hdf5(fid,'plasma_mass',ier,DBLVAR=plasma_mass,ATT='Plasma Mass [kg]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'plasma_mass',ier)
               CALL write_scalar_hdf5(fid,'plasma_Zmean',ier,DBLVAR=plasma_Zmean,ATT='Plasma [Z]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'plasma_Zmean',ier)
               CALL write_var_hdf5(fid,'raxis',nr,ier,DBLVAR=raxis,ATT='Radial Axis [m]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'raxis',ier)
               CALL write_var_hdf5(fid,'phiaxis',nphi,ier,DBLVAR=phiaxis,ATT='Toroidal Axis [rad]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'phiaxis',ier)
               CALL write_var_hdf5(fid,'zaxis',nz,ier,DBLVAR=zaxis,ATT='Vertical Axis [rad]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'zaxis',ier)
               CALL write_var_hdf5(fid,'B_R',nr,nphi,nz,ier,DBLVAR=B_R,ATT='Radial Trajectory Eq. (BR)',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'B_R',ier)
               CALL write_var_hdf5(fid,'B_Z',nr,nphi,nz,ier,DBLVAR=B_Z,ATT='Vertical Trajectory Eq. (BZ)',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'B_Z',ier)
               CALL write_var_hdf5(fid,'B_PHI',nr,nphi,nz,ier,DBLVAR=B_PHI,ATT='Toroidal Trajectory (BPHI)',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'B_PHI',ier)
               CALL write_var_hdf5(fid,'S_ARR',nr,nphi,nz,ier,DBLVAR=S_ARR,ATT='Normalized Toroidal Flux',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'S_ARR',ier)
               CALL write_var_hdf5(fid,'U_ARR',nr,nphi,nz,ier,DBLVAR=U_ARR,ATT='Equilibrium Poloidal Angle [rad]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'U_ARR',ier)
               CALL write_var_hdf5(fid,'POT_ARR',nr,nphi,nz,ier,DBLVAR=POT_ARR,ATT='Electrostatic Potential [V]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'POT_ARR',ier)
               IF (ALLOCATED(GFactor)) THEN
                  CALL write_var_hdf5(fid,'GFactor',ns_prof1,ier,DBLVAR=GFactor,ATT='(1-l31)/Zeff',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'GFactor',ier)
               END IF
               IF (ASSOCIATED(TE)) THEN
                  CALL write_var_hdf5(fid,'TE',nr,nphi,nz,ier,DBLVAR=TE,ATT='Electron Temperature [eV]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'TE',ier)
               END IF
               IF (ASSOCIATED(NE)) THEN
                  CALL write_var_hdf5(fid,'NE',nr,nphi,nz,ier,DBLVAR=NE,ATT='Electron Density [m^-3]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'NE',ier)
               END IF
               IF (ASSOCIATED(NI)) THEN
                  CALL write_var_hdf5(fid,'NI',nion,nr,nphi,nz,ier,DBLVAR=NI,ATT='Ion Densities [m^-3]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'NI',ier)
               END IF
               IF (ASSOCIATED(TI)) THEN
                  CALL write_var_hdf5(fid,'TI',nr,nphi,nz,ier,DBLVAR=TI,ATT='Ion Temperature [eV]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'TI',ier)
               END IF
               IF (ASSOCIATED(ZEFF_ARR)) THEN
                  CALL write_var_hdf5(fid,'ZEFF_ARR',nr,nphi,nz,ier,DBLVAR=ZEFF_ARR,ATT='Effective Ion Charge',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ZEFF_ARR',ier)
               END IF
               IF (ASSOCIATED(vertex)) THEN
                  CALL write_scalar_hdf5(fid,'nvertex',ier,INTVAR=nvertex,ATT='Number of Wall Vertices',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nvertex',ier)
                  CALL write_var_hdf5(fid,'wall_vertex',nvertex,3,ier,DBLVAR=vertex,ATT='Wall Verticies (x,y,z) [m]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'wall_vertex',ier)
               END IF
               IF (ASSOCIATED(face)) THEN
                  CALL write_scalar_hdf5(fid,'nface',ier,INTVAR=nface,ATT='Number of Wall Faces',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nface',ier)
                  CALL write_var_hdf5(fid,'wall_faces',nface,3,ier,INTVAR=face,ATT='Wall Faces',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'wall_faces',ier)
               END IF
               IF (lkick) THEN
                  CALL write_scalar_hdf5(fid,'B_kick_min',ier,DBLVAR=B_kick_min,ATT='|B|_min Kick Model [T]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'B_kick_min',ier)
                  CALL write_scalar_hdf5(fid,'B_kick_max',ier,DBLVAR=B_kick_max,ATT='|B|_max Kick Model [T]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'B_kick_max',ier)
                  CALL write_scalar_hdf5(fid,'freq_kick',ier,DBLVAR=freq_kick,ATT='Frequency Kick Model [Hz]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'freq_kick',ier)
                  CALL write_scalar_hdf5(fid,'E_kick',ier,DBLVAR=E_kick,ATT='E-Field Kick Model [V/m]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'E_kick',ier)
               END IF
            CASE('TRAJECTORY_PARTIAL')
               CALL open_hdf5('beams3d_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.false.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'beams3d_'//TRIM(id_string)//'.h5',ier)
               CALL write_scalar_hdf5(fid,'nparticles',ier,INTVAR=nparticles,ATT='Number of Trajectories',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nparticles',ier)
               CALL write_scalar_hdf5(fid,'nbeams',ier,INTVAR=nbeams,ATT='Number of Beams',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nbeams',ier)
               CALL write_scalar_hdf5(fid,'nsteps',ier,INTVAR=nsteps+1,ATT='Number of Steps Along Trajectory',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nsteps',ier)
               CALL write_scalar_hdf5(fid,'npoinc',ier,INTVAR=npoinc,ATT='Number of steps per trajectory period',&
                                      ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'npoinc',ier)
               !CALL write_var_hdf5(fid,'t_end',nparticles,ier,DBLVAR=t_end,ATT='Time at End of Trajectory [s]',&
               !                    ATT_NAME='description')
               !IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'t_end',ier)
               CALL write_var_hdf5(fid,'mass',nparticles,ier,DBLVAR=mass,ATT='Particle Mass [kg]',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'mass',ier)
               CALL write_var_hdf5(fid,'charge',nparticles,ier,DBLVAR=charge,ATT='Particle Charge [C]',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'charge',ier)
               CALL write_var_hdf5(fid,'Weight',nparticles,ier,DBLVAR=weight,ATT='Weight',&
                                  ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'Weight',ier)
               CALL write_var_hdf5(fid,'Beam',nparticles,ier,INTVAR=beam,ATT='Beam Number',&
                                      ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'Beam',ier)
               CALL write_var_hdf5(fid,'Zatom',nparticles,ier,DBLVAR=Zatom,ATT='Particle Charge Number',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'Zatom',ier)
               CALL write_var_hdf5(fid,'end_state',nparticles,ier,INTVAR=end_state,ATT='0: Orbiting; 1: Thermalized; 2: Wall Strike; 3: Shine-through; 4: Port-Load',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'end_state',ier)
               CALL write_var_hdf5(fid,'V_NEUT',3,nparticles,ier,DBLVAR=V_NEUT,ATT='Neutral Velocity (Vx, Vy,Vz) [m/s]',&
                                      ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'V_NEUT',ier)
               CALL write_var_hdf5(fid,'Energy',nbeams,ier,DBLVAR=e_beams,ATT='Beam Energy [J]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'E_BEAMS',ier)
               IF (ASSOCIATED(ihit_array)) THEN
                  CALL write_var_hdf5(fid,'wall_strikes',nface,ier,INTVAR=ihit_array,&
                                   ATT='Wall Strikes',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'wall_strikes',ier)
               END IF
               IF (ASSOCIATED(wall_load)) THEN
                  CALL write_var_hdf5(fid,'wall_load',nbeams,nface,ier,DBLVAR=wall_load,&
                                   ATT='Wall Loads [W/m^2]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'wall_load',ier)
               END IF
               IF (ASSOCIATED(wall_shine)) THEN
                  CALL write_var_hdf5(fid,'wall_shine',nbeams,nface,ier,DBLVAR=wall_shine,&
                                   ATT='Neutral Beam Shine-through [W/m^2]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'wall_shine',ier)
               END IF
            CASE('TRAJECTORY_FULL')
               CALL open_hdf5('beams3d_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.false.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'beams3d_'//TRIM(id_string)//'.h5',ier)
               CALL write_scalar_hdf5(fid,'nparticles',ier,INTVAR=nparticles,ATT='Number of Trajectories',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nparticles',ier)
               CALL write_scalar_hdf5(fid,'nbeams',ier,INTVAR=nbeams,ATT='Number of Beams',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nbeams',ier)
               CALL write_scalar_hdf5(fid,'nsteps',ier,INTVAR=nsteps+1,ATT='Number of Steps Along Trajectory',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nsteps',ier)
               CALL write_scalar_hdf5(fid,'npoinc',ier,INTVAR=npoinc,ATT='Number of steps per trajectory period',&
                                      ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'npoinc',ier)
               CALL write_var_hdf5(fid,'t_end',nparticles,ier,DBLVAR=t_end,ATT='Time at End of Trajectory [s]',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'t_end',ier)
               CALL write_var_hdf5(fid,'mass',nparticles,ier,DBLVAR=mass,ATT='Particle Mass [kg]',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'mass',ier)
               CALL write_var_hdf5(fid,'charge',nparticles,ier,DBLVAR=charge,ATT='Particle Charge [C]',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'charge',ier)
               CALL write_var_hdf5(fid,'Weight',nparticles,ier,DBLVAR=weight,ATT='Weight',&
                                  ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'Weight',ier)
               CALL write_var_hdf5(fid,'Beam',nparticles,ier,INTVAR=beam,ATT='Beam Number',&
                                      ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'Beam',ier)
               CALL write_var_hdf5(fid,'Zatom',nparticles,ier,DBLVAR=Zatom,ATT='Particle Charge Number',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'Zatom',ier)
               CALL write_var_hdf5(fid,'end_state',nparticles,ier,INTVAR=end_state,ATT='0: Orbiting; 1: Thermalized; 2: Wall Strike; 3: Shienthrough',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'end_state',ier)
               CALL write_var_hdf5(fid,'V_NEUT',3,nparticles,ier,DBLVAR=V_NEUT,ATT='Neutral Velocity (Vx, Vy,Vz) [m/s]',&
                                      ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'V_NEUT',ier)
               CALL write_var_hdf5(fid,'Energy',nbeams,ier,DBLVAR=e_beams,ATT='Beam Energy [J]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'E_BEAMS',ier)
               CALL write_var_hdf5(fid,'R_lines',npoinc+1,nparticles,ier,DBLVAR=R_lines,ATT='Cylindrical R of Trajectory [m]',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'R_lines',ier)
               CALL write_var_hdf5(fid,'Z_lines',npoinc+1,nparticles,ier,DBLVAR=Z_lines,ATT='Cylindrical Z of Trajectory [m]',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'Z_lines',ier)
               CALL write_var_hdf5(fid,'PHI_lines',npoinc+1,nparticles,ier,DBLVAR=PHI_lines,ATT='Cylindrical Phi of Trajectory [rad]',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'PHI_lines',ier)
               CALL write_var_hdf5(fid,'vll_lines',npoinc+1,nparticles,ier,DBLVAR=vll_lines,ATT='Parallel Particle Velocity [m/s]',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'vll_lines',ier)
               CALL write_var_hdf5(fid,'neut_lines',npoinc+1,nparticles,ier,BOOVAR=neut_lines,ATT='Neutral Indicator [1=Neut.]',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'neut_lines',ier)
               CALL write_var_hdf5(fid,'moment_lines',npoinc+1,nparticles,ier,DBLVAR=moment_lines,&
                                   ATT='Magnetic Moment [kg m^2 /s^2 T ]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'moment_lines',ier)
               CALL write_var_hdf5(fid,'S_lines',npoinc+1,nparticles,ier,DBLVAR=S_lines,ATT='Toroidal Flux Coordinate',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'S_lines',ier)
               CALL write_var_hdf5(fid,'U_lines',npoinc+1,nparticles,ier,DBLVAR=U_lines,ATT='U Flux Coordinate',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'U_lines',ier)
               CALL write_var_hdf5(fid,'B_lines',npoinc+1,nparticles,ier,DBLVAR=B_lines,ATT='|B| along Fieldline',&
                                   ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'B_lines',ier)
               IF (ASSOCIATED(ihit_array)) THEN
                  CALL write_var_hdf5(fid,'wall_strikes',nface,ier,INTVAR=ihit_array,&
                                   ATT='Wall Strikes',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'wall_strikes',ier)
               END IF
               IF (ASSOCIATED(wall_load)) THEN
                  CALL write_var_hdf5(fid,'wall_load',nbeams,nface,ier,DBLVAR=wall_load,&
                                   ATT='Wall Loads [W/m^2]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'wall_strikes',ier)
               END IF
               IF (ASSOCIATED(wall_shine)) THEN
                  CALL write_var_hdf5(fid,'wall_shine',nbeams,nface,ier,DBLVAR=wall_shine,&
                                   ATT='Neutral Beam Shine-through [W/m^2]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'wall_shine',ier)
               END IF
            CASE('DIAG')
               CALL open_hdf5('beams3d_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.false.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'beams3d_'//TRIM(id_string)//'.h5',ier)
               CALL write_scalar_hdf5(fid,'partvmax',ier,DBLVAR=partvmax,&
                                      ATT='Maximum velocity of dist func [m/s]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'partvmax',ier)
               CALL write_scalar_hdf5(fid,'ns_prof1',ier,INTVAR=ns_prof1,&
                                   ATT='Rho Dist. Grid Points [0,1]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ns_prof1',ier)
               CALL write_scalar_hdf5(fid,'ns_prof2',ier,INTVAR=ns_prof2,&
                                   ATT='U Dist. Grid Points [0,2pi]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ns_prof2',ier)
               CALL write_scalar_hdf5(fid,'ns_prof3',ier,INTVAR=ns_prof3,&
                                   ATT='PHI Dist. Grid Points [0,2pi]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ns_prof3',ier)
               CALL write_scalar_hdf5(fid,'ns_prof4',ier,INTVAR=ns_prof4,&
                                   ATT='VLL Dist. Grid Points[-vmax,vmax]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ns_prof4',ier)
               CALL write_scalar_hdf5(fid,'ns_prof5',ier,INTVAR=ns_prof5,&
                                   ATT='Vperp Dist. Grid Points [0, vmax]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ns_prof5',ier)
               IF (ASSOCIATED(ndot_prof)) THEN
                  CALL write_var_hdf5(fid,'ndot_prof',nbeams,ns_prof1,ier,DBLVAR=ndot_prof,&
                                      ATT='Fast Ion Source [m^-3/s]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ndot_prof',ier)
               END IF
               IF (ASSOCIATED(epower_prof)) THEN
                  CALL write_var_hdf5(fid,'epower_prof',nbeams,ns_prof1,ier,DBLVAR=epower_prof,&
                                      ATT='Electron Power Deposition [W*m^-3]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'epower_prof',ier)
               END IF
               IF (ASSOCIATED(ipower_prof)) THEN
                  CALL write_var_hdf5(fid,'ipower_prof',nbeams,ns_prof1,ier,DBLVAR=ipower_prof,&
                                      ATT='Ion Power Deposition [W*m^-3]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ipower_prof',ier)
               END IF
               IF (ASSOCIATED(j_prof)) THEN
                  CALL write_var_hdf5(fid,'j_prof',nbeams,ns_prof1,ier,DBLVAR=j_prof,&
                                      ATT='Total Beam Current Density [A*m^-2]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'j_prof',ier)
               END IF
               IF (ASSOCIATED(dense_prof)) THEN
                  CALL write_var_hdf5(fid,'dense_prof',nbeams,ns_prof1,ier,DBLVAR=dense_prof,&
                                      ATT='Fast Ion Density [m^-3]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'dense_prof',ier)
               END IF
               IF (ASSOCIATED(dist5d_prof)) THEN
                  CALL write_var_hdf5(fid,'dist_prof',nbeams,ns_prof1,ns_prof2,ns_prof3,ns_prof4,ns_prof5,ier,DBLVAR=dist5d_prof,&
                                      ATT='Distribution Function [part/(m^6/s^3)] (nbeam,nrho,npol,ntor,nvll,nvperp)',ATT_NAME='description') !its not volume normalized, so shouldnt the units be [part]?
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'dist_prof',ier)
               END IF
               IF (lbeam) THEN
                  CALL write_var_hdf5(fid,'Shinethrough',nbeams,ier,DBLVAR=shine_through,&
                                   ATT='Total Beam Shine Through [%]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'Shinethrough',ier)
                  CALL write_var_hdf5(fid,'Shineport',nbeams,ier,DBLVAR=shine_port,&
                                   ATT='Loss to Port [%]',ATT_NAME='description')
                  IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'shine_port',ier)
               END IF
         END SELECT
         CALL close_hdf5(fid,ier)
         IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'beams3d_'//TRIM(id_string)//'.h5',ier)
      END IF
#else
      WRITE(6,'(A)')  '   FILE: '//'beams3d_'//TRIM(id_string)//'.bin'
      CALL safe_open(iunit,ier,'beams3d_'//TRIM(id_string)//'.bin','replace','unformatted')
      WRITE(iunit,*) BEAMS3D_VERSION
      WRITE(iunit,*) lvmec,lpies,lspec,lcoil,lmgrid,lvessel,lvac,lbeam_simple
      WRITE(iunit,*) nparticles,nsteps,npoinc,nbeams
      WRITE(iunit,*) weight
      WRITE(iunit,*) beam
      WRITE(iunit,*) t_end
      WRITE(iunit,*) R_lines
      WRITE(iunit,*) Z_lines
      WRITE(iunit,*) PHI_lines
      WRITE(iunit,*) vll_lines
      WRITE(iunit,*) S_lines
      WRITE(iunit,*) U_lines
      WRITE(iunit,*) V_lines
      IF (lbeam) THEN
         WRITE(iunit,*) weight
         WRITE(iunit,*) beam
         WRITE(iunit,*) e_beams
         IF (.not.ldepo) THEN
            WRITE(iunit,*) shine_through
            WRITE(iunit,*) ndot_prof
            WRITE(iunit,*) epower_prof
            WRITE(iunit,*) ipower_prof
            WRITE(iunit,*) j_prof
         END IF
      END IF
      WRITE(iunit,*) nr,nphi,nz
      WRITE(iunit,*) raxis
      WRITE(iunit,*) phiaxis
      WRITE(iunit,*) zaxis
      WRITE(iunit,*) B_R
      WRITE(iunit,*) B_Z
      WRITE(iunit,*) B_PHI
      CLOSE(iunit)
#endif

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_write
