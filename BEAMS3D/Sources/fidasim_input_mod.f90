

MODULE fidasim_input_mod

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
      nsh_prof4,  r_h, p_h, z_h, e_h, pi_h, h1_prof
   USE beams3d_grid, ONLY: nr, nphi, nz, B_R, B_PHI, B_Z, raxis, &
      zaxis, phiaxis, POT_ARR, &
      TE, TI, NE, npot, nte, nti, &
      POT4D, NE4D, TE4D, TI4D, ZEFF4D, &
      BR4D, BPHI4D, BZ4D, VTOR4D,&
      hr, hp, hz, hri, hpi, hzi, U4D, &
      rmin, rmax,  phimin, phimax, &
      rmin_fida, rmax_fida, zmin_fida, zmax_fida, phimin_fida, phimax_fida, &
      raxis_fida, zaxis_fida, phiaxis_fida, nr_fida, nphi_fida, nz_fida, &
      nenergy_fida, npitch_fida, energy_fida, pitch_fida, t_fida,&
      nne, nte, nti, nzeff, nvtor
   USE beams3d_runtime
   ! , ONLY: id_string, nbeams, beam, lverb, handle_err, &
   !    HDF5_OPEN_ERR,HDF5_WRITE_ERR,HDF5_CLOSE_ERR, BEAMS3D_VERSION, weight, &
   !    charge, mass,pi, pi2, mass_beams, lfidasim_cyl, lsplit,TE_AUX_S, TE_AUX_F, NE_AUX_S, NE_AUX_F, TI_AUX_S, TI_AUX_F,&
   !    POT_AUX_S, POT_AUX_F, ZEFF_AUX_S, ZEFF_AUX_F
   USE safe_open_mod, ONLY: safe_open
   USE beams3d_write_par
   USE beams3d_physics_mod, ONLY: beams3d_SUFLX
   USE mpi_params
   USE mpi_inc



   !!!! Namelist
   INTEGER, parameter :: MAXCHAN = 512
   INTEGER :: ier, iunit,istat
   INTEGER(HID_T) ::  qid_gid,temp_gid
   LOGICAL :: lexist
   CHARACTER(LEN=1000) :: line
   CHARACTER(128) :: comment,nbi_data_source, spec_data_source,runid,device,name,system,&
      tables_file,equilibrium_file,geometry_file,distribution_file,neutrals_file,result_dir
   CHARACTER(20), DIMENSION(MAXCHAN) :: id
   DOUBLE PRECISION, DIMENSION(3) :: origin,current_fractions,src,axis_nbi,  divy,   divz, q
   DOUBLE PRECISION, DIMENSION(3,MAXCHAN) :: axis_spec,  lens
   DOUBLE PRECISION, DIMENSION(MAXCHAN) :: sigma_pi,spot_size,radius
   DOUBLE PRECISION  :: time, lambdamin,  lambdamax,   alpha,   beta,   gamma,      xminf,&
      xmax,   ymin,   ymax,   zmin,   zmax, ab,   ai,   pinj,   einj,&
      focy,focz,widz,widy,emax_wght,lambdamin_wght, lambdamax_wght
   INTEGER ::  nlambda, nx_beam, ny_beam,nz_beam, nz_tmp,impurity_charge, ne_wght, np_wght, nphi_wght, nlambda_wght,& !!!! nz, 
      calc_npa,   calc_fida,   calc_pnpa,   calc_pfida,   calc_bes,   calc_dcx,   calc_halo, &
      calc_cold,   calc_brems,   calc_birth,   calc_fida_wght,   calc_npa_wght,   calc_neutron,&
      shape, load_neutrals,verbose,flr,naperture,nchan,namelist_present
   INTEGER, DIMENSION(10) :: ashape
   DOUBLE PRECISION, DIMENSION(10) :: awidy,awidz, aoffy, aoffz, adist
   INTEGER(HID_T) :: shot,  n_fida,   n_nbi,   n_pfida,   n_pnpa,   n_dcx,   n_npa,   n_halo,   n_birth,&
      seed


   NAMELIST /fidasim_inputs_b3d/ comment,tables_file, equilibrium_file,geometry_file,distribution_file,&
      result_dir,neutrals_file,shot, time, runid, device, nlambda,seed,load_neutrals,verbose,flr,&
      lambdamin, lambdamax, nx, ny, nz, alpha, beta, gamma, origin,&
      xmin, xmax, ymin, ymax, zmin, zmax, ab, ai, current_fractions,&
      pinj, einj, impurity_charge, n_fida, n_nbi, n_pfida, n_pnpa,&
      n_dcx, n_npa, n_halo, n_birth, ne_wght, np_wght, nphi_wght,&
      emax_wght, nlambda_wght, lambdamin_wght, lambdamax_wght,&
      calc_npa, calc_fida, calc_pnpa, calc_pfida, calc_bes, calc_dcx,&
      calc_halo, calc_cold, calc_brems, calc_birth, calc_fida_wght, calc_npa_wght,&
      calc_neutron,&
      nbi_data_source,name,shape,src,axis_nbi,divy,divz,focy,focz,widz,widy,naperture,ashape,awidy,awidz, aoffy, aoffz, adist,&
      spec_data_source,nchan,system,id,lens,axis_spec,sigma_pi,spot_size,radius

!    NAMELIST /fidasim_inputs/ tables_file, equilibrium_file,geometry_file,distribution_file,&
!       result_dir,neutrals_file,shot, time, runid,  nlambda,seed,load_neutrals,verbose,flr,&
!       lambdamin, lambdamax, nx, ny, nz, alpha, beta, gamma, origin,&
!       xmin, xmax, ymin, ymax, zmin, zmax, ab, ai, current_fractions,&
!       pinj, einj, impurity_charge, n_fida, n_nbi, n_pfida, n_pnpa,&
!       n_dcx, n_npa, n_halo, n_birth, ne_wght, np_wght, nphi_wght,&
!       emax_wght, nlambda_wght, lambdamin_wght, lambdamax_wght,&
!       calc_npa, calc_fida, calc_pnpa, calc_pfida, calc_bes, calc_dcx,&
!       calc_halo, calc_cold, calc_brems, calc_birth, calc_fida_wght, calc_npa_wght,&
!       calc_neutron!comment,device,,id


      CONTAINS
!-----------------------------------------------------------------------
!     Module:        beams3d_write_fidasim
!     Authors:       D. Kulla (david.kulla@ipp.mpg.de)
!     Date:          04/01/2021
!     Description:   This subroutine outputs simulation data for the
!                    FIDASIM code.
!-----------------------------------------------------------------------
SUBROUTINE beams3d_write_fidasim(write_type)

#ifdef LHDF5
   USE hdf5
   USE ez_hdf5
#endif

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
   INTEGER :: ier, iunit,istat, i, j, d1, d2, d3, k, k1, k2, kmax ,ider, &
      l, m, n, b, i3, j3, k3
   !INTEGER(HID_T) ::    qid_gid, temp_gid
   INTEGER(HID_T) :: qid_gid2
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: mask

   REAL*8 :: fvalE(1,3), fval(1), fval2(1), xparam, yparam, zparam
   REAL(rprec) :: jac, v_parr, v_perp, pitch, v
   REAL(rprec), DIMENSION(4) :: rt,zt,pt
   REAL(rprec), DIMENSION(:,:,:,:,:), POINTER :: dist5d_temp

   DOUBLE PRECISION         :: x0, y0, z0, vol
   DOUBLE PRECISION, ALLOCATABLE :: rtemp(:,:,:), rtemp2(:,:,:), rtemp3(:,:,:), rtemp4(:,:,:), r1dtemp(:), r2dtemp(:,:), r4dtemp(:,:,:,:)

   CHARACTER(LEN=8) :: temp_str8

   INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/), ictE(8)=(/0,1,1,1,0,0,0,0/)
   DOUBLE PRECISION, PARAMETER :: e_charge      = 1.602176565e-19 !e_c
   DOUBLE PRECISION, PARAMETER :: zero          = 0.0D0 ! 0.0
   DOUBLE PRECISION, PARAMETER :: t_min          = 1.0D-3 !

   !-----------------------------------------------------------------------
   !     Begin Subroutine
   !-----------------------------------------------------------------------
   IF (myworkid == master) THEN
      SELECT CASE (TRIM(write_type))
         CASE('INIT')
            
            CALL init_fidasim_input
            nz_tmp=nz !This is important for grid interpolation and necessary to prevent renaming the namelist
            ! Read namelist
            istat=0
            iunit=12
            INQUIRE(FILE=TRIM('input.' // TRIM(id_string)),EXIST=lexist)
            IF (.not.lexist) stop 'Could not find input file'
            CALL safe_open(iunit,istat,'input.' // TRIM(id_string),'old','formatted')
            IF (istat /= 0) CALL handle_err(NAMELIST_READ_ERR,'fidasim_inputs_b3d in: ' //'input.' //TRIM(id_string),istat)
            READ(iunit,NML=fidasim_inputs_b3d,IOSTAT=istat)
            IF (istat /= 0) THEN
               backspace(iunit)
               read(iunit,fmt='(A)') line
               write(6,'(A)') 'Invalid line in namelist: '//TRIM(line)
               CALL handle_err(NAMELIST_READ_ERR,'fidasim_inputs_b3d in: '//'input.' //TRIM(id_string),istat)
            END IF
            CLOSE(iunit)
            IF (equilibrium_file .eq. 'none')   equilibrium_file = TRIM(result_dir) // TRIM(runid) //  '_equilibrium.h5'   !! File containing plasma parameters and fields
            IF (geometry_file .eq. 'none')  geometry_file = TRIM(result_dir) // TRIM(runid) //  '_geometry.h5'      !! File containing NBI and diagnostic geometry
            IF (distribution_file .eq. 'none') distribution_file = TRIM(result_dir) // TRIM(runid) //  '_distribution.h5'      !! File containing fast-ion distribution
            IF (neutrals_file .eq. 'none') neutrals_file = TRIM(result_dir) // TRIM(runid) //  '_neutrals.h5'

            ! Open the inputs.dat file
            iunit = 10
            CALL safe_open(iunit,istat,'fidasim_'//TRIM(id_string)//'_inputs.dat','replace','formatted')
            !WRITE(iunit,nml=fidasim_inputs)
            CALL write_fidasim_namelist(iunit, istat)
            CALL FLUSH(iunit)
            CLOSE(iunit)
   
            nz=nz_tmp !This is important for grid interpolation and necessary to prevent renaming the namelist


            CALL write_fidasim_geometry

            CALL open_hdf5('fidasim_'//TRIM(id_string)//'_distribution.h5',fid,ier,LCREATE=.true.)
            CALL h5gopen_f(fid,'/', qid_gid, ier)
            CALL write_att_hdf5(qid_gid,'data_source','Data initialized from BEAMS3D',ier)
            CALL write_var_hdf5(qid_gid,'type',ier,INTVAR=1)
            CALL h5dopen_f(fid, 'type', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'description','Distribution type: 1="Guiding Center Density Function", 2="Guiding Center Monte Carlo", 3="Full Orbit Monte Carlo"',ier)
            CALL h5dclose_f(temp_gid,ier)


            ! FIDASIM GRID
            CALL write_fidasim_grid(qid_gid)

            !Energy grid is only in distribution
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

            ! Close file
            CALL h5gclose_f(qid_gid, ier)
            CALL close_hdf5(fid,ier)
            IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'fidasim_'//TRIM(id_string)//'_distribution.h5',ier)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            CALL write_fidasim_equilibrium
       

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         CASE('DENF')
            ! Distribution needs to be in density [part/m^3] here - called from distnorm (before velocity space normalization)
            !ALLOCATE(r4dtemp(nbeams,nr_fida,nphi_fida,nz_fida))
            ALLOCATE(rtemp(nr_fida,nz_fida,nphi_fida))
            rtemp = 0.0
            IF (lfidasim_cyl) THEN
               DO i=1,nr_fida
                  vol =(raxis_fida(i) + 1 / 2.0 / r_h) / r_h / z_h / p_h
                  dist5d_fida(i,:,:,:,:) = dist5d_fida(i,:,:,:,:) / vol
                  DO j = 1, nz_fida
                     DO k=1,nphi_fida
                        rtemp(i,j,k) = SUM(dist5d_fida(i,j,k,:,:))
                     END DO
                  END DO
               END DO
            ELSE
               !convert to r z phi
               DO i=1,nr_fida
                  DO k = 1, nz_fida
                     DO j=1,nphi_fida
                        !Get S and U
                        q(1) = raxis_fida(i)
                        q(2) = phiaxis_fida(j)
                        q(3) = zaxis_fida(k)
                        CALL beams3d_SUFLX(q,y0,z0)

                        ! Skip if S > 1
                        IF (y0 > 1.0) CYCLE

                        ! Keep U positive
                        IF (z0 < 0) z0 = z0 + pi2

                        ! Create x0 (index over full toridal angle)
                        x0 = MOD(phiaxis_fida(j),pi2)
                        IF (x0<0) x0 = x0 + pi2

                        ! Calc dist func bins
                        l = MAX(MIN(CEILING(SQRT(y0)*h1_prof      ), ns_prof1), 1) ! Rho Bin
                        m = MAX(MIN(CEILING( z0*h2_prof           ), ns_prof2), 1) ! U Bin
                        n = MAX(MIN(CEILING( x0*h3_prof           ), ns_prof3), 1) ! V Bin
                        rtemp(i,k,j) = SUM(dist5d_prof(:,l,m,n,:,:))!output in r-z-phi

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
            ! Do phase space change of coordinates
            !Allocate with Radial-like dimensions for clean transfer and to avoid explicitly looping over every element
            IF (lfidasim_cyl) THEN
               ALLOCATE(dist5d_temp(nenergy_fida, npitch_fida,nr_fida,nz_fida,nphi_fida)) !need temp as velocity bins are in vll/vperp initially
               dist5d_temp = 0
            ELSE
               ALLOCATE(dist5d_fida(ns_prof1, ns_prof2, ns_prof3, nenergy_fida, npitch_fida)) !nenergy and npitch are always aligned to distribution
               dist5d_fida = 0
            END IF

            IF (lfidasim_cyl) THEN
               DO d1 = 1, nenergy_fida
                  DO d2 = 1, npitch_fida
                     dist5d_temp(d1,d2,:,:,:) = dist5d_fida(:,:,:,d1,d2)*e_h*pi_h*REAL(1.0e-6)
                  END DO
               END DO
            ELSE
               DO b = 1,nbeams
                  DO d1 = 1, nenergy_fida
                     DO d2 = 1, npitch_fida
                        v = SQRT(2 * energy_fida(d1) *1000.0 * e_charge / mass_beams(b))
                        pitch = pitch_fida(d2)
                        v_parr = pitch * v
                        v_perp = SQRT(1- pitch * pitch) * v
                        !determine beams3d-grid indices (velocity space)
                        i3 = MAX(MIN(1+nsh_prof4+FLOOR(h4_prof*v_parr), ns_prof4), 1) ! vll
                        j3 = MAX(MIN(CEILING(v_perp*h5_prof         ), ns_prof5), 1) ! Vperp
                        jac = pi2 * v / mass_beams(b) * e_charge / REAL(1000) ! * pi2
                        dist5d_fida(:,:,:,d1,d2) = dist5d_fida(:,:,:,d1,d2) + dist5d_prof(b,:,:,:,i3,j3) * jac ! conversion to final grid comes in next steps
                     END DO
                  END DO
               END DO
            END IF

            IF (.not. lfidasim_cyl) THEN
               ALLOCATE(dist5d_temp(ns_prof1, ns_prof2, ns_prof3, nenergy_fida, npitch_fida))
               dist5d_temp(:,:,:,:,:) = dist5d_fida(:,:,:,:,:)
               DEALLOCATE(dist5d_fida)
               !Interpolate rho u v to r phi z distribution function (nearest neighbor at the moment)
               !Now allocate with correct dimensions
               ALLOCATE(dist5d_fida(nenergy_fida,npitch_fida,nr_fida,nz_fida,nphi_fida))
               DO i=1,nr_fida
                  DO k = 1, nz_fida
                     DO j=1,nphi_fida
                        !Get S and U
                        q(1) = raxis_fida(i)
                        q(2) = phiaxis_fida(j)
                        q(3) = zaxis_fida(k)
                        CALL beams3d_SUFLX(q,y0,z0)

                        ! Skip if S > 1
                        IF (y0 > 1.0) CYCLE

                        ! Handle negative y0
                        IF (y0 < 0) THEN
                           z0 = z0 + pi2*0.5
                           y0 = -y0
                        END IF

                        ! Create x0 (index over full toridal angle)
                        x0 = MOD(phiaxis_fida(j),pi2)
                        IF (x0<0) x0 = x0 + pi2

                        ! Calc dist func bins
                        l = MAX(MIN(CEILING(SQRT(y0)*h1_prof     ), ns_prof1), 1) ! Rho Bin
                        m = MAX(MIN(CEILING( z0*h2_prof           ), ns_prof2), 1) ! U Bin
                        n = MAX(MIN(CEILING( x0*h3_prof           ), ns_prof3), 1) ! V Bin
                        dist5d_fida(:,:,i,k,j) = dist5d_temp(l,m,n,:,:) !output in r-z-phi

                     END DO
                  END DO
               END DO
            END IF

            CALL open_hdf5('fidasim_'//TRIM(id_string)//'_distribution.h5',fid,ier,LCREATE=.false.)
            IF (ASSOCIATED(dist5d_fida)) THEN
               IF (lfidasim_cyl) THEN
                  CALL write_var_hdf5(fid,'f',nenergy_fida,npitch_fida,nr_fida,nz_fida,nphi_fida,ier,DBLVAR=dist5d_temp)
               ELSE
                  CALL write_var_hdf5(fid,'f',nenergy_fida,npitch_fida,nr_fida,nz_fida,nphi_fida,ier,DBLVAR=dist5d_fida)
               END IF
               CALL h5dopen_f(fid, '/f', temp_gid, ier)
               CALL write_att_hdf5(temp_gid,'description','Distribution Function (nenergy_fida,npitch_fida,nr_fida,nz_fida,nphi_fida)',ier)
               CALL write_att_hdf5(temp_gid,'units','part/(cm^3 keV)',ier)
               CALL h5dclose_f(temp_gid,ier)
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'dist_fida',ier)
            END IF
            CALL close_hdf5(fid,ier)
            IF (.not. lfidasim_cyl) THEN
               DEALLOCATE(dist5d_fida)
            END IF
            DEALLOCATE(dist5d_temp)


         CASE('DISTRIBUTION_GC_MC')
         CASE('DISTRIBUTION_FO')

      END SELECT
   END IF

   RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
END SUBROUTINE beams3d_write_fidasim

SUBROUTINE init_fidasim_input

   !Default values
   shot = 0    !! Shot Number
   time = 0.0    !! Time [s]
   runid = 'hard-coded value!'    !! runID
   result_dir = 'hard-coded value!'   !! Result Directory
   comment = 'none'
   device =  'none'
   !! Input Files
   tables_file = 'none'  !! Atomic Tables File
   equilibrium_file = 'none'    !! File containing plasma parameters and fields
   geometry_file = 'none'    !! File containing NBI and diagnostic geometry
   distribution_file = 'none'    !! File containing fast-ion distribution
   neutrals_file = 'none'
   !! Simulation Switches
   calc_bes = 0    !! Calculate NBI Spectra
   calc_dcx = 0   !! Calculate Direct CX Spectra
   calc_halo = 0    !! Calculate Halo Spectra
   calc_cold = 0    !! Calculate Cold D-alpha Spectra
   calc_brems = 0    !! Calculate Bremsstrahlung
   calc_fida = 0    !! Calculate FIDA Spectra
   calc_npa = 0   !! Calculate NPA
   calc_pfida = 0    !! Calculate Passive FIDA Spectra
   calc_pnpa = 0   !! Calculate Passive NPA
   calc_neutron = 0   !! Calculate B-T Neutron Rate
   calc_birth = 0    !! Calculate Birth Profile
   calc_fida_wght = 0    !! Calculate FIDA weights
   calc_npa_wght = 0    !! Calculate NPA weights
   !! Debugging Switches
   seed = -1    !! RNG Seed. If seed is negative a random seed is used
   flr = 1    !! Turn on Finite Larmor Radius corrections
   load_neutrals = 0    !! Load neutrals from neutrals file
   verbose = 1    !! Verbose
   !! Monte Carlo Settings
   n_fida = 5000000    !! Number of FIDA mc particles
   n_npa = 5000000    !! Number of NPA mc particles
   n_pfida = 50000000    !! Number of Passive FIDA mc particles
   n_pnpa = 50000000    !! Number of Passive NPA mc particles
   n_nbi = 50000    !! Number of NBI mc particles
   n_halo = 50000    !! Number of HALO mc particles
   n_dcx = 500000     !! Number of DCX mc particles
   n_birth = 10000    !! Number of BIRTH mc particles
   !! Neutral Beam Settings
   ab = 0.000000     !! Beam Species mass [amu]
   pinj = 0.0     !! Beam Power [MW]
   einj = 0.0     !! Beam Energy [keV]
   current_fractions(1) = 0.0 !! Current Fractions (Full component)
   current_fractions(2) = 0.0 !! Current Fractions (Half component)
   current_fractions(3) = 0.0 !! Current Fractions (Third component)
   !! Plasma Settings
   ai = 0.000000     !! Ion Species mass [amu]
   impurity_charge = 5     !! Impurity Charge
   !! Beam Grid Settings
   nx_beam = 90    !! Number of cells in X direction (Into Plasma)
   ny_beam = 40    !! Number of cells in Y direction
   nz_beam = 40    !! Number of cells in Z direction
   xmin = 0.000000     !! Minimum X value [cm]
   xmax = 180.000000     !! Maximum X value [cm]
   ymin = -40.000000     !! Minimum Y value [cm]
   ymax = 40.000000     !! Maximum Y value [cm]
   zmin = -40.000000     !! Minimum Z value [cm]
   zmax = 40.000000     !! Maximum Z value [cm]
   !! Tait-Bryan Angles for z-y`-x`` rotation
   alpha = 0.0     !! Rotation about z-axis [rad]
   beta  = 0.0     !! Rotation about y`-axis [rad]
   gamma = 0.000000     !! Rotation about x``-axis [rad]
   !! Beam Grid origin in machine coordinates (cartesian)
   origin(1) = 0.0    !! U value [cm]
   origin(2) = 0.0    !! V value [cm]
   origin(3) = 0.0     !! W value [cm]
   !! Wavelength Grid Settings
   nlambda = 1024    !! Number of Wavelengths
   lambdamin = 647.000000    !! Minimum Wavelength [nm]
   lambdamax = 669.000000    !! Maximum Wavelength [nm]
   !! Weight Function Settings
   ne_wght = 10    !! Number of Energies for Weights
   np_wght = 10    !! Number of Pitches for Weights
   nphi_wght = 8    !! Number of Gyro-angles for Weights
   emax_wght = 150.000000    !! Maximum Energy for Weights [keV]
   nlambda_wght = 10    !! Number of Wavelengths for Weights
   lambdamin_wght = 647.000000    !! Minimum Wavelength for Weights [nm]
   lambdamax_wght = 669.000000    !! Maximum Wavelength for Weights [nm]

   nbi_data_source = 'none'
   name = 'none'
   shape = 1
   src = (/0.0,0.0,0.0/)
   axis_nbi = (/0.0,0.0,0.0/)
   divy = (/0.0,0.0,0.0/)
   divz = (/0.0,0.0,0.0/)
   focy= -1.0
   focz= -1.0
   widz= -1.0
   widy= -1.0
   spec_data_source = 'none'
   nchan=0
   system ='none'
   id = SPREAD('none',1,MAXCHAN)
   lens=-1.0
   axis_spec=-1.0
   sigma_pi = -1.0
   spot_size = -1.0
   radius = -1.0
   naperture = 1
   ashape = -1
   awidy = -1.0
   awidz = -1.0
   aoffy = -1.0
   aoffz = -1.0
   adist = -1.0

END SUBROUTINE init_fidasim_input

SUBROUTINE write_fidasim_namelist(iunit_out, istat)
      INTEGER, INTENT(in) :: iunit_out
      !INTEGER, INTENT(out) :: istat
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
      WRITE(iunit_out,'(A)') '&FIDASIM_INPUTS'
      WRITE(iunit_out,'(A)') '!---------- STRING VARIABLES ------------'
      WRITE(iunit_out,outstr) 'TABLES_FILE',TRIM(tables_file)
      WRITE(iunit_out,outstr) 'EQUILIBRIUM_FILE',TRIM(equilibrium_file)
      WRITE(iunit_out,outstr) 'GEOMETRY_FILE',TRIM(geometry_file)
      WRITE(iunit_out,outstr) 'DISTRIBUTION_FILE',TRIM(distribution_file)
      WRITE(iunit_out,outstr) 'NEUTRALS_FILE',TRIM(neutrals_file)
      WRITE(iunit_out,outstr) 'RESULT_DIR',TRIM(result_dir)
      WRITE(iunit_out,outstr) 'RUNID',TRIM(runid)
      WRITE(iunit_out,'(A)') '!---------- INTEGER VARIABLES ------------'
      WRITE(iunit_out,outint) 'SHOT',shot
      WRITE(iunit_out,outint) 'SEED',seed
      WRITE(iunit_out,outint) 'FLR',flr
      WRITE(iunit_out,outint) 'VERBOSE',verbose
      WRITE(iunit_out,outint) 'NLAMBDA',nlambda
      WRITE(iunit_out,outint) 'LOAD_NEUTRALS',load_neutrals
      WRITE(iunit_out,outint) 'NX',nx_beam
      WRITE(iunit_out,outint) 'NY',ny_beam
      WRITE(iunit_out,outint) 'NZ',nz_beam
      WRITE(iunit_out,outint) 'N_NBI',n_nbi
      WRITE(iunit_out,outint) 'N_DCX',n_dcx
      WRITE(iunit_out,outint) 'N_HALO',n_halo
      WRITE(iunit_out,outint) 'N_BIRTH',n_birth
      WRITE(iunit_out,outint) 'N_FIDA',n_fida
      WRITE(iunit_out,outint) 'N_PFIDA',n_pfida
      WRITE(iunit_out,outint) 'N_NPA',n_npa
      WRITE(iunit_out,outint) 'NE_WGHT',ne_wght
      WRITE(iunit_out,outint) 'NP_WGHT',np_wght
      WRITE(iunit_out,outint) 'NPHI_WGHT',nphi_wght
      WRITE(iunit_out,outint) 'NLAMBDA_WGHT',nlambda_wght
      WRITE(iunit_out,outint) 'CALC_BIRTH',calc_birth
      WRITE(iunit_out,outint) 'CALC_BREMS',calc_brems
      WRITE(iunit_out,outint) 'CALC_COLD',calc_cold
      WRITE(iunit_out,outint) 'CALC_NEUTRON',calc_neutron
      WRITE(iunit_out,outint) 'CALC_DCX',calc_dcx
      WRITE(iunit_out,outint) 'CALC_HALO',calc_halo
      WRITE(iunit_out,outint) 'CALC_BES',calc_bes
      WRITE(iunit_out,outint) 'CALC_FIDA',calc_fida
      WRITE(iunit_out,outint) 'CALC_PFIDA',calc_pfida
      WRITE(iunit_out,outint) 'CALC_NPA',calc_npa
      WRITE(iunit_out,outint) 'CALC_FIDA_WGHT',calc_fida_wght
      WRITE(iunit_out,outint) 'CALC_NPA_WGHT',calc_npa_wght
      WRITE(iunit_out,outint) 'IMPURITY_CHARGE',impurity_charge
      WRITE(iunit_out,'(A)') '!---------- FLOAT VARIABLES ------------'
       WRITE(iunit_out,outflt) 'XMIN',xmin  
       WRITE(iunit_out,outflt) 'XMAX',xmax  
       WRITE(iunit_out,outflt) 'YMIN',ymin  
       WRITE(iunit_out,outflt) 'YMAX',ymax         
       WRITE(iunit_out,outflt) 'ZMIN',zmin  
       WRITE(iunit_out,outflt) 'ZMAX',zmax         
       WRITE(iunit_out,outflt) 'PINJ',pinj
       WRITE(iunit_out,outflt) 'EINJ',einj
       WRITE(iunit_out,outflt) 'AI',ai
       WRITE(iunit_out,outflt) 'ab',ab
       WRITE(iunit_out,outflt) 'ALPHA',alpha
       WRITE(iunit_out,outflt) 'BETA',beta
       WRITE(iunit_out,outflt) 'GAMMA',gamma
       WRITE(iunit_out,outflt) 'LAMBDAMIN',lambdamin
       WRITE(iunit_out,outflt) 'LAMBDAMAX',lambdamax
       WRITE(iunit_out,outflt) 'LAMBDAMIN_WGHT',lambdamin_wght
       WRITE(iunit_out,outflt) 'LAMBDAMAX_WGHT',lambdamax_wght
       WRITE(iunit_out,outflt) 'EMAX_WGHT',emax_wght
       WRITE(iunit_out,'(A)') '!---------- ARRAY VARIABLES ------------'
       WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'CURRENT_FRACTIONS',(current_fractions(n), n=1,3)
       WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'ORIGIN',(origin(n), n=1,3)
      WRITE(iunit_out,'(A)') '/'

END SUBROUTINE write_fidasim_namelist



SUBROUTINE write_fidasim_geometry

   CALL open_hdf5('fidasim_'//TRIM(id_string)//'_geometry.h5',fid,ier,LCREATE=.true.)
   IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'fidasim_'//TRIM(id_string)//'_geometry.h5',ier)
   CALL h5gopen_f(fid,'/', qid_gid, ier)
   CALL write_att_hdf5(qid_gid,'description','Geometric quantities for FIDASIM',ier)
   CALL h5gclose_f(qid_gid, ier)

   !--------------------------------------------------------------
   !           NBI
   !--------------------------------------------------------------

   CALL h5gcreate_f(fid,'nbi', qid_gid, ier)
   CALL write_att_hdf5(qid_gid,'coordinate_system','Right-handed cartesian',ier)
   CALL write_att_hdf5(qid_gid,'description','Neutral Beam Geometry',ier)
   CALL write_var_hdf5(qid_gid,'data_source',ier,STRVAR=TRIM(nbi_data_source))
   CALL write_var_hdf5(qid_gid,'name',ier,STRVAR=TRIM(name))

   CALL write_var_hdf5(qid_gid,'shape',ier,INTVAR=shape)
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'shape',ier)
   CALL h5dopen_f(qid_gid, 'shape', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','-',ier)
   CALL write_att_hdf5(temp_gid,'description','Shape of the beam source grid (1 or 2)',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'src',3,ier,DBLVAR=DBLE(src))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'src',ier)
   CALL h5dopen_f(qid_gid, 'src', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Position of the source grid in machine coordinates',ier)
   CALL h5dclose_f(temp_gid,ier)


   CALL write_var_hdf5(qid_gid,'axis',3,ier,DBLVAR=DBLE(axis_nbi))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'axis',ier)
   CALL h5dopen_f(qid_gid, 'axis', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','-',ier)
   CALL write_att_hdf5(temp_gid,'description','Axis of the beam centerline: Centerline(t) = src + axis*t ',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'widy',ier,DBLVAR=DBLE(widy))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'widy',ier)
   CALL h5dopen_f(qid_gid, 'widy', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Source grid half-width in the horizontal direction',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'widz',ier,DBLVAR=DBLE(widz))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'widz',ier)
   CALL h5dopen_f(qid_gid, 'widz', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Source grid half-height in the vertical direction',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'divy',3,ier,DBLVAR=DBLE(divy))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'divy',ier)
   CALL h5dopen_f(qid_gid, 'divy', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','rad',ier)
   CALL write_att_hdf5(temp_gid,'description','Horizontal beam divergence',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'divz',3,ier,DBLVAR=DBLE(divz))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'divz',ier)
   CALL h5dopen_f(qid_gid, 'divz', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','rad',ier)
   CALL write_att_hdf5(temp_gid,'description','Vertical beam divergence',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'focy',ier,DBLVAR=DBLE(focy))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'focy',ier)
   CALL h5dopen_f(qid_gid, 'focy', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Horizontal focal length',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'focz',ier,DBLVAR=DBLE(focz))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'focz',ier)
   CALL h5dopen_f(qid_gid, 'focz', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Vertical focal length',ier)
   CALL h5dclose_f(temp_gid,ier)


   CALL write_var_hdf5(qid_gid,'naperture',ier,INTVAR=naperture)
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'naperture',ier)
   CALL h5dopen_f(qid_gid, 'naperture', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','-',ier)
   CALL write_att_hdf5(temp_gid,'description','Number of apertures',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'ashape',naperture,ier,INTVAR=ashape)
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ashape',ier)
   CALL h5dopen_f(qid_gid, 'ashape', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','-',ier)
   CALL write_att_hdf5(temp_gid,'description','Number of apertures',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'awidy',naperture,ier,DBLVAR=awidy)
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'awidy',ier)
   CALL h5dopen_f(qid_gid, 'awidy', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Half-width of the aperture(s)',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'awidz',naperture,ier,DBLVAR=awidz)
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'awidz',ier)
   CALL h5dopen_f(qid_gid, 'awidz', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Half-width of the aperture(s)',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'aoffy',naperture,ier,DBLVAR=aoffy)
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'aoffy',ier)
   CALL h5dopen_f(qid_gid, 'aoffy', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Horizontal (y) offset of the aperture(s) relative to the +x aligned beam centerline',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'aoffz',naperture,ier,DBLVAR=aoffz)
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'aoffz',ier)
   CALL h5dopen_f(qid_gid, 'aoffz', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Vertical (z) offset of the aperture(s) relative to the +x aligned beam centerline',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'adist',naperture,ier,DBLVAR=adist)
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'adist',ier)
   CALL h5dopen_f(qid_gid, 'adist', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Distance from the center of the beam source grid to the aperture(s) plane',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL h5gclose_f(qid_gid, ier)

   !--------------------------------------------------------------
   !           SPEC
   !--------------------------------------------------------------
   CALL h5gcreate_f(fid,'spec', qid_gid, ier)
   !CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
   CALL write_att_hdf5(qid_gid,'coordinate_system','Right-handed cartesian',ier)
   CALL write_att_hdf5(qid_gid,'description','FIDA/BES Chord Geometry',ier)


   CALL write_var_hdf5(qid_gid,'data_source',ier,STRVAR=TRIM(spec_data_source))
   CALL write_var_hdf5(qid_gid,'system',ier,STRVAR=TRIM(system))

   CALL write_var_hdf5(qid_gid,'adist',naperture,ier,DBLVAR=adist)
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'adist',ier)
   CALL h5dopen_f(qid_gid, 'adist', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Distance from the center of the beam source grid to the aperture(s) plane',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'nchan',ier,INTVAR=nchan)
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nchan',ier)
   CALL h5dopen_f(qid_gid, 'nchan', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','-',ier)
   CALL write_att_hdf5(temp_gid,'description','Number of channels',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'radius',nchan,ier,DBLVAR=radius(1:nchan))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'radius',ier)
   CALL h5dopen_f(qid_gid, 'radius', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Line of sight radius at midplane or tangency point',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'lens',3,nchan,ier,DBLVAR=lens(:,1:nchan))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lens',ier)
   CALL h5dopen_f(qid_gid, 'lens', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Lens location in machine coordinates',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'axis',3,nchan,ier,DBLVAR=axis_spec(:,1:nchan))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'axis',ier)
   CALL h5dopen_f(qid_gid, 'axis', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','-',ier)
   CALL write_att_hdf5(temp_gid,'description','Optical axis/direction of the lines of sight',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'spot_size',nchan,ier,DBLVAR=spot_size(1:nchan))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'spot_size',ier)
   CALL h5dopen_f(qid_gid, 'spot_size', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','cm',ier)
   CALL write_att_hdf5(temp_gid,'description','Radius of the collecting volume',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'sigma_pi',nchan,ier,DBLVAR=sigma_pi(1:nchan))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'sigma_pi',ier)
   CALL h5dopen_f(qid_gid, 'sigma_pi', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','-',ier)
   CALL write_att_hdf5(temp_gid,'description','Ratio of the intensities of the sigma and pi stark lines',ier)
   CALL h5dclose_f(temp_gid,ier)

   CALL write_var_hdf5(qid_gid,'id',nchan,ier,STRVAR=id(1:nchan))
   IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'id',ier)
   CALL h5dopen_f(qid_gid, 'id', temp_gid, ier)
   CALL write_att_hdf5(temp_gid,'units','-',ier)
   CALL write_att_hdf5(temp_gid,'description','Channel ID',ier)
   CALL h5dclose_f(temp_gid,ier)


   !Close file
   CALL h5gclose_f(qid_gid, ier)
   CALL close_hdf5(fid,ier)
   IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'fidasim_'//TRIM(id_string)//'_geometry.h5',ier)

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
END SUBROUTINE write_fidasim_geometry


SUBROUTINE write_fidasim_equilibrium

        !INTEGER :: ier
        INTEGER :: i, j, k, l, m, n
        !INTEGER(HID_T) ::  qid_gid, temp_gid
        INTEGER(HID_T) ::  qid_gid2
        INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: mask

        REAL*8 :: fvalE(1,3), fval(1), xparam, yparam, zparam

        DOUBLE PRECISION         :: x0, y0, z0, vol
        DOUBLE PRECISION, ALLOCATABLE :: rtemp(:,:,:), rtemp2(:,:,:), rtemp3(:,:,:), rtemp4(:,:,:), rtemp5(:,:,:), r1dtemp(:), r2dtemp(:,:), r4dtemp(:,:,:,:)

        CHARACTER(LEN=8) :: temp_str8

        INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/), ictE(8)=(/0,1,1,1,0,0,0,0/)
        DOUBLE PRECISION, PARAMETER :: e_charge      = 1.602176565e-19 !e_c
        DOUBLE PRECISION, PARAMETER :: zero          = 0.0D0 ! 0.0
        DOUBLE PRECISION, PARAMETER :: t_min          = 1.0D-3 !        
        !--------------------------------------------------------------
        !           EQUILIBRIUM - FIELDS
        !--------------------------------------------------------------

        CALL open_hdf5('fidasim_'//TRIM(id_string)//'_equilibrium.h5',fid,ier,LCREATE=.true.)
        IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'fidasim_'//TRIM(id_string)//'_equilibrium.h5',ier)


        CALL h5gcreate_f(fid,'fields', qid_gid, ier)
        !CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
        CALL write_att_hdf5(qid_gid,'data_source','Data initialized from BEAMS3D',ier)

        CALL write_fidasim_grid(qid_gid)

        !           MASK
        ALLOCATE(mask(nr_fida,nz_fida,nphi_fida))
        mask = 1
        CALL write_var_hdf5(qid_gid,'mask',nr_fida,nz_fida,nphi_fida,ier,INTVAR=mask) !PLACEHOLDER, should be "Boolean mask that indicates where the fields are well defined", Dim: [nr,nz,nphi]
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'mask',ier)
        CALL h5dopen_f(qid_gid, 'mask', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'description','Boolean mask that indicates where the fields are well defined',ier)
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'mask',ier)
        CALL h5dclose_f(temp_gid,ier)
        DEALLOCATE(mask)

        ALLOCATE(rtemp(nr_fida,nz_fida, nphi_fida))
        ALLOCATE(rtemp2(nr_fida,nz_fida, nphi_fida))
        ALLOCATE(rtemp3(nr_fida,nz_fida, nphi_fida))
        DO l = 1,nr_fida
        DO n = 1,nz_fida
            DO m = 1,nphi_fida
                ! Eval Spline
                x0 = MOD(phiaxis_fida(m), phimax)
                IF (x0 < 0) x0 = x0 + phimax
                i = MIN(MAX(COUNT(raxis < raxis_fida(l)),1),nr-1)
                j = MIN(MAX(COUNT(phiaxis < x0),1),nphi-1)
                k = MIN(MAX(COUNT(zaxis < zaxis_fida(n)),1),nz-1)
                xparam = (raxis_fida(l) - raxis(i)) * hri(i)
                yparam = (x0 - phiaxis(j)) * hpi(j)
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
        ALLOCATE(r4dtemp(nr_fida,nz_fida,nphi_fida,3))
        ALLOCATE(r1dtemp(3))
        r1dtemp = 1
        DO l = 1,nr_fida
            DO n = 1,nz_fida
                DO m = 1,nphi_fida
                    ! Eval Spline
                    ! Get the gridpoint info (this is possible since all grids are the same)
                    x0 = MOD(phiaxis_fida(m), phimax)
                    IF (x0 < 0) x0 = x0 + phimax
                    i = MIN(MAX(COUNT(raxis < raxis_fida(l)),1),nr-1)
                    j = MIN(MAX(COUNT(phiaxis < x0),1),nphi-1)
                    k = MIN(MAX(COUNT(zaxis < zaxis_fida(n)),1),nz-1)
                    xparam = (raxis_fida(l) - raxis(i)) * hri(i)
                    yparam = (x0 - phiaxis(j)) * hpi(j)
                    zparam = (zaxis_fida(n) - zaxis(k)) * hzi(k)
                    ! Evaluate the Splines
                    CALL R8HERM3FCN(ictE,1,1,fvalE,i,j,k,xparam,yparam,zparam,& !evaluate at grid points
                    hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                    POT4D(1,1,1,1),nr,nphi,nz)
                    r1dtemp(1:3) =-fvalE(1,1:3)
                    r1dtemp(2)=r1dtemp(2)/(raxis_fida(l)+1/r_h/2)
                    r4dtemp(l,n,m,1:3) = r1dtemp(1:3)
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

        CALL write_fidasim_grid(qid_gid)

        !           MASK
        ALLOCATE(mask(nr_fida,nz_fida,nphi_fida))
        mask = 1
        CALL write_var_hdf5(qid_gid,'mask',nr_fida,nz_fida,nphi_fida,ier,INTVAR=mask) !PLACEHOLDER, should be "Boolean mask that indicates where the fields are well defined", Dim: [nr,nz,nphi]
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
        CALL write_var_hdf5(qid_gid,'vz',nr_fida,nz_fida, nphi_fida,ier,DBLVAR=rtemp)
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'vz',ier)
        CALL h5dopen_f(qid_gid, 'vz', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'description','Bulk plasma flow in the z-direction: Vz(r,z,phi)',ier)
        CALL write_att_hdf5(temp_gid,'units','cm/s',ier)
        CALL h5dclose_f(temp_gid,ier)


        ALLOCATE(rtemp2(nr_fida,nz_fida, nphi_fida))
        ALLOCATE(rtemp3(nr_fida,nz_fida, nphi_fida))
        ALLOCATE(rtemp4(nr_fida,nz_fida, nphi_fida))
        IF (nvtor>0) THEN
            ALLOCATE(rtemp5(nr_fida,nz_fida, nphi_fida))
        ELSE
            CALL write_var_hdf5(qid_gid,'vt',nr_fida,nz_fida, nphi_fida,ier,DBLVAR=rtemp)
            IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'vt',ier)
            CALL h5dopen_f(qid_gid, 'vt', temp_gid, ier)
            CALL write_att_hdf5(temp_gid,'description','Bulk plasma flow in the toroidal phi-direction: Vphi(r,z,phi)',ier)
            CALL write_att_hdf5(temp_gid,'units','cm/s',ier)
            CALL h5dclose_f(temp_gid,ier)         
        END IF

        DO l = 1,nr_fida
        DO n = 1,nz_fida
            DO m = 1,nphi_fida
                ! Eval Spline
                ! Get the gridpoint info (this is possible since all grids are the same)
                x0 = MOD(phiaxis_fida(m), phimax)
                IF (x0 < 0) x0 = x0 + phimax
                i = MIN(MAX(COUNT(raxis < raxis_fida(l)),1),nr-1)
                j = MIN(MAX(COUNT(phiaxis < x0),1),nphi-1)
                k = MIN(MAX(COUNT(zaxis < zaxis_fida(n)),1),nz-1)
                xparam = (raxis_fida(l) - raxis(i)) * hri(i)
                yparam = (x0 - phiaxis(j)) * hpi(j)
                zparam = (zaxis_fida(n) - zaxis(k)) * hzi(k)
                ! Evaluate the Splines
                CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                    hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                    TE4D(1,1,1,1),nr,nphi,nz)
                rtemp(l,n,m) = max(fval(1),t_min)
                CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                    hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                    NE4D(1,1,1,1),nr,nphi,nz)
                rtemp2(l,n,m) = max(fval(1),zero)
                CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                    hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                    TI4D(1,1,1,1),nr,nphi,nz)
                rtemp3(l,n,m) = max(fval(1),t_min)
                CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                    hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                    ZEFF4D(1,1,1,1),nr,nphi,nz)
                rtemp4(l,n,m) = max(fval(1),one)
                IF (nvtor>0) THEN
                  CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                     hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                     VTOR4D(1,1,1,1),nr,nphi,nz)
                  rtemp5(l,n,m) = max(fval(1),zero)        
                END IF
            END DO
        END DO
        END DO


        IF (nvtor>0) THEN
         CALL write_var_hdf5(qid_gid,'vt',nr_fida,nz_fida, nphi_fida,ier,DBLVAR=rtemp5)
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'vt',ier)
         CALL h5dopen_f(qid_gid, 'vt', temp_gid, ier)
         CALL write_att_hdf5(temp_gid,'description','Bulk plasma flow in the toroidal phi-direction: Vphi(r,z,phi)',ier)
         CALL write_att_hdf5(temp_gid,'units','cm/s',ier)
         CALL h5dclose_f(temp_gid,ier)
        END IF
        

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
        IF (nvtor>0) DEALLOCATE(rtemp5)

        !--------------------------------------------------------------
        !           Profiles
        !--------------------------------------------------------------
        CALL h5gcreate_f(qid_gid,'profiles', qid_gid2, ier)
        CALL write_att_hdf5(qid_gid2,'data_source','Data initialized from BEAMS3D ',ier)
        CALL write_att_hdf5(qid_gid2,'description','no bulk plasma rotation/flow',ier)

        CALL write_var_hdf5(qid_gid2,'rho',nne, ier,DBLVAR=DBLE(SQRT(NE_AUX_S(1:nne))))
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'rho',ier)
        CALL h5dopen_f(qid_gid2, 'rho', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','-',ier)
        CALL write_att_hdf5(temp_gid,'description','sqrt(s), nne',ier)
        CALL h5dclose_f(temp_gid,ier)


        CALL write_var_hdf5(qid_gid2,'dene',nne, ier,DBLVAR=DBLE(NE_AUX_F(1:nne)*1.0E-6))
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'dene',ier)
        CALL h5dopen_f(qid_gid2, 'dene', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','[m^-3]',ier)
        CALL write_att_hdf5(temp_gid,'description','Electron Density',ier)
        CALL h5dclose_f(temp_gid,ier)

        CALL write_var_hdf5(qid_gid2,'te',nte, ier,DBLVAR=DBLE(TE_AUX_F(1:nte)*1.0E-3))
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'te',ier)
        CALL h5dopen_f(qid_gid2, 'te', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','[eV]',ier)
        CALL write_att_hdf5(temp_gid,'description','Electron Temperature',ier)
        CALL h5dclose_f(temp_gid,ier)


        CALL write_var_hdf5(qid_gid2,'ti',nti, ier,DBLVAR=DBLE(TI_AUX_F(1:nti)*1.0E-3))
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ti',ier)
        CALL h5dopen_f(qid_gid2, 'ti', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','[eV]',ier)
        CALL write_att_hdf5(temp_gid,'description','Ion Temperature',ier)
        CALL h5dclose_f(temp_gid,ier)


        CALL write_var_hdf5(qid_gid2,'zeff',nzeff, ier,DBLVAR=DBLE(ZEFF_AUX_F(1:nzeff)))
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'zeff',ier)
        CALL h5dopen_f(qid_gid2, 'zeff', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','[-]',ier)
        CALL write_att_hdf5(temp_gid,'description','Effective Charge',ier)
        CALL h5dclose_f(temp_gid,ier)

        CALL h5gclose_f(qid_gid2, ier)

        CALL h5gclose_f(qid_gid, ier)
        ! Close file
        CALL close_hdf5(fid,ier)
        IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'fidasim_'//TRIM(id_string)//'_equilibrium.h5',ier)
END SUBROUTINE write_fidasim_equilibrium


SUBROUTINE write_fidasim_grid(qid_gid_in)
         INTEGER(HID_T), INTENT(IN)::  qid_gid_in
       !  INTEGER :: ier
      !   INTEGER(HID_T) ::  temp_gid

        DOUBLE PRECISION, ALLOCATABLE :: r2dtemp(:,:)

               !           GRID
        CALL write_var_hdf5(qid_gid_in,'nr',ier,INTVAR=nr_fida)
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nr',ier)
        CALL h5dopen_f(qid_gid_in, 'nr', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','-',ier)
        CALL write_att_hdf5(temp_gid,'description','Number of R values',ier)
        CALL h5dclose_f(temp_gid,ier)
        CALL write_var_hdf5(qid_gid_in,'nphi',ier,INTVAR=nphi_fida)
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nphi',ier)
        CALL h5dopen_f(qid_gid_in, 'nphi', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','-',ier)
        CALL write_att_hdf5(temp_gid,'description','Number of Phi values',ier)
        CALL h5dclose_f(temp_gid,ier)
        CALL write_var_hdf5(qid_gid_in,'nz',ier,INTVAR=nz_fida)
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nz',ier)
        CALL h5dopen_f(qid_gid_in, 'nz', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','-',ier)
        CALL write_att_hdf5(temp_gid,'description','Number of Z values',ier)
        CALL h5dclose_f(temp_gid,ier)
        CALL write_var_hdf5(qid_gid_in,'time',ier,DBLVAR=DBLE(0)) !!!!Assumes steady state/dummy
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'time',ier)
        CALL h5dopen_f(qid_gid_in, 'time', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','s',ier)
        CALL write_att_hdf5(temp_gid,'description','Distribution time',ier)
        CALL h5dclose_f(temp_gid,ier)

        CALL write_var_hdf5(qid_gid_in,'r',nr_fida, ier,DBLVAR=DBLE(raxis_fida*100)) !convert from m to cm
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'r',ier)
        CALL h5dopen_f(qid_gid_in, 'r', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','cm',ier)
        CALL write_att_hdf5(temp_gid,'description','Radius',ier)
        CALL h5dclose_f(temp_gid,ier)
        CALL write_var_hdf5(qid_gid_in,'phi',nphi_fida, ier,DBLVAR=phiaxis_fida)
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'phi',ier)
        CALL h5dopen_f(qid_gid_in, 'phi', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','rad',ier)
        CALL write_att_hdf5(temp_gid,'description','Toroidal angle',ier)
        CALL h5dclose_f(temp_gid,ier)
        CALL write_var_hdf5(qid_gid_in,'z',nz_fida, ier,DBLVAR=DBLE(zaxis_fida*100)) !convert from m to cm
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'z',ier)
        CALL h5dopen_f(qid_gid_in, 'z', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','cm',ier)
        CALL write_att_hdf5(temp_gid,'description','Z',ier)
        CALL h5dclose_f(temp_gid,ier)


        ALLOCATE(r2dtemp(nr_fida,nz_fida))
        r2dtemp = SPREAD(raxis_fida, 2, nz_fida)
        CALL write_var_hdf5(qid_gid_in,'r2d',nr_fida, nz_fida, ier,DBLVAR=DBLE(r2dtemp*100))
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'r2d',ier)
        CALL h5dopen_f(qid_gid_in, 'r2d', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','cm',ier)
        CALL write_att_hdf5(temp_gid,'description','Radius grid: R(r,z)',ier)
        CALL h5dclose_f(temp_gid,ier)
        r2dtemp = SPREAD(zaxis_fida, 1, nr_fida)
        CALL write_var_hdf5(qid_gid_in,'z2d',nr_fida, nz_fida, ier,DBLVAR=DBLE(r2dtemp*100))
        IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'z2d',ier)
        CALL h5dopen_f(qid_gid_in, 'z2d', temp_gid, ier)
        CALL write_att_hdf5(temp_gid,'units','cm',ier)
        CALL write_att_hdf5(temp_gid,'description','Z grid: Z(r,z)',ier)
        CALL h5dclose_f(temp_gid,ier)
        DEALLOCATE(r2dtemp)
END SUBROUTINE write_fidasim_grid

END MODULE fidasim_input_mod