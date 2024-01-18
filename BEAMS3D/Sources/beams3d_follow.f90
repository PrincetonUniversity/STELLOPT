!-----------------------------------------------------------------------
!     Module:        beams3d_follow
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          10/19/2021
!     Description:   This subroutine organizes running the particle
!                    following parts of the code.
!                    
!-----------------------------------------------------------------------
SUBROUTINE beams3d_follow
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE beams3d_runtime
    USE beams3d_lines
    USE beams3d_grid, ONLY: tmin, tmax, delta_t, BR_spl, BZ_spl, BPHI_spl, &
                            MODB_spl, S_spl, U_spl, TE_spl, NE_spl, TI_spl, &
                            TE_spl, TI_spl, wall_load, wall_shine, rho_fullorbit, &
                            plasma_mass, plasma_Zmean, therm_factor, &
                            nr_fida, nphi_fida, nz_fida, nenergy_fida, &
                            npitch_fida, BEAM_DENSITY
    USE mpi_params ! MPI
    USE beams3d_physics_mod
    USE beams3d_write_par
    USE safe_open_mod, ONLY: safe_open
    USE wall_mod, ONLY: wall_free, ihit_array, nface
    USE mpi_inc
    !-----------------------------------------------------------------------
    !     Local Variables
    !          status       MPI stats indicator
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: MPI_COMM_LOCAL
    INTEGER :: i, j, l, ier, mystart, mypace
    INTEGER, ALLOCATABLE :: mnum(:), moffsets(:)
    INTEGER, ALLOCATABLE :: itemp(:,:)
    REAL :: dist
    REAL(rprec) :: tf_max, vel_max, weight_save
    DOUBLE PRECISION :: tf_nag, t_nag
    DOUBLE PRECISION, ALLOCATABLE :: q(:)

    DOUBLE PRECISION, PARAMETER :: electron_mass = 9.10938356D-31 !m_e
    DOUBLE PRECISION, PARAMETER :: sqrt_pi       = 1.7724538509   !pi^(1/2)
    DOUBLE PRECISION, PARAMETER :: e_charge      = 1.60217662E-19 !e_c
    !-----------------------------------------------------------------------
    !     Begin Subroutine
    !-----------------------------------------------------------------------

    ! Calc max time to follow particles
    i = MAXLOC(ABS(t_end),1)
    tf_max = t_end(i)

    ! Calculate timestep for integration
    vel_max = MAX(MAXVAL(ABS(vll_start)),1E6)
    dt = SIGN(MAX(lendt_m/vel_max,1D-9),tf_max)

    ! Calculate number of integration timesteps per output timestep
    ndt_max = MAX(CEILING(tf_max/(dt*NPOINC)),1)

    ! Adjust dt to match ndt_max  (note dt doesn't really matter)
    dt = tf_max/(ndt_max*NPOINC)
    
    ! Break up the work
    CALL MPI_CALC_MYRANGE(MPI_COMM_BEAMS, 1, nparticles, mystart, myend)

    ! Save mystart and myend
    mystart_save = mystart
    myend_save = myend

    IF (lhitonly) THEN
        npoinc = 2
    END IF

#if defined(MPI_OPT)
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
      ALLOCATE(mnum(nprocs_beams), moffsets(nprocs_beams))
      CALL MPI_ALLGATHER((myend-mystart+1)*(npoinc+1),1,MPI_INTEGER,mnum,1,MPI_INTEGER,MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_ALLGATHER((mystart-1)*(npoinc+1),1,MPI_INTEGER,moffsets,1,MPI_INTEGER,MPI_COMM_BEAMS,ierr_mpi)
#endif

    ! Deallocations
    IF (ALLOCATED(q)) DEALLOCATE(q)
    IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
    IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
    IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
    IF (ALLOCATED(vll_lines)) DEALLOCATE(vll_lines)
    IF (ALLOCATED(S_lines)) DEALLOCATE(S_lines)
    IF (ALLOCATED(U_lines)) DEALLOCATE(U_lines)
    IF (ALLOCATED(B_lines)) DEALLOCATE(B_lines)
    IF (ALLOCATED(vr_lines)) DEALLOCATE(vr_lines)
    IF (ALLOCATED(vphi_lines)) DEALLOCATE(vphi_lines)
    IF (ALLOCATED(vz_lines)) DEALLOCATE(vz_lines)
    IF (ALLOCATED(moment_lines)) DEALLOCATE(moment_lines)
    IF (ALLOCATED(neut_lines)) DEALLOCATE(neut_lines)
    ! Allocations
    ALLOCATE(q(4), STAT = ier)
    IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'Q', ier)
    ALLOCATE(R_lines(0:npoinc, mystart:myend), Z_lines(0:npoinc, mystart:myend), &
             PHI_lines(0:npoinc, mystart:myend), vll_lines(0:npoinc, mystart:myend), moment_lines(0:npoinc, mystart:myend), &
             neut_lines(0:npoinc, mystart:myend), S_lines(0:npoinc, mystart:myend), U_lines(0:npoinc, mystart:myend), &
             vr_lines(0:npoinc, mystart:myend), vphi_lines(0:npoinc, mystart:myend), vz_lines(0:npoinc, mystart:myend), &
              B_lines(0:npoinc, mystart:myend), STAT = ier)
    IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'R_LINES, PHI_LINES, Z_LINES', ier)
    ALLOCATE(t_last(mystart:myend), STAT = ier)
    IF (ier /= 0) CALL handle_err(ALLOC_ERR, 't_last', ier)

    ! Set lbeam to false if doing a box simulation
    lbeam = (lbeam .and. (.not.lboxsim))

    ! Initializations
    lbeam = (lbeam .and. (.not.lboxsim))
    R_lines = 0.0; Z_lines = 0.0; PHI_lines = -1.0
    vll_lines = 0.0; moment_lines = 0.0
    S_lines = 1.5; U_lines = 0.0; B_lines = -1.0
    t_last = 0.0
    R_lines(0, mystart:myend)      = R_start(mystart:myend)
    Z_lines(0, mystart:myend)      = Z_start(mystart:myend)
    PHI_lines(0, mystart:myend)    = phi_start(mystart:myend)
    vll_lines(0, mystart:myend)    = vll_start(mystart:myend)
    moment_lines(0, mystart:myend) = mu_start(mystart:myend)
    vr_lines(0, mystart:myend)     = vr_start(mystart:myend)
    vphi_lines(0, mystart:myend)   = vphi_start(mystart:myend)
    vz_lines(0, mystart:myend)     = vz_start(mystart:myend)
    neut_lines(0, mystart:myend)   = .FALSE.
    IF (lbeam) neut_lines(0, mystart:myend) = .TRUE.

    ! Some helpers
    fact_vsound = 1.5*sqrt(e_charge/plasma_mass)*therm_factor
    fact_crit = SQRT(2*e_charge/plasma_mass)*(0.75*sqrt_pi*sqrt(plasma_mass/electron_mass))**(1.0/3.0) ! Wesson pg 226 5.4.9
    !fact_kick = pi2*2*SQRT(pi*1E-7*plasma_mass)*E_kick*freq_kick !old
    !fact_kick = 2*freq_kick*E_kick

    ! Handle the Beam defaults
    IF (lbeam) THEN
       ! Set VLL_START to Vtotal (may not need to do this)
       vll_start(mystart:myend) = sqrt(  vr_start(mystart:myend) * vr_start(mystart:myend) &
                                       + vphi_start(mystart:myend) * vphi_start(mystart:myend) &
                                       + vz_start(mystart:myend) * vz_start(mystart:myend) )
    END IF

    ! Initialize the particles
    DO i = mystart, myend
       lneut = lbeam
       ltherm = .FALSE.
       q(1) = R_start(i)
       q(2) = phi_start(i)
       q(3) = Z_start(i)
       q(4) = vll_start(i)
       xlast = q(1)*cos(q(2))
       ylast = q(1)*sin(q(2))
       zlast = q(3)
       tf_nag = 0.0
       mycharge = charge(i)
       myZ = Zatom(i)
       mymass = mass(i)
	   E_by_v=mymass*0.5d-3/e_charge
       mybeam = Beam(i)
       moment = mu_start(i)
       fact_pa   = plasma_mass/(mymass*plasma_Zmean)
       fact_coul = myZ*(mymass+plasma_mass)/(mymass*plasma_mass*6.02214076208E+26)
       ! Save the IC of the neutral
       my_end = t_end(i)
       myline = i
       mytdex = 0
       ndt_max = 1
       ndt     = 0
       CALL out_beams3d_nag(tf_nag,q)
       t_last(i) = tf_nag-dt
    END DO

    ! Beam Deposition
    IF (lbeam) THEN
      IF (lverb) WRITE(6, '(A)') '----- FOLLOWING NEUTRAL TRAJECTORIES -----'
       DO i = mystart, myend
          lneut = lbeam
          ltherm = .FALSE.
          q(1) = R_start(i)
          q(2) = phi_start(i)
          q(3) = Z_start(i)
          q(4) = vll_start(i)
          xlast = q(1)*cos(q(2))
          ylast = q(1)*sin(q(2))
          zlast = q(3)
          mycharge = charge(i)
          myZ = Zatom(i)
          mymass = mass(i)
		  E_by_v=mymass*0.5d-3/e_charge
          mybeam = Beam(i)
          moment = mu_start(i)
          my_end = t_end(i)
          myline = i
          mytdex = 1
          fact_pa   = plasma_mass/(mymass*plasma_Zmean)
          fact_coul = myZ*(mymass+plasma_mass)/(mymass*plasma_mass*6.02214076208E+26)
          ! Define neutral trajectory
          myv_neut(1) = vr_start(i)*cos(phi_start(i)) - vphi_start(i)*sin(phi_start(i))
          myv_neut(2) = vr_start(i)*sin(phi_start(i)) + vphi_start(i)*cos(phi_start(i))
          myv_neut(3) = vz_start(i)
          ! This is a trick which works because V_neut hasn't changed
          vr_lines(1:2,myline)   = vr_start(i)
          vphi_lines(1:2,myline) = vphi_start(i)
          vz_lines(1:2,myline)   = vz_start(i)
          lcollision = .FALSE.
          ! Inject Beam
          tf_nag = 0.0
          CALL beams3d_follow_neut(tf_nag,q)
          ! Save the point to index 1 with weight set to 0
          mytdex = 1; ndt=0; ndt_max = 1 ! Save the point
          t_last(i) = tf_nag ! Save timestep after follow_neut
          CALL out_beams3d_nag(tf_nag,q)
          IF (tf_nag > t_end(i)) CYCLE
          ! Step to gyrocenter
          CALL beams3d_ionize(q)
          ! Save to index 2
          ! no weight becasue it just jumped to this point
          mytdex = 2; ndt=0; ndt_max = 1 ! Save point
          weight_save = weight(myline)
          weight(myline) = 0
          CALL out_beams3d_nag(tf_nag,q)
          weight(myline) = weight_save
       END DO
    END IF

    DEALLOCATE(q)

    ! Calcualte Beam Density
#if defined(MPI_OPT)
    CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)
#endif
    IF (lbeam .and. lbeamdensity) CALL beams3d_beam_density

    ! Follow Trajectories
    IF (.not.ldepo) THEN
        CALL beams3d_follow_gc
        IF (rho_fullorbit < 100) CALL beams3d_follow_fo
    END IF

    ! Fix U_lines
    CALL beams3d_fix_poloidal

    ! First reduce the cumulative arrays over shared memory groups then allreduce between shared memeory groups
#if defined(MPI_OPT)
    IF (myid_sharmem == master) THEN
       CALL MPI_REDUCE(MPI_IN_PLACE,   end_state,     nparticles,          MPI_INTEGER, MPI_MAX, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(MPI_IN_PLACE, epower_prof, nbeams*ns_prof1, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(MPI_IN_PLACE, ipower_prof, nbeams*ns_prof1, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(MPI_IN_PLACE,   ndot_prof, nbeams*ns_prof1, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
    ELSE
       CALL MPI_REDUCE(end_state,     end_state,     nparticles,          MPI_INTEGER, MPI_MAX, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(epower_prof, epower_prof, nbeams*ns_prof1, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(ipower_prof, ipower_prof, nbeams*ns_prof1, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(ndot_prof,     ndot_prof, nbeams*ns_prof1, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
    END IF

    i = MPI_UNDEFINED
    IF (myid_sharmem == master) i = 0
    CALL MPI_COMM_SPLIT( MPI_COMM_BEAMS,i,myworkid,MPI_COMM_LOCAL,ierr_mpi)
    IF (myid_sharmem == master) THEN
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,   end_state,     nparticles,          MPI_INTEGER, MPI_MAX, MPI_COMM_LOCAL, ierr_mpi)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, epower_prof, nbeams*ns_prof1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_LOCAL, ierr_mpi)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, ipower_prof, nbeams*ns_prof1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_LOCAL, ierr_mpi)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,   ndot_prof, nbeams*ns_prof1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_LOCAL, ierr_mpi)
       ! This only works becasue of how FORTRAN orders things.
       DO l = 1, ns_prof5
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, dist5d_prof(:,:,:,:,:,l), nbeams*ns_prof1*ns_prof2*ns_prof3*ns_prof4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_LOCAL, ierr_mpi)
       END DO
       IF (lfidasim_cyl) THEN
          DO l = 1, npitch_fida
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, dist5d_fida(:,:,:,:,l), nr_fida*nphi_fida*nz_fida*nenergy_fida, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_LOCAL, ierr_mpi)
          END DO
       END IF
       IF (ASSOCIATED(ihit_array)) THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,ihit_array,nface,MPI_INTEGER,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
       END IF
       IF (ASSOCIATED(wall_load)) THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,wall_load,nface*nbeams,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
       END IF
       IF (ASSOCIATED(wall_shine)) THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,wall_shine,nface*nbeams,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
       END IF
       IF (ASSOCIATED(BEAM_DENSITY)) THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, BEAM_DENSITY, nbeams*nr*nphi*nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_LOCAL, ierr_mpi)
       END IF
       CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
    END IF
    CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)

    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart_save, myend_save,      'R_lines', DBLVAR=R_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart_save, myend_save,    'PHI_lines', DBLVAR=PHI_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart_save, myend_save,      'Z_lines', DBLVAR=Z_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart_save, myend_save,    'vll_lines', DBLVAR=vll_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart_save, myend_save, 'moment_lines', DBLVAR=moment_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart_save, myend_save,     'vr_lines', DBLVAR=vr_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart_save, myend_save,   'vphi_lines', DBLVAR=vphi_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart_save, myend_save,     'vz_lines', DBLVAR=vz_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart_save, myend_save,      'S_lines', DBLVAR=S_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart_save, myend_save,      'U_lines', DBLVAR=U_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart_save, myend_save,      'B_lines', DBLVAR=B_lines)
    CALL beams3d_write1d_parhdf5(         1, nparticles, mystart_save, myend_save,      't_end',   DBLVAR=t_last,FILENAME='beams3d_'//TRIM(id_string))
    ALLOCATE(itemp(0:npoinc,mystart_save:myend_save))
    itemp = 0; WHERE(neut_lines) itemp=1
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart_save, myend_save,   'neut_lines', INTVAR=itemp)
    DEALLOCATE(itemp)
    IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
    IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
    CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= 0) CALL handle_err(MPI_BARRIER_ERR, 'beams3d_follow', ierr_mpi)
#endif

    ! Adjust T_END back to values of T_last
    t_end(mystart_save:myend_save) = t_last(mystart_save:myend_save)
    IF (ALLOCATED(t_last)) DEALLOCATE(t_last)

    RETURN
    !-----------------------------------------------------------------------
    !     End Subroutine
    !-----------------------------------------------------------------------
END SUBROUTINE beams3d_follow
