!-----------------------------------------------------------------------
!     Module:        beams3d_follow_fo
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          10/19/2021
!     Description:   This subroutine follows the particles through
!                    the grid.  The general ODE which must be solved
!                    can be written:
!                        dX        
!                       ----   = fgc(t,q,qdot)
!                        dt        
!                    where X=(R,phi,Z).  This allows particle trajectories
!                    to be parametrized in terms of time.
!                    
!-----------------------------------------------------------------------
SUBROUTINE beams3d_follow_fo
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE beams3d_runtime
    USE beams3d_lines
    USE beams3d_grid, ONLY: tmin, tmax, delta_t, BR_spl, BZ_spl, BPHI_spl, &
                            MODB_spl, S_spl, U_spl, TE_spl, NE_spl, TI_spl, &
                            TE_spl, TI_spl, wall_load, wall_shine, &
                            plasma_mass, plasma_Zavg, plasma_Zmean, therm_factor
    USE mpi_params ! MPI
    USE beams3d_physics_mod
    USE beams3d_write_par
    USE safe_open_mod, ONLY: safe_open
    USE wall_mod, ONLY: wall_free, ihit_array, nface
    USE mpi_inc
    !-----------------------------------------------------------------------
    !     Local Variables
    !          status       MPI stats indicator
    !          mystart      Starting fieldline for each processor
    !          mypace       Number of fieldlines each processor should follow
    !          i            Dummy index
    !          j            Dummy index
    !          sender       Dummy index
    !          ier          Error flag
    !          l            Dummy index
    !          neqs_nag     Number of ODEs to solve (not limited to NAG routines)
    !          l2           Index helper
    !          itol         LSODE tolerance flag
    !          itask        LSODE flag (overshoot and interpolate)
    !          istate       LSODE restart flag
    !-----------------------------------------------------------------------
    IMPLICIT NONE
#if defined(MPI_OPT)
    INTEGER :: status(MPI_STATUS_size) !mpi stuff
    INTEGER :: mystart, mypace, i, j, sender
    INTEGER,ALLOCATABLE :: mnum(:), moffsets(:)
    INTEGER :: MPI_COMM_LOCAL
#endif
    INTEGER :: ier, l, neqs_nag, l2, itol, itask, &
               istate, iopt, lrw, liw, mf, out, iunit
    INTEGER, ALLOCATABLE :: iwork(:), itemp(:,:)
    REAL :: dist
    REAL(rprec) :: tf_max, vel_max, dt_out
    DOUBLE PRECISION, ALLOCATABLE :: w(:), q(:), t_last(:)
    DOUBLE PRECISION :: tf_nag, eps_temp, t_nag, t1_nag, &
                        tol_nag, rtol
    DOUBLE PRECISION :: atol(4), rwork(84)
    DOUBLE PRECISION :: rkh_work(4, 2)
    DOUBLE PRECISION :: qdot1(4)
    CHARACTER*1 :: relab

    DOUBLE PRECISION, PARAMETER :: electron_mass = 9.10938356D-31 !m_e
    DOUBLE PRECISION, PARAMETER :: sqrt_pi       = 1.7724538509   !pi^(1/2)
    DOUBLE PRECISION, PARAMETER :: e_charge      = 1.60217662E-19 !e_c

    !-----------------------------------------------------------------------
    !     External Functions
    !          fgc_nag            RHS of ODE integrator (for NAG)    for BEAMS3D
    !          D02CJF               NAG ODE Solver
    !          D02CJW               NAG Dummy function
    !          jacobian_lsode       Jacobian function (for LSODE, not currently utilized)
    !-----------------------------------------------------------------------
    EXTERNAL D02CJF, D02CJW, fgc_nag, D02CJX, out_beams3d_nag
    EXTERNAL fgc_lsode, jacobian_lsode
    EXTERNAL fgc_rkh68
    !-----------------------------------------------------------------------
    !     Begin Subroutine
    !-----------------------------------------------------------------------
    ! Initializations
    ier = 0
    tol_nag = follow_tol
    neqs_nag = 6
    relab = "M"
    mf = 10

    ! Calc max time to follow particles
    i = MAXLOC(ABS(t_end),1)
    tf_max = t_end(i)

    ! Calculate timestep for integration (here vll is Vtotal)
    vel_max = MAX(MAXVAL(ABS(vll_start)),1E6)
    dt = SIGN(MAX(lendt_m/vel_max,1D-9),tf_max)

    ! Calculate number of integration timesteps per output timestep
    ndt_max = MAX(CEILING(tf_max/(dt*NPOINC)),1)

    ! Adjust dt to match ndt_max
    nsteps = ndt_max * NPOINC
    dt = tf_max/(ndt_max*NPOINC)
    dt_out = tf_max/NPOINC
    
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
    IF (ALLOCATED(w)) DEALLOCATE(w)
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
    IF (ALLOCATED(neut_lines)) DEALLOCATE(neut_lines)
    
    ! Output some stuff
    IF (lverb) THEN
       WRITE(6, '(A)') '----- FOLLOWING PARTICLE TRAJECTORIES -----'
       WRITE(6, '(A,A)')          '      Method: ', TRIM(int_type)
       WRITE(6, '(A,I9)')         '   Particles: ', nparticles
       WRITE(6, '(A,I9,A,EN12.3)') '       Steps: ', nsteps, '   Delta-t: ', dt
       WRITE(6, '(A,I9,A,EN12.3)') '      NPOINC: ', npoinc, '    dt_out: ', dt_out
       SELECT CASE(TRIM(int_type))
          CASE("NAG")
             WRITE(6, '(A,EN12.3,A,A1)') '         Tol: ', follow_tol, '  Type: ', relab
          CASE("LSODE")
             WRITE(6, '(A,EN12.3,A,I2)') '         Tol: ', follow_tol, '  Type: ', mf
       END SELECT
       WRITE(6, '(5X,A,I3,A)', ADVANCE = 'no') 'Trajectory Calculation [', 0, ']%'
       CALL FLUSH(6)
    END IF
    ! Allocations
    ALLOCATE(q(neqs_nag), STAT = ier)
    IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'Q', ier)
    ALLOCATE(R_lines(0:npoinc, mystart:myend), Z_lines(0:npoinc, mystart:myend), &
             PHI_lines(0:npoinc, mystart:myend), &
             neut_lines(0:npoinc, mystart:myend), S_lines(0:npoinc, mystart:myend), U_lines(0:npoinc, mystart:myend), &
              B_lines(0:npoinc, mystart:myend), &
              vr_lines(0:npoinc, mystart:myend), vphi_lines(0:npoinc, mystart:myend), vz_lines(0:npoinc, mystart:myend),STAT = ier)
    IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'R_LINES, PHI_LINES, Z_LINES', ier)
    ALLOCATE(t_last(mystart:myend), STAT = ier)
    IF (ier /= 0) CALL handle_err(ALLOC_ERR, 't_last', ier)

    ! Initializations
    R_lines = 0.0; Z_lines = 0.0; PHI_lines = -1.0
    vr_lines = 0.0; vphi_lines = 0.0; vz_lines = 0.0
    S_lines = 1.5; U_lines = 0.0; B_lines = -1.0
    t_last = 0.0
    R_lines(0, mystart:myend) = R_start(mystart:myend)
    Z_lines(0, mystart:myend) = Z_start(mystart:myend)
    PHI_lines(0, mystart:myend) = phi_start(mystart:myend)
    vr_lines(0, mystart:myend) = vr_start(mystart:myend)
    vphi_lines(0, mystart:myend) = vphi_start(mystart:myend)
    vz_lines(0, mystart:myend) = vz_start(mystart:myend)
    neut_lines(0, mystart:myend) = .FALSE.
    IF (lbeam) neut_lines(0, mystart:myend) = .TRUE.

    ! Some helpers
    fact_vsound = 1.5*sqrt(e_charge/plasma_mass)*therm_factor
    fact_crit = SQRT(2*e_charge/plasma_mass)*(0.75*sqrt_pi*sqrt(plasma_mass/electron_mass))**(1.0/3.0) ! Wesson pg 226 5.4.9
    fact_kick = pi2*2*SQRT(pi*1E-7*plasma_mass)*E_kick*freq_kick

    ! Follow Trajectories
    IF (mystart <= nparticles) THEN
        SELECT CASE (TRIM(int_type))
            CASE ("NAG")
#if defined(NAG)
                ALLOCATE(w(neqs_nag * 21 + 28), STAT = ier)
                IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'W', ier)
                DO l = mystart, myend
                    ltherm = .false.
                    lneut  = .false.
                    ndt    = ndt_max-1 ! we only record first point
                    q(1) = R_start(l)
                    q(2) = phi_start(l)
                    q(3) = Z_start(l)
                    q(4) = vr_start(l)
                    q(5) = vphi_start(l)
                    q(6) = vz_start(l)
                    t_nag = 0.0
                    tf_nag = 0.0
                    mycharge = charge(l)
                    myZ = Zatom(l)
                    mymass = mass(l)
                    mybeam = Beam(l)
                    my_end = t_end(l)
                    fact_pa   = plasma_mass/(mymass*plasma_Zmean)
                    fact_coul = myZ*(mymass+plasma_mass)/(mymass*plasma_mass*6.02214076208E+26)
                    myv_neut(1) = vr_start(l)*cos(phi_start(l)) - vphi_start(l)*sin(phi_start(l))
                    myv_neut(2) = vr_start(l)*sin(phi_start(l)) + vphi_start(l)*cos(phi_start(l))
                    myv_neut(3) = vz_start(l)
                    IF (lbeam) lneut = .TRUE.
                    CALL out_beams3d_nag(tf_nag,q)
                    IF (lbeam) THEN
                       lcollision = .FALSE.
                       ! Follow into plasma
                       CALL beams3d_follow_neut(t_nag,q)
                       mytdex = 1; ndt =1
                       tf_nag = t_nag
                       CALL out_beams3d_nag(tf_nag,q)
                       IF (tf_nag > t_end(l)) CYCLE  ! Detect end shinethrough particle
                       ! Ionize
                       CALL beams3d_ionize(tf_nag,q)
                       mytdex = 2; ndt =1
                       CALL out_beams3d_nag(tf_nag,q)
                       ltherm = .FALSE.
                       lcollision = .TRUE.
                       mytdex = 3; ndt =1
                       tf_nag = tf_nag - dt  ! Because out advances t by dt
                       ! Adjust timestep timestep
                       !CALL beams3d_calc_dt(q,moment,mymass,dt)
                    END IF
                    IF (ldepo) CYCLE
                    DO ! Must do it this way becasue lbeam changes q(4) values
                       CALL D02CJF(t_nag,tf_nag,neqs_nag,q,fpart_nag,tol_nag,relab,out_beams3d_nag,D02CJW,w,ier)
                       IF (ier < 0) CALL handle_err(D02CJF_ERR, 'beams3d_follow', ier)
                       t_last(l) = tf_nag ! Save the value here in case out_beams3d changes it
                       CALL out_beams3d_nag(tf_nag,q)
                       IF (ABS(tf_nag) > ABS(my_end)) EXIT
                    END DO
                END DO
#else
                ier = -1
                CALL handle_err(NAG_ERR, 'beams3d_follow', ier)
#endif
            CASE ("RKH68")
                ier = 0
                DO l = mystart, myend
                    myline = l
                    mytdex = 0
                    ndt    = ndt_max-1 ! we only record first point
                    ltherm = .false.
                    lneut  = .false.
                    q(1) = R_start(l)
                    q(2) = phi_start(l)
                    q(3) = Z_start(l)
                    q(4) = vr_start(l)
                    q(5) = vphi_start(l)
                    q(6) = vz_start(l)
                    xlast = q(1)*cos(q(2))
                    ylast = q(1)*sin(q(2))
                    zlast = q(3)
                    t_nag = 0.0
                    tf_nag = 0.0
                    mycharge = charge(l)
                    myZ = Zatom(l)
                    mymass = mass(l)
                    mybeam = Beam(l)
                    my_end = t_end(l)
                    fact_pa   = plasma_mass/(mymass*plasma_Zmean)
                    fact_coul = myZ*(mymass+plasma_mass)/(mymass*plasma_mass*6.02214076208E+26)
                    myv_neut(1) = vr_start(l)*cos(phi_start(l)) - vphi_start(l)*sin(phi_start(l))
                    myv_neut(2) = vr_start(l)*sin(phi_start(l)) + vphi_start(l)*cos(phi_start(l))
                    myv_neut(3) = vz_start(l)
                    IF (lbeam) lneut = .TRUE.
                    CALL out_beams3d_nag(tf_nag,q)
                    IF (lbeam) THEN
                       lcollision = .FALSE.
                       ! Follow into plasma
                       CALL beams3d_follow_neut(t_nag,q)
                       mytdex = 1; ndt =1
                       tf_nag = t_nag
                       CALL out_beams3d_nag(tf_nag,q)
                       IF (tf_nag > t_end(l)) CYCLE  ! Detect end shinethrough particle
                       ! Ionize
                       CALL beams3d_ionize(tf_nag,q)
                       mytdex = 2; ndt =1
                       CALL out_beams3d_nag(tf_nag,q)
                       ltherm = .FALSE.
                       lcollision = .TRUE.
                       mytdex = 3; ndt =1
                       tf_nag = tf_nag - dt  ! Because out advances t by dt
                       ! Adjust timestep timestep
                       !CALL beams3d_calc_dt(q,moment,mymass,dt)
                    END IF
                    IF (ldepo) CYCLE
                    DO
                        CALL drkhvg(t_nag, q, neqs_nag, dt, 2, fpart_rkh68, rkh_work, iopt, ier)
                        IF (ier < 0) CALL handle_err(RKH68_ERR, 'beams3d_follow', ier)
                        q(1)=rkh_work(1,2)
                        q(2)=rkh_work(2,2)
                        q(3)=rkh_work(3,2)
                        q(4)=rkh_work(4,2)
                        t_nag = t_nag+dt
                        tf_nag = tf_nag+dt
                        t_last(l) = tf_nag ! Save the value here in case out_beams3d changes it
                        CALL out_beams3d_nag(tf_nag,q)
                        IF ((istate == -1) .or. (istate ==-2) .or. (ABS(tf_nag) > ABS(my_end)) ) EXIT
                    END DO
                END DO
            CASE ("LSODE","DLSODE")
                ! Open a file
                iunit = myworkid+400
                ier = 0
                IF (ldebug) THEN
                   !CALL safe_open(iunit,ier,'lsode_err_'//TRIM(id_string),'replace','formatted')
                   iunit = 6
                   CALL xsetun(iunit)
                   iunit = 1
                   CALL xsetf(iunit)
                ELSE
                   iunit = 0
                   CALL xsetf(iunit)
                END IF
                ! Allocate helpers
                lrw = 20 + 16 * neqs_nag
                liw = 20
                ALLOCATE(w(lrw))
                ALLOCATE(iwork(liw))
                ier = 0
                DO l = mystart, myend
                    ! Setup LSODE parameters
                    iopt = 0 ! No optional output
                    w = 0; iwork = 0; itask = 1; istate = 1;
                    itol = 2; rtol = follow_tol; atol(:) = follow_tol
                    myline = l
                    mytdex = 0
                    ndt    = ndt_max-1 ! we only record first point
                    ! Initialize the calculation
                    ltherm = .false.
                    lneut  = .false.
                    q(1) = R_start(l)
                    q(2) = phi_start(l)
                    q(3) = Z_start(l)
                    q(4) = vr_start(l)
                    q(5) = vphi_start(l)
                    q(6) = vz_start(l)
                    xlast = q(1)*cos(q(2))
                    ylast = q(1)*sin(q(2))
                    zlast = q(3)
                    t_nag = 0.0
                    tf_nag = 0.0
                    mycharge = charge(l)
                    myZ = Zatom(l)
                    mymass = mass(l)
                    mybeam = Beam(l)
                    moment = mu_start(l)
                    my_end = t_end(l)
                    fact_pa   = plasma_mass/(mymass*plasma_Zmean)
                    fact_coul = myZ*(mymass+plasma_mass)/(mymass*plasma_mass*6.02214076208E+26)
                    myv_neut(1) = vr_start(l)*cos(phi_start(l)) - vphi_start(l)*sin(phi_start(l))
                    myv_neut(2) = vr_start(l)*sin(phi_start(l)) + vphi_start(l)*cos(phi_start(l))
                    myv_neut(3) = vz_start(l)
                    ! Setup timestep
                    !CALL beams3d_calc_dt(q,moment,mymass,dt)
                    ! Begin handling particle.
                    IF (lbeam) lneut = .TRUE.
                    CALL out_beams3d_nag(tf_nag,q)
                    IF (lbeam) THEN
                       lcollision = .FALSE.
                       ! Follow into plasma
                       CALL beams3d_follow_neut(t_nag,q)
                       mytdex = 1
                       tf_nag = t_nag
                       CALL out_beams3d_nag(tf_nag,q)
                       IF (tf_nag > t_end(l)) CYCLE  ! Detect end shinethrough particle
                       ! Ionize
                       CALL beams3d_ionize(tf_nag,q)
                       mytdex = 2
                       CALL out_beams3d_nag(tf_nag,q)
                       ltherm = .FALSE.
                       lcollision = .TRUE.
                       mytdex = 3; ndt =1
                       tf_nag = tf_nag - dt  ! Because out advances t by dt
                       ! Adjust timestep timestep
                       !CALL beams3d_calc_dt(q,moment,mymass,dt)
                    END IF
                    IF (ldepo) CYCLE
                    DO
                        IF (lcollision) istate = 1
                        CALL FLUSH(6)
                        CALL DLSODE(fpart_lsode, neqs_nag, q, t_nag, tf_nag, itol, rtol, atol, itask, istate, &
                                   iopt, w, lrw, iwork, liw, jacobian_lsode, mf)
                        IF ((istate == -3) .or. (istate == -4)) THEN
                           ! BIG  DEBUG MESSAGE
                           CALL FLUSH(6)
                           WRITE(6,*) '------------------'
                           WRITE(6,*) '     ',myworkid, l, myline,mf
                           WRITE(6,*) '     ',myworkid, neqs_nag, t_nag, tf_nag
                           WRITE(6,*) '     ',myworkid, q, moment
                           WRITE(6,*) '     ',myworkid, itol, rtol, atol
                           WRITE(6,*) '     ',myworkid, itask,istate,iopt
                           WRITE(6,*) '     ',myworkid, w(5),w(6),w(7)
                           WRITE(6,*) '     ',myworkid, iwork(5),iwork(6),iwork(7)
                           WRITE(6,*) '     ',myworkid, w(11:14)
                           WRITE(6,*) '     ',myworkid, iwork(11:18)
                           WRITE(6,*) '------------------'
                           CALL FLUSH(6)
                           CALL handle_err(LSODE_ERR, 'beams3d_follow', istate)
                        END IF
                        iwork(11) = 0; iwork(12) = 0; iwork(13) = 0
                        t_last(l) = tf_nag ! Save the value here in case out_beams3d changes it
                        CALL out_beams3d_nag(tf_nag,q)
                        IF ((istate == -1) .or. (istate ==-2) .or. (ABS(tf_nag) > ABS(my_end)) ) EXIT
                    END DO
                END DO
                IF (ldebug) CLOSE(iunit)
             CASE ('DEBUG')
                DO l = 0, npoinc
                   R_lines(l,mystart:myend) = REAL(l)
                END DO
                B_lines(0:npoinc,mystart:myend) = REAL(myworkid)
        END SELECT
    END IF
      
    IF (lverb) THEN
       CALL backspace_out(6,36)
       CALL FLUSH(6)
       WRITE(6,'(36X)',ADVANCE='no')
       CALL FLUSH(6)
       CALL backspace_out(6,36)
       WRITE(6,*)
       CALL FLUSH(6)
    END IF    

    ! Check for crash
    istate = 0
    CALL handle_err(MPI_CHECK,'beams3d_follow',istate)

    !Deallocations
    IF (ALLOCATED(q)) DEALLOCATE(q)
    IF (ALLOCATED(w)) DEALLOCATE(w)
    IF (ALLOCATED(iwork)) DEALLOCATE(iwork)

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
       IF (ASSOCIATED(ihit_array)) THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,ihit_array,nface,MPI_INTEGER,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
       END IF
       IF (ASSOCIATED(wall_load)) THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,wall_load,nface*nbeams,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
       END IF
       IF (ASSOCIATED(wall_shine)) THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,wall_shine,nface*nbeams,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
       END IF
       CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
    END IF
    CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)

    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,      'R_lines', DBLVAR=R_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,    'PHI_lines', DBLVAR=PHI_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,      'Z_lines', DBLVAR=Z_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,     'vr_lines', DBLVAR=vr_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,   'vphi_lines', DBLVAR=vphi_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,     'vz_lines', DBLVAR=vz_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,      'S_lines', DBLVAR=S_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,      'U_lines', DBLVAR=U_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,      'B_lines', DBLVAR=B_lines)
    CALL beams3d_write1d_parhdf5(         1, nparticles, mystart, myend,      't_end',   DBLVAR=t_last,FILENAME='beams3d_'//TRIM(id_string))
    ALLOCATE(itemp(0:npoinc,mystart:myend))
    itemp = 0; WHERE(neut_lines) itemp=1;
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,   'neut_lines', INTVAR=itemp)
    DEALLOCATE(itemp)
    IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
    IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
    IF (ALLOCATED(t_last)) DEALLOCATE(t_last)
    CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= 0) CALL handle_err(MPI_BARRIER_ERR, 'beams3d_follow', ierr_mpi)
#endif

    RETURN
    !-----------------------------------------------------------------------
    !     End Subroutine
    !-----------------------------------------------------------------------
END SUBROUTINE beams3d_follow_fo
