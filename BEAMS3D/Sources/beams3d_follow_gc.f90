!-----------------------------------------------------------------------
!     Module:        beams3d_follow_gc
!     Authors:       S. Lazerson (lazerson@pppl.gov), M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          06/20/2012
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
SUBROUTINE beams3d_follow_gc
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE beams3d_runtime
    USE beams3d_lines
    USE beams3d_grid, ONLY: tmin, tmax, delta_t, BR_spl, BZ_spl, BPHI_spl, &
                            MODB_spl, S_spl, U_spl, TE_spl, NE_spl, TI_spl, &
                            TE_spl, TI_spl, wall_load, wall_shine, &
                            plasma_mass, plasma_Zmean, therm_factor, &
                            rho_fullorbit, rho_help, E_kick, freq_kick, &
                            nr_fida, nphi_fida, nz_fida, nenergy_fida, npitch_fida,raxis
    USE mpi_params ! MPI
    USE beams3d_write_par
    USE beams3d_physics_mod, ONLY: beams3d_calc_dt
    USE safe_open_mod, ONLY: safe_open
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
    INTEGER :: i, j
    INTEGER :: ier, l, neqs_nag, itol, itask, &
               istate, iopt, lrw, liw, mf, out, iunit
    INTEGER, ALLOCATABLE :: iwork(:), itemp(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: w(:), q(:)
    DOUBLE PRECISION :: tf_nag, eps_temp, t_nag, &
                        tol_nag, rtol, s_fullorbit, &
                        dtmin, dtmax
    DOUBLE PRECISION :: atol(4)
    DOUBLE PRECISION :: rkh_work(4, 2)
    CHARACTER*1 :: relab
    DOUBLE PRECISION, PARAMETER :: e_charge      = 1.60217662E-19 !e_c

    !-----------------------------------------------------------------------
    !     External Functions
    !          fgc_eom            RHS of ODE integrator (for NAG)    for BEAMS3D
    !          D02CJF               NAG ODE Solver
    !          D02CJW               NAG Dummy function
    !          jacobian_lsode       Jacobian function (for LSODE, not currently utilized)
    !-----------------------------------------------------------------------
    EXTERNAL D02CJF, D02CJW, fgc_eom, D02CJX, out_beams3d_gc
    EXTERNAL fgc_lsode, jacobian_lsode
    EXTERNAL fgc_rkh68
    !-----------------------------------------------------------------------
    !     Begin Subroutine
    !-----------------------------------------------------------------------
    ! Initializations
    s_fullorbit = SIGN(rho_fullorbit*rho_fullorbit,rho_fullorbit)
    ier = 0
    tol_nag = follow_tol
    neqs_nag = 4
    relab = "M"
    mf = 10
    ALLOCATE(q(neqs_nag))

! Screen output so we know what's happening
    IF (lverb) THEN
       ! Do a calculation of nstep and delta-t
       myline = MAXLOC(ABS(vll_start),1)
       my_end = t_end(myline)
       CALL beams3d_calc_dt(1,q(1),q(2),q(3),dtmin)
       myline = MINLOC(ABS(vll_start),1)
       my_end = t_end(myline)
       CALL beams3d_calc_dt(1,q(1),q(2),q(3),dtmax)
       ! Screen output
       WRITE(6, '(A)') '----- FOLLOWING GYROCENTER TRAJECTORIES -----'
       WRITE(6, '(A,A)')          '       Method: ', TRIM(int_type)
       WRITE(6, '(A,I9)')          '   Particles: ', nparticles
       WRITE(6, '(A,I9,2(A,EN12.3))') '       Steps: ', ndt_max*NPOINC, '   dt_min: ', dtmin,'   dt_max: ', dtmax
       WRITE(6, '(A,I9)')          '      NPOINC: ', NPOINC
       SELECT CASE(TRIM(int_type))
          CASE("NAG")
             WRITE(6, '(A,EN12.3,A,A1)') '         Tol: ', follow_tol, '  Type: ', relab
          CASE("LSODE")
             WRITE(6, '(A,EN12.3,A,I2)') '         Tol: ', follow_tol, '  Type: ', mf
       END SELECT
       WRITE(6, '(5X,A,I3,A)', ADVANCE = 'no') 'Trajectory Calculation [', 0, ']%'
       CALL FLUSH(6)
    END IF

    ! Follow Trajectories
    IF (mystart_save <= nparticles) THEN
        SELECT CASE (TRIM(int_type))
            CASE ("NAG")
                ALLOCATE(w(neqs_nag * 21 + 28), STAT = ier)
                IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'W', ier)
                DO l = mystart_save, myend_save
                    t_nag = t_last(l)
                    ! Don't do particle if stopped or beyond the full_orbit limit
                    IF (t_nag>t_end(l)) CYCLE
                    ! Particle indicies
                    myline = l
                    mytdex = 1; ndt = 1
                    IF (lbeam) mytdex = 3
                    ! Don't do full_orbit particles
                    IF (sqrt(S_lines(mytdex-1,l))>rho_fullorbit) CYCLE
                    ! Particle Parameters
                    q(1) = R_lines(mytdex-1,l)
                    q(2) = PHI_lines(mytdex-1,l)
                    q(3) = Z_lines(mytdex-1,l)
                    q(4) = vll_lines(mytdex-1,l)
                    xlast = q(1)*cos(q(2))
                    ylast = q(1)*sin(q(2))
                    zlast = q(3)
                    moment = moment_lines(mytdex-1,l)
                    mycharge = charge(l)
                    myZ = Zatom(l)
                    mymass = mass(l)
                    E_by_v=mymass*0.5d-3/e_charge
                    mybeam = Beam(l)
                    my_end = t_end(l)
                    ltherm = .false.
                    lneut  = .false.
                    lcollision = lbeam
                    ! Collision parameters
                    fact_kick = 2*E_kick*mycharge/(mymass*pi2*pi2*freq_kick*freq_kick*SQRT(pi*1E-7*plasma_mass))
                    fact_pa   = plasma_mass/(mymass*plasma_Zmean)
                    fact_coul = myZ*(mymass+plasma_mass)/(mymass*plasma_mass*6.02214076208E+26)

                    ! Now calc dt
                    CALL beams3d_calc_dt(1,q(1),q(2),q(3),dt)
                    tf_nag = t_nag+dt
                    ndt = 1
                    DO ! Must do it this way becasue lbeam changes q(4) values
#if defined(NAG)
                       CALL D02CJF(t_nag,tf_nag,neqs_nag,q,fgc_eom,tol_nag,relab,out_beams3d_gc,D02CJW,w,ier)
#endif
                       IF (ier < 0) CALL handle_err(D02CJF_ERR, 'beams3d_follow', ier)
                       t_last(l) = tf_nag ! Save the value here in case out_beams3d changes it
                       CALL out_beams3d_gc(tf_nag,q)
                       IF (ABS(tf_nag) > ABS(my_end)) EXIT
                    END DO
                END DO
            CASE ("RKH68")
                ier = 0
                DO l = mystart_save, myend_save
                    tf_nag = t_last(l)
                    ! Don't do particle if stopped
                    IF (tf_nag>t_end(l)) CYCLE
                    ! Particle indicies
                    myline = l
                    mytdex = 1; ndt = 1
                    IF (lbeam) mytdex = 3
                    ! Don't do full_orbit particles
                    IF (sqrt(S_lines(mytdex-1,l))>rho_fullorbit) CYCLE
                    ! Particle Parameters
                    q(1) = R_lines(mytdex-1,l)
                    q(2) = PHI_lines(mytdex-1,l)
                    q(3) = Z_lines(mytdex-1,l)
                    q(4) = vll_lines(mytdex-1,l)
                    xlast = q(1)*cos(q(2))
                    ylast = q(1)*sin(q(2))
                    zlast = q(3)
                    moment = moment_lines(mytdex-1,l)
                    t_nag = tf_nag - dt
                    mycharge = charge(l)
                    myZ = Zatom(l)
                    mymass = mass(l)
                    E_by_v=mymass*0.5d-3/e_charge
                    mybeam = Beam(l)
                    my_end = t_end(l)
                    ltherm = .false.
                    lneut  = .false.
                    lcollision = lbeam
                    ! Collision parameters
                    fact_kick = 2*E_kick*mycharge/(mymass*pi2*pi2*freq_kick*freq_kick*SQRT(pi*1E-7*plasma_mass))
                    fact_pa   = plasma_mass/(mymass*plasma_Zmean)
                    fact_coul = myZ*(mymass+plasma_mass)/(mymass*plasma_mass*6.02214076208E+26)
					
                    ! Now calc dt
                    CALL beams3d_calc_dt(1,q(1),q(2),q(3),dt)
                    tf_nag = t_nag+dt
                    ndt = 1
                    ! Setup DRKHVG parameters
                    iopt = 0 
                    DO
                        CALL drkhvg(t_nag, q, neqs_nag, dt, 2, fgc_rkh68, rkh_work, iopt, ier)
                        IF (ier < 0) CALL handle_err(RKH68_ERR, 'beams3d_follow', ier)
                        q(1)=rkh_work(1,2)
                        q(2)=rkh_work(2,2)
                        q(3)=rkh_work(3,2)
                        q(4)=rkh_work(4,2)
                        t_nag = t_nag+dt
                        tf_nag = tf_nag+dt
                        t_last(l) = tf_nag ! Save the value here in case out_beams3d changes it
                        CALL out_beams3d_gc(tf_nag,q)
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
                DO l = mystart_save, myend_save
                    t_nag = t_last(l)
                    ! Don't do particle if stopped
                    IF (t_nag>t_end(l)) CYCLE
                    ! Particle indicies
                    myline = l
                    mytdex = 1
                    IF (lbeam) mytdex = 3
                    ! Don't do full_orbit particles
                    IF (S_lines(mytdex-1,l)>=s_fullorbit) CYCLE
                    ! Particle Parameters
                    q(1) = R_lines(mytdex-1,l)
                    q(2) = PHI_lines(mytdex-1,l)
                    q(3) = Z_lines(mytdex-1,l)
                    q(4) = vll_lines(mytdex-1,l)
                    xlast = q(1)*cos(q(2))
                    ylast = q(1)*sin(q(2))
                    zlast = q(3)
                    moment = moment_lines(mytdex-1,l)
                    mycharge = charge(l)
                    myZ = Zatom(l)
                    mymass = mass(l)
                    E_by_v=mymass*0.5d-3/e_charge
                    mybeam = Beam(l)
                    my_end = t_end(l)
                    ltherm = .false.
                    lneut  = .false.
                    lcollision = lbeam .or. lcollision
                    ! Collision parameters
                    fact_kick = 2*E_kick*mycharge/(mymass*pi2*pi2*freq_kick*freq_kick*SQRT(pi*1E-7*plasma_mass))
                    fact_pa   = plasma_mass/(mymass*plasma_Zmean)
                    fact_coul = myZ*(mymass+plasma_mass)/(mymass*plasma_mass*6.02214076208E+26)
					
                    ! Now calc dt
                    CALL beams3d_calc_dt(1,q(1),q(2),q(3),dt)
                    tf_nag = t_nag+dt
                    ndt = 1
                    ! Setup LSODE parameters
                    iopt = 0 ! No optional output
                    w = 0; iwork = 0; itask = 1; istate = 1;
                    itol = 2; rtol = follow_tol; atol(:) = follow_tol
                    DO
                        IF (lcollision) istate = 1
                        CALL FLUSH(6)
                        CALL DLSODE(fgc_lsode, neqs_nag, q, t_nag, tf_nag, itol, rtol, atol, itask, istate, &
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
                           CALL handle_err(LSODE_ERR, 'beams3d_follow_gc', istate)
                        END IF
                        iwork(11) = 0; iwork(12) = 0; iwork(13) = 0
                        t_last(l) = tf_nag ! Save the value here in case out_beams3d changes it
                        CALL out_beams3d_gc(tf_nag,q)
                        IF ( (istate == -1) .or. (istate ==-2) &
                                            .or. (ABS(tf_nag) > ABS(my_end)) &
                                            .or. (rho_help > rho_fullorbit)) EXIT
                    END DO
                END DO
                IF (ldebug) CLOSE(iunit)
             CASE ('DEBUG')
                DO l = 0, npoinc
                   R_lines(l,mystart_save:myend_save) = REAL(l)
                END DO
                B_lines(0:npoinc,mystart_save:myend_save) = REAL(myworkid)
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
    RETURN
    !-----------------------------------------------------------------------
    !     End Subroutine
    !-----------------------------------------------------------------------
END SUBROUTINE beams3d_follow_gc
