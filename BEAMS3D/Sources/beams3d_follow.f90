!-----------------------------------------------------------------------
!     Module:        beams3d_follow
!     Authors:       S. Lazerson (lazerson@pppl.gov), M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          06/20/2012
!     Description:   This subroutine follows the particles through
!                    the grid.  The general ODE which must be solved
!                    can be written:
!                        dX        
!                       ----   = fpart(t,q,qdot)
!                        dt        
!                    where X=(R,phi,Z).  This allows particle trajectories
!                    to be parametrized in terms of time.
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
                            TE_spl, TI_spl
    USE mpi_params ! MPI
    USE beams3d_physics_mod
    USE beams3d_write_par
    USE safe_open_mod, ONLY: safe_open
    USE wall_mod, ONLY: wall_free, ihit_array, nface
!DEC$ IF DEFINED (MPI_OPT)
    USE mpi
!DEC$ ENDIF
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
!DEC$ IF DEFINED (MPI_OPT)
    INTEGER :: status(MPI_STATUS_size) !mpi stuff
    INTEGER :: mystart, mypace, i, j, sender
    INTEGER,ALLOCATABLE :: mnum(:), moffsets(:)
    INTEGER :: MPI_COMM_LOCAL
!DEC$ ENDIF
    INTEGER :: ier, l, neqs_nag, l2, itol, itask, &
               istate, iopt, lrw, liw, mf, out, iunit
    INTEGER, ALLOCATABLE :: iwork(:), itemp(:,:)
    REAL :: dist
    REAL(rprec) :: tf_max, vel_max
    DOUBLE PRECISION, ALLOCATABLE :: w(:), q(:)
    DOUBLE PRECISION :: tf_nag, eps_temp, t_nag, t1_nag, &
                        tol_nag, rtol
    DOUBLE PRECISION :: atol(4), rwork(84)
    DOUBLE PRECISION :: rkh_work(4, 2)
    DOUBLE PRECISION :: qdot1(4)
    CHARACTER*1 :: relab

    DOUBLE PRECISION, PARAMETER :: lendt_m = 0.05  ! Maximum length to travel before updating physics [m]
    !-----------------------------------------------------------------------
    !     External Functions
    !          fpart_nag            RHS of ODE integrator (for NAG)    for BEAMS3D
    !          D02CJF               NAG ODE Solver
    !          D02CJW               NAG Dummy function
    !          jacobian_lsode       Jacobian function (for LSODE, not currently utilized)
    !-----------------------------------------------------------------------
    EXTERNAL D02CJF, D02CJW, fpart_nag, D02CJX, out_beams3d_nag
    EXTERNAL fpart_lsode, jacobian_lsode
    EXTERNAL fpart_rkh68
    !-----------------------------------------------------------------------
    !     Begin Subroutine
    !-----------------------------------------------------------------------
    ! Initializations
    ier = 0
    i = MAXLOC(ABS(t_end),1)
    tf_max = t_end(i)
    dt_out = tf_max/npoinc
    vel_max = MAXVAL(vll_start)
    dt = SIGN(lendt_m/vel_max,tf_max)      ! Keep this here so we print out max(dt)
    IF (ABS(dt) < 1E-9) dt = SIGN(1D-9,tf_max)  ! This is a limiter term for STELLOPT
    nsteps = FLOOR(tf_max/dt)
    tol_nag = follow_tol
    neqs_nag = 4
    relab = "M"
    mf = 10
    ! Do a quick check
    IF (nsteps < npoinc) THEN
        npoinc = nsteps
    END IF

    ! Break up the work
    mypace  = nparticles
    mystart = 1
    myend   = nparticles
    IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
    ALLOCATE(mnum(nprocs_beams))
    mnum = 0
    DO
       IF (SUM(mnum) == nparticles) EXIT
       DO i = 1, nprocs_beams
          mnum(i) = mnum(i) + 1
          IF (SUM(mnum) == nparticles) EXIT
       END DO
       IF (SUM(mnum) == nparticles) EXIT
    END DO
    mypace  = mnum(myworkid+1)
    IF (myworkid == 0) THEN
       mystart = 1
    ELSE
       mystart = SUM(mnum(1:myworkid))+1
    END IF
    myend   = mystart+mypace-1
    DEALLOCATE(mnum)

    ! Save mystart and myend
    mystart_save = mystart
    myend_save = myend

    IF (lhitonly) THEN
        npoinc = 2
        dt_out = 0
    END IF

!DEC$ IF DEFINED (MPI_OPT)
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
      ALLOCATE(mnum(nprocs_beams), moffsets(nprocs_beams))
      CALL MPI_ALLGATHER((myend-mystart+1)*(npoinc+1),1,MPI_INTEGER,mnum,1,MPI_INTEGER,MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_ALLGATHER((mystart-1)*(npoinc+1),1,MPI_INTEGER,moffsets,1,MPI_INTEGER,MPI_COMM_BEAMS,ierr_mpi)
!DEC$ ENDIF

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
    IF (ALLOCATED(moment_lines)) DEALLOCATE(moment_lines)
    IF (ALLOCATED(neut_lines)) DEALLOCATE(neut_lines)
    IF (ALLOCATED(lost_lines)) DEALLOCATE(lost_lines)
    IF (ALLOCATED(PE_lines)) DEALLOCATE(PE_lines)
    IF (ALLOCATED(PI_lines)) DEALLOCATE(PI_lines)
    IF (ALLOCATED(j_lines)) DEALLOCATE(j_lines)
    
    ! Output some stuff
    IF (lverb) THEN
       WRITE(6, '(A)') '----- FOLLOWING PARTICLE TRAJECTORIES -----'
       WRITE(6, '(A,A)')          '      Method: ', TRIM(int_type)
       WRITE(6, '(A,I9)')         '   Particles: ', nparticles
       WRITE(6, '(A,I9,A,E11.4)') '       Steps: ', nsteps, '   Delta-t: ', dt
       WRITE(6, '(A,I9,A,E11.4)') '      NPOINC: ', npoinc, '    dt_out: ', dt_out
       SELECT CASE(TRIM(int_type))
          CASE("NAG")
             WRITE(6, '(A,E11.4,A,A1)') '         Tol: ', follow_tol, '  Type: ', relab
          CASE("LSODE")
             WRITE(6, '(A,E11.4,A,I2)') '         Tol: ', follow_tol, '  Type: ', mf
       END SELECT
       WRITE(6, '(5X,A,I3,A)', ADVANCE = 'no') 'Trajectory Calculation [', 0, ']%'
       CALL FLUSH(6)
    END IF
    ! Allocations
    ALLOCATE(q(neqs_nag), STAT = ier)
    IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'Q', ier)
    ALLOCATE(R_lines(0:npoinc, mystart:myend), Z_lines(0:npoinc, mystart:myend), &
             PHI_lines(0:npoinc, mystart:myend), vll_lines(0:npoinc, mystart:myend), moment_lines(0:npoinc, mystart:myend), &
             neut_lines(0:npoinc, mystart:myend),PE_lines(0:npoinc, mystart:myend),PI_lines(0:npoinc, mystart:myend), &
             S_lines(0:npoinc, mystart:myend), U_lines(0:npoinc, mystart:myend), B_lines(0:npoinc, mystart:myend), &
             j_lines(0:npoinc, mystart:myend), STAT = ier)
    IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'R_LINES, PHI_LINES, Z_LINES', ier)
    ALLOCATE(lost_lines(mystart:myend), STAT = ier)
    IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'LOST_LINES', ier)

    ! Initializations
    R_lines = 0.0; Z_lines = 0.0; PHI_lines = -1.0
    vll_lines = 0.0; moment_lines = 0.0
    PE_lines = 0.0; PI_lines = 0.0; j_lines = 0.0
    lost_lines = .FALSE.
    S_lines = 1.5; U_lines = 0.0; B_lines = -1.0
    R_lines(0, mystart:myend) = R_start(mystart:myend)
    Z_lines(0, mystart:myend) = Z_start(mystart:myend)
    PHI_lines(0, mystart:myend) = phi_start(mystart:myend)
    vll_lines(0, mystart:myend) = vll_start(mystart:myend)
    moment_lines(0, mystart:myend) = mu_start(mystart:myend)
    neut_lines(0, mystart:myend) = .FALSE.
    IF (lbeam) neut_lines(0, mystart:myend) = .TRUE.

    ! Follow Trajectories
    IF (mystart <= nparticles) THEN
        SELECT CASE (TRIM(int_type))
            CASE ("NAG")
!DEC$ IF DEFINED (NAG)
                ALLOCATE(w(neqs_nag * 21 + 28), STAT = ier)
                IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'W', ier)
                DO l = mystart, myend
                    ltherm = .false.
                    lneut  = .false.
                    q(1) = R_start(l)
                    q(2) = phi_start(l)
                    q(3) = Z_start(l)
                    q(4) = vll_start(l)
                    B_temp(:) = 1.0
                    t_nag = 0.0
                    tf_nag = 0.0
                    mycharge = charge(l)
                    myZ = Zatom(l)
                    mymass = mass(l)
                    moment = mu_start(l)
                    myv_neut(:) = v_neut(:,myline)
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
                       IF (ldepo) CYCLE
                       ltherm = .FALSE.
                       lcollision = .TRUE.
                       mytdex = 3
                       tf_nag = tf_nag - dt  ! Because out advances t by dt
                       dt_out = (t_end(l) - t_nag)/(npoinc-2) ! Adjust dt slightly to keep indexing correct.
                    END IF
                    DO ! Must do it this way becasue lbeam changes q(4) values
                       CALL D02CJF(t_nag,tf_nag,neqs_nag,q,fpart_nag,tol_nag,relab,out_beams3d_nag,D02CJW,w,ier)
                       IF (ier < 0) CALL handle_err(D02CJF_ERR, 'beams3d_follow', ier)
                       CALL out_beams3d_nag(tf_nag,q)
                       IF (tf_nag > t_end(l)) EXIT
                    END DO
                END DO
!DEC$ ELSE
                ier = -1
                CALL handle_err(NAG_ERR, 'beams3d_follow', ier)
!DEC$ ENDIF
            CASE ("RKH68")
                ier = 0
                DO l = mystart, myend
                    ltherm = .false.
                    lneut  = .false.
                    q(1) = R_start(l)
                    q(2) = phi_start(l)
                    q(3) = Z_start(l)
                    q(4) = vll_start(l)
                    B_temp(:) = 1.0
                    t_nag = 0.0
                    tf_nag = 0.0
                    mycharge = charge(l)
                    myZ = Zatom(l)
                    mymass = mass(l)
                    moment = mu_start(l)
                    myv_neut(:) = v_neut(:,myline)
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
                       IF (ldepo) CYCLE
                       ltherm = .FALSE.
                       lcollision = .TRUE.
                       mytdex = 3
                       tf_nag = tf_nag - dt  ! Because out advances t by dt
                       dt_out = (t_end(l) - t_nag)/(npoinc-2) ! Adjust dt slightly to keep indexing correct.
                    END IF
                    DO
                        CALL drkhvg(t_nag, q, neqs_nag, dt, 2, fpart_rkh68, rkh_work, iopt, ier)
                        IF (ier < 0) CALL handle_err(RKH68_ERR, 'beams3d_follow', ier)
                        q(1)=rkh_work(1,2)
                        q(2)=rkh_work(2,2)
                        q(3)=rkh_work(3,2)
                        q(4)=rkh_work(4,2)
                        t_nag = t_nag+dt
                        tf_nag = tf_nag+dt
                        CALL out_beams3d_nag(tf_nag,q)
                        IF (tf_nag > t_end(l)) EXIT
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
                    ! Initialize the calculation
                    ltherm = .false.
                    lneut  = .false.
                    q(1) = R_start(l)
                    q(2) = phi_start(l)
                    q(3) = Z_start(l)
                    q(4) = vll_start(l)
                    xlast = q(1)*cos(q(2))
                    ylast = q(1)*sin(q(2))
                    zlast = q(3)
                    B_temp(:) = 1.0
                    t_nag = 0.0
                    tf_nag = 0.0
                    mycharge = charge(l)
                    myZ = Zatom(l)
                    mymass = mass(l)
                    mybeam = Beam(l)
                    moment = mu_start(l)
                    myv_neut(:) = v_neut(:,myline)
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
                       IF (ldepo) CYCLE
                       ltherm = .FALSE.
                       lcollision = .TRUE.
                       mytdex = 3
                       tf_nag = tf_nag - dt  ! Because out advances t by dt
                       dt_out = (t_end(l) - t_nag)/(npoinc-2) ! Adjust dt slightly to keep indexing correct.
                       CALL FLUSH(6)
                    END IF
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
                        CALL out_beams3d_nag(tf_nag,q)
                        IF ((istate == -1) .or. (istate ==-2) .or. (ABS(tf_nag) > ABS(t_end(l))) ) EXIT
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

!       Deallocations
    IF (ALLOCATED(q)) DEALLOCATE(q)
    IF (ALLOCATED(w)) DEALLOCATE(w)
    IF (ALLOCATED(iwork)) DEALLOCATE(iwork)

    ! Adjust PE_lines and PI_lines to be in real units
!    DO i = mystart, myend
!       PE_lines(:,i) = PE_lines(:,i)*weight(i)
!       PI_lines(:,i) = PI_lines(:,i)*weight(i)
!       j_lines(:,i)  = j_lines(:,i)*weight(i)
!    END DO

    ! First reduce the cumulative arrays over shared memory groups then allreduce between shared memeory groups
!DEC$ IF DEFINED (MPI_OPT)
    IF (myid_sharmem == master) THEN
       CALL MPI_REDUCE(MPI_IN_PLACE, epower_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(MPI_IN_PLACE, ipower_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(MPI_IN_PLACE,   ndot_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(MPI_IN_PLACE,      j_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(MPI_IN_PLACE,   dist_prof, nbeams*ns_prof*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
    ELSE
       CALL MPI_REDUCE(epower_prof, epower_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(ipower_prof, ipower_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(ipower_prof,   ndot_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(ipower_prof,      j_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
       CALL MPI_REDUCE(ipower_prof,   dist_prof, nbeams*ns_prof*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_SHARMEM, ierr_mpi)
    END IF
!DEC$ ENDIF

!DEC$ IF DEFINED (MPI_OPT)
    i = MPI_UNDEFINED
    IF (myid_sharmem == master) i = 0
    CALL MPI_COMM_SPLIT( MPI_COMM_BEAMS,i,myworkid,MPI_COMM_LOCAL,ierr_mpi)
    IF (myid_sharmem == master) THEN
!       partvmax = MAXVAL(MAXVAL(ABS(vll_lines),DIM=2),DIM=1)
!       CALL MPI_ALLREDUCE(MPI_IN_PLACE, partvmax,                 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_LOCAL, ierr_mpi)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, epower_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_LOCAL, ierr_mpi)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, ipower_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_LOCAL, ierr_mpi)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,   ndot_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_LOCAL, ierr_mpi)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,      j_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_LOCAL, ierr_mpi)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,   dist_prof, nbeams*ns_prof*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_LOCAL, ierr_mpi)
       IF (ASSOCIATED(ihit_array)) THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,ihit_array,nface,MPI_INTEGER,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
       END IF
       CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
    END IF
    CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)
!DEC$ ELSE
!    partvmax = MAXVAL(MAXVAL(ABS(vll_lines),DIM=2),DIM=1)
!DEC$ ENDIF


!DEC$ IF DEFINED (MPI_OPT)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,      'R_lines', DBLVAR=R_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,    'PHI_lines', DBLVAR=PHI_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,      'Z_lines', DBLVAR=Z_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,    'vll_lines', DBLVAR=vll_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend, 'moment_lines', DBLVAR=moment_lines)
    !CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,     'PE_lines', DBLVAR=PE_lines)
    !CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,     'PI_lines', DBLVAR=PI_lines)
    !CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,      'j_lines', DBLVAR=j_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,      'S_lines', DBLVAR=S_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,      'U_lines', DBLVAR=U_lines)
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,      'B_lines', DBLVAR=B_lines)
    ALLOCATE(itemp(0:npoinc,mystart:myend))
    itemp = 0; WHERE(neut_lines) itemp=1;
    CALL beams3d_write_parhdf5(0, npoinc, 1, nparticles, mystart, myend,   'neut_lines', INTVAR=itemp)
    DEALLOCATE(itemp)
    IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
    IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
    CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= 0) CALL handle_err(MPI_BARRIER_ERR, 'beams3d_follow', ierr_mpi)

!DEC$ ENDIF
    RETURN
    !-----------------------------------------------------------------------
    !     End Subroutine
    !-----------------------------------------------------------------------
END SUBROUTINE beams3d_follow
