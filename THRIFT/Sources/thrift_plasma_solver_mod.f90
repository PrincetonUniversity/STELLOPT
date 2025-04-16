!-----------------------------------------------------------------------
!     Module:        thrift_plasma_solver_mod
!     Authors:       A. J. Coelho
!     Date:          04/XX/2025
!     Description:   This module contains the variables and subroutines
!                    concerning evolution of density and pressure equations
!-----------------------------------------------------------------------
MODULE thrift_plasma_solver_mod
    !-------------------------------------------------------------------
    !     Libraries
    !-------------------------------------------------------------------
    USE thrift_runtime
    USE safe_open_mod
    USE thrift_profiles_mod
    USE thrift_globals
    USE thrift_vars
    USE thrift_equil, ONLY : eq_Aminor, vp_spl
    !-------------------------------------------------------------------
    !     Module Variables
    !          lverb         Logical to control screen output
    !-------------------------------------------------------------------
    IMPLICIT NONE
    REAL(rprec) :: drho_plasma_solver, dr_plasma_solver
    REAL(rprec), DIMENSION(:), ALLOCATABLE :: rho_plasma_grid, r_plasma_grid, time_plasma_grid
    REAL(rprec), DIMENSioN(:), POINTER :: raxis_source, taxis_source
    REAL(rprec), DIMENSION(:,:,:,:), POINTER :: SE4D, Sn4D
    INTEGER :: ilogplasma, win_SE4D, win_Sn4D, win_rho_plasma_grid, num_species, nt_source, nrho_source, &
    win_raxis_source, win_taxis_source, win_r_plasma_grid
    REAL(rprec), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: plasma_N, plasma_T, plasma_P
    REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: plasma_N_keep, plasma_T_keep
    INTEGER, PRIVATE :: mytimestep_plasma_solver
    INTEGER :: N_plasma_steps_per_THRIFT_step, Nt_total_plasma_solver
    !           
    !REAL(rprec), PARAMETER ::
!-----------------------------------------------------------------------
!     Input Namelists
!         NONE
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!     Subroutines
!         evolve_plasma_equations: ...
!         read_external_plasma_sources:  ....
!-----------------------------------------------------------------------
    PUBLIC  :: evolve_plasma_equations, initialize_plasma_solver, &
    update_splines
    PRIVATE :: write_header_plasma_solver_logfile, &
    write_to_plasma_solver_logfile, update_pressure_and_temperature, &
    get_LHS_density_ions, solve_tridiagonal_system, &
    solve_sparse_nontridiag_system
    
    CONTAINS

    SUBROUTINE evolve_plasma_equations

        IMPLICIT NONE
        INTEGER :: istat, subiter, i, idx, irho, j, ispecies, ier, plasma_iteration
        REAL(rprec) :: t_current, t_old, t_new, rho, delta_p, delta_n
        REAL(rprec), DIMENSION(:), ALLOCATABLE :: pressure_total, pressure_total_old, ne_old
        REAL(rprec), DIMENSION(:), ALLOCATABLE :: RHS_density, lower_diag, upper_diag, main_diag, RHS_pressure
        real(rprec), DIMENSION(:,:), ALLOCATABLE :: LHS_pressure

        ! Allocate
        ! Module arrays
        IF( .NOT. ALLOCATED(plasma_N)) ALLOCATE(plasma_N(num_species,Nr_plasma_solver))
        IF( .NOT. ALLOCATED(plasma_T)) ALLOCATE(plasma_T(num_species,Nr_plasma_solver))
        IF( .NOT. ALLOCATED(plasma_P)) ALLOCATE(plasma_P(num_species,Nr_plasma_solver))
        IF( .NOT. ALLOCATED(plasma_N_keep)) ALLOCATE(plasma_N_keep(num_species,Nt_total_plasma_solver,Nr_plasma_solver))
        IF( .NOT. ALLOCATED(plasma_T_keep)) ALLOCATE(plasma_T_keep(num_species,Nt_total_plasma_solver,Nr_plasma_solver))
        ! Local arrays
        IF( .NOT. ALLOCATED(pressure_total)) ALLOCATE(pressure_total(num_species*Nr_plasma_solver))
        IF( .NOT. ALLOCATED(pressure_total_old)) ALLOCATE(pressure_total_old(num_species*Nr_plasma_solver))
        IF( .NOT. ALLOCATED(ne_old)) ALLOCATE(ne_old(Nr_plasma_solver))
        IF( .NOT. ALLOCATED(lower_diag)) ALLOCATE(lower_diag(Nr_plasma_solver-1))
        IF( .NOT. ALLOCATED(upper_diag)) ALLOCATE(upper_diag(Nr_plasma_solver-1))
        IF( .NOT. ALLOCATED(main_diag)) ALLOCATE(main_diag(Nr_plasma_solver))
        IF( .NOT. ALLOCATED(RHS_density)) ALLOCATE(RHS_density(Nr_plasma_solver))
        IF( .NOT. ALLOCATED(LHS_pressure)) ALLOCATE(LHS_pressure(Nr_plasma_solver*num_species,Nr_plasma_solver*num_species))
        IF( .NOT. ALLOCATED(RHS_pressure)) ALLOCATE(RHS_pressure(Nr_plasma_solver*num_species))


        IF(mytimestep .eq. 1) THEN

            mytimestep_plasma_solver = 1

            WRITE(6,*) 'Setting initial profiles...'
            CALL set_initial_profiles ! Use iniial profiles, that are either pre-defined or come from restart
            CALL update_plasma_keep
            WRITE(6,*) 'Initial profiles set!'

            ! Create log file where the info of the plasma solver will be written
            ilogplasma = 30
            CALL safe_open(ilogplasma, istat, 'plasma_solver.log', 'replace', 'formatted')
            IF (istat .ne. 0) STOP 'Error opening log file for plasma solver'

            ! Write header of log file
            CALL write_header_plasma_solver_logfile

            ! Write info of t=0 in log file
            CALL write_to_plasma_solver_logfile(time_plasma_grid(mytimestep_plasma_solver),0, &
            plasma_T_keep(1,mytimestep_plasma_solver,1), plasma_N_keep(1,mytimestep_plasma_solver,1), &
            plasma_T_keep(2,mytimestep_plasma_solver,1), plasma_N_keep(2,mytimestep_plasma_solver,1), &
            0.0_rprec,0.0_rprec)

            ! Update splines and exit routine
            CALL update_splines

            mytimestep_plasma_solver = mytimestep_plasma_solver + 1
            RETURN
        ENDIF
        
        r_plasma_grid = rho_plasma_grid * eq_Aminor
        dr_plasma_solver = r_plasma_grid(2) - r_plasma_grid(1)

        !
        t_old = THRIFT_T(mytimestep-1)
        t_new = THRIFT_T(mytimestep)

        ! Fill in plasma_N, plasma_P and plasma_T with previous time step
        plasma_N(:,:) = plasma_N_keep(:,mytimestep_plasma_solver-1,:)
        plasma_T(:,:) = plasma_T_keep(:,mytimestep_plasma_solver-1,:)
        plasma_P = plasma_N * plasma_T * e_charge

        ! Flatten plasma_P into plasma_total_old
        idx=1
        DO ispecies=1,num_species
            DO j=1,Nr_plasma_solver
                pressure_total_old(idx) = plasma_P(ispecies,j)
                idx = idx+1
            END DO
        END DO
        ! ne_old from previous time step
        ne_old = plasma_N(1,:)

        t_current = t_old
        DO plasma_iteration = 1,N_plasma_steps_per_THRIFT_step

            ! IN PYTHON WE MAKE THIS UPDATE HERE; BUT IN HERE I DON'T THINK WE NED IT
            ! CAUSE WE ARE NOT GOING TO SAVE IT; 
            ! IN CASE WANTS TO SAVE IT, HERE IS THE GOOD PLACE TO DO IT!
            !the first subiter corresponds to the last time step []
            ! plasma_N = 
            ! plasma_P =

            ! SUBCYCLE
            delta_p = 10*tol_plasma_solver
            delta_n = 10*tol_plasma_solver
            subiter=1
            DO WHILE ( (delta_p > tol_plasma_solver .OR. delta_n > tol_plasma_solver) .AND. (subiter .LE. max_subiter_plasma_solver) )

                ! Run PENTA if NEO fluxes are to be added
                IF(add_NEO) THEN
                    
                    ! ier = 0
                    ! CALL thrift_penta(.FALSE.,ier)
                    
                    ! PROBABLY MORE CORRECT TO DO THIS ??
                    CALL thrift_paraexe('penta',proc_string,lscreen_subcodes)

                    IF (ier /= 0) STOP 'Error running PENTA inside plasma solver'
                END IF
  
                DO i=1,nion_prof
                    CALL get_LHS_density_ions(i,lower_diag,main_diag,upper_diag)
                    CALL get_RHS_density_ions(i,RHS_density)
                    CALL solve_tridiagonal_system(lower_diag,main_diag,upper_diag, RHS_density, plasma_N(1+i,:))
                END DO
                ! electrons from quasi neutrality
                DO j=1,Nr_plasma_solver
                    plasma_N(1,j) = SUM(plasma_N(2:,j)*Zatom_prof)
                END DO

                CALL get_LHS_pressure(LHS_pressure)
                CALL get_RHS_pressure(RHS_pressure)
                CALL solve_sparse_nontridiag_system(LHS_pressure,RHS_pressure,pressure_total)

                ! update 'plasma_P' and 'plasma_T' with 'pressure_total' and 'plasma_N'
                CALL update_pressure_and_temperature(pressure_total)

                ! updates plasma_N_keep and plasma_T_keep
                CALL update_plasma_keep

                !compute errors delta_p and delta_n
                delta_p = MAXVAL(ABS(pressure_total - pressure_total_old)/pressure_total_old)
                delta_n = MAXVAL(ABS(ne_old - plasma_N(1,:))/ne_old) !only with ne since ne contains info from all densities which are actually evolved

                ! update olds
                pressure_total_old = pressure_total
                ne_old = plasma_N(1,:)

                ! update splines if add_NEO or beurskens
                IF(add_NEO .OR. beurskens_ions) THEN
                    STOP 'THIS IS NOT DONE YET... NEET TO UPDATE SPLINES TO USE PENTA PROPERLY; BEURSKENS FOR THE PROPER FIT OF dTd/dr I guess?'
                    CALL update_splines
                    ! measure time of updating splines at every subiter... measure its impact
                END IF

                ! write to plasma solver logfile
                CALL write_to_plasma_solver_logfile(time_plasma_grid(mytimestep_plasma_solver),subiter, &
                plasma_T_keep(1,mytimestep_plasma_solver,1), plasma_N_keep(1,mytimestep_plasma_solver,1), &
                plasma_T_keep(2,mytimestep_plasma_solver,1), plasma_N_keep(2,mytimestep_plasma_solver,1), &
                delta_p,delta_n)

                subiter = subiter+1

            END DO

            t_current = t_current + dt_plasma_solver
            mytimestep_plasma_solver = mytimestep_plasma_solver + 1
        END DO


        CALL update_splines()

        ! Deallocate local arrays
        DEALLOCATE(pressure_total)
        DEALLOCATE(pressure_total_old)
        DEALLOCATE(ne_old)
        DEALLOCATE(lower_diag)
        DEALLOCATE(upper_diag)
        DEALLOCATE(main_diag)
        DEALLOCATE(RHS_density)
        DEALLOCATE(LHS_pressure)
        DEALLOCATE(RHS_pressure)

        RETURN

    END SUBROUTINE evolve_plasma_equations


    SUBROUTINE initialize_plasma_solver(filename)
        ! sets external sources from external file; sets nion_prof, Zatom_prof, Matom_prof;
        ! and sets rho_plasma_grid and r_plasma_grid from Nr_plasma_solver
        ! and  allocates Nx3D, Tx3D, P3D
        USE mpi_inc
        USE mpi_params
        USE mpi_sharmem
        USE EZspline
        USE EZspline_obj
#if defined(LHDF5)
        USE ez_hdf5
#endif
        IMPLICIT NONE
        CHARACTER(*), INTENT(in) :: filename
        INTEGER :: i, ier
        INTEGER :: bcs0(2)
        TYPE(EZspline2_r8) :: temp_spl2d
        REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: S_energy, S_particle
        bcs0=(/ 0, 0/)
        ierr_mpi = 0

        IF (lverb) THEN
            WRITE(6,'(A)')  '----- Reading Sources File -----'
            WRITE(6,'(A)')  '   FILE: '//TRIM(filename)
        END IF

        IF (myid_sharmem == master) THEN
            CALL open_hdf5(TRIM(filename),fid,ier,LCREATE=.false.)
            IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,TRIM(filename),ier)
            
            CALL read_scalar_hdf5(fid,'nrho',ier,INTVAR=nrho_source)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nrho_source',ier)
            !
            CALL read_scalar_hdf5(fid,'nt',ier,INTVAR=nt_source)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nt_source',ier)
            
            CALL read_scalar_hdf5(fid,'nion',ier,INTVAR=nion_prof)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nion_prof',ier)
        END IF

        CALL MPI_BARRIER(MPI_COMM_SHARMEM,ierr_mpi)

        ! Broadcast nrho_source, nt_source, nion_prof
        CALL MPI_BCAST(nrho_source,1,MPI_INTEGER,master,MPI_COMM_SHARMEM,ierr_mpi)
        IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'read_thrift_profh5: nrho_source',ierr_mpi)
        CALL MPI_BCAST(nt_source,1,MPI_INTEGER,master,MPI_COMM_SHARMEM,ierr_mpi)
        IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'read_thrift_profh5: nt_source',ierr_mpi)
        CALL MPI_BCAST(nion_prof,1,MPI_INTEGER,master,MPI_COMM_SHARMEM,ierr_mpi)
        IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'bcasting: nion_prof',ierr_mpi)

        IF(nion_prof .GT. nions_max) STOP 'Number of ions larger that nions_max'

        num_species = nion_prof + 1

        ! Allocate the shared memory objects
        CALL mpialloc(raxis_source, nrho_source, myid_sharmem, 0, MPI_COMM_SHARMEM, win_raxis_source)
        CALL mpialloc(taxis_source, nt_source,   myid_sharmem, 0, MPI_COMM_SHARMEM, win_taxis_source)
        CALL mpialloc(Zatom_prof, nion_prof,   myid_sharmem, 0, MPI_COMM_SHARMEM, win_Zatom_prof)
        CALL mpialloc(Matom_prof, nion_prof,   myid_sharmem, 0, MPI_COMM_SHARMEM, win_Matom_prof)
        CALL mpialloc(SE4D, 4, nt_source, nrho_source, num_species, myid_sharmem, 0, MPI_COMM_SHARMEM, win_SE4D)
        CALL mpialloc(Sn4D, 4, nt_source, nrho_source, num_species, myid_sharmem, 0, MPI_COMM_SHARMEM, win_Sn4D)

        IF (myid_sharmem == master) THEN
            ! Read 
            CALL read_var_hdf5(fid,'Z_prof',nion_prof,ier,INTVAR=Zatom_prof)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'Zatom_prof',ier)
            CALL read_var_hdf5(fid,'mass_prof',nion_prof,ier,DBLVAR=Matom_prof)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'Matom_prof',ier)

            ! Get the axis arrays
            CALL read_var_hdf5(fid,'raxis_source',nrho_source,ier,DBLVAR=raxis_source)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'raxis_source',ier)
            CALL read_var_hdf5(fid,'taxis_source',nt_source,ier,DBLVAR=taxis_source)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'taxis_source',ier)

            ALLOCATE(S_energy(num_species,nt_source,nrho_source),S_particle(num_species,nt_source,nrho_source))

            !
            CALL read_var_hdf5(fid,'S_energy',num_species,nt_source,nrho_source,ier,DBLVAR=S_energy)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'S_energy',ier)
            CALL read_var_hdf5(fid,'S_particle',num_species,nt_source,nrho_source,ier,DBLVAR=S_particle)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'S_particle',ier)

            !
            DO i = 1,num_species
                ! Energy Source
                CALL EZspline_init(temp_spl2d,nt_source,nrho_source,bcs0,bcs0,ier)
                IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'init: S_energy',ier)
                temp_spl2d%x1          = taxis_source
                temp_spl2d%x2          = raxis_source
                temp_spl2d%isHermite   = 1
                CALL EZspline_setup(temp_spl2d,S_energy(i,:,:),ier,EXACT_DIM=.true.)
                IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'setup: S_energy',ier)
                SE4D(:,:,:,i) = temp_spl2d%fspl
                CALL EZspline_free(temp_spl2d,ier)
                ! Particle Source
                CALL EZspline_init(temp_spl2d,nt_source,nrho_source,bcs0,bcs0,ier)
                IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'init: S_particle',ier)
                temp_spl2d%x1          = taxis_source
                temp_spl2d%x2          = raxis_source
                temp_spl2d%isHermite   = 1
                CALL EZspline_setup(temp_spl2d,S_particle(i,:,:),ier,EXACT_DIM=.true.)
                IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'setup: S_particle',ier)
                Sn4D(:,:,:,i) = temp_spl2d%fspl
                CALL EZspline_free(temp_spl2d,ier)
            END DO

            DEALLOCATE(S_energy,S_particle)

            ! Close the HDF5 file
            CALL close_hdf5(fid,ier)
            IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,TRIM(filename),ier)

        END IF

        ! Plasma spatial grid (rho and r)
        ALLOCATE(rho_plasma_grid(Nr_plasma_solver))
        ALLOCATE(r_plasma_grid(Nr_plasma_solver))
        !
        FORALL(i = 1:Nr_plasma_solver)  rho_plasma_grid(i)  = DBLE(i-1)/DBLE(Nr_plasma_solver-1)
        !
        drho_plasma_solver = rho_plasma_grid(2) - rho_plasma_grid(1)

        CALL MPI_BARRIER(MPI_COMM_SHARMEM,ierr_mpi)

        IF (lverb) WRITE(6,*) 'Allocating Splines for plasma solver...'

        nrho_prof = SIZE(rho_plasma_grid) 
        nt_prof = Nt_total_plasma_solver

        ! Broadcast the helpers
        CALL MPI_BCAST(nrho_prof,1,MPI_INTEGER,master,MPI_COMM_SHARMEM,ierr_mpi)
        IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'read_thrift_profh5: nrho_prof',ierr_mpi)
        CALL MPI_BCAST(nt_prof,1,MPI_INTEGER,master,MPI_COMM_SHARMEM,ierr_mpi)
        IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'read_thrift_profh5: nt_prof',ierr_mpi)
        
        ! Allocate the shared memory objects
        CALL mpialloc(raxis_prof, nrho_prof, myid_sharmem, 0, MPI_COMM_SHARMEM, win_raxis_prof)
        CALL mpialloc(taxis_prof, nt_prof,   myid_sharmem, 0, MPI_COMM_SHARMEM, win_taxis_prof)

        CALL mpialloc(NE3D, 4, nt_prof, nrho_prof, myid_sharmem, 0, MPI_COMM_SHARMEM, win_NE3D)
        CALL mpialloc(TE3D, 4, nt_prof, nrho_prof, myid_sharmem, 0, MPI_COMM_SHARMEM, win_TE3D)
        CALL mpialloc(P3D,  4, nt_prof, nrho_prof, myid_sharmem, 0, MPI_COMM_SHARMEM, win_P3D)
        CALL mpialloc(NI4D, 4, nt_prof, nrho_prof, nion_prof, myid_sharmem, 0, MPI_COMM_SHARMEM, win_NI4D)
        CALL mpialloc(TI4D, 4, nt_prof, nrho_prof, nion_prof, myid_sharmem, 0, MPI_COMM_SHARMEM, win_TI4D)

        raxis_prof = rho_plasma_grid
        taxis_prof = time_plasma_grid ! Be careful, not defined yet; need to move this somewhere else!

        CALL setup_grids

        IF (lverb) WRITE(6,*) 'Splines Allocated!'

        RETURN

    END SUBROUTINE initialize_plasma_solver

    SUBROUTINE update_plasma_keep
        IMPLICIT NONE
        INTEGER :: ispecies
        DO ispecies=1,num_species
            plasma_N_keep(ispecies,mytimestep_plasma_solver,:) = plasma_N(ispecies,:)
            plasma_T_keep(ispecies,mytimestep_plasma_solver,:) = plasma_T(ispecies,:)  
        END DO
        RETURN
    END SUBROUTINE update_plasma_keep

    SUBROUTINE get_LHS_density_ions(iion,lower_diag,main_diag,upper_diag)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: iion
        REAL(rprec), DIMENSION(:), INTENT(INOUT) :: lower_diag, main_diag, upper_diag
        REAL(rprec) :: Dn_turb, cn_turb, dr, dr2, dt, rho, temp
        REAL(rprec) :: Vp_plus, Vp_minus, VDplus, VDminus, cplus, cminus, Dn_plus, Dn_minus
        REAL(rprec), DIMENSION(:), ALLOCATABLE :: Dn, cn, Vp
        INTEGER :: ir, ier, Nr

        Nr = Nr_plasma_solver

        lower_diag = 0.0_rprec
        main_diag = 0.0_rprec
        upper_diag = 0.0_rprec

        Dn_turb = Dn_ions(iion)
        cn_turb = cn_ions(iion)

        ALLOCATE(Dn(Nr),cn(Nr),Vp(Nr))

        IF(add_NEO) THEN
            STOP 'Not implemented yet...'
            !D_NEO = ... ! read from spline made in thrift_penta
            !Dn = D_NEO + Dn_turb
            !cn = c_NEO + cn_turb
        ELSE
            Dn = Dn_turb
            cn = cn_turb
        END IF

        dr = dr_plasma_solver
        dr2 = dr*dr
        dt = dt_plasma_solver

        ! Vp = dV/dr
        DO ir=1,Nr
            rho = rho_plasma_grid(ir)
            ier = 0
            CALL EZspline_interp(vp_spl, rho, temp, ier) ! temp = dV/dPhi
            Vp(ir) = 2.0_rprec * rho * THRIFT_PHIEDGE(mytimestep-1) * temp / eq_Aminor
        END Do

        ! r=0
        main_diag(1) = one + dt*4.0_rprec*Dn(1)/dr2 + 2.0_rprec*dt*cn(1)/dr
        upper_diag(1) = -4.0_rprec*dt*Dn(1)/dr2

        ! 0<r<a
        DO ir=2,Nr-1
            Vp_plus  = ( Vp(ir)+Vp(ir+1) ) / 2.0_rprec
            Vp_minus = ( Vp(ir)+Vp(ir-1) ) / 2.0_rprec

            Dn_plus  = ( Dn(ir)+Dn(ir+1) ) / 2.0_rprec
            Dn_minus = ( Dn(ir)+Dn(ir-1) ) / 2.0_rprec

            VDplus  = Vp_plus*Dn_plus / (Vp(ir)*dr2)
            VDminus = Vp_minus*Dn_minus / (Vp(ir)*dr2)

            cplus  = cn(ir+1)*Vp(ir+1) / (2*Vp(ir)*dr)
            cminus = cn(ir-1)*Vp(ir-1) / (2*Vp(ir)*dr)

            main_diag(ir) = one + dt*(VDplus + VDminus)
            upper_diag(ir) = dt*(-VDplus + cplus)
            lower_diag(ir-1) = dt*(-VDminus - cminus)
        END DO

        ! r=a
        main_diag(Nr) = one
        lower_diag(Nr-1) = 0.0_rprec

        DEALLOCATE(Dn,cn,Vp)

        RETURN
    END SUBROUTINE get_LHS_density_ions

    SUBROUTINE get_RHS_density_ions(iion,RHS_density)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: iion
        REAL(rprec), DIMENSION(:), INTENT(INOUT) :: RHS_density
        INTEGER :: Nr, ir, ispecies
        REAL(rprec) :: rho, t, t_prev, explicit_source, ni_previous

        Nr = Nr_plasma_solver
        ispecies = 1 + iion
        t = time_plasma_grid(mytimestep_plasma_solver)
        t_prev = time_plasma_grid(mytimestep_plasma_solver-1)

        DO ir=1,Nr
            rho = rho_plasma_grid(ir)
            CALL get_S_particle(rho,t,ispecies,explicit_source)
            !
            ni_previous = plasma_N_keep(ispecies,mytimestep_plasma_solver-1,ir) !CALL get_prof_ni(rho,t_prev,iion,ni_previous)
            RHS_density(ir) = ni_previous + dt_plasma_solver*explicit_source
        END DO

        ! Boundary condition
        RHS_density(Nr) = plasma_N(ispecies,Nr)

        RETURN
    END SUBROUTINE get_RHS_density_ions

    SUBROUTINE get_LHS_pressure(LHS_pressure)
        USE fusion_mod
        IMPLICIT NONE
        REAL(rprec), DIMENSION(:,:), INTENT(INOUT) :: LHS_pressure
        INTEGER :: ier, ir, Nr, ispecies, kk, row
        REAL(rprec) :: dt_fact, dr, dr2, Dp_turb, cp_turb, rho
        REAL(rprec) :: Vp_plus, Vp_minus, VDplus, VDminus, cplus, cminus, Dp_plus, Dp_minus
        REAL(rprec), DIMENSION(:), ALLOCATABLE :: Dp, cp, Vp
        REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: LHS_coll_heat_exchange

        LHS_pressure = 0.0_rprec

        Nr = Nr_plasma_solver

        dr = dr_plasma_solver
        dr2 = dr*dr
        dt_fact = (2.0_rprec/3.0_rprec)*dt_plasma_solver

        ALLOCATE(Dp(Nr),cp(Nr),Vp(Nr),LHS_coll_heat_exchange(Nr*num_species,Nr*num_species))

        ! Vp = dV/dr
        DO ir=1,Nr
            rho = rho_plasma_grid(ir)
            ier = 0
            ! CALL EZspline_interp(vp_spl, rho, temp, ier) ! temp = dV/dPhi
            ! Vp(ir) = 2.0_rprec * rho * THRIFT_PHIEDGE(mytimestep-1) * temp / eq_Aminor
            Vp(ir) = 4.0_rprec * pi * pi * 20.0_rprec * rho * eq_Aminor
        END DO

        kk = 1
        DO ispecies=1,num_species

            Dp_turb = chi_all(ispecies)
            cp_turb = cp_all(ispecies)

            IF(add_NEO) THEN
                STOP 'Not implemented yet...'
                !D_NEO = ... ! read from spline made in thrift_penta
                !Dp = D_NEO + Dp_turb
                !cp = c_NEO + cp_turb
            ! IF(beurskens_model) ...
            ELSE
                Dp = Dp_turb
                cp = cp_turb
            END IF

            ! r=0
            ! main_diag(1) = one + dt*4.0_rprec*Dp(1)/dr2 + 2.0_rprec*dt*cp(1)/dr
            LHS_pressure(kk,kk) = one + dt_fact*4.0_rprec*Dp(1)/dr2 + 2.0_rprec*dt_fact*cp(1)/dr
            ! upper_diag(1) = -4.0_rprec*dt*Dp(1)/dr2
            LHS_pressure(kk,kk+1) = -4.0_rprec*dt_fact*Dp(1)/dr2
            kk = kk+1

            ! 0<r<a
            DO ir=2,Nr-1
                Vp_plus  = ( Vp(ir)+Vp(ir+1) ) / 2.0_rprec
                Vp_minus = ( Vp(ir)+Vp(ir-1) ) / 2.0_rprec

                Dp_plus  = ( Dp(ir)+Dp(ir+1) ) / 2.0_rprec
                Dp_minus = ( Dp(ir)+Dp(ir-1) ) / 2.0_rprec

                VDplus  = Vp_plus*Dp_plus / (Vp(ir)*dr2)
                VDminus = Vp_minus*Dp_minus / (Vp(ir)*dr2)

                cplus  = cp(ir+1)*Vp(ir+1) / (2*Vp(ir)*dr)
                cminus = cp(ir-1)*Vp(ir-1) / (2*Vp(ir)*dr)

                ! main_diag(ir) = one + dt*(VDplus + VDminus)
                LHS_pressure(kk,kk) = one + dt_fact*(VDplus + VDminus)
                ! upper_diag(ir) = dt*(-VDplus + cplus)
                LHS_pressure(kk,kk+1) = dt_fact*(-VDplus + cplus)
                ! lower_diag(ir-1) = dt*(-VDminus - cminus)
                LHS_pressure(kk,kk-1) = dt_fact*(-VDminus - cminus)

                kk = kk+1
            END DO

            ! r=a (do not set yet BC's)
            kk = kk+1
        END DO

        ! Add collisional heat exchange
        CALL get_collisional_heat_exchange_matrix(LHS_coll_heat_exchange)
        LHS_pressure = LHS_pressure + dt_fact*LHS_coll_heat_exchange

        ! Impose Dirichlet Boundary Conditions
        DO ispecies=1,num_species
            row = ispecies * Nr
            ! Set the whole row to 0.0
            LHS_pressure(row,:) = 0.0_rprec
            ! Set diagonal entry to 1.0
            LHS_pressure(row,row) = one
        END DO

        DEALLOCATE(Dp,cp,Vp,LHS_coll_heat_exchange)

        RETURN
    END SUBROUTINE

    SUBROUTINE get_RHS_pressure(RHS_pressure)
        IMPLICIT NONE
        REAL(rprec), DIMENSION(:), INTENT(INOUT) :: RHS_pressure
        INTEGER :: Nr, ir, iion, offset
        REAL(rprec) :: rho, t, t_prev, explicit_source, n_previous, T_previous, p_previous

        Nr = Nr_plasma_solver

        t = time_plasma_grid(mytimestep_plasma_solver)
        t_prev = time_plasma_grid(mytimestep_plasma_solver-1)

        ! Electrons
        DO ir=1,Nr
            rho = rho_plasma_grid(ir)
            CALL get_S_energy(rho,t,1,explicit_source)
            n_previous = plasma_N_keep(1,mytimestep_plasma_solver-1,ir) !CALL get_prof_ne(rho,t_prev,n_previous)
            T_previous = plasma_T_keep(1,mytimestep_plasma_solver-1,ir) !CALL get_prof_Te(rho,t_prev,T_previous)
            p_previous = n_previous * T_previous * e_charge
            RHS_pressure(ir) = p_previous + (2.0_rprec/3.0_rprec)*dt_plasma_solver*explicit_source
            !Add Bremsstrahlung ...
            ! ... TBD
            ! Add alpha power ... ! ... detect which one is nD and nT with Zatom !...
            ! ... TBD
        END DO
        ! Boundary condition
        RHS_pressure(Nr) = plasma_P(1,Nr)

        ! Ions
        DO iion=1,nion_prof
            offset = Nr*iion
            DO ir=1,Nr-1
                rho = rho_plasma_grid(ir)
                CALL get_S_energy(rho,t,1+iion,explicit_source)
                !
                n_previous = plasma_N_keep(1+iion,mytimestep_plasma_solver-1,ir) !CALL get_prof_ni(rho,t_prev,iion,n_previous)
                T_previous = plasma_T_keep(1+iion,mytimestep_plasma_solver-1,ir) !CALL get_prof_Ti(rho,t_prev,iion,T_previous)
                p_previous = n_previous * T_previous * e_charge
                RHS_pressure(offset+ir) = p_previous + (2.0_rprec/3.0_rprec)*dt_plasma_solver*explicit_source
                ! Add alpha power ... ! ... detect which one is nD and nT with Zatom !...
            ! ... TBD 
            END DO
            ! Boundary condition
            RHS_pressure(offset+Nr) = plasma_P(1+iion,Nr)
        END DO

        RETURN
    END SUBROUTINE get_RHS_pressure

    SUBROUTINE update_splines
        ! updates splines NE3D, NI3D, TE3D, TI3D and P3D
        USE EZspline
        USE EZspline_obj
        IMPLICIT NONE
        INTEGER :: ispecies, i, ier, it
        REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: plasma_N_spline, plasma_T_spline
        REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: plasma_P_spline
        REAL(rprec), DIMENSION(:), ALLOCATABLE :: taxis
        INTEGER :: bcs0(2)
        TYPE(EZspline2_r8) :: temp_spl2d
        bcs0=(/ 0, 0/)

        ! Prepare data to update splies
        ALLOCATE(plasma_N_spline(num_species,nt_prof,nrho_prof))
        ALLOCATE(plasma_T_spline(num_species,nt_prof,nrho_prof))
        ALLOCATE(plasma_P_spline(nt_prof,nrho_prof))

        plasma_N_spline = plasma_N_keep 
        plasma_T_spline = plasma_T_keep 
        ! To make sure the spline doesn't mess-up if later asking value at t = current_time + eps
        ! we set all the fields at t>current_time equal to current_time
        DO it = mytimestep_plasma_solver+1, Nt_total_plasma_solver
            plasma_N_spline(:,it,:) = plasma_N
            plasma_T_spline(:,it,:) = plasma_T
        END DO
        
        plasma_P_spline = SUM(plasma_N_spline*plasma_T_spline*e_charge,dim=1)

        ! Create splines and re-write NE3D, TE3D, NI4D, TI4D and P3D
        ! NE
        CALL EZspline_init(temp_spl2d,nt_prof,nrho_prof,bcs0,bcs0,ier)
        temp_spl2d%x1          = taxis_prof
        temp_spl2d%x2          = raxis_prof
        temp_spl2d%isHermite   = 1
        CALL EZspline_setup(temp_spl2d,plasma_N_spline(1,:,:),ier,EXACT_DIM=.true.)
        NE3D = temp_spl2d%fspl
        CALL EZspline_free(temp_spl2d,ier)

        ! TE
        CALL EZspline_init(temp_spl2d,nt_prof,nrho_prof,bcs0,bcs0,ier)
        temp_spl2d%x1          = taxis_prof
        temp_spl2d%x2          = raxis_prof
        temp_spl2d%isHermite   = 1
        CALL EZspline_setup(temp_spl2d,plasma_T_spline(1,:,:),ier,EXACT_DIM=.true.)
        TE3D = temp_spl2d%fspl
        CALL EZspline_free(temp_spl2d,ier)

        ! NI
        DO i = 1, nion_prof
            CALL EZspline_init(temp_spl2d,nt_prof,nrho_prof,bcs0,bcs0,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'init: ni_prof',ier)
            temp_spl2d%x1          = taxis_prof
            temp_spl2d%x2          = raxis_prof
            temp_spl2d%isHermite   = 1
            CALL EZspline_setup(temp_spl2d,plasma_N_spline(1+i,:,:),ier,EXACT_DIM=.true.)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'setup: ni_prof',ier)
            NI4D(:,:,:,i) = temp_spl2d%fspl
            CALL EZspline_free(temp_spl2d,ier)
        END DO

        ! TI
        DO i = 1, nion_prof
           CALL EZspline_init(temp_spl2d,nt_prof,nrho_prof,bcs0,bcs0,ier)
           IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'init: ti_prof',ier)
           temp_spl2d%x1          = taxis_prof
           temp_spl2d%x2          = raxis_prof
           temp_spl2d%isHermite   = 1
           CALL EZspline_setup(temp_spl2d,plasma_T_spline(1+i,:,:),ier,EXACT_DIM=.true.)
           IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'setup: ti_prof',ier)
           TI4D(:,:,:,i) = temp_spl2d%fspl
           CALL EZspline_free(temp_spl2d,ier)
        END DO

        ! P
        CALL EZspline_init(temp_spl2d,nt_prof,nrho_prof,bcs0,bcs0,ier)
        temp_spl2d%x1          = taxis_prof
        temp_spl2d%x2          = raxis_prof
        temp_spl2d%isHermite   = 1
        CALL EZspline_setup(temp_spl2d,plasma_P_spline,ier,EXACT_DIM=.true.)
        P3D = temp_spl2d%fspl
        CALL EZspline_free(temp_spl2d,ier)

        DEALLOCATE(plasma_N_spline,plasma_T_spline,plasma_P_spline)
        RETURN

    END SUBROUTINE update_splines

    SUBROUTINE update_pressure_and_temperature(press_total)
        IMPLICIT NONE
        REAL(rprec), INTENT(IN), DIMENSION(:) :: press_total

        plasma_P = RESHAPE(press_total, SHAPE=(/num_species,Nr_plasma_solver/), ORDER=(/2,1/))
        plasma_T = plasma_P / (plasma_N * e_charge)
        RETURN
    END SUBROUTINE

    SUBROUTINE get_S_energy(rho_val,t_val,ispecies,val)
        USE EZspline_type
        IMPLICIT NONE
        REAL(rprec), INTENT(IN) :: rho_val, t_val
        REAL(rprec), INTENT(out) :: val
        INTEGER, INTENT(IN) :: ispecies
        REAL(rprec) :: rhomin_source, rhomax_source, tmin_source, tmax_source, eps1_source, eps2_source, t
        REAL(rprec), DIMENSION(:), ALLOCATABLE :: hr_source, ht_source, hri_source, hti_source
        INTEGER :: i, j
        REAL*8, parameter :: small = 1.e-10_ezspline_r8
        INTEGER, parameter :: ict(4)=(/1,0,0,0/)
        REAL*8  :: xparam, yparam
        REAL*8 :: fval(1)

        ALLOCATE(hr_source(nrho_source),ht_source(nt_source),hri_source(nrho_source),hti_source(nt_source))

        FORALL(i = 1:nrho_source-1) hr_source(i) = raxis_source(i+1) - raxis_source(i)
        FORALL(i = 1:nt_source-1)   ht_source(i) = taxis_source(i+1) - taxis_source(i)
        hri_source = one / hr_source
        hti_source = one / ht_source
        rhomin_source = MINVAL(raxis_source)
        rhomax_source = MAXVAL(raxis_source)
        tmin_source = MINVAL(taxis_source)
        tmax_source = MAXVAL(taxis_source)
        eps1_source = (rhomax_source-rhomin_source)*small
        eps2_source = (tmax_source-tmin_source)*small

        t = MIN(t_val,tmax_source)
        IF ((rho_val >= rhomin_source-eps1_source) .and. (rho_val <= rhomax_source+eps1_source) .and. &
         (t   >= tmin_source-eps2_source)   .and. (t   <= tmax_source+eps2_source)) THEN
            i = MIN(MAX(COUNT(taxis_source < t),1),nt_source-1)
            j = MIN(MAX(COUNT(raxis_source < rho_val),1),nrho_source-1)
            xparam = (t - taxis_source(i)) * hti_source(i)
            yparam = (rho_val   - raxis_source(j)) * hri_source(j)
            CALL R8HERM2FCN(ict,1,1,fval,i,j,xparam,yparam,&
                            ht_source(i),hti_source(i),hr_source(j),hri_source(j),&
                            SE4D(1,1,1,ispecies),nt_source,nrho_source)
            val = fval(1)
        ELSE
            STOP 'Why asking for energy sources outside domain?...'
        END IF

        DEALLOCATE(hr_source,ht_source,hri_source,hti_source)

        RETURN
    END SUBROUTINE

    SUBROUTINE get_S_particle(rho_val,t_val,ispecies,val)
        USE EZspline_type
        IMPLICIT NONE
        REAL(rprec), INTENT(IN) :: rho_val, t_val
        REAL(rprec), INTENT(out) :: val
        INTEGER, INTENT(IN) :: ispecies
        REAL(rprec) :: rhomin_source, rhomax_source, tmin_source, tmax_source, eps1_source, eps2_source, t
        REAL(rprec), DIMENSION(:), ALLOCATABLE :: hr_source, ht_source, hri_source, hti_source
        INTEGER :: i, j
        REAL*8, parameter :: small = 1.e-10_ezspline_r8
        INTEGER, parameter :: ict(4)=(/1,0,0,0/)
        REAL*8  :: xparam, yparam
        REAL*8 :: fval(1)

        ALLOCATE(hr_source(nrho_source),ht_source(nt_source),hri_source(nrho_source),hti_source(nt_source))

        FORALL(i = 1:nrho_source-1) hr_source(i) = raxis_source(i+1) - raxis_source(i)
        FORALL(i = 1:nt_source-1)   ht_source(i) = taxis_source(i+1) - taxis_source(i)
        hri_source = one / hr_source
        hti_source = one / ht_source
        rhomin_source = MINVAL(raxis_source)
        rhomax_source = MAXVAL(raxis_source)
        tmin_source = MINVAL(taxis_source)
        tmax_source = MAXVAL(taxis_source)
        eps1_source = (rhomax_source-rhomin_source)*small
        eps2_source = (tmax_source-tmin_source)*small

        t = MIN(t_val,tmax_source)
        IF ((rho_val >= rhomin_source-eps1_source) .and. (rho_val <= rhomax_source+eps1_source) .and. &
         (t   >= tmin_source-eps2_source)   .and. (t   <= tmax_source+eps2_source)) THEN
            i = MIN(MAX(COUNT(taxis_source < t),1),nt_source-1)
            j = MIN(MAX(COUNT(raxis_source < rho_val),1),nrho_source-1)
            xparam = (t - taxis_source(i)) * hti_source(i)
            yparam = (rho_val   - raxis_source(j)) * hri_source(j)
            CALL R8HERM2FCN(ict,1,1,fval,i,j,xparam,yparam,&
                            ht_source(i),hti_source(i),hr_source(j),hri_source(j),&
                            Sn4D(1,1,1,ispecies),nt_source,nrho_source)
            val = fval(1)
        ELSE
            STOP 'Why asking for energy sources outside domain?...'
        END IF

        DEALLOCATE(hr_source,ht_source,hri_source,hti_source)
        RETURN
    END SUBROUTINE

    SUBROUTINE set_initial_profiles
        IMPLICIT NONE
        INTEGER :: i, j

        IF(lrestart_from_file) THEN
            ! read here from restart ...
            PRINT *, 'reading plasma profiles from restart'
            STOP 'not implemented yet!!'

        ELSE
            ! ions: ni = 5E19 on axis, 0.8*5E-19 on edge, and quadratic decay
            DO i=1,nion_prof
                plasma_N(1+i,:) = 5.0E19_rprec * (0.8_rprec + 0.2_rprec*(1.0_rprec-rho_plasma_grid*rho_plasma_grid))
            END DO
            ! electrons from quasi neutrality
            DO j=1,Nr_plasma_solver
                plasma_N(1,j) = SUM(plasma_N(2:,j)*Zatom_prof)
            END DO

            ! T=200eV on axis, 0.8*200eV on edge, and quadratic decay (for all species)
            DO i=1,num_species
                plasma_T(i,:) = 200.0_rprec * (0.8_rprec + 0.2_rprec*(1.0_rprec-rho_plasma_grid*rho_plasma_grid))
            END DO

            plasma_P = plasma_N * plasma_T * e_charge

        END IF

        RETURN
    END SUBROUTINE set_initial_profiles

    SUBROUTINE solve_tridiagonal_system(lower_diag,main_diag,upper_diag,RHS_vec,result)
        IMPLICIT NONE
        REAL(rprec), DIMENSION(:), INTENT(IN) :: lower_diag, main_diag, upper_diag, RHS_vec
        REAL(rprec), DIMENSION(:), INTENT(INOUT) :: result
        INTEGER :: ier
        ier = 0
        CALL DGTSV(Nr_plasma_solver, 1, lower_diag, main_diag, upper_diag, RHS_vec, Nr_plasma_solver, ier)
        IF(ier/=0) STOP 'ERROR ON DENSITY SOLVER SOLUTION'
        result = RHS_vec
        ! Check result has non NaNs
        IF(ANY(ISNAN(result))) THEN
            PRINT *, 'results=', result
            STOP 'NaN values found on density. Exiting program...'
        END IF
        ! Look for negative values
        IF (ANY(result < 0.0)) THEN
            PRINT *, 'results=', result
            STOP 'Negative values found on density. Exiting program...'
        END IF
        RETURN
    END SUBROUTINE solve_tridiagonal_system

    SUBROUTINE solve_sparse_nontridiag_system(LHS_matrix,RHS_vec,result)
        IMPLICIT NONE
        REAL(rprec), DIMENSION(:,:), INTENT(IN) :: LHS_matrix
        REAL(rprec), DIMENSION(:), INTENT(IN) :: RHS_vec
        REAL(rprec), DIMENSION(:), INTENT(INOUT) :: result
        INTEGER :: ier, mat_size
        INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
        ier = 0
        mat_size = Nr_plasma_solver * num_species
        ALLOCATE(ipiv(mat_size))
        ! MIGHT WANT TO CHANGE THE WAY IPIV IS COMPUTED (USING SAMUEL ROUTINE)
        CALL DGESV(mat_size,1,LHS_matrix,mat_size,ipiv,RHS_vec,mat_size,ier)
        IF(ier/=0) STOP 'ERROR ON DENSITY SOLVER SOLUTION'
        result = RHS_vec
        ! Check result has non NaNs
        IF(ANY(ISNAN(result))) STOP 'NaN values found on pressure. Exiting program...'
        DEALLOCATE(ipiv)
        RETURN
    END SUBROUTINE solve_sparse_nontridiag_system

    SUBROUTINE get_collisional_heat_exchange_matrix(LHS_heat_exchange_matrix)
        USE collision_operators
        IMPLICIT NONE
        REAL(rprec), DIMENSION(:,:), INTENT(INOUT) :: LHS_heat_exchange_matrix
        REAL(rprec), DIMENSION(:), ALLOCATABLE :: mass_all, Z_all
        REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: W_s1_s2, aux_B
        INTEGER :: is1, is2, j, p, ir1, ir2, ir
        REAL(rprec) :: m1,n1,T1,m2,n2,T2,clog,const,vth_s1_sqr,vth_s2_sqr
        REAL(rprec) :: den, gamma, Z1, Z2

        LHS_heat_exchange_matrix = 0.0_rprec

        ALLOCATE(W_s1_s2(num_species,num_species,Nr_plasma_solver))
        ALLOCATE(aux_B(num_species,num_species,Nr_plasma_solver))
        ALLOCATE(mass_all(num_species),Z_all(num_species))

        mass_all(1) = electron_mass
        mass_all(2:) = Matom_prof
        Z_all(1) = -1.0_rprec
        Z_all(2:) = REAL(Zatom_prof, kind=rprec) ! need to be reals, because clog functions only accept reals

        DO is1=1,num_species
            m1 = mass_all(is1)
            Z1 = Z_all(is1)

            DO is2=1,num_species
                m2 = mass_all(is2)
                Z2 = Z_all(is2)

                DO ir=1,Nr_plasma_solver  
                    n1 = plasma_N(is1,ir)
                    T1 = plasma_T(is1,ir)
                    !
                    n2 = plasma_N(is2,ir)
                    T2 = plasma_T(is2,ir)

                    ! get Coulomb logarithm
                    IF (Z1>0 .AND. Z2>0) THEN
                        clog = COULOMB_LOG_NRL_II(m1,Z1,n1,T1,m2,Z2,n2,T2)
                    ELSE IF (Z1>0 .AND. Z2<0) THEN
                        clog = COULOMB_LOG_NRL_IE(n2,T2,m1,Z1,n1,T1)
                    ELSE IF (Z1<0 .AND. Z2>0) THEN
                        clog = COULOMB_LOG_NRL_IE(n1,T1,m2,Z2,n2,T2)
                    ELSE
                        clog = 0.0
                    END IF

                    const = (8/SQRT(pi))*(Z1*Z2*e_charge*e_charge)**2 * clog / (8*pi*EPS0**2)

                    vth_s1_sqr = 2*e_charge*T1/m1
                    vth_s2_sqr = 2*e_charge*T2/m2
                    
                    den = m1 * m2 * (vth_s1_sqr + vth_s2_sqr)**1.5_rprec
                    
                    gamma = const / den

                    W_s1_s2(is1,is2,ir) = gamma*n1   
                    aux_B(is1,is2,ir) = gamma*n2
                END DO
            END DO
        END DO

        ! Add aux_B matrix
        DO is1=1,num_species
            DO ir=1,Nr_plasma_solver
                W_s1_s2(is1,is1,ir) = W_s1_s2(is1,is1,ir) - SUM(aux_B(is1,:,ir))
            END DO
        END DO
        
        ! Fill LHS matrix
        j=1
        DO is1=1,num_species
            DO ir1=1,Nr_plasma_solver
                p=1
                DO is2=1,num_species
                    DO ir2=1,Nr_plasma_solver
                        ! Note the minus sign; this is to have it LHS
                        IF(ir1 .EQ. ir2) LHS_heat_exchange_matrix(j,p) = -W_s1_s2(is1,is2,ir2)
                        p = p+1
                    END DO
                END DO
                j = j+1
            END DO
        END DO

        DEALLOCATE(W_s1_s2,aux_B,mass_all,Z_all)
        RETURN
    END SUBROUTINE get_collisional_heat_exchange_matrix

    SUBROUTINE write_header_plasma_solver_logfile

        IMPLICIT NONE
        CHARACTER(len = 200) :: header_str

        header_str = '  T       NSUB    TE_AXIS [keV]    NE_AXIS[m-3]      TI1_AXIS [keV]    NI1_AXIS [m-3]      MAX(dp/p_old)      MAX(dn/n_old)'
               
        WRITE(ilogplasma,'(A)')' '
        WRITE(ilogplasma,'(A)') TRIM(header_str)
        WRITE(ilogplasma,'(A)')'  ========================================================================================================================'

    END SUBROUTINE write_header_plasma_solver_logfile

    SUBROUTINE write_to_plasma_solver_logfile(t,nsub,te_eV,ne,ti_eV,ni,max_dp,max_dn)

        IMPLICIT NONE
        INTEGER, INTENT(in) :: nsub
        REAL(rprec), INTENT(in) :: t,te_eV,ne,ti_eV,ni,max_dp,max_dn
        CHARACTER(len = 200) :: progress_str

        WRITE(progress_str,'(1X,F6.3,3X,I2,6X,F7.3,11X,ES8.2,11X,F7.3,13X,ES8.2,12X,ES8.2,11X,ES8.2)') &
                  t,nsub,te_eV/1000,ne,ti_eV/1000,ni,max_dp,max_dn
        WRITE(ilogplasma,'(A)') TRIM(progress_str)

    END SUBROUTINE write_to_plasma_solver_logfile


END MODULE thrift_plasma_solver_mod