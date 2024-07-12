      PROGRAM coilopt
      USE boundary
      USE bnorm_mod
      USE modular_coils
      USE saddle_coils
      USE saddle_surface
      USE bcoils_mod
      USE vf_coils
      USE tf_coils
      USE Vcoilpts
      USE Vname
      USE Vwire
      USE coils
      USE control_mod
      USE safe_open_mod
      USE mpi_params                                         !mpi stuff
      USE biotsavart, zerob => zero
      IMPLICIT NONE
#if defined(MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
#endif
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
#if defined(WIN32)
      CHARACTER(LEN=*), PARAMETER :: start_path = ".\"
#else
      !CHARACTER(LEN=*), PARAMETER :: start_path = "./"
      CHARACTER(LEN=*), PARAMETER :: start_path = ""
#endif
      INTEGER, PARAMETER :: nxc=1000
      INTEGER :: istat, ifunc, nvar, numargs, nfev
      INTEGER :: index_dot, index_for
      INTEGER :: i, n, iunit, color, key
      REAL(rprec) :: xc(nxc), bnorm1
      CHARACTER(LEN=100)  :: arg1
      REAL(rprec) :: start_time, finish_time, del_time, del_time_max
!-----------------------------------------------
#if defined(GEOM_ONLY)
#undef MPI_OPT
      lgeom_only = .true.
#endif

!     CALL MPI Initialization routines:

#if defined(MPI_OPT)
      CALL MPI_INIT( ierr_mpi )                                  !mpi stuff
      color = 0
      CALL MPI_COMM_SPLIT( MPI_COMM_WORLD,color,myid,
     1                     MPI_COMM_STEL,ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_STEL, myid, ierr_mpi )      
      CALL MPI_COMM_SIZE( MPI_COMM_STEL, numprocs, ierr_mpi )
#endif

      CALL second0(start_time)
!-----------------------------------------------
!
!     REF: D. J. STRICKLER, L. A. BERRY, S. P. HIRSHMAN (Fusion Science and
!          Technology, Vol 41, 2002, p 107.)
!

!     CALLING SYNTAX:

!     Xcoilopt [PATH][FILENAME.]EXTENSION

!     THE RECOMMENDED WAY TO RUN THIS IS FROM THE DIRECTORY WHERE
!     THE INPUT (FOR05 or INPUT_DATA_FILE) FILES EXIST. AN OPTIONAL
!     PATH IS ALLOWED, BUT CONFUSING...
!     THE FILES <INPUT or FOR05.EXTENSION>, <WOUT.EXTENSION>
!     MUST EXIST IN PATH. IF THEY DO NOT, THEN THIS IS BEING
!     RUN IN STELLOPT OPTIMIZER, AND THE NAMELIST IS LOOKED FOR
!     IN THE <INPUT_DATA_FILE.EXTENSION>
!
!     COILSNAMIN.F IN LIBSTELL\MODULES GIVES DEFINITIONS OF VARIABLES CONTROLLING
!     THE COILS AND THEIR OPTIMIZATION
!
!     RUN XBNORM TO GET MATCHING BNORM ON THE LCFS WHEN RUNNING IN COILOPT (NOT GEOM_ONLY) MODE
!     COEFFICIENTS ARE STORED IN BNORM OUTPUT FILE AND THE NESCIN FILE GIVES AN ESTIMATE
!     FOR THE WINDING SURFACE IN NESCOIL (COILOPT) COORDINATES, WITH TRIG ARGS = mu+nv (NOT
!     VMEC mu-nv)
!

!     READ COMMAND LINE EXTENSION

      CALL getcarg (1, arg1, numargs)
      IF( numargs .le. 0 ) THEN
         IF( myid .eq. master ) THEN                         !mpi stuff
            PRINT *,' MUST ENTER FILE SUFFIX ON COMMAND LINE'
         END IF     !IF( myid .eq. master )                  !mpi stuff
        STOP 01
      ENDIF

      index_dot = INDEX(arg1, '.', BACK=.TRUE.)
      IF (index_dot .gt. 2) THEN
         extension = arg1(index_dot+1:)      
      ELSE
         extension = TRIM(arg1)
      END IF
      for05_file = 'for05.' // TRIM(extension)
      for05_file_new = TRIM(for05_file) // '.new'
      input_data_file = 'input.' // TRIM(extension)
      
      path = start_path
      index_for = INDEX(arg1, "for05.")
      IF (index_for .le. 0) index_for = INDEX(arg1, "input.")
      IF (index_for .gt. 1) path = arg1(1:index_for-1)

      xc = zero
      nvar = -1

      CALL initialize_opt(xc, nxc, nvar, numargs)


!     READ COEFFICIENTS AND EVALUATE BNORM (NOTE: LBNORM ONLY USED IN FOLLOWING BLOCK OF CODE)
      IF (lplasma) THEN
         bnormal_match = zero
         IF (.not.lgeom_only) THEN
            IF (lbnorm) THEN
               IF (myid .eq. master) THEN
                  WRITE(6,*) '----- BNORM Information -----'
               END IF
               !IF (myid .eq. master) PRINT *, 'READ ', TRIM(bnorm_file)
               CALL read_bnorm_coefs
               CALL evaluate_bnorm
            ELSE IF (myid .eq. master) THEN
               PRINT *, 'LBNORM = FALSE: BNORM FILE WILL BE IGNORED'
            END IF
         END IF
      END IF

!     INITIALIZE POINTS FOR ACCESS ZONES
      IF (laccess) CALL initialize_access

!     PERFORM OPTIMIZATION

!     ifunc needs to be changed when fvec components are added

      ifunc = nedge + 4*nmod_unique_coils + nmod_coils_per_period
     1      + 9*nsad_unique_coils + n_access + 2*num_vf + 5
     2      +   nsad_unique_coils*nsad_coils
      nfev = niter_opt*nvar
      PRINT *,nvar
      CALL optimize (xc, nvar, ifunc, epsfcn, nfev, istat)
      IF (istat .ge. 9) THEN
         WRITE (6, *)' EXITING OPTIMIZE FROM COILOPT, ISTAT = ', istat
         GOTO 2000
      END IF

      IF (myid .EQ. master) THEN
         iunit = 14
         CALL safe_open(iunit, istat, 'for05', 'unknown',
     1        'formatted')
         IF (istat .EQ. 0) CLOSE(iunit,iostat=istat,status='delete')

!     WRITE b-norm error to file

         IF (.NOT.lgeom_only) THEN

            CALL loadparams (xc, nvar)
            IF (lmodular) CALL eval_modular_coils
            IF (lsaddle) CALL eval_saddle_coils
            IF (lvf) CALL eval_vf_coils
            CALL evaluate_field_error

            iunit = 23
            CALL safe_open(iunit, istat, 'b_norm_error.dat', 'unknown',
     1                     'formatted')
            WRITE (iunit, '(a)')'  NEDGE      BN-MATCH      BN-ERROR'
            WRITE (iunit, *)
            bnorm1 = SQRT(SUM(bnormal_match**2)/nedge)
            DO n = 1, nedge
               WRITE(iunit, 1000) n, bnormal_match(n), b_error(n)*bnorm1
            END DO

            CLOSE (iunit)
         END IF

 1000 FORMAT(i10,1p,2e14.4)

         PRINT 1030, rbphi_avg
 1030    FORMAT(/,' <R*Bphi> [T-m] = ',1pe14.4)


         IF (.NOT.lgeom_only .AND. lmodular) THEN

!        Evaluate max field at coils

            iunit = 29
            CALL safe_open(iunit, istat, 'surf_norm.dat', 'unknown',
     1                     'formatted')
            CALL evaluate_max_field (iunit)
            PRINT 1100
            DO i=1, nmid
               PRINT 1110, i, b_max(i), curmod(i)
            END DO
 1100 FORMAT (" Coil",4x,"BMAX (T)",5x,"Imod (Amp)")
 1110 FORMAT (i4,1p,2e14.4)
               CLOSE (iunit)
         END IF


!        Write coil and surface data to file

         CALL write_coils (extension)
         IF (.not.lgeom_only) THEN
            CALL write_surfaces
            CALL evaluate_bg_field
         END IF

      END IF

 2000 CONTINUE

!     FREE MEMORY

      IF (ALLOCATED(xm_bmn)) DEALLOCATE (xm_bmn, xn_bmn, bmn)
      IF (ALLOCATED(rb)) DEALLOCATE (rb, zb, phib, thetab, x_p, y_p,
     1  z_p, rb_ph, rb_th, zb_ph, zb_th, n_r, n_phi, n_z, d_area,
     2  theta_d, phi_d, bnormal_match, b_error, b_mod, bmn_error,
     3  bsl_error, luv)
      IF (ALLOCATED(rmnc_b)) DEALLOCATE (rmnc_b, zmns_b,
     1  xm_b, xn_b, lmns_b)
      IF (ALLOCATED(nbrho_opt)) DEALLOCATE (nbrho_opt, mbrho_opt,
     1  rbc, zbs, rhobc, nrz0_opt, delta_mn)
      IF (ALLOCATED(bcoil_x)) DEALLOCATE (bcoil_x, bcoil_y, bcoil_z)
      IF (ALLOCATED(x_vf)) DEALLOCATE (x_vf, y_vf, z_vf, vertical)
      IF (ALLOCATED(modular)) DEALLOCATE (modular)
      IF (ALLOCATED(saddle)) DEALLOCATE (saddle)
      IF (ALLOCATED(tfc_x)) DEALLOCATE (tfc_x, tfc_y, tfc_z)
      IF (ALLOCATED(x_mod)) DEALLOCATE (x_mod, y_mod, z_mod)      
      IF (ALLOCATED(rho)) DEALLOCATE (rho, phi, rcoil, zcoil)

      CALL cleanup_biotsavart

      CALL second0(finish_time)
      del_time = finish_time - start_time
#if defined(MPI_OPT)
      CALL MPI_Reduce(del_time, del_time_max, 1, MPI_REAL8,
     1  MPI_max, master, MPI_COMM_WORLD, ierr_mpi)
#else
      del_time_max = del_time
#endif
      IF(myid .eq. master) THEN
         WRITE(*,'("Elapsed time (from rtc) = ",e15.7)') del_time_max
      ENDIF
#if defined(MPI_OPT)
      CALL MPI_FINALIZE(ierr_mpi)                            !mpi stuff
#endif
      END PROGRAM coilopt
