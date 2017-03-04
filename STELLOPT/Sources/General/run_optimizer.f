      SUBROUTINE run_optimizer(xc_opt, var_descript, 
     1                         nopt, nvar, lwa, info)
      USE stel_kinds
      USE optim, ONLY: home_dir, lone_step, lrestart, describe_string,
     1                 constraint_string
      USE optim_params, ONLY: epsfcn, niter_opt, seq_ext,
     1   num_processors, num_levmar_params, nopt_alg
      USE vparams, ONLY: zero
      USE chisq_mod, ONLY: chisq_descript
      USE mpi_params   
      USE safe_open_mod                                                         ! MPI
      IMPLICIT NONE
C-----------------------------------------------
C     D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nopt, nvar, lwa
      REAL(rprec), DIMENSION(nvar) :: xc_opt
      CHARACTER(len=*), DIMENSION(nvar) :: var_descript                         !description of independent variables
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: info_size = 36
      CHARACTER(len=*), PARAMETER :: 
     1   stop_string = 'Allocation error in STELLOPT run-optimizer!'
      INTEGER :: info, niter, iunit_vars, ier
      INTEGER :: mode
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: fvec, diag
      REAL(rprec) :: tol
      CHARACTER(LEN=120), DIMENSION(1:info_size) :: info_array
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL lsfun1
C-----------------------------------------------
      DATA info_array/
     1     'error in read_wout_opt opening wout file',

     2     'error in load_physics_codes in system call',

     3     'error opening indata file in write_indata',

     4     'error opening output file in lsfun',

     5     'error writing new input file in load_params',

     6     'vmec2000 executable file not found',

     7     'i/o error in clean_up routine',

     8     'error reading wout file in call to read_wout_opt',

     9     'allocation error in load_target',

     a     'boozer array dimension mismatch in load_target',                             !10

     b     'nvar and nopt do not match in load_target',

     c     'i/o error opening cobra file in chisq_ballooning',

     d     'i/o error opening bootstrap file in chisq_bootsj',

     e     'i/o error in open_comm_files',

     f     'error in chk_rzmnb',

     g     'boozer transform module - xbooz_xform - not found',

     h     'could not locate executable in load_physics_codes',

     i     'error reading output file in chisq_jinvar subroutine',

     j     'error in external kink computation',

     k     'error running xdkes code in chisq_dkes subroutine',                          !20

     l     'error in vacuum vessel matching subroutine',

     m     'system call to xcoilgeom failed in generate_mgrid',

     n     'coils data file was not produced by xcoilgeom',

     o     'error opening extcur file in generate_mgrid',

     p     'error opening coil_targets file in chisq_coilgeom',

     q     'error reading boozmn file in call to read_boozer_file',

     r     'error opening neo code input file neo_in',

     s     'trouble running eq3d in chisq_orbit',

     t     'trouble running mkjmc in chisq_orbit',

     u     'trouble running orbit in chisq_orbit',                                       !30

     v     'error opening orbsum in chisq_orbit',
     w     'error opening ft79jmc in chisq_dsubr',
     X     'error in chisq_vac_island',

     x     'error running xv3post in chisq_diagnostics',                              
     x     'error in v3post routine',
     y     'error opening v3post output file'
     Z    /

      ALLOCATE (fvec(nopt), diag(nvar), stat = info)
      IF (info .ne. 0) STOP stop_string

      tol = 1.e-9_dp
      mode = 1
      niter = niter_opt
      IF (lone_step) niter = 1

!
!     print description of independent variables at top of screen
!
      IF (myid .eq. master) THEN
         WRITE (6, '(60("="),/,a,/,60("="),/,a)') describe_string,
     1         ' VAR #   TYPE'
         DO info = 1, nvar
            WRITE(6, '(i5,2x,a)') info, var_descript(info)
         END DO
         WRITE (6, *)
         WRITE (6, '(60("="),/,a,/,60("="),/,a)') constraint_string,
     1         ' VAR #   CONSTRAINT'
         DO info = 1, nopt
            WRITE(6, '(i5,2x,a)') info, chisq_descript(info)
         END DO
         WRITE (6, *)
         CALL safe_open(iunit_vars,ier,'var_labels',
     1                  'replace','formatted')
         WRITE (iunit_vars,'(i5)') nvar
         DO info = 1, nvar
            WRITE(iunit_vars, '(i5,2x,a)') info, var_descript(info)
         END DO
         WRITE (iunit_vars, *)
         WRITE (iunit_vars,'(i5)') nopt
         DO info = 1, nopt
            WRITE(iunit_vars, '(i5,2x,a)') info, chisq_descript(info)
         END DO
         WRITE (iunit_vars, *)
         CLOSE(iunit_vars)
      END IF

      CALL flush(6)

      IF(NOPT_ALG .eq. 0) THEN
!DEC$ IF DEFINED (MPI_OPT)
        CALL lmdif1_mp (lsfun1, nopt, nvar, xc_opt, fvec, tol, epsfcn,
     1     niter, diag, mode, info, lwa)
!DEC$ ELSE
        CALL lmdif1 (lsfun1, nopt, nvar, xc_opt, fvec, tol, epsfcn,
     1     niter, diag, mode, info, lwa, num_processors,
     1     num_levmar_params)
!DEC$ ENDIF
      ELSE IF(nopt_alg .eq. 1) THEN
         CALL ga_driver (lsfun1, nopt, nvar, xc_opt, fvec, tol, epsfcn,
     1   niter, num_processors, seq_ext, info, lwa, lrestart )

      ELSE IF(nopt_alg .eq. 2) THEN
         CALL de_driver (lsfun1, nopt, nvar, xc_opt, fvec, tol, epsfcn,
     1   niter, num_processors, seq_ext, info, lwa, lrestart )

      ELSE
         IF (myid .eq. master) 
     1      WRITE(6,*) "NOPT_ALG = ", nopt_alg, "; unable to proceed"
         STOP
      ENDIF
      DEALLOCATE (fvec, diag)

      IF (myid.eq.master .and. info.lt.0) WRITE (*, '(/,1x,a,a)')
     1     'STELLOPT status: ',TRIM(info_array(-info))

      END SUBROUTINE run_optimizer
