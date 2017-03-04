      SUBROUTINE optimize(in_file)
      USE optim
      USE boozer_params, ONLY: xm_bdy, xn_bdy, rmnc_bdy, zmns_bdy,
     1                         rmnc_opt, zmns_opt, lmns_opt
      USE vparams, ONLY: zero
      USE system_mod
      USE legendre_params
      USE mpi_params
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include 'mpif.h'                                                   !MPI
!DEC$ ENDIF
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      CHARACTER(LEN=*) :: in_file
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: k, nvar, nopt, lwa, istat, info, itry
      REAL(rprec), DIMENSION(nxc) :: xc_opt
      REAL(rprec) :: tstart, tstop
      CHARACTER(len=64) :: var_descript(nxc)
      CHARACTER(len=256) :: scratch_dir
      CHARACTER(len=256) :: temp
      LOGICAL :: lprint
      INTEGER, EXTERNAL :: getcwd
C-----------------------------------------------
      iter_min = 0

      CALL get_extension(in_file)
      scratch_dir = "stellopt" // "_" // TRIM(seq_ext)
!
!     DO ALL CALCULATIONS IN SCRATCH DIRECTORY
!
      IF(myid .eq. master) THEN                                          !START MPI
         temp = remd // TRIM(scratch_dir)                                !produces file find error message on PC
         CALL system(temp)
         temp = makd // TRIM(scratch_dir)
         CALL system(temp)
         itry = 1
 321     k = chdir(scratch_dir)
         IF (itry .lt. 200 .and. k .ne. 0) THEN
            WRITE(6,*) 'MYID:',myid,' CHDIR:',TRIM(scratch_dir),itry
            itry = itry + 1
            GOTO 321
         END IF
         IF (k .ne. 0) THEN
            PRINT *,
     1      'ierr = ',k,': Unable to chdir to working directory '
     2       ,scratch_dir
             k = getcwd(scratch_dir)
             PRINT *,'myid = ',myid,'& presently, dir = ',scratch_dir
            STOP
         END IF

!
!     COPY INPUT FILE FROM ORIGINAL DIRECTORY
!
         temp = copy // ".." // sep // TRIM(in_file) // " ." // 
     1                          sep // TRIM(in_file)
         CALL system(temp)
      END IF                                                             !END MPI

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)                         !MPI
      IF (ierr_mpi .ne. 0) 
     1   STOP 'MPI_BARRIER error in STELLOPT OPTIMIZE'
      IF (myid .ne. master) THEN
         itry = 1
 322     k = chdir(scratch_dir)                       !MPI
         IF (itry .lt. 200 .and. k .ne. 0) THEN
            WRITE(6,*) 'MYID:',myid,' CHDIR:',TRIM(scratch_dir),itry
            itry = itry + 1
            GOTO 322
         END IF
         IF (k .ne. 0) THEN
            PRINT *,
     1      'ierr = ',k,': Unable to chdir to working directory '
     2       ,scratch_dir
             k = getcwd(scratch_dir)
             PRINT *,'myid = ',myid,'& presently, dir = ',scratch_dir
            STOP
         END IF
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)                         !MPI
      IF (ierr_mpi .ne. 0) 
     1   STOP 'MPI_BARRIER error in STELLOPT OPTIMIZE'
!DEC$ ENDIF
!
!     Initialize variables, READ input file. On EXIT, in_file is the input file name
!     Set up input, output file extensions
!
      xc_opt = zero
      CALL osetup(in_file)
!
!     SETUP OPTIMIZATION ROUTINE PARAMETERS
!
      CALL initialize_opt (xc_opt, var_descript, nvar, nopt)
      lwa = (nopt+5)*nvar + nopt
!
!     RUN OPTIMIZE ROUTINE
!
      IF (.not.lscale_only) THEN
         CALL second0 (tstart)

         CALL run_optimizer (xc_opt, var_descript, nopt, nvar, 
     1                       lwa, info)
      ELSE
         CALL write_indata (TRIM(min_input_file), info)
      END IF

      CALL second0 (tstop)

      DEALLOCATE (rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy)
      deallocate (rmnc_opt, zmns_opt, lmns_opt)
      DEALLOCATE (ns_surf, ns_booz, lneed_modB)
      DEALLOCATE (ns_ball, ns_neo)                                       !VMECCOBRA, NEO
      IF(l_legendre .and. lcur_prof_opt)                                 !LEGENDRE
     1  DEALLOCATE (a_leg, b_leg, a_leg_inv, b_leg_inv, tc)              !LEGENDRE
      IF(l_legendre .and. liota_prof_opt)                                !LEGENDRE
     1  DEALLOCATE (a_leg, b_leg, a_leg_inv, b_leg_inv, ti)              !LEGENDRE
      IF(l_legendre .and. lpres_prof_opt)                                !LEGENDRE
     1  DEALLOCATE (a_leg, b_leg, a_leg_inv, b_leg_inv, tm)              !LEGENDRE

      IF (myid .eq. master) THEN                                         !START MPI
         lprint = .not.lniter1
         IF (lprint) PRINT 90, chisq_min, iter_min, TRIM(min_input_file)
!DEC$ IF .NOT.DEFINED(MPI_OPT)
!
!        CLEAN UP SCRATCH FILES
!
         temp = remove // "fort.*"
         CALL system(temp)
!DEC$ ENDIF
         PRINT *, 'All input/output files for this run are stored in',
     1          ' the directory: '
         PRINT *, TRIM(scratch_dir)
!
!        APPEND TIME, CHI-SQ TO END OF OUTPUT FILE
!
         OPEN (unit=iunit_opt, file=output_file, status='old',
     1         position='append', action='readwrite', iostat=istat)
         IF (istat .eq. 0) THEN
            WRITE(iunit_opt, '(/,a,i4,a,i6)')' Number of parameters = ',
     1      nvar,' Number of constraints = ',nopt
            IF (lprint) WRITE (iunit_opt, 90)
     1          chisq_min, iter_min, TRIM(min_input_file)
            WRITE (iunit_opt, 100) tstop-tstart
         END IF

         CLOSE (iunit_opt, iostat=istat)

         PRINT  100, tstop - tstart
!
!        CHANGE FILE EXTENSIONS IN INPUT FILES, OUTPUT FILES TO ELIMINATE _OPT
!        (OTHERWISE, INPUT FILES WILL NOT RUN CORRECTLY WITH THE WOUT FILE, WHICH
!         HAS EITHER THE .MIN EXTENSION OR NONE AT ALL)
!
         temp = ""
         IF (lprint) temp = ".min"
         CALL clean_up_extensions (min_ext, temp)

!
!        CONCATENATE INFO FROM SEVERAL FILES INTO A DATA_SUMMARY FILE
!
         CALL data_harvest (min_ext, output_file, var_descript, 
     1                      nvar, nopt)

      ENDIF                                                              !END MPI

 90   FORMAT(/,' Minimum Chi-sq = ',1p,e10.3,' at iteration ',i6,/,
     1  ' Minimum Chi-sq state stored in file ',a)
 100  FORMAT(/,' SECONDS IN OPTIMIZATION CODE:  ',1p,e12.3,' SEC')

      END SUBROUTINE optimize
