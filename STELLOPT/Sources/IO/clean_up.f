      SUBROUTINE clean_up (nvar, lastgo, iflag)
      USE optim
      USE vmec_input, ONLY: raxis_cc, raxis_cs, rbc, zbs,
     1                      zaxis_cs, zaxis_cc, mgrid_file,
     2                      lfreeb
      USE safe_open_mod
      USE system_mod
      USE mpi_params
      USE gade_mod, ONLY: npopsiz
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include 'mpif.h'                                                   !MPI
!DEC$ ENDIF
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: nvar, iflag
      LOGICAL :: lastgo
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: istore = 3
      CHARACTER(len=*), PARAMETER :: all_files="AllFiles"
      INTEGER :: iunit, icount, istat, istat2, itemp, icmax, 
     1           nc_ext, iteration
      REAL(rprec) :: aspect, chisq_tot
      LOGICAL :: lexist, lnew_min
      CHARACTER(len=200) :: temp_input, save_min_ext, testfile
      CHARACTER(len=25) :: ch1, ch2, ch3
!-----------------------------------------------
!
!     Note: be sure to catenate -iflag output here, too
!    (i starts at 0 for that case in do loop below)
!
      iflag = 0
      save_min_ext = " "
      lnew_min = .false.

      ! Only the master process executes cleanup
      IF (myid .ne. master) GOTO 2000

      INQUIRE (file=TRIM(output_file), exist=lexist, iostat=istat)
      IF (lexist) THEN
         OPEN(unit=iunit_opt, file=TRIM(output_file), status='old',
     1        position='append', action='readwrite', iostat=istat)
      ELSE
         OPEN(unit=iunit_opt, file=TRIM(output_file), status='new',
     1        iostat=istat)
      END IF

      IF (istat .ne. 0) THEN
         iflag = -7
         GOTO 2000
      END IF

      iunit = iunit_opt+1
      chisq_tot = chisq_min
      icmax = MAX(nvar, num_levmar_params)
      IF (nopt_alg .gt. 0) icmax = MAX(icmax, npopsiz)

      findmin: DO icount = 0, icmax
         IF (lniter1 .and. icount.ne.0) EXIT
!        BE SURE THIS OPT_EXT IS CONSISTENT WITH WHAT LSFUN1 WRITES         
         WRITE (temp_input,'(i5)') icount
         opt_ext =
     1      TRIM(seq_ext) // '_opt' // TRIM(ADJUSTL(temp_input))
         temp_input = 'output.' // TRIM(opt_ext)
         IF (lniter1) save_min_ext = TRIM(opt_ext)

         CALL safe_open(iunit, istat, temp_input, 'old',
     1             'formatted')

         IF (istat .ne. 0) CYCLE
         
         DO WHILE (istat .eq. 0)
            READ  (iunit,'(a)', iostat=istat) temp_input
            WRITE (iunit_opt, '(a)') TRIM(temp_input)
            IF (INDEX(temp_input,'Aspect Ratio =').ne.0 .and.
     1          INDEX(temp_input,'Chi-Sq =').ne.0) EXIT
         END DO

         IF (istat .eq. 0) THEN
            READ  (temp_input,'(a16,f12.5,a10,f12.5,a23,i6)',
     1             iostat=istat)
     2             ch1, aspect, ch2, chisq_tot, ch3, iteration

!           Save in output file, delete temp file
            DO
               READ  (iunit, '(a)', iostat=istat) temp_input
               IF (istat .ne. 0) EXIT
               WRITE (iunit_opt,'(a)') TRIM(temp_input)
            END DO

!
!        STORE PARAMETERS CORRESPONDING TO MINIMUM CHI-SQ STATE
!        IN INPUT.OPT_EXT.MIN FILE. ALSO STORE WOUT.OPT_EXT.MIN FOR
!        LRESET_OPT=F FAST RESTART OPTION. SAVE ITERATION IN ITER_min
!

            IF (chisq_tot .lt. chisq_min) THEN
               chisq_min = chisq_tot
               iter_min = iteration
               save_min_ext = TRIM(opt_ext)
               lnew_min = .true.
            END IF

         END IF
         
         CLOSE (iunit, iostat=istat) !PPPL

!         CLOSE (iunit, status='delete', iostat=istat) !ORNL

      END DO findmin

!
!     STORE MIN STATE FILES, BUT ONLY IF ONE WAS FOUND
!
      itemp = iunit+1

      IF (lnew_min) THEN
         min_count = min_count + 1
         ch1 = ' '

!
!     Save the input.*_opt?????  file ->  input.*.min
!
         temp_input = copy // "input." // TRIM(save_min_ext) //
     1                " " // TRIM(min_input_file)
         CALL system(temp_input,istat)

!
!     If keeping all the .min files (for run history reconstruction), 
!     store this one in sequence and mark it with label 'min_count'
!

         IF (lkeep_mins) THEN
            ch1 = ' '
!            WRITE(ch1,'(".",i4.4)') min_count
            write(ch1,'(".",i6.6)') iter_min
            ! VMEC INPUT
!            temp_input = copy // "input." // TRIM(save_min_ext) //
!     1                   " " // TRIM(min_input_file) // TRIM(ch1)
            temp_input = copy //"  input."//trim(save_min_ext) // " " 
     1                  // "input." // trim(seq_ext) // trim(ch1)
            CALL system(temp_input,istat)
            ! ERROR Estimate
            temp_input = copy //"  errvals."//trim(save_min_ext)
     1                // " " //"errvals." // trim(seq_ext) // trim(ch1)
            INQUIRE(FILE="errvals."//trim(save_min_ext),
     1              EXIST=lexist)
            IF (lexist) CALL system(temp_input,istat)
!ADDED BY EAL  VMEC WOUT
!DEC$ IF DEFINED (NETCDF)
!            temp_input = copy // "wout_" // TRIM(save_min_ext) //".nc"
!     1               //  " " // TRIM(min_wout_file) // TRIM(ch1) //".nc"
            temp_input = copy // "wout_" // TRIM(save_min_ext) //".nc"
     1              //  " wout_" // TRIM(seq_ext) // TRIM(ch1) //".nc"
            
!DEC$ ELSE
!            temp_input = copy // "wout." // TRIM(save_min_ext) //
!     1                 " " // TRIM(min_wout_file) // TRIM(ch1)
            temp_input = copy // "wout." // TRIM(save_min_ext) //
     1                 " wout." // TRIM(seq_ext) // TRIM(ch1)  
!DEC$ ENDIF
            CALL system(temp_input,istat)
            ! VMEC MERCIER
!DEC$ IF DEFINED (NETCDF)
            temp_input =copy//"  mercier_"//trim(save_min_ext)//".nc  " 
     1             // "mercier_" // trim(seq_ext) // trim(ch1) // ".nc"
            INQUIRE(FILE="mercier_"//trim(save_min_ext)//".nc",
     1              EXIST=lexist)
!DEC$ ELSE
            temp_input = copy //"  mercier."//trim(save_min_ext) // " " 
     1                  // "mercier." // trim(seq_ext) // trim(ch1)
            INQUIRE(FILE="mercier."//trim(save_min_ext),EXIST=lexist)
!DEC$ ENDIF
            IF (lexist) CALL system(temp_input,istat)
            ! VMEC THREED1
            temp_input = copy //"  threed1."//trim(save_min_ext) // " " 
     1                 // "threed1." // trim(seq_ext) // trim(ch1)
            CALL system(temp_input,istat)
            ! VMEC DCON
            temp_input = copy //"  dcon_"//trim(save_min_ext) // ".txt " 
     1                  // "dcon_" // trim(seq_ext) // trim(ch1) //
     2                  ".txt"
            INQUIRE(FILE="dcon_"//trim(save_min_ext) // ".txt",
     1              EXIST=lexist)
            IF (lexist) CALL system(temp_input,istat)
            ! V3 FIT
            IF (lv3post) THEN
                temp_input = 
     1             copy // "chisq_diagno." // TRIM(save_min_ext) //
     1                   " " // TRIM(min_diagno_file) // TRIM(ch1)
                CALL system(temp_input,istat)
            ENDIF
!           DIAGNO MAGNETIC DIAGNOSTIC TARGETING
            IF (ldiagno_opt) THEN
               temp_input = copy //"  diagno_bth."//trim(save_min_ext)
     1        // " " //"diagno_bth." // trim(seq_ext) // trim(ch1)
               INQUIRE(FILE="diagno_bth."//trim(save_min_ext),
     1              EXIST=lexist)
               IF (lexist) CALL system(temp_input,istat)
               temp_input = copy //"  diagno_flux."//trim(save_min_ext)
     1        // " " //"diagno_flux." // trim(seq_ext) // trim(ch1)
               INQUIRE(FILE="diagno_flux."//trim(save_min_ext),
     1              EXIST=lexist)
               IF (lexist) CALL system(temp_input,istat)
               temp_input = copy //"  diagno_seg."//trim(save_min_ext)
     1        // " " //"diagno_seg." // trim(seq_ext) // trim(ch1)
               INQUIRE(FILE="diagno_seg."//trim(save_min_ext),
     1              EXIST=lexist)
               IF (lexist) CALL system(temp_input,istat)
!               temp_input = copy //"  mag_diags."//trim(save_min_ext)
!     1        // " " //"mag_diags." // trim(seq_ext) // trim(ch1)
!               CALL system(temp_input,istat)
            END IF
!           PRESSURE PROFILE TARGETING
            IF ( ANY(sigma_p_prof .lt. bigno)
     1          .or. ANY(sigma_ne_prof < 1.0E30) 
     2          .or. ANY(sigma_te_prof < bigno) 
     2          .or. ANY(sigma_ti_prof < bigno) ) THEN
               temp_input = copy //"  p_prof."//trim(save_min_ext)
     1            // " " //"p_prof." // trim(seq_ext) // trim(ch1)
               CALL system(temp_input,istat)
               IF (isote) THEN
                  temp_input = copy //"  te_prof."//trim(save_min_ext)
     1               // " " //"te_prof." // trim(seq_ext) // trim(ch1)
                  CALL system(temp_input,istat)
               END IF
            END IF
!           MSE DIAGNOSTIC TARGETING
            IF (ANY(sigma_mse_pol .lt. bigno)) THEN
               temp_input = copy //"  mse_prof."//trim(save_min_ext)
     1            // " " //"mse_prof." // trim(seq_ext) // trim(ch1)
               CALL system(temp_input,istat)
            END IF
         ENDIF	!  lkeep_mins

!
!     Construct the .min_fp file containing the new boundary and modified axis
!     This is used to restart a free-boundary (fb) run in fixed-boundary mode
!
         temp_input = 'input.' // TRIM(save_min_ext)

!        Open temp_input (minimum input file), load modules for writing, and delete temp_input file
         CALL read_input_opt (temp_input, istat2)                       

         IF (istat2 .eq. 0) THEN
!           Read in axis and boundary data corresponding to this wout file prior to writing out
!           For free-bdy run, set boundary values to optimized state for subsequent optimization
            CALL read_wout_opt(istore, save_min_ext, istat, istat2)
            CALL write_indata (TRIM(min_input_file) // "_fb", iflag)
            IF (iflag .ne. 0) GOTO 2000

!           Blend old and new axis positions
!           This should be optional (SPH-01/10/01). Does NOT work well to update(blend) axis
!           UNLESS low enough ftol values are used...otherwise turn this option OFF
C            raxis_old(:,1) = (3*raxis_old(:,1) + raxis_cc(:))/4
C            raxis_old(:,2) = (3*raxis_old(:,2) + raxis_cs(:))/4
C            zaxis_old(:,1) = (3*zaxis_old(:,1) + zaxis_cs(:))/4
C            zaxis_old(:,2) = (3*zaxis_old(:,2) + zaxis_cc(:))/4

!           Save (with min_input_file extension) all files for this minimum state
!           Note: wout MAY have the .nc extension at end, so changed *.ext to *ext* 
!           to catch this case in temp_input below
!
            temp_input = list // "*" // TRIM(save_min_ext) 
     1                // "* > min_file"
            CALL system(temp_input, istat)
            CALL safe_open(itemp, istat, "min_file", "old", "formatted")
            DO WHILE (istat .eq. 0)
               READ (itemp, '(a)', END=50, iostat=istat) testfile
               icount = INDEX(testfile, TRIM(save_min_ext)) - 1
               nc_ext = icount + 1 + LEN_TRIM(save_min_ext)              !possible "nc" extension
               IF (nc_ext.le.LEN_TRIM(testfile) .and.
     1             testfile(nc_ext:nc_ext+2).ne.'.nc') CYCLE             !avoid copying,e.g., _opt10 instead of _opt1
               temp_input = move // TRIM(testfile) // " " 
     1                       // testfile(1:icount) // TRIM(min_ext)
     2                       // testfile(nc_ext:)
               CALL system(temp_input, istat)
            END DO
 50         CONTINUE
            CLOSE (itemp, status='delete')
         ELSE IF (save_min_ext(1:1) .ne. " ") THEN
            iflag = -7
            PRINT *,'I/O Error reading file ', TRIM(temp_input),
     1              ' in STELLOPT routine CLEAN_UP, istat = ', istat2
            GOTO 2000
         END IF
         
      END IF              !lnew_min
!
!     Delete input, wout files that are NOT the minimum state (i.e, without .min extension)
!     Keep output file record and mgrid file and coils.dot file
!
      IF (.not.lniter1) THEN
!
!     Remove the _opt files
!        1) rm -f *_optX? >& /dev/null (remove groups of 10)
!  
         do icount = 10, icmax, 10
            ! Get the _opt files
            write (temp_input,'(i5)') icount/10

            temp_input = remove //' *_opt'//trim(adjustl(temp_input))
     1                      //'? >& /dev/null'
            call system(temp_input, istat)
            
            ! Get the _opt.nc files
            write (temp_input,'(i5)') icount/10

            temp_input = remove //' *_opt'//trim(adjustl(temp_input))
     1                      //'?.nc >& /dev/null'
            call system(temp_input, istat)
            
            ! Get the _opt.txt files
            write (temp_input,'(i5)') icount/10

            temp_input = remove //' *_opt'//trim(adjustl(temp_input))
     1                      //'?.txt >& /dev/null'
            call system(temp_input, istat)
         enddo
!        2) rm -f *_opt? > & /dev/null (pick up any stray opts)
         call system('rm -f *_opt? >& /dev/null', istat)
         call system('rm -f *_opt?.nc >& /dev/null',istat)
         call system('rm -f *_opt?.txt >& /dev/null',istat)
         call system('rm -f *_opt?? >& /dev/null', istat)
         call system('rm -f *_opt??.nc >& /dev/null',istat)
         call system('rm -f *_opt??.txt >& /dev/null',istat)
         call system('rm -f *_opt??? >& /dev/null', istat)
         call system('rm -f *_opt???.nc >& /dev/null',istat)
         call system('rm -f *_opt???.txt >& /dev/null',istat)
!
!        Avoid wiping out files necessary for continuing the run
!        
         IF (lastgo) THEN
            temp_input = list // "> " // all_files
            CALL system(temp_input)
            CALL safe_open(itemp, istat, all_files, "old", "formatted")
            DO WHILE (istat .eq. 0)
               READ (itemp, '(a)', END=60) testfile
            IF ((nopt_alg.eq.1 .and. INDEX(testfile,"ga_restart").gt.0)
     1     .or. (nopt_alg.eq.2 .and. INDEX(testfile,"de_restart").gt.0)
     2     .or. (INDEX(testfile, TRIM(output_file)) .gt. 0)
     3     .or. (INDEX(testfile, TRIM(mgrid_file)) .gt. 0)
     4     .or. (INDEX(testfile,"param") .gt. 0)
     5     .or. (TRIM(testfile) .eq. TRIM(all_files))
     6     .or. (INDEX(testfile,".dat") .gt. 0)
     7     .or. (INDEX(testfile,"coils.").gt.0 .and. lvac_opt)
     8     .or. (INDEX(testfile,"ft5") .gt. 0)
     9     .or. (index(testfile,"diagno.control") .gt. 0)) CYCLE

            IF (lv3post .and. ((INDEX(testfile,TRIM(v3post_in)) .gt. 0)
     1     .or. (INDEX(testfile,"v3post.in") .gt. 0)
     1     .or. (INDEX(testfile,"crfun_") .gt. 0)
     2     .or. (INDEX(testfile,".LIST")  .gt. 0)
     2     .or. (INDEX(testfile,"chisq_diagno.")  .gt. 0)
     2     .or. (INDEX(testfile,"diagnostics.") .gt. 0))) CYCLE

               icount = INDEX(TRIM(testfile), TRIM(min_ext),
     1                         back=.true.)
               IF (icount .gt. 0) CYCLE
               temp_input = remove // TRIM(testfile)
               CALL system(temp_input, istat)
            END DO
 60         CLOSE (itemp, status='delete', iostat=istat)
         END IF ! Lastgo
      END IF

      CLOSE (iunit_opt, iostat=istat)

 2000 CONTINUE
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BCAST(chisq_min, 1, MPI_REAL8, master, MPI_COMM_WORLD,
     1     istat)
      IF (istat .ne. 0) STOP 'MPI_BCAST error in STELLOPT CLEAN_UP'
      CALL MPI_BCAST(iflag, 1, MPI_INTEGER, master, MPI_COMM_WORLD,
     1     istat)
      IF (istat .ne. 0) STOP 'MPI_BCAST error in STELLOPT CLEAN_UP'

!     Load boundary values for free-boundary
      IF (lfreeb) THEN
         CALL MPI_BCAST(rbc, SIZE(rbc), MPI_REAL8, master, 
     1                  MPI_COMM_WORLD, istat)
         CALL MPI_BCAST(zbs, SIZE(zbs), MPI_REAL8, master, 
     1                  MPI_COMM_WORLD, istat)
         CALL MPI_BCAST(raxis_cc, SIZE(raxis_cc), MPI_REAL8, master, 
     1                  MPI_COMM_WORLD, istat)
         CALL MPI_BCAST(zaxis_cs, SIZE(zaxis_cs), MPI_REAL8, master, 
     1                  MPI_COMM_WORLD, istat)
         CALL MPI_BCAST(raxis_old, SIZE(raxis_old), MPI_REAL8, master, 
     1                  MPI_COMM_WORLD, istat)
         CALL MPI_BCAST(zaxis_old, SIZE(zaxis_old), MPI_REAL8, master, 
     1                  MPI_COMM_WORLD, istat)
      END IF

!DEC$ ENDIF
      END SUBROUTINE clean_up
