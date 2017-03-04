      SUBROUTINE clean_up (nvar, iflag)
      USE control_mod
      USE Vname
      USE safe_open_mod
      USE system_mod, ONLY: system
      USE mpi_params
      USE gade_mod, ONLY: npopsiz
      IMPLICIT NONE
#if defined(MPI_OPT)
      include 'mpif.h'                                                   !MPI
#endif
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nvar, iflag
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
#if defined(WIN32)
      CHARACTER(LEN=*), PARAMETER :: list = "dir /B "
      CHARACTER(LEN=*), PARAMETER :: move = "ren "
      CHARACTER(LEN=*), PARAMETER :: remove = "del /Q "
#else
      CHARACTER(LEN=*), PARAMETER :: list = "/bin/ls -1 "
      CHARACTER(LEN=*), PARAMETER :: move = "/bin/mv -f "
      CHARACTER(LEN=*), PARAMETER :: remove = "/bin/rm -f "
#endif
      INTEGER :: iunit, icount, istat, itemp, icmax, itemp1
      REAL(rprec) ::  chisq_tot
      LOGICAL :: lnew_min
      CHARACTER(LEN=200) :: temp_input, save_min_ext, testfile, opt_ext
      REAL(rprec) :: rms_err, avg_err, max_err
C-----------------------------------------------

      iflag = 0
      save_min_ext = " "
      lnew_min = .false.

      IF (myid .NE. master) GOTO 2000

      iunit = 14
      chisq_tot = chisq_min
      icmax = MAX(nvar, num_levmar_params)
      IF (nopt_alg .gt. 0) icmax = MAX(icmax, npopsiz)

      findmin: DO icount = 0, icmax
         WRITE (temp_input,'(i5)') icount
         opt_ext =
     1      TRIM(extension)//'_opt'//TRIM(ADJUSTL(temp_input))
         temp_input = 'coil_targets.' // TRIM(opt_ext)

         CALL safe_open(iunit, istat, temp_input, 'old',
     1             'formatted')
        
         IF (istat .NE. 0) CYCLE

         DO WHILE (istat .eq. 0)
            READ  (iunit,'(a)', iostat=istat) temp_input
            IF (INDEX(temp_input,'chisq_tot') .NE. 0 ) EXIT
         END DO

         IF (istat .eq. 0) THEN
            READ (iunit,'(5e16.4)',iostat=istat)
     1             rms_err, avg_err, max_err, chisq_tot

            IF (chisq_tot .lt. chisq_min) THEN
               chisq_min = chisq_tot
               save_min_ext = TRIM(opt_ext)
               lnew_min = .true.
            END IF

         END IF

         CLOSE (iunit)

      END DO findmin

!
!     STORE MIN STATE FILES, BUT ONLY IF ONE WAS FOUND
!
      IF (lnew_min) THEN

         temp_input = list // "*." // TRIM(save_min_ext) // 
     1                " > min_file"
         CALL system(temp_input, istat)
         itemp = iunit+1; itemp1 = itemp+1
         CALL safe_open(itemp, istat, "min_file", "old", "formatted")
         DO WHILE (.TRUE.)
            READ (itemp, '(a)', END=50) testfile
            icount = INDEX(testfile, TRIM(save_min_ext))-1
            temp_input = testfile(1:icount) // TRIM(extension) // ".min"
            CALL safe_open(itemp1,istat,temp_input,"old","formatted")
            IF (istat .EQ. 0) CLOSE(itemp1,status='delete')
            temp_input=move // TRIM(testfile) // " " // TRIM(temp_input)
            CALL system(temp_input, istat)
         END DO
 50      CONTINUE
         CLOSE (itemp, status='delete')

      END IF              !lnew_min

      temp_input = list // " > all_files"
      CALL system(temp_input)
      itemp = iunit + 1
      CALL safe_open(itemp, istat, "all_files", "old", "formatted")
      DO WHILE (.TRUE.)
         READ (itemp, '(a)', END=60) testfile
         temp_input = ADJUSTL(testfile)
!      Avoid wiping out necessary files for continuing run
         IF ((nopt_alg.eq.1 .and. INDEX(testfile,"ga_restart").gt.0)
     1     .or. (nopt_alg.eq.2 .and. INDEX(testfile,"de_restart").gt.0)
     2     .or. (TRIM(testfile) .eq. "all_files")
     3     .or. (INDEX(testfile,".dat") .gt. 0)
     4     .or. (temp_input(1:4) == "wout")
     5     .or. (temp_input(1:5) == "bnorm") 
     6     .or. (INDEX(testfile,".new") .gt. 0)) CYCLE

         icount = INDEX(testfile, ".min", BACK=.TRUE.)
         IF (icount .eq. (LEN_TRIM(testfile) - 3)) CYCLE
         temp_input = remove // TRIM(testfile)
         CALL system(temp_input, istat)
      END DO

 60   CLOSE (itemp, status='delete', iostat=istat)

 2000 CONTINUE

!#if defined(MPI_OPT)
!      CALL MPI_BCAST(chisq_min, 1, MPI_REAL8, master, MPI_COMM_WORLD,
!     1     istat)
!      IF (istat .ne. 0) STOP 'MPI_BCAST error in COILOPT CLEAN_UP'
!      CALL MPI_BCAST(iflag, 1, MPI_INTEGER, master, MPI_COMM_WORLD,
!     1     istat)
!      IF (istat .ne. 0) STOP 'MPI_BCAST error in COILOPT CLEAN_UP'
!#endif
      END SUBROUTINE clean_up
