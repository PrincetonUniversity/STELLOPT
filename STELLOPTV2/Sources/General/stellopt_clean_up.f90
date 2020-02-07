!-----------------------------------------------------------------------
!     Subroutine:    stellopt_clean_up
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/26/2012
!     Description:   This subroutine handles cleaning up the function
!                    calls.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_clean_up(ncnt,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE safe_open_mod, ONLY: safe_open
      USE mpi_inc
      USE vmec_input
      USE fdjac_mod, ONLY: flag_singletask, flag_cleanup, &
                           JAC_CLEANUP => flag_cleanup_jac, &
                           LEV_CLEANUP => flag_cleanup_lev, &
                           LBFGSB_CLEANUP => flag_cleanup_lbfgsb
      USE gade_mod, ONLY: GADE_CLEANUP, PSO_CLEANUP
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in)    :: ncnt
      INTEGER, INTENT(inout) :: iflag
!----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER ::  ier, ik, iunit, ctype, temp_max, m, n
      CHARACTER(len = 256)   :: temp_str
      REAL(rprec), ALLOCATABLE :: fvec_temp(:)
      
      INTEGER, PARAMETER :: JUST_INPUT  = -110
      INTEGER, PARAMETER :: LAST_GO     = -500
      
      LOGICAL, PARAMETER :: lkeep_extra = .TRUE.  ! For debugging the code will keep the _optXXX files
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      ctype = iflag
      iflag = 0
      iunit = 10
      iunit_out = 12
      ier = 0
      IF (ctype == PSO_CLEANUP) THEN
         ! Write Input File
         IF (ncnt /= 1) CALL stellopt_write_inputfile(ncnt,.false.)
         ! If Keeping minimums
         IF (lkeep_mins) THEN
            ! Write EQULIBRIUM File
            WRITE(temp_str,'(i5.5)') ncnt
            proc_string = TRIM(id_string) // '.' // TRIM(ADJUSTL(temp_str))
            CALL stellopt_write_eqfile
            CALL stellopt_write_auxfiles
         END IF
         ! Now open the Output file
         ALLOCATE(fvec_temp(mtargets))
         CALL safe_open(iunit_out,iflag,TRIM('stellopt.'//TRIM(id_string)),'unknown','formatted',ACCESS_IN='APPEND')
         iflag = 1
         IF (ncnt == 0) WRITE(iunit_out,'(A,1X,F5.2)') 'VERSION',STELLOPT_VERSION
         WRITE(iunit_out,'(A,1X,I5.5)') 'ITER',ncnt
         CALL stellopt_load_targets(mtargets,fvec_temp,iflag,ncnt)          ! Count
         !WRITE(iunit_out,'(A,2(2X,I8))') 'TARGETS ',mtargets,1
         !WRITE(iunit_out,'(A)') 'TARGETS'
         !WRITE(iunit_out,'(ES22.12E3)') targets(1:mtargets)
         !WRITE(iunit_out,'(A,2(2X,I8))') 'SIGMAS ',mtargets,1
         !WRITE(iunit_out,'(A)') 'SIGMAS'
         !WRITE(iunit_out,'(ES22.12E3)') sigmas(1:mtargets)
         !WRITE(iunit_out,'(A,2(2X,I8))') 'VALS ',mtargets,1
         !WRITE(iunit_out,'(A)') 'VALUES'
         !WRITE(iunit_out,'(ES22.12E3)') vals(1:mtargets)
         CLOSE(iunit_out)
         DEALLOCATE(fvec_temp)

      ELSE IF ((ctype == LEV_CLEANUP) .or. (ctype == GADE_CLEANUP)) THEN
          IF (ncnt /= 1 .or. ctype == GADE_CLEANUP) CALL stellopt_write_inputfile(ncnt,.false.)
          ! Overwrite the restart file
          proc_string = 'reset_file'
          CALL stellopt_write_eqfile
          ! Keep minimum states
          IF (lkeep_mins) THEN
             WRITE(temp_str,'(i5.5)') ncnt
             proc_string = TRIM(id_string) // '.' // TRIM(ADJUSTL(temp_str))

             CALL stellopt_write_eqfile
             CALL stellopt_write_auxfiles
          END IF
          ! Now open the Output file
          ALLOCATE(fvec_temp(mtargets))

          !WRITE COIL KNOTS, CONTROL POINTS TO FILE
          IF (ANY(lcoil_spline)) THEN
             CALL safe_open(iunit_out,iflag,TRIM('cbspline.'//TRIM(id_string)),'unknown','formatted',ACCESS_IN='APPEND')
             IF (ncnt == 0) THEN
                WRITE(iunit_out,'(A)') 'COIL KNOTS'
                DO n = LBOUND(lcoil_spline,DIM=1), UBOUND(lcoil_spline,DIM=1)
                   IF (ANY(lcoil_spline(n,:))) THEN
                      WRITE(iunit_out,'(2X,A,2X,I5.5)') 'COIL', n
                           ik = coil_nctrl(n) + 4
                      WRITE(iunit_out,"(4X,'k =',4(2X,ES22.12E3))") (coil_splinesx(n,m), m = 1, ik)
                   ENDIF
                END DO !n
             ENDIF
             WRITE(iunit_out,'(A,1X,I5.5)') 'COIL CTRL PTS, ITER',ncnt
             DO n = LBOUND(lcoil_spline,DIM=1), UBOUND(lcoil_spline,DIM=1)
                IF (ANY(lcoil_spline(n,:))) THEN
                   WRITE(iunit_out,'(2X,A,2X,I5.5)') 'COIL', n
                         ik = coil_nctrl(n)
                   IF (lwindsurf) THEN
                      WRITE(iunit_out,"(4X,'u =',4(2X,ES22.12E3))") (coil_splinefx(n,m), m = 1, ik)
                      WRITE(iunit_out,"(4X,'v =',4(2X,ES22.12E3))") (coil_splinefy(n,m), m = 1, ik)
                   ELSE
                      WRITE(iunit_out,"(4X,'x =',4(2X,ES22.12E3))") (coil_splinefx(n,m), m = 1, ik)
                      WRITE(iunit_out,"(4X,'y =',4(2X,ES22.12E3))") (coil_splinefy(n,m), m = 1, ik)
                      WRITE(iunit_out,"(4X,'z =',4(2X,ES22.12E3))") (coil_splinefz(n,m), m = 1, ik)
                   END IF !lwindsurf
                END IF
             END DO !n
             CLOSE(iunit_out)
          END IF !lcoil_spline

          CALL safe_open(iunit_out,iflag,TRIM('stellopt.'//TRIM(id_string)),'unknown','formatted',ACCESS_IN='APPEND')
          iflag = 1
          IF (ncnt == 0) WRITE(iunit_out,'(A,1X,F5.2)') 'VERSION',STELLOPT_VERSION
          WRITE(iunit_out,'(A,1X,I5.5)') 'ITER',ncnt
          CALL stellopt_load_targets(mtargets,fvec_temp,iflag,ncnt)          ! Count
          !WRITE(iunit_out,'(A,2(2X,I8))') 'TARGETS ',mtargets,1 
          !WRITE(iunit_out,'(A)') 'TARGETS'
          !WRITE(iunit_out,'(ES22.12E3)') targets(1:mtargets)
          !WRITE(iunit_out,'(A,2(2X,I8))') 'SIGMAS ',mtargets,1
          !WRITE(iunit_out,'(A)') 'SIGMAS'
          !WRITE(iunit_out,'(ES22.12E3)') sigmas(1:mtargets)
          !WRITE(iunit_out,'(A,2(2X,I8))') 'VALS ',mtargets,1
          !WRITE(iunit_out,'(A)') 'VALUES'
          !WRITE(iunit_out,'(ES22.12E3)') vals(1:mtargets)
          CLOSE(iunit_out)
          DEALLOCATE(fvec_temp)
      ELSE IF (ctype == JAC_CLEANUP) THEN
      ELSE IF (ctype == JUST_INPUT) THEN
         ! Write the input file
         CALL stellopt_write_inputfile(ncnt,.false.)
      ELSE IF (ctype == LAST_GO) THEN
         CALL stellopt_write_inputfile(ncnt,.true.)
         ! Now open the Output file
         ALLOCATE(fvec_temp(mtargets))
         CALL safe_open(iunit_out,iflag,TRIM('stellopt.'//TRIM(id_string)),'unknown','formatted',ACCESS_IN='APPEND')
         iflag = 1
         WRITE(iunit_out,'(A)') 'ITER MIN'
         ik = 999
         CALL stellopt_load_targets(mtargets,fvec_temp,iflag,ik)    ! Count
         !WRITE(iunit_out,'(A,2(2X,I8))') 'TARGETS ',mtargets,1
         !WRITE(iunit_out,'(A)') 'TARGETS'
         !WRITE(iunit_out,'(ES22.12E3)') targets(1:mtargets)
         !WRITE(iunit_out,'(A,2(2X,I8))') 'SIGMAS ',mtargets,1
         !WRITE(iunit_out,'(A)') 'SIGMAS'
         !WRITE(iunit_out,'(ES22.12E3)') sigmas(1:mtargets)
         !WRITE(iunit_out,'(A,2(2X,I8))') 'VALS ',mtargets,1
         !WRITE(iunit_out,'(A)') 'VALUES'
         !WRITE(iunit_out,'(ES22.12E3)') vals(1:mtargets)
         CLOSE(iunit_out)
         DEALLOCATE(fvec_temp)
         IF (.not.lkeep_extra) THEN
            temp_max = max(nvars,numprocs)
            DO ik = 0, temp_max, numprocs-1
               WRITE(temp_str,'(i5)') ik
               proc_string = TRIM(id_string) // '_opt' // TRIM(ADJUSTL(temp_str))
               OPEN(iunit,FILE=TRIM('answers_plot.'//TRIM(proc_string)),STATUS='unknown',IOSTAT=ier)
               IF (ier == 0) CLOSE(iunit,STATUS='delete')
               OPEN(iunit,FILE=TRIM('answers.'//TRIM(proc_string)),STATUS='unknown',IOSTAT=ier)
               IF (ier == 0) CLOSE(iunit,STATUS='delete')
               OPEN(iunit,FILE=TRIM('boozmn_'//TRIM(proc_string)//'.nc'),STATUS='unknown',IOSTAT=ier)
               IF (ier == 0) CLOSE(iunit,STATUS='delete')
               OPEN(iunit,FILE=TRIM('diagno_bth.'//TRIM(proc_string)),STATUS='unknown',IOSTAT=ier)
               IF (ier == 0) CLOSE(iunit,STATUS='delete')
               OPEN(iunit,FILE=TRIM('diagno_flux.'//TRIM(proc_string)),STATUS='unknown',IOSTAT=ier)
               IF (ier == 0) CLOSE(iunit,STATUS='delete')
               OPEN(iunit,FILE=TRIM('jBbs.'//TRIM(proc_string)),STATUS='unknown',IOSTAT=ier)
               IF (ier == 0) CLOSE(iunit,STATUS='delete')
               OPEN(iunit,FILE=TRIM('mercier.'//TRIM(proc_string)),STATUS='unknown',IOSTAT=ier)
               IF (ier == 0) CLOSE(iunit,STATUS='delete')
               OPEN(iunit,FILE=TRIM('jxbout_'//TRIM(proc_string)//'.nc'),STATUS='unknown',IOSTAT=ier)
               IF (ier == 0) CLOSE(iunit,STATUS='delete')
               OPEN(iunit,FILE=TRIM('wout_'//TRIM(proc_string)//'.nc'),STATUS='unknown',IOSTAT=ier)
               IF (ier == 0) CLOSE(iunit,STATUS='delete')
            END DO
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_clean_up
