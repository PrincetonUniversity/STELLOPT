      SUBROUTINE data_harvest (extension, output_file, var_descript,
     1                         nvar, nopt)
      USE stel_kinds
      USE safe_open_mod
      USE vmec_input, ONLY: lfreeb
      USE optim, ONLY: lextcur, nextcur_vmec, describe_string
      USE saddle_coils, ONLY: nsad_unique_coils
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in)          :: nvar, nopt
      CHARACTER(len=*), INTENT(in) :: extension, output_file, 
     1                                var_descript(nvar)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      CHARACTER(LEN=*), PARAMETER :: file_list(8) =  (/
     1  "output  ", "threed1 ", "nescin  ", "nescout ", "dkes_opt",
     1  "threed1 ", "ft6tpr  ", "ft6tpr2 "/)
      CHARACTER(LEN=*), PARAMETER :: header(8) = (/
     1   "OPTIMIZATION PARAMETERS      ",
     2   "GEOMETRIC/MAGNETIC PARAMETERS",
     3   "COIL (NESCOIL) PARAMETERS    ",
     4   "COIL (NESCOIL) PARAMETERS    ",
     5   "TRANSPORT (DKES) PARAMETERS  ",
     6   "COIL PARAMETERS              ",
     7   "KINK (N=1)                   ",
     8   "KINK (N=0)                   "/)
      CHARACTER(len=*), PARAMETER :: parse_coils(6) = (/
     1   "modular currents (amp)    ",
     1   "saddle coil-coil distance ",
     2   "saddle radius of curvature",
     3   "saddle coil lengths (m)   ",
     4   "minimum y - coils (m)     ",
     5   "min plasma-coil distance (" /)
      CHARACTER(LEN=*), PARAMETER :: fmt(5) = (/
     1   "(a)   ", "(a)   ","(1x,a)", "(1x,a)", "(a)   " /)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: icount, iread=8, istat, iscratch, iwrite, jndex,
     1      lcount, istat2, icoil, jcount, ireadc
      REAL (rprec) :: chisq_min, chisq, dum, iota(3), es(3), dum1, dum2
      REAL (rprec), ALLOCATABLE, DIMENSION(:,:) :: coil_info
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: extcur
      CHARACTER(LEN=30),  ALLOCATABLE, DIMENSION(:) :: curlabel
      LOGICAL :: lwrite, lrewind, lheader, lflush
      CHARACTER(LEN=200) :: line, filename, temp
C-----------------------------------------------
!
!     CONCATENATES INFORMATION FROM SEVERAL DIFFERENT FILES INTO A
!     SINGLE "DATA_SUMMARY" FILE
!
      chisq_min = HUGE (chisq_min)
      lwrite = .false.;   lrewind = .false.
      es = -1

      iscratch = iread+1
      jndex = INDEX(extension,".min",BACK=.TRUE.)
      IF (jndex .eq. 0) jndex = LEN_TRIM(extension)+1
      CALL safe_open (iscratch, istat, "NULL", 'scratch', 'formatted')
      IF (istat .ne. 0) RETURN
      iwrite = iscratch+1
      CALL safe_open (iwrite, istat, "data_summary" // "." //
     1      TRIM(extension(1:jndex-1)), 'replace', 'formatted')
      IF (istat .ne. 0) RETURN

      FILE_LOOP: DO icount = 1, SIZE(file_list)
         IF (icount .eq. 1) THEN
            filename = TRIM(output_file)
         ELSE
            filename = TRIM(file_list(icount)) //"."// TRIM(extension)
         END IF
         CALL safe_open(iread, istat, filename, 'old', 'formatted')
         IF (istat .ne. 0) CYCLE

         lcount = 0
         IF (icount.ne.4) lheader = .true.
         PARSE_LOOP: DO WHILE (istat .eq. 0)
            READ (iread, '(a)', END=200) line

!              parse line, depending on which one it is

            SELECT CASE (icount)
            CASE (1)
               jndex = INDEX (line, "Chi-Sq =")
               IF (jndex .ne. 0) THEN
                  lrewind = .true.
                  READ (line(jndex+8:), *) chisq
                  lwrite = (chisq .lt. chisq_min)
                  IF (lwrite) chisq_min = chisq
               END IF
               CASE (2)
               jndex = INDEX (line, "   IOTA   ")
               IF (jndex .ne. 0) THEN
                  DO lcount = 1,3
                     READ (iread, '(a)') line
                  END DO
                  lcount = 1
                  istat2 = 0
                  DO WHILE (istat2 .eq. 0)
                    READ (iread, '(a)') line
                    READ (line, *, iostat=istat2) dum1, dum, dum, dum2
                    IF ((lcount.eq.1 .and. dum1.ge.0._dp) .or.
     1                 (lcount.eq.2 .and. dum1.ge.0.5_dp) .or.
     2                 (lcount.eq.3 .and. dum1.ge.0.9999_dp)) THEN
                       es(lcount) = dum1
                       iota(lcount) = dum2
                       lcount = lcount+1
                    ELSE IF (lcount .gt. 3) THEN
                       EXIT
                    END IF
                  END DO
                  lcount = 0
               END IF
               jndex = INDEX (line, "Aspect Ratio")
               IF (jndex .ne. 0) THEN
                  lwrite = .true.; lrewind = .true.
               END IF
               IF (lwrite) lcount = lcount+1
               IF (lcount .gt. 20) lwrite = .false.             !Read 20 lines from file

               CASE (3)
               jndex = INDEX (line, "Coil-Plasma separation =")
               IF (jndex .ne. 0) THEN
                  READ (line(jndex+24:), *) chisq
                  WRITE(line,'(a,1p,e10.3)')
     1                  "Coil-Plasma Separation = ", chisq
                  lwrite = .true.;   lrewind = .true.
               END IF
               IF (lwrite) lcount = lcount+1
               IF (lcount .gt. 1) lwrite = .false.

               CASE (4)
               jndex = INDEX (line, "Complexity")
               IF (jndex .ne. 0) lwrite = .true.
               IF (lwrite) lcount = lcount+1
               IF (lcount .gt. 8) lwrite = .false.

               CASE (5)
               lwrite = .true.
               lcount = lcount+1
               IF (lcount .eq. 1) lrewind = .true.

               CASE (6)
               jndex = INDEX (line, "EXTERNAL CURRENTS")
               IF (lfreeb .and.  jndex.ne.0) THEN
                  jcount = 5*(nextcur_vmec/5) + 5      !!5 to a line in threed1
                  READ (iread, '(a)') line             !!ignore blank line
                  IF (.not. ALLOCATED(extcur))
     1                ALLOCATE(extcur(jcount), curlabel(jcount))
                  jcount = 1
                  istat2 = 0
                  DO WHILE (jcount.lt.nextcur_vmec .and. istat2.eq.0)
                     READ (iread, '(a)') line
                     READ (line, *, iostat=istat2)
     1                    (curlabel(icoil),icoil=jcount,jcount+4)
                     READ (iread, '(a)') line
                     READ (line, *, iostat=istat2)
     1                    (extcur(icoil),icoil=jcount,jcount+4)
                     jcount = jcount+5                 !!5 currents/labels per line
                  END DO !  DO WHILE
                  ! 09/08/11 - SAL
                  ! Added this section to properly read the last line of
                  ! EXTCUR values
                  jcount = MOD(nextcur_vmec,5)
                  !jcount = nextcur_vmec-5*(nextcur_vmec/5)  ! Get remaining datapoints
                  IF (jcount .gt. 0) THEN
                     READ (iread, '(a)') line
                     READ (line, *, iostat=istat2)
     1                    (curlabel(icoil),icoil=nextcur_vmec-jcount+1,
     2                     nextcur_vmec)
                     READ (iread, '(a)') line
                     READ (line, *, iostat=istat2)
     1                    (extcur(icoil),icoil=nextcur_vmec-jcount+1,
     2                     nextcur_vmec)
                  END IF
                  lwrite = .true.
                  istat = -1
               ENDIF !  IF(lfreeb) THEN

               CASE (7:8)
               jndex = INDEX (line, "EIGENVALUE FROM WP / WK")
               IF (jndex .ne. 0) THEN
                  istat=INDEX(line,'=')
                  IF(istat.ne.0) WRITE (iscratch, '(a,1x,a)')
     .             TRIM(line(1:istat)),TRIM(line(istat+1:))
               ENDIF


            END SELECT

            IF (lwrite) THEN
               IF (lrewind) THEN
                  REWIND (iscratch, iostat=istat)
                  lrewind = .false.
               END IF
               IF ((icount .le. 3) .or.
     1             (icount .eq. 4 .and. (lcount==1 .or. lcount==4 .or.
     2              lcount==5 .or. lcount==8)) .or.
     3             (icount .eq. 5 .and. (lcount>2)))
     4             WRITE (iscratch, fmt(icount)) TRIM(line)
            END IF

         END DO PARSE_LOOP

 200     CONTINUE

         CLOSE (iread)
         lwrite = .false.
         REWIND (iscratch, iostat=istat)
         IF (istat .eq. 0) THEN
            IF (lheader) WRITE (iwrite, '(1x,a)') header(icount)
!
!     WRITE NO. VARIABLES, PARAMETERS TO HEAD OF FILE
!     LOGIC HERE PACKS MULTIPLE COEFFIENCTS OF THE SAME TYPE INTO ONE LINE
!
            IF (icount .eq. 1) THEN
               WRITE (iwrite, '(a)') ' Case: ' // TRIM(extension)
               WRITE (iwrite, '(a,i4,a,i6,/)')
     1            ' Number of parameters: ', nvar,
     2            ' Number of constraints: ',nopt
               WRITE (iwrite, '(60("="),/,a,/,60("="),/,a)') 
     1            describe_string, ' VARIABLE ID    VARIABLE TYPE'
               temp = " "
               lflush = .false.
               DO jcount = 1, nvar
                 jndex = INDEX(var_descript(jcount),'(i=')
                 IF (jndex .gt. 1) THEN
                    line = var_descript(jcount)(1:jndex-1)
                    lflush = (temp .ne. " ")
                    IF (.not.lflush) THEN
                       temp = line
                       lcount = 0
                    END IF
                    IF (TRIM(line) .eq. TRIM(temp)) THEN
                       lcount = lcount + 1
                       IF (jcount .ne. nvar) CYCLE
                    END IF
                 END IF
                 IF (lflush) THEN               !Flush accumulated data
                    istat2 = jcount
                    IF ((jcount.eq.nvar) .and. (jndex.gt.0)) 
     1                     istat2 = jcount+1
                    WRITE (iwrite,'(i5,a,i5,6x,2a,i1,a,i5,a)') 
     1                  istat2-lcount,'-',istat2-1,
     2                  TRIM(temp),'(',1,'-',lcount,')'
                       lflush = .false.
                       temp = line
                       lcount = 1
                 END IF
                 IF (jndex .le. 1) THEN
                    WRITE(iwrite, '(i5,12x,a)') 
     1                    jcount,var_descript(jcount)
                    temp = " "
                 END IF
                 
               END DO
               WRITE (iwrite, *)
                 
            ELSE IF (icount .eq. 3) THEN
               lheader = .false.
            ELSE IF (icount .eq. 6) THEN
               GOTO 400
            END IF

              
            DO WHILE (istat .eq. 0)
               READ (iscratch, '(a)', iostat=istat, END=300) line
               IF (icount.eq.1 .and. LEN_TRIM(line).eq.0) EXIT
               WRITE (iwrite, '(a)') TRIM(line)
            END DO
 
 300        CONTINUE
 
            IF (icount .eq. 2) THEN
               WRITE (iwrite,'(a,2(f3.2,a),f4.2,a,1p,3e10.3)')
     1         ' iota(s=',es(1), ',', es(2), ',', es(3), ')  =      ',
     2         iota(1), iota(2), iota(3)
            END IF
            
            IF (icount .ne. 3) WRITE (iwrite, *)

 400        CONTINUE

            IF (icount.eq.6 .and. lfreeb) THEN
               IF (.not.ALLOCATED (extcur)) EXIT
               icoil = MAXVAL(LEN_TRIM(curlabel(1:nextcur_vmec)))
               icoil = MAX(icoil+4,8)                 !8 = 4+LEN("COIL")

               WRITE (line, '(a,i2.2,a,i2.2,a)') "(1x,'ID',",
     1                icoil-4,"x,'COIL',9x,'CURRENT',3x,'OPTIMIZE',
     2                /,", icoil+30,"('='))"
               WRITE (iwrite, line)
               WRITE (line, '(a,i2.2,a)') "(1x,i2.2,a",icoil,
     1                                    ",3x,1p,e13.6,1x,l6)"
               DO icoil = 1,nextcur_vmec
                  WRITE (iwrite, line) icoil, TRIM(curlabel(icoil)),
     1                                 extcur(icoil), lextcur(icoil)
               ENDDO

               ireadc = iread+1
               temp = "coil_targets." // TRIM(extension)
               CALL safe_open(ireadc, istat2, temp, 'old', 'formatted')
               IF (istat2 .eq. 0) THEN
                  ALLOCATE(coil_info(6,nsad_unique_coils))
                  coil_info = 0
                  jcount = 1
                  DO WHILE (jcount .le. 6)
                     READ (ireadc, '(a)', end=1020) line
                     jndex = INDEX(line, TRIM(parse_coils(jcount)))
                     IF (jndex .gt. 0) THEN
                        READ(ireadc, '(a)', end=1020) line
                        READ(line, *, iostat=istat2) 
     1                        coil_info(jcount, 1:nsad_unique_coils)
                        jcount = jcount+1
                     END IF
                  END DO
 1020             CONTINUE                  
                  CLOSE (ireadc)
                  WRITE (iwrite, '(/,a,f10.3,/)')
     1                  ' Plasma-coil minimum separation (m): ',
     2                  coil_info(6,1)
                  WRITE (iwrite, 1025)
 1025  FORMAT(1x,'SADDLE COIL',6x,'CURRENT (kA)',3x,'CURVATURE (m)',
     1        3x,'CC-SEP (m)',5x,'LENGTH (m)',5x,'MIN-Y (m)',/,95('='))
                  DO icoil = 1, nsad_unique_coils
                     WRITE (iwrite, 1030) icoil, 
     1                      coil_info(1,icoil)*1.E-3,
     1                      coil_info(3,icoil), coil_info(2,icoil),
     2                      coil_info(4,icoil), coil_info(5,icoil)
                  ENDDO
 1030 FORMAT(5x,i4,4x,5(3x,1pe12.3))
                  DEALLOCATE (coil_info)
               END IF
               
            ENDIF                     !IF (lfreeb) THEN

            REWIND (iscratch)
         END IF                    !IF istat.eq.0

      END DO FILE_LOOP

      IF (ALLOCATED(extcur)) DEALLOCATE(extcur, curlabel)
      CLOSE (iscratch)
      CLOSE (iwrite)

      END SUBROUTINE data_harvest
