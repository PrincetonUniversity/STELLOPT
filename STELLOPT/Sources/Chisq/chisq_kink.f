      SUBROUTINE chisq_kink (num, nopt, iflag,
     1    extension, lscreen)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: home_dir, chisq_min, bigno, link
      USE optim_params, ONLY: lkink_opt, sigma=>sigma_kink,
     1                        target_loc=>target_kink
      USE system_mod, ONLY: system
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nopt, iflag, num
      CHARACTER(LEN=*) :: extension
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: unit_kink = 27
      INTEGER :: iunit, ifail, istat, n
      REAL(rprec), PARAMETER :: zero = 0
      REAL(rprec) :: growth_rate
      CHARACTER(len=LEN_TRIM(home_dir)+30) :: version, ft5_ver
      LOGICAL :: ex
      CHARACTER(LEN=1) :: tag
C-----------------------------------------------
      version = TRIM(home_dir) // 'xtprp'
      iflag = 0

      IF (nopt > 0) THEN
         DO n=1, SIZE(sigma)
            IF( sigma(n) >= bigno) CYCLE
            IF( n == 1) THEN
               tag = ' '
            ELSE
               WRITE(tag,'(i1)') n
            ENDIF

            iunit = unit_kink
            iflag = 0

            IF (lscreen) WRITE(6,*)'Running TERPSICHORE('//tag//
     1                 ') for kink stability'
            CALL load_physics_codes(TRIM(version)//tag, " ", "-b " //
     1         TRIM(extension),'ft79tpr'//tag, extension, iunit, iflag)

            IF (iflag .ne. 0) THEN
               WRITE(6,*) 'Problems with load_physics(terp.), iflag=',
     1                     iflag
               RETURN
            ENDIF

            ifail = 0
            READ(iunit, *, iostat=istat) ifail, growth_rate

            IF (istat.ne.0 .or. ifail.eq.1) THEN
c             IF (lscreen)
                 WRITE(6,*) 'Error running Terpsichore'//tag//
     1               '/kink, status=',istat, ' Ifail = ', ifail
              IF( ifail .eq. 1) THEN

c        Terpsichore could not construct vacuum region.
c             IF (lscreen)
                 WRITE(6,*) 'Terpsichore'//tag//
     1                  ' error construction vacuum, ',
     1                  'Deprecate this direction!'

                 growth_rate = -100*SQRT(chisq_min*sigma(n))
                 iflag = 0
              ELSE
                 iflag = -19
                 RETURN
              ENDIF
           ENDIF

           num = num + 1
           index_array(num) = ivar_kink
           chisq_match(num) = 0

           chisq_target(num) = target_loc(n)
!          chisq_target(num) = sigma/2

           chisq_match(num) = MIN (chisq_target(num), growth_rate)
 
           CLOSE(iunit, status='delete')
           wegt(num) = sigma(n)
         ENDDO

      ELSE
        
         DO n=1, SIZE(sigma)
            IF( sigma(n) >= bigno) CYCLE
            IF( n == 1) THEN
               tag = ' '
            ELSE
               WRITE(tag,'(i1)') n
            ENDIF

            INQUIRE(file=TRIM(version)//tag, exist=ex, iostat=istat)
            IF (istat.ne.0 .or. .not.ex) THEN
               IF(lscreen) THEN
                  PRINT *, 'xtprp'//tag//' file not found in '//
     1                TRIM(home_dir)
                  PRINT *, '*** Kink evaluation disabled !'
               ENDIF
               lkink_opt = .false.
            ELSE
               INQUIRE(file='../ft5tpr'//tag, exist=ex, iostat=istat)
               IF (istat.ne.0 .or. .not.ex) THEN
                  IF(lscreen) THEN
                     PRINT *,' TERPSICHORE control file ft5tpr'//tag//
     1                    ' is missing'
                     PRINT *, '*** Kink evaluation disabled !'
                  ENDIF
                  lkink_opt = .false.
               ELSE
                  ifail = 0
                  ft5_ver = link // "../ft5tpr"//tag//" ./ft5tpr"
     1                      //tag
                  CALL system(TRIM(ft5_ver), ifail)
                  IF (ifail.lt.127 .and. ifail.ne.0) THEN
                     IF(lscreen) THEN
                        PRINT *,' TERPSICHORE control file ft5tpr'
     1                        //tag//' unlinked, status=', ifail
                        PRINT *, '*** Kink evaluation disabled !'
                     ENDIF
                     lkink_opt = .false.
                  ELSE
                     num = num + 1
                     IF (nopt .eq. -2) chisq_descript(num) = 
     1                                 descript(ivar_kink)
                  END IF
               END IF
            END IF
         ENDDO
      
      ENDIF

      END SUBROUTINE chisq_kink
