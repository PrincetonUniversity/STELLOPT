      SUBROUTINE chisq_dsubr (sigma, ivar, num,
     1                      nopt, iflag, extension, lscreen)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: home_dir, nrad
      USE optim_params, ONLY: ldsubr_opt
      USE system_mod, ONLY: system
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ivar, nopt, iflag, num
      REAL(rprec), DIMENSION(*) :: sigma
      CHARACTER(LEN=*) :: extension
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: unit_jmc = 79
      INTEGER :: iunit, ierr, istat, i, ns_jmc_max
      INTEGER :: nsval
      REAL(rprec), PARAMETER :: zero = 0
      REAL(rprec), DIMENSION(nrad)  :: dsubr
      INTEGER, DIMENSION(nrad) :: ns_jmc
      CHARACTER(LEN=200) :: version
      LOGICAL :: ex
C-----------------------------------------------
      version = TRIM(home_dir) // 'xjmc'

      IF (nopt > 0) THEN

         iunit = unit_jmc

         IF (lscreen) THEN
           WRITE(6,*)
     1           'Running JMC for D_sub_r calculations'
         ENDIF

         ierr = 0
         version = TRIM(version) // ' ' // extension
         CALL system(version,ierr)
         IF( ierr .lt. 127 .and. ierr .ne. 0 ) THEN
             PRINT *, "Trouble running JMC, ierr = ", ierr
             RETURN
         ENDIF

         CALL safe_open(iunit, istat, 'ft79jmc.'//TRIM(extension),
     1                  'old', 'formatted')
         IF( istat .ne. 0 ) THEN
            iflag = -32
            RETURN
         ELSE
            READ(unit_jmc,*) ns_jmc_max
            READ(unit_jmc,*) (ns_jmc(i),dsubr(i),i=1,ns_jmc_max)
            IF( ns_jmc_max .ne. nrad - 3 ) STOP "Error in JMC"
            CLOSE(unit_jmc)
         ENDIF
!
!         calculated on all surfaces from 2 to nrad-2
!
          DO i=1, ns_jmc_max
             num = num + 1
             index_array(num) = ivar
             nsval = ns_jmc(i)
             chisq_match(num) =  -dsubr(i)
             chisq_target(num) =  MIN(chisq_match(num), zero)
             wegt(num) = sigma(nsval)
          ENDDO

      ELSE
         INQUIRE(file=version, exist=ex, iostat=istat)
         IF (istat.ne.0 .or. .not.ex) THEN
            PRINT *, 'xjmc file not found in ' // TRIM(home_dir)
            ldsubr_opt = .false.
         ELSE
            DO i=1, ns_jmc_max
               num = num + 1
               IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
            END DO
         END IF
      ENDIF

      END SUBROUTINE chisq_dsubr
