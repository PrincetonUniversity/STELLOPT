      SUBROUTINE chisq_vac_island (ivar, num, nopt, iflag, extension,
     1   lscreen)

      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: home_dir, chisq_min, bigno, link
      USE optim_params, ONLY: sigma=>sigma_vac_island,
     1                        n_vac_island, m_vac_island
      USE system_mod, ONLY: system
      USE safe_open_mod
      USE vmec_input, ONLY: mgrid_file

      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ivar, nopt
      INTEGER :: num, iflag
      LOGICAL :: lscreen
      CHARACTER(LEN=*) :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(kind=rprec), PARAMETER :: zero = 0, one = 1,
     1   sjmax = one, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: unit_vac = 27
      INTEGER :: jcount, kcount, ierr, istat
      INTEGER :: j, k, iunit, m, n
      REAL(kind=rprec) :: ds, bnavg, ds_sum
      CHARACTER(len=LEN_TRIM(home_dir)+50) :: version
      CHARACTER(LEN=200) :: tmp
      CHARACTER(LEN=6), PARAMETER :: mgrid_prefix='mgrid_'
      LOGICAL :: ex
C-----------------------------------------------

      IF (.not. ANY(sigma(:) < bigno)) RETURN

      jcount = count(sigma(:) < bigno  .and.
     1               n_vac_island(:) /= 0
     1         .and. m_vac_island(:) > 0)

      version=TRIM(home_dir) // 'xvacisld'
      iunit=unit_vac

      IF (nopt > 0) THEN
!
!         If coils.extension does not exist, copy from parent
!         using extension from mgrid

          tmp='coils.' // TRIM(extension)
          INQUIRE(file=tmp, exist=ex, iostat=istat)
          IF(istat.ne.0 .or. .not.ex) THEN
             j = INDEX(mgrid_file,mgrid_prefix)
             IF (j .eq. 0 ) THEN
                iflag = -33
                RETURN
             ELSE
                j = LEN(mgrid_prefix) + 1
                k = LEN_TRIM(mgrid_file)
                IF (mgrid_file(k-2:k) == ".nc") k = k-3
                IF (mgrid_file(k-3:k) == ".bin") k = k-4
                tmp = link // "../coils." // TRIM(mgrid_file(j:k)) //
     1               " ./coils." // TRIM(extension)
                PRINT *, tmp
             ENDIF
             CALL system(tmp, ierr)
             IF (ierr .lt. 127 .and. ierr .ne. 0) THEN
                iflag = -33
                RETURN
             ENDIF
          ENDIF

          version = TRIM(version) // ' ' // extension
          CALL system(version, ierr)

          IF( ierr .lt. 127 .and. ierr .ne. 0 ) THEN
             PRINT *, "Trouble running XVACISLD, ierr = ", ierr
             iflag = -33     ! temporary flag setting
             RETURN
          ENDIF

          CALL safe_open(iunit, istat,'visldout.'//TRIM(extension),
     1                  'unknown', 'formatted')
          IF( istat .ne. 0 ) THEN
            PRINT *, "In Open VISLDOUT istat = ", istat
            iflag = -33
            RETURN
          ELSE

            IF( lscreen) THEN
                WRITE(6,*) ' Vac-Island Estimate'
                WRITE(6,*) ' n  m   ds     bnavg'
            ENDIF


            DO j=1, jcount

              ds_sum = zero
              READ(iunit,*) kcount
              READ(iunit,*)
              READ(iunit,*)
              DO k=1, kcount
                READ(iunit,*) m, n, ds, bnavg
                IF(m .eq. m_vac_island(j) .and.
     1             n .eq. n_vac_island(j) ) THEN

                   ds_sum=ds_sum+ds

                   IF( lscreen) WRITE(6,'(2i3, 2(1p,e11.2))')
     1                 n, m, ds, bnavg

                ENDIF
              ENDDO
              num = num + 1
              wegt(num) = sigma(j)
              index_array(num) = ivar
              chisq_match(num) = ds_sum
              chisq_target(num) = zero
              REWIND iunit
            ENDDO

            CLOSE(iunit)
          ENDIF

      ELSE

          INQUIRE(file=version, exist=ex, iostat=istat)
          IF(istat.ne.0 .or. .not.ex) THEN
            PRINT *, 'xvacisld not found in ' // TRIM(home_dir)
            sigma(:) = bigno
          ENDIF
       
          DO j=1, jcount
             num = num + 1
             IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
          ENDDO

          CALL safe_open(iunit, istat,'vacisld.param',
     1                  'unknown', 'formatted')
          WRITE(iunit,*) jcount
          DO j=1, jcount
            WRITE(iunit,*) m_vac_island(j), n_vac_island(j)
          ENDDO
          CLOSE(iunit)

      END IF

      END SUBROUTINE chisq_vac_island
