      SUBROUTINE write_rbzb(iunit, istat)
      USE optim
      USE vmec_input
      USE write_array_generic
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: iunit, istat
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      CHARACTER(LEN=*), DIMENSION(7), PARAMETER :: outfmt =
     1   (/ "(2(a6,i1,a1,i1,a4,1p,e22.14,3x))",
     2      "(2(a6,i1,a1,i2,a4,1p,e21.14,3x))",
     3      "(2(a6,i2,a1,i1,a4,1p,e21.14,3x))",
     4      "(2(a6,i2,a1,i2,a4,1p,e21.14,3x))",
     5      "(2(a6,i3,a1,i1,a4,1p,e21.14,3x))",
     6      "(2(a6,i3,a1,i2,a4,1p,e21.14,3x))",
     7      "(2(a6,i ,a1,i ,a4,1p,e21.14,3x))" /)
      INTEGER, PARAMETER :: n_max_leg = 6                               !! LEGENDRE (uses 6+1coeff.)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n, m, i
      INTEGER :: js                                                     !! LEGENDRE
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: ac_leg, ai_leg          !! LEGENDRE
      REAL(rprec), DIMENSION(nrad) :: x_rad                             !! LEGENDRE
      CHARACTER(len=LEN(outfmt(1))) :: outcfmt
C-----------------------------------------------
!
!     CONVERT MESH QUANTITIES (iota_opt, jcurv_opt) TO S-EXPANSION COEFFICIENTS
!

      IF (lcur_prof_opt .and. nrad.gt.3) THEN
        x_rad(2:nrad) = (/( (REAL(js, rprec)-0.5_dp)
     1    /REAL(nrad-1, rprec), js=1,nrad-1)/)
        ALLOCATE(ai_leg(0:n_max_leg), stat=istat)
        IF (istat .ne. 0) GOTO 2000
        CALL point_wise_to_power(nrad-1, x_rad(2), iota_opt(2),
     1    n_max_leg, ai_leg, 1, 1000)

        ai = 0
        ai(0:n_max_leg) = ai_leg(0:n_max_leg)
        DEALLOCATE(ai_leg)

      ELSE IF (liota_prof_opt .and. nrad.gt.3) THEN                     !(JCURV_OPT on full-mesh, NOT integrated in s)

        x_rad(1:nrad) = (/( REAL(js-1, rprec)/(nrad-1), js=1,nrad)/)
        ALLOCATE(ac_leg(0:n_max_leg), stat=istat)
        IF (istat .ne. 0) GOTO 2000
        CALL point_wise_to_power(nrad, x_rad, jcurv_opt,
     1    n_max_leg, ac_leg, 1, 1000)

        ac = 0
        ac(0:n_max_leg) = ac_leg(0:n_max_leg)
        DEALLOCATE(ac_leg)
      ENDIF

      WRITE (iunit, '(2x,3a)') "PMASS_TYPE = '",TRIM(pmass_type),"'"
      WRITE (iunit, '(2x,a,1p,e14.6)') 'GAMMA = ', gamma
      WRITE (iunit, '(2x,a,1p,e14.6)') 'BLOAT = ', bloat
      WRITE (iunit, 100) '  SPRES_PED = ',spres_ped
      WRITE (iunit, 100) '  PRES_SCALE = ',pres_scale
      WRITE (iunit, 100, err=2000) '  AM = ', (am(n-1), n=1,SIZE(am))
      i = minloc(am_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit, 100, err=2000) '  AM_AUX_S = ',
     1         (am_aux_s(n), n=1,i)
         WRITE (iunit, 100, err=2000) '  AM_AUX_F = ', 
     1         (am_aux_f(n), n=1,i)
      END IF
      WRITE (iunit, '(2x,3a)') "PIOTA_TYPE = '",TRIM(piota_type),"'"
      WRITE (iunit, '(2x,3a)') "PCURR_TYPE = '",TRIM(pcurr_type),"'"
      WRITE (iunit, '(2x,a,i4)') 'NCURR = ', ncurr
      WRITE (iunit, 100) '  CURTOR = ',curtor
      WRITE (iunit, 100, err=2000) '  AI = ', (ai(n-1), n=1,SIZE(ai))
      i	= minloc(ai_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit, 100, err=2000) '  AI_AUX_S = ', 
     1         (ai_aux_s(n), n=1,i)
         WRITE (iunit, 100, err=2000) '  AI_AUX_F = ', 
     1         (ai_aux_f(n), n=1,i)
      END IF
      WRITE (iunit, 100, err=2000) '  AC = ', (ac(n-1), n=1,SIZE(ac))
      i	= minloc(ac_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit, 100, err=2000) '  AC_AUX_S = ', 
     1         (ac_aux_s(n), n=1,i)
         WRITE (iunit, 100, err=2000) '  AC_AUX_F = ',
     1         (ac_aux_f(n), n=1,i)
      END IF
      IF (lanimec) THEN
         WRITE (iunit, 100) '  BCRIT = ',bcrit
         WRITE (iunit, 100, err=2000) '  AH = ', (ah(n-1), n=1,SIZE(ah))
         WRITE (iunit, 100, err=2000) '  AT = ', (at(n-1), n=1,SIZE(at))
      END IF
         
!     IF (ANY(aphi .ne. zero))
!    1  WRITE (iunit, 100, err=2000) '  APHI = ', (aphi(n), n=0,10)
!
!     WRITE out both the new axis and the original one, so the user can easily
!     try either in case of convergence problems.  The original axis is written
!     second, so it will have precedence. NOTE: RBC, ZBS WILL BE REASSIGNED TO
!     THE BOUNDARY VALUES IN read_wout_opt IF LFREEB = T
!
!     BE SURE TO WRITE OUT ORIGINAL RAXIS, ZAXIS (THEY MAY BE OVERRIDDEN IN read_wout_opt, clean_up)
!     THIS WILL PRESERVE CHI-SQ IN .MIN FILE
!
      CALL write_array(iunit, 'RAXIS_CC', raxis_cc, ntor+1, low_index=0)
      CALL write_array(iunit, 'ZAXIS_CS', zaxis_cs, ntor+1, low_index=0)
      CALL write_array(iunit, 'RAXIS_CS', raxis_cs, ntor+1, low_index=0)
      CALL write_array(iunit, 'ZAXIS_CC', zaxis_cc, ntor+1, low_index=0)
      CALL write_array(iunit, 'RAXIS_CC', raxis_old(0:ntor,1), ntor+1, 
     1                 low_index=0)
      CALL write_array(iunit, 'ZAXIS_CS', zaxis_old(0:ntor,1), ntor+1, 
     1                 low_index=0)
      CALL write_array(iunit, 'RAXIS_CS', raxis_old(0:ntor,2), ntor+1, 
     1                 low_index=0)
      CALL write_array(iunit, 'ZAXIS_CC', zaxis_old(0:ntor,2), ntor+1, 
     1                 low_index=0)
      DO m = 0, mpol-1
         DO n = -ntor, ntor
            IF ((rbc(n,m).ne.zero) .or. (zbs(n,m).ne.zero)) THEN
!
!     handle formatting for up to 2 digit n by 3 digit m.
!     while this is probably overkill, we at least have to handle up
!     thru 2 digit n by 2 digit m to handle known cases.
!
               outcfmt = outfmt(1)
               IF( m > 9) outcfmt = outfmt(2)
               IF(( n>-10 .and. n<0 ) .or. (n > 9 .and. n < 100)) THEN
                  outcfmt = outfmt(3)
                  IF( m > 9) outcfmt = outfmt(4)
               ELSE IF( n>-100 .and. n< -9 ) THEN
                  outcfmt = outfmt(5)
                  IF( m > 9) outcfmt = outfmt(6)
               ENDIF

               IF( n>= 100 .or. n <= -100 .or. m >=100) THEN
                  outcfmt = outfmt(7)
               ENDIF

               WRITE (iunit, outcfmt, err=2000)
     1            '  RBC(', n, ',', m, ') = ', rbc(n,m),
     2            '  ZBS(', n, ',', m, ') = ', zbs(n,m)
               IF (lasym) WRITE (iunit, outcfmt, err=2000)
     1            '  RBS(', n, ',', m, ') = ', rbs(n,m),
     2            '  ZBC(', n, ',', m, ') = ', zbc(n,m)
            ENDIF
         END DO
      END DO

 100  FORMAT(a,(1p,4e22.14))

      istat = 0
      RETURN

 2000 istat = -5

      END SUBROUTINE write_rbzb
