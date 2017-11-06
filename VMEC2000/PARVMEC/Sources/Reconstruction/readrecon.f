      SUBROUTINE readrecon
      USE vmec_main
      USE vsvd
      USE mgrid_mod, ONLY: nbfldn, nobser, nobd, nbcoilsn
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      CHARACTER(LEN=50), DIMENSION(2), SAVE :: raxis_message
      CHARACTER(LEN=50), DIMENSION(4), SAVE :: phiedge_message
C-----------------------------------------------
      data raxis_message/'Magnetic axis position fixed',
     1   'Magnetic axis position optimized'/
      data phiedge_message/
     1   'Phiedge determined to match pressure minor radius',
     2   'Phiedge matched to input value',
     3   'Phiedge determined by limiter position',
     4   'Phiedge determined by Ip'/

      iphidiam = 0
      IF (imse > nmse) STOP 'IMSE>NMSE'
      IF (itse > ntse) STOP 'ITSE>NTSE'
      lrecon = ((imse>0 .or. nflxs>0 .or. nbfldn>0) .and. itse.ge.0)
      IF (lrecon) iresidue = 0

      IF (.not.lrecon) RETURN                   !No reconstruction matching

      ncurr = 0                                  !Just to be safe
      IF (sigma_current .ge. cbig) STOP 'SIGMA_CURRENT missing'
      IF (sigma_delphid .ge. cbig) PRINT *, ' SIGMA_DELPHID missing'
      IF (sigma_current < zero) sigma_current = ABS(sigma_current*
     1   curtor)
      IF (sigma_delphid < zero) sigma_delphid = ABS(sigma_delphid*
     1   phidiam)
      WRITE (nthreed, 150)
  150 FORMAT(/' DATA MATCHING PARAMETERS: ',/,1x,35('-'))
      WRITE (nthreed, 155) imse, itse, nflxs, nobser, nobd, nbfldn,
     1   nbcoilsn, sigma_current, 1.e3_dp*sigma_delphid, tensp, tensi,
     2   tensi2, fpolyi, mseangle_offset, presfac, pres_offset, lpofr
      WRITE (nthreed, 152) mseangle_offsetm
  152 FORMAT('mse-angleM offset',/,f13.3)
  155 FORMAT('   imse       itse      nflxs     nobser       nobd',
     1'     nbfldn    nbcoils  sigma_current(A)   sigma_delphid(mWb)',/,
     2   i7,6i11,3x,1p,e15.3,4x,0p,f17.3,/,
     3   '    tension(p)   tension(i)  tension2(i)  fpolyi  ',
     4   'mse-angle offset  pres scale factor pressure offset  lpofr',/,
     5   3f13.3,f9.3,f18.3,f19.3,f16.3,6x,l1)
      WRITE (nthreed, 200)
  200 FORMAT(/,' LEGEND',/,1x,6('-'))

      IF (curtor < cbig) THEN
         WRITE (nthreed, 210) 1.e-6_dp*curtor
      ELSE
         WRITE (nthreed, *) 'Need toroidal plasma current'
         STOP 15
      ENDIF
  210 FORMAT(' Matching to toroidal current = ',f10.3,' [MA]')
      sigma_current = mu0*sigma_current
      IF (nflxs > 0) THEN
         WRITE (nthreed, *) 'Fitting ', nflxs,
     1      ' EXTERNAL flux loop measurements'
      ELSE
         WRITE (nthreed, *)
     1      'Not fitting EXTERNAL flux loop measurements.'
      ENDIF
      IF (phidiam<cbig .and. sigma_delphid<cbig) THEN
         iphidiam = 1
         WRITE (nthreed, 220) 1.e3_dp*phidiam
      ELSE
         WRITE (nthreed, *) 'No fit to diamagnetic flux'
      ENDIF
  220 FORMAT(' Fitting diamagnetic flux     = ',f10.3,' [mWb]')
      WRITE (nthreed, *) raxis_message(iopt_raxis+1)
      WRITE (nthreed, *) phiedge_message(imatch_phiedge+1)

      END SUBROUTINE readrecon
