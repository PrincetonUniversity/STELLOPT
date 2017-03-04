      SUBROUTINE magnetics_data
      USE vmec_main
#ifdef _VACUUM2
      USE vac2_vacmod
#else
      USE vacmod
#endif
      USE vsvd
      USE mgrid_mod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iobset, j, i, n, ibfld
      REAL(rprec) :: denwgt, dsitotal,
     1   chierr, wght0, abs_flux, chiwgt, raddeg, denbwgt, btotal,
     2   relerr_b, wghtb, chisqx, chibwgt
C-----------------------------------------------
!
!       ALGEBRAIC SUM OF FLUXES FOR GROUPS OF OBSERVER POSITIONS
!
      IF (nobd .ne. 0) THEN
         flmwgt = 0
         denwgt = 0
         WRITE (nthreed, 400)
  400    FORMAT(/' External Poloidal Flux Loop Measurements (Wb)',/,
     1      '   I     Connection          Coils       Plasma    ',
     2      '    Total     Measured   Chi-Sq Err        Sigma',
     3      '   Designation',/,1x,3('-'),1x,16('-'),7(1x,12('-')),'-')

         DO iobset = 1, nobd
            dsitotal = dsiext(iobset) + plflux(iobset)
            chierr = dsiobt(iobset) - dsitotal
            wght0 = cbig + one
            IF (ANY(iobset.eq.indxflx(:nflxs)))
     1          wght0 = sigma_flux(iobset)
            flmwgt = flmwgt + (chierr/wght0)**2
            denwgt = denwgt + (dsitotal/wght0)**2
            chierr = (chierr/wght0)**2
            IF (wght0 .lt. cbig) THEN
               WRITE (nthreed, 410) iobset, (iconnect(j,iobset),j=1,4),
     1            dsiext(iobset), plflux(iobset), dsitotal, dsiobt(
     2            iobset), chierr, wght0, dsilabel(iobset)
            ELSE
               WRITE (nthreed, 420) iobset, (iconnect(j,iobset),j=1,4),
     1            dsiext(iobset), plflux(iobset), dsitotal, dsiobt(
     2            iobset), chierr, dsilabel(iobset)
            ENDIF
         END DO

         WRITE (nthreed, 600)
         DO i = 1, nobser
            abs_flux = plflux(nobd+i)
            WRITE (nthreed, 610) i, xobser(i), zobser(i), psiext(i),
     1         abs_flux, psiext(i) + ABS_flux
         END DO
  600    FORMAT(//,' Poloidal Flux at Individual Coil Locations',/,
     1      '   I         R [M]        Z [M]        Coils       Plasma',
     2      '        Total',/,2x,3('-'),5(1x,12('-')))
  610    FORMAT(i4,1x,5f13.4)
      ENDIF

      nchisaddle = nflxs
      total_saddle_chi = 0.
      IF (nflxs .GT. 0) THEN
         chiwgt = flmwgt/(nflxs)
         flmwgt = SQRT(flmwgt/denwgt)
         total_saddle_chi = (nflxs)*chiwgt
         WRITE (nthreed, 430) chiwgt, flmwgt
      ENDIF
  410 FORMAT(i4,1x,4i4,6f13.4,7x,a8)
  420 FORMAT(i4,1x,4i4,5f13.4,5x,'No Match',7x,a8)
  430 FORMAT(/' MEAN CHI-SQ ERROR IN SADDLE LOOP MATCH: ',1p,e10.3,/,
     1   ' RMS ERROR IN SADDLE LOOP MATCH: ',e10.3)

!
!       BR, BZ FIELD COMPONENTS FOR MATCHING AT OBSERVER POSITIONS
!
      total_b_chi = 0.
      IF (nbfldn .ne. 0) THEN
  500    FORMAT(//' BFIELD GROUP: ',a,' // ',
     1      'External B Loop Magnetic Field Measurements (T)',/,
     2      '  I  Rcoil[M]  Zcoil[M]    Degree    B[coil]',
     3      '  B[plasma]   B[Total]  ','  B[Meas]    Sigma  Chi-Sq Err',
     4      /,1x,2('-'),3(1x,9('-')),4(1x,10('-')),1(1x,8('-')),1(1x,11(
     5      '-')))

         raddeg = 360.0_dp/twopi
         DO n = 1, nbsets                        ! over all allowed sets
            WRITE (nthreed, 500) bloopnames(n)
            bcwgt = 0.
            denbwgt = 0.
            b_chi(n) = 0.
            ibfld = 0
            DO j = 1, nbcoils(n)
               ibfld = ibfld + 1
               btotal = bcoil(j,n) + plbfld(j,n)
               relerr_b = btotal - bbc(j,n)
               wghtb = cbig + 1.0_dp
c                                                ! over NAMELIST sets
               IF (ANY(j.eq.indxbfld(:nbfld(n),n))) wghtb = sigma_b(j,n)
               bcwgt = bcwgt + (relerr_b/wghtb)**2
               denbwgt = denbwgt + (btotal/wghtb)**2
               chisqx = (relerr_b/wghtb)**2
               IF (wghtb .LT. cbig) THEN
                  WRITE (nthreed, 520) ibfld, rbcoil(j,n), zbcoil(j,n),
     1               abcoil(j,n)*raddeg, bcoil(j,n), plbfld(j,n), btotal
     2               , bbc(j,n), wghtb, chisqx
               ELSE
                  WRITE (nthreed, 530) ibfld, rbcoil(j,n), zbcoil(j,n),
     1               abcoil(j,n)*raddeg, bcoil(j,n), plbfld(j,n), btotal
     2               , bbc(j,n)
               ENDIF
            END DO
            IF (nbfld(n) .gt. 0) THEN
               chibwgt = bcwgt/(nbfld(n))
               bcwgt = SQRT(bcwgt/denbwgt)
               b_chi(n) = (nbfld(n))*chibwgt
               WRITE (nthreed, 540) bloopnames(n), chibwgt,
     1            bloopnames(n), bcwgt
               total_b_chi = total_b_chi + b_chi(n)
            ENDIF
         END DO

      ENDIF
  520 FORMAT(i3,3f10.3,4f11.3,f9.4,f12.4)
  530 FORMAT(i3,3f10.3,4f11.3,3x,'No Match')
  540 FORMAT(/' MEAN CHI-SQ ERROR IN B-COIL GROUP ',a,': ',1p,e10.3,/,
     1   ' RMS ERROR IN B-COIL GROUP ',a,': ',e10.3)

      END SUBROUTINE magnetics_data
