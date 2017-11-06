      SUBROUTINE newprofil(phipog)
      USE vmec_main
      USE vacmod
      USE realspace
      USE vforces, lu => czmn, lv => crmn
      USE vsvd
      USE vspline
      USE xstuff
      USE timer_sub
      USE mgrid_mod, ONLY: nbfldn
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(*) :: phipog
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: icount, inodes, js
      INTEGER, SAVE :: idata(jchix), isize(jchix)       !INCREMENTAL UPDATES
      REAL(rprec), DIMENSION(ipnodes + isnodes) :: datalsq, wten
      REAL(rprec), DIMENSION(1+nrzt) :: r12sqr
      REAL(rprec) :: amat_lsq(isnodes+ipnodes,isnodes+ipnodes),
     1   djac_p(ipnodes,nmeasurements), djac_i(isnodes,nmeasurements)
      REAL(rprec), DIMENSION(nmeasurements) :: datainput
      REAL(rprec) :: delt1, pfac0, aminout, aminin, ymin, pfactor

C-----------------------------------------------
!     IDATA: GIVES STARTING INDEX IN DATAINPUT ARRAY
!            FOR STORING EACH OF THE DATA TYPES
!     ISIZE: GIVES NUMBER OF DATA STARTING AT IDATA
!     INDEXING OF DJAC ARRAY (ASSUMES K STARTS AT 0)

      inodes = isnodes + ipnodes

      CALL second0 (treconon)

!     Unfreeze magnetic axis
      IF (iresidue.eq.0 .and. fsq*1.e6_dp.lt.one) iresidue = 1
      delt1 = one/REAL(ipedsvd,rprec)
      IF (iresidue .eq. 0) delt1 = one

!
!       COMPUTE AVERAGE RADIAL FORCE BALANCE CONSTRAINT
!
      CALL radfor (pfac0)
      pfac = pfac + delt1*(pfac0 - pfac)
!
!       UPDATE PHI SCALE FACTOR (PHIFAC)
!       SCALE TOROIDAL FLUX TO MATCH PRESSURE WIDTH OR LIMITER
!
      aminout = MAX(rthommax,rstarkmax) - r00
      aminin = r00 - MIN(rthommin,rstarkmin)
      apres = (aminin*ipreSIN + aminout*ipresout)/(ipreSIN +ipresout)
      aminor = ((r00 - rinner)*ipreSIN + (router - r00)*ipresout)/
     1       (ipreSIN + ipresout)
      IF (imatch_phiedge.ne.1 .and. ivac.gt.1 .or. imatch_phiedge.eq.3
     1   .and. (.not.lfreeb)) THEN
         CALL newphi (phipog)
         CALL gettflux
      ENDIF

      icount = 0

      IF (.not.(MOD(iter2 - iter1,ipedsvd).ne.0 .and. iequi.eq.0
     1     .and. iresidue.gt.0)) THEN

!
!       SETUP COMMON BLOCKS FOR FLUX-MATCHING ROUTINES
!

         IF (iphidiam + nflxs + nbfldn.gt.0 .or. iequi.gt.0) THEN
            r12sqr(2:nrzt) = SQRT(armn_o(2:nrzt))
            CALL flux_init (phipog)
         ENDIF

!
!       COMPUTE MATRIX ELEMENTS FOR THOMPSON SCATTERING DATA
!
         idata(ithom0) = icount + 1
         CALL getthom(djac_i(1,idata(ITHOM0)), djac_p(1,idata(ITHOM0)),
     1     datainput(idata(ITHOM0)), r1(1:,0), r1(1:,1), isize(ITHOM0))
         icount = icount + isize(ithom0)

!
!       COMPUTE MOTIONAL STARK EFFECT. THIS CALL ALSO INITIALIZES
!       THE ALSQ, DATALSQ ARRAYS AND SETS UP THE SPLINE NODES.
!
         idata(istark0) = icount + 1
         CALL getmse(djac_i(1,idata(ISTARK0)), djac_p(1,idata(ISTARK0)),
     1     datainput(idata(ISTARK0)), r1(1:,0), r1(1:,1),lu,
     2     lu(1+nrzt), zu(1:,0), zu(1:,1), phipog, isize(ISTARK0))
         icount = icount + isize(istark0)

!
!       COMPUTE MATRIX ELEMENTS FOR DIAMAGNETIC FLUX LOOP
!
         idata(idiam0) = icount + 1
         CALL getdiam (djac_i(1,idata(idiam0)), djac_p(1,idata(idiam0))
     1      , datainput(idata(idiam0)), isize(idiam0))
         icount = icount + isize(idiam0)

!
!       COMPUTE MATRIX ELEMENTS FOR EXTERNAL POLOIDAL FLUXES
!
         idata(iflxs0) = icount + 1
         CALL getflux (djac_i(1,idata(iflxs0)), djac_p(1,idata(iflxs0))
     1      , datainput(idata(iflxs0)), r12sqr, clmn_e(1), clmn_o(1),
     2      blmn_o(1), armn_o(1), blmn_e(1), azmn_o(1), isize(iflxs0))
         icount = icount + isize(iflxs0)

!
!       COMPUTE MATRIX ELEMENTS FOR EXTERNAL MAGNETIC FIELD MATCHING
!
         idata(ibrzfld) = icount + 1
         CALL getbfld (djac_i(1,idata(ibrzfld)), djac_p(1,idata(ibrzfld)
     1      ), datainput(idata(ibrzfld)), r12sqr, azmn_o(1), blmn_o(1),
     2      clmn_e(1), clmn_o(1), armn_o(1), blmn_e(1), isize(ibrzfld))
         icount = icount + isize(ibrzfld)

!
!       SQUARE DATA MATRIX ELEMENTS AND STORE IN ALSQ
!

         IF (icount .gt. nmeasurements) STOP 'icount>nmeasurements'
         IF (iequi .eq. 0) THEN

            CALL sgemvmm (djac_i, djac_p, amat_lsq, datainput, datalsq,
     1         wten, icount, isnodes, ipnodes, inodes)

!
!       COMPUTE IOTA, PRESSURE SPLINE COEFFICIENTS
!
            CALL set_dual (datalsq, hstark, ystark, y2stark, hthom,
     1         ythom, y2thom, wten, amat_lsq, isnodes, ipnodes, inodes)

            IF (.not.lpprof) THEN
               ymin = MINVAL(ythom(1:ipnodes))
               ythom(:ipnodes) = ythom(:ipnodes) - ymin
            ENDIF

!
!       COMPUTE IOTA, PRESSURE AT R(js) FROM SPLINED INPUT
!       DATA ALONG THE MIDPLANE
!
            DO js = 1, ns
               CALL splint (sknots, ystark, y2stark, isnodes, 
     1              sqrts(js), isplinef(js))
               CALL splint (pknots, ythom, y2thom, ipnodes, 
     1              sqrts(js), psplinef(js))
               IF (js .eq. 1) CYCLE
               CALL splint (sknots, ystark, y2stark, isnodes, 
     1              shalf(js),isplineh(js))
               CALL splint (pknots, ythom, y2thom, ipnodes, 
     1              shalf(js), psplineh(js))
            END DO

            pfactor = mu0*pthommax             !!*pfac moved to getthom
            DO js = 1,ns
              isplinef(js) = isplinef(js) - iotaf(js)
              isplineh(js) = isplineh(js) - iotas(js)
              psplinef(js) = pfactor*psplinef(js) - presf(js)
              psplineh(js) = pfactor*psplineh(js) - mass(js)
            END DO

         ENDIF                     ! iequi>0
!
!       COMPUTE CHISQ
!
         CALL chisq (djac_i, djac_p, datainput, idata, isize, icount)

      ENDIF                        ! MOD(iter2-iter1,ipedsvd) == 0

      IF (iequi .eq. 0) THEN
!
!       UPDATE PRESSURE SPLINE AND ESTABLISH INTERNAL
!       (CODE) UNITS OF PRESSURE. P(real) = mu0 * pthommax * P(splined)
!       WHERE P(code-units) = mu0 * P(real)
!       SMOOTH TIME VARIATION OF PROFILES
!
         DO js = 1,ns
           iotaf(js) = iotaf(js) + delt1*isplinef(js)
           iotas(js) = iotas(js) + delt1*isplineh(js)
           presf(js) = presf(js) + delt1*psplinef(js)
           mass(js)  = mass(js)  + delt1*psplineh(js)
         END DO
      ENDIF

!
!     STORE CHISQ
!
      CALL store_chisq

!
!       OPTIMIZE MAGNETIC AXIS POSITION BY MINIMIZING RMS ERROR
!       RETURNS RAXMSE AS UPDATED GUESS FOR NEW AXIS POSITION
!       TO BE USED IN SUBROUTINE RESIDUE
!
      CALL axisopt (fsq, r00, iresidue, ivac)

!       Compute force to fix axis at RAXMSE
      grmse = -0.05*(r00 - raxmse)

      CALL second0 (treconoff)
      timer(trecon) = timer(trecon) + (treconoff - treconon)

      END SUBROUTINE newprofil
