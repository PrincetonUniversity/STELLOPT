      SUBROUTINE fixrecon(ier)
      USE vmec_main
      USE vsvd
      USE vspline
      USE mgrid_mod, ONLY: nobser, nobd, nbsets, needflx, iconnect,
     1                     xobser, zobser, dsiext, psiext, plflux,
     2                     nbcoils, needbfld, rbcoil, zbcoil,
     3                     abcoil, plbfld, bcoil
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER ier
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: angle_variance  = 0.3_dp
      REAL(rprec), PARAMETER :: radial_variance = 0.3_dp
      REAL(rprec), PARAMETER :: p_threshold = 1.e-3_dp
      REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
      INTEGER, PARAMETER :: inode_max = 15
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(itse) :: isortp
      INTEGER, DIMENSION(1) :: isamax
      INTEGER :: istat1 = 0, ii, i, ineg,
     1   ipos, ipmax, ioff, ileft, n, n1, icount, index1, ind, m, m1
      REAL(rprec) :: delstark, datapos, dataneg
      REAL(rprec), DIMENSION(imse+1) :: datalsq_s, stark_temp
      REAL(rprec) :: datamax, datamin, t1, rneg,
     1   rpos, datanorm, presmax, presin, presout, presmin, tsign
      REAL(rprec), DIMENSION(itse) ::
     1   datalsq_p, ythom0, y2thom0, qtemp
      LOGICAL :: l1v(imse), l4v(itse)
C-----------------------------------------------

      CALL free_mem_recon

      ier = 0
!
!       ONLY QUANTITIES INDEPENDENT OF RADIAL MESH CAN GO HERE
!
!       STARK-DATA CONSISTENCY CHECK
!
      delstark = TAN(dcon*angle_variance)

!
!       SORT STARK DATA IN ASCENDING ORDER IN R-SPACE
!       AND RE-index1 DATASTARK, SIGMA_STARK ARRAYS ON ASCENDING RSTARK ARRAY
!
!       SCALE MOTIONAL STARK DATA TO RADIANS
!       RSTARK = R position along Z=0 plane of measurement
!       DATASTARK = ARCTAN ( Bpol/Btor ) at RSTARK, in degrees
!       AND IS CONVERTED TO Bpol/Btor
!

      ii = 0
      l1v(:imse) = (rstark(:imse) .gt. zero)
      DO i = 1, imse
         IF (l1v(i)) ii = ii + 1
      END DO
      IF (ii .ne. imse) THEN
         PRINT *, 'There is a zero in the RSTARK array ?!'
         ier = 1
         RETURN
      ENDIF

      ALLOCATE (indexs1(imse+3), indexu1(imse+3), isortr(imse+3),
     1   isorts(imse+3), delse1(imse+3), delso1(imse+3),
     2   starkcal(imse+3), qmeas(imse+3), qcalc(imse+3),
     3   fpsical(imse+3), stark_weight(imse+3),
     4   rsort(imse+3), rsort0(imse+3), stat=istat1)
      IF (istat1 .ne. 0) THEN
         PRINT *, ' ISTAT1 = ', istat1, ' ALLOCATING index1S1'
         STOP
      ENDIF

      datalsq_s(:imse) = datastark(:imse)
      stark_temp(:imse) = sigma_stark(:imse)
      CALL sort_data (rstark, isorts, imse)
      DO i = 1, imse
         datastark(i) = datalsq_s(isorts(i))
         sigma_stark(i) = stark_temp(isorts(i))
         IF (sigma_stark(i) .ge. cbig) THEN
            PRINT *, 'SIGMA_STARK missing'
            ier = 1
            RETURN
         ENDIF
         IF (sigma_stark(i) .lt. zero) sigma_stark(i) =
     1     ABS(sigma_stark(i)*datastark(i))
          !CONVERT TO Bpol/Btor ... applying profile offsets along the way!
         datastark(i) = TAN(dcon*(datastark(i)+mseangle_offset+
     1      mseangle_offsetm*mseprof(i)))
         sigma_stark(i) = TAN(dcon*sigma_stark(i))
         IF (sigma_stark(i) .eq. zero) sigma_stark(i) = 0.1*delstark
      END DO
      rstarkMIN = rstark(1)
      rstarkMAX = rstark(imse)

!
!     NEED FOR SCALING EDGE IOTA
!     HERE, SINCE RSTARK IS ORDERED, DATAMAX -> RIGHT OF AXIS
!     AND DATAMIN -> LEFT OF AXIS
!
!       GET RC0MSE (ESTIMATE FOR MAGNETIC AXIS)
!
      rwidth = radial_variance*(rbc(0,1)+ABS(rbs(0,1)))

      IF (imse .gt. 0) THEN
         datamax = datastark(imse)
         datamin = datastark(1)
         rc0mse = 0.
         ineg = 1
c                                     !NO ZERO CROSSING: FIND MIN ANYHOW
          IF (datamax * datamin.gt.zero) THEN      !NO ZERO CROSSING: FIND MIN ANYHOW
            datamin = ABS(datamin)
            DO i = 2,imse
              t1 = ABS(datastark(i))
              IF (t1.lt.datamin) THEN
                datamin = t1
                ineg = i
              ENDIF
            END DO
             IF (ineg.eq.imse) ineg = ineg-1
            ipos = ineg + 1
            GOTO 310
          ELSE IF (( datamax*signiota.lt.zero ) .or.
     >             (datamin*signiota.gt.zero)) THEN
            datastark = -datastark
          ENDIF

!
!       ALLOW FOR POSSIBLE MULTIPLE ZERO CROSSINGS (WIGGLES) IN DATA
!

         DO i = 1, imse
            IF (datastark(i)*signiota .le. zero) THEN
               ineg = i                          !LEFT OF MAGNETIC AXIS
            ELSE
               EXIT                              !RIGHT OF MAGNETIC AXIS
            ENDIF
         END DO
         DO i = imse, ineg + 1, -1
            IF (datastark(i)*signiota .le. zero) THEN
               EXIT                              !LEFT OF MAGNETIC AXIS
            ELSE
               ipos = i                          !RIGHT OF MAGNETIC AXIS
            ENDIF
         END DO

  310    CONTINUE

         rneg = rstark(ineg)
         rpos = rstark(ipos)
         dataneg = datastark(ineg)
         datapos = datastark(ipos)
      ENDIF                                      !End of IF(imse>0)
      IF (datapos .ne. dataneg) THEN
         rc0mse = (datapos*rneg - dataneg*rpos)/(datapos - dataneg)
         rwidth=delstark*ABS((rneg-rpos)/(datapos-dataneg))+rwidth
      ENDIF
      IF (ipos .gt. ineg + 1) rc0mse = 0.5_dp*(rpos + rneg)

!
!       ESTIMATE MAGNETIC AXIS FROM RAXIS
!
      raxmse = rc0mse
      IF (rc0mse.eq.0.0_dp .or. iopt_raxis.ne.1) raxmse = raxis_cc(0)
      rstepx0 = 0.005_dp*rwidth

!
!       COMPUTE SPLINES IN R-SPACE FOR MATCHING IOTA(0)
!

      datanorm = zero
      delse1(:imse) = one
      qcalc(:imse) = one/sigma_stark(:imse)**2
      datalsq_s(:imse) = datastark(:imse)*qcalc(:imse)
      scstark = SUM(ABS(datastark(:imse)/sigma_stark(:imse)))
      datanorm = SUM(ABS(datastark(:imse)))
      scstark = datanorm/scstark
c04-96        CALL setspline(rstark,qcalc,datalsq_s,stark_temp,ystark0,
c04-96     >  y2stark0,delse1,0.1*tensi/scstark**2,imse,NATUR)

!
!       DETERMINE NUMBER OF IOTA KNOTS IN SQRT(S)-SPACE
!
      IF (isnodes .le. 0) THEN
         isnodes = MIN(inode_max,imse + 1)
         isnodes = MAX(5,isnodes)               !At least 5 spline knots
      ENDIF
      IF (isnodes .lt. 5) STOP 'MUST PICK ISNODES > 4'
      WRITE (nthreed, *)
     1   'Number of iota-spline knots (in s-space):     ', isnodes

      ALLOCATE (nk_ia(isnodes), nk_ib(isnodes), hstark(isnodes),
     1  y2stark(isnodes), ystark(isnodes), sknots(isnodes), stat=istat1)
      IF (istat1 .ne. 0) THEN
         PRINT *, ' ISTAT1 = ', istat1, ' ALLOCATING NK_IA'
         STOP
      ENDIF

!
!       COMPUTES NODES IN SQRT(S) SPACE FOR SPLINES
!       THIS ASSUMES A FIXED NO - ISNODES - OF KNOTS
!       THIS MAY NOT PRESERVE ISNODES.
!       ALSO, IT IS NOT NECESSARY TO TAKE EQUALLY SPACED KNOTS IN
!       SQRT(S) SPACE. INDEED, THE FOLLOWING CHOICES ARE POSSIBLE:
!
!       SKNOTS(I) = HNODES*(I-1)   .eq.>   EQUAL-SPACED IN SQRT(S)
!
!       SKNOTS(I) = SQRT(HNODES*(I-1)) .eq.>  EQUAL-SPACED IN S
!
!       DO NOT - UNDER ANY CIRCUMSTANCES - CHANGE THE ARGUMENTS TO
!       THE SPLINT, GETSPLINE, SETUP_int ROUTINES FROM SQRTS,SHALF
!       TO SQRTS**2, SHALF**2 TO DO S-INTERPOLATION. RATHER, CHANGE
!       SKNOTS (AND PKNOTS) ACCORDING TO THE ABOVE FORMULA. THIS IS
!       ABSOLUTELY CRUCIAL, SINCE ONLY IN SQRT(S) SPACE DO THE
!       FIRST DERIVATIVE BOUNDARY CONDITIONS, d IOTA/d SQRT(S) = 0
!       (SIMILAR FOR P) APPLY AT THE AXIS, S=0.
!
!
      DO i = 1, isnodes
         sknots(i) = REAL(i - 1,rprec)/(isnodes - 1)
      END DO

      hstark(:isnodes-1) = sknots(2:isnodes) - sknots(:isnodes-1)


!
!       SET UP DATA ARRAY SCALE FACTORS
!       ACTUAL PRESSURE = PRESPEAK * PFAC * P(INTERNAL)
!       IF (LPOFR) THEN DATA ARE INPUT vs R (REAL SPACE)
!       IF (.NOT.LPOFR),DATA ARE INPUT vs S (FLUX SPACE)
!
      IF (itse .eq. 0) CALL getpresprofile         !!Simulate 'data'
      IF (itse .gt. 0) THEN

         ALLOCATE (sthom(itse), delse2(itse), delso2(itse), pcalc(itse),
     1      indexs2(itse), indexu2(itse), stat=istat1)
         IF (istat1 .ne. 0) THEN
            PRINT *, ' ISTAT1 = ', istat1, ' ALLOCATING STHOM'
            STOP
         ENDIF


         datathom(:itse) = datathom(:itse)*presfac
         presmax = MAXVAL(datathom(:itse))

!
!       SORT DATA IN ASCENDING ORDER IN R-SPACE (LPOFR) OR S-SPACE(.NOT.LPOFR)
!       AND RE-index1 DATATHOM, SIGMA_THOM ARRAYS ON ASCENDING RTHOM ARRAY

         datalsq_p(:itse) = datathom(:itse)
         qtemp(:itse) = sigma_thom(:itse)
         CALL sort_data (rthom, isortp, itse)
         IF (lpofr) THEN
            DO i = 1, itse
               datathom(i) = datalsq_p(isortp(i))
               rthom(i) = rthom(i) + pres_offset
               IF (rthom(i) .le. zero) THEN
                  PRINT *, 'Check units of PRES_OFFSET: rthom < 0!'
                  ier = 1
                  RETURN
               ENDIF
               IF (datathom(i) .eq. presmax) ipMAX = i
               sigma_thom(i) = qtemp(isortp(i))
               IF (sigma_thom(i) .ge. cbig) THEN
                  PRINT *, 'SIGMA_THOM missing'
                  ier = 1
                  RETURN
               ENDIF
               IF (sigma_thom(i) .lt. zero) THEN
                  sigma_thom(i) = ABS(sigma_thom(i)*datathom(i))
               ELSE
                  IF (sigma_thom(i) .gt. zero) THEN
                     sigma_thom(i) = presfac*sigma_thom(i)
                  ELSE
                     sigma_thom(i) = p_threshold*presmax
                  ENDIF
               ENDIF
            END DO
         ELSE
            DO i = 1, itse
               datathom(i) = datalsq_p(isortp(i))
               sthom(i) = rthom(i)
               IF (datathom(i) .eq. presmax) ipMAX = i
               sigma_thom(i) = qtemp(isortp(i))
               IF (sigma_thom(i) .ge. cbig) THEN
                  PRINT *, 'SIGMA_THOM missing'
                  ier = 1
                  RETURN
               ENDIF
               IF (sigma_thom(i) .lt. zero) THEN
                  sigma_thom(i) = ABS(sigma_thom(i)*datathom(i))
               ELSE
                  IF (sigma_thom(i) .gt. zero) THEN
                     sigma_thom(i) = presfac*sigma_thom(i)
                  ELSE
                     sigma_thom(i) = p_threshold*presmax
                  ENDIF
               ENDIF
            END DO
         ENDIF

!
!       THROW AWAY NOISY (SMALL) PRESSURE DATA BELOW P_THRESHOLD
!       STARTING FROM PEAK WORKING TO LARGER, SMALLER R
!
         ineg = ipmax
         ipos = ipmax
         DO WHILE(ineg.gt.1 .and.
     1            datathom(ineg-1).ge.p_threshold*presmax)
            ineg = ineg - 1
         END DO
         DO WHILE(ipos.lt.itse .and.
     1            datathom(ipos+1).ge.p_threshold*presmax)
            ipos = ipos + 1
         END DO
         itse = ipos - ineg + 1
         ioff = ineg - 1
         DO i = 1, itse
            datathom(i) = datathom(ioff+i)
         END DO
         DO i = 1, itse
            rthom(i) = rthom(ioff+i)
         END DO
         DO i = 1, itse
            sigma_thom(i) = sigma_thom(ioff+i)
         END DO
!
!       COMPUTE PRESSURE AND 1/SIGMA SPLINES IN R-SPACE (OR S-SPACE)
!       a. PRESSURE SPLINE
!
         datanorm = zero
         delse2(:itse) = one
         pcalc(:itse) = one/sigma_thom(:itse)**2
         datalsq_p(:itse) = datathom(:itse)*pcalc(:itse)
         scthom = SUM(ABS(datathom(:itse)/sigma_thom(:itse)))
         datanorm = SUM(ABS(datathom(:itse)))
         scthom = datanorm/scthom
         CALL setspline (rthom, pcalc, datalsq_p, qtemp, ythom0,
     1      y2thom0, delse2, 0.1*tensp/scthom**2, itse, natur)


!
!       FIND PRESSURE PEAK USING SMOOTHED DATA
!
         isamax = MAXLOC(ythom0(:itse))
         i = isamax(1)
         pthommax = ythom0(i)
         rthompeak = rthom(i)

         ileft = 0                    !Count data points to left of peak
         l4v(:itse) = rthom(:itse) < rthompeak
         DO i = 1, itse
            IF (l4v(i)) ileft = ileft + 1
         END DO

         IF (ipnodes .le. 0) THEN
            ipnodes = MAX(ileft + 1,itse - ileft)
            ipnodes = MIN(inode_max,ipnodes)
            ipnodes = MAX(5,ipnodes)            !At least 5 spline knots
            IF (.not.lpprof) ipnodes = 7
         ENDIF
         IF (ipnodes < 5) STOP 'MUST PICK IPNODES > 4'
         WRITE (nthreed, *)
     1      'Number of pressure-spline knots (in s-space): ', ipnodes

         ALLOCATE( nk_pa(ipnodes), nk_pb(ipnodes), ythom(ipnodes),
     1      y2thom(ipnodes), hthom(ipnodes), pknots(ipnodes) )
         IF (istat1 .ne. 0) THEN
            PRINT *, ' ISTAT1 = ', istat1, ' ALLOCATION NK_PA'
            STOP
         ENDIF

!
!       COMPUTE NODES IN SQRT(S) SPACE FOR SPLINES
!       (SEE COMMENTS ABOVE PRECEDING SKNOTS(I) CALCULATION)
!
         DO i = 1, ipnodes
            pknots(i) = REAL(i - 1,rprec)/(ipnodes - 1)
         END DO
         hthom(:ipnodes-1) = pknots(2:ipnodes) - pknots(:ipnodes-1)
!
!       COMPUTE MINOR RADII FOR DETERMINING PHIEDGE
!
         IF (lpofr) THEN
            rthommax = rthom(itse)
            rthommin = rthom(1)
            preSIN = datathom(1)
            presout = datathom(itse)
            presmin = MIN(presin,presout)
            ipreSIN = 0
            ipresout = 0
            IF (preSIN .eq. presmin) THEN
               ipreSIN = 1
               IF (presout.le.c1p5*presmin .or. presout<=0.1_dp*presmax)
     1            ipresout = 1
            ELSE
               ipresout = 1
               IF(presin.le.c1p5*presmin .or. presin<=0.1_dp*presmax)
     1            ipresin=1
            ENDIF
         ELSE
            ipreSIN = 0                        !Only USE theta=0 in pofs
            ipresout = 1
         ENDIF
      ENDIF                                   !End of IF(itse.gt.0) test

!
!       COMPUTE INDICES OF FLUX LOOP MEASUREMENTS NEEDED FOR
!       LOOP SIGNAL MATCHING
!       ALSO COMPUTE CONNECTED EXTERNAL POLOIDAL FLUX ARRAY

      IF (.not.lfreeb) THEN
         nflxs = 0
         nobser = 0
         nobd  = 0
         nbsets = 0
      END IF

      nmeasurements = imse + itse + 2 + nflxs   !!1 for diamag, 1 for edge MSE
      DO n = 1, nobser
         needflx(n) = idontneed
         iconnect(1,nobd+n) = n            !For outputting absolute flux

         IF (lasym) CYCLE

!        Save Index of up-down symmetric spatial observation points
         DO n1 = 1,n-1
           IF ((xobser(n1).eq. xobser(n)) .and.
     >         (zobser(n1).eq.-zobser(n))) needflx(n) = n1
         ENDDO
      END DO

      DO n = 1, nobd
         dsiext(n) = zero
         IF (sigma_flux(n) .lt. zero) sigma_flux(n) =
     1     ABS(sigma_flux(n)*dsiobt(n))
         IF (sigma_flux(n) .eq. zero) sigma_flux(n) = 0.0001
         DO icount = 1, 4
            index1 = iconnect(icount,n)
            tsign = SIGN(1,index1)
            IF (index1 .ne. 0) dsiext(n) = 
     1                         dsiext(n)+psiext(ABS(index1))*tsign
         END DO
      END DO

      DO n = 1,nflxs
        index1 = indxflx(n)
        IF (index1.gt.0) THEN
          plflux(index1) = dsiobt(index1) - dsiext(index1)        !n-th connected PLASMA flux
          DO icount = 1,4
            ind = ABS(iconnect(icount,index1))
            IF ((ind.gt.0).and.(needflx(ind).eq.IDONTNEED))
     1      needflx(ind) = NEEDIT
          ENDDO
        ENDIF
      ENDDO

!
!       COMPUTE INDICES OF EXTERNAL BFIELD MEASUREMENTS NEEDED
!       FOR SIGNAL MATCHING
!       FOR MULTIPLE ANGLES AT THE SAME R,Z,PHI LOCATION, IT IS
!       ASSUMED THE LOOP DATA ARE CONSECUTIVELY ORDERED
!
      DO n = 1, nbsets
         nmeasurements = nmeasurements + nbfld(n)
         DO m = 1,nbcoils(n)
            needbfld(m,n) = IDONTNEED
            IF( (m.gt.1).and.(rbcoil(m,n).eq.rbcoil(m-1,n)).and.
     1        (zbcoil(m,n).eq.zbcoil(m-1,n)) )needbfld(m,n)=ISAMECOIL
            IF( sigma_b(m,n).lt.zero )sigma_b(m,n) =
     1      ABS(sigma_b(m,n) * bbc(m,n))
            IF( sigma_b(m,n).eq.zero )sigma_b(m,n) = 0.0001

            IF (lasym) CYCLE
!       CHECK FOR ANTISYMMETRIC SITUATED COIL FOR M1 < M, N <= NSETS
            DO n1 = 1, nbsets
               DO m1 = 1, m-1
               IF( (rbcoil(m1,n1).eq.rbcoil(m,n)) .and.
     1            (zbcoil(m1,n1).eq.-zbcoil(m,n)).and.
     2            (abcoil(m1,n1).eq.abcoil(m,n)) )
     3         needbfld(m,n) = n1 + nbsets*(m1-1)
               ENDDO
            ENDDO

         ENDDO
      ENDDO

      DO n1 = 1, nbsets
        DO m1 = 1,nbfld(n1)
          index1 = indxbfld(m1,n1)
          IF( index1.gt.0 )then
!m-th PLASMA B-field
            plbfld(index1,n1) = bbc(index1,n1) - bcoil(index1,n1)
            IF( needbfld(index1,n1).eq.IDONTNEED )
     1      needbfld(index1,n1) = NEEDIT
          ENDIF
        ENDDO
      ENDDO

      END SUBROUTINE fixrecon
