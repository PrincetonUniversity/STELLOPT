      SUBROUTINE getbfld(amat_i, amat_p, data_array, r12sqr,
     1  gsqrt, orsq, gobsr1, gobsz1, r12, z12, kcbfld)
      USE vmec_main
      USE vmec_params, ONLY: signgs
      USE vsvd
      USE realspace, ONLY: wint
      USE vspline, ONLY: hthom, hstark
      USE mgrid_mod, ONLY: nbfldn, nbsets, nbcoils, needbfld, abcoil,
     1                     plbfld
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER kcbfld
      REAL(rprec), DIMENSION(isnodes,*) :: amat_i
      REAL(rprec), DIMENSION(ipnodes,*) :: amat_p
      REAL(rprec), DIMENSION(*) :: data_array
      REAL(rprec), DIMENSION(nrzt) :: r12sqr,
     1  gsqrt, orsq, gobsr1, gobsz1, r12, z12
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, n1, m1, iobs, isym, iobsmax, iloop,
     1   msym, nsym, indexbfld, l, lk
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: gobsz2, gobsr2
      REAL(rprec):: tpisg, sumir, sumiz, wscaleb, sumpr,
     1    sumpz, deltab, plasbfld, coscoil, sincoil, t2, t1
C----------------------------------------------

!     IRESIDUE > 0, OTHERWISE FLUX_INIT NOT CALLED YET
      kcbfld = 0
      IF (iresidue.le.0 .or. (nbfldn.eq.0 .and. iequi.eq.0))return

      ALLOCATE( gobsz2(nrzt), gobsr2(nrzt), stat=l)
      IF (l .ne. 0) STOP 'allocation problem in getbfld'
!
!       COMPUTE "GREEN'S FUNCTION" KERNEL ONLY FOR NEEDED OBSERVATION POINTS
!       (OR FOR FINAL OUTPUT IF IEQUI=1)
!       GOBSR1 = BR, GOBSZ1 = BZ
!       IF A (BR,BZ) OR (BRHO,BTHETA) PAIR, JUST CALL GRNBFLD ONCE
!
      IF (ANY(rm2(2:ns) .eq. zero)) STOP 'rm2 = 0'
      tpisg = twopi * signgs                !Positive volume integral
      DO n1 = 1, nbsets
         DO m1 = 1, nbcoils(n1)
            isym = needbfld(m1,n1)
            IF (isym.eq.needit .or. iequi.eq.1) THEN
               CALL grnbfld(r12sqr,r12,z12,gobsr1,gobsz1,nrzt,m1,n1)
               DO l = 2,nrzt
                 gobsr2(l) = gobsr1(l)*orsq(l)
                 gobsr1(l) = gobsr1(l)*gsqrt(l)
                 gobsz2(l) = gobsz1(l)*orsq(l)
                 gobsz1(l) = gobsz1(l)*gsqrt(l)
               END DO
!
!       DO INTEGRAL OVER ANGLES (ALL INTEGRALS ARE FROM THETA=0,TWOPI)
!
               DO js = 2, ns
                  sumir = zero
                  sumiz = zero
                  sumpr = zero
                  sumpz = zero
                  DO lk = js,nrzt,ns
                    sumir = sumir + gobsr2(lk)*wint(lk)
                    sumpr = sumpr + gobsr1(lk)*wint(lk)
                    sumiz = sumiz + gobsz2(lk)*wint(lk)
                    sumpz = sumpz + gobsz1(lk)*wint(lk)
                  ENDDO
                  imb(js,m1,n1,1) = tpisg*sumir
                  pmb(js,m1,n1,1) = (-tpisg*sumpr) + imb(js,m1,n1,1)
     1               /rm2(js)
                  imb(js,m1,n1,2) = tpisg*sumiz
                  pmb(js,m1,n1,2) = (-tpisg*sumpz) + imb(js,m1,n1,2)
     1               /rm2(js)
               END DO

            ELSE IF (isym .eq. ISAMECOIL) THEN       !Same coil position as previous coil
              DO js = 2,ns
                imb(js,m1,n1,1) = imb(js,m1-1,n1,1)
                pmb(js,m1,n1,1) = pmb(js,m1-1,n1,1)
                imb(js,m1,n1,2) = imb(js,m1-1,n1,2)
                pmb(js,m1,n1,2) = pmb(js,m1-1,n1,2)
              ENDDO
            ENDIF
          ENDDO   !m1
        ENDDO   !n1

!
!       CHECK FOR SYMMETRIC COIL (MAY BE IN DIFFERENT COIL SET,
!       SO HAD TO MOVE OUT OF M1,N1 LOOP ABOVE)
!
        DO n1 = 1,nbsets
          DO m1 = 1,nbcoils(n1)
            isym = needbfld(m1,n1)
            IF (isym .ge. ISYMCOIL) THEN
              msym = 1 + (isym-1)/nbsets
              nsym = 1 + MOD(isym-1,nbsets)
              DO js = 2,ns                     !BR(-Z) = -BR(Z), BZ(-Z) = BZ(Z)
                imb(js,m1,n1,1) =-imb(js,msym,nsym,1)
                pmb(js,m1,n1,1) =-pmb(js,msym,nsym,1)
                imb(js,m1,n1,2) = imb(js,msym,nsym,2)
                pmb(js,m1,n1,2) = pmb(js,msym,nsym,2)
              ENDDO
            ENDIF
          ENDDO
        ENDDO

!
!       COMPUTE SPLINE MATRIX ELEMENTS BY INTEGRATING OVER RADIUS
!
        DO 2000 iloop = 0,iequi                !iequi = 0 normally, = 1 at END
          DO n1 = 1, nbsets
            iobsmax = nbfld(n1)
            IF (iloop .eq. 1) iobsmax = nbcoils(n1)
            IF (iobsmax .gt. 0) THEN
              DO 1000 iobs = 1, iobsmax
                indexbfld = indxbfld(iobs,n1)
                IF (iloop .eq. 1) indexbfld = iobs
                IF (indexbfld .le. 0) GOTO 1000
                coscoil = COS( abcoil(indexbfld,n1) )
                sincoil = SIN( abcoil(indexbfld,n1) )
                DO js = 2,ns
                  pmb(js,0,n1,1) = ochip(js) *
     >            (pmb(js,indexbfld,n1,1)*coscoil +
     >             pmb(js,indexbfld,n1,2)*sincoil)
                  imb(js,0,n1,1) = ovrm2(js) *
     >            (imb(js,indexbfld,n1,1)*coscoil +
     >             imb(js,indexbfld,n1,2)*sincoil)
                END DO

                     IF (iloop .eq. 0) THEN
                        deltab = plbfld(indexbfld,n1)
                        kcbfld = kcbfld + 1

                        CALL splinint (imb(1,0,n1,1), current,
     1                     amat_i(1,kcbfld), hstark, u_ib, u1_ib,
     2                     w_ib, w1_ib, nk_ib, isnodes, intder, ns)

                        CALL splinint (pmb(1,0,n1,1), presint,
     1                     amat_p(1,kcbfld), hthom, u_pb, u1_pb, w_pb,
     2                     w1_pb, nk_pb, ipnodes, intder, ns)

                        wscaleb = one/sigma_b(indexbfld,n1)
                        data_array(kcbfld) = wscaleb*deltab
                        t2 = mu0*pthommax      !!*pfac moved to getthom

                        amat_i(:,kcbfld) = wscaleb*amat_i(:,kcbfld)
                        wscaleb = wscaleb*t2
                        amat_p(:,kcbfld) = wscaleb*amat_p(:,kcbfld)

                     ELSE      !Store plasma fluxes in EXTFLX for output
                        plasbfld = zero
                        DO js = 2, ns
                           t1 = current(js)*iotaf(js) - current(js-1)*
     1                        iotaf(js - 1)
                           plasbfld = plasbfld + pmb(js,0,n1,1)*
     1                        presph(js) + imb(js,0,n1,1)*t1
                        END DO
                        plbfld(iobs,n1) = plasbfld
                     ENDIF
 1000         CONTINUE
            ENDIF
          ENDDO       !n1
 2000   CONTINUE

      DEALLOCATE( gobsz2, gobsr2, stat=l)

      END SUBROUTINE getbfld
