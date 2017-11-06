      SUBROUTINE getflux(amat_i, amat_p, data_array, r12sqr,
     1  gobser1, gobser2, orsq, r12, z12, gsqrt, kcflux)
      USE vmec_main
      USE vmec_params, ONLY: signgs
      USE realspace
      USE vsvd
      USE vspline
      USE vparams, ONLY: zero
      USE mgrid_mod, ONLY: nobser, nobd, needflx, iconnect, plflux
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER kcflux
      REAL(rprec), DIMENSION(isnodes,*) :: amat_i
      REAL(rprec), DIMENSION(ipnodes,*) :: amat_p
      REAL(rprec), DIMENSION(*) :: data_array
      REAL(rprec), DIMENSION(nrzt) :: r12sqr,
     1  gobser1, gobser2, orsq, r12, z12, gsqrt
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, l, iloop, iobsmax, iobs, index1, indexflx,
     1   isym, n1, lk
      REAL(rprec) :: t1, t2, tpisg, sign0, sumi, sump, delta
C-----------------------------------------------

!       IRESIDUE > 0, OTHERWISE FLUX_INIT NOT CALLED YET!
        kcflux = 0
        IF (iresidue.le.0 .or. (nflxs.eq.0 .and. iequi.eq.0)) RETURN

!
!       COMPUTES MATRIX ELEMENTS NEEDED TO RELATE OBSERVED FLUX
!       MEASUREMENTS TO CURRENT AND PRESSURE EXPANSION COEFFICIENTS
!       R12,Z12 ARE THE PLASMA R,Z COORDINATES AT THE HALF
!       RADIAL NODE POINTS
!


!
!       COMPUTE SYMMETRIZED PSI(R,Z)+PSI(R,-Z) FLUX "GREEN'S FUNCTION"
!
!
!       COMPUTE "GREEN'S FUNCTION" KERNEL ONLY FOR NEEDED OBSERVATION POINTS
!       (OR FOR FINAL OUTPUT IF IEQUI=1)
!
      tpisg = twopi*signgs                     !Positive volume integral
      DO n1 = 1, nobser
         isym = needflx(n1)
         IF (isym.eq.needit .or. iequi.eq.1) THEN
            CALL grnflx (r12sqr, r12, z12, gobser1, nrzt, n1)
            DO l = 2,nrzt
              gobser2(l) = gobser1(l)*orsq(l)
              gobser1(l) = gobser1(l)*gsqrt(l)
            END DO
!
!       DO INTEGRAL OVER ANGLES (ALL INTEGRALS ARE FROM THETA=0,TWOPI)
!       IM = <G/R**2>, PM = <G(1/R**2/<R**-2> - 1)>
!
            DO js = 2, ns
              sumi = zero
              sump = zero
              DO lk = js ,nrzt, ns
                sumi = sumi + gobser2(lk)*wint(lk)
                sump = sump + gobser1(lk)*wint(lk)
              ENDDO
              im(js,n1) = tpisg*sumi
              pm(js,n1) = (-tpisg*sump) + im(js,n1)/rm2(js)
            END DO

          ELSE IF( isym.ge.ISYMCOIL )THEN    !ONLY for up-down symmetric plasmas
            DO js = 2,ns
              im(js,n1) = im(js,isym)
              pm(js,n1) = pm(js,isym)
            ENDDO
          ENDIF
        ENDDO    !n1 loop
!
!       COMPUTE SPLINE MATRIX ELEMENTS BY INTEGRATING OVER RADIUS
!
        DO 2000 iloop = 0,iequi                !iequi = 0 normally, = 1 at END
          iobsmax = nflxs
          IF( iloop.eq.1 )iobsmax = nobd + nobser
          DO 1000 iobs = 1,iobsmax
            indexflx = indxflx(iobs)
            IF( iloop.eq.1 )indexflx = iobs
            IF( indexflx.le.0 )GOTO 1000
            DO js = 2,ns
              pm(js,0) = zero
              im(js,0) = zero
            ENDDO

            DO l = 1,4                !This could be halved by using symmetry
              index1 = iconnect(l,indexflx)
              IF( index1.ne.0 )then
                sign0 = 1.0
                IF( index1.lt.0 )then
                  sign0 = -sign0
                  index1 = -index1
                ENDIF
                DO js = 2,ns
                  pm(js,0) = pm(js,0) + sign0*pm(js,index1)
                  im(js,0) = im(js,0) + sign0*im(js,index1)
                ENDDO
              ENDIF
            ENDDO

            DO js = 2,ns
              pm(js,0) = pm(js,0)*ochip(js)
              im(js,0) = im(js,0)*ovrm2(js)
            ENDDO

               IF (iloop .eq. 0) THEN
                  kcflux = kcflux + 1
                  delta = plflux(indexflx)

                  CALL splinint (im, current, amat_i(1,kcflux), hstark,
     1               u_ib, u1_ib, w_ib, w1_ib, nk_ib, isnodes, intder,
     2               ns)

                  CALL splinint (pm, presint, amat_p(1,kcflux), hthom,
     1               u_pb, u1_pb, w_pb, w1_pb, nk_pb, ipnodes, intder,
     2               ns)

                  t1 = one/sigma_flux(indexflx)
                  data_array(kcflux) = t1*delta
                  amat_i(:,kcflux) = t1*amat_i(:,kcflux)
                  t2 = t1*mu0*pthommax         !!*pfac moved to getthom
                  amat_p(:,kcflux) = t2*amat_p(:,kcflux)


               ELSE            !Store plasma fluxes in PLFLUX for output
                  plflux(indexflx) = zero
                  DO js = 2, ns
                     plflux(indexflx) = plflux(indexflx) + pm(js,0)*
     1                  presph(js) + im(js,0)*(current(js)*iotaf(js)-
     2                  current(js-1)*iotaf(js-1))
                  END DO
               ENDIF
 1000     CONTINUE
 2000   CONTINUE

      END SUBROUTINE getflux
