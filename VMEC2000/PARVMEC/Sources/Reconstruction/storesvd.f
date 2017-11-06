      SUBROUTINE storesvd(re, ro, lue, luo, lve, phipog, zu00)
      USE vmec_main, fpsi => bvco
      USE realspace
      USE vsvd
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,nzeta,*) :: re, ro, lue
      REAL(rprec), DIMENSION(ns,*) :: luo
      REAL(rprec), DIMENSION(*) :: lve, phipog
      REAL(rprec), DIMENSION(ns,nzeta,*) :: zu00
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: lk, l, noff, jmin, lt, js, j, ks, j1, j2
      REAL(rprec), DIMENSION(ns,ntheta2) :: luef, luof
      REAL(rprec) :: signi, phipog0, fpsie, w1,
     1    fpsi2, fpsi1, diota
C-----------------------------------------------
!
!     COMPUTE FULL MESH VALUES OF LAMBDA AT THETA=0,PI
!
      CALL lamfull (luef, luof, lue, luo)
!
!     COMPUTE PHIP*GUU*(IOTA-LV)/GSQRT AS S->0 FOR CURRENT(0)
!
        DO lk = 1,nznt
          l = 2+ns*(lk-1)
          phipog0 = 1.5_dp*phipog(l) - 0.5_dp*phipog(l+1)
          bsubu0(lk) = phipog0*( (iotas(2) + lve(l))*guu(l)
     1    + shalf(2)*luo(2,lk)*guv(l) )
        ENDDO

      IF (.not.lrecon) RETURN

      signi = SIGN(one,iotaf(2))
      noff = ntheta2
      jmin = 2
      DO lt = 1, 2
         j2   = ns*nzeta*(noff-1)
         IF (lt .eq. 1) THEN
            l = ns+1-jmin
            imid(l:1:(-1)) = (/(j1,j1=jmin+j2,ns+j2)/)
            rmid(l:1:(-1)) = re(jmin:ns,1,noff) + sqrts(jmin:ns)
     1         *ro(jmin:ns,1,noff)
            datamse(l:1:(-1)) = ATAN(iotaf(jmin:ns)*zu00(jmin:ns
     1         ,1,noff)/(rmid(l:1:(-1))*(luef(jmin:ns,noff)+
     2         sqrts(jmin:ns)*luof(jmin:ns,noff))))/dcon
            qmid(l:1:(-1)) = signi/iotaf(jmin:ns)
            presmid(l:1:(-1)) = presf(jmin:ns)/mu0
         ELSE
            l = ns+jmin-1
            imid(l:2*ns-1)=(/(j1,j1=jmin+j2,ns+j2)/)
            rmid(l:2*ns-1) = re(jmin:ns,1,noff) + sqrts(jmin:ns)
     1         *ro(jmin:ns,1,noff)
            datamse(l:2*ns-1) = ATAN(iotaf(jmin:ns)*zu00(jmin:ns
     1         ,1,noff)/(rmid(l:2*ns-1)*(luef(jmin:ns,noff)+
     2         sqrts(jmin:ns)*luof(jmin:ns,noff))))/dcon
            qmid(l:2*ns-1) = signi/iotaf(jmin:ns)
            presmid(l:2*ns-1) = presf(jmin:ns)/mu0
         ENDIF
         noff = 1
         jmin = 1
      END DO

      CALL findphi (re, ro, rthom, delse2, delso2, rmid, indexs2,
     1   indexu2, indexr, itse)
      pcalc(:itse) = (presf(indexs2(:itse))) + ABS(delse2(:itse))*(presf
     1   (indexs2(:itse)+1)-(presf(indexs2(:itse))))

!
!       SORT ON RSORT(KS) ARRAY
!       INCLUDE EDGE PITCH MATCH TO TOTAL CURRENT
!
      fpsie = 1.5_dp*fpsi(ns) - 0.5_dp*fpsi(ns1)

      DO j = 1, imse2
         ks = isorts(j)
         js = indexs1(ks)
         w1 = delse1(ks)
         IF (js .eq. ns) THEN
            fpsi2 = fpsie
            fpsi1 = fpsie
         ELSE IF (js .eq. ns1) THEN
            fpsi2 = fpsie
         ELSE
            fpsi2 = 0.5_dp*(fpsi(js+1)+fpsi(js+2))
         ENDIF
         IF (js.lt.ns .and. js.gt.1) fpsi1 = .5_dp*(fpsi(js)+fpsi(js+1))
         IF (js .eq. 1) fpsi1 = fpsi(2)
         diota = (one - w1)*iotaf(js) + w1*iotaf(js+1)
         fpsical(ks) = (one - w1)*fpsi1 + w1*fpsi2
         IF (stark_weight(ks) .ne. zero) qmeas(ks) = datastark(ks)/
     1      stark_weight(ks)
         qcalc(ks) = diota
         starkcal(ks) = diota*stark_weight(ks)
         rsort0(ks) = SQRT(hs*(js - 1 + w1))
      END DO

      END SUBROUTINE storesvd
