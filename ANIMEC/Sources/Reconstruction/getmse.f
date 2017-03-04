      SUBROUTINE getmse(amat_i, amat_p, data_array, re, ro, lue, luo,
     1   zue, zuo, phipog, kcstark)
      USE vmec_main
      USE vmec_params, ONLY: signgs
      USE realspace
      USE vsvd
      USE vspline
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER kcstark
      REAL(rprec), DIMENSION(isnodes,*) :: amat_i
      REAL(rprec), DIMENSION(ipnodes,*) :: amat_p
      REAL(rprec), DIMENSION(*) :: data_array
      REAL(rprec), DIMENSION(ns,nzeta,*) ::
     1   re, ro, lue, luo, zue, zuo
      REAL(rprec), DIMENSION(*) :: phipog
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: lt, i, js, ks, j, irnodes
      REAL(rprec), DIMENSION(ns,ntheta2) :: luef, luof
      REAL(rprec) :: dlu, dzu, guu1, edgeiota, guu2, edgefactor,
     1   wedge, t1
      REAL(rprec), SAVE :: facedge
C-----------------------------------------------

!
!       THIS SUBROUTINE COMPUTES THE LEAST-SQUARES AMATRIX AND DATA
!       ARRAYS FOR MATCHING TO THE MSE DATA AT EQUALLY-SPACED KNOTS
!       IN SQRT(PHI-FLUX) SPACE (ISNODES, INCLUDING S=0 AND S=1)
!
!       THE RANGE OF MSE DATA IS EXTENDED TO INCLUDE S=1 BY USING THE
!       CURRENT MATCHING CONDITION TO CONSTRAIN IOTA(S=1)
!
!       COMING INTO THIS ROUTINE, THE RSTARK, DATASTARK HAVE BEEN
!       PREVIOUSLY ORDERED SO RSTARK(1) < RSTARK(2) < ....
!

!
!       COMPUTE FULL MESH VALUES OF LAMBDA
!
      CALL lamfull (luef, luof, lue, luo)
!
!       COMPUTE OUTER EDGE PITCH (IOTA) TO MATCH TOTAL CURRENT
!       IOTA(EDGE) = MU0 * IPLASMA/ 2 * PI *< Guu/sqrtg > PHIP(s=1)
!       NEED THIS TO SPLINE IOTA OVER FULL S-RANGE ( 0 .le. S .le.1 )
!
      guu1 = DOT_PRODUCT(c1p5*wint(ns:nrzt:ns)*guu(ns:nrzt:ns),
     1   phipog(ns:nrzt:ns))
      guu2 = DOT_PRODUCT(cp5*wint(ns-1:nrzt-1:ns)*guu(ns-1:nrzt-1:ns),
     1   phipog(ns-1:nrzt-1:ns))

      IF (iresidue.eq.0 .or. iotaf(ns).eq.zero) THEN
         facedge = one
      ELSE IF (MOD(iter2 - iter1,ipedsvd) .eq. 0) THEN
         facedge = (guu1*iotas(ns) - guu2*iotas(ns1))/(iotaf(ns)*(guu1
     1       - guu2))
      ENDIF
      edgefactor = facedge*(guu1 - guu2)*signgs*twopi
      edgeiota = currv/edgefactor

      irnodes = MAX(0,imse) + 1
      lt = 1                                     !Outer R edge
      dlu = luef(ns,lt) + luof(ns,lt)
      wedge = (zue(ns,1,lt)+zuo(ns,1,lt))/(dlu*router)
      rstark(irnodes) = router
      datastark(irnodes) = wedge*edgeiota        !Edge pitch
      sigma_stark(irnodes) = ABS(sigma_current*wedge/edgefactor)

!
!       THROW AWAY POINTS OUTSIDE GRID
!       NOTE: IF ONLY OUTER POINT KEPT, THE CALL TO SORT IS UNNECESSARY
!
      rsort0(:irnodes) = rstark(:irnodes)
      CALL sort_data (rsort0,isortr,irnodes)
      kcstark = 0
      DO i = 1,irnodes
        j = isortr(i)
        IF( ((rsort0(i).gt.rinner) .and.
     >       (rsort0(i).le.router)) .or. (iequi.ne.0) )then
           kcstark = kcstark+1
           rsort(kcstark) = rsort0(i)                        !kcstark <= i
           starkcal(kcstark) = datastark(j)                 !sorted data array
           qcalc(kcstark) = one/sigma_stark(j)                !qcalc = sorted weight array
        ENDIF
      ENDDO

!
!       COMPUTE IOTA(0) FROM SLOPE AT RSTARK=R00
!
c04-96        kcstark = kcstark+1
c04-96        rsort(kcstark) = r00                !Magnetic axis (s=0)
c04-96        qcalc(kcstark) = 1.0/scstark
c04-96
c04-96        IF( imse.gt.0 )then
c04-96        slope0 = 1.0
c04-96        CALL splint(rstark,ystark0,y2stark0,
c04-96     >  imse,r00,dum,slope0,1)
c04-96        starkcal(kcstark) = r00*slope0*luef(1,1)/dkappa
c04-96        ELSE
c04-96c       EXTEND BOUNDARY POINTS TO INCLUDE FULL RANGE IN THETA
c04-96        starkcal(kcstark) = ai(0)
c04-96        ENDIF

!
!       FIND S,THETA INDICES FOR RSORT ARRAY
!
      CALL findphi (re, ro, rsort, delse1, delso1, rmid, indexs1,
     1   indexu1, indexr, kcstark)

!
!       COMPUTE MATRIX ELEMENTS FOR IOTA SPLINE NODES CORRESPONDING
!       TO ORDERED RSORT ARRAY ( = RSORT S )
!
      IF (kcstark .gt. nmse) STOP 'kcstark>nmse'
      CALL getspline (amat_i, sknots, hstark, delse1, hs, indexs1,
     1   isorts, kcstark, isnodes)

!
!       MATCH TO STARK MSE DATA ACCORDING TO THE FORMULA:
!
!       Bz/Btor = IOTA*Zu/[ R*(1+LAMu) ]
!
!       NOTE: QCALC, DATA = STARKCAL CORRESPOND TO RSORT_S(I)
!       WITH INDEX KS = ISORTS(I) (INDEXED ON RSORT BEFORE IT WAS SORTED)
!       SAME IS TRUE FOR DELSE,O1, INDEXS1, INDEXU1
!

      islope = 0
      DO i = 1, kcstark
c                     !Index BEFORE sorting on sknots (INDEXed on rsort)
         ks = isorts(i)
         js = indexs1(ks)
         lt = indexu1(ks)
!
!       COMPUTE WEIGHT(J) = Zu / (R * [1 + LAMu]), WHICH IS THE FACTOR
!       RELATING MSE PITCH = WEIGHT(J) * IOTA(J) TO ROTATIONAL TRANSFORM.
!
!       ON ENTRY INTO THIS LOOP,
!       QCALC = 1/SIGMA_STARK
!
         dlu = (one - delse1(ks))*luef(js,lt) + delse1(ks)*luef(js+1,lt
     1      ) + (one - delso1(ks))*sqrts(js)*luof(js,lt) + delso1(ks)*
     2      sqrts(js+1)*luof(js+1,lt)
         dzu = (one - delse1(ks))*zue(js,1,lt) + delse1(ks)*zue(js+1,1,
     1      lt) + (one - delso1(ks))*sqrts(js)*zuo(js,1,lt) + delso1(ks
     2      )*sqrts(js+1)*zuo(js+1,1,lt)
         stark_weight(ks) = dzu/(rsort(ks)*dlu)

         IF (rsort(ks) .eq. router) icurrout = i

c04-96        IF( rsort(ks).eq.r00 )THEN                        !IOTA(0)
c04-96          islope = i
c04-96          stark_weight(ks) = ABS(wedge)                  !Need in
c04-96          starkcal(ks) = weight(ks)*starkcal(ks)        !Need for
c04-96          data_array(i) = starkcal(ks) * qcalc(ks)
c04-96          amat_i(1,i) = amat_i(1,i) * stark_weight(ks) * qcalc(ks)
c04-96         ELSE
         data_array(i) = starkcal(ks)*qcalc(ks)
         t1 = stark_weight(ks)*qcalc(ks)
         amat_i(:isnodes,i) = t1*amat_i(:isnodes,i)
c04-96        ENDIF
         IF (i.eq.icurrout .and. qcalc(ks).eq.zero) STOP 'CURR ERR'
      END DO

      imse2 = kcstark
      amat_p(:,:kcstark) = zero

      END SUBROUTINE getmse
