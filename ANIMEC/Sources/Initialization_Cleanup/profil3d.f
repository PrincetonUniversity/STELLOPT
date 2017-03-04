      SUBROUTINE profil3d(rmn, zmn, lreset, linterp)
      USE vmec_main
      USE vmec_params
      USE vsvd
      USE vspline, ONLY: sknots, pknots, hstark, hthom
      USE realspace
      USE xstuff
#ifdef _HBANGLE
      USE angle_constraints, ONLY: store_init_array
#endif
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax),
     1    INTENT(inout) ::  rmn, zmn
      LOGICAL, INTENT(in) :: lreset, linterp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: js, l, lk, lt, lz, ntype, m, n, mn
      REAL(rprec), DIMENSION(0:ntor,ntmax) :: rold, zold
      REAL(rprec) :: sm0, t1, facj, si, rax1, zax1
      INTEGER :: jcount, jk, k
!-----------------------------------------------

!
!                INDEX OF LOCAL VARIABLES
!
!        phip     radial derivative of phi/(2*pi) on half-grid
!        chip     radial derivative of chi/(2*pi) on half-grid
!        shalf    sqrt(s) ,two-dimensional array on half-grid
!        sqrts    sqrt(s), two-dimensional array on full-grid
!        wint     two-dimensional array for normalizing angle integrations
!        ireflect two-dimensional array for computing 2pi-v angle

      DO js = 1, ns
         phip(js:nrzt:ns) = phips(js)
         chip(js:nrzt:ns) = chips(js)
      END DO

      phip(nrzt+1) = 0
      faclam = 0
      sigma_an = 1         !INITIALIZE sigma FOR ANISOTROPIC PLASMA
      wint(1:nrzt:ns) = 0

      lk = 0
      DO lt = 1, ntheta3
         DO lz = 1, nzeta
            lk = lk + 1
            DO js=2,ns
               wint(js+ns*(lk-1)) = cosmui3(lt,0)/mscale(0)
            END DO
         END DO
      END DO

!
!     COMPUTE ARRAY FOR REFLECTING v = -v (ONLY needed for lasym)
!
      jcount = 0
      DO k = 1, nzeta
         jk = nzeta + 2 - k
         IF (k .eq. 1) jk = 1
         DO js = 1,ns
           jcount = jcount+1
           ireflect(jcount) = js+ns*(jk-1)           !Index for -zeta[k]
         ENDDO
      END DO

!     INDEX FOR u = -u (need for lasym integration in wrout)
      lk = 0
      IF (.NOT.ALLOCATED(uminus)) ALLOCATE (uminus(nznt))
      DO lt = 1, ntheta2
         k = ntheta1-lt+2                  
         IF (lt .eq. 1) k = 1             !u=-0 => u=0
         DO lz = 1, nzeta
            lk = lk + 1
            uminus(lk) = k                !(-u), for u = 0,pi
         END DO
      END DO

!
!     COMPUTE INITIAL R AND Z FOURIER COEFFICIENTS,
!     FROM SCALED BOUNDARY VALUES, AND SCALXC ARRAY
!     (1/SQRTS FACTOR FOR ODD M VALUES)
!

      DO js = 1, ns
         si = sqrts(js)*sqrts(js)
         sm0 = one - si
         DO ntype = 1, ntmax
            DO m = 0, mpol1
               DO n = 0, ntor
                  t1 = one/(mscale(m)*nscale(n))
                  mn = n + ntor1*m
                  l = js + ns*mn + (ntype - 1)*mns
                  IF (MOD(m,2) .eq. 0) THEN
                     scalxc(l) = one
                  ELSE
                     scalxc(l) = one/MAX(sqrts(js),sqrts(2))
                  ENDIF
!                 Do not overwrite r,z if read in from wout file AND in free bdy mode
!                 For fixed boundary, edge values MAY have been perturbed, so must execute this loop
                  IF (.not.lreset .and. lfreeb) CYCLE        
                  IF (m .eq. 0) THEN
                     IF (.not.lreset) CYCLE        !Freeze axis if read in from wout file
                     rmn(js,n,m,ntype) = rmn(js,n,m,ntype) + si*
     1                  (rmn_bdy(n,m,ntype)*t1 - rmn(ns,n,m,ntype))
                     zmn(js,n,m,ntype) = zmn(js,n,m,ntype) + si*
     1                  (zmn_bdy(n,m,ntype)*t1 - zmn(ns,n,m,ntype))
                     IF (js .eq. 1) THEN
                        rold(n,ntype) = rmn(1,n,0,ntype)
                        zold(n,ntype) = zmn(1,n,0,ntype)
                     ENDIF
                     IF (ntype .eq. rcc) rax1 = raxis_cc(n)
                     IF (ntype .eq. zcs) zax1 =-zaxis_cs(n)
                     IF (ntype .eq. rcs) rax1 =-raxis_cs(n)
                     IF (ntype .eq. zcc) zax1 = zaxis_cc(n)
                     IF (ntype.eq.rcc .or. ntype.eq.rcs) THEN
                        rmn(js,n,m,ntype) = rmn(js,n,m,ntype) + sm0*
     1                                     (rax1*t1 - rold(n,ntype))
                     END IF
                     IF (ntype.eq.zcs .or. ntype.eq.zcc) THEN
                        zmn(js,n,m,ntype) = zmn(js,n,m,ntype) + sm0*
     1                                     (zax1*t1 - zold(n,ntype))
                     END IF
                  ELSE
                     facj = sqrts(js)**m        !!TURN OFF NEXT 3 LINES IF THIS ONE ACTIVATED
!                       IF (MOD(m,2) .eq. 0) THEN
!                           facj = sqrts(js)*sqrts(js)
!                       ELSE IF (MOD(m,2) .eq. 1) THEN
!                          facj = sqrts(js)**MIN(m,3)
!                       END IF
                     rmn(js,n,m,ntype) = rmn(js,n,m,ntype) + (rmn_bdy
     1                  (n,m,ntype)*t1 - rmn(ns,n,m,ntype))*facj
                     zmn(js,n,m,ntype) = zmn(js,n,m,ntype) + (zmn_bdy
     1                  (n,m,ntype)*t1 - zmn(ns,n,m,ntype))*facj
                  ENDIF
               END DO
            END DO
         END DO
      END DO

      scalxc(1+irzloff:2*irzloff)   = scalxc(:irzloff)              !Z-components
      scalxc(1+2*irzloff:3*irzloff) = scalxc(:irzloff)              !Lamda-components

!
!     STORE PHIFAC IN XC(NEQS1) ARRAY ELEMENT
!     STORE DEL-RMSE IN XC(NEQS2) ARRAY ELEMENT
!
      xc(neqs1) = phifac;  xc(neqs2) = 0
      scalxc(neqs1) = 1
      scalxc(neqs2) = 1

#ifdef _HBANGLE
      IF (.NOT.linterp) CALL store_init_array(xc)
#endif

      IF (lrecon) THEN
!
!       COMPUTE SPLINE SEGMENTS FOR INTEGRALS OF FUNCTIONS IN SPLININT
!       THESE ARE FIXED DURING THE ITERATION SEQUENCE FOR FIXED SK,PK-NOTS,
!       BUT MAY CHANGE IF SKNOTS, PKNOTS CHANGE FOR DIFFERENT DATA FILES
!
        CALL setup_int (sknots, shalf(2), hstark, w_ia, w1_ia, u_ia,
     1       u1_ia, nk_ia, isnodes, ns1)
c-08-96 CALL setup_int(pknots,shalf(2),hthom,w_pa,w1_pa,u_pa,u1_pa,
c-08-96 >  nk_pa,ipnodes,ns1)
!
!       COMPUTE SPLINE SEGMENTS FOR INTEGRALS OF DERIVATIVES IN SPLININT
!
        CALL setup_int (sknots, sqrts, hstark, w_ib, w1_ib, u_ib,
     1       u1_ib, nk_ib, isnodes, ns)
        CALL setup_int (pknots, sqrts, hthom, w_pb, w1_pb, u_pb,
     1       u1_pb, nk_pb, ipnodes, ns)
      ENDIF

      END SUBROUTINE profil3d
