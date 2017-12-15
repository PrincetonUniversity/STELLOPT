      SUBROUTINE profil3d_par(rmn, zmn, lreset, linterp)
      USE vmec_main
      USE vmec_params
      USE vsvd
      USE vspline, ONLY: sknots, pknots, hstark, hthom
      USE realspace
      USE xstuff
#ifdef _HBANGLE
      USE angle_constraints, ONLY: store_init_array
#endif
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,ntmax),
     1    INTENT(inout) ::  rmn, zmn
      LOGICAL, INTENT(in) :: lreset, linterp
#if defined(SKS)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: js, l, lk, lt, lz, ntype, m, n, mn
      REAL(dp), DIMENSION(0:ntor,ntmax) :: rold, zold
      REAL(dp) :: sm0, t1, facj, si, rax1, zax1
      INTEGER :: jcount, jk, k
      INTEGER :: i, j, nsmin, nsmax, lpar
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: bcast_buf
      REAL(dp) :: tprofon, tprofoff, tbroadon, tbroadoff

!-----------------------------------------------
!                INDEX OF LOCAL VARIABLES
!
!     phip     radial derivative of phi/(2*pi) on half-grid
!     chip     radial derivative of chi/(2*pi) on half-grid
!     shalf    sqrt(s) ,two-dimensional array on half-grid
!     sqrts    sqrt(s), two-dimensional array on full-grid
!     wint     two-dimensional array for normalizing angle integrations
!     ireflect two-dimensional array for computing 2pi-v angle
!-----------------------------------------------
      CALL second0(tprofon)

      nsmin=t1lglob; nsmax=t1rglob
      DO js = nsmin, nsmax
        pphip(:,js) = phips(js)
        pchip(:,js) = chips(js)
      END DO

      faclam = 0
      sigma_an = 1         !INITIALIZE sigma FOR ANISOTROPIC PLASMA
      pwint(:,1) = 0
      pfaclam(:,:,nsmin:nsmax,:)=0
!
!     COMPUTE ARRAY FOR REFLECTING v = -v (ONLY needed for lasym)
!
      DO k = 1, nzeta
        jk = nzeta + 2 - k
        IF (k .eq. 1) jk = 1
        ireflect_par(k) = jk
      END DO

      lk = 0
      DO lt = 1, ntheta3
        DO lz = 1, nzeta
          lk = lk + 1
          pwint_ns(lk)=cosmui3(lt,0)/mscale(0)
          !DO js=MAX(2,nsmin), nsmax
          DO js=MAX(2,t1lglob), t1rglob
            !pwint(lk,js)=cosmui3(lt,0)/mscale(0)
            pwint(lk,js)=pwint_ns(lk)
          END DO
        END DO
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

      nsmin=t1lglob; nsmax=t1rglob

      rold(0:ntor,1:ntmax) = rmn(0:ntor,0,1,1:ntmax)
      zold(0:ntor,1:ntmax) = zmn(0:ntor,0,1,1:ntmax)

      IF(nranks.GT.1) THEN
        ALLOCATE(bcast_buf(0:2*ntor+1,1:ntmax))
        bcast_buf(0:ntor,1:ntmax)=rold(0:ntor,1:ntmax)
        bcast_buf(ntor+1:2*ntor+1,1:ntmax)=zold(0:ntor,1:ntmax)
        CALL second0(tbroadon)
        CALL MPI_Bcast(bcast_buf,2*(ntor+1)*ntmax,MPI_REAL8,0,
     1                NS_COMM,MPI_ERR)
        CALL second0(tbroadoff)
        broadcast_time = broadcast_time + (tbroadoff - tbroadon)
        rold(0:ntor,1:ntmax)=bcast_buf(0:ntor,1:ntmax)
        zold(0:ntor,1:ntmax)=bcast_buf(ntor+1:2*ntor+1,1:ntmax)
      END IF

      nsmin=t1lglob; nsmax=t1rglob
      DO js = nsmin, nsmax
        si = psqrts(1,js)*psqrts(1,js)
        sm0 = one - si
        DO ntype = 1, ntmax
          DO m = 0, mpol1
            DO n = 0, ntor
              t1 = one/(mscale(m)*nscale(n))
              mn = n + ntor1*m
              lpar = mn+mnsize*(js-1)+ (ntype - 1)*mns+1
              IF (MOD(m,2) .eq. 0) THEN
                pscalxc(lpar) = one
              ELSE
                pscalxc(lpar) = one/psqrts(1,MAX(2,js))
              ENDIF

              pscalxc(lpar+irzloff)=pscalxc(lpar)
              pscalxc(lpar+2*irzloff)=pscalxc(lpar)

!             Do not overwrite r,z if read in from wout file AND in free bdy mode
!             For fixed boundary, edge values MAY have been perturbed, so must execute this loop
              IF (.not.lreset .and. lfreeb) CYCLE        
              IF (m .eq. 0) THEN
                IF (.not.lreset) CYCLE        !Freeze axis if read in from wout file

                rmn(n,m,js,ntype) = rmn(n,m,js,ntype) + si*
     1                  (rmn_bdy(n,m,ntype)*t1 - rmn(n,m,ns,ntype))
                zmn(n,m,js,ntype) = zmn(n,m,js,ntype) + si*
     1                  (zmn_bdy(n,m,ntype)*t1 - zmn(n,m,ns,ntype))

                IF (ntype .eq. rcc) rax1 = raxis_cc(n)
                IF (ntype .eq. zcs) zax1 =-zaxis_cs(n)
                IF (ntype .eq. rcs) rax1 =-raxis_cs(n)
                IF (ntype .eq. zcc) zax1 = zaxis_cc(n)

                IF (ntype.eq.rcc .or. ntype.eq.rcs) THEN
                  rmn(n,m,js,ntype) = rmn(n,m,js,ntype) + sm0*
     1                                (rax1*t1 - rold(n,ntype))
                END IF
                IF (ntype.eq.zcs .or. ntype.eq.zcc) THEN
                  zmn(n,m,js,ntype) = zmn(n,m,js,ntype) + sm0*
     1                                (zax1*t1 - zold(n,ntype))
                END IF
              ELSE
                facj = psqrts(1,js)**m        !!TURN OFF NEXT 3 LINES IF THIS ONE ACTIVATED
                rmn(n,m,js,ntype) = rmn(n,m,js,ntype) + (rmn_bdy
     1                  (n,m,ntype)*t1 - rmn(n,m,ns,ntype))*facj
                zmn(n,m,js,ntype) = zmn(n,m,js,ntype) + (zmn_bdy
     1                  (n,m,ntype)*t1 - zmn(n,m,ns,ntype))*facj
              ENDIF
            END DO
          END DO
        END DO
      END DO

!     STORE PHIFAC IN XC(NEQS1) ARRAY ELEMENT
!     STORE DEL-RMSE IN XC(NEQS2) ARRAY ELEMENT

      pxc(neqs1) = phifac;  pxc(neqs2) = 0
      pscalxc(neqs1) = 1
      pscalxc(neqs2) = 1

#ifdef _HBANGLE
      IF (.NOT.linterp) CALL store_init_array(xc)
#endif

      IF (lrecon) THEN
        STOP 'Check sqrts passed to setup_int; should be psqrts'
!
!       COMPUTE SPLINE SEGMENTS FOR INTEGRALS OF FUNCTIONS IN SPLININT
!       THESE ARE FIXED DURING THE ITERATION SEQUENCE FOR FIXED SK,PK-NOTS,
!       BUT MAY CHANGE IF SKNOTS, PKNOTS CHANGE FOR DIFFERENT DATA FILES
!
        CALL setup_int (sknots, shalf(2), hstark, w_ia, w1_ia, u_ia,
     1       u1_ia, nk_ia, isnodes, ns1)
!
!       COMPUTE SPLINE SEGMENTS FOR INTEGRALS OF DERIVATIVES IN SPLININT
!
        CALL setup_int (sknots, sqrts, hstark, w_ib, w1_ib, u_ib,
     1       u1_ib, nk_ib, isnodes, ns)
        CALL setup_int (pknots, sqrts, hthom, w_pb, w1_pb, u_pb,
     1       u1_pb, nk_pb, ipnodes, ns)
      ENDIF

      CALL second0(tprofoff)
      profile3d_time = profile3d_time + (tprofoff - tprofon)

#endif
      END SUBROUTINE profil3d_par

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
#if defined(SKS)
      USE parallel_include_module
#endif      
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor,0:mpol1,ntmax),
     1    INTENT(inout) ::  rmn, zmn
      LOGICAL, INTENT(in) :: lreset, linterp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: js, l, lk, lt, lz, ntype, m, n, mn
      REAL(dp), DIMENSION(0:ntor,ntmax) :: rold, zold
      REAL(dp) :: sm0, t1, facj, si, rax1, zax1
      INTEGER :: jcount, jk, k
#if defined(SKS)
      INTEGER :: i, j, nsmin, nsmax
#endif      
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
