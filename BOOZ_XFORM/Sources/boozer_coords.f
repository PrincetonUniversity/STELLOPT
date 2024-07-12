      SUBROUTINE boozer_coords(jrad,jsurf)
      USE booz_params
      USE booz_persistent
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: jrad, jsurf
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nparity, istat1, nv2_b, i1, nrep                        
!      INTEGER, SAVE :: jsurf = 0                                       ! moved to function call for multi-processing
      INTEGER :: ier_arr(50)
      REAL(rprec) ::  bmodv(4), bmodb(4), err(4), jacfac
      REAL(rprec) :: u_b(4), v_b(4), piu, piv
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1   r1, z1, rodd, zodd, r12, z12, p1, q1, xjac, 
     1   lt, lz, lam, wt, wz, wp
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1   cosmm, sinmm, cosnn, sinnn
C-----------------------------------------------
c       jrad         radial point where Boozer coords. are needed
c       ns           number of vmec radial grid points
c       nu_boz       number of theta points in integration
c       nv_boz       number of zeta points in integration
c       mpol         number of theta harmonics from vmec for r,z,l
c       ntor         number of zeta harmonics from vmec (no. zeta modes = 2*ntor+1) for r,z,l
c       mpol_nyq     number of theta harmonics from vmec for bsubumn, bsubvmn
c       ntor_nyq     number of zeta harmonics from vmec (no. zeta modes = 2*ntor+1) for bsubumn, bsubvmn
c       mboz         number of boozer theta harmonics
c       nboz         number of boozer zeta harmonics
c


      IF (jsurf .eq. 0) THEN
!
!        ALLOCATE GLOBAL ARRAYS
!

         CALL setup_booz (ntorsum, ns, mnmax, ohs, xmb, xnb, 
     1      sfull, scl, mboz, nboz, mnboz, nu2_b, nu_boz, 
     2      nv_boz, nfp, lasym_b)

!
!        SET UP FIXED ANGLE ARRAYS
!

         IF (lasym_b) THEN
            nu3_b = nu_boz
         ELSE
            nu3_b = nu2_b                            !!ONLY need top half of theta mesh for symmetric plasma
         END IF

         nunv = nu3_b*nv_boz

         CALL foranl (nu3_b, nv_boz, nfp, nunv, lasym_b)

         IF (lscreen) WRITE(6, 50) mboz-1, -nboz, nboz, nu_boz, nv_boz
  50     FORMAT('  0 <= mboz <= ',i4,3x,i4,' <= nboz <= ',i4,/,
     1          '  nu_boz = ',i5,' nv_boz = ',i5,//,
     1         13x,'OUTBOARD (u=0)',14x,'JS',10x,'INBOARD (u=pi)'
     2         /,77('-')/,'  v     |B|vmec    |B|booz    Error',13x,
     3         '|B|vmec    |B|booz    Error'/)

      ENDIF

      jsurf = jsurf + 1
!
!     ALLOCATE LOCAL MEMORY
!
      ier_arr = 0
      IF (ALLOCATED(r12)) DEALLOCATE(r12);
      ALLOCATE(r12(nunv),stat=ier_arr(1))
      IF (ALLOCATED(z12)) DEALLOCATE(z12); 
      ALLOCATE(z12(nunv),stat=ier_arr(2))
      IF (ALLOCATED(r1)) DEALLOCATE(r1); 
      ALLOCATE(r1(nunv),stat=ier_arr(3))
      IF (ALLOCATED(rodd)) DEALLOCATE(rodd);
      ALLOCATE(rodd(nunv),stat=ier_arr(4))
      IF (ALLOCATED(z1)) DEALLOCATE(z1); 
      ALLOCATE(z1(nunv),stat=ier_arr(5))
      IF (ALLOCATED(zodd)) DEALLOCATE(zodd);
      ALLOCATE(zodd(nunv),stat=ier_arr(6))
      IF (ALLOCATED(lt)) DEALLOCATE(lt); 
      ALLOCATE(lt(nunv),stat=ier_arr(7))
      IF (ALLOCATED(lz)) DEALLOCATE(lz); 
      ALLOCATE(lz(nunv),stat=ier_arr(8))
      IF (ALLOCATED(p1)) DEALLOCATE(p1); 
      ALLOCATE(p1(nunv),stat=ier_arr(9))
      IF (ALLOCATED(q1)) DEALLOCATE(q1); 
      ALLOCATE(q1(nunv),stat=ier_arr(10))
      IF (ALLOCATED(xjac)) DEALLOCATE(xjac); 
      ALLOCATE(xjac(nunv),stat=ier_arr(11))
      IF(ANY(ier_arr .ne. 0))STOP 'Allocation error #1 in boozer_coords'
      
!      ALLOCATE (r12(nunv), z12(nunv), r1(nunv), rodd(nunv), z1(nunv), 
!     1          zodd(nunv), lt(nunv), lz(nunv), p1(nunv), q1(nunv), 
!     2          xjac(nunv), stat=istat1 )
!      IF (istat1 .ne. 0) STOP 'Allocation error #1 in boozer_coords'

!
!     COMPUTE FOURIER COEFFICIENTS (pmn) OF THE "SOURCE" CONTRIBUTIONS (right of Eq.10)
!     OF THE BOOZER-TO-VMEC TRANSFORMATION FUNCTION P:
!
!     Theta-Booz = Theta-VMEC + Lambda + Iota*p
!     Zeta-Booz  = Zeta-VMEC  + p
!
      CALL transpmn (pmns, bsubumnc(1,jrad), bsubvmnc(1,jrad), 
     1               pmnc, bsubumns(1,jrad), bsubvmns(1,jrad), 
     2               xm_nyq, xn_nyq, gpsi, ipsi, mnmax_nyq, jrad, 
     3               lasym_b)

!
!     BEGIN CALCULATION OF BOOZER QUANTITIES AT HALF-RADIAL
!     MESH POINT JRAD
!     (ALL TRANSFORMED QUANTITIES MUST BE ON HALF-MESH FOR ACCURACY)
!
      ier_arr = 0
      IF (ALLOCATED(lam)) DEALLOCATE(lam);
      ALLOCATE(lam(nunv),stat=ier_arr(1))
      IF (ALLOCATED(wt)) DEALLOCATE(wt);
      ALLOCATE(wt(nunv),stat=ier_arr(2))
      IF (ALLOCATED(wz)) DEALLOCATE(wz);
      ALLOCATE(wz(nunv),stat=ier_arr(3))
      IF (ALLOCATED(wp)) DEALLOCATE(wp); 
      ALLOCATE(wp(nunv),stat=ier_arr(4))
      IF(ANY(ier_arr .ne. 0))STOP 'Allocation error #2 in boozer_coords'
      !ALLOCATE (lam(nunv), wt(nunv), wz(nunv), wp(nunv), stat=istat1)
      !IF (istat1 .ne. 0) STOP 'Allocation error #2 in boozer_coords'
!
!     COMPUTE EVEN (in poloidal mode number) AND ODD COMPONENTS
!     OF R,Z, LAMDA IN REAL SPACE (VMEC) COORDINATES
!
      nparity = 0
      CALL vcoords_rz (rmnc, zmns, lmns, rmns, zmnc, lmnc, xm, xn, 
     1   ntorsum, ns, jrad, mnmax, r1, z1, lt, lz, lam, sfull, 
     2   nparity, nunv, nfp, lasym_b)

      nparity = 1
      CALL vcoords_rz (rmnc, zmns, lmns, rmns, zmnc, lmnc, xm, xn, 
     1   ntorsum, ns, jrad, mnmax, rodd, zodd, lt, lz, lam, sfull, 
     2   nparity, nunv, nfp, lasym_b)

!     COMPUTE "SOURCE" PART OF TRANSFORMATION FUNCTION p==wp (RIGHT-SIDE OF EQ.10), ITS DERIVATIVES,
!     AND |B| ALL IN VMEC COORDINATES 
      CALL vcoords_w (bmodmnc(1,jrad), bmodmns(1,jrad), pmns, pmnc, 
     1                xm_nyq, xn_nyq, jrad, mnmax_nyq, bmod_b, wt, 
     2                wz, wp, nunv, nfp, lasym_b)

!
!     COMPUTE MAPPING FUNCTIONS P1 = lambda+iota*p, Q1 = p, AND MAPPING JACOBIAN (XJAC)
!     FOR DOING FOURIER INTEGRALS IN REAL SPACE (VMEC) COORDINATES
!
      CALL harfun (jacfac, hiota, gpsi, ipsi, jrad, nunv,
     1             lt, lz, lam, wt, wz, wp, p1, q1, xjac)
      DEALLOCATE (lam, wt, wz, wp, stat=istat1)
      IF (istat1 .ne. 0) STOP 'Deallocation error in boozer_coords'

!
!     COMPUTE R12, Z12 ON RADIAL HALF-GRID (IN ORIGINAL VMEC COORDINATES)
!
      nrep = 1
      CALL booz_rzhalf(r1, z1, rodd, zodd, r12, z12, ohs, 
     1                 jrad, nunv, nrep)
! 
!     Store VMEC-Space fixed point values for |B| for checking accuracy later
!
      nv2_b = nv_boz/2+1               !Index of v=pi (for non-axisymetry)
      bmodv(1) = bmod_b(1,1)           !(v=0,u=0)
      bmodv(2) = bmod_b(1,nu2_b)       !(0,pi)
      bmodv(3) = bmod_b(nv2_b,1)       !(pi,0)
      bmodv(4) = bmod_b(nv2_b,nu2_b)   !(pi,pi)

      DEALLOCATE (r1, rodd, z1, zodd, lt, lz, stat=istat1)
      IF (istat1 .ne. 0) STOP 'Deallocation error in boozer_coords'

!
!     COMPUTE BOOZER-SPACE FOURIER COEFFICIENTS FOR R,Z,P, AND |B|
!
      ier_arr = 0
      IF (ALLOCATED(cosmm)) DEALLOCATE(cosmm); 
      ALLOCATE(cosmm(nunv,0:mboz),stat=ier_arr(1))
      IF (ALLOCATED(sinmm)) DEALLOCATE(sinmm); 
      ALLOCATE(sinmm(nunv,0:mboz),stat=ier_arr(2))
      IF (ALLOCATED(cosnn)) DEALLOCATE(cosnn); 
      ALLOCATE(cosnn(nunv,0:nboz),stat=ier_arr(3))
      IF (ALLOCATED(sinnn)) DEALLOCATE(sinnn); 
      ALLOCATE(sinnn(nunv,0:nboz),stat=ier_arr(4))
      IF(ANY(ier_arr .ne. 0))STOP 'Allocation error #3 in boozer_coords'
!      ALLOCATE (cosmm(nunv,0:mboz), sinmm(nunv,0:mboz),
!     1   cosnn(nunv,0:nboz), sinnn(nunv,0:nboz), stat=istat1)
!      IF (istat1 .ne. 0) STOP 'Deallocation error in boozer_coords'

      CALL boozer (thgrd, ztgrd, bmod_b, r12, z12, xmb, xnb,
     1   bmncb(1,jsurf), rmncb(1,jsurf), zmnsb(1,jsurf),
     2   pmnsb(1,jsurf), gmncb(1,jsurf), bmnsb(1,jsurf), 
     3   rmnsb(1,jsurf), zmncb(1,jsurf), pmncb(1,jsurf), 
     4   gmnsb(1,jsurf), scl, p1, q1, xjac,
     5   cosmm, sinmm, cosnn, sinnn, mnboz, nunv, mboz, nboz,
     6   nfp, nu2_b, nv_boz, jacfac)

!
!     STORE BOOZER ANGLES CORRESPONDING TO VMEC ANGLES (u,v) OF 0 and pi
!
      piv = ztgrd(nv2_b)
      u_b(1) = p1(1);          v_b(1) = q1(1)            !(u=0,v=0)
      u_b(3) = p1(nv2_b);      v_b(3) = piv+q1(nv2_b)    !(u=0,v=pi)
      i1 = 1+nv_boz*(nu2_b-1)
      piu = thgrd(i1)
      u_b(2) = piu+p1(i1);     v_b(2) = q1(i1)           !(u=pi,v=0)
      i1 = nv2_b+nv_boz*(nu2_b-1)
      u_b(4) = piu+p1(i1);     v_b(4) = piv+q1(i1)       !(u=pi,v=pi)

      DEALLOCATE (p1, q1, xjac, stat=istat1)
      IF (istat1 .ne. 0) STOP 'Deallocation error in boozer_coords'

!
!     COMPUTE BOOZER-SPACE MOD-B FOR GENERAL CASE (LASYM = TRUE, CAN'T USE
!     FIXED-POINT SYMMETRY ANYMORE)
!
      CALL modbooz(bmncb(1,jsurf), bmnsb(1,jsurf), 
     1   bmodb, xmb, xnb, u_b, v_b, cosmm, sinmm, cosnn, sinnn,
     2   mnboz, mboz, nboz, nfp, lasym_b)

      err = ABS(bmodb - bmodv)/MAX(bmodb,bmodv)
      IF (lscreen) WRITE(6,100)
     1          bmodv(1),bmodb(1),err(1),jrad,bmodv(2),bmodb(2),err(2),
     2          bmodv(3),bmodb(3),err(3),bmodv(4),bmodb(4),err(4)
 100  FORMAT('  0  ',1p,3e11.3,i5,2x,3e11.3,/' pi  ',2(3e11.3,7x) )

!SPH ADDED 102309 - Write out flux surface shapes in Boozer coordinates (at nvplane=1 for now...)
      IF (lscreen) THEN
         ier_arr=0
         IF (ALLOCATED(r1)) DEALLOCATE(r1)
         ALLOCATE(r1(nunv), stat=ier_arr(1))
         IF (ALLOCATED(rodd)) DEALLOCATE(rodd)
         ALLOCATE(rodd(nunv), stat=ier_arr(2))
         IF (ALLOCATED(z1)) DEALLOCATE(z1)
         ALLOCATE(z1(nunv), stat=ier_arr(3))
         IF (ALLOCATED(zodd)) DEALLOCATE(zodd)
         ALLOCATE(zodd(nunv), stat=ier_arr(4))
         IF (ANY(ier_arr.ne.0))
     1        STOP 'Allocation error #3 in boozer_coords'
!         ALLOCATE (r1(nunv), rodd(nunv), z1(nunv), zodd(nunv), 
!     1          stat=istat1 )
!         IF (istat1 .ne. 0) STOP 'Allocation error #3 in boozer_coords'

         CALL vcoords_rzb (rmncb, zmnsb, rmnsb, zmncb, xmb, xnb, 
     1                  cosmm, sinmm, cosnn, sinnn, mboz, nboz,
     2                  mnboz, jsurf, ns, r1, z1, nunv, nfp, lasym_b)

         nrep = 2                 !1 for vmec, 2 for booz
         CALL booz_rzhalf(r1, z1, rodd, zodd, r12, z12, ohs, 
     1                 jrad, nunv, nrep)

         DEALLOCATE (cosmm, sinmm, cosnn, sinnn, 
     1            r1, rodd, z1, zodd, r12, z12, stat=istat1)
         IF (istat1 .ne. 0) STOP 'Deallocation error in boozer_coords'
      END IF

      END SUBROUTINE boozer_coords
