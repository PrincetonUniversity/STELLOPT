#if defined(SKS)
      SUBROUTINE scalpot_par(bvec, amatrix, wint, ivacskip)
      USE vacmod
      USE parallel_include_module
      USE vmec_main, ONLY: nznt
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ivacskip
      REAL(rprec), INTENT(out) :: bvec(mnpd2), amatrix(mnpd2*mnpd2)
      REAL(rprec), INTENT(in) :: wint(nuv3)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ip, istore, istart, istore_max, ndim
      REAL(rprec), ALLOCATABLE :: grpmn(:), green(:), gstore(:)
      REAL(rprec), ALLOCATABLE :: greenp(:,:)

      INTEGER :: i, j, k, l, q, qq, cartx, carty, numops
      INTEGER, DIMENSION(2) :: cart_coord
      REAL(rprec) :: ton, toff, skston, skstoff

      INTEGER :: blksize
C-----------------------------------------------

      CALL second0(skston)
      IF (.not.ALLOCATED(amatsav))
     1   STOP 'AMATSAV: Allocation error in scalpot'

      ALLOCATE (grpmn(nuv3*mnpd2), stat=ip)
      IF (ip .ne. 0) STOP 'GRPMN: Allocation error in scalpot'

!
!     COMPUTE TRANFORM OF ANALYTIC SOURCE AND KERNEL
!     ON EXIT, BVEC CONTAINS THE TRANSFORM OF THE ANALYTIC SOURCE
!     AND GRPMN CONTAINS TRANSFORM OF NORMAL DERIVATIVE
!     OF THE GREENS FUNCTION [PKM, EQ.(2.15)]
!
!     FOR ivacskip != 0, USE PREVIOUSLY COMPUTED bvecsav FOR SPEED
!
      ndim = mnpd2/mnpd
      CALL second0(ton)
      CALL analyt_par (grpmn, bvec, ivacskip, ndim)
      CALL second0(toff)
      timer_vac(tanal) = timer_vac(tanal) + (toff-ton)

      IF (ivacskip .NE. 0) THEN
         bvec = bvec + bvecsav
       ELSE

         istore_max = MIN(64,nuv3)

         ALLOCATE (green(nuv), gstore(nuv), greenp(nuv,istore_max),
     1             stat=ip)
         IF (ip .NE. 0) STOP 'Allocation error in scalpot'
!         bvecsav = bvec   !Save in fouri now
         gstore  = 0
!
!        COMPUTE SURFACE INTEGRALS OF SOURCE AND GREENS FUNCTION NEEDED
!        FOR SPECTRAL DECOMPOSITION OF POTENTIAL INTEGRAL EQUATION
!        NOTE: SOURCE IS THE RHS OF EQ.(3.2), KERNEL IS THE LHS OF EQ (3.2).
!        IP IS THE INDEX OF THE PRIMED VARIABLE MESH.
!

!         istart = lcol(cartx)-1
!         istart = 0
         istart = nuv3min-1 

         CALL second0(ton)

         !SKS: Have to parallelize over ip since arrays computed in
         !surface.f like rzb2, z1b, etc but used in green.f are
         !known only within [nuv3min, nuv3max].

         DO ip = nuv3min, nuv3max
           istore = 1 + MOD(ip-nuv3min,istore_max)
!
!        COMPUTE DIFFERENCE BETWEEN THE EXACT AND ANALYTIC GREENS FUNCTION AND GRADIENT 
!        [FIRST TERMS IN EQ.(2.14, 2.16)].
!
!        BECAUSE OF THE LARGE SIZES INVOLVED (nuv3*nuv3), THIS IS DONE HERE BY STORING
!        THESE QUANTITIES - FOR ALL VALUES OF THE UNPRIMED U,V COORDINATE MESH - ON A
!        LIMITED SUBSET (ISTORE_max) OF PRIMED MESH VALUES (INDEXED BY ISTORE~IP), 
!        THE FOURIER TRANSFORM OVER THE PRIMED MESH IS "BUILT-UP" BY MULTIPLE CALLS TO FOURP
!        WITHIN THIS LOOP. 
!
           CALL second0(ton)
           CALL greenf_par (green, greenp(1,istore), ip)
           CALL second0(toff)
           timer_vac(tgreenf) = timer_vac(tgreenf) + (toff-ton)


!        PERFORM INTEGRAL (SUM) OVER PRIMED MESH OF NON-SINGULAR SOURCE TERM 
!        [(h-hsing)(u,v,u',v') == bexni(ip)*green(u,v; ip) in Eq. 2.16]
!        AND STORE IT - FOR UNPRIMED MESH VALUES - IN GSTORE

           gstore = gstore + bexni(ip)*green

!
!        PERFORM FOURIER INTEGRAL OF GRADIENT KERNEL (GREENP) OVER THE UNPRIMED MESH
!        (FOR )
!        AND STORE IN GRPMN (NOTE THAT GRPMN IS ADDED TO THE ANALYTIC PIECE IN EQ. 2.14,
!        - COMPUTED IN ANALYT - WHICH HAS THE APPROPRIATE SIN, COS FACTORS ALREADY)
!

           IF (istore.eq.istore_max .or. ip.eq.nuv3max) THEN
             CALL second0(ton)
             CALL fourp_par (grpmn, greenp, istore, istart, ip, ndim)
             CALL second0(toff)
             timer_vac(tfourp) = timer_vac(tfourp) + (toff-ton)
           END IF

         END DO

         CALL second0(ton)
         IF (vlactive) THEN
           CALL MPI_Allreduce(MPI_IN_PLACE,gstore,SIZE(gstore),
     1                        MPI_REAL8,MPI_SUM,VAC_COMM, MPI_ERR)
         END IF
         CALL second0(toff)
         allreduce_time = allreduce_time + (toff - ton)
         timer_vac(tallr) = timer_vac(tallr) + (toff-ton)
!
!        COMPUTE FOURIER INTEGRAL OF GRADIENT (GRPMN) OVER PRIMED MESH IN EQ. 2.14
!        AND SOURCE (GSTORE) OVER UNPRIMED MESH IN EQ. 2.16
!

         CALL second0(ton)
         CALL fouri_par (grpmn, gstore, amatrix, amatsav, bvec, 
     1                bvecsav, ndim)
         CALL second0(toff)
         timer_vac(tfouri) = timer_vac(tfouri) + (toff-ton)
         DEALLOCATE (green, greenp, gstore)

       ENDIF

       DEALLOCATE (grpmn, stat=ip)

       amatrix = amatsav

!
!     FINAL REDUCTION OVER VAC_COMM DONE ONCE HERE
!
       CALL second0(ton)
       IF (vlactive) THEN
         CALL MPI_Allreduce(MPI_IN_PLACE,bvec,SIZE(bvec),MPI_REAL8,
     1                      MPI_SUM,VAC_COMM,MPI_ERR)
       END IF
       CALL second0(toff)
       allreduce_time = allreduce_time + (toff - ton)
       timer_vac(tanar) = timer_vac(tanar) + (toff-ton)

!      IF (ivacskip .EQ. 0) THEN
!
!     SAVE NON-SINGULAR CONTRIBUTION TO BVEC (IN BVECSAV) - DONE IN FOURI NOW
!
!      bvecsav(:mnpd2) = bvec - bvecsav(:mnpd2)
!      END IF

       CALL second0(skstoff)
       scalpot_time = scalpot_time + (skstoff - skston)

      END SUBROUTINE scalpot_par
#endif

      SUBROUTINE scalpot(bvec, amatrix, wint, ns, ivacskip)
      USE vacmod
#if defined(SKS)
      USE vmec_dim, ONLY: nrzt 
      USE parallel_include_module
#endif
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ns, ivacskip
      REAL(rprec), INTENT(out) :: bvec(mnpd2), amatrix(mnpd2*mnpd2)
      REAL(rprec), INTENT(in) :: wint(ns,nuv3)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, ip, istore, istart, istore_max, ndim
      REAL(rprec), ALLOCATABLE :: grpmn(:), green(:), gstore(:)
      REAL(rprec), ALLOCATABLE :: greenp(:,:), tmpwint(:)

#if defined(SKS)
      INTEGER :: j, k, l, q, cartx, carty, numops
      INTEGER, DIMENSION(2) :: cart_coord
      REAL(rprec) :: skston, skstoff, ton, toff
#endif
C-----------------------------------------------
#if defined(SKS)
      CALL second0(skston)
#endif
      IF (.not.ALLOCATED(amatsav))
     1   STOP 'AMATSAV: Allocation error in scalpot'

      ALLOCATE (grpmn(nuv3*mnpd2), stat=ip)
      IF (ip .ne. 0) STOP 'GRPMN: Allocation error in scalpot'

!
!     COMPUTE TRANFORM OF ANALYTIC SOURCE AND KERNEL
!     ON EXIT, BVEC CONTAINS THE TRANSFORM OF THE ANALYTIC SOURCE
!     AND GRPMN CONTAINS TRANSFORM OF NORMAL DERIVATIVE
!     OF THE GREENS FUNCTION [PKM, EQ.(2.15)]
!
!     FOR ivacskip != 0, USE PREVIOUSLY COMPUTED bvecsav FOR SPEED
!
      ndim = mnpd2/mnpd
      CALL analyt (grpmn, bvec, ivacskip, ndim)

      IF (ivacskip .ne. 0) THEN
         bvec = bvec + bvecsav
      ELSE

         istore_max = MIN(64,nuv3)
         !istore_max = MIN(32,nuv3)

         ALLOCATE (green(nuv), gstore(nuv), greenp(nuv,istore_max))
         bvecsav = bvec
         gstore  = 0
!
!        COMPUTE SURFACE INTEGRALS OF SOURCE AND GREENS FUNCTION NEEDED
!        FOR SPECTRAL DECOMPOSITION OF POTENTIAL INTEGRAL EQUATION
!        NOTE: SOURCE IS THE RHS OF EQ.(3.2), KERNEL IS THE LHS OF EQ (3.2).
!        IP IS THE INDEX OF THE PRIMED VARIABLE MESH.
!
         istart = 0

         DO ip = 1, nuv3
            istore = 1 + MOD(ip-1,istore_max)
!
!        COMPUTE DIFFERENCE BETWEEN THE EXACT AND ANALYTIC GREENS FUNCTION AND GRADIENT 
!        [FIRST TERMS IN EQ.(2.14, 2.16)].
!
!        BECAUSE OF THE LARGE SIZES INVOLVED (nuv3*nuv3), THIS IS DONE HERE BY STORING
!        THESE QUANTITIES - FOR ALL VALUES OF THE UNPRIMED U,V COORDINATE MESH - ON A
!        LIMITED SUBSET (ISTORE_max) OF PRIMED MESH VALUES (INDEXED BY ISTORE~IP), 
!        THE FOURIER TRANSFORM OVER THE PRIMED MESH IS "BUILT-UP" BY MULTIPLE CALLS TO FOURP
!        WITHIN THIS LOOP. 
!
            CALL greenf (green, greenp(1,istore), ip)

!
!        PERFORM INTEGRAL (SUM) OVER PRIMED MESH OF NON-SINGULAR SOURCE TERM 
!        [(h-hsing)(u,v,u',v') == bexni(ip)*green(u,v; ip) in Eq. 2.16]
!        AND STORE IT - FOR UNPRIMED MESH VALUES - IN GSTORE
!
            gstore = gstore + bexni(ip)*green

!
!        PERFORM FOURIER INTEGRAL OF GRADIENT KERNEL (GREENP) OVER THE UNPRIMED MESH
!        (FOR )
!        AND STORE IN GRPMN (NOTE THAT GRPMN IS ADDED TO THE ANALYTIC PIECE IN EQ. 2.14,
!        - COMPUTED IN ANALYT - WHICH HAS THE APPROPRIATE SIN, COS FACTORS ALREADY)
!

            IF (istore.eq.istore_max .or. ip.eq.nuv3) THEN
              CALL fourp (grpmn, greenp, istore, istart, ip, ndim)
            END IF
          END DO

!
!        COMPUTE FOURIER INTEGRAL OF GRADIENT (GRPMN) OVER PRIMED MESH IN EQ. 2.14
!        AND SOURCE (GSTORE) OVER UNPRIMED MESH IN EQ. 2.16
!

         ALLOCATE (tmpwint(nuv3))
         DO i=1, nuv3
           tmpwint(i)=wint(ns,i)
         END DO

         CALL fouri (grpmn, gstore, amatrix, amatsav, bvec, tmpwint, 
     1               ndim, ns)


         DEALLOCATE (green, greenp, gstore, tmpwint)
       ENDIF
 50    FORMAT(I4,1P,F20.10)

       DEALLOCATE (grpmn)
      amatrix = amatsav

      IF (ivacskip .ne. 0) RETURN
!
!     SAVE NON-SINGULAR CONTRIBUTION TO BVEC (IN BVECSAV)
!
      bvecsav(:mnpd2) = bvec - bvecsav(:mnpd2)

#if defined(SKS)
      CALL second0(skstoff)
      s_scalpot_time =s_ scalpot_time + (skstoff - skston)
#endif

      END SUBROUTINE scalpot
