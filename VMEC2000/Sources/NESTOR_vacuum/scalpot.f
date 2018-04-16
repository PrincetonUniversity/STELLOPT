      SUBROUTINE scalpot(bvec, amatrix, ivacskip)
      USE vacmod
      USE parallel_include_module
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ivacskip
      REAL(dp), INTENT(out) :: bvec(mnpd2), amatrix(mnpd2*mnpd2)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ip, istore, istart, istore_max, ndim
      REAL(dp), ALLOCATABLE :: grpmn(:), green(:), gstore(:)
      REAL(dp), ALLOCATABLE :: greenp(:,:)
      REAL(dp) :: ton, toff, tonscal
C-----------------------------------------------
      CALL second0(tonscal)

      IF (.NOT.ALLOCATED(amatsav))
     1   STOP 'AMATSAV: Allocation error in scalpot'

      ALLOCATE (grpmn(nuv3*mnpd2), stat=ip)
      IF (ip .NE. 0) STOP 'GRPMN: Allocation error in scalpot'

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
         istart = nuv3min-1 

         !SKS: Have to parallelize over ip since arrays computed in
         !surface.f like rzb2, z1b, etc but used in green.f are
         !known only within [nuv3min, nuv3max].

         PRIMED: DO ip = nuv3min, nuv3max
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
           CALL greenf (green, greenp(1,istore), ip)


!        PERFORM INTEGRAL (SUM) OVER PRIMED MESH OF NON-SINGULAR SOURCE TERM 
!        [(h-hsing)(u,v,u',v') == bexni(ip)*green(u,v; ip) in Eq. 2.16]
!        AND STORE IT - FOR UNPRIMED MESH VALUES - IN GSTORE

           gstore = gstore + bexni(ip)*green

!
!        PERFORM FOURIER INTEGRAL OF GRADIENT KERNEL (GREENP) OVER THE UNPRIMED MESH
!        AND STORE IN GRPMN (NOTE THAT GRPMN IS ADDED TO THE ANALYTIC PIECE IN EQ. 2.14,
!        - COMPUTED IN ANALYT - WHICH HAS THE APPROPRIATE SIN, COS FACTORS ALREADY)
!

           IF (istore.EQ.istore_max .OR. ip.EQ.nuv3max)
     1        CALL fourp (grpmn, greenp, istore, istart, ip, ndim)

         END DO PRIMED

         CALL second0(ton)
#if defined(SKS)
         IF (vlactive) THEN
           CALL MPI_Allreduce(MPI_IN_PLACE,gstore,SIZE(gstore),
     1                        MPI_REAL8,MPI_SUM,VAC_COMM, MPI_ERR)
         END IF
#endif
         CALL second0(toff)
         allreduce_time = allreduce_time + (toff - ton)
         timer_vac(tallr) = timer_vac(tallr) + (toff-ton)
!
!        COMPUTE FOURIER INTEGRAL OF GRADIENT (GRPMN) OVER PRIMED MESH IN EQ. 2.14
!        AND SOURCE (GSTORE) OVER UNPRIMED MESH IN EQ. 2.16
!
         CALL fouri (grpmn, gstore, amatrix, amatsav, bvec, 
     1               bvecsav, ndim)
         DEALLOCATE (green, greenp, gstore)

       ENDIF

       DEALLOCATE (grpmn, stat=ip)

       amatrix = amatsav

!
!     FINAL REDUCTION OVER VAC_COMM DONE ONCE HERE
!
       CALL second0(ton)
#if defined(SKS)
       IF (vlactive) THEN
         CALL MPI_Allreduce(MPI_IN_PLACE,bvec,SIZE(bvec),MPI_REAL8,
     1                      MPI_SUM,VAC_COMM,MPI_ERR)
       END IF
#endif
       CALL second0(toff)
       allreduce_time = allreduce_time + (toff - ton)
       timer_vac(tanar) = timer_vac(tanar) + (toff-ton)

       scalpot_time = scalpot_time + (tonscal - toff)

      END SUBROUTINE scalpot

