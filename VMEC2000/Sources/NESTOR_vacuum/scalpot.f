      SUBROUTINE scalpot(bvec, amatrix, wint, ns, ivacskip)
      USE vacmod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ns, ivacskip
      REAL(rprec), INTENT(out) :: bvec(mnpd2), amatrix(mnpd2*mnpd2)
      REAL(rprec), INTENT(in) :: wint(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ip, istore, istart, istore_max, ndim
      REAL(rprec), ALLOCATABLE :: grpmn(:), green(:), gstore(:)
      REAL(rprec), ALLOCATABLE :: greenp(:,:)
C-----------------------------------------------
      IF (.not.ALLOCATED(amatsav))
     1   STOP 'AMATSAV: Allocation error in scalpot'

      ALLOCATE (grpmn(nuv2*mnpd2), stat=ip)
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

         istore_max = MIN(64,nuv2)

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
         DO ip = 1, nuv2
            istore = 1 + MOD(ip-1,istore_max)
!
!        COMPUTE DIFFERENCE BETWEEN THE EXACT AND ANALYTIC GREENS FUNCTION AND GRADIENT 
!        [FIRST TERMS IN EQ.(2.14, 2.16)].
!
!        BECAUSE OF THE LARGE SIZES INVOLVED (NUV2*NUV2), THIS IS DONE HERE BY STORING
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
            IF (istore.eq.istore_max .or. ip.eq.nuv2)
     1        CALL fourp (grpmn, greenp, istore, istart, ip, ndim)

         END DO

!
!        COMPUTE FOURIER INTEGRAL OF GRADIENT (GRPMN) OVER PRIMED MESH IN EQ. 2.14
!        AND SOURCE (GSTORE) OVER UNPRIMED MESH IN EQ. 2.16
!
         CALL fouri (grpmn, gstore, amatrix, amatsav, bvec, wint, 
     1               ndim, ns)
         DEALLOCATE (green, greenp, gstore)
      ENDIF

      DEALLOCATE (grpmn)
      amatrix = amatsav

      IF (ivacskip .ne. 0) RETURN
!
!     SAVE NON-SINGULAR CONTRIBUTION TO BVEC (IN BVECSAV)
!
      bvecsav(:mnpd2) = bvec - bvecsav(:mnpd2)

      END SUBROUTINE scalpot
