      SUBROUTINE fouri (grpmn, gsource, amatrix, amatsq, bvec, 
     1                  bvecNS, ndim)
      USE vacmod
      USE parallel_include_module
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: ndim
      REAL(dp), DIMENSION(mnpd2,nuv3), INTENT(IN) :: grpmn
      REAL(dp), DIMENSION(nuv), INTENT(IN) :: gsource
      REAL(dp), DIMENSION(mnpd,mnpd,ndim**2), INTENT(OUT) :: amatrix
      REAL(dp), DIMENSION(mnpd2,mnpd2), INTENT(OUT) :: amatsq
      REAL(dp), DIMENSION(mnpd,ndim), INTENT(INOUT) :: bvec, bvecNS
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
C     interior  (int_ext=-1), exterior  (int_ext=+1)  neumann problem
      REAL(dp), PARAMETER :: int_ext = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, k, m, mn, mn0, n
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: source
      REAL(dp) :: ton, toff, tfourion, tfourioff
!-----------------------------------------------
!
!     AMATRIX(,1) = A(Sin)(Sin');  AMATRIX(,2) = A(Sin)(Cos');
!     AMATRIX(,3) = A(Cos)(Sin');  AMATRIX(,4) = A(Cos)(Cos')
!
!     ARG OF TRIG FUNCTIONS: m'u' - n'v' CORRESPONDING TO THE PRIMED MESH IN
!     PKM EQ.(2.14), AND mu  - nv (UNPRIMED MESH) IN EQ. (2.16)
!
!     ON ENTRY, GRPMN(MN,...) HAS BEEN FOURIER-TRANSFORMED OVER THE UNPRIMED
!     COORDINATES (IN FOURP), SO THE FIRST INDEX OF GRPMN CORRESPONDS TO THE FIRST 
!     INDEX OF AMATRIX. THUS, THE FOURIER TRANSFORMS OF GRPMN HERE ARE OVER THE PRIMED MESH.
!
!     IN CONTRAST, THE INTEGRAL OF THE SOURCE TERM OVER THE PRIMED MESH WAS ALREADY
!     DONE (IN SCALPOT), SO HERE THE FT ARE OVER THE UNPRIMED MESH FOR THE SOURCE. 
!  
      CALL second0(tfourion)

      ALLOCATE (source(nuv3), stat=i)
      IF (i .NE. 0) STOP 'Allocation error in fouri'
!
!     STELLARATOR-SYMMETRIZE SOURCE TERMS (with respect to u,v and -u,-v)
!     INDEX (1) IS ANTI-SYMMETRIC, INDEX (2) IS SYMMETRIC
!
!     GSOURCE = - (2pi)**2 B * dS (h - hsing) * WINT
!
!     WINT: needed to normalize integral over angles to unity 
!
      IF (lasym) THEN
         source(nuv3min:nuv3max) = onp*gsource(nuv3min:nuv3max)
      ELSE
         k = 0
         DO i = 1, nu2
            DO j = 1, nv
               k = k + 1
               IF (nuv3min.LE.k .AND. k.LE.nuv3max) THEN
                  source(k) = p5*onp*(gsource(k) - gsource(imirr(k)))
               END IF
            END DO
         END DO
      END IF

!
!     INITIALIZE RUNNING-SUM ARRAYS
!
      bvecNS = 0
      amatrix = 0
    
!
!     PERFORM M,N TRANSFORMS
!
      mn = 0
      NLOOP2: DO n = -nf, nf
      MLOOP: DO m = 0, mf
         mn = mn+1
         j = 0
         IF (m.EQ.0 .AND. n.LT.0) CYCLE
         DO i = nuv3min, nuv3max
            j = j+1
            bvecNS(mn,1) = bvecNS(mn,1)
     1                   + sinmni(j,mn)*source(i)
 
            amatrix(:,mn,1) = amatrix(:,mn,1) 
     1                      + sinmni(j,mn)*grpmn(:mnpd,i)

            IF (.NOT.lasym) CYCLE

            bvecNS(mn,2) = bvecNS(mn,2) 
     1                   + cosmni(j,mn)*source(i)

            amatrix(:,mn,2) = amatrix(:,mn,2) 
     1                      + cosmni(j,mn)*grpmn(:mnpd,i)
            amatrix(:,mn,3) = amatrix(:,mn,3) 
     1                      + sinmni(j,mn)*grpmn(mnpd+1:,i)
            amatrix(:,mn,4) = amatrix(:,mn,4) 
     1                      + cosmni(j,mn)*grpmn(mnpd+1:,i)
         END DO
      END DO MLOOP
      END DO NLOOP2

#if defined(SKS)
      IF (vlactive) THEN
        CALL second0(ton)
        CALL MPI_Allreduce(MPI_IN_PLACE,amatrix,SIZE(amatrix),MPI_REAL8,
     1                     MPI_SUM,VAC_COMM,MPI_ERR)
        CALL second0(toff)
        allreduce_time = allreduce_time + (toff - ton)
      END IF
#endif
!
!     ADD (still not reduced) ANALYTIC AND NON-SINGULAR PARTS 
!
      bvec = bvec + bvecNS

      DEALLOCATE (source, stat=i)

      mn0 = 1+mf1*nf                                            !Index of m,n=(0,0)

!SANITY CHECKS
!      IF (ANY(bvec(1:mn0-mf1:mf1,:) .NE. 0._dp)) STOP 'BVEC != 0'
!      IF (ANY(amatrix(:,1:mn0-mf1:mf1,:) .ne. 0._dp)) STOP 'AMAT1 != 0'
!      IF (ANY(amatrix(1:mn0-mf1:mf1,:,:) .ne. 0._dp)) STOP 'AMAT2 != 0'
!
!     ADD DIAGONAL TERMS TO AMATRIX [THE FIRST TERM IN EQ(3.2) OF PKM]
!
      DO mn = 1, mnpd
         amatrix(mn,mn,1) = amatrix(mn,mn,1) + pi3*int_ext
      END DO

      IF (lasym) THEN
         DO mn = 1, mnpd
            amatrix(mn,mn,4) = amatrix(mn,mn,4) + pi3*int_ext
         END DO
         amatrix(mn0,mn0,4) = amatrix(mn0,mn0,4) + pi3*int_ext     !!COS(0-0) mode *2
      END IF

!
!     PUT ELEMENTS INTO SQUARE MATRIX
!
      amatsq(:mnpd,:mnpd) = amatrix(:,:,1)                      !Sin-Sin'

      IF (lasym) THEN

      amatsq(:mnpd,1+mnpd:mnpd2) = amatrix(:,:,2)               !Sin-Cos'
      amatsq(1+mnpd:mnpd2,:mnpd) = amatrix(:,:,3)               !Cos-Sin'
      amatsq(1+mnpd:mnpd2,1+mnpd:mnpd2) = amatrix(:,:,4)        !Cos-Cos'

      END IF

      CALL second0(tfourioff)
      timer_vac(tfouri) = timer_vac(tfouri) + (tfourioff-tfourion)
      fouri_time = timer_vac(tfouri)

      END SUBROUTINE fouri

