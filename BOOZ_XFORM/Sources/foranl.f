      SUBROUTINE foranl(nu, nv, nfp, nunv, lasym)
      USE stel_kinds
      USE booz_persistent, ONLY: cosm_b, cosn_b, sinm_b, sinn_b,
     1    cosm_nyq, cosn_nyq, sinm_nyq, sinn_nyq, thgrd, ztgrd
      USE booz_params, ONLY: mpol1, ntor, mpol_nyq, ntor_nyq
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nu, nv, nfp, nunv
      LOGICAL, INTENT(in) :: lasym
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: lk, lt, lz, istat=0
      INTEGER :: ier_arr(16)
      REAL(rprec) :: dth, dzt, twopi
C-----------------------------------------------
      ier_arr = 0
      IF (ALLOCATED(cosm_b)) DEALLOCATE(cosm_b);
      ALLOCATE(cosm_b(nunv,0:mpol1),stat=ier_arr(1))
      IF (ALLOCATED(sinm_b)) DEALLOCATE(sinm_b);
      ALLOCATE(sinm_b(nunv,0:mpol1),stat=ier_arr(2))
      IF (ALLOCATED(cosn_b)) DEALLOCATE(cosn_b);
      ALLOCATE(cosn_b(nunv,0:ntor),stat=ier_arr(3))
      IF (ALLOCATED(sinn_b)) DEALLOCATE(sinn_b);
      ALLOCATE(sinn_b(nunv,0:ntor),stat=ier_arr(4))
      IF (ALLOCATED(cosm_nyq)) DEALLOCATE(cosm_nyq);
      ALLOCATE(cosm_nyq(nunv,0:mpol_nyq),stat=ier_arr(5))
      IF (ALLOCATED(sinm_nyq)) DEALLOCATE(sinm_nyq);
      ALLOCATE(sinm_nyq(nunv,0:mpol_nyq),stat=ier_arr(6))
      IF (ALLOCATED(cosn_nyq)) DEALLOCATE(cosn_nyq);
      ALLOCATE(cosn_nyq(nunv,0:ntor_nyq),stat=ier_arr(7))
      IF (ALLOCATED(sinn_nyq)) DEALLOCATE(sinn_nyq);
      ALLOCATE(sinn_nyq(nunv,0:ntor_nyq),stat=ier_arr(8))
      IF (ALLOCATED(thgrd)) DEALLOCATE(thgrd); 
      ALLOCATE(thgrd(nunv),stat=ier_arr(9))
      IF (ALLOCATED(ztgrd)) DEALLOCATE(ztgrd); 
      ALLOCATE(ztgrd(nunv),stat=ier_arr(10))
      IF (ANY(ier_arr .ne. 0)) STOP 'Allocation error in foranl'
!      IF (.not.ALLOCATED(cosm_b)) ALLOCATE (
!     1    cosm_b(nunv,0:mpol1), sinm_b(nunv,0:mpol1),
!     2    cosn_b(nunv,0:ntor),  sinn_b(nunv,0:ntor),
!     3    cosm_nyq(nunv,0:mpol_nyq), sinm_nyq(nunv,0:mpol_nyq),
!     4    cosn_nyq(nunv,0:ntor_nyq),  sinn_nyq(nunv,0:ntor_nyq),
!     6    thgrd(nunv), ztgrd(nunv), stat=istat)
!      IF (istat .ne. 0) STOP 'Allocation error in foranl'

      twopi = 8*ATAN(1.0_dp)
!
!     COMPUTE POLOIDAL (thgrd) AND TOROIDAL (ztgrd) ANGLES
!
      IF (lasym) THEN
         dth = twopi/nu                  !USE THIS FOR FULL 2-pi
      ELSE
         dth = twopi/(2*(nu-1))          !Half-around in theta
      END IF

      dzt = twopi/(nv*nfp)
      lk = 0

      DO lt = 1, nu
         DO lz=1, nv
           lk = lk + 1
           thgrd(lk) = (lt-1)*dth
           ztgrd(lk) = (lz-1)*dzt
          END DO
      END DO

      CALL trigfunc (thgrd, ztgrd, cosm_b, sinm_b, cosn_b, sinn_b, 
     1               mpol1, ntor, nunv)

      CALL trigfunc (thgrd, ztgrd, cosm_nyq, sinm_nyq, cosn_nyq, 
     1               sinn_nyq, mpol_nyq, ntor_nyq, nunv)

      END SUBROUTINE foranl
