      SUBROUTINE belicu(bx, by, bz, cos1, sin1, rp, zp)
      USE vacmod
      USE biotsavart
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(nuv3), INTENT(in)  :: cos1, sin1, rp, zp
      REAL(dp), DIMENSION(nuv3), INTENT(out) :: bx, by, bz
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: j
      REAL(dp), DIMENSION(3) :: xpt, bvec
      REAL(dp) :: tbelon, tbeloff
!-----------------------------------------------
      CALL second0(tbelon)

      DO j = nuv3min, nuv3max
         xpt(1) = rp(j) * cos1(j)
         xpt(2) = rp(j) * sin1(j)
         xpt(3) = zp(j) 
         CALL bsc_b (single_coil, xpt, bvec)
         bx(j) = bvec(1);  by(j) = bvec(2);  bz(j) = bvec(3)
      END DO
      CALL cleanup_biotsavart

      CALL second0(tbeloff)
      belicu_time = belicu_time + (tbeloff - tbelon)

      END SUBROUTINE belicu
