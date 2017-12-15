      SUBROUTINE bextern(plascur, wint)
      USE vacmod
      USE mgrid_mod, ONLY: bvac
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), INTENT(IN) :: plascur
      REAL(dp), INTENT(IN) :: wint(nuv3)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER  :: i, k
      REAL(dp), ALLOCATABLE :: brad(:), bphi(:), bz(:)
      REAL(dp) :: tbexon, tbexoff
C-----------------------------------------------
c
c  exterior Neumann problem
c
      CALL second0(tbexon)

      IF (.not.ALLOCATED(bvac)) STOP 'BVAC unallocated in bextern'
      ALLOCATE (brad(nuv3), bphi(nuv3), bz(nuv3), stat=i)
      IF (i .ne. 0) STOP 'allocation error in bextern'

!
!     THIS ROUTINE COMPUTES THE B DOT DS ARISING FROM EXTERNAL COILS AND INTERNAL PLASMA CURRENT
!     NOTE THAT BEXN = - BEX * DS IS THE EFFECTIVE SOURCE TERM
!
!     COMPUTE B FROM COILS ON THE PLASMA BOUNDARY
!

      CALL becoil (r1b,z1b,brad,bphi,bz,bvac(1,1),bvac(1,2),bvac(1,3))

!
!     COMPUTE B (ON PLASMA BOUNDARY) FROM NET TOROIDAL PLASMA CURRENT
!     THE NET CURRENT IS MODELLED AS A WIRE AT THE MAGNETIC AXIS, AND THE
!     BIOT-SAVART LAW IS USED TO COMPUTE THE FIELD AT THE PLASMA SURFACE
!
!     USE BEXU, BEXV, BEXN AS TEMPORARY STORAGE FOR BX, BY, BZ
!
      CALL tolicu (plascur)
      CALL belicu (bexu, bexv, bexn, cosuv, sinuv, r1b, z1b)

      DO i = nuv3min, nuv3max
         brad(i) = brad(i) + bexu(i)*cosuv(i) + bexv(i)*sinuv(i)
         bphi(i) = bphi(i) - bexu(i)*sinuv(i) + bexv(i)*cosuv(i)
         bz(i) = bz(i) + bexn(i)
      END DO

!
!     COMPUTE COVARIANT COMPONENTS OF EXTERNAL FIELD: BEXU = B0 dot dx/du, 
!     BEXV = B0 dot dx/dv. HERE, BEXN = -B0*SURF_NORM CORRESPONDS TO THE
!     "exterior Neumann problem" convention of PKM (sign flipped as noted in PKM)
!     THUS, THE UNIT NORMAL SHOULD POINT INTO THE PLASMA (OUTWARD FROM VACUUM),
!     WHICH IT DOES FOR A NEGATIVE JACOBIAN (SIGNGS) SYSTEM
!
      DO i = nuv3min, nuv3max
        bexu(i) = rub(i)*brad(i) + zub(i)*bz(i)
        bexv(i) = rvb(i)*brad(i) + zvb(i)*bz(i) + r1b(i)*bphi(i)
        bexn(i) =-(brad(i)*snr(i) + bphi(i)*snv(i) + bz(i)*snz(i))
!
!     COMPUTE NORMALIZED [(2*pi)**2], READY-TO-INTEGRATE (WINT FACTOR) SOURCE TERM
!
!     NOTE: BEXN == NP*F = -B0 dot [Xu cross Xv] NP        (see PKM, Eq. 2.13)
        bexni(i) = bexn(i)*wint(i)*pi2*pi2
      END DO

      DEALLOCATE (brad, bphi, bz)

      CALL second0(tbexoff)
      bextern_time = bextern_time + (tbexoff - tbexon)

      END SUBROUTINE bextern

