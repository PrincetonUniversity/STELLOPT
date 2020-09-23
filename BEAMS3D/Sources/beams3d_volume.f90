!-----------------------------------------------------------------------
!     Module:        beams3d_volume
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          09/23/2020
!     Description:   This subroutine computes the plasma volume based
!                    on the background grid.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_volume
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE beams3d_grid
      USE beams3d_runtime, ONLY: pi2
      USE beams3d_lines, ONLY: ns_prof1
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input Variables
!        NONE
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!     Local Variables
!        i           Helper index
!-----------------------------------------------------------------------
      INTEGER ::  i, nfp, ns_local, ier
      REAL(rprec) :: dR, dZ, dphi, s1, s2, ds
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: saxis, dVds
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: dV3d, temp
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Cylindrical coordinates dV=R*dR*dZ*dPHI
      ns_local = ns_prof1*4
      nfp = NINT(pi2/phimax)
      dR = (rmax-rmin)/nr
      dZ = (zmax-zmin)/nz
      dphi = (phimax-phimin)/nphi
      ALLOCATE(dV3d(nr,nphi,nz),temp(nr,nphi,nz))
      DO i = 1, nr
         dV3d(i,:,:) = raxis(i)*dR*dZ*dphi
      END DO

      ! Calculated dVds      
      ALLOCATE(dVds(ns_local))
      dVds = 0
      DO i = 1, ns_local
         s1 = REAL(i-1)/REAL(ns_local)
         s2 = REAL(i)/REAL(ns_local)
         ds   = s2-s1
         temp = 0
         !WHERE((S_ARR =< s2) .and. (S_ARR >s1)) temp = dV3d
         WHERE(S_ARR<=s2) temp = dV3d
         WHERE(S_ARR<=s1) temp = 0
         dVds(i) = SUM(SUM(SUM(temp,DIM=3),DIM=2),DIM=1)/ds
      END DO
      dVds = dVds*nfp ! Correct for [0,2*pi/nfp] range
      DEALLOCATE(dV3d, temp)

      ! Basically now we're done becasue we have dVds as in VMEC

      ! TEST
      DO i = 1, ns_local
         s1 = REAL(i-1)/REAL(ns_local)
         s2 = REAL(i)/REAL(ns_local)
         s1 = s1 + 0.5*(s2-s1)
         CALL EZspline_interp(Vp_spl_s,s1,ds,ier)
         PRINT *,i,dVds(i),ds
      END DO

      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_volume