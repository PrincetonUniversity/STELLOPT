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
      USE beams3d_runtime, ONLY: pi2, handle_err, EZSPLINE_ERR
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
      REAL(rprec) :: dR, dZ, dphi, s1, s2, ds, a1, a2, detinv
      REAL(rprec), DIMENSION(4) :: fmat
      REAL(rprec), DIMENSION(4,4) :: A, B
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: slocal, dVds, Vtemp
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: dV3d, temp
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Cylindrical coordinates dV=R*dR*dZ*dPHI
      ns_local = nr
      nfp = NINT(pi2/phimax)
      dR = (rmax-rmin)/nr
      dZ = (zmax-zmin)/nz
      dphi = (phimax-phimin)/nphi
      ALLOCATE(dV3d(nr,nphi,nz),temp(nr,nphi,nz))
      DO i = 1, nr
         dV3d(i,:,:) = raxis(i)*dR*dZ*dphi
      END DO

      ! Create an Saxis
      ALLOCATE(slocal(ns_local))
      DO i = 1, ns_local
         slocal(i) = REAL(i-1)/REAL(ns_local-1)
      END DO

      ! Calculated Volume profile     
      ALLOCATE(Vtemp(ns_local))
      !PRINT *,'VOLUME'
      Vtemp = 0
      DO i = 2, ns_local
         s1 = slocal(i-1)
         s2 = slocal(i)
         ds   = s2-s1
         temp = 0
         WHERE(S_ARR<=s2) temp = dV3d
         Vtemp(i) = SUM(SUM(SUM(temp,DIM=3),DIM=2),DIM=1)*nfp
         !PRINT *,slocal(i),Vtemp(i)
      END DO
      DEALLOCATE(dV3d, temp)

      ! We fit a 3rd order polynomial to the volume
      ! f  = a0 + a1*s + a2*s*s + a3*s*s*s
      ! f' =      a1   + 2*a2*s + 3*a3*s*s
      i=2
      A(1,2) = slocal(i)
      fmat(1)   = Vtemp(i)
      i  = ns_local/4
      A(2,2) = slocal(i)
      fmat(2)   = Vtemp(i)
      i  = ns_local/2
      A(3,2) = slocal(i)
      fmat(3)   = Vtemp(i)
      i  = ns_local
      A(4,2) = slocal(i)
      fmat(4)   = Vtemp(i)
      A(:,1) = 1
      A(:,3) = A(:,2)*A(:,2)
      A(:,4) = A(:,2)*A(:,2)*A(:,2)
!      PRINT *,'MATRIX'
!      PRINT *,A
!      PRINT *,fmat

      detinv = &
      1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
       - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
       + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
       - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

      B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
      B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
      B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
      B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
      B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
      B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
      B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
      B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
      B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
      B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
      B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
      B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
      B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
      B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
      B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
      B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))

      fmat=MATMUL(B,fmat)
      !PRINT *,'COEFS'
      !PRINT *,fmat
      DEALLOCATE(Vtemp)

      ! Now create functional form of dVds
      ALLOCATE(dVds(ns_local))
      dVds = 0
      dVds(1) = fmat(2)
      DO i = 2, ns_local
         s1 = REAL(i-1)/REAL(ns_local)
         s2 = REAL(i)/REAL(ns_local)
         ds   = s2-s1
         s1 = s1 + 0.5*ds
         dVds(i) = fmat(2) + 2*fmat(3)*s1 + 3*fmat(4)*s1*s1
      END DO


      ! TEST
      !PRINT *,'============'
      !DO i = 1, ns_local
      !   s1 = REAL(i-1)/REAL(ns_local)
      !   s2 = REAL(i)/REAL(ns_local)
      !   s1 = s1 + 0.5*(s2-s1)
      !   CALL EZspline_interp(Vp_spl_s,s1,ds,ier)
      !   PRINT *,i,dVds(i),ds
      !END DO

      IF (EZspline_allocated(Vp_spl_s))   CALL EZspline_free(Vp_spl_s,ier)
      !bcs1_s=(/ 0, 0 /)
      CALL EZspline_init(Vp_spl_s,ns_local,(/ 0, 0 /),ier)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_volume',ier)
      Vp_spl_s%isHermite   = 1
      CALL EZspline_setup(Vp_spl_s,dVds,ier,EXACT_DIM=.true.)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_volume',ier)

      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_volume