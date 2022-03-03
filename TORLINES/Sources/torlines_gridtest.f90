!-----------------------------------------------------------------------
!     Function:      torlines_gridtest
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/3/2015
!     Description:   This is designed to test the various coordinate
!                    transforms in the code.
!-----------------------------------------------------------------------
      SUBROUTINE torlines_gridtest
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE torlines_realspace, ONLY: brreal, bphireal, bzreal, bsreal, &
                                    bureal, bvreal, rs, ru, rv, zs, zu, &
                                    zv, rreal, nrho, nu, nv
      USE torlines_runtime, ONLY: vsurf
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, ik, k1, k2,u,v
      REAL(rprec) :: z
      REAL(rprec), ALLOCATABLE :: brreal2(:,:,:),bphireal2(:,:,:),bzreal2(:,:,:)
      REAL(rprec), ALLOCATABLE :: bsreal2(:,:,:),bureal2(:,:,:),bvreal2(:,:,:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      PRINT *,'IN GRIDTEST'
      ! Test metric elements (convert to realspace) (init_wout)
      !ALLOCATE(brreal(1:k,1:nu,1:nv),bphireal(1:k,1:nu,1:nv),bzreal(1:k,1:nu,1:nv))
      !ALLOCATE(bsreal2(1:k,1:nu,1:nv),bureal2(1:k,1:nu,1:nv),bvreal2(1:k,1:nu,1:nv))
      !brreal = bureal*ru + bvreal*rv
      !bphireal = rreal*bvreal
      !bzreal = bureal*zu + bvreal*zv

      !z = 1.0
      !CALL cyl2suv(2,k,nu,nv,rreal(2:k,:,:),&
      !             brreal(2:k,:,:),bphireal(2:k,:,:),bzreal(2:k,:,:),&
      !             bsreal2(2:k,:,:),bureal2(2:k,:,:),bvreal2(2:k,:,:),z)
      !bsreal2(1,:,:) = bsreal(1,:,:) ! becasue cyl2suv doesn't work for ns=1
      !bureal2(1,:,:) = bureal(1,:,:)
      !bvreal2(1,:,:) = bvreal(1,:,:)
      !DO ik = 1, vsurf*nu*nv
      !   WRITE(327,'(I,2E20.10)') ik, bsreal2(ik),bsreal(ik)
      !   WRITE(328,'(I,2E20.10)') ik, bureal2(ik),bureal(ik)
      !   WRITE(329,'(I,2E20.10)') ik, bvreal2(ik),bvreal(ik)
      !END DO
      !PRINT *,'BS_ERROR: ',SUM(SUM(SUM(bsreal2(2:k,:,:)-bsreal(2:k,:,:),DIM=3),DIM=2))
      !PRINT *,'BU_ERROR: ',SUM(SUM(SUM(bureal2(2:k,:,:)-bureal(2:k,:,:),DIM=3),DIM=2))
      !PRINT *,'BV_ERROR: ',SUM(SUM(SUM(bvreal2(2:k,:,:)-bvreal(2:k,:,:),DIM=3),DIM=2))
      !DEALLOCATE(brreal,bphireal,bzreal)
      !DEALLOCATE(bsreal2,bureal2,bvreal2)
      ! Test metric elements (convert to flux and back) (init_external)
      k1 = vsurf+1
      k2 = nrho
      ALLOCATE(brreal2(k1:k2,1:nu,1:nv),bphireal2(k1:k2,1:nu,1:nv),bzreal2(k1:k2,1:nu,1:nv))
      ALLOCATE(bsreal2(k1:k2,1:nu,1:nv),bureal2(k1:k2,1:nu,1:nv),bvreal2(k1:k2,1:nu,1:nv))
      z = 1.0
      CALL cyl2suv(k1,k2,nu,nv,rreal(k1:k2,:,:),&
                   brreal(k1:k2,:,:),bphireal(k1:k2,:,:),bzreal(k1:k2,:,:),&
                   bsreal2(k1:k2,:,:),bureal2(k1:k2,:,:),bvreal2(k1:k2,:,:),z)
      brreal2(k1:k2,1:nu,1:nv) = bsreal2(k1:k2,1:nu,1:nv)*rs(k1:k2,1:nu,1:nv) + bureal2(k1:k2,1:nu,1:nv)*ru(k1:k2,1:nu,1:nv) + bvreal2(k1:k2,1:nu,1:nv)*rv(k1:k2,1:nu,1:nv)
      bphireal2(k1:k2,1:nu,1:nv) = rreal(k1:k2,1:nu,1:nv)*bvreal2(k1:k2,1:nu,1:nv)
      bzreal2(k1:k2,1:nu,1:nv) = bsreal2(k1:k2,1:nu,1:nv)*zs(k1:k2,1:nu,1:nv) + bureal2(k1:k2,1:nu,1:nv)*zu(k1:k2,1:nu,1:nv) + bvreal2(k1:k2,1:nu,1:nv)*zv(k1:k2,1:nu,1:nv)
      PRINT *,'got here'
      DO ik = k1,k2
         DO u = 1, nu
            DO v = 1, nv
         !WRITE(326,'(I,2E20.10)') ik, rs(ik),zs(ik)
               WRITE(327,'(I5,2E20.10)') ik, brreal2(ik,u,v),brreal(ik,u,v)
               WRITE(328,'(I5,2E20.10)') ik, bphireal2(ik,u,v),bphireal(ik,u,v)
               WRITE(329,'(I5,2E20.10)') ik, bzreal2(ik,u,v),bzreal(ik,u,v)
            END DO
         END DO
      END DO
      CALL FLUSH(327)
      CALL FLUSH(328)
      CALL FLUSH(329)
      PRINT *,'BR_ERROR: ',SUM(SUM(SUM(ABS(brreal2(k1:k2,:,:))-ABS(brreal(k1:k2,:,:)),DIM=3),DIM=2))
      PRINT *,'BP_ERROR: ',SUM(SUM(SUM(ABS(bphireal2(k1:k2,:,:))-ABS(bphireal(k1:k2,:,:)),DIM=3),DIM=2))
      PRINT *,'BZ_ERROR: ',SUM(SUM(SUM(ABS(bzreal2(k1:k2,:,:))-ABS(bzreal(k1:k2,:,:)),DIM=3),DIM=2))
      DEALLOCATE(brreal2,bphireal2,bzreal2)
      DEALLOCATE(bsreal2,bureal2,bvreal2)
      
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE torlines_gridtest
