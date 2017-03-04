!-----------------------------------------------------------------------
!     Subroutine:    virtual_casing
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This subroutine calculates the plasma response in
!                    real space coordinates through a virtual casing
!                    principle.
!-----------------------------------------------------------------------
      SUBROUTINE virtual_casing(surf)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_runtime
      USE pies_realspace
      USE pies_background
      USE virtual_casing_mod, pi2_vc => pi2
!-----------------------------------------------------------------------
!     Input Parameters
!          surf        Radial index of the plasma surface
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: surf
!-----------------------------------------------------------------------
!     Local Variables
!          pi2          2 * PI
!-----------------------------------------------------------------------
      INTEGER     :: ik,u,v, i, ntot
      REAL(rprec) :: pi2, xp, yp, zp, bx, by, bz, br, bphi
      REAL(rprec), ALLOCATABLE ::  br_temp(:,:,:), bphi_temp(:,:,:), bz_temp(:,:,:) 
      REAL(rprec), ALLOCATABLE ::  bs_temp(:,:,:), bu_temp(:,:,:), bv_temp(:,:,:) 
      DOUBLE PRECISION, ALLOCATABLE :: rmnc_temp(:,:), zmns_temp(:,:), bumnc_temp(:,:),bvmnc_temp(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: rmns_temp(:,:), zmnc_temp(:,:), bumns_temp(:,:),bvmns_temp(:,:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!----------------------------------------------------------------------- 
      WRITE(6,*) '-----Applying Virtual Casing-----'
      pi2 = 8 * ATAN(1._rprec)
      IF (nu_vc < 0.0) nu_vc = nu*10
      IF (nv_vc < 0.0) nv_vc = nv*10
      ntot = nu*nv*(k-surf)
      WRITE(6,*) '  k_sheet: ',surf,'  field_points:',ntot
      WRITE(6,*) ' nu_sheet: ',nu_vc,' nv_sheet:',nv_vc*nfp
      WRITE(6,'(a)',ADVANCE='no') 'Applying Virtual Casing '
      ALLOCATE(br_temp(0:k,1:nu,1:nv),bphi_temp(0:k,1:nu,1:nv),bz_temp(0:k,1:nu,1:nv))
      ALLOCATE(bs_temp(0:k,1:nu,1:nv),bu_temp(0:k,1:nu,1:nv),bv_temp(0:k,1:nu,1:nv))
      br_temp = 0.0; bphi_temp = 0.0; bz_temp = 0.0;
      ! Initialize the virtual casing principle
      ALLOCATE(rmnc_temp(1:mnmax,1),zmns_temp(1:mnmax,1),bumnc_temp(1:mnmax,1),bvmnc_temp(1:mnmax,1))
      rmnc_temp(1:mnmax,1)  = rmnc(1:mnmax,surf)
      zmns_temp(1:mnmax,1)  = zmns(1:mnmax,surf)
      bumnc_temp(1:mnmax,1) = bumnc(1:mnmax,surf)
      bvmnc_temp(1:mnmax,1) = bvmnc(1:mnmax,surf)
      IF (lasym) THEN
         ALLOCATE(rmns_temp(1:mnmax,1),zmnc_temp(1:mnmax,1),bumns_temp(1:mnmax,1),bvmns_temp(1:mnmax,1))
         rmns_temp(1:mnmax,1)  = rmns(1:mnmax,surf)
         zmnc_temp(1:mnmax,1)  = zmnc(1:mnmax,surf)
         bumns_temp(1:mnmax,1) = bumns(1:mnmax,surf)
         bvmns_temp(1:mnmax,1) = bvmns(1:mnmax,surf)
         CALL init_virtual_casing(mnmax,nu,nv,xm,xn,&
                                rmnc_temp,zmns_temp,nfp,&
                                BUMNC=bumnc_temp,BVMNC=bvmnc_temp,&
                                RMNS=rmns_temp,ZMNC=zmnc_temp,&
                                BUMNS=bumns_temp,BVMNS=bvmns_temp)
         !CALL init_virtual_casing(mnmax,nu_vc,nv_vc,xm,xn,&
         !                      rmnc_temp,zmns_temp,&
         !                      bumnc_temp,bvmnc_temp,nfp,&
         !                      RMNS=rmns_temp,ZMNC=zmnc_temp,&
         !                      BUMNS=bumns_temp,BVMNS=bvmns_temp)
         DEALLOCATE(rmns_temp,zmnc_temp,bumns_temp,bvmns_temp)
      ELSE
         CALL init_virtual_casing(mnmax,nu,nv,xm,xn,&
                                rmnc_temp,zmns_temp,nfp,&
                                BUMNC=bumnc_temp,BVMNC=bvmnc_temp)
         !CALL init_virtual_casing(mnmax,nu_vc,nv_vc,xm,xn,&
         !                      rmnc_temp,zmns_temp,&
         !                      bumnc_temp,bvmnc_temp,nfp)
      END IF
      DEALLOCATE(rmnc_temp,zmns_temp,bumnc_temp,bvmnc_temp)
      ! Calculate the field
      i = 1
      WRITE(6,'(a,i3,a)',ADVANCE='no') '[',0,']%'
      DO ik = surf+1, k
         DO u = 1, nu
            DO v = 1, nv
               xp = rreal(ik,u,v)*cos(pi2*(v-1)/(nfp*nv))
               yp = rreal(ik,u,v)*sin(pi2*(v-1)/(nfp*nv))
               zp = zreal(ik,u,v)
               CALL bfield_virtual_casing(xp,yp,zp,bx,by,bz)
               br   = bx * cos(pi2*(v-1)/(nfp*nv)) + by * sin(pi2*(v-1)/(nfp*nv))
               bphi = by * cos(pi2*(v-1)/(nfp*nv)) - bx * sin(pi2*(v-1)/(nfp*nv))
               br_temp(ik,u,v) = br
               bphi_temp(ik,u,v) = bphi
               bz_temp(ik,u,v) = bz
               WRITE(28,*) xp,yp,zp,bx,by,bz
               i = i + 1
            END DO
            CALL backspace_out(6,6)
            WRITE(6,'(a,i3,a)',ADVANCE='no') '[',INT(100.*REAL(i)/REAL(ntot)),']%'
            CALL FLUSH(6)
         END DO
      END DO
      CALL backspace_out(6,30)
      WRITE(6,'(30X)')
      CALL FLUSH(6)
      ! Convert from cylindrical(R,phi,Z) to toroidal coordiantes (B^s,B^u,B^v)
      CALL cyl2suv(0,k,nu,nv,rreal,br_temp,bphi_temp,bz_temp,bs_temp,bu_temp,bv_temp,1._rprec)
      bsreal(surf+1:k,:,:) = bsreal(surf+1:k,:,:) + bs_temp(surf+1:k,:,:)
      bureal(surf+1:k,:,:) = bureal(surf+1:k,:,:) + bu_temp(surf+1:k,:,:)
      bvreal(surf+1:k,:,:) = bvreal(surf+1:k,:,:) + bv_temp(surf+1:k,:,:)
      DEALLOCATE(br_temp,bphi_temp,bz_temp)
      DEALLOCATE(bs_temp,bu_temp,bv_temp)
      CALL free_virtual_casing
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE virtual_casing
