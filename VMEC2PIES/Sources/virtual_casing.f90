!!-----------------------------------------------------------------------
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
      USE virtual_casing_mod, pi_vc => pi2
!-----------------------------------------------------------------------
!     Input Parameters
!          surf        Radial index of the plasma surface
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: surf
!-----------------------------------------------------------------------
!     Local Variables
!          ik           Radial dummy index
!          u            Poloidal dummy index
!          v            Toroidal dummy index
!          i            Dummy index
!          ntot         Total number of surface gridpoints
!          pi2          2 * PI
!          xp           X at desired field location
!          yp           Y at desired field location
!          zp           Z at desired field location
!          bx           Bx
!          by           By
!          bz           Bz
!          br           Br
!          bphi         Bphi
!          br_temp      Br in (s,u,v)
!          bphi_temp    Bphi in (s,u,v)
!          bz_temp      Bz in (s,u,v)
!          bs_temp      Bs in (s,u,v) (contravariant B^s)
!          bu_temp      Bu in (s,u,v) (contravariant B^u)
!          bv_temp      Bv in (s,u,v) (contravariant B^v)
!-----------------------------------------------------------------------
      INTEGER     :: ik,u,v, i, ntot, ier
      REAL        :: xp, yp, zp, bx, by, bz, br, bphi
      REAL        :: ax, ay, az, ax2, ay2, az2  ! Can be removed if not testing SPEC
      REAL(rprec) :: pi2
      REAL, ALLOCATABLE:: rmnc_temp(:,:),rmns_temp(:,:),zmnc_temp(:,:),zmns_temp(:,:)
      REAL, ALLOCATABLE :: bumnc_temp(:,:),bumns_temp(:,:),bvmnc_temp(:,:),bvmns_temp(:,:)
      REAL(rprec), ALLOCATABLE ::  br_temp(:,:,:), bphi_temp(:,:,:), bz_temp(:,:,:) 
      REAL(rprec), ALLOCATABLE ::  bs_temp(:,:,:), bu_temp(:,:,:), bv_temp(:,:,:) 
!-----------------------------------------------------------------------
!     Begin Subroutine
!----------------------------------------------------------------------- 
      pi2 = 8 * ATAN(1._rprec)
      nu_vc = nu*4
      nv_vc = nv*4
      !IF (nu_vc < 360) nu_vc = 360
      !IF (nv_vc < 3600/nfp) nv = 3600/nfp
      ntot = nu*nv*(k-surf)
      !IF (lverb) THEN
      !   WRITE(6,*) '-----Applying Virtual Casing-----'
      !   WRITE(6,*) '  k_sheet: ',surf,'  field_points:',ntot
      !   WRITE(6,*) ' nu_sheet: ',nu_vc,' nv_sheet:',nv_vc*nfp
      !   WRITE(6,'(a)',ADVANCE='no') 'Applying Virtual Casing '
      !END IF
      ALLOCATE(br_temp(0:k,1:nu,1:nv),bphi_temp(0:k,1:nu,1:nv),bz_temp(0:k,1:nu,1:nv))
      ALLOCATE(bs_temp(0:k,1:nu,1:nv),bu_temp(0:k,1:nu,1:nv),bv_temp(0:k,1:nu,1:nv))
      br_temp = 0.0; bphi_temp = 0.0; bz_temp = 0.0;
      ! Initialize the virtual casing principle
      IF (lasym) THEN
         ALLOCATE(rmnc_temp(1:mnmax,0:k),zmns_temp(1:mnmax,0:k))
         ALLOCATE(rmns_temp(1:mnmax,0:k),zmnc_temp(1:mnmax,0:k))
         ALLOCATE(bumnc_temp(1:mnmax,0:k),bvmnc_temp(1:mnmax,0:k))
         ALLOCATE(bumns_temp(1:mnmax,0:k),bvmns_temp(1:mnmax,0:k))
         rmnc_temp  = rmnc
         rmns_temp  = rmns
         zmnc_temp  = zmnc
         zmns_temp  = zmns
         bumnc_temp = bumnc
         bumns_temp = bumns
         bvmnc_temp = bvmnc
         bvmns_temp = bvmns
         CALL init_virtual_casing_flt(mnmax,nu_vc,nv_vc,xm,xn,&
                               rmnc_temp(1:mnmax,surf),zmns_temp(1:mnmax,surf),nfp,&
                               BUMNC_FLT=bumnc_temp(1:mnmax,surf),BVMNC_FLT=bvmnc_temp(1:mnmax,surf),&
                               RMNS_FLT=rmns_temp(1:mnmax,surf),ZMNC_FLT=zmnc_temp(1:mnmax,surf),&
                               BUMNS_FLT=bumns_temp(1:mnmax,surf),BVMNS_FLT=bvmns_temp(1:mnmax,surf))
         !CALL init_virtual_casing_flt(mnmax,nu_vc,nv_vc,xm,xn,&
         !                      rmnc_temp(1:mnmax,surf),zmns_temp(1:mnmax,surf),&
         !                      bumnc_temp(1:mnmax,surf),bvmnc_temp(1:mnmax,surf),nfp,&
         !                      RMNS_FLT=rmns_temp(1:mnmax,surf),ZMNC_FLT=zmnc_temp(1:mnmax,surf),&
         !                      BUMNS_FLT=bumns_temp(1:mnmax,surf),BVMNS_FLT=bvmns_temp(1:mnmax,surf))
         DEALLOCATE(rmnc_temp,zmns_temp)
         DEALLOCATE(rmns_temp,zmnc_temp)
         DEALLOCATE(bumnc_temp,bvmnc_temp)
         DEALLOCATE(bumns_temp,bvmns_temp)
      ELSE
         ALLOCATE(rmnc_temp(1:mnmax,0:k),zmns_temp(1:mnmax,0:k))
         ALLOCATE(bumnc_temp(1:mnmax,0:k),bvmnc_temp(1:mnmax,0:k))
         rmnc_temp  = rmnc
         zmns_temp  = zmns
         bumnc_temp = bumnc
         bvmnc_temp = bvmnc
         CALL init_virtual_casing_flt(mnmax,nu_vc,nv_vc,xm,xn,&
                               rmnc_temp(1:mnmax,surf),zmns_temp(1:mnmax,surf),nfp,&
                               BUMNC_FLT=bumnc_temp(1:mnmax,surf),BVMNC_FLT=bvmnc_temp(1:mnmax,surf))
         !CALL init_virtual_casing_flt(mnmax,nu_vc,nv_vc,xm,xn,&
         !                      rmnc_temp(1:mnmax,surf),zmns_temp(1:mnmax,surf),&
         !                      bumnc_temp(1:mnmax,surf),bvmnc_temp(1:mnmax,surf),nfp)
         DEALLOCATE(rmnc_temp,zmns_temp)
         DEALLOCATE(bumnc_temp,bvmnc_temp)
      END IF
      adapt_tol = 1.0E-5
      IF (lverb) CALL virtual_casing_info(6)
      ! Calculate the field
      i = 1
      IF (lverb) WRITE(6,'(a,i3,a)',ADVANCE='no') '   CALCULATING [',0,']%'
      CALL FLUSH(6)
      DO ik = k, surf+1, -1
         DO u = 1, nu
            DO v = 1, nv
               xp = rreal(ik,u,v)*cos(pi2*(v-1)/(nfp*nv))
               yp = rreal(ik,u,v)*sin(pi2*(v-1)/(nfp*nv))
               zp = zreal(ik,u,v)
               ier = 1
               !IF (virtual_casing_dist(xp,yp,zp) < 3.*min_delta_x) THEN
                  CALL bfield_virtual_casing_adapt(xp,yp,zp,bx,by,bz,ier)
               !ELSE
               !   CALL bfield_virtual_casing(xp,yp,zp,bx,by,bz)
               !END IF
               br   = bx * cos(pi2*(v-1)/(nfp*nv)) + by * sin(pi2*(v-1)/(nfp*nv))
               bphi = by * cos(pi2*(v-1)/(nfp*nv)) - bx * sin(pi2*(v-1)/(nfp*nv))
               br_temp(ik,u,v) = br
               bphi_temp(ik,u,v) = bphi
               bz_temp(ik,u,v) = bz
               !WRITE(28,*) xp,yp,zp,bx,by,bz
               i = i + 1
            END DO
            IF (lverb) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(a,i3,a)',ADVANCE='no') '[',INT(100.*REAL(i)/REAL(ntot)),']%'
               CALL FLUSH(6)
            END IF
         END DO
      END DO
      IF (lverb) THEN
         CALL backspace_out(6,30)
         WRITE(6,'(30X)')
         CALL FLUSH(6)
      END IF
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
