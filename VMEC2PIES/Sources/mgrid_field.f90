!-----------------------------------------------------------------------
!     Subroutine:    mgrid_field
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/3/2012
!     Description:   This subroutine reads the mgrid file and
!                    interpolates over the grid to add the vacuum field
!                    to the pies realspace grid.
!
!     Note:          Should probably modularize this so other codes can
!                    use the mgrid file for fields.
!-----------------------------------------------------------------------
      SUBROUTINE mgrid_field(surf,load)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE mgrid_mod
      USE pies_runtime
      USE pies_realspace
      USE pies_background
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Input Parameters
!          surf       First surface to begin vacuum field calculation.
!          load       Load the vacuum field into splines.
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: surf
      INTEGER, INTENT(in) :: load
!-----------------------------------------------------------------------
!     Local Variables
!          mn           
!-----------------------------------------------------------------------
      INTEGER :: ier, i, j, u, v, ik, nv_mgrid
      INTEGER :: bcs1(2), bcs2(2), bcs3(2)
      REAL(rprec) :: pi2, x_int, y_int, z_int
      REAL(rprec) :: br(0:k,1:nu,1:nv)
      REAL(rprec) :: bphi(0:k,1:nu,1:nv)
      REAL(rprec) :: bz(0:k,1:nu,1:nv)
      REAL(rprec) :: bs_temp(0:k,1:nu,1:nv)
      REAL(rprec) :: bsreal2(0:k,1:nu,1:nv)
      REAL(rprec) :: bureal2(0:k,1:nu,1:nv)
      REAL(rprec) :: bvreal2(0:k,1:nu,1:nv)
      REAL(rprec), ALLOCATABLE :: br_vac(:,:,:), bphi_vac(:,:,:), bz_vac(:,:,:)
      REAL(rprec), ALLOCATABLE :: r_vac(:), phi_vac(:), z_vac(:)
      TYPE(EZspline3_r8), SAVE :: br_spl, bz_spl, bphi_spl
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi2 = 8 * ATAN(1._rprec)
      ! Setup splines if load == 1
      IF (load == 1) THEN
         IF (lverb) write(6,*) '-----Vacuum Field Parameters-----'
      	 nv_mgrid = nv ! This should probably be supplied by the user.
         CALL read_mgrid(TRIM(mgrid_file),extcur,nv_mgrid,nfp,.false.,ier)
         IF (ier == 9) THEN
            nv_mgrid = np0b
            nfp      = nfper0
            IF (lverb) WRITE(6,*) '!!!!! Setting nv_mgrid to',nv_mgrid,' !!!!!'
            ier = 0
            CALL read_mgrid(TRIM(mgrid_file),extcur,nv_mgrid,nfp,.false.,ier)
         END IF
         IF (ier /= 0) CALL handle_err(MGRID_ERR,mgrid_file,ier)
         IF (lverb) THEN
            write(6,'(A,A)') '     file:',TRIM(mgrid_file)
            write(6,'(A,I3,A,F5.2,A,F5.2,A)') '       nr:',nr0b,'   R=[',rminb,',',rmaxb,'] (m)'
            write(6,'(A,I3,A,F5.2,A)') '     nphi:',np0b,' PHI=[ 0.00,',pi2/nfp,'] (rad)'
            write(6,'(A,I3,A,F5.2,A,F5.2,A)') '       nz:',nz0b,'   Z=[',zminb,',',zmaxb,'] (m)'
            CALL FLUSH(6)
         END IF
         ! Recompose vacuum field for splines
         IF (np0b == 1) THEN
            ALLOCATE(br_vac(1:nr0b,1:2,1:nz0b),bphi_vac(1:nr0b,1:2,1:nz0b),&
                     bz_vac(1:nr0b,1:2,1:nz0b),STAT=ier)
            IF (ier /= 0) CALL handle_err(ALLOC_ERR,'Vacuum Fields',ier)
            ALLOCATE(r_vac(1:nr0b),phi_vac(1:2),z_vac(1:nz0b),STAT=ier)
            IF (ier /= 0) CALL handle_err(ALLOC_ERR,'Vacuum Grids',ier)
         ELSE
            ALLOCATE(br_vac(1:nr0b,1:np0b,1:nz0b),bphi_vac(1:nr0b,1:np0b,1:nz0b),&
                     bz_vac(1:nr0b,1:np0b,1:nz0b),STAT=ier)
            IF (ier /= 0) CALL handle_err(ALLOC_ERR,'Vacuum Fields',ier)
            ALLOCATE(r_vac(1:nr0b),phi_vac(1:np0b),z_vac(1:nz0b),STAT=ier)
            IF (ier /= 0) CALL handle_err(ALLOC_ERR,'Vacuum Grids',ier)
         END IF
         v = 1
         DO j = 1, np0b
            DO ik = 1, nz0b
               DO i = 1, nr0b
                  br_vac(i,j,ik) = bvac(v,1)
                  bphi_vac(i,j,ik) = bvac(v,2)
                  bz_vac(i,j,ik) = bvac(v,3)
                  v = v + 1
               END DO
            END DO
         END DO
         FORALL(i = 1:nr0b) r_vac(i) = (i-1)*(rmaxb-rminb)/(nr0b-1) + rminb
         FORALL(i = 1:nz0b) z_vac(i) = (i-1)*(zmaxb-zminb)/(nz0b-1) + zminb
         IF (np0b == 1) THEN
            phi_vac(1) = 0.0
            phi_vac(2) = pi2
            br_vac(:,2,:) = br_vac(:,1,:)
            bphi_vac(:,2,:) = bphi_vac(:,1,:)
            bz_vac(:,2,:) = bz_vac(:,1,:)
            np0b = 2
         ELSE
            FORALL(i = 1:np0b) phi_vac(i) = (i-1)*(pi2/nfp)/(np0b-1)
         END IF
         ! Now setup spline
         bcs1=(/ 0, 0/)
         bcs2=(/-1,-1/)
         bcs3=(/ 0, 0/)
         !BR_SPL
         CALL EZspline_init(br_spl,nr0b,np0b,nz0b,bcs1,bcs2,bcs3,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'BR_SPL',ier)
         br_spl%isHermite = 1
         br_spl%x1 = r_vac
         br_spl%x2 = phi_vac
         br_spl%x3 = z_vac
         CALL EZspline_setup(br_spl,br_vac,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'BR_SPL',ier)
         !BPHI_SPL
         CALL EZspline_init(bphi_spl,nr0b,np0b,nz0b,bcs1,bcs2,bcs3,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'BPHI_SPL',ier)
         bphi_spl%isHermite = 1
         bphi_spl%x1 = r_vac
         bphi_spl%x2 = phi_vac
         bphi_spl%x3 = z_vac
         CALL EZspline_setup(bphi_spl,bphi_vac,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'BPHI_SPL',ier)
         !BZ_SPL
         CALL EZspline_init(bz_spl,nr0b,np0b,nz0b,bcs1,bcs2,bcs3,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'BZ_SPL',ier)
         bz_spl%isHermite = 1
         bz_spl%x1 = r_vac
         bz_spl%x2 = phi_vac
         bz_spl%x3 = z_vac
         CALL EZspline_setup(bz_spl,bz_vac,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'BZ_SPL',ier)
         DEALLOCATE(br_vac,bphi_vac,bz_vac)
         DEALLOCATE(r_vac,phi_vac,z_vac)
         CALL free_mgrid(ier)
         IF (ier /=0) CALL handle_err(MGRID_ERR,mgrid_file,ier)
      END IF
      DO ik = surf, k
         DO u = 1, nu
            DO v = 1, nv
               x_int = rreal(ik,u,v)
               y_int = pi2*(v-1)/(nv*nfp)
               z_int = zreal(ik,u,v)
               call EZspline_interp(br_spl, x_int, y_int, z_int, br(ik,u,v), ier)
               call EZspline_interp(bphi_spl, x_int, y_int, z_int, bphi(ik,u,v), ier)
               call EZspline_interp(bz_spl, x_int, y_int, z_int, bz(ik,u,v), ier)
            END DO
         END DO
      END DO
      CALL cyl2suv(0,k,nu,nv,rreal(0:k,:,:),br(0:k,:,:),bphi(0:k,:,:),bz(0:k,:,:),bsreal2(0:k,:,:),bureal2(0:k,:,:),bvreal2(0:k,:,:),1._rprec)
      bsreal(surf:k,:,:) = bsreal(surf:k,:,:) + bsreal2(surf:k,:,:)
      bureal(surf:k,:,:) = bureal(surf:k,:,:) + bureal2(surf:k,:,:)
      bvreal(surf:k,:,:) = bvreal(surf:k,:,:) + bvreal2(surf:k,:,:)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE mgrid_field
