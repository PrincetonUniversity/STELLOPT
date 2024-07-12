!-----------------------------------------------------------------------
!     Module:        fieldlines_write_emc3
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/22/2012
!     Description:   This subroutine outputs the fieldline data to an
!                    HDF5 file or binary file.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_write_emc3
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_runtime, ONLY: id_string, lverb, lvac, pi2, &
                                    lafield_only, lbfield_only
      USE fieldlines_grid, ONLY: nr, nphi, nz, B_R, B_PHI, B_Z, raxis, &
                                 zaxis, phiaxis, nfp_m
      USE read_wout_mod, ONLY: ns, mnmax_nyq, xm_nyq, xn_nyq,rmnc, zmns,&
                               nfp
      USE safe_open_mod, ONLY: safe_open
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          iunit        File ID
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: u,v,ier, iunit, mn, i, j, ik, k, nexternal, mn0
      INTEGER :: bcs1(2)
      INTEGER, ALLOCATABLE :: xm_temp(:), xn_temp(:)
      REAL(rprec) :: dr, f0_temp
      REAL(rprec), ALLOCATABLE :: rho(:), rho_vmec(:), rho_ext(:), &
                                  ftemp(:),fmn(:)
      REAL(rprec), ALLOCATABLE :: rmnc1(:),zmns1(:),rmnc2(:),zmns2(:)
      REAL(rprec), ALLOCATABLE :: rmns1(:),zmnc1(:),rmns2(:),zmnc2(:)
      DOUBLE PRECISION, ALLOCATABLE :: xu(:),xv(:)
      DOUBLE PRECISION, ALLOCATABLE :: r(:,:,:),z(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: rmnc_sav(:,:),zmns_sav(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: rmns_sav(:,:),zmnc_sav(:,:)
      TYPE(EZspline1_r8) :: f_spl

      INTEGER, PARAMETER :: NU = 360
      INTEGER, PARAMETER :: NV = 12
      INTEGER, PARAMETER :: N_EXT = 20
      REAL(rprec), PARAMETER :: bound_separation = 1.25
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (lverb) THEN
         WRITE(6,'(A)')  '----- EMC3-EIRENE Output -----'
      END IF
      ! Write the vmec information
      IF (.not. lvac) THEN
         ! Calculate gridpoints
         k = ns + N_EXT
         bcs1=(/0,0/)
         ! Setup Splines
         CALL EZspline_init(f_spl,ns,bcs1,ier)
         ALLOCATE(xu(nu),xv(nv),r(NU,NV,k),z(NU,NV,k))
         ALLOCATE(rho(k),rho_vmec(ns),fmn(k),ftemp(ns))
         ALLOCATE(rmnc_sav(mnmax_nyq,k),zmns_sav(mnmax_nyq,k))
         ALLOCATE(xm_temp(mnmax_nyq),xn_temp(mnmax_nyq))
         ! Setup RHO Arrays
         FORALL(ik = 1:ns) rho_vmec(ik) =  SQRT(REAL(ik-1) / REAL(ns-1))
         f_spl%isHermite = 1
         f_spl%x1 = rho_vmec
         FORALL(ik = 1:k) rho(ik) =  REAL(ik-1) / REAL(k-1)
         rho = rho * bound_separation  ! Scale rho
         j=k
         DO ik = 1, k
            IF (rho(ik) > 1.0) THEN
               rho(ik) = 1.0
               j = ik
               nexternal = k-j
               EXIT
            END IF
         END DO
         FORALL(ik = 1:k) rho(ik) =  REAL(ik-1) / REAL(j-1)
         ! Setup fourier index arrays
         xm_temp = INT(xm_nyq)
         xn_temp = INT(xn_nyq)
         nfp_m = nfp
         ! Spline to rho axis (all quantities now on full grid)
         DO mn = 1, mnmax_nyq
            IF (xm_nyq(mn) == 0 .and. xn_nyq(mn) == 0) mn0 = mn
            ! RMNC
            f0_temp = rmnc(mn,1)
            ftemp(1:ns) = (rmnc(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm_nyq(mn))/2.+1)
            ftemp(1) = 0.0
            CALL EZspline_setup(f_spl,ftemp,ier)
            !IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup',ier)
            CALL EZspline_interp(f_spl,k,rho,fmn,ier)
            rmnc_sav(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm_nyq(mn))-2)
            ! ZMNS
            f0_temp = zmns(mn,1)
            ftemp(1:ns) = (zmns(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm_nyq(mn))/2.+1)
            ftemp(1) = 0.0
            CALL EZspline_setup(f_spl,ftemp,ier)
            !IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup',ier)
            CALL EZspline_interp(f_spl,k,rho,fmn,ier)
            zmns_sav(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm_nyq(mn))-2)
         END DO
         DEALLOCATE(ftemp,fmn)
         dr = 1./(ns-1)
         ! Calculate external boundary
         ALLOCATE(rmnc1(mnmax_nyq),zmns1(mnmax_nyq),rmnc2(mnmax_nyq),zmns2(mnmax_nyq))
         ALLOCATE(rmns1(mnmax_nyq),zmnc1(mnmax_nyq),rmns2(mnmax_nyq),zmnc2(mnmax_nyq))
         DO ik = 1, k
            IF (rho(ik) > 1.0) THEN
               rmnc_sav(:,ik) = 0.0
               zmns_sav(:,ik) = 0.0
               dr = rho(ik) - rho(ik-1)
               rmnc1(1:mnmax_nyq) = rmnc_sav(1:mnmax_nyq,ik-1)
               zmns1(1:mnmax_nyq) = zmns_sav(1:mnmax_nyq,ik-1)
               rmns1 = 0.0; zmnc1 = 0.0
               !IF (lasym) rmns1(1:mnmax_nyq) = rmns_sav(1:mnmax_nyq,ik-1)
               !IF (lasym) zmnc1(1:mnmax_nyq) = zmnc_sav(1:mnmax_nyq,ik-1)
               CALL scaleup_boundary(dr,mnmax_nyq,xm_temp,xn_temp,rmnc1,zmns1,rmnc2,zmns2,rmns1,zmnc1,rmns2,zmnc2)
               rmnc_sav(:,ik) = rmnc2
               zmns_sav(:,ik) = zmns2
               !IF (lasym) rmns_sav(:,ik) = rmns2
               !IF (lasym) zmnc_sav(:,ik) = zmnc2
               nexternal = nexternal - 1
            END IF
         END DO
         DEALLOCATE(rmnc1,zmns1,rmnc2,zmns2)
         DEALLOCATE(rmns1,zmnc1,rmns2,zmnc2)
         


         ! Transform to real space
         xn_temp = -xn_temp/nfp
         FORALL(u = 1:nu) xu(u) = DBLE(u-1)/DBLE(NU)
         FORALL(v = 1:nv) xv(v) = DBLE(v-1)/DBLE(NV)
         STOP "NEED TO FIX FIELDLINES_WRITE_EMC3 TO SUPPORT mntouv_local NOW PRIVATE"
         !CALL mntouv_local(1,k,mnmax_nyq,NU,NV,xu,xv,rmnc_sav,xm_temp,xn_temp,r,0,1)
         !CALL mntouv_local(1,k,mnmax_nyq,NU,NV,xu,xv,zmns_sav,xm_temp,xn_temp,z,1,0)
         ! Output to file
         WRITE(6,'(A)')  '   FILE: '//'emc3_equil_'//TRIM(id_string)//'.txt'
         iunit = 100
         CALL safe_open(iunit,ier,'emc3_equil_'//TRIM(id_string)//'.txt','replace','formatted')
         WRITE(iunit,'(4(2X,I4.4))') NU,NV,INT(k*.5),k
         DO ik = INT(k*.5), k
            WRITE(iunit,'(I4.4)') ik
            DO v = 1, NV
               DO u = 1, NU
                  WRITE(iunit,'(2(E23.16,2X))') r(u,v,ik),z(u,v,ik)
               END DO
            END DO
         END DO
         CLOSE(iunit)
         DEALLOCATE(xm_temp,xn_temp,rmnc_sav,zmns_sav)
         DEALLOCATE(xu,xv,r,z)
      END IF

      ! Write the Field information
      IF (lbfield_only) THEN
         CALL safe_open(iunit,ier,'layout','replace','formatted')
         WRITE(iunit,*) nr,nz,nphi,INT(pi2/phiaxis(nphi)),raxis(1),raxis(nr),zaxis(1),zaxis(nz)
         CLOSE(iunit)
         CALL safe_open(iunit,ier,'B_R','replace','formatted')
         WRITE(iunit,*) (((B_R(i,k,j), k=1,nphi-1), j=1,nz), i=1,nr)
         CLOSE(iunit)
         CALL safe_open(iunit,ier,'B_PHI','replace','formatted')
         WRITE(iunit,*) (((B_PHI(i,k,j), k=1,nphi-1), j=1,nz), i=1,nr)
         CLOSE(iunit)
         CALL safe_open(iunit,ier,'B_Z','replace','formatted')
         WRITE(iunit,*) (((B_Z(i,k,j), k=1,nphi-1), j=1,nz), i=1,nr)
         CLOSE(iunit)
      ELSE IF (lafield_only) THEN
         CALL safe_open(iunit,ier,'layout','replace','formatted')
         WRITE(iunit,*) nr,nz,nphi,INT(pi2/phiaxis(nphi)),raxis(1),raxis(nr),zaxis(1),zaxis(nz)
         CLOSE(iunit)
         CALL safe_open(iunit,ier,'A_R','replace','formatted')
         WRITE(iunit,*) (((B_R(i,k,j), k=1,nphi-1), j=1,nz), i=1,nr)
         CLOSE(iunit)
         CALL safe_open(iunit,ier,'A_PHI','replace','formatted')
         WRITE(iunit,*) (((B_PHI(i,k,j), k=1,nphi-1), j=1,nz), i=1,nr)
         CLOSE(iunit)
         CALL safe_open(iunit,ier,'A_Z','replace','formatted')
         WRITE(iunit,*) (((B_Z(i,k,j), k=1,nphi-1), j=1,nz), i=1,nr)
         CLOSE(iunit)
      END IF

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_write_emc3
