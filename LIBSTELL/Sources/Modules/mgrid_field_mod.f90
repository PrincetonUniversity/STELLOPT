!-----------------------------------------------------------------------
!     Module:    mgrid_field_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This module contains routines for working with
!                    mgrid files.  
!-----------------------------------------------------------------------
      MODULE mgrid_field_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE mgrid_mod, nextcur_mgrid => nextcur
!      USE mgrid_mod, ONLY: read_mgrid, np0b, nfper0, rminb, rmaxb,&
!                           zminb, zmaxb, nr0b, np0b, nz0b, bvac,&
!                           nextcur_mgrid => nextcur, free_mgrid, &
!                           curlabel, mgrid_path_old
      USE EZspline_obj
      USE EZspline
      
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER            :: nr_mgrid, nphi_mgrid, nz_mgrid
      REAL(rprec)        :: pi2, rmin_mgrid, rmax_mgrid, zmin_mgrid,&
                            zmax_mgrid, phimin_mgrid, phimax_mgrid
      CHARACTER(256)     :: mgrid_filename
      TYPE(EZspline3_r8) :: brm_spl, bzm_spl, bphim_spl
      REAL(rprec), DIMENSION(:),       POINTER, PRIVATE :: r_vac, phi_vac, z_vac
      REAL(rprec), DIMENSION(:,:,:,:), POINTER, PRIVATE :: BRV4D, BPV4D, BZV4D
      INTEGER, PRIVATE   :: win_BRV4D, win_BPV4D, win_BZV4D, eps1, eps2, eps3, &
                            win_rvac, win_phivac, win_zvac
      REAL*8, parameter, PRIVATE :: small = 1.e-10_ezspline_r8
!-----------------------------------------------------------------------
!     Subroutines
!         mgrid_load:   Reads an mgrid file and loads a splines
!         mgrid_b:      Calculates the field at a point in space
!         mgrid_free:   Free's allocated variables
!         mgrid_info:   Displays information about the mgrid file
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE mgrid_load(filename,extcur,ncur,nv,nfp,istat,ithread,comm)
      USE mpi_sharmem
      USE mpi_params 
      USE mpi_inc
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(inout)        :: ncur, nv, nfp
      REAL(rprec), INTENT(in)       :: extcur(ncur)
      INTEGER, INTENT(inout)        :: istat
      INTEGER, INTENT(in)           :: ithread
      INTEGER, INTENT(in), OPTIONAL :: comm
      INTEGER :: i,j,v,ik, nv_temp, nfp_temp
      INTEGER :: shar_comm, shar_rank, shar_size, local_master, mylocalid
      INTEGER :: bcs1(2), bcs2(2), bcs3(2)
      REAL(rprec), ALLOCATABLE :: br_vac(:,:,:), bphi_vac(:,:,:), bz_vac(:,:,:)
      LOGICAL, PARAMETER :: lbug=.false.
      !pi2 = 2.0*3.14159265358979323846
      pi2 = 8 * ATAN(1._rprec)
      istat = 0
      local_master = 0
      mgrid_filename=TRIM(filename)
      nv_temp = nv
      nfp_temp = nfp
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
         CALL MPI_COMM_RANK(shar_comm, shar_rank, istat)
         CALL MPI_COMM_SIZE(shar_comm, shar_size, istat)
      ELSE
#endif
         shar_rank = 0; shar_size = 1
#if defined(MPI_OPT)
      END IF
#endif

      IF (.not.ASSOCIATED(bvac)) THEN
#if defined(MPI_OPT)
         IF (PRESENT(comm)) THEN
            CALL read_mgrid(TRIM(mgrid_filename),extcur,nv_temp,nfp_temp,lbug,istat,comm)
         ELSE
#endif
            CALL read_mgrid(TRIM(mgrid_filename),extcur,nv_temp,nfp_temp,lbug,istat)
#if defined(MPI_OPT)
         END IF
#endif
         nv_temp = np0b
         nfp_temp = nfper0
         IF (istat == 9) THEN
            IF (ASSOCIATED(curlabel)) DEALLOCATE(curlabel)
            mgrid_path_old = " "
            nv_temp  = np0b
            nfp_temp = nfper0
            istat = 0
#if defined(MPI_OPT)
            IF (PRESENT(comm)) THEN
               CALL read_mgrid(TRIM(mgrid_filename),extcur,nv_temp,nfp_temp,lbug,istat,comm)
            ELSE
#endif
               CALL read_mgrid(TRIM(mgrid_filename),extcur,nv_temp,nfp_temp,lbug,istat)
#if defined(MPI_OPT)
            END IF
#endif
            IF (istat .ne. 0) stop 'ERROR Reading mgrid_file'
         END IF
      ELSE
         nv_temp=np0b
         nfp_temp=nfper0
      END IF
      nv  = nv_temp
      nfp = nfp_temp
      rmin_mgrid   = rminb
      rmax_mgrid   = rmaxb
      zmin_mgrid   = zminb
      zmax_mgrid   = zmaxb
      phimin_mgrid = 0.0
      phimax_mgrid = pi2/nfp_temp
      nr_mgrid     = nr0b
      nphi_mgrid   = np0b+1
      nz_mgrid     = nz0b
      eps1 = (rmax_mgrid-rmin_mgrid)*small
      eps2 = (phimax_mgrid-phimin_mgrid)*small
      eps3 = (zmax_mgrid-zmin_mgrid)*small
      ! Recompose vacuum field for splines
      IF (EZspline_allocated(brm_spl)) CALL EZspline_free(brm_spl,istat)
      IF (EZspline_allocated(bphim_spl)) CALL EZspline_free(bphim_spl,istat)
      IF (EZspline_allocated(bzm_spl)) CALL EZspline_free(bzm_spl,istat)
!      IF (np0b == 1) THEN
!         ALLOCATE(br_vac(1:nr0b,1:2,1:nz0b),bphi_vac(1:nr0b,1:2,1:nz0b),&
!                  bz_vac(1:nr0b,1:2,1:nz0b),STAT=istat)
!         ALLOCATE(r_vac(1:nr0b),phi_vac(1:2),z_vac(1:nz0b),STAT=istat)
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPIALLOC(BRV4D,8,nr0b,nphi_mgrid,nz0b,shar_rank,0,shar_comm,win_BRV4D)
         CALL MPIALLOC(BPV4D,8,nr0b,nphi_mgrid,nz0b,shar_rank,0,shar_comm,win_BPV4D)
         CALL MPIALLOC(BZV4D,8,nr0b,nphi_mgrid,nz0b,shar_rank,0,shar_comm,win_BZV4D)
         CALL MPIALLOC(r_vac,nr0b,shar_rank,0,shar_comm,win_rvac)
         CALL MPIALLOC(phi_vac,nphi_mgrid,shar_rank,0,shar_comm,win_phivac)
         CALL MPIALLOC(z_vac,nz0b,shar_rank,0,shar_comm,win_zvac)
      ELSE
#endif
         ALLOCATE(BRV4D(8,nr0b,nphi_mgrid,nz0b), STAT=istat)
         ALLOCATE(BPV4D(8,nr0b,nphi_mgrid,nz0b), STAT=istat)
         ALLOCATE(BZV4D(8,nr0b,nphi_mgrid,nz0b), STAT=istat)
         ALLOCATE(r_vac(1:nr0b),phi_vac(1:nphi_mgrid),z_vac(1:nz0b),STAT=istat)
#if defined(MPI_OPT)
      END IF
#endif
      IF (shar_rank == 0) THEN
         ALLOCATE(br_vac(1:nr0b,1:nphi_mgrid,1:nz0b),bphi_vac(1:nr0b,1:nphi_mgrid,1:nz0b),&
                  bz_vac(1:nr0b,1:nphi_mgrid,1:nz0b),STAT=istat)
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
         FORALL(i = 1:nphi_mgrid) phi_vac(i) = (i-1)*(pi2/nfp)/(nphi_mgrid-1)
         br_vac(:,nphi_mgrid,:)   = br_vac(:,1,:)
         bphi_vac(:,nphi_mgrid,:) = bphi_vac(:,1,:)
         bz_vac(:,nphi_mgrid,:)   = bz_vac(:,1,:)
         bcs1=(/ 0, 0/)
         bcs2=(/-1,-1/)
         bcs3=(/ 0, 0/)
         CALL EZspline_init(brm_spl,nr_mgrid,nphi_mgrid,nz_mgrid,bcs1,bcs2,bcs3,istat)
         brm_spl%isHermite = 1
         brm_spl%x1 = r_vac
         brm_spl%x2 = phi_vac
         brm_spl%x3 = z_vac
         CALL EZspline_setup(brm_spl,br_vac,istat,EXACT_DIM=.true.)
         CALL EZspline_init(bphim_spl,nr_mgrid,nphi_mgrid,nz_mgrid,bcs1,bcs2,bcs3,istat)
         bphim_spl%isHermite = 1
         bphim_spl%x1 = r_vac
         bphim_spl%x2 = phi_vac
         bphim_spl%x3 = z_vac
         CALL EZspline_setup(bphim_spl,bphi_vac,istat,EXACT_DIM=.true.)
         CALL EZspline_init(bzm_spl,nr_mgrid,nphi_mgrid,nz_mgrid,bcs1,bcs2,bcs3,istat)
         bzm_spl%isHermite = 1
         bzm_spl%x1 = r_vac
         bzm_spl%x2 = phi_vac
         bzm_spl%x3 = z_vac
         CALL EZspline_setup(bzm_spl,bz_vac,istat,EXACT_DIM=.true.)
         DEALLOCATE(br_vac,bphi_vac,bz_vac)
         DEALLOCATE(r_vac,phi_vac,z_vac)
         BRV4D = brm_spl%fspl
         BPV4D = bphim_spl%fspl
         BZV4D = bzm_spl%fspl
         CALL EZspline_free(brm_spl,istat)
         CALL EZspline_free(bphim_spl,istat)
         CALL EZspline_free(bzm_spl,istat)
      END IF
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_BARRIER(comm,istat)
         CALL MPI_COMM_FREE(shar_comm,istat) 
      END IF
#endif
      phimax_mgrid = MAXVAL(phi_vac)
      phimin_mgrid = MINVAL(phi_vac) 
      RETURN
      END SUBROUTINE mgrid_load
      
      SUBROUTINE mgrid_bcart(xp,yp,zp,bx,by,bz,istat)
      IMPLICIT NONE
      REAL(rprec), INTENT(in)  :: xp, yp, zp
      REAL(rprec), INTENT(out) :: bx, by, bz
      INTEGER, INTENT(inout)   :: istat
      REAL(rprec) :: rp,phip,br, bphi
      rp = SQRT(xp*xp+yp*yp)
      phip = ATAN2(yp,xp)
      IF (phip < 0.0) phip = phip + pi2
      CALL mgrid_bcyl(rp, phip, zp, br, bphi, bz, istat)
      bx = br * cos(phip) - bphi * sin(phip)
      by = br * sin(phip) + bphi * cos(phip)
      RETURN
      END SUBROUTINE mgrid_bcart
      
      SUBROUTINE mgrid_bcyl(rp,phip,zp,br,bphi,bz,istat)
      IMPLICIT NONE
      REAL(rprec), INTENT(in)  :: rp, phip, zp
      REAL(rprec), INTENT(out) :: br, bphi, bz
      INTEGER, INTENT(inout)   :: istat
      REAL(rprec) :: phi2
      ! For splines
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(1)
      INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
      REAL*8, PARAMETER :: one = 1
      phi2 = phip
      IF (phi2 > phimax_mgrid) phi2 = MOD(phi2,phimax_mgrid)
      !CALL EZspline_isInDomain(brm_spl,rp,phi2,zp,istat)
      IF ((rp >= rmin_mgrid-eps1) .and. (rp <= rmax_mgrid+eps1) .and. &
          (phi2 >= phimin_mgrid-eps2) .and. (phi2 <= phimax_mgrid+eps2) .and. &
          (zp >= zmin_mgrid-eps3) .and. (zp <= zmax_mgrid+eps3)) THEN
      !IF (istat == 0) THEN
         i = MIN(MAX(COUNT(r_vac < rp),1),nr_mgrid-1)
         j = MIN(MAX(COUNT(phi_vac < phi2),1),nphi_mgrid-1)
         k = MIN(MAX(COUNT(z_vac < zp),1),nz_mgrid-1)
         hx     = r_vac(i+1) - r_vac(i)
         hy     = phi_vac(j+1) - phi_vac(j)
         hz     = z_vac(k+1) - z_vac(k)
         hxi    = one / hx
         hyi    = one / hy
         hzi    = one / hz
         xparam = (rp - r_vac(i)) * hxi
         yparam = (phi2 - phi_vac(j)) * hyi
         zparam = (zp - z_vac(k)) * hzi
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BRV4D(1,1,1,1),nr_mgrid,nphi_mgrid,nz_mgrid)
         br = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BPV4D(1,1,1,1),nr_mgrid,nphi_mgrid,nz_mgrid)
         bphi = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BZV4D(1,1,1,1),nr_mgrid,nphi_mgrid,nz_mgrid)
         bz = fval(1)
         !CALL EZspline_interp(brm_spl, rp, phi2, zp, br, istat)
         !CALL EZspline_interp(bphim_spl, rp, phi2, zp, bphi, istat)
         !CALL EZspline_interp(bzm_spl, rp, phi2, zp, bz, istat)
      END IF
      RETURN
      END SUBROUTINE mgrid_bcyl
      
      SUBROUTINE mgrid_free(istat,comm)
      USE mpi_sharmem
      IMPLICIT NONE
      INTEGER, INTENT(out) :: istat
      INTEGER, INTENT(in), OPTIONAL :: comm
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL mpidealloc(BRV4D,win_BRV4D)
         CALL mpidealloc(BPV4D,win_BPV4D)
         CALL mpidealloc(BZV4D,win_BZV4D)
         CALL mpidealloc(r_vac,win_rvac)
         CALL mpidealloc(phi_vac,win_phivac)
         CALL mpidealloc(z_vac,win_zvac)
         CALL free_mgrid(istat,comm)
      ELSE
#endif
         DEALLOCATE(BRV4D, BPV4D, BZV4D, stat=istat)
         DEALLOCATE(r_vac, phi_vac, z_vac, stat=istat)
         CALL free_mgrid(istat)
#if defined(MPI_OPT)
      END IF
#endif
      IF (ASSOCIATED(curlabel)) DEALLOCATE(curlabel)
      mgrid_path_old = " "
      END SUBROUTINE
      
      SUBROUTINE mgrid_info(iunit)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: iunit
      WRITE(iunit,'(A)')   '----- MGRID Information -----'
      WRITE(iunit,'(A,A)') '   FILE:',TRIM(mgrid_filename)
      WRITE(iunit,'(A,F8.5,A,F8.5,A,I4)')   '   R   = [',rmin_mgrid,',',rmax_mgrid,'];  NR   = ',nr_mgrid
      WRITE(iunit,'(A,F8.5,A,F8.5,A,I4)')   '   PHI = [',phimin_mgrid,',',phimax_mgrid,'];  NPHI = ',nphi_mgrid
      WRITE(iunit,'(A,F8.5,A,F8.5,A,I4)')   '   Z   = [',zmin_mgrid,',',zmax_mgrid,'];  NZ   = ',nz_mgrid
      CALL FLUSH(iunit)
      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE mgrid_field_mod
