!-----------------------------------------------------------------------
!     Module:        read_nescoil_mod
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          06/07/2024
!     Description:   This module stores routines for reading
!                    NESCOIL data.
!                    pot(u,v) = sum_mn Phi_mn sin(mu+nv) 
!                               - Ip*v/nfp - It*u
!                    Eq 4. Merkel et al. NF 27, 5 (1987)
!-----------------------------------------------------------------------
      MODULE read_nescoil_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE safe_open_mod
      USE EZspline_obj
      USE EZspline

!-----------------------------------------------------------------------
!     Module Variables
!           NTITLE:   
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL :: lasym
      INTEGER :: nu, nv, nu1, nv1, mpol, ntor, mf, nf, md, nd, np, &
                 ibex, mstrt, mstep, mkeep, mdspw, w_psurf, w_csurf, &
                 w_bnuv, w_jsurf, w_xerr, w_svd, mnmax_plasma, &
                 mnmax_surface, nmax, mnd, nuv, nuv1, nuvh, nuvh1, &
                 nvp, mnmax_pot
      REAL(rprec) :: iota_edge, phip_edge, curpol, cut, cup, curwt, &
                     trgwt, fnuv, alp
      INTEGER, DIMENSION(:), ALLOCATABLE ::xm_plasma, xn_plasma,      & 
                                          xm_surface, xn_surface,     &
                                          xm_pot, xn_pot 
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::   &
                                          rmnc_plasma, zmns_plasma,   &
                                          rmns_plasma, zmnc_plasma,   &
                                          lmnc_plasma, lmns_plasma,       &
                                          rmnc_surface, zmns_surface, &
                                          rmns_surface, zmnc_surface, &
                                          potmns_surface
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: x_plasma, y_plasma,   &
                                                z_plasma, r_plasma,   &
                                          nx_plasma, ny_plasma,       &
                                          nz_plasma,   dsur_plasma,   &
                                          dxdu_plasma, dydu_plasma,   &
                                          dxdv_plasma, dydv_plasma,   &
                                          bn_plasma,                  & 
                                          x_surface, y_surface,       &
                                          z_surface, r_surface,       &
                                          nx_surface, ny_surface,     &
                                          nz_surface, dsur_surface,   &
                                          xcur_surface, ycur_surface, &
                                          zcur_surface,               &
                                          db_normal, Babs
      ! FILE ACCESS
      CHARACTER(LEN=256), PRIVATE :: line
      ! For the info routine
      REAL(rprec), PRIVATE :: surf_area
      ! SHARED MEMORY
      INTEGER, PRIVATE :: win_X3D, win_Y3D, win_Z3D, win_KX3D, win_KY3D, win_KZ3D, win_x1, win_x2
      DOUBLE PRECISION, DIMENSION(:), POINTER, PRIVATE :: x1, x2
      DOUBLE PRECISION, DIMENSION(:,:,:),   POINTER, PRIVATE :: X3D, Y3D, Z3D, KX3D, KY3D, KZ3D
      ! Integration Related
      INTEGER, PRIVATE :: nu_int, nx1, nx2, u1,v1
      REAL(rprec), PRIVATE :: norm
      REAL*8, PRIVATE :: eps1, eps2, x1_min, x1_max, x2_min, x2_max
      DOUBLE PRECISION, PRIVATE :: x_f, y_f, z_f, norm_fsub


      ! CONSTANTS
      REAL*8, parameter, PRIVATE :: small = 1.e-10_ezspline_r8
      DOUBLE PRECISION, PARAMETER, PRIVATE :: pi2 = 6.283185482025146D+00
      DOUBLE PRECISION, PARAMETER, PRIVATE :: zero = 0.0D+00
      DOUBLE PRECISION, PARAMETER, PRIVATE :: one  = 1.0D+00

!-----------------------------------------------------------------------
!     Private Subroutines
!-----------------------------------------------------------------------
      PRIVATE :: mntouv_local, lookupgrid2d, funsub_b, evalbi2D, &
                 read_nescout_plasmaboundary, read_nescout_coilsurface, &
                 read_nescout_plasmaxyz, read_nescout_plasmanormal, &
                 read_nescout_plasmaderiv, read_nescout_surfacexyz, &
                 read_nescout_surfacecurrent, read_nescout_accuracy

!-----------------------------------------------------------------------
!     Module Subroutines
!          nescoil_bfield        B-Field using direct integration
!          nescoil_bfield_adapt  B-Field using adaptive integration
!-----------------------------------------------------------------------
      INTERFACE nescoil_bfield
         MODULE PROCEDURE nescoil_bfield_dbl, nescoil_bfield_flt
      END INTERFACE

      INTERFACE nescoil_bfield_adapt
         MODULE PROCEDURE nescoil_bfield_adapt_dbl, nescoil_bfield_adapt_flt
      END INTERFACE

      CONTAINS

!-----------------------------------------------------------------------
!     Module Subroutines
!           read_nescout:  For reading NESCOIL output files
!           
!-----------------------------------------------------------------------
      SUBROUTINE read_nescout(filename,istat)
         IMPLICIT NONE
         CHARACTER(*), INTENT(in) :: filename
         INTEGER, INTENT(inout) :: istat
         LOGICAL :: lexist
         INTEGER :: iunit, n, nlines, mtemp,ntemp

         ! Check for the file
         INQUIRE(FILE=TRIM(filename),EXIST=lexist)
         IF (.not.lexist) THEN
            istat=-1
            WRITE(6,*) " ERROR: Could not find file: "//TRIM(filename)
            RETURN
         END IF

         ! Open the file
         CALL safe_open (iunit, istat, filename, 'old', 'formatted')
         IF (istat /= 0) THEN
            WRITE(6,*) " ERROR: Could not open file: "//TRIM(filename)
            RETURN
         END IF

         CALL read_nescout_deallocate

         ! Read through file
         DO
            READ(iunit, '(A)', iostat=istat) line
            IF (istat /=0) EXIT
            ! Handle each option separately
            SELECT CASE (TRIM(line))
               CASE("----- Grid Spatial Dimensions -----")
                  READ(iunit, '(A)', iostat=istat) line
                  READ(iunit, '(6i6,l)') nu, nv, nu1, nv1, mpol, ntor, lasym
               CASE("----- Fourier Dimensions -----")
                  READ(iunit, '(A)', iostat=istat) line
                  READ(iunit, '(4i6)') mf, nf, md, nd
                  nmax  = 3*nv*(nu+2)+1
                  mnd   = (md + 1)*(2*nd + 1)
                  nuv   = nu*nv
                  nuv1  = nu1*nv1
                  nuvh  = nuv/2 + nu
                  nuvh1 = nuv1/2 + nu1
                  fnuv  = 1.0 / nuv
               CASE("----- Plasma information from VMEC -----")
                  READ(iunit, '(A)', iostat=istat) line
                  READ(iunit, '(i6,3g25.16)') np, iota_edge, phip_edge, curpol
                  alp   = pi2/np
               CASE("----- Current Controls -----")
                  READ(iunit, '(A)', iostat=istat) line
                  READ(iunit, '(2g25.16,i6)') cut, cup, ibex
               CASE("----- SVD controls -----")
                  READ(iunit, '(A)', iostat=istat) line
                  READ(iunit, '(4i6,2g25.16)') mstrt, mstep, mkeep, mdspw, curwt, trgwt
               CASE("----- Output controls -----")
                  READ(iunit, '(A)', iostat=istat) line
                  READ(iunit, '(6i6)') w_psurf, w_csurf, w_bnuv, w_jsurf, w_xerr, w_svd
               CASE("----- Plasma Surface -----")
                  READ(iunit, '(A)', iostat=istat) line
                  READ(iunit, *) mnmax_plasma
               CASE("----- Plasma boundary fourier coefficients  -----")
                  READ(iunit, '(A)', iostat=istat) line
                  CALL read_nescout_plasmaboundary(iunit,istat)
               CASE("----- Coil Surface -----")
                  READ(iunit, '(A)', iostat=istat) line
                  READ(iunit, *) mnmax_surface
               CASE("----- Coil surface fourier coefficients -----")
                  READ(iunit, '(A)', iostat=istat) line
                  CALL read_nescout_coilsurface(iunit,istat)
               CASE("----- Calling Surface_Plasma -----")
                  IF (w_psurf>0 .or. w_psurf==-1) THEN
                     n = MAXVAL(ABS(xn_plasma))
                     nlines = 6 + 2*n+1 + 5 +2*n+1
                     DO n = 1, nlines
                        READ(iunit, '(A)', iostat=istat) line
                     END DO
                  END IF
                  IF (w_psurf>1 .or. w_psurf==-2) CALL read_nescout_plasmaxyz(iunit,istat)
                  IF (w_psurf>2 .or. w_psurf==-3) CALL read_nescout_plasmanormal(iunit,istat)
                  IF (w_psurf>2 .or. w_psurf==-3) CALL read_nescout_plasmaderiv(iunit,istat)
               CASE("Calculated bn_ext(u,v) on plasma surface")
                  IF ( w_bnuv > 1 .or. w_bnuv == -2 ) THEN
                     READ(iunit, '(A)', iostat=istat) line
                     ALLOCATE(bn_plasma(nuvh1))
                     READ(iunit,*) bn_plasma
                  END IF
               CASE("----- Calling Surface_Coil -----")
                  IF (w_csurf>0 .or. w_csurf==-1) THEN
                     n = MAXVAL(ABS(xn_surface))
                     nlines = 5 + 2*n+1 + 5 +2*n+1
                     DO n = 1, nlines
                        READ(iunit, '(A)', iostat=istat) line
                     END DO
                  END IF
                  IF (w_csurf>1 .or. w_csurf==-2) CALL read_nescout_surfacexyz(iunit,istat)
                  IF (w_csurf>2 .or. w_csurf==-3) CALL read_nescout_surfacenormal(iunit,istat)
                  IF (w_csurf>3 .or. w_csurf==-4) CALL read_nescout_surfacecurrent(iunit,istat)
               CASE("---- Phi(m,n) for least squares ---")
                  mnmax_pot = (mf+1)*(2*nf+1)
                  ALLOCATE(xm_pot(mnmax_pot),xn_pot(mnmax_pot))
                  ALLOCATE(potmns_surface(mnmax_pot))
                  DO n = 1, mnmax_pot
                     READ(iunit,'(i3,2x,i3,2x,g25.16)') xm_pot(n),xn_pot(n),potmns_surface(n)
                  END DO
               CASE("----- Calling Surfcur_Diag -----")
                  READ(iunit, '(A)', iostat=istat) line
                  READ(iunit, '(A)', iostat=istat) line
                  IF (w_jsurf>0 .or. w_jsurf==-1) THEN
                     READ(iunit, '(A)', iostat=istat) line
                  END IF
                  IF (w_jsurf>1 .or. w_jsurf==-2) THEN
                     READ(iunit, '(A)', iostat=istat) line
                  END IF
               CASE("----- Calling Accuracy -----")
                  IF (w_bnuv>0 .or. w_bnuv==-1) CALL read_nescout_accuracy(iunit,istat)


            END SELECT
         END DO
         istat = 0

         ! Close file
         CLOSE(iunit)

         RETURN
      END SUBROUTINE read_nescout

      SUBROUTINE read_nescout_plasmaboundary(iunit,istat)
         IMPLICIT NONE
         INTEGER, INTENT(inout) :: iunit
         INTEGER, INTENT(inout) :: istat
         INTEGER :: mn
         ALLOCATE(xm_plasma(mnmax_plasma),xn_plasma(mnmax_plasma))
         ALLOCATE(rmnc_plasma(mnmax_plasma),zmns_plasma(mnmax_plasma))
         ALLOCATE(rmns_plasma(mnmax_plasma),zmnc_plasma(mnmax_plasma))
         ALLOCATE(lmnc_plasma(mnmax_plasma),lmns_plasma(mnmax_plasma))
         DO mn = 1, mnmax_plasma
            READ(iunit,'(2i4,6g20.10)') xm_plasma(mn), xn_plasma(mn), &
                  rmnc_plasma(mn), zmns_plasma(mn), lmns_plasma(mn),  &
                  rmns_plasma(mn), zmnc_plasma(mn), lmnc_plasma(mn)
         END DO
         
         RETURN
      END SUBROUTINE read_nescout_plasmaboundary

      SUBROUTINE read_nescout_coilsurface(iunit,istat)
         IMPLICIT NONE
         INTEGER, INTENT(inout) :: iunit
         INTEGER, INTENT(inout) :: istat
         INTEGER :: mn
         ALLOCATE(xm_surface(mnmax_surface),xn_surface(mnmax_surface))
         ALLOCATE(rmnc_surface(mnmax_surface),zmns_surface(mnmax_surface))
         ALLOCATE(rmns_surface(mnmax_surface),zmnc_surface(mnmax_surface))
         DO mn = 1, mnmax_surface
            READ(iunit,'(2i4,4g20.10)') xm_surface(mn), xn_surface(mn), &
                                    rmnc_surface(mn), zmns_surface(mn), &
                                    rmns_surface(mn), zmnc_surface(mn)
         END DO
         RETURN
      END SUBROUTINE read_nescout_coilsurface

      SUBROUTINE read_nescout_plasmaxyz(iunit,istat)
         IMPLICIT NONE
         INTEGER, INTENT(inout) :: iunit
         INTEGER, INTENT(inout) :: istat
         INTEGER :: mn
         READ(iunit, '(A)', iostat=istat) line
         ALLOCATE(x_plasma(nuvh1),y_plasma(nuvh1))
         ALLOCATE(z_plasma(nuvh1),r_plasma(nuvh1))
         DO mn = 1, nuvh1
            READ(iunit,'(4g16.6)') x_plasma(mn), y_plasma(mn), &
                  z_plasma(mn), r_plasma(mn)
         END DO
         READ(iunit, '(A)', iostat=istat) line
         RETURN
      END SUBROUTINE read_nescout_plasmaxyz

      SUBROUTINE read_nescout_plasmanormal(iunit,istat)
         IMPLICIT NONE
         INTEGER, INTENT(inout) :: iunit
         INTEGER, INTENT(inout) :: istat
         INTEGER :: mn
         READ(iunit, '(A)', iostat=istat) line
         ALLOCATE(nx_plasma(nuvh1),ny_plasma(nuvh1))
         ALLOCATE(nz_plasma(nuvh1),dsur_plasma(nuvh1))
         DO mn = 1, nuvh1
            READ(iunit,'(4g16.6)') dsur_plasma(mn), nx_plasma(mn), &
                                   ny_plasma(mn), nz_plasma(mn)
         END DO
         READ(iunit, '(A)', iostat=istat) line
         RETURN
      END SUBROUTINE read_nescout_plasmanormal

      SUBROUTINE read_nescout_plasmaderiv(iunit,istat)
         IMPLICIT NONE
         INTEGER, INTENT(inout) :: iunit
         INTEGER, INTENT(inout) :: istat
         INTEGER :: mn
         READ(iunit, '(A)', iostat=istat) line
         ALLOCATE(dxdu_plasma(nuvh1),dydu_plasma(nuvh1))
         ALLOCATE(dxdv_plasma(nuvh1),dydv_plasma(nuvh1))
         DO mn = 1, nuvh1
            READ(iunit,'(4g16.6)') dxdu_plasma(mn), dydu_plasma(mn), &
                                   dxdv_plasma(mn), dydv_plasma(mn)
         END DO
         READ(iunit, '(A)', iostat=istat) line
         RETURN
      END SUBROUTINE read_nescout_plasmaderiv

      SUBROUTINE read_nescout_surfacexyz(iunit,istat)
         IMPLICIT NONE
         INTEGER, INTENT(inout) :: iunit
         INTEGER, INTENT(inout) :: istat
         INTEGER :: mn
         READ(iunit, '(A)', iostat=istat) line
         ALLOCATE(x_surface(nuvh),y_surface(nuvh))
         ALLOCATE(z_surface(nuvh),r_surface(nuvh))
         DO mn = 1, nuvh
            READ(iunit,'(4g16.6)') x_surface(mn), y_surface(mn), &
                  z_surface(mn), r_surface(mn)
         END DO
         READ(iunit, '(A)', iostat=istat) line
         RETURN
      END SUBROUTINE read_nescout_surfacexyz

      SUBROUTINE read_nescout_surfacenormal(iunit,istat)
         IMPLICIT NONE
         INTEGER, INTENT(inout) :: iunit
         INTEGER, INTENT(inout) :: istat
         INTEGER :: mn
         READ(iunit, '(A)', iostat=istat) line
         ALLOCATE(nx_surface(nuvh),ny_surface(nuvh))
         ALLOCATE(nz_surface(nuvh),dsur_surface(nuvh))
         DO mn = 1, nuvh
            READ(iunit,'(4g16.6)') dsur_surface(mn), nx_surface(mn), &
                                   ny_surface(mn), nz_surface(mn)
         END DO
         READ(iunit, '(A)', iostat=istat) line
         RETURN
      END SUBROUTINE read_nescout_surfacenormal

      SUBROUTINE read_nescout_surfacecurrent(iunit,istat)
         IMPLICIT NONE
         INTEGER, INTENT(inout) :: iunit
         INTEGER, INTENT(inout) :: istat
         INTEGER :: mn
         READ(iunit, '(A)', iostat=istat) line
         ALLOCATE(xcur_surface(nuvh),ycur_surface(nuvh))
         ALLOCATE(zcur_surface(nuvh))
         DO mn = 1, nuvh
            READ(iunit,'(3g16.6)') xcur_surface(mn), ycur_surface(mn), &
                                   zcur_surface(mn)
         END DO
         READ(iunit, '(A)', iostat=istat) line
         RETURN
      END SUBROUTINE read_nescout_surfacecurrent

      SUBROUTINE read_nescout_accuracy(iunit,istat)
         IMPLICIT NONE
         INTEGER, INTENT(inout) :: iunit
         INTEGER, INTENT(inout) :: istat
         INTEGER :: mn
         REAL(rprec) :: dsur_local
         READ(iunit, '(A)', iostat=istat) line
         ALLOCATE(db_normal(nuvh1), Babs(nuvh1))
         DO mn = 1, nuvh1
            READ(iunit,*) db_normal(mn), Babs(mn), &
                                   dsur_local
         END DO
         RETURN
      END SUBROUTINE read_nescout_accuracy

      SUBROUTINE nescoil_bfield_init_ctypes
         IMPLICIT NONE
         CALL nescoil_bfield_init(64,64)
      END SUBROUTINE nescoil_bfield_init_ctypes

      SUBROUTINE nescoil_info(iunit)
         IMPLICIT NONE
         ! INPUT VARIABLES
         INTEGER, INTENT(in)  :: iunit
         ! LOCAL VARIABLES
         ! BEGIN SUBROUTINE
         WRITE(iunit,'(A)')                    '----- NESCOIL Current Surface -----'
         WRITE(iunit,'(A,ES11.4,A)')           '   Surface Area: ',surf_area,' [m]'
         WRITE(iunit,'(A,ES11.4,A)')           '   Poloidal Current: ',curpol*np,' [A]'
         !WRITE(iunit,'(A,I4,A,I4,A,I4,A,I3)')   '   NR = ',nr_vc,';   NU = ',nu_vc,';  NV = ',nv_vc,';  NFP = ',nvp/nv_vc
         !WRITE(iunit,'(A,I6)')                  '   NUVP = ',nuvp
         !WRITE(iunit,'(A,I6,A,I8,A)')                  '   MIN_CLS = ',MIN_CLS,'   (',IWRK,')'
         CALL FLUSH(iunit)
      END SUBROUTINE nescoil_info

      SUBROUTINE nescoil_bfield_init(nu_local, nv_local, comm)
         USE mpi_sharmem
         USE mpi_inc
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: nu_local
         INTEGER, INTENT(IN) :: nv_local
         INTEGER, INTENT(in), OPTIONAL :: comm
         INTEGER :: mn, shar_rank, shar_size,  shar_comm, nu1, u, v, ier, i1, i2
         DOUBLE PRECISION, ALLOCATABLE :: xu(:), xv(:),           &
               fmn_temp(:), yu(:), yv(:), cop(:), sip(:)
         DOUBLE PRECISION, ALLOCATABLE :: rreal(:,:), zreal(:,:), &
               xreal(:,:), yreal(:,:), rureal(:,:), rvreal(:,:),  &
               zureal(:,:), zvreal(:,:), sxreal(:,:),             &
               syreal(:,:), szreal(:,:), potu(:,:), potv(:,:),    &
               potr(:,:), potp(:,:), potz(:,:), potx(:,:),        &
               poty(:,:)
         TYPE(EZspline2_r8)   :: f_spl
         INTEGER, PARAMETER, DIMENSION(2) :: bcs1 = (/-1,-1/)
         INTEGER, PARAMETER, DIMENSION(2) :: bcs2 = (/-1,-1/)

#if defined(MPI_OPT)
         IF (PRESENT(comm)) CALL MPI_BCAST(np,1,MPI_INTEGER,0,comm,ier)
#endif
         ! Factors
         nu_int = nu_local
         nvp    = (nv_local-1)*np + 1 ! Because of stellarator symmetry
         !norm   = DBLE(np) / DBLE(pi2*pi2*nu_local*nvp/2)
         u1 = nu_local-1
         v1 = nvp - 1
         !norm = 1E-7/(u1*v1)
         norm   = DBLE(np) / DBLE(pi2*2*u1*v1)
         norm_fsub = DBLE(np) / (2*pi2)
         ! These must be consistent with splines below
         nx1    = nu_int;  nx2    = nvp
         x1_min = 0; x2_min = 0
         x1_max = 1; x2_max = 1
         eps1 = (x1_max-x1_min)*small
         eps2 = (x2_max-x2_min)*small
         ! Allocate the background info
#if defined(MPI_OPT)
         IF (PRESENT(comm)) THEN
            ! Get rank
            ier = 0
            CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, ier)
            CALL MPI_COMM_RANK(shar_comm, shar_rank, ier)
            CALL MPI_COMM_SIZE(shar_comm, shar_size, ier)
            ! Deallocate
            IF (ASSOCIATED(x1))   CALL mpidealloc(x1,     win_x1)
            IF (ASSOCIATED(x2))   CALL mpidealloc(x1,     win_x2)
            IF (ASSOCIATED(X3D))  CALL mpidealloc(X3D,    win_X3D)
            IF (ASSOCIATED(Y3D))  CALL mpidealloc(Y3D,    win_Y3D)
            IF (ASSOCIATED(Z3D))  CALL mpidealloc(Z3D,    win_Z3D)
            IF (ASSOCIATED(KX3D)) CALL mpidealloc(X3D,    win_KX3D)
            IF (ASSOCIATED(KY3D)) CALL mpidealloc(Y3D,    win_KY3D)
            IF (ASSOCIATED(KZ3D)) CALL mpidealloc(Z3D,    win_KZ3D)
            ! Allocate
            CALL mpialloc(x1, nx1,                 shar_rank, 0, shar_comm, win_x1)
            CALL mpialloc(x2, nx2,                 shar_rank, 0, shar_comm, win_x2)
            CALL mpialloc(X3D,  4,  nu_local, nvp, shar_rank, 0, shar_comm, win_X3D)
            CALL mpialloc(Y3D,  4,  nu_local, nvp, shar_rank, 0, shar_comm, win_Y3D)
            CALL mpialloc(Z3D,  4,  nu_local, nvp, shar_rank, 0, shar_comm, win_Z3D)
            CALL mpialloc(KX3D, 4,  nu_local, nvp, shar_rank, 0, shar_comm, win_KX3D)
            CALL mpialloc(KY3D, 4,  nu_local, nvp, shar_rank, 0, shar_comm, win_KY3D)
            CALL mpialloc(KZ3D, 4,  nu_local, nvp, shar_rank, 0, shar_comm, win_KZ3D)
         ELSE
#endif
            shar_rank = 0; shar_size = 1
            IF (ASSOCIATED(x1))   DEALLOCATE(x1)
            IF (ASSOCIATED(x2))   DEALLOCATE(x2)
            IF (ASSOCIATED(X3D))  DEALLOCATE(X3D)
            IF (ASSOCIATED(Y3D))  DEALLOCATE(Y3D)
            IF (ASSOCIATED(Z3D))  DEALLOCATE(Z3D)
            IF (ASSOCIATED(KX3D))  DEALLOCATE(KX3D)
            IF (ASSOCIATED(KY3D))  DEALLOCATE(KY3D)
            IF (ASSOCIATED(KZ3D))  DEALLOCATE(KZ3D)
            ALLOCATE(x1(nx1))
            ALLOCATE(x2(nx2))
            ALLOCATE(X3D(4,nu_local,nvp))
            ALLOCATE(Y3D(4,nu_local,nvp))
            ALLOCATE(Z3D(4,nu_local,nvp))
            ALLOCATE(KX3D(4,nu_local,nvp))
            ALLOCATE(KY3D(4,nu_local,nvp))
            ALLOCATE(KZ3D(4,nu_local,nvp))
#if defined(MPI_OPT)
         ENDIF
#endif
         IF (shar_rank == 0) THEN
            ALLOCATE(xu(nu_local),xv(nv_local))
            ALLOCATE(rreal(nu_local,nv_local), &
                     zreal(nu_local,nv_local), &
                     rureal(nu_local,nv_local), &
                     rvreal(nu_local,nv_local), &
                     zureal(nu_local,nv_local), &
                     zvreal(nu_local,nv_local), &
                     sxreal(nu_local,nv_local), &
                     syreal(nu_local,nv_local), &
                     szreal(nu_local,nv_local), &
                     potu(nu_local,nv_local), &
                     potv(nu_local,nv_local), &
                     potr(nu_local,nv_local), &
                     potp(nu_local,nv_local), &
                     potz(nu_local,nv_local), &
                     potx(nu_local,nv_local), &
                     poty(nu_local,nv_local))
            ALLOCATE(fmn_temp(mnmax_surface))
            FORALL(u=1:nu_local) xu(u) = DBLE(u-1)/DBLE(nu_local-1)
            FORALL(v=1:nv_local) xv(v) = DBLE(v-1)/DBLE(nv_local-1)
            rreal = zero; rureal = zero; rvreal = zero
            zreal = zero; zureal = zero; zvreal = zero
            potu = zero; potv = zero
            CALL mntouv_local(mnmax_surface,nu_local,nv_local,xu,xv,            &
                              rmnc_surface,xm_surface,xn_surface,  &
                              rreal,0,1)
            CALL mntouv_local(mnmax_surface,nu_local,nv_local,xu,xv,            &
                              zmns_surface,xm_surface,xn_surface,  &
                              zreal,1,0)
            fmn_temp = -rmnc_surface*xm_surface
            CALL mntouv_local(mnmax_surface,nu_local,nv_local,xu,xv,fmn_temp,xm_surface,xn_surface,rureal,1,0)
            fmn_temp = -rmnc_surface*xn_surface
            CALL mntouv_local(mnmax_surface,nu_local,nv_local,xu,xv,fmn_temp,xm_surface,xn_surface,rvreal,1,0)  
            fmn_temp =  zmns_surface*xm_surface
            CALL mntouv_local(mnmax_surface,nu_local,nv_local,xu,xv,fmn_temp,xm_surface,xn_surface,zureal,0,0) 
            fmn_temp =  zmns_surface*xn_surface
            CALL mntouv_local(mnmax_surface,nu_local,nv_local,xu,xv,fmn_temp,xm_surface,xn_surface,zvreal,0,0)
            DEALLOCATE(fmn_temp)
            ALLOCATE(fmn_temp(mnmax_pot))
            fmn_temp =  potmns_surface*xm_pot
            CALL mntouv_local(mnmax_pot,nu_local,nv_local,xu,xv,fmn_temp,xm_pot,xn_pot,potu,0,1)
            fmn_temp =  potmns_surface*xn_pot
            CALL mntouv_local(mnmax_pot,nu_local,nv_local,xu,xv,fmn_temp,xm_pot,xn_pot,potv,0,0)
            DEALLOCATE(xu,xv,fmn_temp)



            ! Correct derivatives for missing pi2
            rureal = pi2*rureal
            rvreal = pi2*rvreal
            zureal = pi2*zureal
            zvreal = pi2*zvreal
            potu   = potu
            potv   = potv

            ! Add secular pieces
            potu = potu - cut
            potv = potv - cup

            ! Calculate surface coords and normals
            ALLOCATE(xu(nv_local),xv(nv_local),yu(nv_local),yv(nv_local),cop(nv_local),sip(nv_local))
            FORALL(v=1:nv_local) cop(v) = COS(alp*DBLE(v-1)/DBLE(nv_local-1))
            FORALL(v=1:nv_local) sip(v) = SIN(alp*DBLE(v-1)/DBLE(nv_local-1))
            DO u = 1, nu_local
               X3D(1,u,1:nv_local) = rreal(u,:)*cop
               Y3D(1,u,1:nv_local) = rreal(u,:)*sip
               xu       = rureal(u,:)*cop
               yu       = rureal(u,:)*sip
               xv       = rvreal(u,:)*cop - rreal(u,:)*sip*alp
               yv       = rvreal(u,:)*sip + rreal(u,:)*cop*alp
               sxreal(u,:) = -yu(:)*zvreal(u,:) + zureal(u,:)*yv(:)
               syreal(u,:) = -xv(:)*zureal(u,:) + zvreal(u,:)*xu(:)
               szreal(u,:) = -xu(:)*yv(:)       + yu(:)*xv(:)
               ! Potential
               !potr(u,:) = potu(u,:)*rureal(u,:) + potv(u,:)*rvreal(u,:)
               !potp(u,:) = potv(u,:)*rreal(u,:)*alp
               !potx(u,:) = potr(u,:)*cop - potp(u,:)*sip
               !poty(u,:) = potr(u,:)*sip + potp(u,:)*cop
               potx(u,:) = potu(u,:)*xu          + potv(u,:)*xv
               poty(u,:) = potu(u,:)*yu          + potv(u,:)*yv
               potz(u,:) = potu(u,:)*zureal(u,:) - potv(u,:)*zvreal(u,:)
            END DO
            WRITE(327,*) X3D(1,:,1:nv_local)
            WRITE(328,*) Y3D(1,:,1:nv_local)
            WRITE(329,*) zreal
            u = nu_local - 1
            v = nv_local - 1
            surf_area = np*SUM(SQRT( sxreal(1:u,1:v)**2+syreal(1:u,1:v)**2+szreal(1:u,1:v)**2))/(u*v)
            PRINT *,surf_area
            Z3D(1,:,1:nv_local) = zreal
            KX3D(1,:,1:nv_local) = syreal*potz - szreal*poty
            KY3D(1,:,1:nv_local) = szreal*potx - sxreal*potz
            KZ3D(1,:,1:nv_local) = sxreal*poty - syreal*potx
            WRITE(337,*) KX3D(1,:,1:nv_local)
            WRITE(338,*) KY3D(1,:,1:nv_local)
            WRITE(339,*) KZ3D(1,:,1:nv_local)
            DEALLOCATE(xu,xv,yu,yv,cop,sip)
            DEALLOCATE(sxreal,syreal,szreal,rureal,rvreal,zureal,zvreal,potu,potv,potr,potp,potz,potx,poty)

            ! Now extend to more field periods
            ALLOCATE(cop(np),sip(np))
            FORALL(v=1:np) cop(v) = COS(alp*DBLE(v-1))
            FORALL(v=1:np) sip(v) = SIN(alp*DBLE(v-1))
            i1 = nv_local
            i2 = i1 + nv_local - 1
            DO v=2,np
               X3D(1,:,i1:i2)  =  X3D(1,:,1:nv_local)*cop(v) -  Y3D(1,:,1:nv_local)*sip(v)
               Y3D(1,:,i1:i2)  =  X3D(1,:,1:nv_local)*sip(v) +  Y3D(1,:,1:nv_local)*cop(v)
               Z3D(1,:,i1:i2)  =  Z3D(1,:,1:nv_local)
               KX3D(1,:,i1:i2) = KX3D(1,:,1:nv_local)*cop(v) - KY3D(1,:,1:nv_local)*sip(v)
               KY3D(1,:,i1:i2) = KX3D(1,:,1:nv_local)*sip(v) + KY3D(1,:,1:nv_local)*cop(v)
               KZ3D(1,:,i1:i2) = KZ3D(1,:,1:nv_local)
               i1 = i2
               i2 = i1 + nv_local - 1
            END DO
            DEALLOCATE(cop,sip)

            ! Create splines
            CALL EZspline_init(f_spl,nu_local,nvp,bcs1,bcs2,ier)
            f_spl%isHermite  = 0
            CALL EZspline_setup(f_spl,X3D(1,:,:),ier)
            X3D  = f_spl%fspl
            x1   = f_spl%x1
            x2   = f_spl%x2
            CALL EZspline_free(f_spl,ier)
            CALL EZspline_init(f_spl,nu_local,nvp,bcs1,bcs2,ier)
            f_spl%isHermite  = 0
            CALL EZspline_setup(f_spl,Y3D(1,:,:),ier)
            Y3D  = f_spl%fspl
            CALL EZspline_free(f_spl,ier)
            CALL EZspline_init(f_spl,nu_local,nvp,bcs1,bcs2,ier)
            f_spl%isHermite  = 0
            CALL EZspline_setup(f_spl,Z3D(1,:,:),ier)
            Z3D  = f_spl%fspl
            CALL EZspline_free(f_spl,ier)
            CALL EZspline_init(f_spl,nu_local,nvp,bcs1,bcs2,ier)
            f_spl%isHermite  = 0
            CALL EZspline_setup(f_spl,KX3D(1,:,:),ier)
            KX3D  = f_spl%fspl
            CALL EZspline_free(f_spl,ier)
            CALL EZspline_init(f_spl,nu_local,nvp,bcs1,bcs2,ier)
            f_spl%isHermite  = 0
            CALL EZspline_setup(f_spl,KY3D(1,:,:),ier)
            KY3D  = f_spl%fspl
            CALL EZspline_free(f_spl,ier)
            CALL EZspline_init(f_spl,nu_local,nvp,bcs1,bcs2,ier)
            f_spl%isHermite  = 0
            CALL EZspline_setup(f_spl,KZ3D(1,:,:),ier)
            KZ3D  = f_spl%fspl
            CALL EZspline_free(f_spl,ier)
         END IF
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_BARRIER(shar_comm,ier)
         CALL MPI_COMM_FREE(shar_comm,ier)
      END IF
#endif
      RETURN
      END SUBROUTINE nescoil_bfield_init

      SUBROUTINE nescoil_bfield_dbl(x,y,z,bx,by,bz)
         ! Based on BESURFCUR
         IMPLICIT NONE
         REAL(rprec), INTENT(in) :: x,y,z
         REAL(rprec), INTENT(out) :: bx,by,bz 
         ! LOCAL VARIABLES
         DOUBLE PRECISION ::  gf(nu_int-1,nvp-1), gf3(nu_int-1,nvp-1)

         gf = one / SQRT(   (x - X3D(1,1:u1,1:v1))**2 &
                          + (y - Y3D(1,1:u1,1:v1))**2 &
                          + (z - Z3D(1,1:u1,1:v1))**2 )
         gf3 = gf*gf*gf*norm
         bx  = SUM((   KY3D(1,1:u1,1:v1) * (z - Z3D(1,1:u1,1:v1)) &
                            - KZ3D(1,1:u1,1:v1) * (y - Y3D(1,1:u1,1:v1)) ) *gf3 )
         by  = SUM((   KZ3D(1,1:u1,1:v1) * (x - X3D(1,1:u1,1:v1)) &
                            - KX3D(1,1:u1,1:v1) * (z - Z3D(1,1:u1,1:v1)) ) *gf3 )
         bz  = SUM((   KX3D(1,1:u1,1:v1) * (y - Y3D(1,1:u1,1:v1)) &
                            - KY3D(1,1:u1,1:v1) * (x - X3D(1,1:u1,1:v1)) ) *gf3 )
         RETURN
      END SUBROUTINE nescoil_bfield_dbl

      SUBROUTINE nescoil_bfield_flt(x_flt,y_flt,z_flt,bx_flt,by_flt,bz_flt)
         IMPLICIT NONE
         ! INPUT VARIABLES
         REAL, INTENT(in)  :: x_flt, y_flt, z_flt
         REAL, INTENT(out) :: bx_flt, by_flt, bz_flt
         ! LOCAL VARIABLES
         DOUBLE PRECISION :: xt,yt,zt,bxt,byt,bzt
         ! BEGIN SUBROUTINE
         xt  = x_flt
         yt  = y_flt
         zt  = z_flt
         bxt = zero
         byt = zero
         bzt = zero
         CALL nescoil_bfield_dbl(xt,yt,zt,bxt,byt,bzt)
         bx_flt = bxt
         by_flt = byt
         bz_flt = bzt
         RETURN
      END SUBROUTINE nescoil_bfield_flt

      SUBROUTINE nescoil_bfield_adapt_dbl(x,y,z,bx,by,bz,istat)
         IMPLICIT NONE
         ! INPUT VARIABLES
         DOUBLE PRECISION, INTENT(in)  :: x, y, z
         DOUBLE PRECISION, INTENT(out) :: bx, by, bz
         INTEGER, INTENT(inout) :: istat
         ! LOCAL VARIABLES
         LOGICAL            :: adapt_rerun
         INTEGER(KIND=8), PARAMETER :: ndim = 2 ! theta,zeta
         INTEGER(KIND=8), PARAMETER :: nfun = 3 ! Bx, By, Bz
         INTEGER(KIND=8) :: maxcls,mincls, restar, wrklen, funcls
         DOUBLE PRECISION :: absreq, relreq
         DOUBLE PRECISION :: a(ndim), b(ndim), &
                              finest(nfun), absest(nfun)
         DOUBLE PRECISION, ALLOCATABLE :: vrtwrk(:)

         EXTERNAL :: dcuhre

         ! BEGIN SUBROUTINE
         IF (istat == -327) THEN
            CALL nescoil_bfield(x,y,z,bx,by,bz)
            RETURN
         END IF

         a(1:2) = zero
         b(1:2) = one
         mincls = 0
         maxcls = 16777216
         absreq = 1.0E-6
         relreq = 1.0E-2
         finest = zero
         absest = zero
         x_f      = x
         y_f      = y
         z_f      = z
         adapt_rerun = .true.
         restar = 0
         DO WHILE (adapt_rerun) 
            IF (.not.ALLOCATED(vrtwrk)) THEN
               wrklen = ((maxcls-ndim)/(2*ndim) + 1)*(2*ndim+2*nfun+2) + 17*nfun + 1
               ALLOCATE(vrtwrk(wrklen),STAT=istat)
               IF (istat .ne. 0) THEN
                  WRITE(6,*) ' ALLOCATION ERROR IN: nescoil_bfield_adapt_dbl'
                  WRITE(6,*) '   VARIABLE: VRTWRK, SIZE: ',wrklen
                  WRITE(6,*) '   ISTAT: ',istat
                  RETURN
               END IF
            END IF
            CALL dcuhre(ndim,nfun,a,b,mincls,maxcls,funsub_b,absreq,&
                        relreq,0,wrklen,restar,finest,absest,funcls,istat,vrtwrk)
            IF (istat == 1) THEN
               ! For now we don't try to restart and just live with the result
               bx = finest(1)
               by = finest(2)
               bz = finest(3)
               adapt_rerun=.false.
               DEALLOCATE(vrtwrk)
            ELSE IF (istat > 1) THEN
               bx = zero
               by = zero
               bz = zero
               adapt_rerun=.false.
               DEALLOCATE(vrtwrk)
            ELSE
               bx = finest(1)
               by = finest(2)
               bz = finest(3)
               adapt_rerun=.false.
               DEALLOCATE(vrtwrk)
            END IF
         END DO
         RETURN
      END SUBROUTINE nescoil_bfield_adapt_dbl

      SUBROUTINE nescoil_bfield_adapt_flt(x_flt,y_flt,z_flt,bx_flt,by_flt,bz_flt,istat)
         IMPLICIT NONE
         ! INPUT VARIABLES
         REAL, INTENT(in)  :: x_flt, y_flt, z_flt
         REAL, INTENT(out) :: bx_flt, by_flt, bz_flt
         INTEGER, INTENT(inout) :: istat
         ! LOCAL VARIABLES
         DOUBLE PRECISION :: xt,yt,zt,bxt,byt,bzt
         ! BEGIN SUBROUTINE
         xt  = x_flt
         yt  = y_flt
         zt  = z_flt
         bxt = zero
         byt = zero
         bzt = zero
         CALL nescoil_bfield_adapt_dbl(xt,yt,zt,bxt,byt,bzt,istat)
         bx_flt = bxt
         by_flt = byt
         bz_flt = bzt
         RETURN
      END SUBROUTINE nescoil_bfield_adapt_flt
         
      SUBROUTINE funsub_b(ndim, vec, nfun, f)
         IMPLICIT NONE
         ! INPUT VARIABLES
         INTEGER, INTENT(in) :: ndim, nfun
         DOUBLE PRECISION, INTENT(in) :: vec(ndim)
         DOUBLE PRECISION, INTENT(out) :: f(nfun)
         ! LOCAL VARIABLES
         INTEGER :: ier
         DOUBLE PRECISION :: bn, xs, ys, zs, gf, gf3, nx, ny, &
                             nz, kx, ky, kz
         INTEGER :: i,j
         REAL*8 :: xparam, yparam, hx, hy, hxi, hyi
         REAL*8 :: xpi, xp2, xpi2, ypi, yp2, ypi2
         REAL*8 :: cx,cxi,hx2,cy,cyi,hy2 ! Non Hermite quantities
         ! BEGIN SUBROUTINE
         CALL lookupgrid2d(vec(1),vec(2),i,j,hx,hy,hxi,hyi,xparam,yparam)
         xpi  = one - xparam;    ypi  = one - yparam
         xp2  = xparam * xparam; yp2  = yparam * yparam
         xpi2 = xpi * xpi;       ypi2 = ypi * ypi
         ! non Hermite Quatitites
         cx = xparam*(xp2-1); cxi = xpi*(xpi2-1); hx2 = hx*hx
         cy = yparam*(yp2-1); cyi = ypi*(ypi2-1); hy2 = hy*hy
         xs  =  evalbi2D(xparam, xpi, xp2, xpi2, cx, cxi, hx2, yparam, ypi, yp2, ypi2, cy, cyi, hy2, nx1, nx2,  X3D, i, j)
         ys  =  evalbi2D(xparam, xpi, xp2, xpi2, cx, cxi, hx2, yparam, ypi, yp2, ypi2, cy, cyi, hy2, nx1, nx2,  Y3D, i, j)
         zs  =  evalbi2D(xparam, xpi, xp2, xpi2, cx, cxi, hx2, yparam, ypi, yp2, ypi2, cy, cyi, hy2, nx1, nx2,  Z3D, i, j)
         kx  =  evalbi2D(xparam, xpi, xp2, xpi2, cx, cxi, hx2, yparam, ypi, yp2, ypi2, cy, cyi, hy2, nx1, nx2, KX3D, i, j)
         ky  =  evalbi2D(xparam, xpi, xp2, xpi2, cx, cxi, hx2, yparam, ypi, yp2, ypi2, cy, cyi, hy2, nx1, nx2, KY3D, i, j)
         kz  =  evalbi2D(xparam, xpi, xp2, xpi2, cx, cxi, hx2, yparam, ypi, yp2, ypi2, cy, cyi, hy2, nx1, nx2, KZ3D, i, j)

         gf   = one/DSQRT((x_f-xs)*(x_f-xs)+(y_f-ys)*(y_f-ys)+(z_f-zs)*(z_f-zs))
         gf3  = norm_fsub*gf*gf*gf
         f(1) = (ky*(z_f-zs)-kz*(y_f-ys))*gf3
         f(2) = (kz*(x_f-xs)-kx*(z_f-zs))*gf3
         f(3) = (kx*(y_f-ys)-ky*(x_f-xs))*gf3
         !WRITE(427,*) vec(1),vec(2),xs,ys,zs
         RETURN
         ! END SUBROUTINE
      END SUBROUTINE funsub_b

      SUBROUTINE read_nescout_deallocate(comm)
         USE mpi_sharmem
         USE mpi_inc
         IMPLICIT NONE
         INTEGER, INTENT(in), OPTIONAL :: comm
         INTEGER :: ier
         IF (ALLOCATED(xm_plasma)) DEALLOCATE(xm_plasma)
         IF (ALLOCATED(xn_plasma)) DEALLOCATE(xn_plasma)
         IF (ALLOCATED(rmnc_plasma)) DEALLOCATE(rmnc_plasma)
         IF (ALLOCATED(zmns_plasma)) DEALLOCATE(zmns_plasma)
         IF (ALLOCATED(rmns_plasma)) DEALLOCATE(rmns_plasma)
         IF (ALLOCATED(zmnc_plasma)) DEALLOCATE(zmnc_plasma)
         IF (ALLOCATED(lmnc_plasma)) DEALLOCATE(lmnc_plasma)
         IF (ALLOCATED(lmns_plasma)) DEALLOCATE(lmns_plasma)
         IF (ALLOCATED(xm_surface)) DEALLOCATE(xm_surface)
         IF (ALLOCATED(xn_surface)) DEALLOCATE(xn_surface)
         IF (ALLOCATED(rmnc_surface)) DEALLOCATE(rmnc_surface)
         IF (ALLOCATED(zmns_surface)) DEALLOCATE(zmns_surface)
         IF (ALLOCATED(rmns_surface)) DEALLOCATE(rmns_surface)
         IF (ALLOCATED(zmnc_surface)) DEALLOCATE(zmnc_surface)
         IF (ALLOCATED(potmns_surface)) DEALLOCATE(potmns_surface)
         IF (ALLOCATED(xm_pot)) DEALLOCATE(xm_pot)
         IF (ALLOCATED(xn_pot)) DEALLOCATE(xn_pot)
         IF (ALLOCATED(x_plasma)) DEALLOCATE(x_plasma)
         IF (ALLOCATED(y_plasma)) DEALLOCATE(y_plasma)
         IF (ALLOCATED(z_plasma)) DEALLOCATE(z_plasma)
         IF (ALLOCATED(dsur_plasma)) DEALLOCATE(dsur_plasma)
         IF (ALLOCATED(nx_plasma)) DEALLOCATE(nx_plasma)
         IF (ALLOCATED(ny_plasma)) DEALLOCATE(ny_plasma)
         IF (ALLOCATED(nz_plasma)) DEALLOCATE(nz_plasma)
         IF (ALLOCATED(dxdu_plasma)) DEALLOCATE(dxdu_plasma)
         IF (ALLOCATED(dydu_plasma)) DEALLOCATE(dydu_plasma)
         IF (ALLOCATED(dxdv_plasma)) DEALLOCATE(dxdv_plasma)
         IF (ALLOCATED(dydv_plasma)) DEALLOCATE(dydv_plasma)
         IF (ALLOCATED(x_surface)) DEALLOCATE(x_surface)
         IF (ALLOCATED(y_surface)) DEALLOCATE(y_surface)
         IF (ALLOCATED(z_surface)) DEALLOCATE(z_surface)
         IF (ALLOCATED(dsur_surface)) DEALLOCATE(dsur_surface)
         IF (ALLOCATED(nx_surface)) DEALLOCATE(nx_surface)
         IF (ALLOCATED(ny_surface)) DEALLOCATE(ny_surface)
         IF (ALLOCATED(nz_surface)) DEALLOCATE(nz_surface)
         IF (ALLOCATED(db_normal)) DEALLOCATE(db_normal)
         IF (ALLOCATED(Babs)) DEALLOCATE(Babs)
#if defined(MPI_OPT)
         IF (PRESENT(comm)) THEN
            CALL MPI_BARRIER(comm,ier)
            IF (ASSOCIATED(X3D))  CALL mpidealloc(X3D,  win_X3D)
            IF (ASSOCIATED(Y3D))  CALL mpidealloc(Y3D,  win_Y3D)
            IF (ASSOCIATED(Z3D))  CALL mpidealloc(Z3D,  win_Z3D)
            IF (ASSOCIATED(KX3D)) CALL mpidealloc(KX3D,  win_KX3D)
            IF (ASSOCIATED(KY3D)) CALL mpidealloc(KY3D,  win_KY3D)
            IF (ASSOCIATED(KZ3D)) CALL mpidealloc(KZ3D,  win_KZ3D)
            IF (ASSOCIATED(x1))   CALL mpidealloc(x1,  win_x1)
            IF (ASSOCIATED(x2))   CALL mpidealloc(x2,  win_x2)
            CALL MPI_BARRIER(comm,ier)
         ELSE
#endif
            IF (ASSOCIATED(X3D)) DEALLOCATE(X3D)
            IF (ASSOCIATED(Y3D)) DEALLOCATE(Y3D)
            IF (ASSOCIATED(Z3D)) DEALLOCATE(Z3D)
            IF (ASSOCIATED(KX3D)) DEALLOCATE(KX3D)
            IF (ASSOCIATED(KY3D)) DEALLOCATE(KY3D)
            IF (ASSOCIATED(KZ3D)) DEALLOCATE(KZ3D)
            IF (ASSOCIATED(x1)) DEALLOCATE(x1)
            IF (ASSOCIATED(x2)) DEALLOCATE(x2)
#if defined(MPI_OPT)
      END IF
#endif         
         RETURN
      END SUBROUTINE read_nescout_deallocate
         
      SUBROUTINE mntouv_local(mnmax,nul,nvl,xu,xv,fmn,xm,xn,f,signs,calc_trig)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in) :: mnmax
      INTEGER, INTENT(in) :: nul
      INTEGER, INTENT(in) :: nvl
      DOUBLE PRECISION, INTENT(in) :: xu(1:nul)
      DOUBLE PRECISION, INTENT(in) :: xv(1:nvl)           
      DOUBLE PRECISION, INTENT(in) :: fmn(1:mnmax)
      INTEGER, INTENT(in) :: xm(1:mnmax)
      INTEGER, INTENT(in) :: xn(1:mnmax)
      DOUBLE PRECISION, INTENT(inout) :: f(1:nul,1:nvl)
      INTEGER, INTENT(in) :: signs
      INTEGER, INTENT(in) :: calc_trig
      ! LOCAL VARIABLES
      INTEGER     :: mn, i, ier, ik
      DOUBLE PRECISION :: xm_temp(1:mnmax,1)
      DOUBLE PRECISION :: xn_temp(1:mnmax,1)
      DOUBLE PRECISION :: pi2_l
      DOUBLE PRECISION :: mt(1:mnmax,1:nul)
      DOUBLE PRECISION :: nz(1:mnmax,1:nvl)
      DOUBLE PRECISION :: fmn_temp(1:mnmax,1:nul)
      DOUBLE PRECISION :: xu_temp(1,1:nul)
      DOUBLE PRECISION :: xv_temp(1,1:nvl)
      DOUBLE PRECISION :: fmn_help(1:mnmax)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: cosmt(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: sinmt(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: cosnz(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: sinnz(:,:)
      ! BEGIN SUBROUTINE
      pi2_l = 8.0D+0 * ATAN(1.)
      IF (calc_trig == 1) THEN
         IF (ALLOCATED(cosmt)) DEALLOCATE(cosmt)
         IF (ALLOCATED(sinmt)) DEALLOCATE(sinmt)
         IF (ALLOCATED(cosnz)) DEALLOCATE(cosnz)
         IF (ALLOCATED(sinnz)) DEALLOCATE(sinnz)
         ALLOCATE(cosmt(1:mnmax,1:nul),sinmt(1:mnmax,1:nul),&
                  cosnz(1:mnmax,1:nvl),sinnz(1:mnmax,1:nvl),STAT=ier)
         FORALL(i=1:mnmax) xm_temp(i,1)=DBLE(xm(i))
         FORALL(i=1:mnmax) xn_temp(i,1)=DBLE(xn(i))
         FORALL(i=1:nul) xu_temp(1,i)=xu(i)
         FORALL(i=1:nvl) xv_temp(1,i)=xv(i)
         mt = MATMUL(xm_temp,xu_temp)
         nz = MATMUL(xn_temp,xv_temp)
         FORALL(mn=1:mnmax,i=1:nul) cosmt(mn,i) = dcos(pi2_l*mt(mn,i))
         FORALL(mn=1:mnmax,i=1:nul) sinmt(mn,i) = dsin(pi2_l*mt(mn,i))
         FORALL(mn=1:mnmax,i=1:nvl) cosnz(mn,i) = dcos(pi2_l*nz(mn,i))
         FORALL(mn=1:mnmax,i=1:nvl) sinnz(mn,i) = dsin(pi2_l*nz(mn,i))
      END IF
      fmn_temp=SPREAD(fmn,2,nul)
      IF (SIGNS == 0) THEN
         f(1:nul,1:nvl) = f(1:nul,1:nvl) + MATMUL(TRANSPOSE((fmn_temp*cosmt)),cosnz) &
                                         - MATMUL(TRANSPOSE((fmn_temp*sinmt)),sinnz)
      ELSE IF (SIGNS == 1) THEN
         f(1:nul,1:nvl) = f(1:nul,1:nvl) + MATMUL(TRANSPOSE((fmn_temp*sinmt)),cosnz) &
                                         + MATMUL(TRANSPOSE((fmn_temp*cosmt)),sinnz)
      END IF
      ! END SUBROUTINE
      END SUBROUTINE mntouv_local

      SUBROUTINE lookupgrid2d(x1_in,x2_in,i,j,hx,hy,hxi,hyi,xparam,yparam)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x1_in, x2_in
      INTEGER, INTENT(out) :: i,j
      REAL*8, INTENT(out) :: hx, hy, hxi, hyi, xparam, yparam
      i = MIN(MAX(COUNT(x1 < x1_in),1),nx1-1)
      j = MIN(MAX(COUNT(x2 < x2_in),1),nx2-1)
      hx     = x1(i+1) - x1(i)
      hy     = x2(j+1) - x2(j)
      hxi    = one / hx
      hyi    = one / hy
      xparam = (x1_in - x1(i)) * hxi
      yparam = (x2_in - x2(j)) * hyi
      RETURN
      END SUBROUTINE lookupgrid2d

      REAL*8 FUNCTION evalbi2D(xp,xpi,xp2,xpi2,cx,cxi,hx2,yp,ypi,yp2,ypi2,cy,cyi,hy2, n1, n2, F, i, j)
      IMPLICIT NONE
      REAL*8, INTENT(in) :: xp, xpi, xp2, xpi2, cx, cxi, hx2
      REAL*8, INTENT(in) :: yp, ypi, yp2, ypi2, cy, cyi, hy2
      INTEGER, INTENT(in) :: n1, n2, i, j
      REAL*8, INTENT(in) :: F(4,n1,n2)
      REAL*8, PARAMETER :: sixth = 0.166666666666666667
      REAL*8, PARAMETER :: z36th = 0.027777777777777776

      evalbi2D = xpi*(ypi*F(1,i,j)  +yp*F(1,i,j+1))+ &
                  xp*(ypi*F(1,i+1,j)+yp*F(1,i+1,j+1)) &
                +sixth*hx2*( &
                      cxi*(ypi*F(2,i,j)  +yp*F(2,i,j+1))+ &
                      cx*(ypi*F(2,i+1,j)+yp*F(2,i+1,j+1)) ) &
                +sixth*hy2*( &
                      xpi*(cyi*F(3,i,j)  +cy*F(3,i,j+1))+ &
                       xp*(cyi*F(3,i+1,j)+cy*F(3,i+1,j+1))) &
                     +z36th*hx2*hy2*( &
                              cxi*(cyi*F(4,i,j)  +cy*F(4,i,j+1))+ &
                               cx*(cyi*F(4,i+1,j)+cy*F(4,i+1,j+1)) ) 
      RETURN
      END FUNCTION evalbi2D

      END MODULE read_nescoil_mod