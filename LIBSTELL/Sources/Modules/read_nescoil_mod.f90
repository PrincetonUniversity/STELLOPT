!-----------------------------------------------------------------------
!     Module:        read_nescoil_mod
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          06/07/2024
!     Description:   This module stores routines for reading
!                    NESCOIL data.
!-----------------------------------------------------------------------
      MODULE read_nescoil_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE safe_open_mod

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
                 mnmax_surface, nmax, mnd, nuv, nuv1, nuvh, nuvh1
      REAL(rprec) :: iota_edge, phip_edge, curpol, cut, cup, curwt, &
                     trgwt
      INTEGER, DIMENSION(:), ALLOCATABLE ::xm_plasma, xn_plasma,      & 
                                          xm_surface, xn_surface 
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::   &
                                          rmnc_plasma, zmns_plasma,   &
                                          rmns_plasma, zmnc_plasma,   &
                                          lmnc_plasma, lmns_plasma,       &
                                          rmnc_surface, zmns_surface, &
                                          rmns_surface, zmnc_surface, &
                                          potmnc_surface
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
                                          dxdu_surface, dydu_surface, &
                                          dxdv_surface, dydv_surface, &
                                          db_normal, Babs
      CHARACTER(LEN=256), PRIVATE :: line


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
               CASE("----- Plasma information from VMEC -----")
                  READ(iunit, '(A)', iostat=istat) line
                  READ(iunit, '(i6,3g25.16)') np, iota_edge, phip_edge, curpol
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
                  READ(iunit, '(A)', iostat=istat) line
                  ALLOCATE(bn_plasma(nuvh1))
                  READ(iunit,*) bn_plasma
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
                  IF (w_csurf>2 .or. w_csurf==-3) CALL read_nescout_surfacederiv(iunit,istat)
               CASE("---- Phi(m,n) for least squares ---")
                  ALLOCATE(potmnc_surface(mnmax_surface))
                  DO n = 1, mnmax_surface
                     READ(iunit,'(i3,2x,i3,2x,g25.16)') mtemp,ntemp,potmnc_surface(n)
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

      SUBROUTINE read_nescout_surfacederiv(iunit,istat)
         IMPLICIT NONE
         INTEGER, INTENT(inout) :: iunit
         INTEGER, INTENT(inout) :: istat
         INTEGER :: mn
         READ(iunit, '(A)', iostat=istat) line
         ALLOCATE(dxdu_surface(nuvh),dydu_surface(nuvh))
         ALLOCATE(dxdv_surface(nuvh),dydv_surface(nuvh))
         DO mn = 1, nuvh
            READ(iunit,'(4g16.6)') dxdu_surface(mn), dydu_surface(mn), &
                                   dxdv_surface(mn), dydv_surface(mn)
         END DO
         READ(iunit, '(A)', iostat=istat) line
         RETURN
      END SUBROUTINE read_nescout_surfacederiv

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

      SUBROUTINE read_nescout_deallocate
         IMPLICIT NONE
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
         IF (ALLOCATED(potmnc_surface)) DEALLOCATE(potmnc_surface)
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
         IF (ALLOCATED(dxdu_surface)) DEALLOCATE(dxdu_surface)
         IF (ALLOCATED(dydu_surface)) DEALLOCATE(dydu_surface)
         IF (ALLOCATED(dxdv_surface)) DEALLOCATE(dxdv_surface)
         IF (ALLOCATED(dydv_surface)) DEALLOCATE(dydv_surface)
         IF (ALLOCATED(db_normal)) DEALLOCATE(db_normal)
         IF (ALLOCATED(Babs)) DEALLOCATE(Babs)
         RETURN
      END SUBROUTINE read_nescout_deallocate

      END MODULE read_nescoil_mod