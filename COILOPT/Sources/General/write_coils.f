      SUBROUTINE write_coils (extension)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE boundary
      USE modular_coils
      USE saddle_coils
      USE tf_coils
      USE vf_coils
      USE coiltypes
      USE coils
      USE bcoils_mod
      USE Vwire
      USE safe_open_mod
      USE Vname, ONLY: lgeom_only
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*) :: extension
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: cunit=15, munit, gunit, sunit, tunit, ierr
      INTEGER :: n, i, j, js, nc, ncm, ncp, mbw, nmax, ncoil_count
      INTEGER :: nz
      REAL(rprec) :: ph_w, th_w, zz, sfils, rps, zps, plascur, dph_w
      REAL(rprec), DIMENSION (nwdim1) :: xw, yw, zw
      CHARACTER(LEN=7) :: status
      CHARACTER(LEN=10) :: ch_sad
!-----------------------------------------------
      IF (lgeom_only) THEN
         status = 'scratch'
      ELSE
         status = 'replace'
      END IF

      CALL safe_open(cunit, ierr, 'coilxyz.dat', status, 'formatted')
      IF (ierr .ne. 0) STOP 'Error opening coilxyz.dat in write_coils'

      munit = cunit+1
      CALL safe_open(munit, ierr, 'coils.' // TRIM(extension),
     1    'replace', 'formatted')
      IF (ierr .ne. 0)STOP 'Error opening coils-dot file in write_coils'

      gunit = munit+1
      CALL safe_open(gunit, ierr, 'coilgnu.dat', status, 'formatted')
      IF (ierr .ne. 0) STOP 'Error opening coilgnu.dat in write_coils'

      sunit = gunit+1
      CALL safe_open(sunit, ierr, 'coilsad.dat', status, 'formatted')
      IF (ierr .ne. 0) STOP 'Error opening coilsad.dat in write_coils'

      IF (lspline .and. lctrlpt) THEN
        tunit = sunit+1
        CALL safe_open(tunit, ierr, 'ctrlpts.dat', status, 'formatted')
        IF (ierr .ne. 0) STOP 'Error opening ctrlpts.dat in write_coils'
      END IF

      zz = 0
      nz = 0

      WRITE (munit,300) nfp
      WRITE (munit,310)
      WRITE (munit,320)
  300 FORMAT("periods ",i2)
  310 FORMAT("begin filament")
  320 FORMAT("mirror NUL")

      nmax = 0
      ncoil_count = nmax
      ncp = 0

      IF (lmodular) THEN
!     WRITE modular coils, currents
        DO n = 1, nmod_coils
           WRITE(cunit,'(1p,e22.14,2i6)') curcon(n), nwire, nz
           WRITE(gunit,'(/)')
           ncm = 1 + MOD(n-1,nmod_coils_per_period)
           ncp = ncp + 1
           IF (ncp .eq. nmid) ncp = 0
!          compute coil group number
           IF (lmodcur) THEN
!             ONLY good for lsymm = t, nodd = 1 for now
              nc = ncoil_count + ncm
              IF (ncm .ge. nmid) nc = ncoil_count + nmid - ncp
              IF (ncm .eq. nmod_coils_per_period) ncp = 0
              nmax = MAX(nmax, nc)
           ELSE
              nc = ncoil_count + 1
              nmax = nc
           END IF

           DO i = 1,nwire1
             WRITE(cunit,'(1p,3e22.14)')
     1          x_mod(i,1,n), y_mod(i,1,n), z_mod(i,1,n)
             IF (i .eq. nwire1) THEN
                WRITE(munit,'(1p,4e22.14,i4,a9)')
     1             x_mod(1,1,n), y_mod(1,1,n), z_mod(1,1,n), zero, nc,
     2             "  modular"
             ELSE
                WRITE(munit,'(1p,4e22.14)')
     1             x_mod(i,1,n), y_mod(i,1,n), z_mod(i,1,n), curcon(n)
             END IF
             ph_w = modular(nc)%phi(i) + phi_full(n)
             th_w = twopi*(i-1)/nwire
             WRITE(gunit,'(3f14.6,i6,2e14.6)')
     1         x_mod(i,1,n), y_mod(i,1,n), z_mod(i,1,n),i,
     2         ph_w/(twopi/nfper), th_w/twopi
           END DO
        END DO
      END IF

      IF (lsaddle) THEN
        ncoil_count = nmax
        ncp = 0
        IF (lsmod) THEN
           ch_sad = ' modular'
        ELSE
           ch_sad = ' Saddle'
        END IF

        js = 1
        IF (nfils .gt. 4) THEN
          js = 2
          sfils = 4.0_dp
        ELSE IF (nfils .gt. 1) THEN
          sfils = 3.0_dp
        ELSE
          sfils = 1.0_dp
        END IF

!       Write saddle coils, currents
        DO n = 1, nsad_coils
           ncm = 1 + MOD(n-1,nsad_coils_per_period)
           ncp = 1 + MOD(n-1,nsmid)
!          compute coil group number
           nc = ncoil_count + nsad_group(ncp)
           IF (ncm .gt. nsmid) nc = ncoil_count
     1        + nsad_group(nsmid + 1 - ncp)
           nmax = MAX(nmax, nc)

           DO j = js, nfils
              WRITE(cunit,'(1p,e22.14,2i6)') c_sad(n)/sfils, nwire, nz
              WRITE(sunit,'(/)')
              DO i = 1,nwire1
                 WRITE(cunit,'(1p,3e22.14)')
     1              x_sad(i,n,j), y_sad(i,n,j), z_sad(i,n,j)
                 IF (i .eq. nwire1) THEN
                    WRITE(munit, '(1p,4e22.14,i4,a8)')
     1              x_sad(1,n,j), y_sad(1,n,j), z_sad(1,n,j), zz, nc,
     2              ch_sad
                 ELSE
                    WRITE(munit,'(1p,4e22.14)')
     1              x_sad(i,n,j), y_sad(i,n,j), z_sad(i,n,j),
     2              c_sad(n)/sfils
                 END IF
                 WRITE(sunit,'(3f14.6,2e14.6)')
     1              x_sad(i,n,j), y_sad(i,n,j), z_sad(i,n,j),
     2              v_sad(i,n), u_sad(i,n)
              END DO
           END DO
        END DO

        IF (lspline .and. lctrlpt) THEN
           DO n = 1, nsmid
              DO j = 1, nsad_v
                 WRITE(tunit,'(1p,2e22.14)') sad_v_c(n,j), sad_u_c(n,j)
              END DO
              WRITE(tunit,'(/)')
           END DO
        END IF

      END IF

      IF (lvf) THEN
      ncoil_count = nmax

!     WRITE vf coils, currents
        DO n = 1, nvf
         WRITE(gunit,'(/)')
!        compute coil group number
         nc = ncoil_count + (n+1)/2
         nmax = MAX(nmax, nc)

         WRITE(cunit,'(1p,e22.14,2i6)') cvf(n), nwire, nz
         DO i = 1,nwire1
           WRITE(cunit,'(1p,3e22.14)')
     1        x_vf(i,1,n), y_vf(i,1,n), z_vf(i,1,n)
           IF (i .eq. nwire1) THEN
              WRITE(munit, '(1p,4e22.14,i4,a4)')
     1          x_vf(i,1,n), y_vf(i,1,n), z_vf(i,1,n), zz, nc, " VF"
           ELSE
              WRITE(munit,'(1p,4e22.14)')
     1          x_vf(i,1,n), y_vf(i,1,n), z_vf(i,1,n), cvf(n)
           END IF
           WRITE(gunit,'(3f14.6,i6,2e14.6)')
     1      x_vf(i,1,n), y_vf(i,1,n), z_vf(i,1,n), i, zz, zz
          END DO
        END DO
      END IF

      IF (ltfc) THEN
        ncoil_count = nmax

!       compute coil group number
        nc = ncoil_count + 1
        nmax = MAX(nmax, nc)

!       WRITE tf coils, currents
        DO n = 1, mtfcoil
          WRITE(gunit,'(/)')
          WRITE(cunit,'(1p,e22.14,2i6)') tfc_cur(n), mtfwire-1, nz
          DO i = 1,mtfwire
            WRITE(cunit,'(1p,3e22.14)')
     1        tfc_x(n,i), tfc_y(n,i), tfc_z(n,i)
            IF (i .eq. mtfwire) THEN
              WRITE(munit, '(1p,4e22.14,i4,a4)')
     1          tfc_x(n,i), tfc_y(n,i), tfc_z(n,i), zz, nc, " TF"
            ELSE
              WRITE(munit,'(1p,4e22.14)')
     1          tfc_x(n,i), tfc_y(n,i), tfc_z(n,i), tfc_cur(n)
            END IF
            WRITE(gunit,'(3f14.6,i6,2e14.6)')
     1       tfc_x(n,i), tfc_y(n,i), tfc_z(n,i), i, zz, zz
          END DO
        END DO
      END IF

      IF (lbcoil) THEN
      ncoil_count = nmax

!     compute coil group number
      nc = ncoil_count + 1
      nmax = MAX(nmax, nc)
!     WRITE background coils, currents
        DO n = 1, mbcoils
         mbw = mbwires(n) + 1
         WRITE(cunit,'(1p,e22.14,2i6)') bcoil_cur(n), mbwires(n),
     1      mc_bg(n)
         WRITE(gunit,'(/)')
         DO i = 1,mbw
           WRITE(cunit,'(1p,3e22.14)')
     1        bcoil_x(n,i), bcoil_y(n,i), bcoil_z(n,i)
           IF (i .eq. mbw) THEN
              WRITE(munit, '(1p,4e22.14,i4,a4)')
     1          bcoil_x(n,i), bcoil_y(n,i), bcoil_z(n,i), zz,
     2          nc + mc_bg(n), " BG"
           ELSE
              WRITE(munit,'(1p,4e22.14)')
     1          bcoil_x(n,i), bcoil_y(n,i), bcoil_z(n,i), bcoil_cur(n)
           END IF
           WRITE(gunit,'(3f14.6,i6,2e14.6)')
     1       bcoil_x(n,i), bcoil_y(n,i), bcoil_z(n,i), i, zz, zz
         END DO
        END DO
      END IF

!     Compute and write magnetic axis coordinates
      IF (laxis) THEN
         dph_w = twopi/nwire
         DO i = 1, nwire1
            ph_w = (i - 1)*dph_w
            rps = 0
            zps = 0
            DO n = 1, ntor
               rps = rps + rmnc_a(n)*cos( xn_b(n)*ph_w)
               zps = zps + zmns_a(n)*sin(-xn_b(n)*ph_w)
            END DO
            xw(i) = rps*cos(ph_w)
            yw(i) = rps*sin(ph_w)
            zw(i) = zps
         END DO
         plascur = 0
         ch_sad = ' Plasma'
         nc = nmax + 1
         DO i = 1,nwire1
            IF (i .eq. nwire1) THEN
               WRITE(munit, '(1p,4e22.14,i4,a8)')
     1            xw(1), yw(1), zw(1), zz, nc, ch_sad
            ELSE
               WRITE(munit,'(1p,4e22.14)')
     1            xw(i), yw(i), zw(i), plascur
            END IF
         END DO
      END IF

      IF (laccess) THEN
!     WRITE access zones to coils file
        DO n = 1, n_access
         WRITE(gunit,'(/)')
         DO i = 1, n_access_pts
          WRITE(gunit,'(3f14.6,i6,2e14.6)')
     1      x_access(n,i), y_access(n,i), z_access(n,i), i, zz, zz
         END DO
        END DO
      END IF

      WRITE (munit, '(a3)') "end"

      CLOSE(munit)
      CLOSE(cunit)
      CLOSE(gunit)
      CLOSE(sunit)

      END SUBROUTINE write_coils
