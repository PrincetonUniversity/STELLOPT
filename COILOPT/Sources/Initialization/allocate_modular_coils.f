      SUBROUTINE allocate_modular_coils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE boundary, ONLY: nfp
      USE modular_coils
      USE Vcoilpts
      USE Vwire
      USE coils
      USE mpi_params
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nc, i, status, modes
      LOGICAL :: lmaster
!-----------------------------------------------
      lmaster = (myid .EQ. master)
!     ALLOCATE modular coil arrays

      nc = nmod_coils_per_period
      IF (lmaster .AND. nc.LE.0) 
     1  STOP 'nmod_coils_per_period must > 0!'

      ALLOCATE (modular(nc), stat=status)
      ALLOCATE (x_mod(nwdim1,3,ncdim), y_mod(nwdim1,3,ncdim), 
     1          z_mod(nwdim1,3,ncdim), rho(ncdim,nwdim1) , 
     2          phi(ncdim,nwdim1), rcoil(ncdim,nwdim), 
     3          zcoil(ncdim,nwdim), stat=status)
      
      IF (lmaster .AND. status.NE.0)
     1   STOP "Cannot ALLOCATE modular_coils"

      nfper = nfp
      nmod_coils = nc * nfper
      IF (lmaster .AND. nmod_coils.GT.ncdim)
     1   STOP 'nmod_coils > ncdim'

      nmid = (nc+1)/2
      nodd = MOD(nc,2)
      IF (lmaster) THEN
        PRINT *, 'no. modular coils (per period) = ', nc
      END IF

      IF (nmod_coils .LE. 0) RETURN

!     Initialize arrays to values of unique coil parameters

!     First consider the CASE with coils on both symmetry planes at
!     phi = 0 and phi = pi/nfp (lsymm = F). This implies that the number
!     of coils per field period must be even (nodd = 0).

      IF ((nodd.EQ.0) .AND. (.NOT.lsymm)) THEN
        nmid = nmid + 1

        IF (lmaster) THEN
          PRINT *, 'nmid = ', nmid, 'nodd = ', nodd
        END IF

!     Allocate only those coil components which are needed. There
!     is no need to store the Fourier coefficients for coils whose
!     locations in real space will be computed from stellarator
!     symmetry.

!     ALLOCATE coil components for symmetry coil at phi = 0.
        i = 1
        ALLOCATE(modular(i)%rhoc(0:nf_rho), modular(i)%rhos(0:nf_rho),
     1           modular(i)%phic(0:nf_phi), modular(i)%phis(0:nf_phi),
     2           modular(i)%phi(nwire1),    modular(i)%rho(nwire1),
     3           stat=status)
        IF (lmaster .AND. status.NE.0)
     1     STOP "Cannot ALLOCATE center coil components"

!     ALLOCATE coil components for coils 2 through nmid-1.
        DO i = 2, nmid-1
          ALLOCATE(modular(i)%rhoc(0:nf_rho), modular(i)%rhos(0:nf_rho),
     1             modular(i)%phic(0:nf_phi), modular(i)%phis(0:nf_phi),
     2             modular(i)%phi(nwire1),    modular(i)%rho(nwire1),
     3             stat = status)
        IF (lmaster .AND. status .NE. 0) 
     1     STOP "Cannot ALLOCATE modular components"

        END DO               !! DO i = ... loop over half of coils

!     ALLOCATE coil components for symmetry coil at phi = pi/nfp.
        i = nmid

        ALLOCATE(modular(i)%rhoc(0:nf_rho), modular(i)%rhos(0:nf_rho),
     1           modular(i)%phic(0:nf_phi), modular(i)%phis(0:nf_phi),
     2           modular(i)%phi(nwire1),    modular(i)%rho(nwire1),
     3           stat=status)
        IF (lmaster .AND. status.NE.0)
     1     STOP "Cannot ALLOCATE center coil components"

!     ALLOCATE coil components phi, rhi for coils nmid + 1 through nc.
        DO i = 1, nmid-2

          ALLOCATE(modular(nc-i+1)%phi(nwire1),
     1             modular(nc-i+1)%rho(nwire1),
     2             stat=status)
          IF (lmaster .AND. status.NE.0)
     1       STOP "Cannot ALLOCATE modular components"

        END DO               !! DO i = ... loop over half of coils

!     Symmetry coil at phi = 0.
        i = 1
        modes = 0
        phic(i,modes) = 0
        modular(i)%phic(modes) = 0
        modular(i)%phis(modes) = phis(i,modes)
        DO modes = 1, nf_phi
          phic(i,modes) = 0
          modular(i)%phic(modes) = 0
          modular(i)%phis(modes) = phis(i,modes)
        END DO

        modes = 0
        DO modes = 0, nf_rho
          modular(i)%rhos(modes) = rhos(i,modes)
          modular(i)%rhoc(modes) = rhoc(i,modes)
        END DO
        modular(i)%rhos(modes) = 0

!     Coils 2 through nmid-1.
        DO i = 2, nmid-1
          modes = 0
          modular(i)%phic(modes) = phic(i,modes)
          modular(i)%phis(modes) = 0
          DO modes = 1,nf_phi
            modular(i)%phic(modes) = phic(i,modes)
            modular(i)%phis(modes) = phis(i,modes)
          END DO

          modes = 0
          DO modes = 0,nf_rho
            modular(i)%rhos(modes) = rhos(i,modes)
            modular(i)%rhoc(modes) = rhoc(i,modes)
          END DO
          modular(i)%rhos(modes) = 0
        END DO

!     Symmetry coil at phi = pi/nfp.
        i = nmid
        modes = 0
        phic(i,modes) = 0
        modular(i)%phic(modes) = 0
        modular(i)%phis(modes) = phis(i,modes)
        DO modes = 1, nf_phi
          phic(i,modes) = 0
          modular(i)%phic(modes) = 0
          modular(i)%phis(modes) = phis(i,modes)
        END DO

        modes = 0
        DO modes = 1, nf_rho
          modular(i)%rhos(modes) = rhos(i,modes)
          modular(i)%rhoc(modes) = rhoc(i,modes)
        END DO
        modular(i)%rhos(modes) = 0

!     Set currents from input

        DO i = 1, nmid
          modular(i)%current = curmod(i)
        END DO

!     Impose symmetry on coil-coil penalty weights, exponents

        DO i = 1, nmid-2
          dcc_wgt(nc+1-i) = dcc_wgt(i+1)
          dcc_exp(nc+1-i) = dcc_exp(i+1)
          dcc_tgt(nc+1-i) = dcc_tgt(i+1)
          rc_wgt(nc+1-i) = rc_wgt(i+1)
          rc_exp(nc+1-i) = rc_exp(i+1)
          rc_tgt(nc+1-i) = rc_tgt(i+1)
          r_ext(nc+1-i) = r_ext(i+1)
        END DO

!     Next consider the cases with a coil on phi = 0 (lsymm = t), or
!     on phi = pi/nfp (lsymm = F), but not both. Then there may be an
!     even number of coils per period (nodd = 0) or an odd number of
!     coils per period (nodd = 1).

      ELSE

        IF (lmaster) THEN
          PRINT *, 'nmid = ', nmid
          PRINT *, 'nodd = ', nodd
        END IF

        DO i = 1,nmid-nodd

!     Allocate only those coil components which are needed. There
!     is no need to store the Fourier coefficients for coils whose
!     locations in real space will be computed from stellarator
!     symmetry.

          ALLOCATE(modular(i)%rhoc(0:nf_rho), modular(i)%rhos(0:nf_rho),
     1             modular(i)%phic(0:nf_phi), modular(i)%phis(0:nf_phi),
     2             modular(i)%phi(nwire1),    modular(i)%rho(nwire1),
     3             modular(nmod_coils_per_period-i+1)%phi(nwire1),
     4             modular(nmod_coils_per_period-i+1)%rho(nwire1),
     5             stat=status)
          IF (lmaster .AND. status.NE.0)
     1       STOP "Cannot ALLOCATE modular components"

        END DO               !! DO i = ... loop over half of coils

        IF (nodd .eq. 1) THEN
          i = nmid

          ALLOCATE(modular(i)%rhoc(0:nf_rho), modular(i)%rhos(0:nf_rho),
     1             modular(i)%phic(0:nf_phi), modular(i)%phis(0:nf_phi),
     2             modular(i)%phi(nwire1),    modular(i)%rho(nwire1),
     3             stat = status)
          IF (lmaster .AND. status.NE.0)
     1       STOP "Cannot ALLOCATE center coil components"
        ENDIF

        DO i = 1, nmid-nodd
          modes = 0
          modular(i)%phic(modes) = phic(i,modes)
          modular(i)%phis(modes) = 0
          DO modes = 1,nf_phi
            modular(i)%phic(modes) = phic(i,modes)
            modular(i)%phis(modes) = phis(i,modes)
          END DO

          modes = 0
          modular(i)%rhos(modes) = 0
          DO modes = 1,nf_rho
            modular(i)%rhos(modes) = rhos(i,modes)
          END DO
        END DO

        IF (nodd .EQ. 1) THEN
          i = nmid
          modes = 0
          phic(i,modes) = 0
          modular(i)%phic(modes) = 0
          modular(i)%phis(modes) = phis(i,modes)
          DO modes = 1, nf_phi
            phic(i,modes) = 0
            modular(i)%phic(modes) = 0
            modular(i)%phis(modes) = phis(i,modes)
          END DO

          modes = 0
          modular(i)%rhos(modes) = 0
          DO modes = 1, nf_rho
            modular(i)%rhos(modes) = rhos(i,modes)
          END DO
        END IF

!     Set currents from input

        DO i = 1, nmid
          modular(i)%current = curmod(i)
        END DO

!     Impose symmetry on coil-coil penalty weights, exponents

        DO i = 1,nmid-nodd
          dcc_wgt(nc+1-i) = dcc_wgt(i)
          dcc_exp(nc+1-i) = dcc_exp(i)
          dcc_tgt(nc+1-i) = dcc_tgt(i)
          rc_wgt(nc+1-i) = rc_wgt(i)
          rc_exp(nc+1-i) = rc_exp(i)
          rc_tgt(nc+1-i) = rc_tgt(i)
          r_ext(nc+1-i) = r_ext(i)
        END DO

      END IF  ! END IF ((nodd .eq. 0) .and. (lsymm .eqv. .false.))

      nmod_unique_coils = nmid

      END SUBROUTINE allocate_modular_coils
