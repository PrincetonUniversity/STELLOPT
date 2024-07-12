      SUBROUTINE init_modular_coils (nvariables, xvariables, nfp)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE modular_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nc, i, n, nfp
      INTEGER :: nvariables, modes
      REAL(rprec) :: xvariables(*)
!-----------------------------------------------

      nc = nmod_coils_per_period

      nfper = nfp
      nmod_coils = nc*nfper
      nmid = (nc+1)/2
      nodd = MOD(nc,2)
      IF ((nodd.EQ.0) .AND. (.NOT.lsymm)) nmid = nmid + 1

      nvariables = 0
      nmod_coeffs = 0

      IF (nc .LE. 0) RETURN

!     Initialize the variables to values of unique coil parameters
!     and count the number of variables

      n = 0

!     First consider the CASE with coils on both symmetry planes at
!     phi = 0 and phi = pi/nfp (lsymm = F). This implies that the number
!     of coils per field period must be even (nodd = 0).

      IF ((nodd .EQ. 0) .AND. (.NOT.lsymm)) THEN

!     Symmetry coil at phi = 0.

        i = 1
        DO modes = 1, nf_phi
          n = n + 1
          xvariables(n) = phis(i,modes)
        END DO

        DO modes = 1, nf_rho
          n = n + 1
          xvariables(n) = rhos(i,modes)
        END DO

!     Coils 2 through nmid-1.
        DO i = 2, nmid-1
          n = n + 1
          xvariables(n) = phic(i,0)
          DO modes = 1,nf_phi
            n = n + 1
            xvariables(n) = phic(i,modes)
            n = n + 1
            xvariables(n) = phis(i,modes)
          END DO

          n = n + 1
          xvariables(n) = rhoc(i,0)
          DO modes = 1,nf_rho
            n = n + 1
            xvariables(n) = rhoc(i,modes)
            n = n + 1
            xvariables(n) = rhos(i,modes)
          END DO
        END DO

!     Symmetry coil at phi = pi/nfp.
        i = nmid
        DO modes = 1, nf_phi
          n = n + 1
          xvariables(n) = phis(i,modes)
        END DO

        DO modes = 1, nf_rho
          n = n + 1
          xvariables(n) = rhos(i,modes)
        END DO

!     Next consider the cases with a coil on phi = 0 (lsymm = T), or
!     on phi = pi/nfp (lsymm = F), but not both. Then there may be an
!     even number of coils per period (nodd = 0) or an odd number of
!     coils per period (nodd = 1).

      ELSE

        DO i = 1, nmid-nodd
          n = n + 1
          xvariables(n) = phic(i,0)
          DO modes = 1,nf_phi
            n = n + 1
            xvariables(n) = phic(i,modes)
            n = n + 1
            xvariables(n) = phis(i,modes)
          END DO

          n = n + 1
          xvariables(n) = rhoc(i,0)
          DO modes = 1,nf_rho
            n = n + 1
            xvariables(n) = rhoc(i,modes)
            n = n + 1
            xvariables(n) = rhos(i,modes)
          END DO
        END DO

        IF (nodd .EQ. 1) THEN
          i = nmid
          DO modes = 1, nf_phi
            n = n + 1
            xvariables(n) = phis(i,modes)
          END DO

          DO modes = 1, nf_rho
            n = n + 1
            xvariables(n) = rhos(i,modes)
          END DO
        END IF

      END IF  ! END IF ((nodd .eq. 0) .and. (lsymm .eqv. .false.))

      nmod_coeffs = n
      nvariables = n

      END SUBROUTINE init_modular_coils
