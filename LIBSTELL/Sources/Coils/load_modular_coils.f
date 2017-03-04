      SUBROUTINE load_modular_coils (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE modular_coils
      USE tf_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, n
      INTEGER :: nvariables, modes
      REAL(rprec) :: xvariables(*)
!-----------------------------------------------

!     load the unique coil parameters with values from variables in
!     optimization

      nvariables = 0
      n = 0

!     Case with coils on both symmetry planes phi = 0 and phi = pi/nfp
!     (lsymm = F). No. of coils per field period must be even (nodd = 0).

      IF ((nodd.EQ.0) .AND. (.NOT.lsymm)) THEN

!        Symmetry coil at phi = 0.
         i = 1
         phis(i,0) = 0
         DO modes = 1, nf_phi
            n = n + 1
            phis(i,modes) = xvariables(n)
         END DO

         rhos(i,0) = 0
         DO modes = 1, nf_rho
            n = n + 1
            rhos(i,modes) = xvariables(n)
         END DO

!        Coils 2 through nmid - 1
         DO i = 2, nmid-1
            n = n + 1
            phic(i,modes) = xvariables(n)
            phis(i,0) = 0
            DO modes = 1,nf_phi
               n = n + 1
               phic(i,modes) = xvariables(n)
               n = n + 1
               phis(i,modes) = xvariables(n)
            END DO

            n = n + 1
            rhoc(i,0) = xvariables(n)
            rhos(i,0) = 0
            DO modes = 1,nf_rho
               n = n + 1
               rhos(i,modes) = xvariables(n)
            END DO
         END DO

!        Symmetry coil at phi = pi/nfp
         i = nmid
         phis(i,0) = 0
         DO modes = 1, nf_phi
            n = n + 1
            phis(i,modes) = xvariables(n)
         END DO

         rhos(i,0) = 0
         DO modes = 1, nf_rho
            n = n + 1
            rhos(i,modes) = xvariables(n)
         END DO

!     Cases with a coil on phi = 0 (lsymm = T), or on phi = pi/nfp
!     (lsymm = F), but not both. There may be an even number of coils
!     per period (nodd = 0), or an odd number of coils per period
!     (nodd = 1).

      ELSE

         DO i = 1, nmid-nodd
            n = n + 1
            phic(i,0) = xvariables(n)
            phis(i,0) = 0
            DO modes = 1,nf_phi
               n = n + 1
               phic(i,modes) = xvariables(n)
               n = n + 1
               phis(i,modes) = xvariables(n)
            END DO

            n = n + 1
            rhoc(i,0) = xvariables(n)
            rhos(i,0) = 0
            DO modes = 1,nf_rho
               n = n + 1
               rhoc(i,modes) = xvariables(n)
               n = n + 1
               rhos(i,modes) = xvariables(n)
            END DO
         END DO
         IF (nodd .EQ. 1) THEN
            i = nmid
            phis(i,0) = 0
            DO modes = 1, nf_phi
               n = n + 1
               phis(i,modes) = xvariables(n)
            END DO

            rhos(i,0) = 0
            DO modes = 1, nf_rho
               n = n + 1
               rhos(i,modes) = xvariables(n)
            END DO
         END IF

      END IF

      nvariables = n

      END SUBROUTINE load_modular_coils
