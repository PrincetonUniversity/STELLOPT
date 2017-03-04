      SUBROUTINE load_modular_structures
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE boundary, ONLY: nfp
      USE modular_coils
      USE tf_coils
      USE Vcoilpts
      USE Vwire
      USE coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, modes
!-----------------------------------------------

!     load the unique coil parameters with values from variables in
!     optimization

!     First consider the CASE with coils on both symmetry planes at
!     phi = 0 and phi = pi/nfp (lsymm = F). This implies that the number
!     of coils per field period must be even (nodd = 0).

      IF ((nodd.EQ.0) .AND. .NOT.lsymm) THEN

!        Symmetry coil at phi = 0.
         i = 1
         DO modes = 1, nf_phi
            modular(i)%phis(modes) = phis(i,modes)
         END DO

         rhos(i,0) = 0
         DO modes = 1, nf_rho
            modular(i)%rhos(modes) = rhos(i,modes)
         END DO

!        Coils 2 through nmid - 1
         DO i = 2, nmid-1
            phis(i,0) = 0
            DO modes = 0,nf_phi
               modular(i)%phic(modes) = phic(i,modes)
               modular(i)%phis(modes) = phis(i,modes)
            END DO

            rhos(i,0) = 0
            DO modes = 0,nf_rho
               modular(i)%rhoc(modes) = rhoc(i,modes)
               modular(i)%rhos(modes) = rhos(i,modes)
            END DO
         END DO

!        Symmetry coil at phi = pi/nfp
         i = nmid
         phis(i,0) = 0
         DO modes = 0, nf_phi
            modular(i)%phis(modes) = phis(i,modes)
         END DO

         rhos(i,0) = 0
         DO modes = 0, nf_rho
            modular(i)%rhos(modes) = rhos(i,modes)
         END DO

         DO i = 1, nmid
            modular(i)%current = curmod(i)
         END DO

!     Next consider the cases with a coil on phi = 0 (lsymm = T), or
!     on phi = pi/nfp (lsymm = F), but not both. THEN there may be an
!     even number of coils per period (nodd = 0) or an odd number of
!     coils per period (nodd = 1).

      ELSE

         DO i = 1, nmid-nodd
            phis(i,0) = 0
            DO modes = 0,nf_phi
               modular(i)%phic(modes) = phic(i,modes)
               modular(i)%phis(modes) = phis(i,modes)
            END DO

            rhos(i,0) = 0
            DO modes = 0,nf_rho
               modular(i)%rhoc(modes) = rhoc(i,modes)
               modular(i)%rhos(modes) = rhos(i,modes)
            END DO
         END DO
         IF (nodd .eq. 1) THEN
            i = nmid
            phis(i,0) = 0
            DO modes = 0, nf_phi
               modular(i)%phis(modes) = phis(i,modes)
            END DO

            rhos(i,0) = 0
            DO modes = 0, nf_rho
               modular(i)%rhos(modes) = rhos(i,modes)
            END DO
         END IF

         DO i = 1, nmid
            modular(i)%current = curmod(i)
         END DO

      ENDIF   ! END IF ((nodd .eq. 0) .and. (lsymm .eqv. .false.))

!     MOD current scale factor for coil spacing, curvature penalties

      DO i = 1, nmid
         cmod_scl(i) = ABS(i_pol)/(one + ABS(curmod(i)))
      END DO

      END SUBROUTINE load_modular_structures
