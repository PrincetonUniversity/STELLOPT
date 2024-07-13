      SUBROUTINE eval_poloidal_currents
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE boundary
      USE modular_coils
      USE saddle_coils
      USE bcoils_mod
      USE tf_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i
      REAL (rprec) :: SUM0, csum0
!-----------------------------------------------

!     Evaluate the SUM0 of the poloidal currents

      csum0 = 0

      IF (lmodular) THEN
         SUM0 = 0
         IF ((nodd .eq. 0) .and. (.not. lsymm)) THEN
            DO i = 2, nmid - 1
               SUM0 = SUM0 + curmod(i)
            END DO
            SUM0 = 2*sum0 + curmod(1) + curmod(nmid)
         ELSE
            IF (nodd .eq. 0) THEN
               DO i = 1, nmid
                  SUM0 = SUM0 + curmod(i)
               END DO
               SUM0 = 2*sum0
            END IF
            IF (lsymm) THEN
               DO i = 1, nmid - 1
                  SUM0 = SUM0 + curmod(i)
               END DO
               SUM0 = 2*sum0 + curmod(nmid)
            ELSE
               DO i = 2, nmid
                  SUM0 = SUM0 + curmod(i)
               END DO
               SUM0 = 2*sum0 + curmod(1)
            END IF
         END IF
         csum0 = csum0 + nfp*sum0
      END IF

      IF (lsaddle .and. lsmod) THEN
         SUM0 = 0
         DO i = 1, nsmid
            SUM0 = SUM0 + cursad(nsad_group(i))*csad_scl(i)
         END DO
         csum0 = csum0 + 2*nfp*sum0
      END IF

      IF (lbcoil) THEN
         SUM0 = 0
         DO i = 1, mbcoils
            IF (lp_bg(i)) SUM0 = SUM0 + bcoil_cur(i)
         END DO
         csum0 = csum0 + SUM0
      END IF

      IF (ltfc) THEN
         csum0 = csum0 + i_tfc
      END IF

      pol_cur = csum0

      END SUBROUTINE eval_poloidal_currents
