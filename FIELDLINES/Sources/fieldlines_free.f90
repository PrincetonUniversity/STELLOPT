!-----------------------------------------------------------------------
!     Module:        fieldlines_free
!     Authors:       S. Lazerson (samuel.lazerson@guass-fusion.com)
!     Date:          06/29/2025
!     Description:   Deallocate and free all arrays. 
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_free(IN_COMM)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE fieldlines_runtime
      USE fieldlines_grid
      USE fieldlines_lines
      USE mpi_sharmem
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier,i
      INTEGER, INTENT(INOUT), OPTIONAL :: IN_COMM
!-----------------------------------------------------------------------
!     External Functions
!          A00ADF               NAG Detection
!-----------------------------------------------------------------------
!      EXTERNAL A00ADF
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ier = 0
      IF (PRESENT(IN_COMM)) THEN
         IF (ASSOCIATED(raxis))    CALL mpidealloc(raxis,win_raxis)
         IF (ASSOCIATED(phiaxis))  CALL mpidealloc(phiaxis,win_phiaxis)
         IF (ASSOCIATED(zaxis))    CALL mpidealloc(zaxis,win_zaxis)
         IF (ASSOCIATED(B_R))      CALL mpidealloc(B_R,win_B_R)
         IF (ASSOCIATED(B_PHI))    CALL mpidealloc(B_PHI,win_B_PHI)
         IF (ASSOCIATED(B_Z))      CALL mpidealloc(B_Z,win_B_Z)
         IF (ASSOCIATED(PRES_G))      CALL mpidealloc(PRES_G,win_PRES)
         IF (ASSOCIATED(MU3D))      CALL mpidealloc(MU3D,win_MU)
         IF (ASSOCIATED(MODB4D))      CALL mpidealloc(MODB4D,win_MODB4D)
         IF (ASSOCIATED(MU4D))      CALL mpidealloc(MU4D,win_MU4D)
         IF (ASSOCIATED(BR4D))     CALL mpidealloc(BR4D,win_BR4D)
         IF (ASSOCIATED(BZ4D))     CALL mpidealloc(BZ4D,win_BZ4D)
      ELSE
         IF (ASSOCIATED(raxis))    DEALLOCATE(raxis)
         IF (ASSOCIATED(phiaxis))  DEALLOCATE(phiaxis)
         IF (ASSOCIATED(zaxis))    DEALLOCATE(zaxis)
         IF (ASSOCIATED(B_R))      DEALLOCATE(B_R)
         IF (ASSOCIATED(B_PHI))    DEALLOCATE(B_PHI)
         IF (ASSOCIATED(B_Z))      DEALLOCATE(B_Z)
         IF (ASSOCIATED(PRES_G))      DEALLOCATE(PRES_G)
         IF (ASSOCIATED(MU3D))      DEALLOCATE(MU3D)
         IF (ASSOCIATED(MODB4D))      DEALLOCATE(MODB4D)
         IF (ASSOCIATED(MU4D))      DEALLOCATE(MU4D)
         IF (ASSOCIATED(BR4D))     DEALLOCATE(BR4D)
         IF (ASSOCIATED(BZ4D))     DEALLOCATE(BZ4D) 
      ENDIF
      IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
      IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
      IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
      IF (ALLOCATED(Rhc_lines)) DEALLOCATE(Rhc_lines)
      IF (ALLOCATED(Zhc_lines)) DEALLOCATE(Zhc_lines)
      IF (ALLOCATED(B_lines)) DEALLOCATE(B_lines)
      IF (EZspline_allocated(BR_spl)) CALL EZspline_free(BR_spl,ier)
      IF (EZspline_allocated(BZ_spl)) CALL EZspline_free(BZ_spl,ier)
      IF (EZspline_allocated(MU_spl)) CALL EZspline_free(MU_spl,ier)
      IF (EZspline_allocated(MODB_spl)) CALL EZspline_free(MODB_spl,ier)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_free