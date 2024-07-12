      SUBROUTINE free_mem_funct3d_par
      USE vmec_main
      USE realspace
      USE vforces
      USE vacmod

      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1 = 0
C-----------------------------------------------
      IF (ALLOCATED(parmn))
     &   DEALLOCATE(parmn, pazmn, pbrmn, pbzmn, pcrmn, pczmn, pblmn,
     &              pclmn, pr1, pru, prv, pz1, pzu, pzv, pgcon, prcon,
     &              pzcon, prcon0, pzcon0, pguu, pguv, pgvv,
     &              pru0, pzu0, stat=istat1)
      IF (istat1 .ne. 0) THEN
         STOP 'deallocation error#1 in funct3d'
      END IF

      IF (ALLOCATED(pextra1)) THEN
         DEALLOCATE (pextra1, pextra2, pextra3, pextra4, stat=istat1)
      END IF
      IF (istat1 .ne. 0) THEN
         STOP 'deallocation error#3 in funct3d'
      END IF

      END SUBROUTINE free_mem_funct3d_par

      SUBROUTINE free_mem_funct3d
      USE vmec_main
      USE realspace
      USE vforces
      USE vacmod

      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1 = 0
C-----------------------------------------------

      IF (ALLOCATED(armn)) THEN
         DEALLOCATE (armn, azmn, brmn, bzmn, crmn, czmn, blmn, clmn,
     &               r1, ru, rv, z1, zu, zv, gcon, rcon, zcon, ru0, zu0,
     &               rcon0, zcon0, guu, guv, gvv, sigma_an, stat=istat1)
      END IF
      IF (istat1 .ne. 0) THEN
         STOP 'deallocation error#1 in funct3d'
      END IF

#ifdef _ANIMEC
      IF (ALLOCATED(pperp)) THEN
         DEALLOCATE (pperp, ppar, onembc, pp1, pp2, pp3,
     &               stat=istat1)
      END IF
      IF (istat1 .ne. 0) THEN
         STOP 'deallocation error#1A in funct3d'
      END IF
#endif

      IF (ALLOCATED(brv)) THEN
         DEALLOCATE(brv, bphiv, bzv, bsqvac,
     &              bsubu_sur, bsubv_sur,
     &              bsupu_sur, bsupv_sur, stat=istat1)
      END IF
      IF (istat1 .ne. 0) THEN
         STOP 'deallocation error#2 in funct3d'
      END IF

      IF (ALLOCATED(bsqvac0)) THEN
         DEALLOCATE (bsqvac0)
      END IF

      IF (ALLOCATED(extra1)) THEN
         DEALLOCATE (extra1, extra2, extra3, extra4, stat=istat1)
      END IF
      IF (istat1 .ne. 0) THEN
         STOP 'deallocation error#3 in funct3d'
      END IF

      END SUBROUTINE free_mem_funct3d
