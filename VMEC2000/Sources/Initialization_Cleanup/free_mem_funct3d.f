      SUBROUTINE free_mem_funct3d_par
      USE vmec_main
      USE realspace
      USE vforces
      USE vacmod
#if defined(SKS)
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1 = 0
C-----------------------------------------------
      IF (ALLOCATED(parmn))
     1   DEALLOCATE (parmn, pazmn, pbrmn, pbzmn, pcrmn, pczmn, pblmn,
     2       pclmn, pr1, pru, prv, pz1, pzu, pzv, pgcon, prcon,
     3       pzcon, prcon0, pzcon0, pguu, pguv, pgvv, 
     4       pru0, pzu0, stat=istat1)
      IF (istat1 .ne. 0) STOP 'deallocation error#1 in funct3d'

      IF (ALLOCATED(pextra1))
     1   DEALLOCATE (pextra1, pextra2, pextra3, pextra4, stat=istat1)
      IF (istat1 .ne. 0) STOP 'deallocation error#3 in funct3d'
#endif
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

      IF (ALLOCATED(armn))
     1   DEALLOCATE (armn, azmn, brmn, bzmn, crmn, czmn, blmn, clmn,
     2   r1, ru, rv, z1, zu, zv, gcon, rcon, zcon, ru0, zu0,
     3   rcon0, zcon0, guu, guv, gvv, sigma_an, stat=istat1)
      IF (istat1 .ne. 0) STOP 'deallocation error#1 in funct3d'

#ifdef _ANIMEC
      IF (ALLOCATED(pperp))
     1   DEALLOCATE (pperp, ppar, onembc, pp1, pp2, pp3,
     2               stat=istat1)
      IF (istat1 .ne. 0) STOP 'deallocation error#1A in funct3d'
#endif

      IF (ALLOCATED(brv))
     1   DEALLOCATE (brv, bphiv, bzv, bsqvac, 
     2               bsubu_sur, bsubv_sur,
     3               bsupu_sur, bsupv_sur, stat=istat1)
      IF (istat1 .ne. 0) STOP 'deallocation error#2 in funct3d'

      IF (ALLOCATED(bsqvac0)) DEALLOCATE (bsqvac0)

      IF (ALLOCATED(extra1))
     1   DEALLOCATE (extra1, extra2, extra3, extra4, stat=istat1)
      IF (istat1 .ne. 0) STOP 'deallocation error#3 in funct3d'

      END SUBROUTINE free_mem_funct3d
