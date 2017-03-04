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
     1   DEALLOCATE (brv, bphiv, bzv, bsqvac, stat=istat1)
      IF (istat1 .ne. 0) STOP 'deallocation error#2 in funct3d'

      IF (ALLOCATED(extra1))
     1   DEALLOCATE (extra1, extra2, extra3, extra4, stat=istat1)
      IF (istat1 .ne. 0) STOP 'deallocation error#3 in funct3d'

      END SUBROUTINE free_mem_funct3d
