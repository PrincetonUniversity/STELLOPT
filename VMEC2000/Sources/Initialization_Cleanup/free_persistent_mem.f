      SUBROUTINE free_persistent_mem
      USE vmec_main
      USE xstuff
      USE mgrid_mod, ONLY: free_mgrid
      USE parallel_vmec_module, ONLY: RUNVMEC_COMM_WORLD
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1 = 0, istat2 = 0
c-----------------------------------------------
      IF (ALLOCATED(xc)) DEALLOCATE (xc, scalxc, stat=istat1)
      CALL free_mgrid (istat2, RUNVMEC_COMM_WORLD)

      IF (istat1.ne.0 .or. istat2.ne.0) THEN
          PRINT *,'problem in free_persistent_mem'
          PRINT *,' istat1 = ',istat1,' istat2 = ',istat2
      ENDIF

      END SUBROUTINE free_persistent_mem
