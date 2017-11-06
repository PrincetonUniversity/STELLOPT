      SUBROUTINE free_mem_nunv
      USE vmec_main
      USE vacmod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1 = 0, istat2 = 0, istat3 = 0
C-----------------------------------------------

      IF (ALLOCATED(bsubu0))
     1    DEALLOCATE (bsubu0, rbsq, dbsq, stat=istat1)
#ifdef _ANIMEC
      IF (ALLOCATED(pperp_ns)) DEALLOCATE(pperp_ns)
#endif
      IF (ALLOCATED(rmn_bdy))
     1    DEALLOCATE (rmn_bdy, zmn_bdy, stat=istat2)

      IF (ALLOCATED(amatsav))
     1    DEALLOCATE (amatsav, bvecsav, potvac, bsqsav, 
     2                raxis_nestor, zaxis_nestor, stat=istat3)

      IF (istat1.ne.0 .or. istat2.ne.0 .or. istat3.ne.0) THEN
          PRINT *,' deallocation problem in free_mem_nunv'
          PRINT *,' istat1 = ',istat1,' istat2 = ',istat2
          PRINT *,' istat3 = ',istat3
      ENDIF

      END SUBROUTINE free_mem_nunv
