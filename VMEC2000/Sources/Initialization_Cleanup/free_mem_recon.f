      SUBROUTINE free_mem_recon
      USE vsvd
      USE vspline
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1 = 0, istat2 = 0, istat3 = 0
C-----------------------------------------------
      IF (ALLOCATED(nk_ia))
     1  DEALLOCATE (nk_ia,nk_ib, nk_pa, nk_pb,
     2      indexs2,indexu2, indexs1, indexu1,
     3      isortr, isorts, stat=istat1)
      IF (ALLOCATED(hthom))
     1  DEALLOCATE( hthom, ythom, y2thom, pknots,
     2  hstark,y2stark,ystark, sknots, stat=istat2 )
      IF (ALLOCATED(sthom))
     1  DEALLOCATE( sthom, delse2, delso2, pcalc,
     2      delse1, delso1, starkcal, qmeas, qcalc,
     3      fpsical, stark_weight, rsort,rsort0, stat=istat3)
      IF ((istat1 .ne. 0) .or. (istat2 .ne. 0)
     1                    .or. (istat3 .ne. 0) )then
          PRINT *,' in free_mem_recon, istat1 = ',istat1
          PRINT *,' istat2 = ',istat2,' istat3 = ',istat3
        END IF

      END SUBROUTINE free_mem_recon
