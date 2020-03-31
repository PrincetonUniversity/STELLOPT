      SUBROUTINE fdjac_parallel(j, fcn)
      USE fdjac_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: j
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iflag
!DEC$ IF DEFINED (CRAY)
      INTEGER :: i
!DEC$ ENDIF
      REAL(rprec) :: temp, temp2, h, enorm
      EXTERNAL fcn, enorm
C-----------------------------------------------
!DEC$ IF .NOT.DEFINED (MPI_OPT)
!
!     THIS ROUTINE IS PASSED TO THE MULTI-PROCESSOR HANDLING
!     ROUTINE

      temp = xp(j)
      h = eps*ABS(temp)
      IF (h .eq. zero) h = eps
      IF( flip(j)) h = -h
      xp(j) = temp + h
      iflag = j

      CALL fcn (m, n, xp, wap, iflag, ncnt)

      temp2 = enorm(m, wap)
      WRITE(6, '(2x,i6,7x,1es12.4)') ncnt+j, temp2**2

!
!     WRITE TO A UNIQUE FILE FOR I/O IN MULTI-PROCESSOR SYSTEM
!
      WRITE (j+1000) j, iflag, h, temp2
!DEC$ IF DEFINED (CRAY)
      DO i = 1,m
         WRITE (j+1000) wap(i)
      END DO
      DO i = 1,n
         WRITE (j+1000) xp(i)
      END DO
!DEC$ ELSE
      WRITE (j+1000) wap
      WRITE (j+1000) xp
!DEC$ ENDIF

      CLOSE (j+1000)                      !!Needed to run correctly in multi-tasking...

      xp(j) = temp
!DEC$ ENDIF
      END SUBROUTINE fdjac_parallel
