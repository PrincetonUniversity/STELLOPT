!Calculate wallclock times on each routine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

!----------------------------------------------------------------------------------------------- 
!Calculate the total time ttotal spent at routine, number of iterations ntotal and average time
!(ignoring the first iteration)
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input  
  CHARACTER*30 routine
  REAL*8 tstart
  !Input/output
  INTEGER ntotal
  REAL*8 ttotal,t0
  !Others
  REAL*8 tfinish

  IF(.NOT.TIME) RETURN
  CALL CPU_TIME(tfinish)
  ttotal=ttotal+(tfinish-tstart)
  ntotal=ntotal+1
  IF(ntotal.EQ.1) t0=ttotal
  IF(.NOT.FAST_IONS) WRITE(iout,'("Time in routine ",a,". Total: ",f10.3," Average: ",f10.3,". Iterations: ",I6)') &
       & routine,ttotal,(ttotal-t0)/(ntotal-1),ntotal
  CALL FLUSH(iout)

END SUBROUTINE CALCULATE_TIME


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
