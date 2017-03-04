!-----------------------------------------------------------------------
!     Subroutine:    stelltran_eqconsistency
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          09/04/2015
!     Description:   Check for equilibrium consistency
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_eqconsistency(itime,lconsistent)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime
      USE stelltran_vars
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(inout) :: lconsistent
      INTEGER, INTENT(in)    :: itime
      
      REAL(rprec) :: favg
      REAL(rprec), DIMENSION(prof_length) :: fval, f1, f2
      CHARACTER(1) :: found_char

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      lconsistent = .TRUE.
      
      ! Electron Temperature
      !f1   = te_old(:,2)
      !f2   = te(:,2)
      !fval = 0
      !WHERE(f2>0) fval = ABS(f1-f2)/ABS(f2)
      !fval(prof_length) = 0
      !IF (ANY(fval > fit_tol)) lconsistent = .FALSE.
      
      ! Ion Temperature
      !f1   = ti_old(:,2)
      !f2   = ti(:,2)
      !fval = 0
      !WHERE(f2>0) fval = ABS(f1-f2)/ABS(f2)
      !fval(prof_length) = 0
      !IF (ANY(fval > fit_tol)) lconsistent = .FALSE.
      
      ! Electron Density
      !f1   = ne_old(:,2)
      !f2   = ne(:,2)
      !fval = 0
      !WHERE(f2>0) fval = ABS(f1-f2)/ABS(f2)
      !fval(prof_length) = 0
      !IF (ANY(fval > fit_tol)) lconsistent = .FALSE.
      
      ! Z_effective
      !f1   = zeff_old(:,2)
      !f2   = zeff(:,2)
      !fval = 0
      !WHERE(f2>0) fval = ABS(f1-f2)/ABS(f2)
      !fval(prof_length) = 0
      !IF (ANY(fval > fit_tol)) lconsistent = .FALSE.
      
      ! Ohmic current
      f1   = johm_old(:,2)
      f2   = johm(:,2)
      fval = 0
      WHERE(f2/=0) fval = ABS(f1-f2)/ABS(f2)
      fval(prof_length) = 0
      favg = SUM(fval,DIM=1)/prof_length
      IF (favg > fit_tol) lconsistent = .FALSE.
      
      ! Bootstrap Current (average)
      f1   = jboot_old(:,2)
      f2   = jboot(:,2)
      fval = 0
      WHERE(f2/=0) fval = ABS(f1-f2)/ABS(f2)
      fval(prof_length) = 0
      favg = SUM(fval,DIM=1)/prof_length
      IF (favg > fit_tol) lconsistent = .FALSE.
      
      ! RF Current
      f1   = jrf_old(:,2)
      f2   = jrf(:,2)
      fval = 0
      WHERE(f2/=0) fval = ABS(f1-f2)/ABS(f2)
      fval(prof_length) = 0
      IF (ANY(fval > fit_tol)) lconsistent = .FALSE.
      
      ! Beam Current
      f1   = jbeam_old(:,2)
      f2   = jbeam(:,2)
      fval = 0
      WHERE(f2/=0) fval = ABS(f1-f2)/ABS(f2)
      fval(prof_length) = 0
      IF (ANY(fval > fit_tol)) lconsistent = .FALSE.

      found_char = ' '
      IF (lconsistent) THEN
         found_char = '*'
         !WRITE(6,'(2X,I6,A)') itime,'  Equilibrium found'
         lveryfirst_pass = .FALSE.
         CALL stelltran_runeq(itime)
         RETURN
      END IF
      
      IF (lverb.and.lfirst_pass) THEN
         WRITE(6,'(4X,I3.3,3X,8(2X,F7.5),A)')  itime,MAXVAL(te(:,2),DIM=1)*1E-3,  MAXVAL(ti(:,2),DIM=1)*1E-3,&
      	                           MAXVAL(ne(:,2),DIM=1)*1E-20,MAXVAL(zeff(:,2),DIM=1),&
      	                           MAXVAL(ABS(johm(:,2)),DIM=1)*1E-6,MAXVAL(ABS(jboot(:,2)),DIM=1)*1E-6,&
      	                           MAXVAL(ABS(jrf(:,2)),DIM=1)*1E-6,MAXVAL(ABS(jbeam(:,2)),DIM=1)*1E-6,&
                                   found_char
      ELSE
         WRITE(6,'(10X,8(2X,F7.5),A)')  MAXVAL(te(:,2),DIM=1)*1E-3,  MAXVAL(ti(:,2),DIM=1)*1E-3,&
      	                           MAXVAL(ne(:,2),DIM=1)*1E-20,MAXVAL(zeff(:,2),DIM=1),&
      	                           MAXVAL(ABS(johm(:,2)),DIM=1)*1E-6,MAXVAL(ABS(jboot(:,2)),DIM=1)*1E-6,&
      	                           MAXVAL(ABS(jrf(:,2)),DIM=1)*1E-6,MAXVAL(ABS(jbeam(:,2)),DIM=1)*1E-6,&
                                   found_char
      END IF
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE stelltran_eqconsistency
