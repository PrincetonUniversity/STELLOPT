!-----------------------------------------------------------------------
!     Subroutine:    stelltran_read_database_predictive
!     Authors:       J. Mittelstaedt (jmittelstaedt@uchicago.edu)
!     Date:          07/27/2016
!     Description:   This routine sets up the arrays of timeslice
!                    profile info and loads the first time slice.
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_read_database_predictive
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE safe_open_mod
      USE stelltran_runtime
      USE stelltran_data
!-----------------------------------------------------------------------
!     Local Variables
!          lexist         Logical flag for file existance
!          ier            Local Error flag
!          iunit          Iunit for database file
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL :: lexist
      INTEGER :: ier, iunit,i, ntemp, ik
      CHARACTER(256) :: temp_str
      
      LOGICAL :: ldebug = .FALSE.
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      
      ! Check for database file
      INQUIRE(FILE=TRIM(db_file),EXIST=lexist)
      IF (.not.lexist) CALL handle_err(FILE_EXIST_ERR,TRIM(db_file),ier)
      
      ! Open the database file
      CALL safe_open(iunit,ier,TRIM(db_file),'old','formatted')
      READ(iunit,*) ntemp ! ASSUMING 3 TIMESTEPS NOW
      
      DO
         READ(iunit,'(A2,I3)',IOSTAT=ier) temp_str,ntemp
         IF (ier < 0) EXIT ! End of file reached
         SELECT CASE(TRIM(temp_str(1:2)))
            CASE('IP')
               IF (ALLOCATED(Ip)) DEALLOCATE(Ip)
               ALLOCATE(Ip(ntimesteps))
               READ(iunit,*) Ip(1:ntimesteps)
            CASE('VL')
               IF (ALLOCATED(Vloop)) DEALLOCATE(Vloop)
               ALLOCATE(Vloop(ntimesteps))
               READ(iunit,*) Vloop(1:ntimesteps)
            CASE('PE')
               necrh = ntemp
               IF (ALLOCATED(P_ecrh)) DEALLOCATE(P_ecrh)
               ALLOCATE(P_ecrh(necrh,ntimesteps))
               DO ik = 1, ntimesteps
                  READ(iunit,*) P_ecrh(1:necrh,ik)
               END DO
            CASE('NE')
               nne = ntemp
               IF (ALLOCATED(ne_s)) DEALLOCATE(ne_s)
               IF (ALLOCATED(ne_f)) DEALLOCATE(ne_f)
               ALLOCATE(ne_s(nne),ne_f(nne,ntimesteps+1))
               READ(iunit,*) ne_s(1:nne)
               READ(iunit,*) ne_f(1:nne,1)
            CASE('TE')
               nte = ntemp
               IF (ALLOCATED(te_s)) DEALLOCATE(te_s)
               IF (ALLOCATED(te_f)) DEALLOCATE(te_f)
               ALLOCATE(te_s(nte),te_f(nte,ntimesteps+1))
               READ(iunit,*) te_s(1:nte)
               READ(iunit,*) te_f(1:nte,1)
            CASE('TI')
               nti = ntemp
               IF (ALLOCATED(ti_s)) DEALLOCATE(ti_s)
               IF (ALLOCATED(ti_f)) DEALLOCATE(ti_f)
               ALLOCATE(ti_s(nti),ti_f(nti,ntimesteps+1))
               READ(iunit,*) ti_s(1:nti)
               READ(iunit,*) ti_f(1:nti,1)
            CASE('ZE')
               nzeff = ntemp
               IF (ALLOCATED(zeff_s)) DEALLOCATE(zeff_s)
               IF (ALLOCATED(zeff_f)) DEALLOCATE(zeff_f)
               ALLOCATE(zeff_s(nzeff),zeff_f(nzeff,ntimesteps+1))
               READ(iunit,*) zeff_s(1:nzeff)
               READ(iunit,*) zeff_f(1:nzeff,1)
         END SELECT
      END DO
      
      ! Close the database file
      CLOSE(iunit)
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE stelltran_read_database_predictive
