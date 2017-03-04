!-----------------------------------------------------------------------
!     Module:        diagno_init_coil
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This subroutine reads the coils file and
!                    initializes the coils fields.
!-----------------------------------------------------------------------
      SUBROUTINE diagno_init_coil
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE diagno_runtime
      USE read_wout_mod, ONLY: extcur_out => extcur, nextcur_out => nextcur
      USE vmec_input,  ONLY: extcur_in => extcur, read_indata_namelist,&
                             nv_in => nzeta, nfp_in => nfp, nigroup
      USE biotsavart
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, iunit, i, j, myphis, myphie, ik, ig, coil_dex
      REAL(rprec)  :: current, current_first
      LOGICAL :: lcoil_open
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      iunit = 11
      ! Read the wout file if extcur was allocated
      IF (lcoil .and. .not. lmut) THEN
         IF (ALLOCATED(extcur)) DEALLOCATE(extcur)
         IF (ALLOCATED(extcur_out) .and. (nextcur_out > 0)) THEN
            nextcur = nextcur_out
            ALLOCATE(extcur(nextcur))
            extcur(1:nextcur) = extcur_out(1:nextcur)
         ELSE
            OPEN(UNIT=iunit, FILE='input.' // TRIM(id_string), STATUS='OLD', IOSTAT=ier)
            IF (ier /= 0) CALL handle_error(FILE_OPEN_ERR,id_string,ier)
            CALL read_indata_namelist(iunit,ier)
            IF (ier /= 0) CALL handle_error(VMEC_INPUT_ERR,id_string,ier)
            CLOSE(iunit)
            DO i = 1, nigroup
               IF (ABS(extcur_in(i)) > 0) nextcur = i
            END DO
            ALLOCATE(extcur(nextcur))
            extcur = 0.0
            extcur(1:nextcur) = extcur_in(1:nextcur)
         END IF
      END IF
      
      ! Read the coils file
      CALL cleanup_biotsavart                  !Do this in case the coil has already been read
      ! Close the coils file as it may be left open
      INQUIRE(FILE=TRIM(coil_string),NUMBER=iunit,OPENED=lcoil_open)
      IF (lcoil_open) CLOSE(iunit)
      
      ! Process the coil file
      CALL parse_coils_file(TRIM(coil_string))
      
      IF (lmut) THEN
         nextcur = SIZE(coil_group)
         IF (ALLOCATED(extcur)) DEALLOCATE(extcur)
         ALLOCATE(extcur(1:nextcur))
         extcur = 1.0
      END IF
      
      
      ! Initiazlize the Current
      !nextcur = SIZE(coil_group) !SAL
      nextcur = MIN(SIZE(coil_group), SIZE(extcur))  !!SPH
      DO ik = 1, nextcur
         DO j = 1, coil_group(ik) % ncoil
            current = coil_group(ik) % coils(j) % current
            IF (j .eq. 1) current_first = current
            IF (current_first .ne. zero) coil_group(ik) % coils(j) % current = (current/current_first)*extcur(ik)
            !IF (current_first .ne. zero) coil_group(ik) % coils(j) % current = (current/current_first)*extcur_in(ik)
         END DO
      END DO
      
      
      ! Output some information
      IF (lverb) THEN
         WRITE(6,'(A)')   '----- COILS Information -----'
         WRITE(6,'(A,A)') '   FILE: ',TRIM(coil_string)
         WRITE(6,'(A,I3)')'   Coil Periodicity: ',nfp_bs
         WRITE(6,'(A,I3)')'   Current Systems: ',nextcur
         DO ik = 1, nextcur
            IF (ABS(extcur(ik)).ge. 1.0E6) THEN
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(ik)%ncoil,'  EXTCUR = ',extcur(ik)/1.0E6,' [MA]'
            ELSE IF (ABS(extcur(ik)).ge. 1.0E3) THEN
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(ik)%ncoil,'  EXTCUR = ',extcur(ik)/1.0E3,' [kA]'
            ELSE
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(ik)%ncoil,'  EXTCUR = ',extcur(ik),' [A]'
            END IF
         END DO
         CALL FLUSH(6)
      END IF 
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_init_coil
