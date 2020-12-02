!-----------------------------------------------------------------------
!     Subroutine:    stellopt_stella
!     Authors:       J.L. Velasco (joseluis.velasco@ciemat.es),
!                    building on stellopt_knosos
!     Date:          last modified, 29/10/2020
!     Description:   This subroutine calculates turbulent flux
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_stella(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
        USE stellopt_runtime, ONLY:  proc_string, bigno
        USE equil_utils, ONLY: rho
        USE stellopt_targets, ONLY: sigma_stella_q1, sigma_stella_q2, lbooz, nsd
!DEC$ IF DEFINED (STELLA_OPT)
        USE stella_stellopt_mod
!DEC$ ENDIF
        USE safe_open_mod
        USE mpi_params
        USE mpi_inc
            
!-----------------------------------------------------------------------
!     Subroutine Parameters
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Parameters
!----------------------------------------------------------------------      
      INTEGER, PARAMETER :: nbx=11
      INTEGER, PARAMETER :: nsx=1000
      INTEGER nbb,ns,REGB(nbx)
      REAL*8 Ab(nbx),Zb(nbx),s(nsx),Zeff
      
!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      INTEGER ik, jk, istat, w_u3, ierr
      CHARACTER(120) :: out_file
      CHARACTER :: dkes_input_file*64, temp_str*64
      REAL*8, ALLOCATABLE :: s_st(:)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      IF (iflag < 0) RETURN
!DEC$ IF DEFINED (STELLA_OPT)
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------    stella CALCULATION     -------------------------'
      ! Enter the main loop
      IF (ALLOCATED(STELLA_Q1)) DEALLOCATE(STELLA_Q1)
      IF (ALLOCATED(STELLA_Q2)) DEALLOCATE(STELLA_Q2)
      ALLOCATE(STELLA_Q1(nsd),STELLA_Q2(nsd))
      STELLA_Q1=0.0; STELLA_Q2=0.0
      ST_Q1=0.0; ST_Q2=0.0

      ns=-1
      DO ik=2,nsd
         IF(.not. lbooz(ik)) CYCLE
         IF(sigma_stella_q1(ik) >= bigno .and. sigma_stella_q2) CYCLE
         ns=ns+1
         IF(myworkid == ns ) WRITE(temp_str,'(i3.3)') ik
      END DO
      ST_EXT='st_log.'//TRIM(proc_string)//'_s'//TRIM(temp_str)      
      !Read input files (simulation parameters, models, flux-surfaces, species...)
      IF(ANY(sigma_stella_q1 < bigno)) ST_STELLOPT(1)=.TRUE.
      IF(ANY(sigma_stella_q2 < bigno)) ST_STELLOPT(2)=.TRUE.
!      CALL READ_INPUT(ns,s,nbb,Zb,Ab,regb,Zeff)      
      !Open files
      out_file ='st_out.'//TRIM(proc_string)
      w_u3 =222
      istat=0
      IF (myworkid == master) call safe_open(w_u3,istat,out_file,'replace','formatted')
      IF (istat .ne. 0) THEN
         IF (myworkid == master) PRINT *,istat,out_file
         IF (myworkid == master) WRITE(6,*) 'Error opening STELLA output file:',TRIM(out_file),istat
         iflag=-1
         RETURN
      END IF
      ! Loop over magnetic surfaces
      IF (ALLOCATED(s_st))       DEALLOCATE(s_st);        ALLOCATE(s_st(nsd))
      IF(myworkid == master) s_st=rho ; CALL MPI_BCAST(s_st,nsd,MPI_REAL8,master,MPI_COMM_MYWORLD,ierr_mpi)
      
      jk=0
      DO ik=2,nsd
         IF(.not. lbooz(ik)) CYCLE
         IF(sigma_stella_q1(ik) >= bigno .and. sigma_stella_q2(ik) >= bigno) CYCLE
         jk=jk+1
         IF(myworkid+1 /= jk ) CYCLE
!         ST_DKESFILE='boozmn_'//TRIM(proc_string)//'.nc'
!         CALL READ_BFIELD(s_st(ik))
         IF(sigma_stella_q1(ik) >= bigno .and. sigma_stella_q2(ik) >= bigno) CYCLE
!         CALL CALC_DATABASE(ik,s_st(ik))
         STELLA_Q1(ik)=ST_Q1
         STELLA_Q2(ik)=ST_Q2
      END DO
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,STELLA_Q1,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,STELLA_Q2,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF      
      IF (myworkid == master) THEN
         DO ik=2,nsd
            IF(.not. lbooz(ik)) CYCLE
            IF(sigma_stella_q1(ik) >= bigno .and. sigma_stella_q2(ik)) CYCLE
            WRITE(w_u3,'(1(1x,i8),6(1x,e17.10))') ik,STELLA_Q1(ik),STELLA_Q2(ik)
         END DO
         IF (lscreen) WRITE(6,'(2X,I8,1(2X,E17.10))') ik,STELLA_Q1(ik),STELLA_Q2(ik)
         CALL FLUSH(6)
         CLOSE(w_u3)
      END IF
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STELLA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (lscreen) WRITE(6,'(a)') ' -----------------  stella CALCULATION (DONE) ----------------'
!DEC$ ELSE
      IF (lscreen) WRITE(6,'(a)') ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      IF (lscreen) WRITE(6,'(a)') ' !! stella based optimization not supported on your machine !!'
      IF (lscreen) WRITE(6,'(a)') ' !! Check for stella directory and STELLA_OPT in your makefile!!'
      IF (lscreen) WRITE(6,'(a)') ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!DEC$ ENDIF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_stella
