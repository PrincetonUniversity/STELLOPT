!-----------------------------------------------------------------------
!     Module:        vmec_imas_module
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          09/22/2021
!     Description:   Provides IMAS interfaces for running VMEC from 
!                    IMAS.  Currently we assume the code runs
!                    standalone.
!-----------------------------------------------------------------------
MODULE VMEC_IMAS_MODULE
#if defined(IMAS)

   
!-----------------------------------------------------------------------
!     SUBROUTINE:        VMEC_IMAS
!     PURPOSE:           RUNS VMEC ASSUMING XML DEFINES RUN LIKE THE
!                        NAMELIST
!-----------------------------------------------------------------------
SUBROUTINE VMEC_IMAS(IDS_EQ_OUT, INDATA, status_code, status_message)
  !---------------------------------------------------------------------
  !     Libraries
  !---------------------------------------------------------------------
  USE ids_schemas
  USE parallel_vmec_module
  USE parallel_include_module

  !---------------------------------------------------------------------
  !     INPUT/OUTPUT VARIABLES
  !        IDS_EQ_OUT : Equilibrium output
  !        INDATA_XML : VMEC INDATA namelist AS XML
  !        STATUS_CODE : STATUS Flag
  !        STATUS_MESSAGE : Text Message
  !---------------------------------------------------------------------
  IMPLICIT NONE
  TYPE(ids_equilibrium), INTENT(OUT) :: IDS_EQ_OUT
  TYPE(ids_parameters_input), INTENT(IN) :: INDATA
  INTEGER, INTENT(OUT) :: status_code
  CHARACTER(LEN=:), POINTER, INTENT(OUT) :: status_message

  !---------------------------------------------------------------------
  !     SUBROUTINE VARIABLES
  !---------------------------------------------------------------------
  LOGICAL :: lmpi_flag, lscreen
  INTEGER :: impi_flag, ivmec_flag
  INTEGER :: ictrl(5)
  CHARACTER(len = 128)    :: reset_string, filename

  !---------------------------------------------------------------------
  !     BEGIN EXECUTION
  !---------------------------------------------------------------------

  !----  MPI initialisation ----
  CALL MPI_initialized(lmpi_flag, impi_flag)
  if (.not. lmpi_flag)   call MPI_INIT(impi_flag)

  !----  MIMIC CALL MyEnvVariables
  PARVMEC = .TRUE.

  !----  MIMIC CALL InitializeParallel
  CALL MPI_Comm_rank(MPI_COMM_WORLD,grank,MPI_ERR)
  CALL MPI_Comm_size(MPI_COMM_WORLD,gnranks,MPI_ERR)
  
  !----  Mimic Read_Indata
  CALL VMEC_INDATA_IMAS(INDATA, status_code, status_message)

  !----  Run VMEC
  ictrl(1) = restart_flag + readimas_flag + timestep_flag + output_flag &
                          + cleanup_flag
  ictrl(2) = 0 ! ier_flag
  ictrl(3) = 0 ! numsteps
  ictrl(4) = 0 ! ns_index
  ictrl(5) = 0 ! Sequence ID
  reset_string = ''
  filename=''
  RVCCALLNUM = 1
  CALL runvmec(ictrl, extension(index_seq), lscreen, RVC_COMM, &
               reset_string)
  ivmec_flag = ictrl(2)

  !----  No need to call finalize

  !---------------------------------------------------------------------
  !     END EXECUTION - Don't touch below here
  !---------------------------------------------------------------------

  !----  MPI Finalisation ----
  call MPI_finalized(lmpi_flag, impi_flag)
  IF (.NOT. lmpi_flag)   CALL MPI_Finalize(impi_flag)
  
END SUBROUTINE VMEC_IMAS

!-----------------------------------------------------------------------
!     SUBROUTINE:        VMEC_INDATA_IMAS
!     PURPOSE:           Handles setting up the indata variables
!-----------------------------------------------------------------------
SUBROUTINE VMEC_INDATA_IMAS(INDATA, status_code, status_message)
  !---------------------------------------------------------------------
  !     Libraries
  !---------------------------------------------------------------------
  USE ids_schemas
  USE vmec_input

  !---------------------------------------------------------------------
  !     INPUT/OUTPUT VARIABLES
  !        INDATA_XML : VMEC INDATA namelist AS XML
  !        STATUS_CODE : STATUS Flag
  !        STATUS_MESSAGE : Text Message
  !---------------------------------------------------------------------
  IMPLICIT NONE
  TYPE(ids_parameters_input), INTENT(IN) :: INDATA
  INTEGER, INTENT(OUT) :: status_code
  CHARACTER(LEN=:), POINTER, INTENT(OUT) :: status_message

  !---------------------------------------------------------------------
  !     SUBROUTINE VARIABLES
  !---------------------------------------------------------------------
  LOGICAL :: lmpi_flag, lscreen
  INTEGER :: impi_flag, ivmec_flag
  INTEGER :: ictrl(5)
  CHARACTER(len = 128)    :: reset_string

  !---------------------------------------------------------------------
  !     BEGIN EXECUTION
  !---------------------------------------------------------------------

  !----  Intitailizations (This is a hack)
  iunit = -1
  CALL read_indata_namelist(iunit,status_code)
  
  !----  Now setup the run based on the xml

END SUBROUTINE VMEC_INDATA_IMAS

#endif
END MODULE VMEC_IMAS_MODULE