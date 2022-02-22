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

CONTAINS   
!-----------------------------------------------------------------------
!     SUBROUTINE:        VMEC_IMAS
!     PURPOSE:           RUNS VMEC ASSUMING XML DEFINES RUN LIKE THE
!                        NAMELIST
!-----------------------------------------------------------------------
SUBROUTINE VMEC_IMAS(IDS_EQ_OUT, INDATA_XML, status_code, status_message)
  !---------------------------------------------------------------------
  !     Libraries
  !---------------------------------------------------------------------
  USE ids_schemas
  USE parallel_vmec_module
  USE parallel_include_module
  USE vmec_params

  !---------------------------------------------------------------------
  !     INPUT/OUTPUT VARIABLES
  !        IDS_EQ_OUT : Equilibrium output
  !        INDATA_XML : VMEC INDATA namelist AS XML
  !        STATUS_CODE : STATUS Flag
  !        STATUS_MESSAGE : Text Message
  !---------------------------------------------------------------------
  IMPLICIT NONE
  TYPE(ids_equilibrium), INTENT(OUT) :: IDS_EQ_OUT
  TYPE(ids_parameters_input), INTENT(IN) :: INDATA_XML
  INTEGER, INTENT(OUT) :: status_code
  CHARACTER(LEN=:), POINTER, INTENT(OUT) :: status_message

  !---------------------------------------------------------------------
  !     SUBROUTINE VARIABLES
  !---------------------------------------------------------------------
  LOGICAL :: lmpi_flag, lscreen
  INTEGER :: impi_flag, ivmec_flag, RVC_COMM
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
  CALL MPI_Comm_rank(MPI_COMM_WORLD,grank,impi_flag)
  CALL MPI_Comm_size(MPI_COMM_WORLD,gnranks,impi_flag)

  !----  Duplicate the communicator
  CALL MPI_COMM_DUP(MPI_COMM_WORLD,RVC_COMM,impi_flag)
  
  !----  Mimic Read_Indata
  CALL VMEC_INDATA_IMAS(INDATA_XML, status_code, status_message)

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
  CALL runvmec(ictrl, filename, lscreen, RVC_COMM, &
               reset_string)
  ivmec_flag = ictrl(2)

  !----  Write out EQUILIBRIUM


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
SUBROUTINE VMEC_INDATA_IMAS(INDATA_XML, status_code, status_message)
  !---------------------------------------------------------------------
  !     Libraries
  !---------------------------------------------------------------------
  USE ids_schemas
  USE xml2eg_mdl, only: xml2eg_parse_memory, xml2eg_get, type_xml2eg_document, xml2eg_free_doc, set_verbose
  USE vmec_input

  !---------------------------------------------------------------------
  !     INPUT/OUTPUT VARIABLES
  !        INDATA_XML : VMEC INDATA namelist AS XML
  !        STATUS_CODE : STATUS Flag
  !        STATUS_MESSAGE : Text Message
  !---------------------------------------------------------------------
  IMPLICIT NONE
  TYPE(ids_parameters_input), INTENT(IN) :: INDATA_XML
  INTEGER, INTENT(OUT) :: status_code
  CHARACTER(LEN=:), POINTER, INTENT(OUT) :: status_message

  !---------------------------------------------------------------------
  !     SUBROUTINE VARIABLES
  !---------------------------------------------------------------------
  LOGICAL :: lmpi_flag, lscreen
  INTEGER :: impi_flag, ivmec_flag, iunit
  INTEGER :: ictrl(5)
  CHARACTER(len = 128)    :: reset_string

  TYPE(type_xml2eg_document) :: doc

  !---------------------------------------------------------------------
  !     BEGIN EXECUTION
  !---------------------------------------------------------------------

  !----  Intitailizations (This is a hack to set defaults)
  iunit = -1
  CALL read_indata_namelist(iunit,status_code)
  
  !----  Now setup the run based on the xml
  CALL xml2eg_parse_memory(INDATA_XML%parameters_value,doc)
  CALL SET_VERBOSE(.TRUE.) ! Only needed if you want to see what's going on in the parsing
  CALL xml2eg_get(doc,'MGRID_FILE',mgrid_file)
  CALL xml2eg_get(doc,'NFP',nfp)
  CALL xml2eg_get(doc,'NCURR',ncurr)
  CALL xml2eg_get(doc,'NSTEP',nstep)
  CALL xml2eg_get(doc,'NVACSKIP',nvacskip)
  CALL xml2eg_get(doc,'DELT',delt)
  CALL xml2eg_get(doc,'GAMMA',gamma)
  CALL xml2eg_get(doc,'AM',am)
  CALL xml2eg_get(doc,'AI',ai)
  CALL xml2eg_get(doc,'AC',ac)
  CALL xml2eg_get(doc,'PCURR_TYPE',pcurr_type)
  CALL xml2eg_get(doc,'PMASS_TYPE',pmass_type)
  CALL xml2eg_get(doc,'PIOTA_TYPE',piota_type)
  CALL xml2eg_get(doc,'AM_AUX_S',am_aux_s)
  CALL xml2eg_get(doc,'AM_AUX_F',am_aux_f)
  CALL xml2eg_get(doc,'AI_AUX_S',ai_aux_s)
  CALL xml2eg_get(doc,'AI_AUX_F',ai_aux_f)
  CALL xml2eg_get(doc,'AC_AUX_S',ac_aux_s)
  CALL xml2eg_get(doc,'AC_AUX_F',ac_aux_f)
  !CALL xml2eg_get(doc,'RBC',rbc)
  !CALL xml2eg_get(doc,'ZBS',zbs)
  !CALL xml2eg_get(doc,'RBS',rbs)
  !CALL xml2eg_get(doc,'ZBC',zbc)
  CALL xml2eg_get(doc,'SPRES_PED',spres_ped)
  CALL xml2eg_get(doc,'PRES_SCALE',pres_scale)
  CALL xml2eg_get(doc,'RAXIS_CC',raxis_cc)
  CALL xml2eg_get(doc,'ZAXIS_CS',zaxis_cs)
  CALL xml2eg_get(doc,'RAXIS_CS',raxis_cs)
  CALL xml2eg_get(doc,'ZAXIS_CC',zaxis_cc)
  CALL xml2eg_get(doc,'MPOL',mpol)
  CALL xml2eg_get(doc,'NTOR',ntor)
  CALL xml2eg_get(doc,'NTHETA',ntheta)
  CALL xml2eg_get(doc,'NZETA',nzeta)
  CALL xml2eg_get(doc,'NITER_ARRAY',niter_array)
  CALL xml2eg_get(doc,'NS_ARRAY',ns_array)
  CALL xml2eg_get(doc,'FTOL_ARRAY',ftol_array)
  CALL xml2eg_get(doc,'TCON0',tcon0)
  CALL xml2eg_get(doc,'CURTOR',curtor)
  CALL xml2eg_get(doc,'EXTCUR',extcur)
  CALL xml2eg_get(doc,'PHIEDGE',phiedge)
  CALL xml2eg_get(doc,'BLOAT',bloat)
  CALL xml2eg_get(doc,'LFREEB',lfreeb)
  CALL xml2eg_get(doc,'LASYM',lasym)

  ! Make sure to clean up after you!!
  ! When calling "xml2eg_parse_memory" memory was allocated in the "doc" object.
  ! This memory is freed by "xml2eg_free_doc(doc)"
  CALL xml2eg_free_doc(doc)

END SUBROUTINE VMEC_INDATA_IMAS

#endif
END MODULE VMEC_IMAS_MODULE
