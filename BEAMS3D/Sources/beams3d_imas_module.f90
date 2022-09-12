!-----------------------------------------------------------------------
!     Module:        BEAMS3D_IMAS_MODULE
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          09/22/2021
!     Description:   Provides IMAS interfaces for running BEAMS3D from 
!                    a variety of sources.
!-----------------------------------------------------------------------
MODULE BEAMS3D_IMAS_MODULE
#if defined(IMAS)

CONTAINS
!-----------------------------------------------------------------------
!     SUBROUTINE:        BEAMS3D_IMAS_INIT
!     PURPOSE:           Handles Initialization of BEAMS3D
!-----------------------------------------------------------------------
SUBROUTINE BEAMS3D_IMAS_INIT(INDATA_XML, status_code, status_message)
  !---------------------------------------------------------------------
  !     Libraries
  !---------------------------------------------------------------------
  USE ids_schemas
  USE beams3d_input_mod, ONLY: init_beams3d_input, read_beams3d_input

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
  !     BEGIN EXECUTION
  !---------------------------------------------------------------------

  status_code = 0
  !status_message = 'OK'

  !----  Set defaults ----
  CALL init_beams3d_input

  !----  Read XML ----
  CALL beams3d_input_imas(INDATA_XML,status_code,status_message)

  !----  Read XML ----
  CALL read_beams3d_input('IMAS',status_code)

  RETURN
END SUBROUTINE BEAMS3D_IMAS_INIT

!-----------------------------------------------------------------------
!     SUBROUTINE:        BEAMS3D_INPUT_IMAS
!     PURPOSE:           Handles setting up the input variables
!-----------------------------------------------------------------------
SUBROUTINE BEAMS3D_INPUT_IMAS(INDATA_XML, status_code, status_message)
  !---------------------------------------------------------------------
  !     Libraries
  !---------------------------------------------------------------------
  USE ids_schemas
  USE xml2eg_mdl, only: xml2eg_parse_memory, xml2eg_get, type_xml2eg_document, xml2eg_free_doc, set_verbose
  USE beams3d_input_mod

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
  REAL(rprec) :: Xmn(0:mpol1d)

  TYPE(type_xml2eg_document) :: doc

  !---------------------------------------------------------------------
  !     BEGIN EXECUTION
  !---------------------------------------------------------------------

  !----  Intitailizations (This is a hack to set defaults)
  CALL vsetup(0)
  iunit = -327
  CALL read_indata_namelist(iunit,status_code)
  
  !----  Now setup the run based on the xml
  CALL xml2eg_parse_memory(INPUT_XML%parameters_value,doc)
  CALL SET_VERBOSE(.TRUE.) ! Only needed if you want to see what's going on in the parsing
  CALL xml2eg_get(doc,'NR',nr)
  CALL xml2eg_get(doc,'NPHI',nphi)
  CALL xml2eg_get(doc,'NZ',nz)
  CALL xml2eg_get(doc,'RMIN',rmin)
  CALL xml2eg_get(doc,'RMAX',rmax)
  CALL xml2eg_get(doc,'ZMIN',zmin)
  CALL xml2eg_get(doc,'ZMAX',zmax)
  CALL xml2eg_get(doc,'PHIMIN',phimin)
  CALL xml2eg_get(doc,'PHIMAX',phimax)
  CALL xml2eg_get(doc,'NPARTICLES_START',nparticles_start)
  CALL xml2eg_get(doc,'R_START_IN',r_start_in)
  CALL xml2eg_get(doc,'PHI_START_IN',phi_start_in)
  CALL xml2eg_get(doc,'Z_START_IN',z_start_in)
  CALL xml2eg_get(doc,'VLL_START_IN',vll_start_in)
  CALL xml2eg_get(doc,'MU_START_IN',mu_start_in)
  CALL xml2eg_get(doc,'T_END_IN',t_end_in)
  CALL xml2eg_get(doc,'CHARGE_IN',charge_in)
  CALL xml2eg_get(doc,'MASS_IN',mass_in)
  CALL xml2eg_get(doc,'ZATOM_IN',zatom_in)
  CALL xml2eg_get(doc,'NPOINC',npoinc)
  CALL xml2eg_get(doc,'FOLLOW_TOL',follow_tol)
  CALL xml2eg_get(doc,'INT_TYPE',int_type)
  CALL xml2eg_get(doc,'VC_ADAPT_TOL',vc_adapt_tol)
  CALL xml2eg_get(doc,'ADIST_BEAMS',adist_beams)
  CALL xml2eg_get(doc,'ASIZE_BEAMS',asize_beams)
  CALL xml2eg_get(doc,'DIV_BEAMS',div_beams)
  CALL xml2eg_get(doc,'E_BEAMS',E_beams)
  CALL xml2eg_get(doc,'P_BEAMS',P_beams)
  CALL xml2eg_get(doc,'DEX_BEAMS',dex_beams)
  CALL xml2eg_get(doc,'MASS_BEAMS',mass_beams)
  CALL xml2eg_get(doc,'CHARGE_BEAMS',charge_beams)
  CALL xml2eg_get(doc,'ZATOM_BEAMS',Zatom_beams)
  CALL xml2eg_get(doc,'R0_BEAMS',r_beams_in(:,1))
  CALL xml2eg_get(doc,'R1_BEAMS',r_beams_in(:,2))
  CALL xml2eg_get(doc,'PHI0_BEAMS',phi_beams_in(:,1))
  CALL xml2eg_get(doc,'PHI1_BEAMS',phi_beams_in(:,2))
  CALL xml2eg_get(doc,'Z0_BEAMS',z_beams_in(:,1))
  CALL xml2eg_get(doc,'Z1_BEAMS',z_beams_in(:,2))
  CALL xml2eg_get(doc,'TE_AUX_S',TE_AUX_S)
  CALL xml2eg_get(doc,'TE_AUX_F',TE_AUX_F)
  CALL xml2eg_get(doc,'NE_AUX_S',NE_AUX_S)
  CALL xml2eg_get(doc,'NE_AUX_F',NE_AUX_F)
  CALL xml2eg_get(doc,'TI_AUX_S',TI_AUX_S)
  CALL xml2eg_get(doc,'TI_AUX_F',TI_AUX_F)
  CALL xml2eg_get(doc,'POT_AUX_S',POT_AUX_S)
  CALL xml2eg_get(doc,'POT_AUX_F',POT_AUX_F)
  CALL xml2eg_get(doc,'ZEFF_AUX_S',ZEFF_AUX_S)
  CALL xml2eg_get(doc,'ZEFF_AUX_F',ZEFF_AUX_F)
  CALL xml2eg_get(doc,'NI_AUX_Z',NI_AUX_Z)
  CALL xml2eg_get(doc,'NI_AUX_M',NI_AUX_M)
  CALL xml2eg_get(doc,'NI_AUX_S',NI_AUX_S)
  CALL xml2eg_get(doc,'NI01_AUX_F',NI_AUX_F(1,:))
  CALL xml2eg_get(doc,'NI02_AUX_F',NI_AUX_F(2,:))
  CALL xml2eg_get(doc,'NI03_AUX_F',NI_AUX_F(3,:))
  CALL xml2eg_get(doc,'NI04_AUX_F',NI_AUX_F(4,:))
  CALL xml2eg_get(doc,'NE_SCALE',NE_SCALE)
  CALL xml2eg_get(doc,'TE_SCALE',TE_SCALE)
  CALL xml2eg_get(doc,'TI_SCALE',TI_SCALE)
  CALL xml2eg_get(doc,'ZEFF_SCALE',ZEFF_SCALE)
  CALL xml2eg_get(doc,'PLASMA_MASS',plasma_mass)
  CALL xml2eg_get(doc,'PLASMA_ZMEAN',plasma_zmean)
  CALL xml2eg_get(doc,'THERM_FACTOR',therm_factor)
  CALL xml2eg_get(doc,'FUSION_SCALE',fusion_scale)

  ! Make sure to clean up after you!!
  ! When calling "xml2eg_parse_memory" memory was allocated in the "doc" object.
  ! This memory is freed by "xml2eg_free_doc(doc)"
  CALL xml2eg_free_doc(doc)

  RETURN
END SUBROUTINE BEAMS3D_INPUT_IMAS

#endif
END MODULE BEAMS3D_IMAS_MODULE