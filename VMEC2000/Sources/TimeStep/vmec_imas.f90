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
!     SUBROUTINE:        VMEC_IMAS_INIT
!     PURPOSE:           Handles Initialization of VMEC
!-----------------------------------------------------------------------
SUBROUTINE VMEC_IMAS_INIT(INDATA_XML, status_code, status_message)
  !---------------------------------------------------------------------
  !     Libraries
  !---------------------------------------------------------------------
  USE ids_schemas

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

END SUBROUTINE VMEC_IMAS_INIT

!-----------------------------------------------------------------------
!     SUBROUTINE:        VMEC_IMAS
!     PURPOSE:           RUNS VMEC ASSUMING XML DEFINES RUN LIKE THE
!                        NAMELIST
!-----------------------------------------------------------------------
SUBROUTINE VMEC_IMAS(IDS_EQ_IN, IDS_EQ_OUT, INDATA_XML, status_code, status_message)
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
  TYPE(ids_equilibrium), INTENT(IN) :: IDS_EQ_IN
  TYPE(ids_equilibrium), INTENT(OUT) :: IDS_EQ_OUT
  TYPE(ids_parameters_input), INTENT(IN) :: INDATA_XML
  INTEGER, INTENT(OUT) :: status_code
  CHARACTER(LEN=:), POINTER, INTENT(OUT) :: status_message

  !---------------------------------------------------------------------
  !     SUBROUTINE VARIABLES
  !---------------------------------------------------------------------
  LOGICAL :: lmpi_flag, lscreen
  INTEGER :: impi_flag, ivmec_flag, RVC_COMM
  INTEGER :: ictrl(5), myseq
  CHARACTER(len = 128)    :: reset_string, filename

  !---------------------------------------------------------------------
  !     BEGIN EXECUTION
  !---------------------------------------------------------------------

  !---- Default outputs
  status_code = 0
  !status_message = 'OK'

  !----  MPI initialisation ----
  CALL MPI_initialized(lmpi_flag, impi_flag)
  if (.not. lmpi_flag)   call MPI_INIT(impi_flag)

  !----  MIMIC CALL MyEnvVariables
  CALL MyEnvVariables

  !----  MIMIC CALL InitializeParallel
  CALL MPI_Comm_rank(MPI_COMM_WORLD,grank,impi_flag)
  CALL MPI_Comm_size(MPI_COMM_WORLD,gnranks,impi_flag)

  !----  Duplicate the communicator
  CALL MPI_COMM_DUP(MPI_COMM_WORLD,RVC_COMM,impi_flag)
  
  !----  Mimic Read_Indata
  CALL VMEC_INDATA_IMAS(INDATA_XML, status_code, status_message)

  !----  Interface with Equilibirum IDS
  CALL VMEC_EQIN_IMAS(IDS_EQ_IN,status_code,status_message)

  !---- Setup VMEC Run
  myseq=0
  ictrl(1) = restart_flag + imasrun_flag + timestep_flag +  &
                          + cleanup_flag
  ictrl(2) = 0 ! ier_flag
  ictrl(3) = -1 ! numsteps
  ictrl(4) = -1 ! ns_index
  ictrl(5) = 0 ! Sequence ID
  reset_string = ''
  filename=''
  RVCCALLNUM = 1
  lscreen = .True.
  NS_RESLTN = 0 ! Need to do this otherwise situations arrise which cause problems.

  !----  Run VMEC
  CALL runvmec(ictrl, filename, lscreen, RVC_COMM, &
               reset_string)
  status_code = ictrl(2)
  STOP

  !----  Write out EQUILIBRIUM


  !----  No need to call finalize
  CALL FinalizeSurfaceComm(NS_COMM)
  CALL FinalizeRunVmec(RUNVMEC_COMM_WORLD)

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
  CALL xml2eg_parse_memory(INDATA_XML%parameters_value,doc)
  CALL SET_VERBOSE(.TRUE.) ! Only needed if you want to see what's going on in the parsing
  CALL xml2eg_get(doc,'MGRID_FILE',mgrid_file)
  CALL xml2eg_get(doc,'NFP',nfp)
  CALL xml2eg_get(doc,'NCURR',ncurr)
  CALL xml2eg_get(doc,'NSTEP',nstep)
  !CALL xml2eg_get(doc,'NVACSKIP',nvacskip)
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
  CALL xml2eg_get(doc,'SPRES_PED',spres_ped)
  CALL xml2eg_get(doc,'PRES_SCALE',pres_scale)
  CALL xml2eg_get(doc,'RAXIS_CC',raxis_cc)
  CALL xml2eg_get(doc,'ZAXIS_CS',zaxis_cs)
  !CALL xml2eg_get(doc,'RAXIS_CS',raxis_cs)
  !CALL xml2eg_get(doc,'ZAXIS_CC',zaxis_cc)
  CALL xml2eg_get(doc,'MPOL',mpol)
  CALL xml2eg_get(doc,'NTOR',ntor)
  !CALL xml2eg_get(doc,'NTHETA',ntheta)
  CALL xml2eg_get(doc,'NZETA',nzeta)
  CALL xml2eg_get(doc,'NITER_ARRAY',niter_array)
  CALL xml2eg_get(doc,'NS_ARRAY',ns_array)
  CALL xml2eg_get(doc,'FTOL_ARRAY',ftol_array)
  CALL xml2eg_get(doc,'TCON0',tcon0)
  !CALL xml2eg_get(doc,'PRECON_TYPE',precon_type)
  !CALL xml2eg_get(doc,'PREC2D_THRESHOLD',prec2d_threshold)
  CALL xml2eg_get(doc,'CURTOR',curtor)
  CALL xml2eg_get(doc,'EXTCUR',extcur)
  CALL xml2eg_get(doc,'PHIEDGE',phiedge)
  CALL xml2eg_get(doc,'BLOAT',bloat)
  !CALL xml2eg_get(doc,'LFORBAL',lforbal)
  CALL xml2eg_get(doc,'LFREEB',lfreeb)
  CALL xml2eg_get(doc,'LASYM',lasym)
  !CALL xml2eg_get(doc,'LRFP',lrfp)
  !CALL xml2eg_get(doc,'LBSUBS',lbsubs)
  !CALL xml2eg_get(doc,'LNYQUIST',lnyquist)
  CALL xml2eg_get(doc,'RBC_np0',Xmn)
  rbc(0,:) = Xmn
  CALL xml2eg_get(doc,'ZBS_np0',Xmn)
  zbs(0,:) = Xmn

  ! Make sure to clean up after you!!
  ! When calling "xml2eg_parse_memory" memory was allocated in the "doc" object.
  ! This memory is freed by "xml2eg_free_doc(doc)"
  CALL xml2eg_free_doc(doc)

  RETURN
END SUBROUTINE VMEC_INDATA_IMAS

!-----------------------------------------------------------------------
!     SUBROUTINE:        VMEC_INDATA_IMAS
!     PURPOSE:           Handles setting up the indata variables
!-----------------------------------------------------------------------
SUBROUTINE VMEC_EQIN_IMAS(IDS_EQ_IN, status_code, status_message)
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
  TYPE(ids_equilibrium), INTENT(IN) :: IDS_EQ_IN
  INTEGER, INTENT(OUT) :: status_code
  CHARACTER(LEN=:), POINTER, INTENT(OUT) :: status_message

  !---------------------------------------------------------------------
  !     SUBROUTINE VARIABLES
  !---------------------------------------------------------------------
  INTEGER :: cocos_index, npts_imas, itime, u, mn
  REAL*8  :: B0, pi2
  REAL*8, ALLOCATABLE, DIMENSION(:) :: R_BND, Z_BND, radius, theta

  !---------------------------------------------------------------------
  !     BEGIN EXECUTION
  !---------------------------------------------------------------------
  status_code = 0
  itime = 1
  pi2 = 8.0*ATAN(1.0)

  !----  Check the IDS for timeslices
  IF (.not. ASSOCIATED(IDS_EQ_IN%time_slice)) THEN
     WRITE(*,*) 'No time slices in this equilibrium'
  END IF

  !----  Check Cocos Index
  cocos_index = IDS_EQ_IN%time_slice(itime)%profiles_2d(1)%grid_type%index
  WRITE(*,*) ' IDS COCOS Convention: ',cocos_index

  !----  Print Grid Type
  IF (ASSOCIATED(IDS_EQ_IN%time_slice(itime)%profiles_2d(1)%grid_type%name)) &
     WRITE(*,*) 'IDS Grid Type: ',IDS_EQ_IN%time_slice(itime)%profiles_2d(1)%grid_type%name

  !----  Pass Shot Information
  time_slice = IDS_EQ_IN%time(itime)


  !----  Get Global quantities
  curtor = IDS_EQ_IN%time_slice(itime)%global_quantities%ip ! Total Toroidal Current
  !PRINT *,'======'
  !PRINT *,curtor

  !----  Get 1D Profiles
  !      Toroidal Flux
  AI_AUX_S(:) = -1
  AM_AUX_S(:) = -1
  AC_AUX_S(:) = -1
  npts_imas = size(IDS_EQ_IN%time_slice(itime)%profiles_1d%phi)
  AM_AUX_S(1:npts_imas)  = IDS_EQ_IN%time_slice(itime)%profiles_1d%phi
  AM_AUX_S(1:npts_imas)  = AM_AUX_S(1:npts_imas)-AM_AUX_S(1)
  AM_AUX_S(1:npts_imas)  = AM_AUX_S(1:npts_imas)/AM_AUX_S(npts_imas)
  PHIEDGE = IDS_EQ_IN%time_slice(itime)%profiles_1d%phi(npts_imas)
  !PRINT *,'======'
  !PRINT *,PHIEDGE
  !     Pressure
  npts_imas = size(IDS_EQ_IN%time_slice(itime)%profiles_1d%pressure)
  AM_AUX_F(1:npts_imas)  = IDS_EQ_IN%time_slice(itime)%profiles_1d%pressure
  PRES_SCALE = 1.0
  PMASS_TYPE = 'Akima_spline'
  !PRINT *,'======'
  !PRINT *,AM_AUX_F(1:npts_imas)
  !     Rotational Transform
  AI_AUX_S  = AM_AUX_S
  npts_imas = size(IDS_EQ_IN%time_slice(itime)%profiles_1d%q)
  AI_AUX_F(1:npts_imas)  = 1.0/IDS_EQ_IN%time_slice(itime)%profiles_1d%q
  PIOTA_TYPE = 'Akima_spline'
  !PRINT *,'======'
  !PRINT *,AI_AUX_F(1:npts_imas)
  !     Current Density
  AC_AUX_S  = AM_AUX_S
  npts_imas = size(IDS_EQ_IN%time_slice(itime)%profiles_1d%j_parallel)
  AC_AUX_F(1:npts_imas)  = IDS_EQ_IN%time_slice(itime)%profiles_1d%j_parallel
  PCURR_TYPE = 'Akima_spline_ip'
  !PRINT *,'======'
  !PRINT *,AC_AUX_F(1:npts_imas)
  NCURR = 0

  !----  Get Magnetic Axis position
  RAXIS_CC = 0; RAXIS_CS = 0; ZAXIS_CC = 0; ZAXIS_CS = 0
  RAXIS_CC(0) = IDS_EQ_IN%time_slice(itime)%global_quantities%magnetic_axis%r
  ZAXIS_CC(0) = IDS_EQ_IN%time_slice(itime)%global_quantities%magnetic_axis%z

  !----  Get Boundary Shape
  MPOL = 24
  npts_imas = SIZE(IDS_EQ_IN%time_slice(itime)%boundary%outline%R)
  ALLOCATE(R_BND(npts_imas), Z_BND(npts_imas), radius(npts_imas), theta(npts_imas))
  R_BND = IDS_EQ_IN%time_slice(itime)%boundary%outline%R
  Z_BND = IDS_EQ_IN%time_slice(itime)%boundary%outline%Z


  
  !theta = ATAN2(Z_BND-ZAXIS_CS(0),R_BND-RAXIS_CC(0))
  DO u = 1, npts_imas-2
     theta(u) = pi2*DBLE(u-1)/DBLE(npts_imas-2)
  END DO
  RBC(:,:) = 0.0; RBS(:,:) = 0.0; ZBC(:,:) = 0.0; ZBS(:,:) = 0.0;
  DO u = 1, npts_imas-2
     DO mn = 0, MPOL-1
       RBC(0,mn) = RBC(0,mn) + R_BND(u)*COS(mn*theta(u)) 
       ZBC(0,mn) = ZBC(0,mn) + Z_BND(u)*COS(mn*theta(u))
       RBS(0,mn) = RBS(0,mn) + R_BND(u)*SIN(mn*theta(u))
       ZBS(0,mn) = ZBS(0,mn) + Z_BND(u)*SIN(mn*theta(u))
     END DO
  END DO
  RBC(0,0) = RBC(0,0)/(npts_imas-2)
  ZBC(0,0) = ZBC(0,0)/(npts_imas-2)
  RBS(0,0) = RBS(0,0)/(npts_imas-2)
  ZBS(0,0) = ZBS(0,0)/(npts_imas-2)
  RBC(0,1:MPOL-1) = 2*RBC(0,1:MPOL-1)/(npts_imas-2)
  ZBC(0,1:MPOL-1) = 2*ZBC(0,1:MPOL-1)/(npts_imas-2)
  RBS(0,1:MPOL-1) = 2*RBS(0,1:MPOL-1)/(npts_imas-2)
  ZBS(0,1:MPOL-1) = 2*ZBS(0,1:MPOL-1)/(npts_imas-2)
  LASYM = .TRUE.

  DEALLOCATE(R_BND, Z_BND, radius, theta)

  RETURN
  

END SUBROUTINE VMEC_EQIN_IMAS

#endif
END MODULE VMEC_IMAS_MODULE
