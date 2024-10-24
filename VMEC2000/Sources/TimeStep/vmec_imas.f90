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
  INTEGER :: ictrl(5), myseq, npts_imas
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
  IF (status_code .ne. 0) RETURN

  !---- Setup VMEC Run
  myseq=0
  ictrl(1) = restart_flag + imasrun_flag + timestep_flag + output_flag
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

  !----  Write out EQUILIBRIUM
  IF (grank == 0) CALL VMEC_EQOUT_IMAS(IDS_EQ_OUT, status_code, status_message)

  !----  Copy stuff from the input to the output
  npts_imas = SIZE(IDS_EQ_IN%time)
  ALlOCATE(IDS_EQ_OUT%time(npts_imas))
  IDS_EQ_OUT%time = IDS_EQ_IN%time


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
!     SUBROUTINE:        VMEC_EQIN_IMAS
!     PURPOSE:           Handles the equilibrium IDS.
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
  REAL*8  :: B0, pi2, s_temp, x
  REAL*8, ALLOCATABLE, DIMENSION(:) :: R_BND, Z_BND, radius, theta, &
                                       s_in, f_in, f2_in, f3_in

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


  !----  Get Toroidal Current
  IF (ids_is_valid(IDS_EQ_IN%time_slice(itime)%global_quantities%ip)) THEN
     curtor = IDS_EQ_IN%time_slice(itime)%global_quantities%ip ! Total Toroidal Current
  ELSE
     status_code = -1
     ALLOCATE(CHARACTER(len=256) :: status_message)
     status_message = 'Toroidal Current (IDS_EQ_IN%time_slice(itime)%global_quantities%ip) not present in Input IDS'
     RETURN
  END IF

  !----  Get Toroidal Flux
  IF (ids_is_valid(IDS_EQ_IN%time_slice(itime)%profiles_1d%phi)) THEN
     npts_imas = size(IDS_EQ_IN%time_slice(itime)%profiles_1d%phi)
     phiedge = IDS_EQ_IN%time_slice(itime)%profiles_1d%phi(npts_imas) ! Total Toroidal Current
  ELSE
     status_code = -1
     ALLOCATE(CHARACTER(len=256) :: status_message)
     status_message = 'Toroidal Flux (IDS_EQ_IN%time_slice(itime)%profiles_1d%phi) not present in Input IDS'
     RETURN
  END IF

  !----  Get 1D Profiles
  !      Note that IMAS stores everything on a grid which may
  !      extend beyond s=[0,1]
  IF (ids_is_valid(IDS_EQ_IN%time_slice(itime)%profiles_1d%rho_tor_norm)) THEN
     npts_imas = size(IDS_EQ_IN%time_slice(itime)%profiles_1d%rho_tor_norm)
     ALLOCATE(s_in(npts_imas), f_in(npts_imas), f2_in(npts_imas), f3_in(npts_imas))
     s_in = IDS_EQ_IN%time_slice(itime)%profiles_1d%rho_tor_norm
  ELSE
     status_code = -1
     ALLOCATE(CHARACTER(len=256) :: status_message)
     status_message = 'Normalzied Toroidal Flux (IDS_EQ_IN%time_slice(itime)%profiles_1d%rho_tor_norm) not present in Input IDS'
     RETURN
  END IF
  ! Pressure
  IF (ids_is_valid(IDS_EQ_IN%time_slice(itime)%profiles_1d%pressure)) THEN
     f_in = IDS_EQ_IN%time_slice(itime)%profiles_1d%pressure
  ELSE
     status_code = -1
     ALLOCATE(CHARACTER(len=256) :: status_message)
     status_message = 'Pressure (IDS_EQ_IN%time_slice(itime)%profiles_1d%pressure) not present in Input IDS'
     RETURN
  END IF
  ! q
  IF (ids_is_valid(IDS_EQ_IN%time_slice(itime)%profiles_1d%q)) THEN
     f2_in = 1.0/IDS_EQ_IN%time_slice(itime)%profiles_1d%q
  ELSE
     status_code = -1
     ALLOCATE(CHARACTER(len=256) :: status_message)
     status_message = 'Safety Factor (IDS_EQ_IN%time_slice(itime)%profiles_1d%q) not present in Input IDS'
     RETURN
  END IF
  ! J_tor
  IF (ids_is_valid(IDS_EQ_IN%time_slice(itime)%profiles_1d%j_tor)) THEN
     f3_in = IDS_EQ_IN%time_slice(itime)%profiles_1d%j_tor
  ELSE
     status_code = -1
     ALLOCATE(CHARACTER(len=256) :: status_message)
     status_message = 'Toroidal Current Density (IDS_EQ_IN%time_slice(itime)%profiles_1d%pressure) not present in Input IDS'
     RETURN
  END IF
  ! Now Interpolate
  AI_AUX_S(:) = -1
  AM_AUX_S(:) = -1
  AC_AUX_S(:) = -1
  DO u = 1, ndatafmax-1
     s_temp = DBLE(u-1)/DBLE(ndatafmax-2)
     AM_AUX_S(u) = s_temp
     AI_AUX_S(u) = s_temp
     AC_AUX_S(u) = s_temp
     mn = MIN(MAX(COUNT(s_in<s_temp),1),npts_imas-1)
     x  = (s_temp - s_in(mn))/(s_in(mn+1)-s_in(mn))
     AM_AUX_F(u) = f_in(mn)*(1.0-x)+f_in(mn+1)*x
     AI_AUX_F(u) = f2_in(mn)*(1.0-x)+f2_in(mn+1)*x
     AC_AUX_F(u) = f3_in(mn)*(1.0-x)+f3_in(mn+1)*x
  END DO
  DEALLOCATE(s_in, f_in, f2_in, f3_in)
  PRES_SCALE = 1.0
  NCURR = 1
  PMASS_TYPE = 'Akima_spline'
  PIOTA_TYPE = 'Akima_spline'
  PCURR_TYPE = 'Akima_spline_ip'

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

!-----------------------------------------------------------------------
!     SUBROUTINE:        VMEC_EQOUT_IMAS
!     PURPOSE:           Handles setting up the indata variables
!-----------------------------------------------------------------------
SUBROUTINE VMEC_EQOUT_IMAS(IDS_EQ_OUT, status_code, status_message)
  !---------------------------------------------------------------------
  !     Libraries
  !---------------------------------------------------------------------
  USE ids_schemas
  USE vmec_input
  USE vmec_params, ONLY: version_
  USE git_version_mod
  USE vmec_main, ONLY: ctor, wp, r00, z00, phi, chi, presf, jcurv, &
                       iotaf, jdotb, chipf, vp
  USE vmec_io, ONLY: betapol, betator, Aminor_p, b0, volume_p, &
                     cross_area_p, surf_area_p, circum_p
  USE vparams, ONLY: mu0
  USE vmec_dim, ONLY: ns, mnmax
  USE vmec_persistent, ONLY: xm, xn
  USE read_wout_mod, ONLY: rmnc,zmns

  !---------------------------------------------------------------------
  !     INPUT/OUTPUT VARIABLES
  !        INDATA_XML : VMEC INDATA namelist AS XML
  !        STATUS_CODE : STATUS Flag
  !        STATUS_MESSAGE : Text Message
  !---------------------------------------------------------------------
  IMPLICIT NONE
  TYPE(ids_equilibrium), INTENT(OUT) :: IDS_EQ_OUT
  INTEGER, INTENT(OUT) :: status_code
  CHARACTER(LEN=:), POINTER, INTENT(OUT) :: status_message

  !---------------------------------------------------------------------
  !     SUBROUTINE VARIABLES
  !---------------------------------------------------------------------
  INTEGER :: itime, i
  CHARACTER(8) ::date
  CHARACTER(10) ::times

  !---- Defaults
  status_code = 0
  itime = 1

  !---- Check Validity
  IF (.NOT. ASSOCIATED(IDS_EQ_OUT%time_slice)) THEN
     ALLOCATE(IDS_EQ_OUT%time_slice(itime))
  END IF
  IF (.NOT. ASSOCIATED(IDS_EQ_OUT%ids_properties%creation_date)) THEN
     ALLOCATE(IDS_EQ_OUT%ids_properties%creation_date(1))
  END IF

  !---- IDS PROPERTIES
  IDS_EQ_OUT%ids_properties%provider='USER'
  CALL DATE_AND_TIME(DATE=date,TIME=times)
  IDS_EQ_OUT%ids_properties%creation_date=date
  IF(IDS_EQ_OUT%ids_properties%homogeneous_time .lt. 0) &
     IDS_EQ_OUT%ids_properties%homogeneous_time=1


  !---- CODE FIELDS
  IDS_EQ_OUT%code%name='VMEC2000'
  IF (.NOT.ASSOCIATED(IDS_EQ_OUT%code%commit))&
     ALLOCATE(IDS_EQ_OUT%code%commit(1))
  IDS_EQ_OUT%code%commit=git_hash
  IF (.NOT.ASSOCIATED(IDS_EQ_OUT%code%version))&
     ALLOCATE(IDS_EQ_OUT%code%version(1))
  PRINT *,version_
  IDS_EQ_OUT%code%version=version_
  IF (.NOT.ASSOCIATED(IDS_EQ_OUT%code%repository))&
     ALLOCATE(IDS_EQ_OUT%code%repository(1))
  IDS_EQ_OUT%code%repository=git_repository

  !---- Vacuum Field Information
  !IDS_EQ_OUT%vacuum_toroidal_field%r0= Rvac
  IF (.NOT. ASSOCIATED(IDS_EQ_OUT%vacuum_toroidal_field%b0) ) &
     ALLOCATE(IDS_EQ_OUT%vacuum_toroidal_field%b0(size(IDS_EQ_OUT%time)))
  !IDS_EQ_OUT%vacuum_toroidal_field%b0(itime)=Bvac*Bvac_sign

  !---- Global quantities
  IDS_EQ_OUT%time_slice(itime)%global_quantities%beta_pol      = betapol
  IDS_EQ_OUT%time_slice(itime)%global_quantities%beta_tor      = betator
  IDS_EQ_OUT%time_slice(itime)%global_quantities%beta_normal   = 100*betator*Aminor_p*b0/(ctor/mu0)
  IDS_EQ_OUT%time_slice(itime)%global_quantities%Ip            = ctor/mu0
  IDS_EQ_OUT%time_slice(itime)%global_quantities%volume        = volume_p
  !IDS_EQ_OUT%time_slice(itime)%global_quantities%li_3          = xli_3
  IDS_EQ_OUT%time_slice(itime)%global_quantities%area          = cross_area_p
  IDS_EQ_OUT%time_slice(itime)%global_quantities%surface       = surf_area_p
  IDS_EQ_OUT%time_slice(itime)%global_quantities%length_pol    = circum_p
  IDS_EQ_OUT%time_slice(itime)%global_quantities%psi_axis      = chipf(1)
  IDS_EQ_OUT%time_slice(itime)%global_quantities%psi_boundary  = chipf(ns)
  IDS_EQ_OUT%time_slice(itime)%global_quantities%q_axis        = 1/iotaf(1)
  i = FLOOR(0.9025*ns)
  IDS_EQ_OUT%time_slice(itime)%global_quantities%q_95          = 1/iotaf(i)
  IDS_EQ_OUT%time_slice(itime)%global_quantities%energy_mhd    = wp
  !IDS_EQ_OUT%time_slice(itime)%global_quantities%psi_external_average    = ??
  !IDS_EQ_OUT%time_slice(itime)%global_quantities%plasma_inductance    = ??
  IDS_EQ_OUT%time_slice(itime)%global_quantities%magnetic_axis%r = r00
  IF (lasym) THEN
     IDS_EQ_OUT%time_slice(itime)%global_quantities%magnetic_axis%z = z00
  ELSE
     IDS_EQ_OUT%time_slice(itime)%global_quantities%magnetic_axis%z = 0.0
  END IF

  !---- PROFILES 1D
  ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%psi(ns))
  IDS_EQ_OUT%time_slice(itime)%profiles_1d%psi      = chi
  ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%phi(ns))
  IDS_EQ_OUT%time_slice(itime)%profiles_1d%phi      = phi
  ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%pressure(ns))
  IDS_EQ_OUT%time_slice(itime)%profiles_1d%pressure = presf/mu0
  ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%j_tor(ns))
  IDS_EQ_OUT%time_slice(itime)%profiles_1d%j_tor = jcurv/mu0
  ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%j_parallel(ns))
  IDS_EQ_OUT%time_slice(itime)%profiles_1d%j_parallel = jdotb
  ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%q(ns))
  IDS_EQ_OUT%time_slice(itime)%profiles_1d%q = 1.0/iotaf
  ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%rho_tor(ns))
  IDS_EQ_OUT%time_slice(itime)%profiles_1d%rho_tor = Aminor_p*(phi-phi(1))/(phi(ns)-phi(1))
  ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%rho_tor_norm(ns))
  IDS_EQ_OUT%time_slice(itime)%profiles_1d%rho_tor_norm = (phi-phi(1))/(phi(ns)-phi(1))
  ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%dpsi_drho_tor(ns))
  IDS_EQ_OUT%time_slice(itime)%profiles_1d%dpsi_drho_tor = chipf
  ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%volume(ns))
  IDS_EQ_OUT%time_slice(itime)%profiles_1d%volume(1) = 0
  DO i = 2, ns
     IDS_EQ_OUT%time_slice(itime)%profiles_1d%volume(i) = &
        IDS_EQ_OUT%time_slice(itime)%profiles_1d%volume(i-1)+vp(i)
  END DO
  ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%rho_volume_norm(ns))
  IDS_EQ_OUT%time_slice(itime)%profiles_1d%rho_volume_norm = &
     SQRT(IDS_EQ_OUT%time_slice(itime)%profiles_1d%volume/volume_p)
  ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%dvolume_drho_tor(ns))
  IDS_EQ_OUT%time_slice(itime)%profiles_1d%dvolume_drho_tor = vp/Aminor_p
  
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%f(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%f = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%dpressure_dpsi(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%dpressure_dpsi = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%f_df_dpsi(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%f_df_dpsi = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%magnetic_shear(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%magnetic_shear = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%r_inboard(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%r_inboard = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%r_outboard(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%r_outboard = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%geometric_axis%r(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%geometric_axis%r = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%geometric_axis%z(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%geometric_axis%z = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%elongation(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%elongation = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%triangularity_upper(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%triangularity_upper = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%triangularity_lower(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%triangularity_lower = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%squareness_upper_inner(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%squareness_upper_inner = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%squareness_upper_outer(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%squareness_upper_outer = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%squareness_lower_inner(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%squareness_lower_inner = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%squareness_lower_outer(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%squareness_lower_outer = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%dvolume_dpsi(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%dvolume_dpsi = vp/chip
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%area(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%area = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%darea_dpsi(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%darea_dpsi = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%darea_drho_tor(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%darea_drho_tor = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%surface(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%surface = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%trapped_fraction(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%trapped_fraction = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm1(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm1 = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm2(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm2 = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm3(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm3 = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm4(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm4 = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm5(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm5 = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm6(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm6 = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm7(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm7 = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm8(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm8 = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm9(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%gm9 = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%b_field_average(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%b_field_average = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%b_field_min(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%b_field_min = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%b_field_max(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%b_field_max = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%beta_pol(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%beta_pol = ?
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_1d%mass_density(ns))
  !IDS_EQ_OUT%time_slice(itime)%profiles_1d%mass_density = ?

  !---- PROFILES 2D
  ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_2d(1))
  IF (.NOT.ASSOCIATED(IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid_type%name))&
     ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid_type%name(1))
  IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid_type%name='inverse_rhotor_polar_fourier' ! inverse_rhotor_polar_fourier
  IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid_type%index=56 ! inverse_rhotor_polar_fourier
  IF (.NOT.ASSOCIATED(IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid_type%description))&
     ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid_type%description(1))
  IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid_type%description='Flux surface type with radial label sqrt[Phi/pi/B0] (dim1), Phi being toroidal flux, and Fourier modes in the polar poloidal angle (dim2)' ! inverse_rhotor_polar_fourier
  IF (.NOT.ASSOCIATED(IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid%dim1))&
     ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid%dim1(ns))
  DO i = 1, ns
     IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid%dim1(i) = DBLE(i-1)/DBLE(ns-1)
  END DO
  IF (.NOT.ASSOCIATED(IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid%dim2))&
     ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid%dim2(mnmax))
  DO i = 1, mnmax
     IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid%dim2(i) = xm(i)
  END DO
  IF (.NOT.ASSOCIATED(IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%r))&
     ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%r(ns,mnmax))
  IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%r = TRANSPOSE(rmnc)
  IF (.NOT.ASSOCIATED(IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%z))&
     ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%z(ns,mnmax))
  IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%z = TRANSPOSE(zmns)
  !allocate(IDS_EQ_OUT%time_slice(itime)%profiles_2d(2))
  !IDS_EQ_OUT%time_slice(itime)%profiles_2d(1)%grid_type%index=1 ! rectangular
  !IDS_EQ_OUT%time_slice(itime)%profiles_2d(2)%grid_type%index=56 ! inverse_rhotor_polar_fourier
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_2d(2)%grid%dim1(ns))
  !ALLOCATE(IDS_EQ_OUT%time_slice(itime)%profiles_2d(2)%grid%dim2(mpol+1))
  !IDS_EQ_OUT%time_slice(itime)%profiles_2d(2)%grid%dim1 ! inverse_rhotor_polar_fourier


  RETURN
END SUBROUTINE VMEC_EQOUT_IMAS

#endif
END MODULE VMEC_IMAS_MODULE
