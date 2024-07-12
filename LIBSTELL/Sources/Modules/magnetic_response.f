!*******************************************************************************
!>  @file magnetic_response.f
!>  @brief Contains module @ref magnetic_response
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref magnetic_response_class.
!*******************************************************************************

      MODULE magnetic_response
      USE compression
      USE stel_kinds, only: rprec
      USE profiler

      IMPLICIT NONE

!*******************************************************************************
!  magnetic_response module parameters
!*******************************************************************************
!  NOTE: When changing this version changes, the select statement in
!  magnetic_construct_netcdf of magnetic.f in the V3FIT needs to be updated. ???
!  FIXME Check if this is true. The versioning should probably be handeled
!  internally here.

!  Currnet version -------------------------------------------------------------
!>  Version for the MDSIG files. This version adds the point diagnostics.
      CHARACTER (len=*), PARAMETER ::                                          &
     &   magnetic_response_current = 'MRC 2015-04-27'
!  Previous versions -----------------------------------------------------------
!>  Version for the MDSIG files. This version adds the conducting shell.
!>  @depreicated
      CHARACTER (len=*), PARAMETER ::                                          &
     &   magnetic_response_20140928 = 'MRC 2014-09-28'

!>  Length for strings.
      INTEGER, PARAMETER :: magnetic_response_len = 80

!>  Bit position for the use coil response flag.
      INTEGER, PARAMETER :: magnetic_response_use_coil_flag   = 1
!>  Bit position for the force coil response flag.
      INTEGER, PARAMETER :: magnetic_response_use_plasma_flag = 2
!>  Bit position for the stellerator symmetry flag.
      INTEGER, PARAMETER :: magnetic_response_stell_sym_flag  = 3
!>  Bit position for the use conducting shell flag.
      INTEGER, PARAMETER :: magnetic_response_use_shell_flag  = 4
!>  Bit position for the use conducting shell flag.
      INTEGER, PARAMETER :: magnetic_response_is_point_flag   = 5

!  NETCDF file variable names --------------------------------------------------
!>  NETCDF configureation flags.
      CHARACTER (len=*), PARAMETER :: nc_flags = 'mddc_mrf_flags'

!  Identification Variables ----------------------------------------------------
!>  NETCDF code name variable.
      CHARACTER (len=*), PARAMETER :: nc_name = 'mddc_mrf_code_name'
!>  NETCDF version variable.
      CHARACTER (len=*), PARAMETER :: nc_version =                             &
     &   'mddc_mrf_code_version'
!>  NETCDF date run variable.
      CHARACTER (len=*), PARAMETER :: nc_date = 'mddc_mrf_date_run'
!>  NETCDF coil identifier variable.
      CHARACTER (len=*), PARAMETER :: nc_coil_id =                             &
     &   'mddc_mrf_field_coils_id'

!  Coil Responce Function Variables --------------------------------------------
!>  NETCDF number of field coils variable.
      CHARACTER (len=*), PARAMETER :: nc_n_field_cg =                          &
     &   'mddc_mrf_n_field_cg'
!>  NETCDF mutual inductance variable.
      CHARACTER (len=*), PARAMETER :: nc_inductance =                          &
     &   'mddc_mrf_rdiag_coilg_1'
!>  NETCDF current scale variable.
      CHARACTER (len=*), PARAMETER :: nc_current_scale =                       &
     &   'mddc_mrf_extcur_mg'

!  Plasma Response Grid Variables ----------------------------------------------
!>  NETCDF number torodial planes variable.
      CHARACTER (len=*), PARAMETER :: nc_num_t = 'mddc_mrf_kp'
!>  NETCDF number radial grid points variable.
      CHARACTER (len=*), PARAMETER :: nc_num_r = 'mddc_mrf_ir'
!>  NETCDF number vertical grid points variable.
      CHARACTER (len=*), PARAMETER :: nc_num_z = 'mddc_mrf_jz'
!>  NETCDF minimum radial grid variable.
      CHARACTER (len=*), PARAMETER :: nc_rmin = 'mddc_mrf_rmin'
!>  NETCDF maximum radial grid variable.
      CHARACTER (len=*), PARAMETER :: nc_rmax = 'mddc_mrf_rmax'
!>  NETCDF minimum vertical grid variable.
      CHARACTER (len=*), PARAMETER :: nc_zmin = 'mddc_mrf_zmin'
!>  NETCDF maximum vertical grid variable.
      CHARACTER (len=*), PARAMETER :: nc_zmax = 'mddc_mrf_zmax'
!>  NETCDF number of field periods variable.
      CHARACTER (len=*), PARAMETER :: nc_n_field_periods =                     &
     &   'mddc_mrf_n_field_periods'

!>  NETCDF stell symmetric variable.
!>  @depricated
      CHARACTER (len=*), PARAMETER :: nc_stell_sym =                           &
     &   'mddc_mrf_lstell_sym'

!>  NETCDF radial response grid variable.
      CHARACTER (len=*), PARAMETER :: nc_a_r = 'mddc_mrf_a_r'
!>  NETCDF toroidal response grid variable.
      CHARACTER (len=*), PARAMETER :: nc_a_f = 'mddc_mrf_a_f'
!>  NETCDF vertical response grid variable.
      CHARACTER (len=*), PARAMETER :: nc_a_z = 'mddc_mrf_a_z'

!  Conducting Shell Response Arrays --------------------------------------------
!>  NETCDF use conducting shell variable.
!>  @depricated
      CHARACTER (len=*), PARAMETER :: nc_use_shell =                           &
     &   'mddc_mrf_use_con_shell'
!>  NETCDF number torodial shell grid points variable.
      CHARACTER (len=*), PARAMETER :: nc_num_t_shell =                         &
     &   'mddc_mrf_kp_shell'

!>  NETCDF radial response shell grid variable.
      CHARACTER (len=*), PARAMETER :: nc_a_s_r = 'mddc_mrf_a_s_r'
!>  NETCDF toroidal response shell grid variable.
      CHARACTER (len=*), PARAMETER :: nc_a_s_f = 'mddc_mrf_a_s_f'
!>  NETCDF vertical response shell grid variable.
      CHARACTER (len=*), PARAMETER :: nc_a_s_z = 'mddc_mrf_a_s_z'

!  Point Diagnostic ------------------------------------------------------------
!>  NETCDF use conducting shell variable.
      CHARACTER (len=*), PARAMETER :: nc_position =                            &
     &   'mddc_mrf_position'
!>  NETCDF number torodial shell grid points variable.
      CHARACTER (len=*), PARAMETER :: nc_direction =                           &
     &   'mddc_mrf_direction'

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) magnetic base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a magnetic signal response function.
!-------------------------------------------------------------------------------
      TYPE magnetic_response_class
!>  Flag variables stored as a bit packed integer.
         INTEGER                               :: flags

!  Identification Variables ----------------------------------------------------
!>  Name of the code which computed the responses.
         CHARACTER (len=magnetic_response_len) :: name
!>  Version of the code.
         CHARACTER (len=magnetic_response_len) :: version
!>  Data and time the code was run.
         CHARACTER (len=magnetic_response_len) :: date
!>  Idenifier of field-coils.
         CHARACTER (len=magnetic_response_len) :: coil_id

!  Coil Responce Function Variables --------------------------------------------
!>  Number of field-coil groups.
         INTEGER                             :: n_field_cg
!>  Array of diagnostic field-coil-group responses.
         REAL (rprec), DIMENSION(:), POINTER :: inductance => null()
!>  Array of external current scales.
         REAL (rprec), DIMENSION(:), POINTER :: current_scale => null()

!  Plasma Response Grid Variables ----------------------------------------------
!>  Number of phi planes per field period in plasma grid.
         INTEGER                            :: num_t
!>  Minimum R for plasma grid.
         REAL (rprec)                       :: rmin
!>  Maximum R for plasma grid.
         REAL (rprec)                       :: rmax
!>  Minimum z for plasma grid.
         REAL (rprec)                       :: zmin
!>  Maximum z for plasma grid.
         REAL (rprec)                       :: zmax

!>  Number of field periods.
         INTEGER                            :: n_field_periods

!>  Number of phi planes per field period in the conducting shell grid.
         INTEGER                            :: num_t_shell

!  Plasma Response Arrays ------------------------------------------------------
!>  R component of plasma response function.
         TYPE (compression_pointer), DIMENSION(:), POINTER ::                  &
     &      a_r => null()
!>  Phi component of plasma response function.
         TYPE (compression_pointer), DIMENSION(:), POINTER ::                  &
     &      a_f => null()
!>  Z component of plasma response function.
         TYPE (compression_pointer), DIMENSION(:), POINTER ::                  &
     &      a_z => null()

!  Conducting Shell Response Arrays --------------------------------------------
!>  R component of the conducting shell response function.
         TYPE (compression_class), POINTER :: a_s_r => null()
!>  Phi component of the conducting shell response function.
         TYPE (compression_class), POINTER :: a_s_f => null()
!>  Z component of the conducting shell response function.
         TYPE (compression_class), POINTER :: a_s_z => null()

!  Point Diagnostic Response ---------------------------------------------------
!  Use the mutual inductance array to store the vacuum field response.
!>  Position of the bfield measurement point.
         REAL (rprec), DIMENSION(3) :: position
!>  Position of the bfield measurement point.
         REAL (rprec), DIMENSION(3) :: direction
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref magnetic_response_class types using
!>  @ref magnetic_response_construct_new,                                      &
!>  @ref magnetic_response_construct_point or                                  &
!>  @ref magnetic_response_construct_netcdf.
!-------------------------------------------------------------------------------
      INTERFACE magnetic_response_construct
         MODULE PROCEDURE magnetic_response_construct_new,                     &
     &                    magnetic_response_construct_point,                   &
     &                    magnetic_response_construct_netcdf
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref magnetic_response_class object.
!>
!>  Allocates memory and initializes a @ref magnetic_response_class object.
!>
!>  @param[in] name            Name of the code which computed the responses.
!>  @param[in] version         Version of the code.
!>  @param[in] date            Data and time the code was run.
!>  @param[in] coil_id         Idenifier of field-coils.
!>  @param[in] inductance      Array of diagnostic field-coil-group responses.
!>  @param[in] current_scale   Array of external current scales.
!>  @param[in] num_t           Number of phi planes per field period in plasma
!>                             grid.
!>  @param[in] num_t_shell     Number of phi planes per field period in the
!>                             conducting shell grid.
!>  @param[in] rmin            Minimum R for plasma grid.
!>  @param[in] rmax            Maximum R for plasma grid.
!>  @param[in] zmin            Minimum z for plasma grid.
!>  @param[in] zmax            Maximum z for plasma grid.
!>  @param[in] n_field_periods Number of field periods.
!>  @param[in] stell_sym       Use stellarator symmetry.
!>  @param[in] a_r             R component of plasma response function.
!>  @param[in] a_f             Phi component of plasma response function.
!>  @param[in] a_z             Z component of plasma response function.
!>  @param[in] a_s_r           R component of the conducting shell response
!>                             function.
!>  @param[in] a_s_f           Phi component of the conducting shell response
!>                             function.
!>  @param[in] a_s_z           Z component of the conducting shell response
!>                             function.
!>  @param[in] svd_cut_off     Cutoff on singular values for data compression.
!>  @returns A pointer to a constructed @ref magnetic_class object.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_response_construct_new(name, date,                     &
     &                                         coil_id, inductance,            &
     &                                         current_scale,                  &
     &                                         num_t, num_t_shell,             &
     &                                         rmin, rmax, zmin, zmax,         &
     &                                         n_field_periods,                &
     &                                         stell_sym, a_r, a_f, a_z,       &
     &                                         a_s_r, a_s_f, a_s_z,            &
     &                                         svd_cut_off)
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_response_class), POINTER                                  &
     &   :: magnetic_response_construct_new
      CHARACTER (len=*), INTENT(in)           :: name
      CHARACTER (len=*), INTENT(in)           :: date
      CHARACTER (len=*), INTENT(in)           :: coil_id
      REAL (rprec), DIMENSION(:), POINTER     :: inductance
      REAL (rprec), DIMENSION(:), POINTER     :: current_scale
      INTEGER, INTENT(in)                     :: num_t
      INTEGER, INTENT(in)                     :: num_t_shell
      REAL (rprec), INTENT(in)                :: rmin
      REAL (rprec), INTENT(in)                :: rmax
      REAL (rprec), INTENT(in)                :: zmin
      REAL (rprec), INTENT(in)                :: zmax
      INTEGER, INTENT(in)                     :: n_field_periods
      LOGICAL, INTENT(in)                     :: stell_sym
      REAL (rprec), DIMENSION(:,:,:), POINTER :: a_r
      REAL (rprec), DIMENSION(:,:,:), POINTER :: a_f
      REAL (rprec), DIMENSION(:,:,:), POINTER :: a_z
      REAL (rprec), DIMENSION(:,:), POINTER   :: a_s_r
      REAL (rprec), DIMENSION(:,:), POINTER   :: a_s_f
      REAL (rprec), DIMENSION(:,:), POINTER   :: a_s_z
      REAL (rprec), INTENT(in)                :: svd_cut_off

!  local variables
      INTEGER                                 :: phi
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(magnetic_response_construct_new)

      magnetic_response_construct_new%flags = 0

!  Identification Variables ----------------------------------------------------
      magnetic_response_construct_new%name = name
      magnetic_response_construct_new%version =                                &
     &   magnetic_response_current
      magnetic_response_construct_new%date = date
      magnetic_response_construct_new%coil_id = coil_id

!  Coil Responce Function Variables --------------------------------------------
      IF (ASSOCIATED(inductance)) THEN
         magnetic_response_construct_new%flags =                               &
     &      IBSET(magnetic_response_construct_new%flags,                       &
     &            magnetic_response_use_coil_flag)

         ALLOCATE(magnetic_response_construct_new%inductance(                  &
     &               SIZE(inductance)))
         magnetic_response_construct_new%inductance = inductance
         magnetic_response_construct_new%n_field_cg = SIZE(inductance)
      END IF

      IF (ASSOCIATED(current_scale)) THEN
         ALLOCATE(magnetic_response_construct_new%current_scale(               &
     &               SIZE(current_scale)))
         magnetic_response_construct_new%current_scale = current_scale
      END IF

!  Plasma Response Grid Variables ----------------------------------------------
      magnetic_response_construct_new%num_t = num_t
      magnetic_response_construct_new%num_t_shell = num_t_shell
      magnetic_response_construct_new%rmin = rmin
      magnetic_response_construct_new%rmax = rmax
      magnetic_response_construct_new%zmin = zmin
      magnetic_response_construct_new%zmax = zmax
      magnetic_response_construct_new%n_field_periods = n_field_periods

      IF (stell_sym) THEN
         magnetic_response_construct_new%flags =                               &
     &      IBSET(magnetic_response_construct_new%flags,                       &
     &            magnetic_response_stell_sym_flag)
      END IF

!  Plasma Response Arrays ------------------------------------------------------
      IF (ASSOCIATED(a_r)) THEN
         magnetic_response_construct_new%flags =                               &
     &      IBSET(magnetic_response_construct_new%flags,                       &
     &            magnetic_response_use_plasma_flag)

         CALL assert(ASSOCIATED(a_f), 'a_f response function not ' //          &
     &                                'allocated')
         CALL assert(ASSOCIATED(a_z), 'a_z response function not ' //          &
     &                                'allocated')

         ALLOCATE(magnetic_response_construct_new%a_r(SIZE(a_r, 3)))
         ALLOCATE(magnetic_response_construct_new%a_f(SIZE(a_r, 3)))
         ALLOCATE(magnetic_response_construct_new%a_z(SIZE(a_r, 3)))

         DO phi = 1, SIZE(a_r, 3)
            magnetic_response_construct_new%a_r(phi)%p =>                      &
     &         compression_construct(a_r(:,:,phi), svd_cut_off)
            magnetic_response_construct_new%a_f(phi)%p =>                      &
     &         compression_construct(a_f(:,:,phi), svd_cut_off)
            magnetic_response_construct_new%a_z(phi)%p =>                      &
     &         compression_construct(a_z(:,:,phi), svd_cut_off)
         END DO
      END IF

!  Conducting Shell Response Arrays --------------------------------------------
      IF (ASSOCIATED(a_s_r)) THEN
         magnetic_response_construct_new%flags =                               &
     &      IBSET(magnetic_response_construct_new%flags,                       &
     &            magnetic_response_use_shell_flag)

         CALL assert(ASSOCIATED(a_s_f), 'a_s_f response function ' //          &
     &                                  'not allocated')
         CALL assert(ASSOCIATED(a_s_z), 'a_s_z response function ' //          &
     &                                  'not allocated')

         magnetic_response_construct_new%a_s_r =>                              &
     &      compression_construct(a_s_r, svd_cut_off)
         magnetic_response_construct_new%a_s_f =>                              &
     &      compression_construct(a_s_f, svd_cut_off)
         magnetic_response_construct_new%a_s_z =>                              &
     &      compression_construct(a_s_z, svd_cut_off)
      END IF

      CALL profiler_set_stop_time('magnetic_response_construct_new',           &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref magnetic_response_class object.
!>
!>  Allocates memory and initializes a @ref magnetic_response_class object.
!>
!>  @param[in] name          Name of the code which computed the responses.
!>  @param[in] version       Version of the code.
!>  @param[in] date          Data and time the code was run.
!>  @param[in] coil_id       Idenifier of field-coils.
!>  @param[in] position      Position of the point measurement.
!>  @param[in] direction     Direction of the point measurement.
!>  @param[in] vacuum        Array of diagnostic vacuum responses.
!>  @param[in] current_scale Array of external current scales.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_response_construct_point(name, date,                   &
     &                                           coil_id, position,            &
     &                                           direction, vacuum,            &
     &                                           current_scale)
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_response_class), POINTER                                  &
     &   :: magnetic_response_construct_point
      CHARACTER (len=*), INTENT(in)       :: name
      CHARACTER (len=*), INTENT(in)       :: date
      CHARACTER (len=*), INTENT(in)       :: coil_id
      REAL (rprec), DIMENSION(3)          :: position
      REAL (rprec), DIMENSION(3)          :: direction
      REAL (rprec), DIMENSION(:), POINTER :: vacuum
      REAL (rprec), DIMENSION(:), POINTER :: current_scale

!  local variables
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(magnetic_response_construct_point)

      magnetic_response_construct_point%flags = 0

      magnetic_response_construct_point%flags =                                &
     &   IBSET(magnetic_response_construct_point%flags,                        &
     &         magnetic_response_is_point_flag)

!  Identification Variables ----------------------------------------------------
      magnetic_response_construct_point%name = name
      magnetic_response_construct_point%version =                              &
     &   magnetic_response_current
      magnetic_response_construct_point%date = date
      magnetic_response_construct_point%coil_id = coil_id

!  Geometric Variables ---------------------------------------------------------
      magnetic_response_construct_point%position = position
      magnetic_response_construct_point%direction = direction

!  Coil Responce Function Variables --------------------------------------------
      IF (ASSOCIATED(vacuum)) THEN
         magnetic_response_construct_point%flags =                             &
     &      IBSET(magnetic_response_construct_point%flags,                     &
     &            magnetic_response_use_coil_flag)

         ALLOCATE(magnetic_response_construct_point%inductance(                &
     &               SIZE(vacuum)))
         magnetic_response_construct_point%inductance = vacuum
         magnetic_response_construct_point%n_field_cg = SIZE(vacuum)
      END IF

      IF (ASSOCIATED(current_scale)) THEN
         ALLOCATE(magnetic_response_construct_point%current_scale(             &
     &               SIZE(current_scale)))
         magnetic_response_construct_point%current_scale = current_scale
      END IF

      CALL profiler_set_stop_time('magnetic_response_construct_point',         &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref magnetic_response_class object.
!>
!>  Allocates memory and initializes a @ref magnetic_class object from a mdsig
!>  file.
!>
!>  @param[in] mdsig_iou   An instance of a the netcdf id of the open mdsig
!>                         file.
!>  @param[in] svd_cut_off Cutoff on singular values for data compression.
!>  @returns A pointer to a constructed @ref magnetic_response_class object.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_response_construct_netcdf(mdsig_iou,                   &
     &                                            svd_cut_off)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_response_class), POINTER ::                               &
     &   magnetic_response_construct_netcdf
      INTEGER, INTENT(in)                        :: mdsig_iou
      REAL (rprec), INTENT(in)                   :: svd_cut_off

!  local variables
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: temp_buffer
      LOGICAL                                    :: temp_logical
      INTEGER, DIMENSION(3)                      :: dim_lengths
      INTEGER                                    :: phi
      REAL (rprec)                               :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(magnetic_response_construct_netcdf)

!  Identification Variables ----------------------------------------------------
      CALL cdf_read(mdsig_iou, nc_name,                                        &
     &              magnetic_response_construct_netcdf%name)
      CALL cdf_read(mdsig_iou, nc_version,                                     &
     &              magnetic_response_construct_netcdf%version)
      CALL cdf_read(mdsig_iou, nc_date,                                        &
     &              magnetic_response_construct_netcdf%date)
      CALL cdf_read(mdsig_iou, nc_coil_id,                                     &
     &              magnetic_response_construct_netcdf%coil_id)

!  Check the version code verison. If the version is 'MRC 2015-04-27' or later,
!  there is a flag variable that indicates the contents of the file.
      magnetic_response_construct_netcdf%flags = 0

      SELECT CASE (TRIM(magnetic_response_construct_netcdf%version))

         CASE (magnetic_response_current)
!  This version adds the new point diagnostic. The format of the netcdf mdsig
!  file has been modified to hold bit packed flags to determine what is or isn't
!  present in the netcdf file.
            CALL cdf_read(mdsig_iou, nc_flags,                                 &
     &                    magnetic_response_construct_netcdf%flags)

         CASE (magnetic_response_20140928)
!  This version added the option for a coducting shell. Options were stored as
!  logicals.
            CALL cdf_read(mdsig_iou, nc_use_shell, temp_logical)
            IF (temp_logical) THEN
               magnetic_response_construct_netcdf%flags =                      &
     &            IBSET(magnetic_response_construct_netcdf%flags,              &
     &                  magnetic_response_use_shell_flag)
            END IF

            CALL cdf_read(mdsig_iou, nc_stell_sym, temp_logical)
            IF (temp_logical) THEN
               magnetic_response_construct_netcdf%flags =                      &
     &            IBSET(magnetic_response_construct_netcdf%flags,              &
     &                  magnetic_response_stell_sym_flag)
            END IF

!  The mutual inductance and plasma response are expected.
            magnetic_response_construct_netcdf%flags =                         &
     &         IBSET(magnetic_response_construct_netcdf%flags,                 &
     &               magnetic_response_use_coil_flag)

            magnetic_response_construct_netcdf%flags =                         &
     &         IBSET(magnetic_response_construct_netcdf%flags,                 &
     &               magnetic_response_use_plasma_flag)

         CASE DEFAULT
!  This is the orginal format. Options were stored as logicals.
            CALL cdf_read(mdsig_iou, nc_stell_sym, temp_logical)
            IF (temp_logical) THEN
               magnetic_response_construct_netcdf%flags =                      &
     &            IBSET(magnetic_response_construct_netcdf%flags,              &
     &                  magnetic_response_stell_sym_flag)
            END IF

!  The mutual inductance and plasma response are expected.
            magnetic_response_construct_netcdf%flags =                         &
     &         IBSET(magnetic_response_construct_netcdf%flags,                 &
     &               magnetic_response_use_coil_flag)

            magnetic_response_construct_netcdf%flags =                         &
     &         IBSET(magnetic_response_construct_netcdf%flags,                 &
     &               magnetic_response_use_plasma_flag)

      END SELECT

!  Coil Responce Function Variables --------------------------------------------
      IF (BTEST(magnetic_response_construct_netcdf%flags,                      &
     &          magnetic_response_use_coil_flag)) THEN
         CALL cdf_read(mdsig_iou, nc_n_field_cg,                               &
     &                 magnetic_response_construct_netcdf%n_field_cg)

         ALLOCATE(magnetic_response_construct_netcdf%inductance(               &
     &               magnetic_response_construct_netcdf%n_field_cg))
         CALL cdf_read(mdsig_iou, nc_inductance,                               &
     &                 magnetic_response_construct_netcdf%inductance)

         ALLOCATE(magnetic_response_construct_netcdf%current_scale(            &
     &               magnetic_response_construct_netcdf%n_field_cg))
         CALL cdf_read(mdsig_iou, nc_current_scale,                            &
     &                 magnetic_response_construct_netcdf%current_scale)
      END IF

!  Plasma Response Grid Variables ----------------------------------------------
      IF (BTEST(magnetic_response_construct_netcdf%flags,                      &
     &          magnetic_response_use_plasma_flag)) THEN
         CALL cdf_read(mdsig_iou, nc_num_t,                                    &
     &                 magnetic_response_construct_netcdf%num_t)
         CALL cdf_read(mdsig_iou, nc_rmin,                                     &
     &                 magnetic_response_construct_netcdf%rmin)
         CALL cdf_read(mdsig_iou, nc_rmax,                                     &
     &                 magnetic_response_construct_netcdf%rmax)
         CALL cdf_read(mdsig_iou, nc_zmin,                                     &
     &                 magnetic_response_construct_netcdf%zmin)
         CALL cdf_read(mdsig_iou, nc_zmax,                                     &
     &                 magnetic_response_construct_netcdf%zmax)
         CALL cdf_read(mdsig_iou, nc_n_field_periods,                          &
     &           magnetic_response_construct_netcdf%n_field_periods)

         CALL cdf_inquire(mdsig_iou, nc_a_r, dim_lengths)
         ALLOCATE(temp_buffer(dim_lengths(1),                                  &
     &                        dim_lengths(2),                                  &
     &                        dim_lengths(3)))

         CALL cdf_read(mdsig_iou, nc_a_r, temp_buffer)
         ALLOCATE(magnetic_response_construct_netcdf%a_r(                      &
     &               dim_lengths(3)))
         DO phi = 1, dim_lengths(3)
            magnetic_response_construct_netcdf%a_r(phi)%p =>                   &
     &         compression_construct(temp_buffer(:,:,phi), svd_cut_off)
         END DO

         CALL cdf_read(mdsig_iou, nc_a_f, temp_buffer)
         ALLOCATE(magnetic_response_construct_netcdf%a_f(                      &
     &               dim_lengths(3)))
         DO phi = 1, dim_lengths(3)
            magnetic_response_construct_netcdf%a_f(phi)%p =>                   &
     &         compression_construct(temp_buffer(:,:,phi), svd_cut_off)
         END DO

         CALL cdf_read(mdsig_iou, nc_a_z, temp_buffer)
         ALLOCATE(magnetic_response_construct_netcdf%a_z(                      &
     &               dim_lengths(3)))
         DO phi = 1, dim_lengths(3)
            magnetic_response_construct_netcdf%a_z(phi)%p =>                   &
     &         compression_construct(temp_buffer(:,:,phi), svd_cut_off)
         END DO

         DEALLOCATE(temp_buffer)
      END IF

!  Conducting Shell Response Arrays --------------------------------------------
      IF (BTEST(magnetic_response_construct_netcdf%flags,                      &
     &          magnetic_response_use_shell_flag)) THEN
         CALL cdf_read(mdsig_iou, nc_num_t_shell,                              &
     &                 magnetic_response_construct_netcdf%num_t_shell)

         CALL cdf_inquire(mdsig_iou, nc_a_s_r, dim_lengths(1:2))
         ALLOCATE(temp_buffer(1,dim_lengths(1),dim_lengths(2)))

         CALL cdf_read(mdsig_iou, nc_a_s_r, temp_buffer(1,:,:))
         magnetic_response_construct_netcdf%a_s_r =>                           &
     &      compression_construct(temp_buffer(1,:,:), svd_cut_off)

         CALL cdf_read(mdsig_iou, nc_a_s_f, temp_buffer(1,:,:))
         magnetic_response_construct_netcdf%a_s_f =>                           &
     &      compression_construct(temp_buffer(1,:,:), svd_cut_off)

         CALL cdf_read(mdsig_iou, nc_a_s_z, temp_buffer(1,:,:))
         magnetic_response_construct_netcdf%a_s_z =>                           &
     &      compression_construct(temp_buffer(1,:,:), svd_cut_off)

         DEALLOCATE(temp_buffer)
      END IF

!  Point Diagnostic Response ---------------------------------------------------
      IF (BTEST(magnetic_response_construct_netcdf%flags,                      &
     &          magnetic_response_is_point_flag)) THEN
         CALL cdf_read(mdsig_iou, nc_position,                                 &
     &                 magnetic_response_construct_netcdf%position)

         CALL cdf_read(mdsig_iou, nc_direction,                                &
     &                 magnetic_response_construct_netcdf%direction)
      END IF

      CALL profiler_set_stop_time('magnetic_response_construct_netcdf',        &
     &                            start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref magnetic_response_class object.
!>
!>  Deallocates memory and uninitializes a @ref magnetic_response_class object.
!>
!>  @param[inout] this A @ref magnetic_response_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE magnetic_response_destruct(this)

!  Declare Arguments
      TYPE (magnetic_response_class), POINTER :: this

!  local variables
      INTEGER                                 :: phi

!  Start of executable code
      IF (ASSOCIATED(this%inductance)) THEN
         DEALLOCATE(this%inductance)
         this%inductance => null()
      END IF

      IF (ASSOCIATED(this%current_scale)) THEN
         DEALLOCATE(this%current_scale)
         this%current_scale => null()
      END IF

      IF (ASSOCIATED(this%a_r)) THEN
         DO phi = 1, SIZE(this%a_r)
            CALL compression_destruct(this%a_r(phi)%p)
         END DO
         DEALLOCATE(this%a_r)
         this%a_r => null()
      END IF

      IF (ASSOCIATED(this%a_f)) THEN
         DO phi = 1, SIZE(this%a_f)
            CALL compression_destruct(this%a_f(phi)%p)
         END DO
         DEALLOCATE(this%a_f)
         this%a_f => null()
      END IF

      IF (ASSOCIATED(this%a_z)) THEN
         DO phi = 1, SIZE(this%a_z)
            CALL compression_destruct(this%a_z(phi)%p)
         END DO
         DEALLOCATE(this%a_z)
         this%a_z => null()
      END IF

      IF (ASSOCIATED(this%a_s_r)) THEN
         CALL compression_destruct(this%a_s_r)
         this%a_s_r => null()
      END IF

      IF (ASSOCIATED(this%a_s_f)) THEN
         CALL compression_destruct(this%a_s_f)
         this%a_s_f => null()
      END IF

      IF (ASSOCIATED(this%a_s_z)) THEN
         CALL compression_destruct(this%a_s_z)
         this%a_s_z => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Checks if the stellarator symmetric flag is set.
!>
!>  This is a convenience method so that outside callers do not need to test
!>  bits directly.
!>
!>  @param[in] this A @ref magnetic_response_class instance.
!>  @returns True if the flag is set.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_response_is_stell_sym(this)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: magnetic_response_is_stell_sym
      TYPE (magnetic_response_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                               :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      magnetic_response_is_stell_sym =                                         &
     &   BTEST(this%flags, magnetic_response_stell_sym_flag)

      CALL profiler_set_stop_time('magnetic_response_is_stell_sym',            &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Checks if the plasma flag is set.
!>
!>  This is a convenience method so that outside callers do not need to test
!>  bits directly.
!>
!>  @param[in] this A @ref magnetic_response_class instance.
!>  @returns True if the flag is set.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_response_use_plasma(this)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: magnetic_response_use_plasma
      TYPE (magnetic_response_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                               :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      magnetic_response_use_plasma =                                           &
     &   BTEST(this%flags, magnetic_response_use_plasma_flag)

      CALL profiler_set_stop_time('magnetic_response_use_plasma',              &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Checks if the conducting shell flag is set.
!>
!>  This is a convenience method so that outside callers do not need to test
!>  bits directly.
!>
!>  @param[in] this A @ref magnetic_response_class instance.
!>  @returns True if the flag is set.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_response_use_shell(this)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: magnetic_response_use_shell
      TYPE (magnetic_response_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                               :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      magnetic_response_use_shell =                                            &
     &   BTEST(this%flags, magnetic_response_use_shell_flag)

      CALL profiler_set_stop_time('magnetic_response_use_shell',               &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Checks if the coil response flag is set.
!>
!>  This is a convenience method so that outside callers do not need to test
!>  bits directly.
!>
!>  @param[in] this A @ref magnetic_response_class instance.
!>  @returns True if the flag is set.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_response_use_coil(this)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: magnetic_response_use_coil
      TYPE (magnetic_response_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                               :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      magnetic_response_use_coil =                                             &
     &   BTEST(this%flags, magnetic_response_use_coil_flag)

      CALL profiler_set_stop_time('magnetic_response_use_coil',                &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Clears the coil response bit.
!>
!>  This is a convenience method so that outside callers do not need to clear
!>  bits directly.
!>
!>  @param[in] this A @ref magnetic_response_class instance.
!>  @returns True if the flag is set.
!-------------------------------------------------------------------------------
      SUBROUTINE magnetic_response_clr_use_coil(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_response_class), INTENT(inout) :: this

!  local variables
      REAL (rprec)                                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      this%flags = IBCLR(this%flags, magnetic_response_use_coil_flag)

      CALL profiler_set_stop_time('magnetic_response_clr_use_coil',            &
     &                            start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Checks if the coil response flag is set.
!>
!>  This is a convenience method so that outside callers do not need to test
!>  bits directly.
!>
!>  @param[in] this A @ref magnetic_response_class instance.
!>  @returns True if the flag is set.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_response_is_point(this)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: magnetic_response_is_point
      TYPE (magnetic_response_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                               :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      magnetic_response_is_point =                                             &
     &   BTEST(this%flags, magnetic_response_is_point_flag)

      CALL profiler_set_stop_time('magnetic_response_is_point',                &
     &                            start_time)

      END FUNCTION

!*******************************************************************************
!  NETCDF SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Defines the variables for the NETCDF file.
!>
!>  Write out the magnetic response to an mdsig file.
!>
!>  @param[in] this      A @ref magnetic_response_class instance.
!>  @param[in] mdsig_iou An instance of a the netcdf id of the open mdsig file.
!-------------------------------------------------------------------------------
      SUBROUTINE magnetic_response_define(this, mdsig_iou)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_response_class), INTENT(in)  :: this
      INTEGER, INTENT(in)                         :: mdsig_iou

!  local variables
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: temp_buffer
      INTEGER                                     :: phi
      REAL (rprec)                                :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Define Variables ------------------------------------------------------------
      CALL cdf_define(mdsig_iou, nc_flags, this%flags)

!  Identification Variables ----------------------------------------------------
      CALL cdf_define(mdsig_iou, nc_name, this%name)
      CALL cdf_define(mdsig_iou, nc_version, this%version)
      CALL cdf_define(mdsig_iou, nc_date, this%date)
      CALL cdf_define(mdsig_iou, nc_coil_id, this%coil_id)

!  Coil Responce Function Variables --------------------------------------------
      IF (BTEST(this%flags, magnetic_response_use_coil_flag)) THEN
         CALL cdf_define(mdsig_iou, nc_n_field_cg, this%n_field_cg)
         CALL cdf_define(mdsig_iou, nc_inductance, this%inductance)
         CALL cdf_define(mdsig_iou, nc_current_scale,                          &
     &                   this%current_scale)
      END IF

!  Plasma Response Grid Variables ----------------------------------------------
      IF (BTEST(this%flags, magnetic_response_use_plasma_flag)) THEN
         CALL cdf_define(mdsig_iou, nc_num_t, this%num_t)
         CALL cdf_define(mdsig_iou, nc_rmin, this%rmin)
         CALL cdf_define(mdsig_iou, nc_rmax, this%rmax)
         CALL cdf_define(mdsig_iou, nc_zmin, this%zmin)
         CALL cdf_define(mdsig_iou, nc_zmax, this%zmax)
         CALL cdf_define(mdsig_iou, nc_n_field_periods,                        &
     &                   this%n_field_periods)

!  Decompress a single response function phi plane to get the r z dimensions.
!  Then define the netcdf variables based off a temp buffer. All response
!  function directions should have the same dimensions.
         ALLOCATE(temp_buffer(compression_get_dimension1(this%a_r(1)%p),       &
     &                        compression_get_dimension2(this%a_r(1)%p),       &
     &                        SIZE(this%a_r)))

         CALL cdf_define(mdsig_iou, nc_a_r, temp_buffer)
         CALL cdf_define(mdsig_iou, nc_a_f, temp_buffer)
         CALL cdf_define(mdsig_iou, nc_a_z, temp_buffer)

         DEALLOCATE(temp_buffer)
      END IF

!  Conducting Shell Response Arrays --------------------------------------------
      IF (BTEST(this%flags, magnetic_response_use_shell_flag)) THEN
         CALL cdf_define(mdsig_iou, nc_num_t_shell, this%num_t_shell)

!  Decompress a single conducting shell response function phi plane to get the
!  data buffer. All shell response function directions should have the same
!  dimensions.
         ALLOCATE(temp_buffer(1,                                               &
     &               compression_get_dimension1(this%a_s_r),                   &
     &               compression_get_dimension2(this%a_s_r)))

         CALL cdf_define(mdsig_iou, nc_a_s_r, temp_buffer(1,:,:))
         CALL cdf_define(mdsig_iou, nc_a_s_f, temp_buffer(1,:,:))
         CALL cdf_define(mdsig_iou, nc_a_s_z, temp_buffer(1,:,:))

         DEALLOCATE(temp_buffer)
      END IF

!  Point Diagnostic Response ---------------------------------------------------
      IF (BTEST(this%flags, magnetic_response_is_point_flag)) THEN
         CALL cdf_define(mdsig_iou, nc_position, this%position)
         CALL cdf_define(mdsig_iou, nc_direction, this%direction)
      END IF

      CALL profiler_set_stop_time('magnetic_response_define',                  &
     &                            start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write variables to the NETCDF file.
!>
!>  Write out the magnetic response to an mdsig file.
!>
!>  @param[in] this      A @ref magnetic_response_class instance.
!>  @param[in] mdsig_iou An instance of a the netcdf id of the open mdsig file.
!-------------------------------------------------------------------------------
      SUBROUTINE magnetic_response_write(this, mdsig_iou)
      USE ezcdf
      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_response_class), INTENT(in)  :: this
      INTEGER, INTENT(in)                         :: mdsig_iou

!  local variables
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: temp_buffer
      INTEGER                                     :: phi
      REAL (rprec)                                :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Define Variables ------------------------------------------------------------
      CALL cdf_write(mdsig_iou, nc_flags, this%flags)

!  Identification Variables ----------------------------------------------------
      CALL cdf_write(mdsig_iou, nc_name, this%name)
      CALL cdf_write(mdsig_iou, nc_version, this%version)
      CALL cdf_write(mdsig_iou, nc_date, this%date)
      CALL cdf_write(mdsig_iou, nc_coil_id, this%coil_id)

!  Coil Responce Function Variables --------------------------------------------
      IF (BTEST(this%flags, magnetic_response_use_coil_flag)) THEN
         CALL cdf_write(mdsig_iou, nc_n_field_cg, this%n_field_cg)
         CALL cdf_write(mdsig_iou, nc_inductance, this%inductance)
         CALL cdf_write(mdsig_iou, nc_current_scale, this%current_scale)
      END IF

!  Plasma Response Grid Variables ----------------------------------------------
      IF (BTEST(this%flags, magnetic_response_use_plasma_flag)) THEN
         CALL cdf_write(mdsig_iou, nc_num_t, this%num_t)
         CALL cdf_write(mdsig_iou, nc_rmin, this%rmin)
         CALL cdf_write(mdsig_iou, nc_rmax, this%rmax)
         CALL cdf_write(mdsig_iou, nc_zmin, this%zmin)
         CALL cdf_write(mdsig_iou, nc_zmax, this%zmax)
         CALL cdf_write(mdsig_iou, nc_n_field_periods,                         &
     &                  this%n_field_periods)

!  Need to decompress all the response function planes to and write the data
!  buffers to a temp buffer. Decompress the r direction first to get the r z
!  dimensions to allocate the temp buffer.
         CALL compression_decompress(this%a_r(1)%p)
         ALLOCATE(temp_buffer(SIZE(this%a_r(1)%p%data_buffer, 1),              &
     &                        SIZE(this%a_r(1)%p%data_buffer, 2),              &
     &                        SIZE(this%a_r)))

         temp_buffer(:,:,1) = this%a_r(1)%p%data_buffer
         CALL compression_cleanup(this%a_r(1)%p)
         DO phi = 2, SIZE(this%a_r)
            CALL compression_decompress(this%a_r(phi)%p)
            temp_buffer(:,:,phi) = this%a_r(phi)%p%data_buffer
            CALL compression_cleanup(this%a_r(phi)%p)
         END DO
         CALL cdf_write(mdsig_iou, nc_a_r, temp_buffer)

         DO phi = 1, SIZE(this%a_f)
            CALL compression_decompress(this%a_f(phi)%p)
            temp_buffer(:,:,phi) = this%a_f(phi)%p%data_buffer
            CALL compression_cleanup(this%a_f(phi)%p)
         END DO
         CALL cdf_write(mdsig_iou, nc_a_f, temp_buffer)

         DO phi = 1, SIZE(this%a_z)
            CALL compression_decompress(this%a_z(phi)%p)
            temp_buffer(:,:,phi) = this%a_z(phi)%p%data_buffer
            CALL compression_cleanup(this%a_z(phi)%p)
         END DO
         CALL cdf_write(mdsig_iou, nc_a_z, temp_buffer)

         DEALLOCATE(temp_buffer)
      END IF

!  Conducting Shell Response Arrays --------------------------------------------
      IF (BTEST(this%flags, magnetic_response_use_shell_flag)) THEN
         CALL cdf_write(mdsig_iou, nc_num_t_shell, this%num_t_shell)

!  Decompress each shell response function direction then write the data buffer.
         CALL compression_decompress(this%a_s_r)
         CALL cdf_write(mdsig_iou, nc_a_s_r, this%a_s_r%data_buffer)
         CALL compression_cleanup(this%a_s_r)

         CALL compression_decompress(this%a_s_f)
         CALL cdf_write(mdsig_iou, nc_a_s_f, this%a_s_f%data_buffer)
         CALL compression_cleanup(this%a_s_f)

         CALL compression_decompress(this%a_s_z)
         CALL cdf_write(mdsig_iou, nc_a_s_z, this%a_s_z%data_buffer)
         CALL compression_cleanup(this%a_s_z)
      END IF

!  Point Diagnostic Response ---------------------------------------------------
      IF (BTEST(this%flags, magnetic_response_is_point_flag)) THEN
         CALL cdf_write(mdsig_iou, nc_position, this%position)
         CALL cdf_write(mdsig_iou, nc_direction, this%direction)
      END IF

      CALL profiler_set_stop_time('magnetic_response_write', start_time)

      END SUBROUTINE

      END MODULE
