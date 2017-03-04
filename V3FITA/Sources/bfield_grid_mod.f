!     SPH: INTEGER(iprec) -> INTEGER
      MODULE bfield_grid_mod
!--------------------------------------------------------------------------
! 
! FUNCTION: Module containing routines necessary to compute the line-integrated
!           Faraday rotation angle and the line-integrated polarization phase shift.
!
! modified by JS 11/8/06 to use VMEC instead of AJAX for the flux coord transformations
!------------------------------------------------------------------------------

      USE stel_kinds
      USE stel_constants
      USE ajax_mod
      USE safe_open_mod
      USE v3_utilities     ! contains subroutine assert_eq
      USE v3read_wout      ! reads equilibrium from a wout file into AJAX (includes READEQ_FROM_WOUT)
      USE math_utilities   ! module containing vector manipulation routines
      USE b_transform_mod  ! module containing B field transformation routines
      USE b_transform_vmec  ! module containing VMEC B field transformation routines
      IMPLICIT NONE

      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: flux_surface     ! array of pts along each flux surface

      INTEGER, PARAMETER :: USE_VMEC = 0
      INTEGER, PARAMETER :: USE_FS = 1

      INTEGER :: transform_select = USE_VMEC      
!      INTEGER :: transform_select = USE_FS      

      CONTAINS
!===================================================================================
      SUBROUTINE B_GRID_DRIVER()
!--------------------------------------------------------------------------------
!
! Function:  Top-level "driver" subroutine for module "bfield_grid_mod"
!----------------------------------------------------------------------------------
      USE faraday_mod   ! contains conversion routine TRACK_TRANS2FLX
      USE b_transform_vmec, ONLY: VMEC_B_INIT, VMEC_CYL2FLX
      USE safe_open_mod
      IMPLICIT NONE

      INTEGER :: i, j, iflag, ierr_vmec
      INTEGER :: istat
      INTEGER :: m, n
      INTEGER :: ngridpts_R      ! number of grid POINTS along R
      INTEGER :: ngridpts_Z      ! number of grid POINTS along Z
      INTEGER :: m_R             ! R dimension for Bcar array
      INTEGER :: n_Z             ! Z dimension for Bcar array
      INTEGER       :: ns_rho        ! number of radial grid pts in VMEC arrays
!      REAL(rprec)       :: rho
      REAL(rprec)       :: stepsize_R
      REAL(rprec)       :: stepsize_Z
      REAL(rprec)       :: Rmin_temp, Rmax_temp
      REAL(rprec)       :: Zmin_temp, Zmax_temp


      REAL(rprec), DIMENSION(3) :: dims
      REAL(rprec), DIMENSION(3) :: cyl_vec
      REAL(rprec), DIMENSION(3) :: x_cyl_temp
      REAL(rprec), DIMENSION(3) :: xflux_ajax         ! position vector in AJAX definition of flux coords
      REAL(rprec), DIMENSION(3) :: xflux_vmec         ! position vector from VMEC calc of flux coords
      REAL(rprec), DIMENSION(3) :: xflux, xflux2
      REAL(rprec), DIMENSION(6) :: g_cyl
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: B_car
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: rho_rz      ! VMEC radial flux coord array
      REAL(rprec), DIMENSION(3) :: B_car_temp
      REAL(rprec), DIMENSION(3) :: B_temp, B_temp2
      REAL(rprec) :: testarray(5)
      CHARACTER(len=80)  :: nml_infile
      CHARACTER(len=80)  :: default_nml_infile
      CHARACTER(len=80)  :: message

      CHARACTER(len=80) :: mychar1
      CHARACTER(len=*), PARAMETER :: subname = 'B_GRID_DRIVER: '
      LOGICAL        :: QUERY_USER_FOR_NML_FILE
      INTEGER :: iou_wout




!............NML variables...................................................!
      INTEGER :: nsegments_R       ! number of grid SEGMENTS along R
      INTEGER :: nsegments_Z       ! number of grid SEGMENTS along Z
      INTEGER :: nsegments_theta   ! number of grid SEGMENTS along poloidal flux coord theta
      REAL(rprec)    :: Rmin, Rmax
      REAL(rprec)    :: Zmin, Zmax
      REAL(rprec)    :: phi
      CHARACTER(len=80) :: wout_infile
      CHARACTER(len=80) :: cdf_outfile
   
      NAMELIST /mylist/ Rmin, Rmax, Zmin, Zmax,                                &
     &                  nsegments_R, nsegments_Z,  phi,                        &
     &                  nsegments_theta, wout_infile, cdf_outfile

      write(*,*) '===================================================='    
      write(*,*) ' '
!      write(*,*) 'Hi.  Now in ', subname

!      QUERY_USER_FOR_NML_FILE = .false.
      QUERY_USER_FOR_NML_FILE = .true.

!..........Create a name for a default nml file that can be used in place of user input........!
      default_nml_infile = 'bgrid_in.nml'

      if ( QUERY_USER_FOR_NML_FILE) then
        message = "Enter nml filename or select 'd' for default file"
        write(*,*) message
        READ(*,*) nml_infile

        if( nml_infile == 'd') then
          write(*,*) 'Using default nml file: ', default_nml_infile
          nml_infile = default_nml_infile
        end if
      else
!......."hard-wire" in the name of the default input nml file containing the other values.........!
        nml_infile = default_nml_infile
      end if


!..........open the input nml file..........................!
      CALL safe_open(iou_wout, istat, TRIM(nml_infile),                 &
     &               'old', 'formatted')
      CALL assert_eq(0,istat,subname //                                       &
     &   ' Safe_open of input NML file failed')

      READ( UNIT = iou_wout, NML = mylist)
      CLOSE(iou_wout)

      write(*,*) subname, 'Rmin = ', Rmin
      write(*,*) subname, 'Rmax = ', Rmax
      write(*,*) subname, 'Zmin = ', Zmin
      write(*,*) subname, 'Zmax = ', Zmax
      write(*,*) subname, 'm = ', nsegments_R
      write(*,*) subname, 'n = ', nsegments_Z
      write(*,*) subname, 'phi = ', phi
      write(*,*) subname, 'nsegments_theta = ', nsegments_theta
      write(*,*) subname, 'input wout = ', wout_infile
      write(*,*) subname, 'cdf_outfile = ', cdf_outfile
 
!...........fill AJAX stuff in memory from the input wout file.................!
      call DEAL_WITH_WOUT(wout_infile, ns_rho)

!........At this point , AJAX stuff for the wout should be in memory................!
!......also note that wout variables are now globally available via a USE call to read_wout_mod...!

!........initialize the VMEC flux coord transformation routines..........!
      call VMEC_B_INIT


!.......fill the array of (R,Z) points that describe each flux surface...........!
      call FILL_FLUX_SURFACE_ARRAY(ns_rho, phi, nsegments_theta)

!!.......write out the flux_surface array to std output to verify that it filled correctly....!
!      do i = 1, ns_rho
!        write(*,*) 'flux suface number = ', i
!        do j = 1, nsegments_theta
!          write(*,*) flux_surface(i,j,:) 
!         end do
!      end do


!.......convert the # of grid SEGMENTS into the # of grid POINTS (note that it is.....
!.......the grid point values (m_R, n_Z) that define the array size)................!
      ngridpts_R = nsegments_R + 1
      ngridpts_Z = nsegments_Z + 1


!.......allocate the B_car array.......................................!
      m_R = nsegments_R
      n_Z = nsegments_Z
      ALLOCATE( B_car(0:m_R, 0:n_Z, 3) )

!!.....DEBUGGING: Hard-wire cyl_vec values that give results inside the plasma. JS 12/16/05
!    (1.602,2.114, -0.117)   and (1.606, 2.104, -0.0738)
!      cyl_vec(1) = 1.602
!      cyl_vec(2) = 2.114
!      cyl_vec(3) = -0.117
!      Rmin_temp = 1.45
!      Rmax_temp = 1.606
!      Zmin_temp = -0.117
!      Zmax_temp = -0.0738


      stepsize_R = ABS(Rmax - Rmin) / REAL(nsegments_R)
      stepsize_Z = ABS(Zmax - Zmin) / REAL(nsegments_Z)

      do m = 0, m_R
        do n = 0, n_Z

!........Note that the x_cyl order in AJAX is the conventional (r,phi,z)...........!
          x_cyl_temp(1) = Rmin + m * stepsize_R
          x_cyl_temp(2) = phi
          x_cyl_temp(3) = Zmin + n * stepsize_Z
          call TRACK_TRANS2FLX(x_cyl_temp, xflux_ajax, iflag, message)

          if (iflag .ne. 0 ) then
            write(*,*) 'ERROR in TRACK_TRANS2FLX! MESSAGE IS: ', message
          end if

!...............use VMEC to compute flux coord position vector for debugging.  JS 12/4/06
          call VMEC_CYL2FLX(x_cyl_temp, xflux_vmec, ierr_vmec)

       
!           write(*,*) 'm = ', m, 'n = ', n
!           write(*,*) subname, 'x_cyl_temp = ', x_cyl_temp
!           write(*,*) subname, 'xflux_ajax = ', xflux_ajax

!!................verify that the current grid point is INSIDE the plasma................!
!!.........(and set upper limit to 0.9801 instead of 1.0 to avoid weird plasma "edge effects")......!
!          if ( xflux_ajax(1) .gt. 0.9801 ) then

!...............replace AJAX plasma in/out check with FS and/or VMEC check.  JS 12/2/06....!
!          if ( INVALID_RHO(x_cyl_temp) ) then
          if ( INVALID_RHO_VMEC(x_cyl_temp) ) then
            write(*,*) subname, 'WARNING: POINT IS OUTSIDE THE PLASMA!',       &
     &                 'xflux_ajax(1) = ',  xflux_ajax(1)

            B_car(m,n,:) = 0.0
          else
!            write(*,*) subname, 'xflux_ajax(1) = ', xflux_ajax(1)


!...............find B field at point x_cyl_temp & save to B_car array.............!
!             write(*,*) subname, 'now calling COMPUTE_B_CAR'
              call COMPUTE_B_CAR( x_cyl_temp, B_temp)
              B_car(m,n,:) = B_temp
          end if
        end do
      end do

!!..........write B_car arrays to std output for debugging (grouping by m [R] values)...............!
!      do i = 0,m_R
!          write(*,*) '----------------------------------------------'
!          write(*,*) subname, 'm value = ', i
!        do j = 0, n_Z
!          write(*,*) subname, 'B_car = ', B_car(i,j,:)
!        end do
!      end do
!
!!..........write B_car arrays to std output for debugging.(grouping by n [Z] values)...............!
!      do j = 0, n_Z
!          write(*,*) '----------------------------------------------'
!          write(*,*) subname, 'n value = ', j
!        do i = 0,m_R
!          write(*,*) subname, 'B_car = ', B_car(i,j,:)
!        end do
!      end do



!......write the results to an output netCDF file.........................!
      call OUTPUT_TO_CDF(Rmin, Rmax, Zmin, Zmax, nsegments_R,                       &
     &                   nsegments_Z, phi, B_car, ns_rho,                           &
     &                   nsegments_theta, cdf_outfile)

!........now read file back in to verify that the data was output to the .nc correctly....!
!      call READ_TEST_FOR_CDF(cdf_outfile)



      RETURN
      END SUBROUTINE B_GRID_DRIVER




!===================================================================================
      SUBROUTINE DEAL_WITH_WOUT(wout_infile, nr_rz)
!--------------------------------------------------------------------------------
!
! Function:  Retrieve necessary info from VMEC wout file
!
!----------------------------------------------------------------------------------
      USE read_wout_mod  ! module that holds VMEC variables in memory
      IMPLICIT NONE

!............dummy variables......................................!
      CHARACTER(len=*), INTENT(IN)  :: wout_infile
      INTEGER, INTENT(OUT) :: nr_rz ! number of radial grid pts in VMEC arrays

!............local variables......................................!
      INTEGER :: i, j, iflag
      INTEGER :: istat
      REAL(rprec), DIMENSION(3) :: dims
      REAL(rprec), DIMENSION(3) :: xflux
      REAL(rprec), DIMENSION(3) :: x_cyl
      CHARACTER(len=80)  :: message
      CHARACTER(len=*), PARAMETER :: subname = 'DEAL_WITH_WOUT: '


!      write(*,*) '===================================================='    
!      write(*,*) ' '
!      write(*,*) 'Hi.  Now in ', subname

 
!...........fill AJAX stuff in memory from the input wout file.................!
      call READEQ_FROM_WOUT(wout_infile)

!........At this point , AJAX stuff for the wout should be in memory.  Also note that wout variables..........!
!......are now globally available to this subroutine due the call to USE read_wout_mod ...!


!.......extract the number of radial grid points from global VMEC variables.............!
      dims(1:2) = shape(rmnc)
      nr_rz = dims(2)

      write(*,*) subname, 'nr_rz = ', nr_rz, 'ns = ', ns

      if ( nr_rz .ne. ns ) then
        message = 'WARNING! MISMATCH IN # OF VMEC RADIAL GRIDPTS!!'
        write(*,*) subname, message
      end if


      RETURN
      END SUBROUTINE DEAL_WITH_WOUT


!===================================================================================
      SUBROUTINE FILL_FLUX_SURFACE_ARRAY(ns_rho, input_phi,                     &
     &                                   nsegments_theta)
!--------------------------------------------------------------------------------
!
! Function:  Fill an array with cylinderical coordinate points along the flux surfaces
!
!----------------------------------------------------------------------------------
      USE ajax_mod          ! contains conversion routine AJAX_FLX2CYL
      USE b_transform_vmec, ONLY: VMEC_FLX2CYL_JS
      IMPLICIT NONE

!...........dummy variables...............................................!
      REAL(rprec), INTENT(IN)       :: input_phi           ! phi value specified in input NML file
      INTEGER, INTENT(IN)    :: ns_rho              ! number of radial grid pts in VMEC arrays
      INTEGER, INTENT(IN)    :: nsegments_theta     ! # of grid segments along the theta path

!...........local variables...............................................!
      INTEGER    :: i, j, iflag
      INTEGER    :: istat
      REAL(rprec)       :: theta        ! poloidal flux coordinate
      REAL(rprec)       :: delta_theta  ! grid step size in the poloidal direction

      REAL(rprec), DIMENSION(3) :: xflux
      REAL(rprec), DIMENSION(3) :: x_cyl
      REAL(rprec), DIMENSION(3) :: x_cyl_temp, xflx_temp
      REAL(rprec), DIMENSION(3) :: B_con_temp
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: rho_rz      ! VMEC radial flux coord array
      CHARACTER(len=80)  :: message

      CHARACTER(len=80) :: mychar1
      CHARACTER(len=*), PARAMETER :: subname='FILL_FLUX_SURFACE_ARRAY: '
      INTEGER :: iou_wout

      write(*,*) '===================================================='    
      write(*,*) ' '
!      write(*,*) 'Hi.  Now in ', subname



      ALLOCATE( rho_rz(1:ns_rho) )

!...........allocate output flux surface array........................!
      ALLOCATE( flux_surface(ns_rho, nsegments_theta, 2) )
 
!......fill rho_rz as an evenly spaced "full-grid", following Houlberg.........!
      do i = 1, ns_rho
        rho_rz(i) = SQRT( REAL(i-1) / REAL(ns_rho - 1) )
      end do

      if ( nsegments_theta .ne. 0 ) then
      delta_theta = (2.0 * PI) / REAL(nsegments_theta) 
      else
        delta_theta = 0.0
      end if



!..........convert each point along the flux surface contour into cylinderical coords..!

      do i = 1, ns_rho
        do j = 1, nsegments_theta

          xflux(1) = rho_rz(i)
          xflux(2) = (j-1) * delta_theta
          xflux(3) = input_phi

!!...........transform flux vector xflux() into cylin coord vector x_cyl() via AJAX
!          call AJAX_FLX2CYL(xflux, x_cyl, iflag, message)

!...........transform flux vector xflux() into cylin coord vector x_cyl() via VMEC
          call VMEC_FLX2CYL_JS(xflux, x_cyl, iflag)


!...........fill the array of (R,Z) points along each flux surface (note: R=1, Z=2)................!
          flux_surface(i,j,1) = x_cyl(1)
          flux_surface(i,j,2) = x_cyl(3)
         end do
      end do

      RETURN
      END SUBROUTINE FILL_FLUX_SURFACE_ARRAY


!===================================================================================
      SUBROUTINE OUTPUT_TO_CDF(Rmin, Rmax, Zmin, Zmax, nsegments_R,               &
     &                         nsegments_Z, phi, B_Car, ns_rho,                   &
     &                         nsegments_theta, cdf_outfile)
!--------------------------------------------------------------------------------
!
! Function:  Writes filled variables to an output netCDF file.
!
! Global inputs: flux_surface
! Global outputs: none
!
! CAVEAT:   There is no mechanism analogous to "safe_open" for the cdf_open subroutine
!         (at least none that I'm aware of...).  Therefore, it is *theoretically* possible
!         that the i/o unit number for the file could conflict with a previously opened
!          i/o unit.
!----------------------------------------------------------------------------------
      USE ezcdf
      IMPLICIT NONE

!.........passed dummy variables...............................!
      INTEGER :: nsegments_R
      INTEGER :: nsegments_Z
      INTEGER :: ns_rho
      INTEGER :: nsegments_theta
      REAL(rprec) :: Rmin, Rmax
      REAL(rprec) :: Zmin, Zmax
      REAL(rprec) :: phi
      REAL(rprec), DIMENSION(:,:,:) :: B_car
      CHARACTER(len=*) :: cdf_outfile

!...........local variables..............................!
      INTEGER           :: counter = 0
      INTEGER           :: i, j, nbeams, n_edge, iflag
      INTEGER           :: n_ip_units, istat
      REAL(rprec) :: testarray(5) 

      CHARACTER(len=80)  :: nml_infile

      CHARACTER(len=80) :: mychar1
      CHARACTER(len=80) :: mychar2, message
      CHARACTER(len=*), PARAMETER :: subname = 'OUTPUT_TO_CDF: '

      INTEGER :: iou_cdf = 132


!      write(*,*) 'Hi.  Now in ', subname

!........enter values for testarray.....!
      testarray(1) = 5.0
      testarray(2) = 6.0
      testarray(3) = 7.0
      testarray(4) = 8.0
      testarray(5) = 9.0

!      write(*,*) subname, 'testarray = ', testarray

!..........open the output netCDF file..........................!

      call cdf_open(iou_cdf, cdf_outfile, 'w', istat)
      CALL cdf_define(iou_cdf, 'Rmin', Rmin)
      CALL cdf_define(iou_cdf, 'Rmax', Rmax)
      CALL cdf_define(iou_cdf, 'Zmin', Zmin)
      CALL cdf_define(iou_cdf, 'Zmax', Zmax)
      CALL cdf_define(iou_cdf, 'phi', phi)
      CALL cdf_define(iou_cdf, 'nsegments_R', nsegments_R)
      CALL cdf_define(iou_cdf, 'nsegments_Z', nsegments_Z)
      CALL cdf_define(iou_cdf, 'ns_rho', ns_rho)
      CALL cdf_define(iou_cdf, 'nsegments_theta', nsegments_theta)
!      CALL cdf_define(iou_cdf, 'testarray', testarray)
      CALL cdf_define(iou_cdf, 'B_car', B_car)
      CALL cdf_define(iou_cdf, 'flux_surface', flux_surface)

      CALL cdf_write(iou_cdf, 'Rmin', Rmin)
      CALL cdf_write(iou_cdf, 'Rmax', Rmax)
      CALL cdf_write(iou_cdf, 'Zmin', Zmin)
      CALL cdf_write(iou_cdf, 'Zmax', Zmax)
      CALL cdf_write(iou_cdf, 'phi', phi)
      CALL cdf_write(iou_cdf, 'nsegments_R', nsegments_R)
      CALL cdf_write(iou_cdf, 'nsegments_Z', nsegments_Z)
      CALL cdf_write(iou_cdf, 'ns_rho', ns_rho)
      CALL cdf_write(iou_cdf, 'nsegments_theta', nsegments_theta)
!      CALL cdf_write(iou_cdf, 'testarray', testarray)
      CALL cdf_write(iou_cdf, 'B_car', B_car)
      CALL cdf_write(iou_cdf, 'flux_surface', flux_surface)
      call cdf_close(iou_cdf)

!      write(*,*) subname, 'Rmin = ', Rmin
!      write(*,*) subname, 'Rmax = ', Rmax
!      write(*,*) subname, 'Zmin = ', Zmin
!      write(*,*) subname, 'Zmax = ', Zmax
!      write(*,*) subname, 'm = ', nsegments_R
!      write(*,*) subname, 'n = ', nsegments_Z
!      write(*,*) subname, 'phi = ', phi
!      write(*,*) subname, 'cdf_outfile = ', cdf_outfile



      RETURN
      END SUBROUTINE OUTPUT_TO_CDF



!===================================================================================
      SUBROUTINE READ_TEST_FOR_CDF(cdf_infile)
!--------------------------------------------------------------------------------
!
! Function:  Test routine to verify that data was written to netCDF correctly
!
! CAVEAT:   There is no mechanism analogous to "safe_open" for the cdf_open subroutine
!         (at least none that I'm aware of...).  Therefore, it is *theoretically* possible
!         that the i/o unit number for the file could conflict with a previously opened
!          i/o unit.
!----------------------------------------------------------------------------------
      USE ezcdf
      IMPLICIT NONE

!.........passed dummy variables...............................!

      CHARACTER(len=*) :: cdf_infile

!...........local variables..............................!
      INTEGER           :: counter = 0
      INTEGER           :: i, j, nbeams, n_edge, iflag
      INTEGER           :: n_ip_units, istat

      INTEGER :: nsegments_R
      INTEGER :: nsegments_Z
      INTEGER :: ns_rho
      INTEGER :: nsegments_theta
      REAL(rprec) :: Rmin, Rmax
      REAL(rprec) :: Zmin, Zmax
      REAL(rprec) :: phi
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: testarray
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: B_car
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: my_fsurface

      CHARACTER(len=80) :: mychar1
      CHARACTER(len=80) :: mychar2, message
      CHARACTER(len=*), PARAMETER :: subname = 'READ_TEST_FOR_CDF: '

      INTEGER :: iou_cdf = 132
      INTEGER ::  arr_size
      INTEGER ::  size1, size2, size3

      INTEGER, DIMENSION(3) :: dim


      write(*,*) 'Hi.  Now in ', subname


!........open the input .nc file and read in the variables..............!
      call cdf_open(iou_cdf, cdf_infile, 'r', istat)

!      call cdf_inquire(iou_cdf, 'testarray', dim)
      call cdf_inquire(iou_cdf, 'B_car', dim)

!      write(*,*) subname, 'dim for testarray = ', dim
      write(*,*) subname, 'dim for B_car = ', dim
 
!!........allocate testarray size based on cdf_inquire results.......!
!      arr_size = dim(1)
!      write(*,*) subname, 'dim for testarray = ', arr_size
!      ALLOCATE( testarray(arr_size) )
!
!!........read in values for testarray...........................!
!      call cdf_read(iou_cdf, 'testarray', testarray)


!!........allocate B_car size based on cdf_inquire results.......!
      size1 = dim(1)
      size2 = dim(2)
      size3 = dim(3)
      write(*,*) subname, 'dim for B_car = ', dim
      ALLOCATE( B_car(0:(size1 -1),0:(size2 -1),size3) )

!........read in values for B_car...........................!
      call cdf_read(iou_cdf, 'B_car', B_car)

!.......repeat allocation process for flux_surface array...............!
      call cdf_inquire(iou_cdf, 'flux_surface', dim)
       size1 = dim(1)
       size2 = dim(2)
       size3 = dim(3)
      write(*,*) subname, 'dim for B_car = ', dim
      ALLOCATE( my_fsurface(size1, size2, size3) )


!........read in flux_surface values into local array "my_fsurface"...........................!
      call cdf_read(iou_cdf, 'flux_surface', my_fsurface)

!..........read in non-array variables.....................!
      CALL cdf_read(iou_cdf, 'Rmin', Rmin)
      CALL cdf_read(iou_cdf, 'Rmax', Rmax)
      CALL cdf_read(iou_cdf, 'Zmin', Zmin)
      CALL cdf_read(iou_cdf, 'Zmax', Zmax)
      CALL cdf_read(iou_cdf, 'phi', phi)
      CALL cdf_read(iou_cdf, 'nsegments_R', nsegments_R)
      CALL cdf_read(iou_cdf, 'nsegments_Z', nsegments_Z)
      CALL cdf_read(iou_cdf, 'ns_rho', ns_rho)
      CALL cdf_read(iou_cdf, 'nsegments_theta', nsegments_theta)
      call cdf_close(iou_cdf)

      write(*,*) subname, 'Rmin = ', Rmin
      write(*,*) subname, 'Rmax = ', Rmax
      write(*,*) subname, 'Zmin = ', Zmin
      write(*,*) subname, 'Zmax = ', Zmax
      write(*,*) subname, 'm = ', nsegments_R
      write(*,*) subname, 'n = ', nsegments_Z
      write(*,*) subname, 'phi = ', phi
      write(*,*) subname, 'ns_rho = ', ns_rho
      write(*,*) subname, 'nsegments_theta = ', nsegments_theta

!!.........test read-in of B_car array............!
!      do i = 0,nsegments_R
!        do j = 0, nsegments_Z
!          write(*,*) subname, 'B_car = ', B_car(i,j,:)
!        end do
!      end do 


!.........test read-in of flux_surface array............!
!      do i = 1, ns_rho
      do i = 30, ns_rho
        write(*,*) subname, 'flux suface number = ', i
        do j = 1, nsegments_theta
          write(*,*) subname, my_fsurface(i,j,:) 
         end do
      end do


      RETURN
      END SUBROUTINE READ_TEST_FOR_CDF


      end module bfield_grid_mod

