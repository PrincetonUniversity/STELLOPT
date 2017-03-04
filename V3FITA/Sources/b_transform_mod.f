!     SPH: INTEGER(iprec) -> INTEGER
      MODULE b_transform_mod
!--------------------------------------------------------------------------
! 
! FUNCTION: Module containing routines that compute coordinate transforms
!           for the  B field components in the bgrid_mod module. 
!
!------------------------------------------------------------------------------

      USE stel_kinds
      USE stel_constants
!      USE safe_open_mod
      USE v3_utilities     ! contains subroutine assert_eq
      USE v3read_wout      ! reads equilibrium from a wout file into AJAX
      USE math_utilities   ! module containing vector manipulation routines
      IMPLICIT NONE


      CONTAINS

!===================================================================================
      SUBROUTINE COMPUTE_B_CAR(x_cyl, B_car, iflag, MOD_B, B_CYLIN)
!--------------------------------------------------------------------------------
!
! Function:  Given an input cylinderical 3-vector x_cyl(r,phi,Z), this subroutine outputs
!            the local CARTESIAN B field 3-vector at that point.  
!
! IMPORTANT NOTE: in order for this routine to function correctly, the point at x_cyl
!                 *MUST* be inside the plasma (rho <= 1.0).  If x_cyl does not fufill
!                 this criteria, the routine should figure it out and exit gracefully.
!
! created 2/13/06 by J. Shields
!
! Modified 11/8/06 by JS: use VMEC routines for flux coord conversions instead of AJAX          
!----------------------------------------------------------------------------------
!      USE ajax_mod
      USE b_transform_vmec, ONLY: VMEC_B, INVALID_RHO_VMEC
      IMPLICIT NONE

!.........dummy variables..................................!
      INTEGER, OPTIONAL, INTENT(OUT)  :: iflag   ! error flag ( 0 = happiness)
      REAL(rprec), DIMENSION(3), INTENT(IN)  :: x_cyl   ! cylin position vector (R, phi, Z)
      REAL(rprec), DIMENSION(3), INTENT(OUT) :: B_car   ! B vector in Cartesian Coords (Bx, By, Bz)
      REAL(rprec), OPTIONAL, INTENT(OUT)     :: MOD_B   ! Magnitude of the B field
      REAL(rprec), DIMENSION(3), OPTIONAL, INTENT(OUT):: B_CYLIN ! B field in Cylin coords

!..........local variables.................................!
      INTEGER           :: i, j
      INTEGER           :: istat, ierr
      REAL(rprec)       :: rho
      REAL(rprec), DIMENSION(3) :: xflux, xflux2
      REAL(rprec), DIMENSION(3) :: x_flx
      REAL(rprec), DIMENSION(3) :: x_cyl_temp
      REAL(rprec), DIMENSION(6) :: g_cyl
      REAL(rprec), DIMENSION(3) :: B_cyl
      REAL(rprec), DIMENSION(3) :: B_cyl_temp
!      REAL(rprec), DIMENSION(3) :: B_cyl_ajax
!      REAL(rprec), DIMENSION(3) :: B_car_ajax
      REAL(rprec), DIMENSION(3) :: B_co
      REAL(rprec), DIMENSION(3) :: B_con
      REAL(rprec), DIMENSION(3) :: B_con_temp
!      REAL(rprec), DIMENSION(3) :: B_cyl_rewt    ! B_cyl reweighted using ratio of AJAX & VMEC B_mods 
!      REAL(rprec), DIMENSION(3) :: B_car_rewt    ! B_car reweighted using ratio of AJAX & VMEC B_mods 
!      REAL(rprec), DIMENSION(3) :: B_co_ratio    ! ratio of b_co_rewt/b_co_wout
!      REAL(rprec)               :: B_mod         ! B_mod from VMEC bmnc array
!      REAL(rprec)               :: B_mod_ajax     ! B_mod from AJAX_B
!      REAL(rprec)               :: B_mod_ratio    ! VMEC B_mod / AJAX B_mod
!      REAL(rprec)               :: mag_B_cyl_wout  ! B_mod from VMEC B_CO ---> B_Cyl comps
      REAL(rprec)               :: mag_B_car     ! magnitude of AJAX B_Car vector
      REAL(rprec)               :: mag_B_cyl     ! magnitude of AJAX B_cyl vector

      CHARACTER(len=60)  :: message
      CHARACTER(len=*), PARAMETER :: subname = 'COMPUTE_B_CAR: '

!      write(*,*) '===================================================='    
!      write(*,*) ' '
!      write(*,*) 'Hi.  Now in ', subname

!......initialize the error flag............................!
      istat = 0

!.........check if point x_cyl is inside the plasma........................!
      if ( .NOT. INVALID_RHO_VMEC(x_cyl) ) then


!.........use VMEC routine to extract B info from the wout.........!
!        call VMEC_B( x_cyl, B_cyl, ierr, R_FLX=x_flx)
        call VMEC_B( x_cyl, B_cyl, ierr, R_FLX=x_flx, B_CON=B_con)
        istat = ierr
      else
        write(*,*) subname, 'ERROR!!! Input coordinates OUTSIDE PLASMA!'
        istat = 1
        go to 100
      end if  ! end if for "inside plasma" check


!!.......compute the invariant dot product of the AJAX flux coords
!      B_dot_B = B_CO(1) * B_CON(1) + B_CO(2) * B_CON(2)                        &
!     &          + B_CO(3)* B_CON(3)

!.........convert cylin B field vector to Cartesian coords.................!
!      call B_CYL_TO_CAR(x_cyl, B_cyl, B_car, MOD_B)
      call B_CYL_TO_CAR(x_cyl, B_cyl, B_car)

!....compute magnitude of B_car & B_cyl vectors...........................!
      mag_B_car = SQRT(B_car(1)**2.0 + B_car(2)**2.0 + B_car(3)**2.0 ) 
      mag_B_cyl = SQRT(B_cyl(1)**2.0 + B_cyl(2)**2.0 + B_cyl(3)**2.0 ) 


      if( PRESENT(MOD_B)) MOD_B = mag_B_cyl
      if( PRESENT(B_CYLIN)) B_CYLIN = B_cyl

!      write(*,*) subname, 'mag_B_car      = ', mag_B_car
!      write(*,*) subname, 'mag_B_cyl      = ', mag_B_cyl
!      write(*,*) subname, 'FINAL B_cyl = ', B_cyl
!      write(*,*) subname, 'FINAL B_car = ', B_car
!      write(*,*) subname, 'FINAL B_mod = ', B_mod
!      write(*,*) '----------------------------------------'

!      write(*,*) 'Hi.  Now EXITING ', subname

 100  continue
      if( PRESENT(iflag) ) iflag = istat

      RETURN
      END SUBROUTINE COMPUTE_B_CAR




!===================================================================================
      SUBROUTINE COMPUTE_B_FROM_WOUT(xflux, mod_B, B_CO, B_CON)
!--------------------------------------------------------------------------------
!
! FUNCTION:  Determines the components of the EXTERNAL magnetic field (in flux coords)
!            at a point specified by the flux coordinate position vector "xflux" via a DIRECT
!            reading of an input wout file.
!
!        **Based on B_CYL_FROM_WOUT and the (now defunct) FLUX_COORD_B_EXTERNAL subroutine (which took 
!        the B arrays from memory instead of directly from the wout file)
!
! LOGIC: Takes values from bsubumnc, bsubumns, etc Fourier arrays that are output from VMEC
!         and then linearly interpolates between them.
!        Arrays are based on a Fourier expansion of form e.g. :
!          B = SUM_mn{ B_mn * cos( mu - nv) } 
!            = SUM_mn{ B_mn * cos(z) }        where z = mu - nv 
! 
!  **note that the order of the flux vector elements follows the standard convention:
!          (rho, theta, phi), or equivalently:  (s,u,v)
!
!  ** also note that there is NO array for the rho component of B, bsups(s,u,v) ! 
!   This is because  B_rho = 0  is part of the *DEFINITION* of the flux coords.
!             (see the paragraph after (3c) in the Hirshman 1983 paper)
!                  
! CAVEAT: there is currently NO safeguard against the case in which the path point
!         is outside the field vessel and hence the database contains no info about
!         the external field at that point.  This type of error will generally
!          cause the module to crash.  JS 7/1/05
!
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
!      USE track_mod
      USE read_wout_mod
      IMPLICIT NONE

!........passed variables....................................................!
      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            xflux            ! flux coordinate POSITION vector
      REAL(rprec), INTENT(OUT) ::                                              &
     &            mod_B            ! magnitude of MAGNETIC FIELD vector at point xflux
      REAL(rprec), DIMENSION(3), OPTIONAL, INTENT(OUT) ::                                &
     &            B_CO             ! COVAR flux coord MAGNETIC FIELD vector at point xflux
      REAL(rprec), DIMENSION(3), OPTIONAL, INTENT(OUT) ::                                &
     &            B_CON            ! CONTRAVAR flux coord MAGNETIC FIELD vector at point xflux


!..................local variables.............................................!
      INTEGER :: si_hi, si_low       ! indices of grid pts above, below "s"
      INTEGER :: mn_size
      INTEGER :: ierr

      REAL(rprec), ALLOCATABLE, DIMENSION(:) ::                                &
     &                 z, cosz, sinz

      REAL(rprec), ALLOCATABLE, DIMENSION(:) ::                                &
     &                 b_subs_hi, b_subs_low, b_subs,                          &
     &                 b_subu_hi, b_subu_low, b_subu,                          &
     &                 b_subv_hi, b_subv_low, b_subv                           &

      REAL(rprec), ALLOCATABLE, DIMENSION(:) ::                                &
     &                 b_supu_hi, b_supu_low, b_supu,                          &
     &                 b_supv_hi, b_supv_low, b_supv                           &

      REAL(rprec), ALLOCATABLE, DIMENSION(:) ::                                &
     &                 b_mag_hi, b_mag_low, b_mag

      REAL(rprec) :: s, u, v
      REAL(rprec) :: sgrid_low, sgrid_hi     ! radial grid points to left/right of "s"
      REAL(rprec) :: delta_sgrid             ! distance btwn radial grid points
      REAL(rprec) :: b_subs_sum, b_subu_sum, b_subv_sum
      REAL(rprec) :: b_supu_sum, b_supv_sum
      REAL(rprec) :: b_mag_sum
      REAL(rprec) :: b_subs_low_sum, b_subs_hi_sum
      REAL(rprec) :: b_subu_low_sum, b_subu_hi_sum
      REAL(rprec) :: b_subv_low_sum, b_subv_hi_sum

      REAL(rprec) :: b_dot_b, sqrt_b_dot_b
      REAL(rprec) :: metric_mod_B, c2c_mod_B
      REAL(rprec), DIMENSION(3) :: B_co_low, B_co_hi
      REAL(rprec), DIMENSION(3) :: B_co_temp

      LOGICAL           :: lasym_local
      CHARACTER(len=60) :: in_woutfile
      CHARACTER(len=60) :: message
      CHARACTER(len=*), PARAMETER  ::                                          &
     &                     subname = 'COMPUTE_B_FROM_WOUT: '

!............begin executable .................................................!

!      lasym_local =  model_a%eqstate%fixp%lasym
      lasym_local =  lasym
     

      s = xflux(1)
      u = xflux(2)
      v = xflux(3)

!.......allocate z and b_subu arrays based on size of bsubumnc mn modes.......!
!      mn_size = SIZE(model_a%eqstate%aux1%bsubumnc , 1)
      mn_size = SIZE(bsubumnc , 1)

      ALLOCATE( z(mn_size), cosz(mn_size), sinz(mn_size) )
      ALLOCATE( b_mag(mn_size) ) 
      ALLOCATE( b_mag_low(mn_size), b_mag_hi(mn_size) ) 

      if( PRESENT(B_CO) ) then
        ALLOCATE( b_subs_low(mn_size), b_subs_hi(mn_size) ) 
        ALLOCATE( b_subu_low(mn_size), b_subu_hi(mn_size) ) 
        ALLOCATE( b_subv_low(mn_size), b_subv_hi(mn_size) )
        ALLOCATE( b_subs(mn_size), b_subu(mn_size), b_subv(mn_size) ) 
      end if

      if( PRESENT(B_CON) ) then
        ALLOCATE( b_supu_low(mn_size), b_supu_hi(mn_size) ) 
        ALLOCATE( b_supv_low(mn_size), b_supv_hi(mn_size) ) 
        ALLOCATE( b_supu(mn_size), b_supv(mn_size) ) 
      end if

!......calculate the argument of the sinusoidal functions...............!
!      z    = model_a%eqstate%aux1%xm_nyq * u -                                       &
!     &       model_a%eqstate%aux1%xn_nyq * v
      z    = xm_nyq * u -  xn_nyq * v
      cosz = COS(z)
      sinz = SIN(z)


!.............determine the indices of the radial grid points that most closely....!
!.............correspond to the value of the input radial variable "s"..............!
      call S_BOUNDS(s, sgrid_low, sgrid_hi, si_low, si_hi)
      delta_sgrid = sgrid_hi - sgrid_low
      if ( ABS(delta_sgrid) == 0.0 ) then
         write(*,*) subname, 'WARNING! delta_sgrid DIVIDE BY 0 ERROR'
      end if 
  

!.......calculate B field at the grid pts to left, right of input value "s"......!
      b_mag_low(:) = bmnc(:,si_low)
      b_mag_hi(:)  = bmnc(:,si_hi)

      if( PRESENT(B_CO) ) then
        b_subs_low(:) = bsubsmns(:,si_low)
        b_subs_hi(:)  = bsubsmns(:,si_hi)

!        b_subu_low(:) = model_a%eqstate%aux1%bsubumnc(:,si_low)
!        b_subu_hi(:)  = model_a%eqstate%aux1%bsubumnc(:,si_hi)
        b_subu_low(:) = bsubumnc(:,si_low)
        b_subu_hi(:)  = bsubumnc(:,si_hi)

!        b_subv_low(:) = model_a%eqstate%aux1%bsubvmnc(:,si_low)
!        b_subv_hi(:)  = model_a%eqstate%aux1%bsubvmnc(:,si_hi)
        b_subv_low(:) = bsubvmnc(:,si_low)
        b_subv_hi(:)  = bsubvmnc(:,si_hi)
      end if

      if( PRESENT(B_CON) ) then
        b_supu_low(:) = bsupumnc(:,si_low)
        b_supu_hi(:)  = bsupumnc(:,si_hi)

        b_supv_low(:) = bsupvmnc(:,si_low)
        b_supv_hi(:)  = bsupvmnc(:,si_hi)
      end if




!.......Find B field coefficient at "s" via linear interpolation btwn sgrid_low and sgrid_hi.....!
      b_mag = ( (sgrid_hi - s) * b_mag_low + (s - sgrid_low) *                 &
     &         b_mag_hi ) / delta_sgrid


      if( PRESENT(B_CO) ) then
        b_subs = ( (sgrid_hi - s) * b_subs_low + (s - sgrid_low) *                 &
     &         b_subs_hi ) / delta_sgrid

        b_subu = ( (sgrid_hi - s) * b_subu_low + (s - sgrid_low) *                 &
     &         b_subu_hi ) / delta_sgrid

        b_subv = ( (sgrid_hi - s) * b_subv_low + (s - sgrid_low) *                 &
     &         b_subv_hi ) / delta_sgrid
      end if


      if( PRESENT(B_CON) ) then
        b_supu = ( (sgrid_hi - s) * b_supu_low + (s - sgrid_low) *                 &
     &         b_supu_hi ) / delta_sgrid

        b_supv = ( (sgrid_hi - s) * b_supv_low + (s - sgrid_low) *                 &
     &         b_supv_hi ) / delta_sgrid
      end if




!......sum over modes to obtain the NET B_u and B_v components...........!
      b_mag_sum = SUM( b_mag * cosz)

      if (PRESENT(B_CO) ) then
        b_subs_sum = SUM( b_subs * sinz)
        b_subu_sum = SUM( b_subu * cosz)
        b_subv_sum = SUM( b_subv * cosz)
      end if

      if (PRESENT(B_CON) ) then
        b_supu_sum = SUM( b_supu * cosz)
        b_supv_sum = SUM( b_supv * cosz)
      end if

!      b_subs_low_sum = SUM( b_subs_low * sinz)
!      b_subs_hi_sum  = SUM( b_subs_hi * sinz)
!      b_subu_low_sum = SUM( b_subu_low * cosz)
!      b_subu_hi_sum  = SUM( b_subu_hi * cosz)
!      b_subv_low_sum = SUM( b_subv_low * cosz)
!      b_subv_hi_sum  = SUM( b_subv_hi * cosz)



!...........include an option for the asymmetric case..............!
!.........(Warning: mod B and contravar B not implemented for asymmetric case....!  
      if (lasym_local) then

        write(*,*) subname, 'WARNING!  lasym_local = TRUE'

        if (PRESENT(B_CO) ) then
!.......calculate B field at the grid pts to left, right of input value "s"......!
!          b_subu_low(:) = model_a%eqstate%aux1%bsubumns(:,si_low)
!          b_subu_hi(:)  = model_a%eqstate%aux1%bsubumns(:,si_hi)
!          b_subv_low(:) = model_a%eqstate%aux1%bsubvmns(:,si_low)
!          b_subv_hi(:)  = model_a%eqstate%aux1%bsubvmns(:,si_hi)
          b_subu_low(:) = bsubumns(:,si_low)
          b_subu_hi(:)  = bsubumns(:,si_hi)
          b_subv_low(:) = bsubvmns(:,si_low)
          b_subv_hi(:)  = bsubvmns(:,si_hi)


!.......Find B field at "s" via linear interpolation btwn sgrid_low and sgrid_hi.....!
          b_subu = ( (sgrid_hi - s) * b_subu_low + (s - sgrid_low) *               &
     &           b_subu_hi ) / delta_sgrid

          b_subv = ( (sgrid_hi - s) * b_subv_low + (s - sgrid_low) *               &
     &           b_subv_hi ) / delta_sgrid

!..........add sine terms to the previous (cosine) sum over modes..............!
          b_subu_sum = b_subu_sum + SUM( b_subu * sinz)
          b_subv_sum = b_subv_sum * SUM( b_subv * sinz)
        end if
      end if      ! end of "asymmetric case" option



!........assign final values to B field flux coords (keeping in mind that flux coords
!........are DEFINED with B_sub_s = 0 )...........................................!
      mod_B = b_mag_sum

      if (PRESENT(B_CO) ) then
        B_co(1)     = b_subs_sum
        B_co(2)     = b_subu_sum
        B_co(3)     = b_subv_sum
      end if

      if (PRESENT(B_CON) ) then
        B_con(1) = 0.0
        B_con(2) = b_supu_sum
        B_con(3) = b_supv_sum
      end if

 
!........deallocate arrays..........................................!   
      DEALLOCATE( z, cosz, sinz )
      DEALLOCATE( b_subu_low, b_subu_hi ) 
      DEALLOCATE( b_subv_low, b_subv_hi ) 
      DEALLOCATE( b_subu, b_subv ) 

!      write(*,*) 'Hi.  Now exiting ', subname

      RETURN
      END SUBROUTINE COMPUTE_B_FROM_WOUT




!==================================================================================
      SUBROUTINE B_COVAR_TO_CYL(x_cyl, g_cyl, B_co, B_cyl)
!----------------------------------------------------------------------------------
!  FUNCTION:  computes the magnetic field components in standard cylinderical coords, given an
!             input covariant B vector.
!
!  created by J. Shields 2/9/06
!
!  LOGIC:  the components transform via the tensor (t^ij) via
!
!          B^i = t^(i1)* B_1 + t^(i2)* B_2 + t^(i3)* B_3 
!
! where t^ij is the INVERSE of the cyl ---> (covar flux) matrix used by Houlberg in AJAX_B
!
!   Specifically, the "forward" cyl --> covar flux transform that Houlberg uses is:
!
!    B_CO(1) = b_cylt(1)*g_cyl(1) + b_cylt(3)*g_cyl(4)
!
!    B_CO(2) = b_cylt(1)*g_cyl(2) + b_cylt(3)*g_cyl(5)
!
!    B_CO(3) = b_cylt(1)*g_cyl(3) + b_cylt(2)*r_cyl(1) + b_cylt(3)*g_cyl(6)
!
!
!  **Note that by "standard cylinderical coordinates", I mean coordinates associated with
!         with the normalized unit vectors  R_hat, phi_hat, Z_hat.
!
!--------------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE

!........passed variables....................................................!
      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            x_cyl              ! cylin coordinate POSITION vector

      REAL(rprec), DIMENSION(6), INTENT(IN) ::                                 &
     &            g_cyl              ! cylin/flux derivatives
                                     !=(R_rho,R_theta,R_zeta,Z_rho,Z_theta,Z_zeta)
                                     ! [m/rho,m,      m,     m/rho,m,      m     ]

      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            B_co          ! COVAR flux coord MAGNETIC FIELD vector at point x_cyl

      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            B_cyl         ! CYLINDERICAL coord MAGNETIC FIELD vector at point x_cyl


!.................local variables...................................!
      INTEGER :: i
      REAL(rprec) :: ns
      REAL(rprec) :: tau                   ! 2D Jacobian
      REAL(rprec) :: rtau                  ! R * (2D Jacobian)
      REAL(rprec), DIMENSION(3,3) :: t     ! inverse transformation matrix

      character(len=*), PARAMETER ::  subname = 'B_COVAR_TO_CYL'

      tau  = g_cyl(2)*g_cyl(4) - g_cyl(1)*g_cyl(5)
      rtau = x_cyl(1) * tau

!.....express the inverse transform tensor in terms of the g_cyl derivatives.........................!
      t(1,1) = -1.0 * ( g_cyl(5)/tau )
      t(1,2) = g_cyl(4) / tau 
      t(1,3) = 0.0
      t(2,1) = ( g_cyl(3) * g_cyl(5) - g_cyl(2) * g_cyl(6) ) / rtau
      t(2,2) = ( g_cyl(1) * g_cyl(6) - g_cyl(3) * g_cyl(4) ) / rtau
      t(2,3) = 1.0 / x_cyl(1)
      t(3,1) = g_cyl(2) / tau 
      t(3,2) = -1.0 * ( g_cyl(1) / tau ) 
      t(3,3) = 0.0

!........compute the cylin components(B_R, B_phi,B_Z) using the input covariant vector B_co.....!
      do  i = 1,3
        B_cyl(i) = t(i,1)*B_co(1) + t(i,2)*B_co(2)  + t(i,3)*B_co(3)
      end do

!      write(*,*) subname, 'final B_CYL = ', B_cyl

      RETURN
      END SUBROUTINE B_COVAR_TO_CYL



!==================================================================================
      SUBROUTINE B_CONTRAVAR_TO_CYL(x_cyl, g_cyl, B_con, B_cyl)
!----------------------------------------------------------------------------------
!  FUNCTION:  computes the magnetic field components in standard cylinderical coords, given an
!             input CONTRAvariant B vector.
!
!  created by J. Shields 2/9/06, based on W. Houlberg's AJAX_B subroutine
!
!
!
!  **Note that by "standard cylinderical coordinates", I mean coordinates associated with
!         with the normalized unit vectors  R_hat, phi_hat, Z_hat.
!
!--------------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE

!........passed variables....................................................!
      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            x_cyl              ! cylin coordinate POSITION vector

      REAL(rprec), DIMENSION(6), INTENT(IN) ::                                 &
     &            g_cyl              ! cylin/flux derivatives
                                     !=(R_rho,R_theta,R_zeta,Z_rho,Z_theta,Z_zeta)
                                     ! [m/rho,m,      m,     m/rho,m,      m     ]

      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            B_con         ! CONTRAVAR flux coord MAGNETIC FIELD vector at point x_cyl

      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            B_cyl         ! CYLINDERICAL coord MAGNETIC FIELD vector at point x_cyl


!.................local variables...................................!
      INTEGER :: i
      REAL(rprec) :: ns

      character(len=*), PARAMETER ::  subname = 'B_CONTRAVAR_TO_CYL'


!.............B_R component........................................!
      B_cyl(1) = g_cyl(2)*B_con(2) + g_cyl(3)*B_con(3)

!.............B_phi component.....................................!
      B_cyl(2) = x_cyl(1)*B_con(3)

!.............B_Z component........................................!
      B_cyl(3) = g_cyl(5)*B_con(2) + g_cyl(6)*B_con(3)



!      write(*,*) subname, 'final B_CYL = ', B_cyl

      RETURN
      END SUBROUTINE B_CONTRAVAR_TO_CYL



!==================================================================================
      SUBROUTINE B_CYL2FLUX(x_cyl, g_cyl, B_cyl, B_con, B_co)
!----------------------------------------------------------------------------------
!  FUNCTION:  computes the magnetic field components in flux coordinates (both
!             covariant & contravariant), given an input B vector expressed in 
!             standard cylinderical coords,
!
!  created by J. Shields 2/9/06, based on W. Houlberg's AJAX_B subroutine
!
!
!
!  **Note that by "standard cylinderical coordinates", I mean coordinates associated with
!         with the normalized unit vectors  R_hat, phi_hat, Z_hat.
!
!--------------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE

!........passed variables....................................................!
      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            x_cyl              ! cylin coordinate POSITION vector

      REAL(rprec), DIMENSION(6), INTENT(IN) ::                                 &
     &            g_cyl              ! cylin/flux derivatives
                                     !=(R_rho,R_theta,R_zeta,Z_rho,Z_theta,Z_zeta)
                                     ! [m/rho,m,      m,     m/rho,m,      m     ]

      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            B_cyl         ! CYLINDERICAL coord MAGNETIC FIELD vector at point x_cyl

      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            B_con         ! CONTRAVAR flux coord MAGNETIC FIELD vector at point x_cyl

      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            B_co          ! COVAR flux coord MAGNETIC FIELD vector at point x_cyl

      
!.................local variables...................................!
      INTEGER :: i
      REAL(rprec) :: ns

      REAL(rprec), DIMENSION(2,2) :: t    ! 2D  (B_Cyl ---> B_Co) transform matrix

      character(len=*), PARAMETER ::  subname = 'B_CYL2FLUX'



!......compute the COVARIANT flux vector B_CO from B_cyl.....................!
!.......(using the same transform matrix as in AJAX_B).............!

!........B_rho component...............................................!
      B_co(1) = B_cyl(1)*g_cyl(1) + B_cyl(3)*g_cyl(4)

!........B_theta component..............................................!
      B_co(2) = B_cyl(1)*g_cyl(2) + B_cyl(3)*g_cyl(5)

!........B_zeta component...............................................!
      B_co(3) = B_cyl(1)*g_cyl(3)+B_cyl(2)*x_cyl(1)+B_cyl(3)*g_cyl(6)




!......compute the CONTRAVARIANT flux vector B_CON from B_cyl.....................!

!........write out the 2D {B^phi,B^Z} ---> (B^theta,B^zeta} transform matrix......!
      t(1,1) = ( -1.0/x_cyl(1) ) * ( g_cyl(6) / g_cyl(5) )
      t(1,2) = 1.0 / g_cyl(5)
      t(2,1) = 1.0 / x_cyl(1)
      t(2,2) = 0.0

!........B^rho always = 0.0 in flux coords, *by construction*................!
      B_con(1) = 0.0

!.....B^theta..............................................................!
      B_con(2) = t(1,1) * B_cyl(2) + t(1,2) * B_cyl(3)

!.......B_zeta.............................................................!
      B_con(3) = t(2,1) * B_cyl(2)


      RETURN
      END SUBROUTINE B_CYL2FLUX



!==================================================================================
      SUBROUTINE B_COVAR_TO_CONTRAVAR(x_cyl, g_cyl, B_co, B_con, MOD_B)
!----------------------------------------------------------------------------------
!  FUNCTION:  computes the contravariant flux coordinate magnetic field, given
!             an input COVARIANT B field vector.
!
!  created by J. Shields 2/10/06
!
!
!
!
!--------------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE

!........passed variables....................................................!
      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            x_cyl             ! cylin coordinate POSITION vector

      REAL(rprec), DIMENSION(6), INTENT(IN) ::                                 &
     &            g_cyl             ! cylin/flux derivatives
                                    !=(R_rho,R_theta,R_zeta,Z_rho,Z_theta,Z_zeta)
                                    ! [m/rho,m,      m,     m/rho,m,      m     ]

      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                &
     &            B_co              ! COVAR flux coord MAGNETIC FIELD vector at point x_cyl

      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            B_con             ! CONTRAVAR flux coord MAGNETIC FIELD vector at point x_cyl

      REAL(rprec), OPTIONAL, INTENT(OUT) :: MOD_B    ! magnitude of B field vector

      
!.................local variables...................................!
      INTEGER :: i
      REAL(rprec), DIMENSION(3) ::                                              &
     &            B_cyl             ! CYLINDERICAL coord magnetic field vector at point x_cyl

      REAL(rprec), DIMENSION(3) ::  B_cov, B_cont
      character(len=*), PARAMETER ::  subname = 'B_COVAR_TO_CONTRAVAR'


!.......convert the covariant B vector into standard cylin coords..............!
      call B_COVAR_TO_CYL(x_cyl, g_cyl, B_co, B_cyl)      

!........convert the cylin B vector into contravariant components................!
      call B_CYL2FLUX(x_cyl, g_cyl, B_cyl, B_cont, B_cov)

      B_con = B_cont

      if ( PRESENT(MOD_B) ) then
!        MOD_B = SQRT( B_cyl(1)**2.0 + B_cyl(2)**2.0 + B_cyl(3)**2.0 )
        MOD_B = SQRT( B_con(2) * B_co(2) +  B_con(3) * B_co(3)   )
      end if

      RETURN
      END SUBROUTINE B_COVAR_TO_CONTRAVAR



!==================================================================================
      SUBROUTINE B_CONTRAVAR_TO_COVAR(x_cyl, g_cyl, B_con, B_co, MOD_B)
!----------------------------------------------------------------------------------
!  FUNCTION:  computes the contravariant flux coordinate magnetic field, given
!             an input COVARIANT B field vector.
!
!  created by J. Shields 2/10/06
!
!
!
!
!--------------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE

!........passed variables....................................................!
      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            x_cyl            ! cylin coordinate POSITION vector

      REAL(rprec), DIMENSION(6), INTENT(IN) ::                                 &
     &            g_cyl            ! cylin/flux derivatives
                                   !=(R_rho,R_theta,R_zeta,Z_rho,Z_theta,Z_zeta)
                                   ! [m/rho,m,      m,     m/rho,m,      m     ]

      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                &
     &            B_con            ! CONTRAVAR flux coord MAGNETIC FIELD vector at point x_cyl

      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            B_co             ! COVAR flux coord MAGNETIC FIELD vector at point x_cyl

      REAL(rprec), OPTIONAL, INTENT(OUT) :: MOD_B
      
!.................local variables...................................!
      INTEGER :: i
      REAL(rprec), DIMENSION(3) ::                                              &
     &            B_cyl         ! CYLINDERICAL coord magnetic field vector at point x_cyl

      REAL(rprec), DIMENSION(3) ::  B_cov, B_cont
      character(len=*), PARAMETER ::  subname = 'B_CONTRAVAR_TO_COVAR'


!.......convert the covariant B vector into standard cylin coords..............!
      call B_CONTRAVAR_TO_CYL(x_cyl, g_cyl, B_con, B_cyl)      

!........convert the cylin B vector into contravariant components................!
      call B_CYL2FLUX(x_cyl, g_cyl, B_cyl, B_cont, B_cov)

      B_co = B_cov

      if ( PRESENT(MOD_B) ) then
!        MOD_B = SQRT( B_cyl(1)**2.0 + B_cyl(2)**2.0 + B_cyl(3)**2.0 )
        MOD_B = SQRT( B_con(2) * B_co(2) +  B_con(3) * B_co(3)   )
      end if

      RETURN
      END SUBROUTINE B_CONTRAVAR_TO_COVAR


!==================================================================================
      SUBROUTINE B_CYL_TO_CAR(x_cyl, B_cyl, B_car, MOD_B)
!----------------------------------------------------------------------------------
!  FUNCTION:  computes the CARTESIAN coordinate magnetic field vector, given
!             an input CYLINDERICAL B field vector.
!
!  created by J. Shields 2/10/06
!
!  LOGIC:  The 2D (R,phi) ---> (x,y) transform matrix is
!
!   [B_x]  =  [ cos_phi     -sin_phi ]  [ B_R ]
!   [B_y]  =  [ sin_phi      cos_phi ]  [B_phi]
!
!--------------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE

!........passed variables....................................................!
      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            x_cyl              ! cylin coordinate POSITION vector

      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                &
     &            B_cyl         !cylin coord MAGNETIC FIELD vector at point x_cyl

      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            B_car          ! Cartesian coord MAGNETIC FIELD vector at point x_cyl

      REAL(rprec), OPTIONAL, INTENT(OUT) :: MOD_B

      
!.................local variables...................................!
      INTEGER :: i
      character(len=*), PARAMETER ::  subname = 'B_CYL_TO_CAR'


!..................B_x component.....................................!
      B_CAR(1) = B_cyl(1)*COS(x_cyl(2)) - B_cyl(2)*SIN(x_cyl(2))

!...................B_y component.....................................!
      B_CAR(2) = B_cyl(1)*SIN(x_cyl(2)) + B_cyl(2)*COS(x_cyl(2))

!...................B_z component.....................................!
      B_CAR(3) = B_cyl(3)


      if ( PRESENT(MOD_B) ) then
        MOD_B = SQRT( B_cyl(1)**2.0 + B_cyl(2)**2.0 + B_cyl(3)**2.0 )
      end if

      RETURN
      END SUBROUTINE B_CYL_TO_CAR

!==================================================================================
      SUBROUTINE B_CAR_TO_CYL(x_cyl, B_car, B_cyl, MOD_B)
!----------------------------------------------------------------------------------
!  FUNCTION:  computes the CYLINDERICAL coordinate magnetic field vector, given
!             an input CARTESIAN B field vector.
!
!  created by J. Shields 3/7/07
!
!  LOGIC:  The 2D  (x,y) --->(R,phi) transform matrix is
!
!   [ B_R ]  =  [  cos_phi     sin_phi ]  [ B_x ]
!   [B_phi]  =  [ -sin_phi     cos_phi ]  [ B_y ]
!
!--------------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE

!............passed variables....................................................!
      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            x_cyl            ! cylin coordinate POSITION vector

      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                &
     &            B_car            ! Cartesian coord MAGNETIC FIELD vector at point x_cyl

      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            B_cyl            ! cylin coord MAGNETIC FIELD vector at point x_cyl

      REAL(rprec), OPTIONAL, INTENT(OUT) :: MOD_B

!..............local variables...................................!
      INTEGER :: i
      character(len=*), PARAMETER ::  subname = 'B_CAR_TO_CYL: '


!..................B_R component.....................................!
      B_cyl(1) = B_car(1)*COS(x_cyl(2)) + B_car(2)*SIN(x_cyl(2))

!...................B_phi component.....................................!
      B_cyl(2) = -1.0*B_car(1)*SIN(x_cyl(2)) + B_car(2)*COS(x_cyl(2))

!...................B_z component.....................................!
      B_cyl(3) = B_car(3)


      if ( PRESENT(MOD_B) ) then
        MOD_B = SQRT( B_car(1)**2.0 + B_car(2)**2.0 + B_car(3)**2.0 )
      end if

      RETURN
      END SUBROUTINE B_CAR_TO_CYL



!===================================================================================
      SUBROUTINE S_BOUNDS(s, sgrid_low, sgrid_hi, si_low, si_hi)
!--------------------------------------------------------------------------------
!
! FUNCTION:  Determines the indices and radial flux coordinates ("rho") for the 
!            (HALF grid) radial grid points that bound point "s" from the left and right. 
!
! LOGIC: The radial variable "s" is defined on a HALF grid:  s = (i - 1.5)/(ns - 1)
!         Therefore  i = s * (ns - 1) + 1.5
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      USE track_mod
      USE eq_T         ! module containing VMEC output info
!      USE v3f_global   ! module containing VMEC output info (i.e. it SAVEs model_a)
      USE read_wout_mod, only : read_wout_file, rmnc, zmns, gmnc,              &
     &   currumnc, currvmnc, bsubumnc, bsubvmnc, lasym, rmns, zmnc,            &
     &   lmns, gmns, currumns, currvmns, bsubumns, bsubvmns, xm, xn,           &
     &   xm_nyq, xn_nyq
      IMPLICIT NONE


      REAL(rprec), INTENT(IN) :: s
      INTEGER, INTENT(OUT) ::                                                  &
     &               si_hi, si_low     ! indices of grid pts above, below "s"
      REAL(rprec), INTENT(OUT)::                                               &
     &               sgrid_hi, sgrid_low  ! s values corresponding to si_hi, si_low

!.................local variables...................................!
      REAL(rprec) :: ns
      character(len=*), PARAMETER ::  subname = 'S_BOUNDS: '
!      write(*,*) 'Hi.  Now in ', subname

!......determine the number of radial grid points (ns)...........!
!      ns = SIZE(  model_a%eqstate%aux1%bsubumnc , 2)
      ns = SIZE(bsubumnc, 2)

!.......calculate indices for the grid points below, above point "s"..........!
      si_low = AINT (s * (ns - 1) + 1.5)
      si_hi  = si_low + 1.0

!.........calculate radial flux coord for points at si_low & si_hi.........!
      sgrid_low = ( REAL(si_low) - 1.5 )  / (ns - 1.0) 
      sgrid_hi  = ( REAL(si_hi) - 1.5 )  / (ns - 1.0) 

!      write(*,*) 'Hi.  Now exiting ', subname

      RETURN
      END SUBROUTINE S_BOUNDS


      END MODULE b_transform_mod

