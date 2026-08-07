!-----------------------------------------------------------------------
!     Module:        surface_extender
!     Authors:       D. Pfefferlé (david.pfefferle@gauss-fusion.com)
!     Date:          11/07/2026
!     Description:   This module contains routines for extending the
!                    coordinate system beyond the last-closed flux-
!                    surface using splined Fourier coefficients from
!                    the surfaceflow python utility.
!-----------------------------------------------------------------------
MODULE surface_extender_mod

  !-----------------------------------------------------------------------
  !     Libraries
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !     Module Variables
  !-----------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE

  INTEGER :: nmodes
  INTEGER :: fit_type
  INTEGER :: degrees_of_freedom
  INTEGER :: nfp

  DOUBLE PRECISION :: rhoscale

  INTEGER, PARAMETER :: FIT_POLYNOMIAL    = 1
  INTEGER, PARAMETER :: FIT_ENDPOINTCUBIC = 2

  INTEGER, ALLOCATABLE :: m(:), n(:)

  ! Polynomial fit
  DOUBLE PRECISION, ALLOCATABLE :: Rc_coeff(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Zc_coeff(:,:)

  ! Endpoint cubic fit
  DOUBLE PRECISION, ALLOCATABLE :: R0(:), R1(:)
  DOUBLE PRECISION, ALLOCATABLE :: Z0(:), Z1(:)
  DOUBLE PRECISION, ALLOCATABLE :: Rc_ab(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Zc_ab(:,:)

  ! Workspace: Fourier coefficients at a fixed radius
  DOUBLE PRECISION, ALLOCATABLE :: Rc(:), Zc(:)
  DOUBLE PRECISION, ALLOCATABLE :: dRc_dt(:),dZc_dt(:)

  PUBLIC :: load_surface_fit
  PUBLIC :: unload_surface_fit
  PUBLIC :: rhothetazeta2xyz
  PUBLIC :: xyz2rhothetazeta
  PUBLIC :: get_rhoscale
  PUBLIC :: get_field_period

CONTAINS

  INTEGER FUNCTION get_field_period()
    IMPLICIT NONE
    get_field_period=nfp
  END FUNCTION get_field_period

  DOUBLE PRECISION FUNCTION get_rhoscale()
    IMPLICIT NONE
    get_rhoscale=rhoscale
  END FUNCTION get_rhoscale


  SUBROUTINE load_surface_fit(filename)
    USE netcdf

    IMPLICIT NONE

    CHARACTER(*), INTENT(in) :: filename

    INTEGER :: ncid
    INTEGER :: varid
    INTEGER :: dimid
    INTEGER :: ierr

    CHARACTER(len=32) :: fit_name

    !--------------------------------------
    ! Deallocate arrays if load_surface_fit called again
    !--------------------------------------

    CALL unload_surface_fit()

    !--------------------------------------
    ! Open file
    !--------------------------------------

    ierr = nf90_open(filename, NF90_NOWRITE, ncid)
    IF (ierr /= nf90_noerr) STOP nf90_strerror(ierr)

    !--------------------------------------
    ! Dimensions
    !--------------------------------------

    ierr = nf90_inq_dimid(ncid, "mode", dimid)
    IF (ierr /= nf90_noerr) STOP "Error dimension mode not found."
    ierr = nf90_inquire_dimension(ncid, dimid, len=nmodes)

    ierr = nf90_inq_dimid(ncid, "dof", dimid)
    IF (ierr /= nf90_noerr) STOP "Error dimension dof not found."
    ierr = nf90_inquire_dimension(ncid, dimid, len=degrees_of_freedom)

    ALLOCATE(m(nmodes))
    ALLOCATE(n(nmodes))

    ALLOCATE(Rc(nmodes))
    ALLOCATE(Zc(nmodes))

    ALLOCATE(dRc_dt(nmodes))
    ALLOCATE(dZc_dt(nmodes))

    !--------------------------------------
    ! Attributes
    !--------------------------------------

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "fit_type", fit_name)
    IF (ierr /= nf90_noerr) STOP "Error fit_type not found."
    ierr = nf90_get_att(ncid, NF90_GLOBAL, "nfp", nfp)
    IF (ierr /= nf90_noerr) STOP "Error nfp not found."
    ierr = nf90_get_att(ncid, NF90_GLOBAL, "rhoscale", rhoscale)
    IF (ierr /= nf90_noerr) STOP "Error rhoscale not found."

    !--------------------------------------
    ! Basis
    !--------------------------------------

    ierr = nf90_inq_varid(ncid, "xm", varid)
    IF (ierr /= nf90_noerr) STOP "Error xm not found."
    ierr = nf90_get_var(ncid, varid, m)
    IF (ierr /= nf90_noerr) STOP "Error reading m."

    ierr = nf90_inq_varid(ncid, "xn", varid)
    IF (ierr /= nf90_noerr) STOP "Error xn not found."
    ierr = nf90_get_var(ncid, varid, n)
    IF (ierr /= nf90_noerr) STOP "Error reading n."

    !--------------------------------------
    ! Polynomial fit
    !--------------------------------------

    IF (TRIM(fit_name) == "polynomial") THEN

       fit_type = FIT_POLYNOMIAL

       ALLOCATE(Rc_coeff(nmodes,degrees_of_freedom))    ! opposite convention to python
       ALLOCATE(Zc_coeff(nmodes,degrees_of_freedom))

       ierr = nf90_inq_varid(ncid,"rmnc",varid)
       ierr = nf90_get_var(ncid,varid,Rc_coeff)

       ierr = nf90_inq_varid(ncid,"zmns",varid)
       ierr = nf90_get_var(ncid,varid,Zc_coeff)

    ELSE IF (TRIM(fit_name) == "endpoint_cubic") THEN

       fit_type = FIT_ENDPOINTCUBIC

       ALLOCATE(R0(nmodes),R1(nmodes),Z0(nmodes),Z1(nmodes))
       ALLOCATE(Rc_ab(nmodes,degrees_of_freedom),Zc_ab(nmodes,degrees_of_freedom))

       ierr = nf90_inq_varid(ncid, "R0", varid)
       ierr = nf90_get_var(ncid, varid, R0)
       IF (ierr /= nf90_noerr) STOP "Error reading R0"

       ierr = nf90_inq_varid(ncid, "R1", varid)
       ierr = nf90_get_var(ncid, varid, R1)
       IF (ierr /= nf90_noerr) STOP "Error reading R1"

       ierr = nf90_inq_varid(ncid, "Z0", varid)
       ierr = nf90_get_var(ncid, varid, Z0)
       IF (ierr /= nf90_noerr) STOP "Error reading Z0"

       ierr = nf90_inq_varid(ncid, "Z1", varid)
       ierr = nf90_get_var(ncid, varid, Z1)
       IF (ierr /= nf90_noerr) STOP "Error reading Z1"

       ierr = nf90_inq_varid(ncid, "Rc_ab", varid)
       ierr = nf90_get_var(ncid, varid, Rc_ab)
       IF (ierr /= nf90_noerr) STOP "Error reading Rc_ab"

       ierr = nf90_inq_varid(ncid, "Zc_ab", varid)
       ierr = nf90_get_var(ncid, varid, Zc_ab)
       IF (ierr /= nf90_noerr) STOP "Error reading Zc_ab"

    ELSE

       STOP "Unknown fit type."

    END IF

    ierr = nf90_close(ncid)

  END SUBROUTINE load_surface_fit

  SUBROUTINE evaluate_coefficients(t)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in) :: t

    INTEGER :: i, k

    SELECT CASE (fit_type)

       !----------------------------------------
       ! Polynomial fit
       !----------------------------------------

    CASE (FIT_POLYNOMIAL)

       DO i = 1, nmodes

          Rc(i) = Rc_coeff(i,1)
          Zc(i) = Zc_coeff(i,1)
          dRc_dt(i) = 0.0
          dZc_dt(i) = 0.0

          DO k = 2, degrees_of_freedom
             dRc_dt(i) = dRc_dt(i)*t + Rc(i)
             dZc_dt(i) = dZc_dt(i)*t + Zc(i)
             Rc(i) = Rc(i)*t + Rc_coeff(i,k)
             Zc(i) = Zc(i)*t + Zc_coeff(i,k)
          END DO

       END DO

       !----------------------------------------
       ! Endpoint cubic
       !----------------------------------------

    CASE (FIT_ENDPOINTCUBIC)

       Rc = (1.0 - t) * R0 + t * R1 &
            + t * (1.0 - t) * (Rc_ab(:,1) * t + Rc_ab(:,2))
       Zc = (1.0 - t) * Z0 + t * Z1 &
            + t * (1.0 - t) * (Zc_ab(:,1) * t + Zc_ab(:,2))
       dRc_dt = R1 - R0 &
            + (1.0 -2.0*t)*(Rc_ab(:,1)*t + Rc_ab(:,2)) &
            + t*(1.0-t)*Rc_ab(:,1)

       dZc_dt = Z1 - Z0 &
            + (1.0 -2.0*t)*(Zc_ab(:,1)*t + Zc_ab(:,2)) &
            + t*(1.0-t)*Zc_ab(:,1)

    CASE default

       STOP "Unknown fit_type"

    END SELECT

  END SUBROUTINE evaluate_coefficients


  SUBROUTINE evaluate_RZ(theta,zeta,R,Z)
    ! This routine must be called after an update of the coefficients (evaluate_coefficients(t)) if the radial coordinate changes
    ! zeta is the NESCOIL angle i.e. zeta=PHI*nfp where PHI is the geometric toroidal angle
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in)  :: theta
    DOUBLE PRECISION, INTENT(in)  :: zeta

    DOUBLE PRECISION, INTENT(out) :: R
    DOUBLE PRECISION, INTENT(out) :: Z

    INTEGER  :: i
    DOUBLE PRECISION :: phase

    R = 0.0
    Z = 0.0

    DO i = 1, nmodes
       phase = m(i)*theta + n(i)*zeta

       R = R + Rc(i)*COS(phase)
       Z = Z + Zc(i)*SIN(phase)

    END DO

  END SUBROUTINE evaluate_RZ

  ! zeta is the NESCOIL angle, ie PHI=zeta/nfp where PHI is the geometric toroidal angle
  ! t is the parameter used in the polynomial fit of the Fourier coefficients. It ranges (typically) from [0,1]
  ! where t=0 returns the LCFS and t=1 returns the winding surface
  SUBROUTINE evaluate_xyz(t,theta,zeta,x,y,z)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in)  :: t, theta, zeta
    DOUBLE PRECISION, INTENT(out) :: x, y, z

    DOUBLE PRECISION :: R

    CALL evaluate_coefficients(t)
    CALL evaluate_RZ(theta, zeta, R, z)

    x = R*COS(zeta/nfp)
    y = R*SIN(zeta/nfp)

  END SUBROUTINE evaluate_xyz

  ! zeta is the NESCOIL angle, ie PHI=zeta/nfp where PHI is the geometric toroidal angle
  ! rho_in is a radial variable off the LCFS, approximately equal to the (euclidean) distance away from the LCFS.
  ! However the Fourier coefficients are fit over the variable t\in [0,1]. Hence we divide by the rhoscale to convert t=rho_in/rhoscale
  ! rho_in=0 will return the LCFS and rho_in=rhoscale will return the winding surface.
  ! rho_in+1 is an extrapolation of the VMEC radial coordinate.
  SUBROUTINE rhothetazeta2xyz(rho_in,theta_in,zeta_in,x_out,y_out,z_out)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: rho_in, theta_in, zeta_in
    DOUBLE PRECISION, INTENT(out) :: x_out, y_out, z_out

    CALL evaluate_xyz(rho_in/rhoscale, theta_in, zeta_in, x_out, y_out, z_out)

  END SUBROUTINE rhothetazeta2xyz


  SUBROUTINE xyz2rhothetazeta(&
       x_in,y_in,z_in,&
       rho_out,theta_out,zeta_out,&
       iflag_out)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: x_in, y_in, z_in
    DOUBLE PRECISION, INTENT(inout) :: rho_out, theta_out
    DOUBLE PRECISION, INTENT(out) :: zeta_out

    INTEGER, INTENT(out), OPTIONAL :: iflag_out

    DOUBLE PRECISION :: dRdrho, dZdrho, dRdtheta, dZdtheta
    DOUBLE PRECISION :: R_in, R_eval, Z_eval, dR, dZ, determinant
    DOUBLE PRECISION :: delta_rho, delta_theta
    DOUBLE PRECISION :: residual_norm

    INTEGER, PARAMETER :: max_iterations = 50

    DOUBLE PRECISION, PARAMETER :: absolute_tolerance = 1.0e-10
    DOUBLE PRECISION, PARAMETER :: relative_tolerance = 1.0e-12
    DOUBLE PRECISION, PARAMETER :: jacobian_tolerance = 1.0e-12
    DOUBLE PRECISION, PARAMETER :: pi2 = 6.2831853071795864769D0

    DOUBLE PRECISION :: convergence_tolerance, geometry_scale, jacobian_scale

    INTEGER :: iteration
    INTEGER :: status

    ! Status codes
    ! ------------
    ! iflag = 0   converged
    ! iflag = 1   iteration limit
    ! iflag = 2   singular/ill-conditioned Jacobian

    status = 1
    residual_norm = HUGE(1.0)

    R_in = SQRT(x_in*x_in + y_in*y_in)
    zeta_out = ATAN2(y_in,x_in)*get_field_period()  ! the field period should come from the surface_extender_mod and (hopefully) match the Biot-Savart field period

    ! initial guess for Newton search
    rho_out = MIN(MAX(rho_out,0.0),rhoscale)
    theta_out = MODULO(theta_out, pi2)
    ! theta_out=0.0 ! surely there is something more clever than this

    geometry_scale = MAX(1.0, SQRT(R_in*R_in + z_in*z_in))

    convergence_tolerance = absolute_tolerance + relative_tolerance*geometry_scale

    ! Newton step
    DO iteration = 0, max_iterations
       ! Evaluate the analytic Jacobian.
       CALL evaluate_RZ_derivatives(rho_out, theta_out, zeta_out, R_eval, Z_eval, dRdrho, dRdtheta, dZdrho, dZdtheta)
       dR = R_eval - R_in
       dZ = z_eval - z_in

       residual_norm = SQRT(dR*dR + dZ*dZ)
       IF (residual_norm <= convergence_tolerance) THEN
          status = 0
          EXIT
       END IF

       ! Do not take update max_iterations + 1.
       IF (iteration == max_iterations) EXIT

       determinant = dRdrho*dZdtheta - dRdtheta*dZdrho
       jacobian_scale = SQRT(dRdrho*dRdrho + dZdrho*dZdrho)*SQRT(dRdtheta*dRdtheta + dZdtheta*dZdtheta)

       IF (jacobian_scale <= TINY(1.0)) THEN
          status = 2
          EXIT
       END IF

       IF (ABS(determinant) <= &
            jacobian_tolerance*jacobian_scale) THEN
          status = 2
          EXIT
       END IF

       ! Solve J*delta = -residual.
       delta_rho = (-dZdtheta*dR + dRdtheta*dZ)/determinant
       delta_theta = (dZdrho*dR - dRdrho*dZ)/determinant

       rho_out = MAX(rho_out + delta_rho, 0.0)
       theta_out = MODULO(theta_out + delta_theta, pi2)

    END DO

    IF (PRESENT(iflag_out)) iflag_out = status

  END SUBROUTINE xyz2rhothetazeta


  SUBROUTINE evaluate_RZ_derivatives(rho_in,theta_in,zeta_in,R,Z,dRdrho,dRdtheta,dZdrho,dZdtheta)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: rho_in, theta_in, zeta_in
    DOUBLE PRECISION, INTENT(out) :: R,Z,dRdrho,dRdtheta,dZdrho,dZdtheta

    DOUBLE PRECISION :: t

    INTEGER  :: i
    DOUBLE PRECISION :: phase, cphase,sphase

    t = rho_in / rhoscale
    CALL evaluate_coefficients(t)

    R = 0.0
    Z = 0.0
    dRdrho   = 0.0
    dZdrho   = 0.0
    dRdtheta = 0.0
    dZdtheta = 0.0

    DO i = 1, nmodes
       phase = m(i)*theta_in + n(i)*zeta_in
       cphase = COS(phase)
       sphase = SIN(phase)

       R = R + Rc(i)*cphase
       Z = Z + Zc(i)*sphase

       dRdrho = dRdrho + dRc_dt(i)*cphase/rhoscale
       dZdrho = dZdrho + dZc_dt(i)*sphase/rhoscale

       dRdtheta = dRdtheta - m(i)*Rc(i)*sphase
       dZdtheta = dZdtheta + m(i)*Zc(i)*cphase
    END DO

  END SUBROUTINE evaluate_RZ_derivatives


  SUBROUTINE unload_surface_fit()

    ! gracefully deallocate arrays

    IF (ALLOCATED(m)) DEALLOCATE(m)
    IF (ALLOCATED(n)) DEALLOCATE(n)

    IF (ALLOCATED(Rc)) DEALLOCATE(Rc)
    IF (ALLOCATED(Zc)) DEALLOCATE(Zc)

    IF (ALLOCATED(dRc_dt)) DEALLOCATE(dRc_dt)
    IF (ALLOCATED(dZc_dt)) DEALLOCATE(dZc_dt)

    IF (ALLOCATED(Rc_coeff)) DEALLOCATE(Rc_coeff)
    IF (ALLOCATED(Zc_coeff)) DEALLOCATE(Zc_coeff)

    IF (ALLOCATED(R0)) DEALLOCATE(R0)
    IF (ALLOCATED(R1)) DEALLOCATE(R1)
    IF (ALLOCATED(Z0)) DEALLOCATE(Z0)
    IF (ALLOCATED(Z1)) DEALLOCATE(Z1)
    IF (ALLOCATED(Rc_ab)) DEALLOCATE(Rc_ab)
    IF (ALLOCATED(Zc_ab)) DEALLOCATE(Zc_ab)

  END SUBROUTINE unload_surface_fit

END MODULE surface_extender_mod
