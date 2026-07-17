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

  PUBLIC :: load_surface_fit
  PUBLIC :: unload_surface_fit
  PUBLIC :: rhothetazeta2xyz
  PUBLIC :: get_field_period
  
CONTAINS

  INTEGER FUNCTION get_field_period()
    IMPLICIT NONE
    get_field_period=nfp
  END FUNCTION get_field_period
  
  
  SUBROUTINE load_surface_fit(filename)
    USE netcdf

    IMPLICIT NONE

    CHARACTER(*), INTENT(in) :: filename

    INTEGER :: ncid
    INTEGER :: varid
    INTEGER :: dimid
    INTEGER :: ierr
    INTEGER :: fit_type_len

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

          DO k = 2, degrees_of_freedom
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

    CASE default

       STOP "Unknown fit_type"

    END SELECT

  END SUBROUTINE evaluate_coefficients
  
      
  SUBROUTINE evaluate_RZ(theta,zeta,R,Z)
    ! This routine must be called after an update of the coefficients (evaluate_coefficients(t)) if the radial coordinate changes    
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
       phase = m(i)*theta + n(i)*zeta*nfp

       R = R + Rc(i)*COS(phase)
       Z = Z + Zc(i)*SIN(phase)

    END DO

  END SUBROUTINE evaluate_RZ
        
  SUBROUTINE evaluate_xyz(t,theta,zeta,x,y,z)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in)  :: t, theta, zeta
    DOUBLE PRECISION, INTENT(out) :: x, y, z

    DOUBLE PRECISION :: R

    CALL evaluate_coefficients(t)
    CALL evaluate_RZ(theta, zeta, R, z)

    x = R*COS(zeta)
    y = R*SIN(zeta)

  END SUBROUTINE evaluate_xyz


  SUBROUTINE rhothetazeta2xyz(rho_in,theta_in,zeta_in,x_out,y_out,z_out)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: rho_in, theta_in, zeta_in
    DOUBLE PRECISION, INTENT(out) :: x_out, y_out, z_out

    CALL evaluate_xyz(rho_in/rhoscale, theta_in, zeta_in, x_out, y_out, z_out)

  END SUBROUTINE rhothetazeta2xyz

  SUBROUTINE unload_surface_fit()
    
    ! gracefully deallocate arrays
    
    IF (ALLOCATED(m)) DEALLOCATE(m)
    IF (ALLOCATED(n)) DEALLOCATE(n)
    
    IF (ALLOCATED(Rc)) DEALLOCATE(Rc)
    IF (ALLOCATED(Zc)) DEALLOCATE(Zc)
    
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
        



