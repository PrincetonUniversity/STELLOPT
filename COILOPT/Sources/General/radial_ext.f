      SUBROUTINE radial_ext (ncoil, theta, dr_ext)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE modular_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ncoil
      REAL(rprec) :: theta, dr_ext

      dr_ext = 0
      IF ((theta.le.pi/2).or.(theta.ge.(3*pi)/2))
     1  dr_ext = r_ext(ncoil)*COS(theta)

      END SUBROUTINE radial_ext
