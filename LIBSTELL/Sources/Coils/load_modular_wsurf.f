      SUBROUTINE load_modular_wsurf (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE modular_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nvariables
      REAL(rprec) :: xvariables(*)
!-----------------------------------------------
      rmn_sf (1:numsurf) = xvariables (1:numsurf)
      zmn_sf (1:numsurf) = xvariables (numsurf+1:2*numsurf)
      nvariables = 2*numsurf

      END SUBROUTINE load_modular_wsurf
