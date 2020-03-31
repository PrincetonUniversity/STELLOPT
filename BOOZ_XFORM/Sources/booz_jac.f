      SUBROUTINE booz_rzhalf(r, z, rodd, zodd, r12, z12, ohs, 
     1                    js, nznt, nrep)
      USE stel_kinds
      USE booz_params, ONLY: nv_boz
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: js, nznt, nrep
      REAL(rprec) :: ohs
      REAL(rprec), DIMENSION(nznt), INTENT(in) :: r, z, rodd, zodd
      REAL(rprec), DIMENSION(nznt), INTENT(out) :: r12, z12
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nvplane
      REAL(rprec) :: shalf, hs
C-----------------------------------------------
      hs = one/ohs
      shalf = SQRT(hs*ABS(js-c1p5))

      IF (nrep .eq. 1) THEN
         r12 = r + shalf*rodd
         z12 = z + shalf*zodd
      ELSE
         r12 = r
         z12 = z
      END IF

!     WRITE R,Z IN SPECIFIED TOROIDAL PLANE IN VMEC COORDINATES 
      nvplane = 1
!	nvplane = nv_boz/4
      !CALL WriteSurface(js, nvplane, nrep, r12, z12)

      END SUBROUTINE booz_rzhalf

      SUBROUTINE WriteSurface(js, nvplane, nrep, r12, z12)
      USE stel_kinds
      USE booz_params, ONLY: nv_boz
      USE booz_persistent, ONLY: nu3_b
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: js, nvplane, nrep
      REAL(rprec), DIMENSION(nv_boz, nu3_b), INTENT(in) :: r12, z12
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i
      CHARACTER*(256) :: file_name
C-----------------------------------------------
      IF (nrep .eq. 1) THEN
         file_name = "RZ_VMEC_VPLANE"
      ELSE 
         file_name = "RZ_BOOZER_VPLANE"
      END IF

      WRITE (file_name,'(a,i2,a)') TRIM(file_name), nvplane, ".txt"

      OPEN(unit=33,FILE=file_name,STATUS='REPLACE')

      DO I=1,nu3_b

        WRITE (33,100) I, r12(nvplane,I), z12(nvplane,I)

      END DO

 100  FORMAT(i5,1p,2E12.4)

      CLOSE (33)

      END SUBROUTINE WriteSurface
