!-----------------------------------------------------------------------
!     Function:      out_fieldlines
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/2/2011
!     Description:   Save output from field line following while running.
!-----------------------------------------------------------------------
      SUBROUTINE out_fieldlines(phi,q)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_background
      USE pies_fieldlines
!-----------------------------------------------------------------------
!     Input Parameters
!          phi          Length of field line in phi
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(inout) :: phi
      DOUBLE PRECISION, INTENT(inout) :: q(ninteqs)
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      INTEGER             :: jint, ik
      DOUBLE PRECISION    :: pi2
      DOUBLE PRECISION    :: rho2(0:k)
!-----------------------------------------------------------------------
!     Begin Function
!-----------------------------------------------------------------------
      jint = INT(phi/dphi+0.5)
      rho2 = 0.0
      DO ik = 0,k-1
         xsiln(ik,jint) = q(2*ik+1)
         etaln(ik,jint) = q(2*ik+2)
         rho2(ik) = xsiln(ik,jint)*xsiln(ik,jint) + etaln(ik,jint)*etaln(ik,jint)
      END DO
      WHERE (rho2 >= 1.0) hitwal = 1
      IF (MOD(jint,1000)==0) THEN
         CALL backspace_out(6,6)
         !WRITE(6,'(a)',ADVANCE='no') CHAR(8)
         !WRITE(6,'(a)',ADVANCE='no') CHAR(8)
         !WRITE(6,'(a)',ADVANCE='no') CHAR(8)
         !WRITE(6,'(a)',ADVANCE='no') CHAR(8)
         !WRITE(6,'(a)',ADVANCE='no') CHAR(8)
         !WRITE(6,'(a)',ADVANCE='no') CHAR(8)
         WRITE(6,'(a,i3,a)',ADVANCE='no') '[',INT(100.*REAL(jint)/nintw),']%'
         CALL FLUSH(6)
      END IF
      phi = (jint +1)*dphi
!-----------------------------------------------------------------------
!     End Function
!-----------------------------------------------------------------------
      END SUBROUTINE out_fieldlines
