      SUBROUTINE chisq_mercier(match, sigma,
     1           num, nrad, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: bigno
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nrad, nopt
      INTEGER :: num
      REAL(rprec) :: match(*), sigma(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jrad
      INTEGER, DIMENSION(nrad):: mercier_flag, pre_flag, post_flag
C-----------------------------------------------
      IF (ALL(ABS(sigma(3:nrad-1)) .ge. bigno)) RETURN

!
!     LOGIC TO AVOID ISOLATED RESONANCES
!
      IF (nopt. gt. 0) THEN
          mercier_flag = 0; pre_flag = 0; post_flag = 0                            !! MERCIERFLAG (RS)
          WHERE(match(3:nrad-1) < zero) mercier_flag(3:nrad-1) = 1                 !! MERCIERFLAG (RS)
          pre_flag(2:nrad-2) = mercier_flag(3:nrad-1)                              !! MERCIERFLAG (RS)
          post_flag(4:nrad) = mercier_flag(3:nrad-1)                               !! MERCIERFLAG (RS)
          mercier_flag = mercier_flag + pre_flag + post_flag                       !! MERCIERFLAG (RS)

         DO jrad = 3, nrad-1
            num = num + 1
            index_array(num) = ivar_mercier
            wegt(num) = sigma(jrad)
            chisq_target(num) = zero
            IF(mercier_flag(jrad) > 1) THEN                                    !! MERCIERFLAG (RS)
              chisq_match(num) = MIN(match(jrad), zero)
            ELSE
              chisq_match(num) = zero                                          !! MERCIERFLAG (RS)
            ENDIF
         ENDDO
      ELSE
         DO jrad = 3, nrad-1
            num = num + 1
            IF (nopt .eq. -2) chisq_descript(num) = 
     1                        descript(ivar_mercier)
         ENDDO
      ENDIF

      END SUBROUTINE chisq_mercier
