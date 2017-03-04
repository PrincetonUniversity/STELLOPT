!-----------------------------------------------------------------------
!     Subroutine:    stelltran_find_coeffs
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          07/25/2016
!     Description:   Determines the transport coefficients for the current
!                    system configuration
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_find_coeffs(itime)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_vars
      USE stelltran_equilutils
!-----------------------------------------------------------------------
!     Local Variables
!        ier          Error flag
!        ik           Indexing variable
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: itime
      INTEGER ::  n, j, ier, useless_n
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: coeff_init, s_init, useless_r
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      OPEN(UNIT=12,FILE='./coeffs/Xe_sfincs.txt',STATUS='OLD',ACTION='READ')
      ! Get rid of the values we don't care about from eariler times
      DO j=1,itime-1
            READ(12,*) useless_n
            ALLOCATE(useless_r(n))
            READ(12,*) useless_r
            READ(12,*) useless_r
            DEALLOCATE(useless_r)
      END DO
      ! read values we care about
      READ(12,*) n
      ALLOCATE(coeff_init(1:n),s_init(1:n))
      READ(12,*) s_init
      READ(12,*) coeff_init
      ! Spline interpolate onto
      DO j=1,prof_length
            CALL eval_prof_spline(n,s_init,coeff_init,Xe(j,1),Xe(j,2),ier)
      END DO
      DEALLOCATE(coeff_init)
      DEALLOCATE(s_init)
      CLOSE(12)
      OPEN(UNIT=12,FILE='./coeffs/Xi_sfincs.txt',STATUS='OLD',ACTION='READ')
      DO j=1,itime-1
            READ(12,*) useless_n
            ALLOCATE(useless_r(n))
            READ(12,*) useless_r
            READ(12,*) useless_r
            DEALLOCATE(useless_r)
      END DO
      READ(12,*) n
      ALLOCATE(coeff_init(1:n),s_init(1:n))
      READ(12,*) s_init
      READ(12,*) coeff_init
      DO j=1,prof_length
            CALL eval_prof_spline(n,s_init,coeff_init,Xi(j,1),Xi(j,2),ier)
      END DO
      DEALLOCATE(coeff_init)
      DEALLOCATE(s_init)
      CLOSE(12)
      OPEN(UNIT=12,FILE='./coeffs/De_sfincs.txt',STATUS='OLD',ACTION='READ')
      DO j=1,itime-1
            READ(12,*) useless_n
            ALLOCATE(useless_r(n))
            READ(12,*) useless_r
            READ(12,*) useless_r
            DEALLOCATE(useless_r)
      END DO
      READ(12,*) n
      ALLOCATE(coeff_init(1:n),s_init(1:n))
      READ(12,*) s_init
      READ(12,*) coeff_init
      DO j=1,prof_length
            CALL eval_prof_spline(n,s_init,coeff_init,De(j,1),De(j,2),ier)
      END DO
      DEALLOCATE(coeff_init)
      DEALLOCATE(s_init)
      CLOSE(12)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE stelltran_find_coeffs
