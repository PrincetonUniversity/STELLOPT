!-----------------------------------------------------------------------
!     Subroutine:    stellopt_balloon
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/06/2012
!     Description:   This subroutine calculates the ballooning growth
!                    rate.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_balloon(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE stellopt_targets, ONLY: sigma_balloon, balloon_theta, balloon_zeta
      USE equil_vals, ONLY: balloon_grate, nrad
      ! COBRA LIBRARIES
      USE readin_data, ONLY: lscreen_cobra=>lscreen
      USE ballooning_data, ONLY: init_theta, init_zeta, l_geom_input, &
                                 l_tokamak_input, k_w, kth
      USE fmesh_quantities, ONLY: radios
      ! VMEC
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER ::  ier, i, j, ik, ntheta, nzeta, nlis
      INTEGER, ALLOCATABLE :: bsurf(:)
      REAL(rprec), ALLOCATABLE :: grate(:)
      REAL(rprec) :: t1, t2
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      l_geom_input = .TRUE.
      l_tokamak_input = .FALSE.
      k_w = 10
      kth = 1
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  BALLOONING CALCULATION  ------------------------'
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','parvmec','paravmec')
            lscreen_cobra = lscreen
            sigma_balloon(1) = bigno  ! Don't do the first surface (aka magaxis)
            ! First we need to initialize the COBRA variables
            ntheta = COUNT(balloon_theta >= 0.0)
            nzeta  = COUNT(balloon_zeta >= 0.0)
            nlis   = COUNT(sigma_balloon < bigno)
            IF (.not. ALLOCATED(bsurf)) ALLOCATE(bsurf(nlis))
            j=1
            DO ik = 1, SIZE(sigma_balloon,DIM=1)
               IF (sigma_balloon(ik) < bigno) THEN
                  bsurf(j) = ik
                  j = j + 1
               END IF
            END DO
            CALL order_input('',nlis,bsurf,iflag)
            IF (iflag .ne. 0) RETURN
            IF (.not. ALLOCATED(grate)) ALLOCATE(grate(nrad))
            IF (.not. ALLOCATED(radios)) ALLOCATE(radios(nrad))
            IF (.not. ALLOCATED(balloon_grate)) ALLOCATE(balloon_grate(nrad,ntheta,nzeta))
            CALL second0(t1)
            DO i = 1, ntheta
               DO j = 1, nzeta
                  init_theta = balloon_theta(i)
                  init_zeta  = balloon_zeta(j)
                  CALL get_ballooning_grate(grate)
                  CALL second0(t2)
                  IF (lscreen) THEN
                     WRITE(6,120) init_zeta, init_theta, t2-t1
                     WRITE(6,48)
                  END IF
                  balloon_grate(:,i,j) = grate(:)
               END DO
            END DO
         CASE('spec')
      END SELECT
      IF (lscreen) WRITE(6,'(a)') ' -------------------------  BALLOONING CALCULATION DONE  ---------------------'
      CALL FLUSH(6)
      RETURN
  48  format('====================================================')
 120  format(/,'ZETA0 = ', 1pe10.3,' THETA0 = ', 1pe10.3,&
               ' TIME IN COBRA CODE:',1pe10.2,' SEC')
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_balloon
