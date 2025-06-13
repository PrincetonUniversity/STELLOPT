!-----------------------------------------------------------------------
!     Subroutine:    stellopt_generate_coils
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          06/13/2025
!     Description:   This subroutine computes a coils file in cartesian
!                    space from a spline representation.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_generate_coils(lscreen)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_vars, ONLY: ncoils_max, nknots_coils_max, &
            rho_coil_kts, theta_coil_kts, zeta_coil_kts, &
            nw_coil, nh_coil, width_coil, height_coil
      USE stellopt_runtime, ONLY: proc_string
      USE read_wout_mod, ONLY: mnmax, ns, xm, xn, rmnc, zmns
      USE spline_coils_mod
      USE biotsavart, ONLY: write_coils_file
      USE stel_kinds, ONLY: rprec

!-----------------------------------------------------------------------
!     Input Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      INTEGER :: i,n,k, numcoilgroups
      INTEGER, PARAMETER :: nscoil = 128
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: tvec
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: rho, theta, zeta

!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------
      !-----------------------------------------------------------------
      !     Compute helpers
      !-----------------------------------------------------------------
      n = MAXVAL(MAXLOC(rho_coil_kts,DIM=2,BACK=.TRUE.))
      numcoilgroups = COUNT(ANY(rho_coil_kts>0,DIM=2))
      k=1
      ALLOCATE(tvec(n))

      !-----------------------------------------------------------------
      !     Screen Output
      !-----------------------------------------------------------------
      IF (lscreen) THEN
         WRITE(6,'(A)')      '--------  COIL GENERATION  ----------'
         WRITE(6,'(A,I3)')   '     CURRENT GROUPS:  ',numcoilgroups
         WRITE(6,'(A,I3)')   '              KNOTS:  ',n
         WRITE(6,'(A,I3)')   '      COIL SEGMENTS:  ',nscoil
         WRITE(6,'(A,I3)')   '    WIDTH FILAMENTS:  ',nw_coil
         WRITE(6,'(A,I3)')   '   HEIGHT FILAMENTS:  ',nh_coil
         WRITE(6,'(A,I7.3)') '         COIL WIDTH:  ',width_coil
         WRITE(6,'(A,I7.3)') '        COIL HEIGHT:  ',height_coil
         WRITE(6,'(A)')      '-------------------------------------'
      END IF

      !-----------------------------------------------------------------
      !     Load Splines
      !-----------------------------------------------------------------
      FORALL(i=1:n) tvec(i) = dble(i-1)/dble(n-1)
      CALL init_spline_coils(nscoil, numcoilgroups, n, n+k, &
                              rho_coil_kts(1:numcoilgroups,1:n), &
                              theta_coil_kts(1:numcoilgroups,1:n), &
                              zeta_coil_kts(1:numcoilgroups,1:n), &
                              tvec)
      !-----------------------------------------------------------------
      !     Load Boundary
      !-----------------------------------------------------------------
      CALL init_boundary_spline_coils(mnmax,xm,-xn, &
                                      rmnc(:,ns),zmns(:,ns))

      !-----------------------------------------------------------------
      !     Create coils
      !-----------------------------------------------------------------
      CALL spline_to_coils

      !-----------------------------------------------------------------
      !     Make multi-filament
      !-----------------------------------------------------------------
      IF (nw_coil > 1 .or. nh_coil > 1) &
            CALL coils_to_multifilament(nw_coil,nh_coil, &
                                          width_coil,height_coil)

      !-----------------------------------------------------------------
      !     Write coils file
      !-----------------------------------------------------------------
      CALL write_coils_file(TRIM(proc_string))

      !-----------------------------------------------------------------
      !     Deallocations
      !-----------------------------------------------------------------
      DEALLOCATE(tvec)

!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE stellopt_generate_coils