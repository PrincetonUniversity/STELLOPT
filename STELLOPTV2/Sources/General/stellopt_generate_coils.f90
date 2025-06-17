!-----------------------------------------------------------------------
!     Subroutine:    stellopt_generate_coils
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          06/13/2025
!     Description:   This subroutine computes a coils file in cartesian
!                    space from a spline representation.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_generate_coils(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_vars, ONLY: ncoils_max, nknots_coils_max, &
            rho_coil_kts, theta_coil_kts, zeta_coil_kts, &
            nw_coil, nh_coil, width_coil, height_coil
      USE stellopt_runtime, ONLY: proc_string
      USE read_wout_mod, ONLY: mnmax, ns, xm, xn, rmnc, zmns, isigng
      USE vmec_input, ONLY: extcur
      USE spline_coils_mod
      USE biotsavart, ONLY: write_coils_file
      USE stel_kinds, ONLY: rprec

!-----------------------------------------------------------------------
!     Input Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      INTEGER :: i,n,k, numcoilgroups
      INTEGER, PARAMETER :: nscoil = 128
      REAL(rprec) :: c1, c2, c3
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: tvec
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: rho, theta, zeta

!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      !-----------------------------------------------------------------
      !     Compute helpers
      !-----------------------------------------------------------------
      DO i = 1, ncoils_max
         DO k = 1, nknots_coils_max
            IF (rho_coil_kts(i,k)>0) n=k
         ENDDO
      ENDDO
      !n = MAXVAL(MAXLOC(rho_coil_kts,DIM=2,BACK=.TRUE.))
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
         WRITE(6,'(A,F7.3)') '         COIL WIDTH:  ',width_coil
         WRITE(6,'(A,F7.3)') '        COIL HEIGHT:  ',height_coil
      END IF

      !-----------------------------------------------------------------
      !     Load Splines
      !-----------------------------------------------------------------
      CALL init_spline_coils(nscoil, numcoilgroups, n, n+k, &
                              rho_coil_kts(1:numcoilgroups,1:n), &
                              theta_coil_kts(1:numcoilgroups,1:n), &
                              zeta_coil_kts(1:numcoilgroups,1:n))
      !-----------------------------------------------------------------
      !     Load Boundary
      !-----------------------------------------------------------------
      CALL init_boundary_spline_coils(mnmax,xm,-xn, &
                                      rmnc(:,ns),zmns(:,ns))

      !-----------------------------------------------------------------
      !     Create coils
      !-----------------------------------------------------------------
      CALL spline_to_coils(isigng)

      !-----------------------------------------------------------------
      !     Set the Current
      !-----------------------------------------------------------------
      CALL set_currents(numcoilgroups,extcur(1:numcoilgroups))

      !-----------------------------------------------------------------
      !     Make multi-filament
      !-----------------------------------------------------------------
      IF (nw_coil > 1 .or. nh_coil > 1) &
            CALL coils_to_multifilament(nw_coil,nh_coil, &
                                          width_coil,height_coil)
      !-----------------------------------------------------------------------
      !     Compute Curvature and Torsion
      !-----------------------------------------------------------------------
      CALL compute_coil_curvature(TRIM(proc_string))
      IF (lscreen) THEN
         CALL get_coil_curvature_avg(c1,c2,c3)
         WRITE(6,'(A)')            '        COIL CURVATURE:  '
         WRITE(6,'(A,3(2X,F7.3))') '              MIN/MEAN/MAX:  ',c3,c1,c2
         CALL get_coil_torsion_avg(c1,c2,c3)
         WRITE(6,'(A)')            '          COIL TORSION:  '
         WRITE(6,'(A,3(2X,F7.3))') '              MIN/MEAN/MAX:  ',c3,c1,c2
         WRITE(6,'(A)')      '-------------------------------------'
         CALL FLUSH(6)
      END IF


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