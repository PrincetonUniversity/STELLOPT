!-----------------------------------------------------------------------
!     Subroutine:    chisq_gs2_ptsm3d
!     Authors:       A. Bader (abader@engr.wisc.edu)
!     Date:          07/2018
!     Description:   This subroutine targets the energy transfer
!           between stable and unstable modes using the ptsm3d module
!           Geometric quantities are calculated in the gs2 framework
!-----------------------------------------------------------------------
      SUBROUTINE chisq_gs2_ptsm3d(target,sigma,niter,iflag)
      !     Libraries
      !-----------------------------------------------------------------------

      USE stellopt_runtime
      USE stellopt_targets  
      USE vmec2gs2_mod
      USE safe_open_mod
      !-----------------------------------------------------------------------
      !     Input/Output Variables
      !
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      INTEGER                    :: nalpha, nzgrid, vmec_option, i, j
      INTEGER                    :: iunit, ierr
      REAL(rprec)                :: s, zeta_center, periods
      LOGICAL                    :: verbose, should_print_gs2
      REAL(rprec)                :: s_used, q, shat, L_ref, B_ref
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: alpha, zeta
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bmag, gradpar, &
                          gds2, gds21, gds22, gbdrift, gbdrift0, cvdrift, &
                          cvdrift0, sqrt_jac

      !----------------BEGIN SUBROUTINE --------------

      !For now hard code the parameter settings - can be moved to input file or
      !calculated here later

      nalpha = 1
      nzgrid = 6400
      periods = 400
      s = 0.5
      zeta_center = 0.0
      vmec_option = 0
      verbose = .false.
      should_print_gs2 = .true.
      iunit = 7276 !Figure out what number I should actually use later

      ALLOCATE(alpha(nalpha))
      ALLOCATE(zeta(-nzgrid:nzgrid))
      ALLOCATE(bmag(nalpha, -nzgrid:nzgrid))
      ALLOCATE(gradpar(nalpha, -nzgrid:nzgrid))
      ALLOCATE(gds2(nalpha, -nzgrid:nzgrid))
      ALLOCATE(gds21(nalpha, -nzgrid:nzgrid))
      ALLOCATE(gds22(nalpha, -nzgrid:nzgrid))
      ALLOCATE(gbdrift(nalpha, -nzgrid:nzgrid))
      ALLOCATE(gbdrift0(nalpha, -nzgrid:nzgrid))
      ALLOCATE(cvdrift(nalpha, -nzgrid:nzgrid))
      ALLOCATE(cvdrift0(nalpha, -nzgrid:nzgrid))
      ALLOCATE(sqrt_jac(nalpha, -nzgrid:nzgrid))


      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'P2     ',1,3
      IF (niter >= 0) THEN

        write (*,*) 'GS2 PTSM3D TEST FUNCTION'
        CALL vmec2gs2('wout_'//TRIM(PROC_STRING)//'.nc', nalpha, nzgrid,&
             zeta_center, periods, s, vmec_option, verbose, s_used, q, shat, &
             L_ref, B_ref, alpha, zeta, bmag, gradpar, gds2, gds21, gds22, & 
             gbdrift, gbdrift0, cvdrift, cvdrift0, sqrt_jac)
        write (*,*) 'Calculated gs2 coordinates'
        IF (should_print_gs2) THEN
            CALL safe_open(iunit, ierr, 'gs2_stell_'//TRIM(PROC_STRING)//'.txt', &
                           'unknown', 'formatted')

            !Write the header
            DO i = 1,nalpha
                WRITE(iunit, "(A)") "&parameters"
                WRITE(iunit, "(A,2F12.7)") "s0 = ",s
                WRITE(iunit, "(A,2F12.7)") "minor_a = ",L_ref
                WRITE(iunit, "(A,5F12.7)") "Bref = ",B_ref
                WRITE(iunit, "(A,F12.7)") "q0 = ",ABS(q)
                WRITE(iunit, "(A,F12.7)") "shat = ",shat
                WRITE(iunit, "(A,I6)") "gridpoints = ",nzgrid
                WRITE(iunit, "(A,F12.7)") "n_pol = ",periods
                WRITE(iunit, "(A)") "/"
                DO j = -nzgrid,nzgrid
                    WRITE(iunit, "(9ES22.12E3)") gds2(i,j), gds21(i,j), &
                          gds22(i,j), bmag(i,j), sqrt_jac(i,j), cvdrift(i,j), &
                          cvdrift0(i,j), gradpar(i,j), zeta(j)
                END DO
            END DO
            CLOSE(iunit)
        END IF
            
             
        mtargets = mtargets + 1
        targets(mtargets) = target
        sigmas(mtargets) = sigma
        vals(mtargets) = 1.0
        IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,1.0

      ELSE
        IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_gs2_ptsm3d
        END IF
      END IF
      RETURN

      DEALLOCATE(alpha)
      DEALLOCATE(zeta)
      DEALLOCATE(bmag)
      DEALLOCATE(gradpar)
      DEALLOCATE(gds2)
      DEALLOCATE(gds21)
      DEALLOCATE(gds22)
      DEALLOCATE(gbdrift)
      DEALLOCATE(gbdrift0)
      DEALLOCATE(cvdrift)
      DEALLOCATE(cvdrift0)
      DEALLOCATE(sqrt_jac)

      END SUBROUTINE chisq_gs2_ptsm3d
