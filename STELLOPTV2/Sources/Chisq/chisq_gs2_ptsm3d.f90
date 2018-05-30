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
      !-----------------------------------------------------------------------
      !     Input/Output Variables
      !
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      INTEGER                    :: nalpha, nzgrid, vmec_option
      REAL(rprec)                :: s, zeta_center, periods
      LOGICAL                    :: verbose
      REAL(rprec)                :: s_used, q, shat, L_ref, B_ref
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: alpha, zeta
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bmag, gradpar, &
                          gds2, gds21, gds22, gbdrift, gbdrift0, cvdrift, &
                          cvdrift0

      !----------------BEGIN SUBROUTINE --------------

      !For now hard code the parameter settings - can be moved to input file or
      !calculated here later

      nalpha = 5
      nzgrid = 7
      periods = 10
      s = 0.5
      zeta_center = 0.0
      vmec_option = 0
      verbose = .false.

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


      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'P2     ',1,3
      IF (niter >= 0) THEN

        write (*,*) 'GS2 PTSM3D TEST FUNCTION'
        CALL vmec2gs2('wout_'//TRIM(PROC_STRING)//'.nc', nalpha, nzgrid,&
             zeta_center, periods, s, vmec_option, verbose, s_used, q, shat, &
             L_ref, B_ref, alpha, zeta, bmag, gradpar, gds2, gds21, gds22, & 
             gbdrift, gbdrift0, cvdrift, cvdrift0)
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

      END SUBROUTINE chisq_gs2_ptsm3d
