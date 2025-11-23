!-----------------------------------------------------------------------
!     Subroutine:    chisq_bnmns
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          07/07/2025
!     Description:   Calculate the difference between the normal field
!                    on the equilibrium surface and zero.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_bnmns(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals, ONLY: bmns_normal_total, im_normal_total, &
            in_normal_total
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(-bnorm_nmax:bnorm_nmax,0:bnorm_mmax)
      REAL(rprec), INTENT(in)    ::  sigma(-bnorm_nmax:bnorm_nmax,0:bnorm_mmax)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(in)    ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: n, m, mn, mnmax
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      n = COUNT(sigma<bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I8.8))') 'BNMNS ',n,5
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  BNORMAL  N  M'
      IF (niter >= 0) THEN
         mnmax = SIZE(im_normal_total)
         DO m = 0, bnorm_mmax
            DO n = -bnorm_nmax, bnorm_nmax
               IF (sigma(n,m) < bigno) THEN
                  DO mn = 1, mnmax
                     IF (im_normal_total(mn) == m .and. in_normal_total(mn) == n) EXIT
                  END DO
                  mtargets = mtargets + 1
                  targets(mtargets) = target(n,m)
                  sigmas(mtargets)  = sigma(n,m)
                  IF ((mn > mnmax) .and. (niter == 0)) &
                     WRITE(6,*) 'WARNING: (',n,',',m,') Exceeds MNMAX in CHISQ_BNMNS'
                  IF (mn > mnmax) THEN
                        vals(mtargets) = target(n,m)
                  ELSE
                        vals(mtargets)    = bmns_normal_total(mn)
                  END IF
                  IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3,2(2X,I3.3))') &
                        target(n,m),sigma(n,m),vals(mtargets),n,m
               END IF
            END DO
         END DO
      ELSE
         lneed_bnormal = .TRUE.
         DO m = 0, bnorm_mmax
            DO n = -bnorm_nmax, bnorm_nmax
               IF (sigma(n,m) < bigno) THEN
                  mtargets = mtargets + 1
                  IF (niter == -2) target_dex(mtargets)=jtarget_bnmns
               END IF
            END DO
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_bnmns
