!-----------------------------------------------------------------------
!     Subroutine:    chisq_rosenbrock
!     Authors:       J. Schmitt (jschmitt@pppl.gov)
!     Date:          2019
!     Description:   This is a chisq subroutine for a Rosenbrock-like
!                    test function.  It does not depend on any equlibrium
!                    values.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_rosenbrock(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(rosenbrock_dim)
      REAL(rprec), INTENT(in)    ::  sigma(rosenbrock_dim)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: ii
      REAL(rprec)  :: Rosenbrock_F
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'Rosenbrock_X',1,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  XVAL'
      IF (niter >= 0) THEN
         ! Now fill in the targets, sigmas and chi_sq
         DO ii = 1, rosenbrock_dim
            !IF (sigma(ii) >= bigno) CYCLE
            IF (sigma(ii) < bigno) THEN
               mtargets = mtargets + 1
               targets(mtargets) = target(ii)
               sigmas(mtargets)  = sigma(ii)
               ! Fill in vals with the Rosenbrock function values
               if (ii .eq. 1) THEN
                  Rosenbrock_F = 1.0 * (1.0 - Rosenbrock_X(ii))
               elseif (ii .gt. 1) THEN
                  ! Rosenbrock_F = 10.0 * (Rosenbrock_X(ii) - (Rosenbrock_X(ii-1))**2)
                  Rosenbrock_F = 10.0 * (Rosenbrock_X(ii) - Rosenbrock_X(ii-1)**2)
               END IF
               vals(mtargets)    = Rosenbrock_F
               IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,vals(mtargets)
            END IF
         END DO
      ELSE
         IF (ANY(sigma < bigno)) THEN
            ! Fill in the targets
            DO ii = 1, rosenbrock_dim
               IF (sigma(ii) < bigno) THEN
                  mtargets = mtargets + 1
                  IF (niter == -2) THEN
                     target_dex(mtargets)=jtarget_Rosenbrock_F
                  END IF
               END IF
            END DO
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_Rosenbrock
