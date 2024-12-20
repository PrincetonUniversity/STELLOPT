!-----------------------------------------------------------------------
!     Subroutine:    chisq_rosenbrock2D
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          2023
!     Description:   This is a 2D Rosenbrock function see:
!                    https://en.wikipedia.org/wiki/Rosenbrock_function
!-----------------------------------------------------------------------
      SUBROUTINE chisq_rosenbrock2d(target,sigma,niter,iflag)
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
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: ii
      REAL(rprec)  :: F
      REAL(rprec),PARAMETER  ::  A = 1.0
      REAL(rprec),PARAMETER  ::  B = 100.0 
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'ROSENBROCK2D',1,5
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  XVAL  A  B'
      IF (niter >= 0) THEN
         F = (A - xval)**2 + B*(yval - xval**2)**2
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         vals(mtargets)    = F
         IF (iflag == 1) WRITE(iunit_out,'(5ES22.12E3)') target,sigma,vals(mtargets),A,B
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_Rosenbrock2D
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_rosenbrock2d
