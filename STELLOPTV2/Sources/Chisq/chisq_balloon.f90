!-----------------------------------------------------------------------
!     Subroutine:    chisq_balloon
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate difference between desired and actual
!                    ballooning growth rate.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_balloon(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals, ONLY: nrad, balloon_grate
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      REAL(rprec), INTENT(in)    ::  target(nsd)
      REAL(rprec), INTENT(in)    ::  sigma(nsd)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(in)    ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: i,j,ik, ntheta, nzeta, nlis
      REAL(rprec) :: grate
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ntheta = COUNT(balloon_theta >= 0.0)
      nzeta  = COUNT(balloon_zeta >= 0.0)
      nlis   = COUNT(sigma < bigno)
      ik = ntheta * nzeta * nlis
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I5.5))') 'BALLOON ',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  CHI  BALLOON_GRATE  THETA  ZETA  #'
      IF (niter >= 0) THEN
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               DO i = 1, ntheta
                  DO j = 1, nzeta
                     mtargets = mtargets + 1
                     targets(mtargets) = target(ik)
                     sigmas(mtargets)  = sigma(ik)
                     IF (ik <= nrad) THEN
                        vals(mtargets)    = 1+tanh(balloon_grate(ik,i,j)/sigma(ik))
                        vals(mtargets)    = MAX(target(ik),balloon_grate(ik,i,j)) ! Assume negative to be stable
                        grate             = balloon_grate(ik,i,j)
                     ELSE
                        vals(mtargets)   = target(ik)
                        grate            = 0.0
                     END IF
                     IF (iflag == 1) WRITE(iunit_out,'(6ES22.12E3,2X,I3.3)') target(ik),sigma(ik),vals(mtargets),grate,balloon_theta(i),balloon_zeta(j),ik
                  END DO
               END DO
            END IF
         END DO
      ELSE
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               DO i = 1, ntheta
                  DO j = 1, nzeta
                     mtargets = mtargets + 1
                     IF (niter == -2) target_dex(mtargets)=jtarget_balloon
                  END DO
               END DO
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_balloon
