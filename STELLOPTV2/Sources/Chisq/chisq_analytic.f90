!-----------------------------------------------------------------------
!     Subroutine:    chisq_analytic
!     Authors:       J.C. Schmitt (Auburn/PPPL) (jcschmitt@auburn.edu)
!     Date:          2018
!     Description:   Chisq routine(s) for ANALYTIC.
!                    This is a 'typical' chisq routine.  In
!                    general all chisq routines should take a target
!                    variable, a sigma variable, an iteration number 
!                    and error flag.
!
!                    IFLAG:
!                    On entry, If iflag is less than 1, the code returns
!                    with no further actions.
!                    On entry, if iflag is set to zero, the code should
!                    operate with no screen output.
!                    On entry, if iflag is set to a 1, the code should
!                    output to screen.
!                    On exit, negative iflag terminates execution,
!                    positive iflag, indicates error but continues, and
!                    zero indicates the code has functioned properly.
!                    
!                    NITER:
!                    On entry, if niter is less than 1 the
!                    code should increment the mtargets value by
!                    the number of sigmas less than bigno.
!                    On entry, if niter is equlal to -2, the value of
!                    target_dex(mtargets) will be set to
!                    jtarget_analytic
!                    On entry, if niter is 0 or larger, then:
!                       increment mtargets, and
!                       assign targets, sigmas, and vals to the
!                       appropriate quantities from the target and
!                       sigma input arrays.
!
!-----------------------------------------------------------------------
      SUBROUTINE chisq_analytic(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------

      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_input_mod
      USE stellopt_vars, ONLY: analytic_fcnt, analytic_eval
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE

      REAL(rprec), INTENT(in)    :: target(9)
      REAL(rprec), INTENT(in)    :: sigma(9)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag

!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: iunit, counter, ii

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN

      IF (iflag == 1) THEN
          counter = 0
          DO ii = 1,INT(analytic_fcnt)
            IF (sigma(ii) < bigno) counter=counter +1
          END DO
          WRITE(iunit_out,'(A,2(2X,I7))') 'ANALYTIC ', counter, 4
          WRITE(iunit_out,'(A)') 'TARGET  SIGMA  UNUSED  CHI'
      END IF

      IF (niter >= 0) THEN
         ! Now fill in the targets, sigmas and chi_sq
         DO ii = 1, INT(analytic_fcnt)
            !IF (sigma(ii) >= bigno) CYCLE
            IF (sigma(ii) < bigno) THEN
              mtargets = mtargets + 1
              targets(mtargets) = target(ii)
              sigmas(mtargets)  = sigma(ii)
              ! The value of the results is in the analytic_eval variable
              vals(mtargets)    = analytic_eval(ii)
              IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3)') target(ii), &
                                    sigma(ii), 0.0, vals(mtargets)
            END IF
         END DO
      ELSE
         IF (ANY(sigma < bigno)) THEN
            ! Fill in the targets
            DO ii = 1, INT(analytic_fcnt)
               IF (sigma(ii) < bigno) THEN
                  mtargets = mtargets + 1
                  IF (niter == -2) THEN
                     target_dex(mtargets)=jtarget_analytic
                  END IF
               END IF
            END DO
            CALL safe_open(iunit, iflag, TRIM('input.'//TRIM(id_string)), 'old', 'formatted')

            !CALL regcoil_read_input(iunit, iflag)
            READ(iunit, nml=analytic_nml, iostat=iflag)

            close(iunit)
            IF (iflag < 0) THEN
               WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
               WRITE(6,*) '  Analytic Namelist not found     '
               WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            END IF
         END IF
      END IF

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_analytic
