!-----------------------------------------------------------------------
!     Subroutine:    chisq_regcoil_lambda
!     Authors:       J.C. Schmitt (Auburn/PPPL) (jcschmitt@auburn.edu)
!     Date:          2017-2018
!     Description:   Chisq routine(s) for REGCOIL.
!                    This is a 'typical' chisq routine.  In
!                    general all chisq routines should take a target
!                    variable, a sigma variable, an iteration number 
!                    and error flag.
!
!                    IFLAG:
!                    On entry, If iflag is less than 0, the code returns
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
!                    jtarget_regcoil_lambda
!                    On entry, if niter is 0 or larger, then:
!                       increment mtargets, and
!                       assign targets, sigmas, and vals to the
!                       appropriate quantities from the target and
!                       sigma input arrays.
!
!-----------------------------------------------------------------------
      SUBROUTINE chisq_regcoil_lambda(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------

      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_input_mod
      USE stellopt_vars, ONLY: regcoil_nlambda, mnprod_x4_rcws
!DEC$ IF DEFINED (REGCOIL)
      USE regcoil_variables, ONLY:  lambda_target, nlambda, regcoil_nml
!DEC$ ENDIF      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE

      REAL(rprec), INTENT(in)    :: target(mnprod_x4_rcws)
      REAL(rprec), INTENT(in)    :: sigma(mnprod_x4_rcws)
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

!DEC$ IF DEFINED (REGCOIL)
      IF (iflag == 1) THEN
          counter = 0
          DO ii = 1,mnprod_x4_rcws
            IF (sigma(ii) < bigno) counter=counter +1
          END DO
          WRITE(iunit_out,'(A,2(2X,I7))') 'REGCOIL_RMS_K ', counter, 4
          WRITE(iunit_out,'(A)') 'TARGET  SIGMA  UNUSED  CHI'
      END IF

      IF (niter >= 0) THEN
         ! Now fill in the targets, sigmas and values
         DO ii = 1, mnprod_x4_rcws
            !IF (sigma(ii) >= bigno) CYCLE
            IF (sigma(ii) < bigno) THEN
              mtargets = mtargets + 1
              targets(mtargets) = target(ii)
              sigmas(mtargets)  = sigma(ii)
              ! The value of the results is in the lambda_target variable
              vals(mtargets)    = lambda_target
              IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3)') target(ii), &
                                    sigma(ii), 0.0, vals(mtargets)
            END IF
         END DO
      ELSE
         IF (ANY(sigma < bigno)) THEN
            ! Fill in the targets
            DO ii = 1, mnprod_x4_rcws
               IF (sigma(ii) < bigno) THEN
                  mtargets = mtargets + 1
                  IF (niter == -2) THEN
                     target_dex(mtargets)=jtarget_regcoil_lambda
                  END IF
               END IF
            END DO
            CALL safe_open(iunit, iflag, TRIM('input.'//TRIM(id_string)), 'old', 'formatted')

            !CALL regcoil_read_input(iunit, iflag)
            READ(iunit, nml=regcoil_nml, iostat=iflag)

            ! save an internal copy of the value of nlambda here (regcoil may
            ! overwrite it)
            regcoil_nlambda = nlambda
            close(iunit)
            IF (iflag < 0) THEN
               WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
               WRITE(6,*) '  REGCOIL Namelist not found     '
               WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            END IF
         END IF
      END IF
!DEC$ ENDIF      

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_regcoil_lambda
