!-----------------------------------------------------------------------
!     Subroutine:    chisq_focus_bn
!     Authors:       J.C. Schmitt (Auburn/PPPL) (jcschmitt@auburn.edu)
!     Date:          2020
!     Description:   Chisq routine(s) for FOCUS
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
!                    On entry, if niter is less than 0 the
!                    code should increment the mtargets value by
!                    the number of sigmas less than bigno.
!                    On entry, if niter is equlal to -2, the value of
!                    target_dex(mtargets) will be set to
!                    jtarget_focus_bn
!                    On entry, if niter is 0 or larger, then:
!                       increment mtargets, and
!                       assign targets, sigmas, and vals to the
!                       appropriate quantities from the target and
!                       sigma input arrays.
!
!-----------------------------------------------------------------------
      SUBROUTINE chisq_focus_bn(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------

      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_input_mod
      USE stellopt_vars, ONLY: mnprod_fds, focusin_filename
!DEC$ IF DEFINED (FOCUS)
      USE focus_globals, ONLY: bn_focus => bn, Nzeta, Nteta, focusin
!DEC$ ENDIF      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE

      REAL(rprec), INTENT(in)    :: target(mnprod_fds)
      REAL(rprec), INTENT(in)    :: sigma(mnprod_fds)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag

!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: iunit=5152, counter, ii
      REAL(rprec), dimension(:), allocatable  :: bn_reshape

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN

!DEC$ IF DEFINED (FOCUS)
      IF (iflag == 1) THEN
          counter = 0
          DO ii = 1,(Nzeta * Nteta)
            IF (sigma(ii) < bigno) counter=counter +1
          END DO
          WRITE(iunit_out,'(A,2(2X,I7))') 'FOCUS_BN ', counter, 4
          WRITE(iunit_out,'(A)') 'TARGET  SIGMA  UNUSED  CHI'
      END IF

      IF (niter >= 0) THEN
         ! Reshape the FOCUS 'bn' (bn_focus) from a 2-D matrix to a 1-D array
         bn_reshape = reshape(bn_focus, (/ Nzeta * Nteta /) )
         ! Now fill in the targets, sigmas and chi_sq
         DO ii = 1, (Nzeta * Nteta)
            !IF (sigma(ii) >= bigno) CYCLE
            IF (sigma(ii) < bigno) THEN
              mtargets = mtargets + 1
              targets(mtargets) = target(ii)
              sigmas(mtargets)  = sigma(ii)
              ! The value of the results is in the Bnormal_total_target variable
              vals(mtargets)    = bn_reshape(ii)
              IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3)') target(ii), &
                                    sigma(ii), 0.0, vals(mtargets)
            END IF
         END DO
      ELSE
         write(6,*) '<----chisq_focus_driver: iflag=',iflag,' niter=',niter
         write(6,*) '<----chisq_focus_driver counting mtargets and/or assigning target_dex'
         IF (ANY(sigma < bigno)) THEN
            write(6,*) '<----chisq_focus_driver Nzeta=',Nzeta, ' Nteta=',Nteta
            write(6,*) '<----chisq_focus_driver mtargets on entry=',mtargets
            ! Fill in the targets
            DO ii = 1, (Nzeta * Nteta)
               IF (sigma(ii) < bigno) THEN
                  mtargets = mtargets + 1
                  IF (niter == -2) THEN
                     target_dex(mtargets)=jtarget_focus_bn
                  END IF
               END IF
            END DO
            write(6,*) '<----chisq_focus_driver mtargets on exit=',mtargets
            write(6,*) '<----chisq_focus_driver: Reading FOCUSIN namelist from ', focusin_filename
            CALL safe_open(iunit, iflag, TRIM(focusin_filename), 'old', 'formatted')

            !CALL regcoil_read_input(iunit, iflag)
            READ(iunit, nml=focusin, iostat=iflag)

            close(iunit)
            IF (iflag < 0) THEN
               WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
               WRITE(6,*) '  FOUCS Namelist not found     '
               WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            END IF
         END IF
      END IF
!DEC$ ENDIF      

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_focus_bn
