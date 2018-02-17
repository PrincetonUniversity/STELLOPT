!-----------------------------------------------------------------------
!     Subroutine:    chisq_regcoil_chi2_b
!     Authors:       J.C. Schmitt (Auburn/PPPL) (jcschmitt@auburn.edu)
!     Date:          2017
!     Description:   Chisq routine(s) for REGCOIL.
!                    More description needed
!                    This is a template for the chisq routines.  In
!                    general all chisq routines should take a target
!                    variable, a sigma variable, and and error flag. On
!                    entry, if niter is less than 1 the
!                    code should simply increment the mtargets value by
!                    the number of sigmas less than bigno.  On entry, if
!                    iflag is set to a positive number the code should
!                    output to screen.  On entry, if iflag is set to
!                    zero the code should operate with no screen output.
!                    On exit, negative iflag terminates execution,
!                    positive iflag, indicates error but continues, and
!                    zero indicates the code has functioned properly.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_regcoil_chi2_b(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------

! JCS TO DO: Verify that all of these are necessary.
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_input_mod
      USE stellopt_vars, ONLY: nlambda_regcoil
      USE equil_vals, ONLY: curtor
!DEC$ IF DEFINED (REGCOIL)
      USE regcoil_input_mod 
      USE regcoil_variables
!DEC$ ENDIF      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      integer :: iunit
 
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
!DEC$ IF DEFINED (REGCOIL)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') &
         'REGCOIL CHI2_B ',1,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  DUMMY  CHI'
      IF (niter >= 0) THEN
        IF (sigma < bigno) THEN
           mtargets = mtargets + 1
           targets(mtargets) = target
           sigmas(mtargets)  = sigma
           vals(mtargets)    = sqrt(chi2_B_target)
           !   targets(mtargets) = 0.0
           !  sigmas(mtargets)  = bigno
           !  vals(mtargets)    = 0.0
           IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,0.0,vals(mtargets)
        ENDIF
      ELSE
         ! IF (sigma < bigno .and. myid == master) THEN
         IF (sigma < bigno) THEN
           write(6,'(a,i12)') '<---- niter=', niter
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_regcoil_chi2_b
           ! Read the regcoil namelist from the input."id_string" file
           ! WRITE(6,'(a,a)') '<---- id_string=', id_string
        
           CALL safe_open(iunit, iflag, TRIM('input.'//TRIM(id_string)), 'old', 'formatted')
           CALL read_regcoil_input(iunit, iflag)
           ! save an internal copy of the value of nlambda here (regcoil may
           ! overwrite it)
           nlambda_regcoil = nlambda
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
      END SUBROUTINE chisq_regcoil_chi2_b
