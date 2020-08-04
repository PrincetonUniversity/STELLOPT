!-----------------------------------------------------------------------
!     Subroutine:    chisq_famus_bn
!     Authors:       J.C. Schmitt (Auburn/PPPL) (jcschmitt@auburn.edu)
!     Date:          2020
!     Description:   Chisq routine(s) for famus
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
!                    jtarget_famus_bn
!                    On entry, if niter is 0 or larger, then:
!                       increment mtargets, and
!                       assign targets, sigmas, and vals to the
!                       appropriate quantities from the target and
!                       sigma input arrays.
!
!-----------------------------------------------------------------------
      SUBROUTINE chisq_famus_bn(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------

      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_input_mod
      USE stellopt_vars, ONLY: mnprod_fps, mnprod_fds, famusin_filename
      USE mpi_params, ONLY: MPI_COMM_MYWORLD
!DEC$ IF DEFINED (FAMUS)
      USE famus_globals, ONLY: bn_famus => bn, Nzeta, Nteta, focusin, &
                  famus_called_from_ext_opt => called_from_ext_opt, &
                  famus_MPI_COMM_MYWORLD => MPI_COMM_MYWORLD, &
                  famus_ext => ext
!DEC$ ENDIF      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE

      REAL(rprec), INTENT(in)    :: target(mnprod_fps)
      REAL(rprec), INTENT(in)    :: sigma(mnprod_fps)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag

!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: iunit=5152, counter, ii, index_dot
      REAL(rprec), dimension(:), allocatable  :: bn_reshape

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      print *, '<----Entered chisq_famus with iflag=',iflag
      IF (iflag < 0) RETURN

!DEC$ IF DEFINED (FAMUS)
      IF (iflag == 1) THEN
          counter = 0
          DO ii = 1,(Nzeta * Nteta)
            IF (sigma(ii) < bigno) counter=counter +1
          END DO
          WRITE(iunit_out,'(A,2(2X,I7))') 'FAMUS_BN ', counter, 4
          WRITE(iunit_out,'(A)') 'TARGET  SIGMA  UNUSED  CHI'
      END IF

      IF (niter >= 0) THEN
         ! Reshape the FAMUS 'bn' (bn_famus) from a 2-D matrix to a 1-D array
         bn_reshape = reshape(bn_famus, (/ Nzeta * Nteta /) )
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
         write(6,*) '<----chisq_famus_bn: iflag=',iflag,' niter=',niter
         write(6,*) '<----chisq_famus_bn counting mtargets and/or assigning target_dex'
         IF (ANY(sigma < bigno)) THEN
            write(6,*) '<----chisq_famus_bn Nzeta=',Nzeta, ' Nteta=',Nteta
            write(6,*) '<----chisq_famus_bn mtargets on entry=',mtargets
            write(6,*) '<----chisq_famus_bn: Reading FAMUSIN namelist from ', famusin_filename

            index_dot = INDEX(famusin_filename,'.input')
            IF (index_dot .gt. 0)  famus_ext = famusin_filename(1:index_dot-1)

            famus_called_from_ext_opt = .True.
            call famus_initialize
            !call famus_read_namelist(famusin_filename)
            write(6,*) '<----chisq_famus_driver post famus_read_namelist Nzeta=',Nzeta, ' Nteta=',Nteta
!            CALL safe_open(iunit, iflag, TRIM(famusin_filename), 'old', 'formatted')
!
!            !CALL regcoil_read_input(iunit, iflag)
!            READ(iunit, nml=focusin, iostat=iflag)
!
!            close(iunit)
            IF (iflag < 0) THEN
               WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!!!!!!'
               WRITE(6,*) '  FAMUS Namelist (focusin) not found     '
               WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            ELSE
               write(6,*) '<---FAMUS namelist (focusin) found, iflag=',iflag
            END IF

            write(6,*) '<----chisq_famus_driver post namelist read: Nzeta=',Nzeta, ' Nteta=',Nteta
            ! Fill in the targets
            DO ii = 1, (Nzeta * Nteta)
               IF (sigma(ii) < bigno) THEN
                  mtargets = mtargets + 1
                  IF (niter == -2) THEN
                     target_dex(mtargets)=jtarget_famus_bn
                  END IF
               END IF
            END DO
            write(6,*) '<----chisq_famus_bn: mtargets on exit=',mtargets
         END IF
      END IF
!DEC$ ENDIF      

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_famus_bn
