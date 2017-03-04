!-----------------------------------------------------------------------
!     Subroutine:    chisq_ne
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate difference between equilbrium electron
!                    density and target electron density
!-----------------------------------------------------------------------
      SUBROUTINE chisq_ne(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils
      USE EZspline_obj
      USE EZspline
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(nprof)
      REAL(rprec), INTENT(in)    ::  sigma(nprof)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!        lreset_s    Gets set to true if using R,PHI,Z Specification
!        ik          Dummy index
!        ier         Error Flag
!        dex         Length of ne_aux_s array
!        ne_val      Holds profile evaulation
!-----------------------------------------------------------------------
      LOGICAL ::  lreset_s = .true.
      INTEGER ::  ik, ier, dex, num_avg
      REAL(rprec) :: ne_val, u_val, s_temp, avg_ne, avg_ne_val
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno_ne)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'NE ',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'R  PHI  Z  S  TARGET  SIGMA  NE'
      dex = MINLOC(ne_aux_s(2:),DIM=1)
      IF (niter >= 0) THEN
         IF (ANY(s_ne > 0)) lreset_s = .false.
         ! Normalize Targets
         IF (ANY(sigma_ne_line < bigno_ne)) THEN
            avg_ne = 0.0
            s_temp = 3.0
            num_avg = 0
            DO ik = 1, nprof
               IF (sigma(ik) >= bigno_ne) CYCLE
               !IF (target(ik) > avg_ne) avg_ne = target(ik)
               ier = 0
               CALL get_equil_s(r_ne(ik),phi_ne(ik),z_ne(ik),s_ne(ik),ier,u_val)
               ! Old Way
               !IF ((s_ne(ik) < s_temp) .and. (s_ne(ik) <= 1.0)) THEN
               !   s_temp = s_ne(ik)
               !   avg_ne = target(ik)
               !END IF
               ! New Way
               IF ((s_ne(ik) < 0.05) .and. (s_ne(ik) <= 1.0)) THEN
                  avg_ne = avg_ne + target(ik)
                  num_avg = num_avg + 1
               END IF
            END DO
            IF (num_avg > 0) THEN
               avg_ne = avg_ne / num_avg
            ELSE
               dex = MINLOC(s_ne,DIM=1,MASK=(s_ne > 0.0))
               avg_ne = target(dex)
            END IF
            ! Normalize profile
            ier = 0
            s_temp = 0.0
            CALL get_equil_ne(s_temp,TRIM(ne_type),avg_ne_val,ier)
         ELSE
            avg_ne = 1.0
            avg_ne_val = 1.0
         END IF
         ! Calc values
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno_ne) CYCLE
            ! GET s if necessary
            IF (lreset_s) THEN
               ier = 0
               CALL get_equil_s(r_ne(ik),phi_ne(ik),z_ne(ik),s_ne(ik),ier,u_val)
            END IF
            IF (s_ne(ik) <= 1.0 .and. s_ne(ik) >= 0.0) THEN
               ier = 0
               CALL get_equil_ne(s_ne(ik),TRIM(ne_type),ne_val,ier)
            ELSE
               ne_val = 0.0
            END IF
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)/avg_ne
            sigmas(mtargets)  = sigma(ik)/avg_ne
            vals(mtargets)    = ne_val/avg_ne_val
            IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3)') r_ne(ik),phi_ne(ik),z_ne(ik),s_ne(ik),target(ik)/avg_ne,sigma(ik)/avg_ne,ne_val/avg_ne_val
         END DO
         IF (lreset_s) s_ne(:) = -1.0
      ELSE
         DO ik = 1, nprof
            IF (sigma(ik) < bigno_ne) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_ne
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_ne
