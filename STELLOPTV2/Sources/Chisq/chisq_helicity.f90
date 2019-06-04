!-----------------------------------------------------------------------
!     Subroutine:    chisq_helicity
!     Authors:       S. Lazerson (lazerson@pppl.gov) modified by AB
!     Date:          06/28/2012 modified 06/2019
!     Description:   Calculate ratio of magnetic energies in modes with
!                    undesirable helicities to the energy in modes with
!                    the desired helicity (bnorm). 
!-----------------------------------------------------------------------
      SUBROUTINE chisq_helicity(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals
      USE read_boozer_mod
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      REAL(rprec), INTENT(in)    ::  target(nsd)
      REAL(rprec), INTENT(in)    ::  sigma(nsd)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!        lreset_s    Gets set to true if using R,PHI,Z Specification
!        ik          Dummy index
!        ier         Error Flag
!        te_val      Holds profile evaulation
!-----------------------------------------------------------------------
      LOGICAL :: lsym
      INTEGER :: dex, ik, mn, n, m, k_heli, l_heli, num0, n1, n2, i_save
      REAL(rprec) :: bnorm, bmax, bmn, rad_sigma, sj, val
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      IF (lasym) STOP 'ERROR: Helicity targeting requires lasym = .FALSE.'
      dex = COUNT(sigma < bigno)
      l_heli = NINT(REAL(helicity))
      k_heli = NINT(AIMAG(helicity))
      IF (niter >= 0) THEN   
         ! Get total number of modes for output
         IF (iflag == 1) THEN
            WRITE(iunit_out,'(A,2(2X,I8))') 'HELICITY ',dex,5
            WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VALS  BNORM  K'
            CALL FLUSH(iunit_out)
         END IF
         ! Now calculate chi_sq
         DO ik = 1, nsd
            val = 0.0_rprec
            IF (sigma(ik) >= bigno) CYCLE
            bmax  = MAXVAL(ABS(bmnc_b(1:mnboz_b,ik)))
            DO mn = 1, mnboz_b
               n = ixn_b(mn)/nfp_b
               m = ixm_b(mn)
               bmn = bmnc_b(mn,ik)
               !m_save(mn) = m
               !n_save(mn) = n
!               ! Target for minimization Bmn-s with helicities other than the one desired
!               ! General Helical Symmetry: mu - nv ~ Y(lu + kv) for integers Y != 0 (n,k in fp units)
!               ! HELICITY = (1,0)  !QA
!               !          = (0,1)  !QP
!               !          = (1,-1) !QH (l,k) (m*k+n*l)
               lsym = .FALSE.
               IF (k_heli == 0) THEN               !!quasi-axisymmetry
                  IF (n == 0) lsym = .TRUE. 
               ELSE IF (l_heli == 0) THEN          !!quasi-poloidal symmetry
                  IF (m == 0) lsym = .TRUE.       
               ELSE IF (MOD(m,l_heli) == 0) THEN   !!quasi-helical symmetry (lu + kv)
                  IF ((m*k_heli+n*l_heli) == 0) lsym = .TRUE.
               END IF
               
               !set normalizing field to b_00 rather than the symmetric modes
               IF ((n==0) .and. (m==0)) THEN
                  bnorm = bmn
               END IF

               !ignore symmetric modes
               IF (lsym) THEN
               !   bnorm = bnorm + bmn*bmn
                  CYCLE
               END IF

               val = val + bmn*bmn
            END DO
            val = sqrt(val)

            IF (bnorm == 0.0) bnorm = bmax
            mtargets = mtargets + 1
            vals(mtargets) = val/bnorm
            sigmas(mtargets)  = ABS(sigma(ik))
            targets(mtargets) = target(ik)
            IF (iflag ==1) THEN
               WRITE(iunit_out,'(4ES22.12E3,3(1X,I5))') targets(ik),sigmas(ik),vals(ik),bnorm,ik
            END IF
         END DO
         
      ELSE
         ! Consistency check
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               lbooz(ik) = .TRUE.
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_helicity
            END IF
         END DO
      END IF

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_helicity
