!-----------------------------------------------------------------------
!     Subroutine:    chisq_helicity
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/28/2012
!     Description:   Calculate ratio of magnetic energies in modes with
!                    undesirable helicities to the energy in modes with
!                    the desired helicity (bnorm). All harmonics target.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_helicity(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals
      USE read_boozer_mod
      USE vmec_input, ONLY: mpol, ntor
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
      LOGICAL :: booz_xform_initialized
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      IF (lasym) STOP 'ERROR: Helicity targeting requires lasym = .FALSE.'
      dex = COUNT(ABS(sigma) < bigno)*mnboz_b
      l_heli = NINT(REAL(helicity))
      k_heli = NINT(AIMAG(helicity))
      IF (niter >= 0) THEN   
         ! Get total number of modes for output
         IF (iflag == 1) THEN
            WRITE(iunit_out,'(A,2(2X,I8))') 'HELICITY_FULL ',dex,7
            WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VALS  BNORM  K  MBOZ  NBOZ'
            CALL FLUSH(iunit_out)
         END IF
         ! Now calculate chi_sq
         booz_xform_initialized = ALLOCATED(bmnc_b)
         DO ik = 1, nsd
            IF (ABS(sigma(ik)) >= bigno) CYCLE
            IF (booz_xform_initialized) THEN
               bmax  = MAXVAL(ABS(bmnc_b(1:mnboz_b,ik)))
            ELSE
               ! This case arises when mango is querying the sigmas and targets before the first equilibrium
               bmax = 1.0 ! Any value is fine.
            END IF
            sj = (real(ik,rprec) - 1.5_dp)/REAL((ns_b-1),rprec)            !!This is correct (SPH)
            bnorm = 0.0
            num0 = mtargets + 1
            DO mn = 1, mnboz_b
               mtargets = mtargets + 1
               targets(mtargets) = target(ik)
               vals(mtargets)    = 0
               sigmas(mtargets)  = bigno
               IF (booz_xform_initialized) THEN
                  n = ixn_b(mn)/nfp_b
                  m = ixm_b(mn)
                  bmn = bmnc_b(mn,ik)
               ELSE
                  ! m and n must match the initialization in setup_booz.f:
                  IF (mn <= (nboz+1)) THEN
                     m = 0
                     n = mn - 1
                  ELSE
                     m = FLOOR((mn + nboz - 1.0) / (nboz*2+1))
                     n = mn - nboz - 2 - (m-1) * (nboz*2+1) - nboz
                  END IF
                  bmn = 0 ! Any value is fine.
               END IF
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

               IF (lsym) THEN
                  bnorm = bnorm + bmn*bmn
                  CYCLE
               END IF

               sigmas(mtargets)  = 1
               vals(mtargets)    = bmn
            END DO

            bnorm = sqrt(bnorm)
            IF (bnorm == 0.0) bnorm = bmax

            IF (sigma(ik) < 0.0) THEN
               rad_sigma = 1
            ELSE IF (m < 3) THEN
               rad_sigma = sj
            ELSE IF (m == 3) THEN
               rad_sigma = sj**1.5
            ELSE
               rad_sigma = sj*sj
            END IF

            vals(num0:mtargets) = vals(num0:mtargets)/bnorm
            sigmas(num0:mtargets)  = ABS(sigma(ik))*rad_sigma*sigmas(num0:mtargets)
            
            IF (iflag ==1) THEN
               DO mn = num0,mtargets
                  WRITE(iunit_out,'(4ES22.12E3,3(1X,I5))') targets(mn),sigmas(mn),vals(mn),SQRT(bnorm),ik,ixm_b(mn-num0+1),ixn_b(mn-num0+1)/nfp_b
               END DO
            END IF
         END DO
      ELSE
         ! Consistency check
         mboz = MAX(6*mpol, 2, mboz)             
         nboz = MAX(2*ntor-1, 0, nboz)     
         ! CALCULATE mnboz_b becasue we don't know it yet (setup_booz.f)
         mnboz_b = (2*nboz+1)*(mboz-1) + (nboz + 1)
         DO ik = 1, nsd
            IF (ABS(sigma(ik)) < bigno) THEN
               lbooz(ik) = .TRUE.
               DO mn = 1, mnboz_b
                  mtargets = mtargets + 1
                  IF (niter == -2) target_dex(mtargets) = jtarget_helicity
               END DO
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_helicity
