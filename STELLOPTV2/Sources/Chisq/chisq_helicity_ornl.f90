!-----------------------------------------------------------------------
!     Subroutine:    chisq_helicity_ornl
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/28/2012
!     Description:   Calculate ratio of magnetic energies in modes with
!                    undesirable helicities to the energy in modes with
!                    the desired helicity (bnorm). By surf
!-----------------------------------------------------------------------
      SUBROUTINE chisq_helicity_ornl(target,sigma,niter,iflag)
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
      INTEGER :: dex, ik, mn, n, m, k_heli, l_heli, num0, n1, n2
      REAL(rprec) :: bnorm, bmax, bmn, rad_sigma, sj, val
      LOGICAL, ALLOCATABLE :: lmask(:), nlmask(:)
      REAL(rprec), ALLOCATABLE :: rad_wegt(:)
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
            WRITE(iunit_out,'(A,2(2X,I8))') 'HELICITY ',dex,4
            WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VALS  BNORM'
            CALL FLUSH(iunit_out)
         END IF
         ! Now calculate chi_sq
         ALLOCATE(lmask(mnboz_b),rad_wegt(mnboz_b),nlmask(mnboz_b))
         lmask = .FALSE.
         IF (k_heli == 0) THEN  !QAS
            WHERE (ixn_b == 0)  lmask = .TRUE.
         ELSE IF (l_heli == 0) THEN !QPS
            WHERE (ixm_b == 0) lmask = .TRUE.
         ELSE ! Helical symmetry
            WHERE (MOD(ixm_b,l_heli) == 0) lmask = .TRUE.
            WHERE ((ixm_b*k_heli + ixn_b*l_heli/nfp) /= 0) lmask = .FALSE.
         END IF
         nlmask = .TRUE.
         WHERE(lmask) nlmask = .FALSE.
         print "(a,i6,a,i7,a,l1)","chisq_helicity_ornl mnboz_b=",mnboz_b," nsd=",nsd," allocated(bmnc_b)=",allocated(bmnc_b)
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            bmax  = MAXVAL(ABS(bmnc_b(1:mnboz_b,ik)))
            sj = (real(ik,rprec) - 1.5_dp)/REAL((ns_b-1),rprec)            !!This is correct (SPH)
            val = 0.0
            bnorm = 0.0
            rad_wegt = 1
            IF (sigma(ik) < 0) THEN
               rad_wegt = sj*sj
               WHERE(ixm_b < 3) rad_wegt = sj
               WHERE(ixm_b == 3) rad_wegt = sj**1.5
            END IF
            bnorm = SUM(bmnc_b(:,ik)**2,MASK=lmask)
            val   = SUM((bmnc_b(:,ik)/rad_wegt)**2,MASK=nlmask)
!            DO mn = 1, mnboz_b
!               n = ixn_b(mn)/nfp_b
!               m = ixm_b(mn)
!               bmn = bmnc_b(mn,ik)
!               ! Target for minimization Bmn-s with helicities other than the one desired
!               ! General Helical Symmetry: mu - nv ~ Y(lu + kv) for integers Y != 0 (n,k in fp units)
!               ! HELICITY = (1,0)  !QA
!               !          = (0,1)  !QP
!               !          = (1,-1) !QH (l,k) (m*k+n*l)
!               lsym = .FALSE.
!               IF (k_heli == 0) THEN               !!quasi-axisymmetry
!                  IF (n == 0) lsym = .TRUE. 
!               ELSE IF (l_heli == 0) THEN          !!quasi-poloidal symmetry
!                  IF (m == 0) lsym = .TRUE.       
!               ELSE IF (MOD(m,l_heli) == 0) THEN   !!quasi-helical symmetry (lu + kv)
!                  IF ((m*k_heli+n*l_heli) == 0) lsym = .TRUE.
!               END IF
!               ! mimic's oak ridge system
!               rad_sigma = 1
!               IF (sigma(ik) < 0.0) THEN
!                  rad_sigma = 1
!               ELSE IF (m < 3) THEN
!                  rad_sigma = sj
!               ELSE IF (m == 3) THEN
!                  rad_sigma = sj**1.5
!               ELSE
!                  rad_sigma = sj*sj
!               END IF
!               IF (lsym) THEN
!                  bnorm = bnorm + bmn*bmn
!               ELSE
!                  val = val + bmn*bmn/rad_sigma
!               END IF
!            END DO
            IF (bnorm == 0.0) bnorm = bmax*bmax
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = ABS(sigma(ik))
            vals(mtargets)    = SQRT(ABS(val/bnorm))
            IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3)') targets(mtargets),sigmas(mtargets),vals(mtargets),SQRT(bnorm)
         END DO
         DEALLOCATE(lmask, rad_wegt, nlmask)
      ELSE
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               lbooz(ik) = .TRUE.
               IF (niter == -2) target_dex(mtargets) = jtarget_helicity
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_helicity_ornl
