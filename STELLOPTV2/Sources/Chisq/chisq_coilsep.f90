!-----------------------------------------------------------------------
!     Subroutine:    chisq_coilsep
!     Authors:       J. Breslau (jbreslau@pppl.gov)
!     Date:          2017
!     Description:   Calculates the distance of closest approach between
!                    each modifiable coil and every other coil.
!     Warning:       This routine currently assumes all coils are modular,
!                    and none are self-symmetric.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coilsep(targ,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars, ONLY : coil_splinesx
      USE vmec_input, ONLY    : nfp, lasym
      IMPLICIT NONE
      INTRINSIC SIN, COS

!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      REAL(rprec), INTENT(in)    ::  targ, sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      REAL(rprec), PARAMETER                     :: twopi = 6.283185307179586476925286766559D0
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: xyzuniq
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE   :: ctarg
      REAL(rprec)                                :: deltaphi, cdp, sdp, rsgn
      INTEGER                                    :: ic, ifp, irefl, itarg, iu, n_uniq, ic0
      LOGICAL, DIMENSION(:), ALLOCATABLE         :: lmod

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'MIN. COIL SEP. '
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'C1 C2  SIGMA  VAL'
      IF (sigma.ge.bigno) RETURN  !No targets if function is switched off.
      IF (.NOT.lcoil_geom) RETURN  !No targets if coils cannot be modified.

      ! Count unique coils
      n_uniq = 0
      DO ic = 1, nigroup
         IF (ANY(coil_splinesx(ic,:) >-1)) n_uniq = n_uniq + 1
      END DO

      IF (niter >= 0) THEN
         ! Generate all unique coils
         ALLOCATE(xyzuniq(npts_csep,3,n_uniq), lmod(n_uniq), ctarg(npts_csep,3))
         iu = 0
         DO ic=1,nigroup
            IF (ANY(coil_splinesx(ic,:) >-1)) THEN
               iu = iu + 1
               ! Reconstruct coil geometry from spline data
               CALL spline_to_coil(ic, npts_csep, &
                    xyzuniq(:,1,iu), xyzuniq(:,2,iu), xyzuniq(:,3,iu), lmod(iu)) 
            END IF
         END DO !ic

         ! Determine distance from each other
         DO itarg=1,n_uniq-1
            ctarg = xyzuniq(:,:,itarg)
            DO ic=itarg+1,n_uniq
               mtargets = mtargets + 1
               CALL get_coil_sep(xyzuniq(:,1,ic), xyzuniq(:,2,ic), xyzuniq(:,3,ic), npts_csep, &
                    ctarg(:,1), ctarg(:,2), ctarg(:,3), npts_csep, vals(mtargets))
               targets(mtargets) = targ
               sigmas(mtargets)  = sigma
               IF (iflag == 1) WRITE(iunit_out,'(2I5,3ES22.12E3)') ic,itarg, targ,sigma,vals(mtargets)
            END DO !ic
         END DO !itarg

         ! Now target rotations & reflections of these coils
         deltaphi = twopi / REAL(nfp, rprec)
         rsgn = 1.0d0
         DO irefl=0,1
            DO itarg=1,n_uniq
               ctarg(:,3) = rsgn*xyzuniq(:,3,itarg)
               DO ifp=1,nfp+irefl-1
                  ! Reflect (iff irefl==1) and rotate to target field period
                  cdp = COS(ifp*deltaphi);  sdp = SIN(ifp*deltaphi)
                  ctarg(:,1) = cdp*xyzuniq(:,1,itarg) - rsgn*sdp*xyzuniq(:,2,itarg)
                  ctarg(:,2) = sdp*xyzuniq(:,1,itarg) + rsgn*cdp*xyzuniq(:,2,itarg)

                  ! Avoid redundant calcs. for reflections w/in same f.p.
                  IF ((irefl.eq.1).AND.(ifp.eq.1)) THEN
                     ic0 = itarg
                  ELSE
                     ic0 = 1
                  END IF

                  DO ic=ic0,n_uniq
                     mtargets = mtargets + 1
                     CALL get_coil_sep(xyzuniq(:,1,ic), xyzuniq(:,2,ic), xyzuniq(:,3,ic), npts_csep, &
                          ctarg(:,1), ctarg(:,2), ctarg(:,3), npts_csep, vals(mtargets))
                     targets(mtargets) = targ
                     sigmas(mtargets)  = sigma
                     IF (iflag == 1) WRITE(iunit_out,'(2I5,3ES22.12E3)')&
                          ic,itarg+2*n_uniq*(irefl*nfp + ifp), targ,sigma,vals(mtargets)
                  END DO !ic
               END DO !ifp
            END DO !itarg
            IF (lasym) EXIT
            rsgn = -1.0D0
         END DO !irefl

         DEALLOCATE(xyzuniq, lmod, ctarg)
      ELSE
         DO itarg=1,n_uniq-1
            DO ic=itarg+1,n_uniq
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_coilsep
            END DO !ic
         END DO !itarg

         DO irefl=0,1
            DO itarg=1,n_uniq
               DO ifp=1,irefl+nfp-1
                  IF ((irefl.eq.1).AND.(ifp.eq.1)) THEN
                     ic0 = itarg
                  ELSE
                     ic0 = 1
                  END IF
                  DO ic=ic0,n_uniq
                     mtargets = mtargets + 1
                     IF (niter == -2) target_dex(mtargets)=jtarget_coilsep
                  END DO !ic
               END DO !ifp
            END DO !itarg
            IF (lasym) EXIT
         END DO !irefl
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coilsep
