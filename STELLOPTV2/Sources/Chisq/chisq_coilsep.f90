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
      REAL(rprec)                                :: deltaphi, cdp, sdp, coilsep
      INTEGER                                    :: ic, itarg, iu, n_uniq
      LOGICAL                                    :: lmod

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
         ALLOCATE(xyzuniq(npts_csep,3,n_uniq), ctarg(npts_csep,3))
         iu = 0
         DO ic=1,nigroup
            IF (ANY(coil_splinesx(ic,:) >-1)) THEN
               iu = iu + 1
               ! Reconstruct coil geometry from spline data
               CALL spline_to_coil(ic, npts_csep, &
                    xyzuniq(:,1,iu), xyzuniq(:,2,iu), xyzuniq(:,3,iu), lmod) 
            END IF
         END DO !ic
         IF (iu .NE. n_uniq) STOP 'ERROR: Bad unique coil count in chisq_coilsep!'

         ! Determine distance of adjacent coils from each other
         DO itarg=1,n_uniq-1
           CALL get_coil_sep(xyzuniq(:,1,itarg), xyzuniq(:,2,itarg), xyzuniq(:,3,itarg), npts_csep, &
                 xyzuniq(:,1,itarg+1), xyzuniq(:,2,itarg+1), xyzuniq(:,3,itarg+1), npts_csep, coilsep)
            mtargets = mtargets + 1
            targets(mtargets) = targ
            sigmas(mtargets)  = sigma
            vals(mtargets) = MIN(coilsep, targ)  ! One-sided barrier
            IF (iflag == 1) WRITE(iunit_out,'(2I5,3ES22.12E3)') itarg,itarg+1,targ,sigma,coilsep
         END DO !itarg

         ! Now target 1-1'
         ctarg(:,1) =  xyzuniq(:,1,1)
         ctarg(:,2) = -xyzuniq(:,2,1)
         ctarg(:,3) = -xyzuniq(:,3,1)
         mtargets = mtargets + 1
         CALL get_coil_sep(xyzuniq(:,1,1), xyzuniq(:,2,1), xyzuniq(:,3,1), npts_csep, &
              ctarg(:,1), ctarg(:,2), ctarg(:,3), npts_csep, coilsep)
         targets(mtargets) = targ
         sigmas(mtargets)  = sigma
         vals(mtargets) = MIN(coilsep, targ)  ! One-sided barrier
         IF (iflag == 1) WRITE(iunit_out,'(2I5,3ES22.12E3)')&
              1,2*n_uniq*nfp,targ,sigma,coilsep

         ! Now target N-N'
         deltaphi = twopi / REAL(nfp, rprec)
         cdp = COS(deltaphi);  sdp = SIN(deltaphi)
         ctarg(:,1) = cdp*xyzuniq(:,1,n_uniq) + sdp*xyzuniq(:,2,n_uniq)
         ctarg(:,2) = sdp*xyzuniq(:,1,n_uniq) - cdp*xyzuniq(:,2,n_uniq)
         ctarg(:,3) = -xyzuniq(:,3,n_uniq)
         mtargets = mtargets + 1
         CALL get_coil_sep(xyzuniq(:,1,n_uniq), xyzuniq(:,2,n_uniq), xyzuniq(:,3,n_uniq), npts_csep, &
              ctarg(:,1), ctarg(:,2), ctarg(:,3), npts_csep, coilsep)
         targets(mtargets) = targ
         sigmas(mtargets)  = sigma
         vals(mtargets) = MIN(coilsep, targ)  ! One-sided barrier
         IF (iflag == 1) WRITE(iunit_out,'(2I5,3ES22.12E3)')&
              n_uniq,n_uniq+1,targ,sigma,coilsep

         DEALLOCATE(xyzuniq, ctarg)
      ELSE ! Just count targets
         IF (niter == -2) target_dex(mtargets+1:mtargets+n_uniq+1) = jtarget_coilsep
         mtargets = mtargets + n_uniq + 1
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coilsep
