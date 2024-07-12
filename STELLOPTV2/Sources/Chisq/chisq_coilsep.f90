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
      USE windingsurface, ONLY : maxwindsurf
      USE stellopt_vars, ONLY : coil_splinesx, lwindsurf, coil_surf
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
      REAL(rprec)                                :: deltaphi, cdp, sdp, coilsep, s1, s2
      INTEGER                                    :: isurf, ic, itarg, n_uniq
      INTEGER, DIMENSION(maxwindsurf)            :: nuniq
      INTEGER, DIMENSION(nigroup,maxwindsurf)    :: iu
      LOGICAL                                    :: lmod

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN

      n_uniq = 0;  nuniq = 0;  iu = 0
      DO isurf=1,maxwindsurf ! For each valid winding surface
         IF (.NOT.lwindsurf(isurf)) CYCLE

         ! Count unique coils on this surface
         DO ic = 1, nigroup
            IF (coil_surf(ic).NE.isurf) CYCLE
            IF (ANY(coil_splinesx(ic,:) >-1)) THEN
               nuniq(isurf) = nuniq(isurf) + 1
               iu(nuniq(isurf),isurf) = ic
            ENDIF
         END DO !ic

         IF (nuniq(isurf).GT.0) n_uniq = n_uniq + nuniq(isurf) + 1
      ENDDO !isurf

      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'MIN_COIL_SEP ',n_uniq,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'C1 C2  TARGET  SIGMA  VAL  SLOC_1  SLOC_2'
      IF (sigma.ge.bigno) RETURN  !No targets if function is switched off.
      IF (.NOT.lcoil_geom) RETURN  !No targets if coils cannot be modified.

10    FORMAT(2I5,3ES22.12E3,2F11.4)

      IF (niter >= 0) THEN
         DO isurf=1,maxwindsurf
            IF (nuniq(isurf).EQ.0) CYCLE

            ! Generate all unique coils on this surface
            ALLOCATE(xyzuniq(npts_csep,3,nuniq(isurf)), ctarg(npts_csep,3))
            !iu = 0
            !DO ic=1,nigroup
            !   IF (ANY(coil_splinesx(ic,:) >-1)) THEN
            !      iu = iu + 1
                  ! Reconstruct coil geometry from spline data
            DO ic=1,nuniq(isurf)
               CALL spline_to_coil(iu(ic,isurf), npts_csep, &
                    xyzuniq(:,1,ic), xyzuniq(:,2,ic), xyzuniq(:,3,ic), lmod)
            ENDDO !ic
            !   END IF
            !END DO !ic
            !IF (iu .NE. n_uniq) STOP 'ERROR: Bad unique coil count in chisq_coilsep!'

            ! Determine distance of adjacent coils from each other
            DO itarg=1,nuniq(isurf)-1
               CALL get_coil_sep(xyzuniq(:,1,itarg), xyzuniq(:,2,itarg), xyzuniq(:,3,itarg), npts_csep, &
                    xyzuniq(:,1,itarg+1), xyzuniq(:,2,itarg+1), xyzuniq(:,3,itarg+1), npts_csep, &
                    coilsep, s1, s2)
               mtargets = mtargets + 1
               targets(mtargets) = targ
               sigmas(mtargets)  = sigma
               vals(mtargets) = MIN(coilsep, targ)  ! One-sided barrier
               !IF (iflag == 1) WRITE(iunit_out,10) itarg,itarg+1,targ,sigma,coilsep,s1,s2
               IF (iflag == 1) WRITE(iunit_out,10) iu(itarg,isurf),iu(itarg+1,isurf),targ,sigma,coilsep,s1,s2
            END DO !itarg

            ! Now target 1-1'
            ctarg(:,1) =  xyzuniq(:,1,1)
            ctarg(:,2) = -xyzuniq(:,2,1)
            ctarg(:,3) = -xyzuniq(:,3,1)
            mtargets = mtargets + 1
            CALL get_coil_sep(xyzuniq(:,1,1), xyzuniq(:,2,1), xyzuniq(:,3,1), npts_csep, &
                 ctarg(:,1), ctarg(:,2), ctarg(:,3), npts_csep, coilsep, s1, s2)
            targets(mtargets) = targ
            sigmas(mtargets)  = sigma
            vals(mtargets) = MIN(coilsep, targ)  ! One-sided barrier
            !IF (iflag == 1) WRITE(iunit_out,10) 1,2*n_uniq*nfp,targ,sigma,coilsep,s1,s2
            IF (iflag == 1) WRITE(iunit_out,10) iu(1,isurf),iu(1,isurf),targ,sigma,coilsep,s1,s2

            ! Now target N-N'
            deltaphi = twopi / REAL(nfp, rprec)
            cdp = COS(deltaphi);  sdp = SIN(deltaphi)
            ctarg(:,1) = cdp*xyzuniq(:,1,nuniq(isurf)) + sdp*xyzuniq(:,2,nuniq(isurf))
            ctarg(:,2) = sdp*xyzuniq(:,1,nuniq(isurf)) - cdp*xyzuniq(:,2,nuniq(isurf))
            ctarg(:,3) = -xyzuniq(:,3,nuniq(isurf))
            mtargets = mtargets + 1
            CALL get_coil_sep(xyzuniq(:,1,nuniq(isurf)), xyzuniq(:,2,nuniq(isurf)), xyzuniq(:,3,nuniq(isurf)), &
                 npts_csep, ctarg(:,1), ctarg(:,2), ctarg(:,3), npts_csep, coilsep, s1, s2)
            targets(mtargets) = targ
            sigmas(mtargets)  = sigma
            vals(mtargets) = MIN(coilsep, targ)  ! One-sided barrier
            !IF (iflag == 1) WRITE(iunit_out,10) n_uniq,n_uniq+1,targ,sigma,coilsep,s1,s2
            IF (iflag == 1) WRITE(iunit_out,10) iu(nuniq(isurf),isurf),iu(nuniq(isurf),isurf),targ,sigma,coilsep,s1,s2

            DEALLOCATE(xyzuniq, ctarg)
         ENDDO !isurf
      ELSE ! Just count targets
         IF (niter == -2) target_dex(mtargets+1:mtargets+n_uniq) = jtarget_coilsep
         mtargets = mtargets + n_uniq
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coilsep
