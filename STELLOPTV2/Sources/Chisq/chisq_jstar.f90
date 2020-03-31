!-----------------------------------------------------------------------
!     Subroutine:    chisq_jstar
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/28/2012
!     Description:   Calculate J_star at various values of ep/mu ranging
!                    from slightly above the trappd-passing boundary
!                    to slightly below the deeply trapped-forbidden
!                    boundary.  The parameters epl and epu determine
!                    the distance to these boundaries.  First we
!                    Fourier transform the Boozer magnetic field
!                    into real space then evaluate J_STAR on each
!                    surface.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_jstar(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils
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
      LOGICAL ::  lreset_s = .true.
      INTEGER ::  u, v, ik, ier, dex, nu_b, nv_b, nJstar
      REAL(rprec) :: bmin_global, bmax_global, sj, avg_Jstar, epsmu
      REAL(rprec), ALLOCATABLE :: xu_booz(:),xv_booz(:), bmin(:), bmax(:),&
                                  trapJs(:)
      REAL(rprec), ALLOCATABLE :: modb_temp(:,:),temp_arr(:,:)
      REAL(rprec), ALLOCATABLE :: modb_booz(:,:,:)
      
      INTEGER, PARAMETER     :: mskip = 8
      REAL(rprec), PARAMETER :: epl = 0.05_dp, epu = 0.05_dp
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      nu_b = 2 * mboz
      nv_b = 4 * nboz
      dex = COUNT(sigma < bigno)
      ik  = dex*NumJstar*(nu_b/mskip)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I8))') 'J_STAR ',ik,8
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VALS  AVG_JSTAR  TRAPJS  XU_BOOZ  #  V'
      IF (niter >= 0) THEN
         ALLOCATE(modb_booz(nu_b,nv_b,ns_b))
         ALLOCATE(xu_booz(nu_b),xv_booz(nv_b))
         ALLOCATE(bmin(nu_b),bmax(nu_b))
         ALLOCATE(trapJs(nu_b))
         ALLOCATE(temp_arr(nu_b,nv_b))
         ALLOCATE(modb_temp(nv_b,nu_b))
         FORALL(u=1:nu_b) xu_booz(u) = REAL(u-1)/REAL(2*(nu_b-1))                     !Fully around in theta
         FORALL(v=1:nv_b) xv_booz(v) = REAL(v-1)/REAL(nv_b)
         CALL mntouv(1,ns_b,mnboz_b,nu_b,nv_b,xu_booz,xv_booz,bmnc_b,&
                     ixm_b,ixn_b/nfp_b,modb_booz,0,1)
         IF (lasym_b) THEN
            CALL mntouv(1,ns_b,mnboz_b,nu_b,nv_b,xu_booz,xv_booz,bmns_b,&
                        ixm_b,ixn_b/nfp_b,modb_booz,1,0)
         END IF
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            ! Now get b_min and b_max
            DO u = 1, nu_b
               bmin(u) = MINVAL(modb_booz(u,:,ik))
               bmax(u) = MAXVAL(modb_booz(u,:,ik))
            END DO
            bmin_global = MINVAL(bmin)
            bmax_global = MAXVAL(bmax)
            temp_arr    = modb_booz(:,:,ik)
            modb_temp   = TRANSPOSE(temp_arr)
            ! Now get J-Star
            DO v = 1, NumJstar
               sj = REAL(v-1)/REAL(NumJstar-1)
               epsmu = bmin_global*(1.0_rprec + epl) + sj*(bmax_global*(1.0_rprec - epu) - bmin_global*(1.0_rprec + epl))
               CALL j_star(modb_temp, bmin, bmax, epsmu, trapJs, nv_b, nu_b)
               avg_Jstar = SUM(trapJs, mask=(trapJs > 0))
               nJstar    = COUNT(trapJs > 0)
               IF (nJstar > 0) avg_Jstar = avg_Jstar / nJstar
               WHERE(trapJs <= 0) trapJs = avg_Jstar
               DO u = 1, nu_b, mskip
                  ! Target all non-zero Jstars to (d Jstar/du) = 0
                  mtargets = mtargets + 1
                  targets(mtargets) = target(ik)
                  sigmas(mtargets)  = sigma(ik)
                  vals(mtargets)    = avg_Jstar - trapJs(u)
                  IF (iflag == 1) WRITE(iunit_out,'(6ES22.12E3,2(2X,I4))') target(ik),sigma(ik),vals(mtargets),avg_Jstar,trapJs(u),xu_booz(u),ik,v
               END DO
            END DO
         END DO
         DEALLOCATE(modb_booz)
         DEALLOCATE(xu_booz,xv_booz)
         DEALLOCATE(bmin,bmax)
         DEALLOCATE(trapJs)
         DEALLOCATE(modb_temp)
         DEALLOCATE(temp_arr)
      ELSE
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               lbooz(ik) = .TRUE.
               DO v = 1, NumJstar
                  DO u = 1, nu_b, mskip
                     mtargets = mtargets + 1
                     IF (niter == -2) target_dex(mtargets) = jtarget_Jstar
                  END DO
               END DO
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_jstar
