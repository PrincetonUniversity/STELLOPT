!-----------------------------------------------------------------------
!     Subroutine:    chisq_resjac
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/12/2013
!     Description:   This subroutine targets resonant jacobians.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_resjac(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE read_boozer_mod
      USE equil_vals
      USE EZspline_obj
      USE EZspline
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(nsd)
      REAL(rprec), INTENT(in)    ::  sigma(nsd)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER            :: ik, mn, j, ier, mn_norm, mn_res
      REAL(rprec)        :: ds, s, s_old, s_temp,&
                            gvalc, gvals, gval, gnorm, gtot, fmin, f_iotares
      TYPE(EZspline1_r8) :: GMNC_spl, GMNS_spl, GMNCn_spl, GMNSn_spl
      
      EXTERNAL fmin, f_iotares
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'RESJAC ',1,6
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  TARGET_IOTA  XM  XN'
      IF (niter >= 0) THEN
         ds = 1./ns_b
         ! Find basis mn
         DO mn = 1, mnboz_b
            IF (ixm_b(mn) == 0 .and. ixn_b(mn) == 0) THEN
               mn_norm = mn
               EXIT
            END IF
         END DO
         ! Allocate Basis Arrays
         IF (EZspline_allocated(GMNCn_spl)) CALL EZspline_free(GMNCn_spl,ier)
         IF (EZspline_allocated(GMNSn_spl)) CALL EZspline_free(GMNSn_spl,ier)
         CALL EZspline_init(GMNCn_spl,ns_b,bcs0,iflag)
         GMNCn_spl%isHermite = 1
         CALL EZspline_setup(GMNCn_spl,gmnc_b(mn_norm,:),ier)
         IF (lasym_b) THEN
            CALL EZspline_init(GMNSn_spl,ns_b,bcs0,iflag)
            GMNSn_spl%isHermite = 1
            CALL EZspline_setup(GMNSn_spl,gmns_b(mn_norm,:),ier)
         END IF
         
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            ! Find resonant mn
            mn_res = 0
            DO mn = 1, mnboz_b
               IF (xm_resjac(ik) == ixm_b(mn) .and. &
                   xn_resjac(ik)*nfp == ixn_b(mn)) THEN
                   mn_res = mn
                   EXIT
               END IF
            END DO
            ! Define resonant target
            iota_res_tgt = REAL(xn_resjac(ik)*nfp,rprec)/REAL(xm_resjac(ik),rprec)
            ! Check to see if we find a resonant surface
            s = fmin(0.0,1.0,f_iotares,0.0)
            ! Define Splines if Found
            gtot = 0.0
            IF (s > 0.0 .and. s < 1.0) THEN
               IF (EZspline_allocated(GMNC_spl)) CALL EZspline_free(GMNC_spl,ier)
               IF (EZspline_allocated(GMNS_spl)) CALL EZspline_free(GMNS_spl,ier)
               CALL EZspline_init(GMNC_spl,ns_b,bcs0,iflag)
               GMNC_spl%isHermite = 1
               CALL EZspline_setup(GMNC_spl,gmnc_b(mn_res,:),ier)
               IF (lasym_b) THEN
                  CALL EZspline_init(GMNS_spl,ns_b,bcs0,iflag)
                  GMNS_spl%isHermite = 1
                  CALL EZspline_setup(GMNS_spl,gmns_b(mn_res,:),ier)
               END IF
               s_old = 0.0
               DO j = 1, ns_b
                  s_temp = (j-1)/(ns_b-1)
                  s = fmin(s_temp,1.0,f_iotares,0.0)
                  IF (s == s_old .or. s == s_temp) CYCLE
                  gvalc = 0.0; gvals = 0.0
                  CALL EZspline_interp(GMNC_spl,s,gvalc,ier)
                  IF (lasym_b) CALL EZspline_interp(GMNS_spl,s,gvals,ier)
                  gval = gvalc + gvals
                  gvalc = 0.0; gvals = 0.0
                  CALL EZspline_interp(GMNCn_spl,s,gvalc,ier)
                  IF (lasym_b) CALL EZspline_interp(GMNSn_spl,s,gvals,ier)
                  gnorm = gvalc + gvals
                  s_old = s
                  gtot = gtot + (gval/gnorm)**2
               END DO
            END IF
            IF (gtot > 0.0) gtot = SQRT(gtot)
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = gtot
            IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3,2(2X,I4.4))') target(ik),sigma(ik),gtot,iota_res_tgt,xm_resjac(ik),xn_resjac(ik)
         END DO
         IF (EZspline_allocated(GMNC_spl)) CALL EZspline_free(GMNC_spl,ier)
         IF (EZspline_allocated(GMNS_spl)) CALL EZspline_free(GMNS_spl,ier)
         IF (EZspline_allocated(GMNCn_spl)) CALL EZspline_free(GMNCn_spl,ier)
         IF (EZspline_allocated(GMNSn_spl)) CALL EZspline_free(GMNSn_spl,ier)
      ELSE
         IF (ANY(sigma < bigno)) lbooz = .TRUE. ! We need all the surfaces (note that only up to ns_vmec will be calculated)
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_resjac
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_resjac
      
      FUNCTION f_iotares(x)
      USE stel_kinds, ONLY: rprec
      USE equil_utils, ONLY : get_equil_iota
      USE equil_vals, ONLY: iota_res_tgt
      IMPLICIT NONE
      REAL(rprec) :: x
      INTEGER     :: ier
      REAL(rprec) :: iota_eq, f_iotares
      
      CALL get_equil_iota(x,iota_eq,ier)
      f_iotares = ABS(iota_eq - iota_res_tgt)
      END FUNCTION
