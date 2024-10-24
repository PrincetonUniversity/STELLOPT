!-----------------------------------------------------------------------
!     Subroutine:    chisq_bmin
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   This subroutine targets the difference between
!                    the minimum value of |B| on a surface in Boozer
!                    coordinates
!-----------------------------------------------------------------------
      SUBROUTINE chisq_bmin(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils
      USE read_boozer_mod
      
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
      INTEGER ::  u, v, ik, ier, nu_b, nv_b
      REAL(rprec) :: bmin_global
      REAL(rprec), ALLOCATABLE :: xu_booz(:),xv_booz(:), bmin(:)
      REAL(rprec), ALLOCATABLE :: modb_booz(:,:,:)

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik = count(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'B_MIN ',ik,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  NS'
      IF (niter >= 0) THEN
         ! First transform to s,u,v space
         nu_b = 8 * mboz
         nv_b = 8 * nboz
         ALLOCATE(bmin(nu_b))
         ALLOCATE(modb_booz(nu_b,nv_b,ns_b))
         ALLOCATE(xu_booz(nu_b),xv_booz(nv_b))
         FORALL(u=1:nu_b) xu_booz(u) = REAL(u-1)/REAL(2*(nu_b-1))                     !Fully around in theta
         FORALL(v=1:nv_b) xv_booz(v) = REAL(v-1)/REAL(nv_b)
         CALL mntouv(1,ns_b,mnboz_b,nu_b,nv_b,xu_booz,xv_booz,bmnc_b,&
                     ixm_b,ixn_b/nfp_b,modb_booz,0,1)
         IF (lasym_b) THEN
            CALL mntouv(1,ns_b,mnboz_b,nu_b,nv_b,xu_booz,xv_booz,bmns_b,&
                        ixm_b,ixn_b/nfp_b,modb_booz,1,0)
         END IF
         ! Now evaluate
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            DO u = 1, nu_b
               bmin(u) = MINVAL(modb_booz(u,:,ik))
            END DO
            bmin_global = MINVAL(bmin)
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = bmin_global
            IF (iflag == 1) WRITE(iunit_out,'(3(ES22.12E3),1X,I5)') targets(mtargets),sigmas(mtargets),vals(mtargets),ik
         END DO
         ! DEALLOCATE
         DEALLOCATE(modb_booz)
         DEALLOCATE(xu_booz,xv_booz)
         DEALLOCATE(bmin)
      ELSE
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               lbooz(ik) = .true.
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_bmin
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_bmin
