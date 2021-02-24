!-----------------------------------------------------------------------
!     Subroutine:    chisq_qsfp
!     Authors:       E. Rodriguez (eduardor@princeton.edu)
!     Date:          20/02/2021
!     Description:   Calculate the deviation of quasisymmetry as defined
!                    by so-called Paul's form of quasisymmetry, f_P.
!
!-----------------------------------------------------------------------
      SUBROUTINE chisq_qsfp(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE equil_utils
      USE stellopt_targets
      USE vmec_main,          ONLY: hs
      USE read_wout_mod,      ONLY: pres, ns, num_per=>nfp, bmnc_vmec=>bmnc, ntor_vmec=>ntor, mpol_vmec=>mpol, mnmode_vmec=>mnmax, &
                                xm_nyq, xn_nyq, gmnc_vmec=>gmnc, vp_vmec=>vp, phi_vmec=>phi, &
                                iotas, bsupumnc, bsupvmnc, bsubumnc, bsubvmnc, isigng, amin=>Aminor, Rmax=>Rmajor

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
!
!-----------------------------------------------------------------------
      INTEGER     :: nu, nv, iu, iv, is
      real(rprec) :: du, dv, bbar, normalization, psi0, G, I, F, l_heli, k_heli, fp_val
      real(rprec), allocatable :: B(:,:), dBdtheta(:,:), dBdphi(:,:), d2Bdtheta2(:,:),  &
                     d2Bdphi2(:,:), d2Bdthetadphi(:,:), Bsuptheta(:,:), Bsupphi(:,:), &
                     Bsubtheta(:,:), Bsubphi(:,:), dBsupthetadtheta(:,:), dBsupthetadphi(:,:), &
                     dBsupphidtheta(:,:), dBsupphidphi(:,:), BdotgradB(:,:), BcrossgradpsidotgradB(:,:), &
                     xu(:), xv(:), J(:,:), angle(:), fp_local(:,:), cos_angle(:), sin_angle(:)

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lasym) STOP 'ERROR: Helicity targeting requires lasym = .FALSE.'
      l_heli = NINT(REAL(helicity))
      k_heli = num_per*NINT(AIMAG(helicity))
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'QSFP',1,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  QSFP'
      IF (niter >= 0) THEN
            nu = 8 * mpol_vmec
            IF (ntor_vmec > 0) THEN
                nv = 8 * ntor_vmec
            ELSE
                nv = 2
            END IF
            psi0 = phi_vmec(ns)/pi2
            ALLOCATE(B(nu,nv))
            ALLOCATE(dBdtheta(nu,nv))
            ALLOCATE(dBdphi(nu,nv))
            ALLOCATE(d2Bdphi2(nu,nv))
            ALLOCATE(d2Bdtheta2(nu,nv))
            ALLOCATE(d2Bdthetadphi(nu,nv))
            ALLOCATE(Bsubphi(nu,nv))
            ALLOCATE(Bsubtheta(nu,nv))
            ALLOCATE(Bsupphi(nu,nv))
            ALLOCATE(Bsuptheta(nu,nv))
            ALLOCATE(dBsupphidphi(nu,nv))
            ALLOCATE(dBsupphidtheta(nu,nv))
            ALLOCATE(dBsupthetadphi(nu,nv))
            ALLOCATE(dBsupthetadtheta(nu,nv))
            ALLOCATE(BdotgradB(nu,nv))
            ALLOCATE(BcrossgradpsidotgradB(nu,nv))
            ALLOCATE(J(nu,nv))
            ALLOCATE(fp_local(nu,nv))
            ALLOCATE(angle(mnmode_vmec))
            ALLOCATE(cos_angle(mnmode_vmec))
            ALLOCATE(sin_angle(mnmode_vmec))
            ALLOCATE(xu(nu),xv(nv))
            FORALL(iu=1:nu) xu(iu) = real(iu-1)/real(nu)
            FORALL(iv=1:nv) xv(iv) = real(iv-1)/real(nv)
            du = xu(2)-xu(1)
            dv = xv(2)-xv(1)
            B = 0_rprec
            J = 0_rprec
            dBdtheta = 0_rprec
            dBdphi = 0_rprec
            d2Bdphi2 = 0_rprec
            d2Bdtheta2 = 0_rprec
            d2Bdthetadphi = 0_rprec
            Bsubphi = 0_rprec
            Bsubtheta = 0_rprec
            Bsupphi = 0_rprec
            Bsuptheta = 0_rprec
            dBsupphidphi = 0_rprec
            dBsupphidtheta = 0_rprec
            dBsupthetadphi = 0_rprec
            dBsupthetadtheta = 0_rprec
            BdotgradB = 0_rprec
            BcrossgradpsidotgradB = 0_rprec
            normalization = 0
            bbar = 0
            DO is=2, nprof
            	IF (sigma(is) >= bigno) CYCLE
                DO iu=1,nu
                    DO iv=1,nv
                        angle = pi2*(real(xm_nyq(1:mnmode_vmec))*xu(iu)-real(xn_nyq(1:mnmode_vmec))*xv(iv))
                        cos_angle = cos(angle)
                        sin_angle = sin(angle)
                        dBdtheta(iu,iv) = sum(-xm_nyq(1:mnmode_vmec)*bmnc_vmec(1:mnmode_vmec,is)*sin_angle)
                        dBdphi(iu,iv) = sum(xn_nyq(1:mnmode_vmec)*bmnc_vmec(1:mnmode_vmec,is)*sin_angle)
                        d2Bdtheta2(iu,iv) = sum(-xm_nyq(1:mnmode_vmec)*xm_nyq(1:mnmode_vmec)*bmnc_vmec(1:mnmode_vmec,is)*cos_angle)
                        d2Bdphi2(iu,iv) = sum(-xn_nyq(1:mnmode_vmec)*xn_nyq(1:mnmode_vmec)*bmnc_vmec(1:mnmode_vmec,is)*cos_angle)
                        d2Bdthetadphi(iu,iv) = sum(xm_nyq(1:mnmode_vmec)*xn_nyq(1:mnmode_vmec)*bmnc_vmec(1:mnmode_vmec,is)*cos_angle)
                        B(iu,iv) = sum(bmnc_vmec(1:mnmode_vmec,is)*cos_angle)
                        Bsuptheta(iu,iv) = sum(bsupumnc(1:mnmode_vmec,is)*cos_angle)
                        Bsupphi(iu,iv) = sum(bsupvmnc(1:mnmode_vmec,is)*cos_angle)
                        Bsubtheta(iu,iv) = sum(bsubumnc(1:mnmode_vmec,is)*cos_angle)
                        Bsubphi(iu,iv) = sum(bsubvmnc(1:mnmode_vmec,is)*cos_angle)
                        dBsupthetadtheta(iu,iv) = sum(-xm_nyq(1:mnmode_vmec)*bsupumnc(1:mnmode_vmec,is)*sin_angle)
                        dBsupthetadphi(iu,iv) = sum(xn_nyq(1:mnmode_vmec)*bsupumnc(1:mnmode_vmec,is)*sin_angle)
                        dBsupphidtheta(iu,iv) = sum(-xm_nyq(1:mnmode_vmec)*bsupvmnc(1:mnmode_vmec,is)*sin_angle)
                        dBsupphidphi(iu,iv) = sum(xn_nyq(1:mnmode_vmec)*bsupvmnc(1:mnmode_vmec,is)*sin_angle)
                        J(iu,iv) = sum(gmnc_vmec(1:mnmode_vmec,is)*isigng*cos_angle)/psi0
                    END DO
                END DO
                G = sum(Bsubphi)/nu/nv
                I = sum(Bsubtheta)/nu/nv
                F = (l_heli*G-k_heli*I)/(l_heli*iotas(is)+k_heli)

                BdotgradB = Bsuptheta*dBdtheta+Bsupphi*dBdphi
                BcrossgradpsidotgradB =(dBdtheta*Bsubphi-dBdphi*Bsubtheta)/J
                fp_local = BcrossgradpsidotgradB-F*BdotgradB;
                
                ! Find appropriate normalisation
                bbar = sum(B*B*J)/sum(J)
                normalization = bbar**1.5
                fp_val = sqrt(sum(fp_local*fp_local*J)/sum(J))/normalization

                ! Output value
                mtargets = mtargets + 1
                targets(mtargets) = target(is)
                sigmas(mtargets)  = sigma(is)
                vals(mtargets)    = fp_val
                IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target(is),sigma(is),fp_val
            end do
      ELSE
        DO is = 2, nprof
            IF (sigma(is) < bigno) THEN
                mtargets = mtargets + 1
                IF (niter == -2) target_dex(mtargets) = jtarget_qsfp
            END IF
        END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_qsfp
