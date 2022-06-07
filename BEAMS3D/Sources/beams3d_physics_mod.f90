!-----------------------------------------------------------------------
!     Module:        beams3d_physics_mod
!     Authors:       M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          07/20/2012
!     Description:   This module contains the physics subroutines for
!                    manipulating particles and finding when to do so.
!-----------------------------------------------------------------------
MODULE beams3d_physics_mod

      !-----------------------------------------------------------------
      !     Libraries
      !-----------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime, ONLY: lneut, pi, pi2, dt, lverb, ADAS_ERR, &
                                 dt_save, lbbnbi, weight, ndt, &
                                 ndt_max, npoinc, lendt_m, te_col_min, &
                                 NION, NI_AUX_M, NI_AUX_Z
      USE beams3d_lines, ONLY: R_lines, Z_lines, PHI_lines, &
                               myline, mytdex, moment, ltherm, &
                               nsteps, nparticles, vll_lines, &
                               moment_lines, mybeam, mycharge, myZ, &
                               mymass, myv_neut, rand_prob, &
                               cum_prob, tau, &
                               epower_prof, ipower_prof, &
                               end_state, fact_crit, fact_pa, &
                               fact_vsound, fact_coul, fact_kick, &
                               ns_prof1, ns_prof2, ns_prof3, ns_prof4, &
                               ns_prof5, my_end
      USE beams3d_grid, ONLY: BR_spl, BZ_spl, delta_t, BPHI_spl, &
                              MODB_spl, MODB4D, &
                              phimax, S4D, X4D, Y4D, TE4D, NE4D, TI4D, ZEFF4D, &
                              nr, nphi, nz, rmax, rmin, zmax, zmin, &
                              phimin, eps1, eps2, eps3, raxis, phiaxis,&
                              zaxis, U4D, &
                              hr, hp, hz, hri, hpi, hzi, &
                              B_kick_min, B_kick_max, E_kick, freq_kick, &
                              plasma_mass, NI5D
      USE EZspline_obj
      USE EZspline
      USE adas_mod_parallel
      USE mpi_params 

      !-----------------------------------------------------------------
      !     Module PARAMETERS
      !-----------------------------------------------------------------

      DOUBLE PRECISION, PRIVATE, PARAMETER :: electron_mass = 9.10938356D-31 !m_e
      DOUBLE PRECISION, PRIVATE, PARAMETER :: e_charge      = 1.60217662E-19 !e_c
      DOUBLE PRECISION, PRIVATE, PARAMETER :: sqrt_pi       = 1.7724538509   !pi^(1/2)
      DOUBLE PRECISION, PRIVATE, PARAMETER :: inv_sqrt2     = 0.7071067812   !1/sqrt(2)
      DOUBLE PRECISION, PRIVATE, PARAMETER :: mpome         = 5.44602984424355D-4 !e_c
      DOUBLE PRECISION, PRIVATE, PARAMETER :: inv_dalton    = 6.02214076208E+26 ! 1./AMU [1/kg]
      DOUBLE PRECISION, PRIVATE, PARAMETER :: inv_cspeed    = 3.3356409520E-09 ! 1./c [s/m]
      DOUBLE PRECISION, PRIVATE, PARAMETER :: zero          = 0.0D0 ! 0.0
      DOUBLE PRECISION, PRIVATE, PARAMETER :: half          = 0.5D0 ! 1/2
      DOUBLE PRECISION, PRIVATE, PARAMETER :: one           = 1.0D0 ! 1.0

      !-----------------------------------------------------------------
      !     SUBROUTINES
      !        BEAMS3D_PHYISCS(t,q) Slowing down and pitch angle
      !        BEAMS3D_FOLLOW_NEUT(t,q) Neutral particle following
      !-----------------------------------------------------------------
      CONTAINS

      !-----------------------------------------------------------------
      !     Function:      beams3d_physics
      !     Authors:       S. Lazerson (lazerson@pppl.gov)
      !                    M. McMillan (matthew.mcmillan@my.wheaton.edu)
      !     Date:          07/13/2012
      !     Description:   Particle slowing down and pitch angle
      !                    scattering.
      !-----------------------------------------------------------------
      SUBROUTINE beams3d_physics(t, q)
         !--------------------------------------------------------------
         !     Input Parameters
         !          t          Location along fieldline in t
         !          q            (q(1),q(2),q(3),q(4)) = (R,phi,Z,vll)
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(inout) :: t
         DOUBLE PRECISION, INTENT(inout) :: q(4)
         !--------------------------------------------------------------
         !     Local Variables
         !        tau_spit The Spitzer characteristic time
         !        tau_spit_inv Helper 1/tau_spit
         !        v_crit       Critical velocity (transition to ion slowing)
         !        coulomb_log  Couloumb Logarithm (see NRL) 
         !--------------------------------------------------------------
         INTEGER        :: ier
         DOUBLE PRECISION    :: r_temp, phi_temp, z_temp, vll, te_temp, ne_temp, ti_temp, speed, newspeed, &
                          zeta, sigma, zeta_mean, zeta_o, v_s, tau_inv, tau_spit_inv, &
                          reduction, dve,dvi, tau_spit, v_crit, coulomb_log, te_cube, &
                          inv_mymass, speed_cube, vcrit_cube, vfrac, modb, s_temp, &
                          vc3_tauinv, vbeta, zeff_temp
         DOUBLE PRECISION :: Ebench  ! for ASCOT Benchmark
         ! For splines
         INTEGER :: i,j,k, l
         REAL*8 :: xparam, yparam, zparam
         INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
         REAL*8 :: fval(1)

         !--------------------------------------------------------------
         !     Begin Subroutine
         !--------------------------------------------------------------
      
         ier      = 0

         ! Setup position in a vll arrays
         r_temp   = q(1)
         phi_temp = MODULO(q(2), phimax)
         IF (phi_temp < 0) phi_temp = phi_temp + phimax
         z_temp   = q(3)
         vll      = q(4)

         ! Initialize values
         te_temp  = 0; ne_temp  = 0; ti_temp  = 0; zeff_temp=1;
         speed = 0; reduction = 0

         tau_spit_inv = 0.0; v_crit   = 0.0; coulomb_log = 15
         tau_inv = 10.0; vcrit_cube = 0.0; vc3_tauinv = 0

         ! Check that we're inside the domain then proceed
         !CALL EZspline_isInDomain(BR_spl,r_temp,phi_temp,z_temp,ier)
         IF ((r_temp >= rmin-eps1) .and. (r_temp <= rmax+eps1) .and. &
             (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
             (z_temp >= zmin-eps3) .and. (z_temp <= zmax+eps3)) THEN
!         IF (ier == 0) THEN
            ! Get the gridpoint info (this is possible since all grids are the same)
            i = MIN(MAX(COUNT(raxis < r_temp),1),nr-1)
            j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
            k = MIN(MAX(COUNT(zaxis < z_temp),1),nz-1)
            xparam = (r_temp - raxis(i)) * hri(i)
            yparam = (phi_temp - phiaxis(j)) * hpi(j)
            zparam = (z_temp - zaxis(k)) * hzi(k)
            ! Evaluate the Splines
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            MODB4D(1,1,1,1),nr,nphi,nz)
            modb = fval(1)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            TE4D(1,1,1,1),nr,nphi,nz)
            te_temp = max(fval(1),zero)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            NE4D(1,1,1,1),nr,nphi,nz)
            ne_temp = max(fval(1),zero)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            TI4D(1,1,1,1),nr,nphi,nz)
            ti_temp = max(fval(1),zero)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            ZEFF4D(1,1,1,1),nr,nphi,nz)
            zeff_temp = max(fval(1),one)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            S4D(1,1,1,1),nr,nphi,nz)
            s_temp = fval(1)

            !-----------------------------------------------------------
            !  Helpers
            !     v_s       Local Sound Speed
            !     speed     Total particle speed
            !-----------------------------------------------------------
            te_cube = te_temp * te_temp * te_temp
            inv_mymass = 1/mymass
            v_s = fact_vsound*sqrt(ti_temp)
            speed = sqrt( vll*vll + 2*moment*modb*inv_mymass )
            vbeta = max(ABS(speed-v_s)*inv_cspeed,1E-6)

            !-----------------------------------------------------------
            !  Calculate Coulomb Logarithm (NRL pg. 35)
            !     te in eV and ne in cm^-3
            !-----------------------------------------------------------
            IF ((te_temp > te_col_min).and.(ne_temp > 0)) THEN
               coulomb_log = 35 - log( zeff_temp*fact_coul*sqrt(ne_temp*1E-6/te_temp)/(vbeta*vbeta))
!               IF (te_temp < 10*myZ*myZ) THEN
!                  coulomb_log = 23 - log( myZ*sqrt(ne_temp*1E-6/(te_cube) )   )
!               ELSE
!                  coulomb_log = 24 - log( sqrt(ne_temp*1E-6)/(te_temp) )
!               END IF
               coulomb_log = max(coulomb_log,one)
               ! Callen Ch2 pg41 eq2.135 (fact*Vtherm; Vtherm = SQRT(2*E/mass) so E in J not eV)
               v_crit = fact_crit*SQRT(te_temp)
               vcrit_cube = v_crit*v_crit*v_crit
               tau_spit = 3.777183D41*mymass*SQRT(te_cube)/(ne_temp*myZ*myZ*coulomb_log)  ! note ne should be in m^-3 here
               tau_spit_inv = (1.0D0)/tau_spit
               vc3_tauinv = vcrit_cube*tau_spit_inv
            END IF

            !-----------------------------------------------------------
            !  Viscouse Velocity Reduction
            !     v_s       Local Sound Speed
            !     speed     Total particle speed
            !     dve       Speed change due to electron slowing down 
            !     dvi       Speed change due to ion slowing down 
            !     reduction Total change in speed
            !     newspeed  New total speed
            !     vfrac     Ratio between new and old speed (helper) 
            !-----------------------------------------------------------
            dve   = speed*tau_spit_inv
            dvi   = vc3_tauinv/(speed*speed)
            reduction = dve + dvi
            newspeed = speed - reduction*dt
            vfrac = newspeed/speed

            !-----------------------------------------------------------
            !  Thermalize particle or adjust vll and moment
            !  Fowler et al. NF 1990 30 (6) 997--1010
            !-----------------------------------------------------------
            IF (newspeed < v_s) THEN  ! Thermalize
               dve       = dve/reduction
               dvi       = dvi/reduction
               reduction = (speed - v_s)/dt  ! Thermalize
               dve       = dve*reduction
               dvi       = dvi*reduction
               newspeed = speed - reduction*dt
               ltherm = .true.
               vfrac = newspeed/speed
               vll = vfrac*vll
               moment = vfrac*vfrac*moment
               q(4) = vll
               RETURN
            END IF
            l = MAX(MIN(CEILING(SQRT(s_temp)*ns_prof1),ns_prof1),1)
            epower_prof(mybeam,l) = epower_prof(mybeam,l) + mymass*dve*dt*speed*weight(myline)
            ipower_prof(mybeam,l) = ipower_prof(mybeam,l) + mymass*dvi*dt*speed*weight(myline)
            vll = vfrac*vll
            moment = vfrac*vfrac*moment
            speed = newspeed

           !------------------------------------------------------------
           !  Pitch Angle Scattering
           !------------------------------------------------------------
           speed_cube = 2*vc3_tauinv*zeff_temp*fact_pa*dt/(speed*speed*speed) ! redefine as inverse
           zeta_o = vll/speed   ! Record the current pitch.
           CALL gauss_rand(1,zeta)  ! A random from a standard normal (1,1)
           sigma = sqrt( ABS((1.0D0-zeta_o*zeta_o)*speed_cube) ) ! The standard deviation.
           zeta_mean = zeta_o *(1.0D0 - speed_cube )  ! The new mean in the distribution.
           zeta = zeta*sigma + zeta_mean  ! The new pitch angle.
           !!!The pitch angle MUST NOT go outside [-1,1] nor be NaN; but could happen accidentally with the distribution.
           zeta = MIN(MAX(zeta,-0.999D+00),0.999D+00)
           !IF (ABS(zeta) >  0.999D+00) zeta =  SIGN(0.999D+00,zeta)
           vll = zeta*speed

           !------------------------------------------------------------
           !  Kick Model Scattering (old)
           !------------------------------------------------------------
           !IF (modb>=B_kick_min .and. modb<=B_kick_max) THEN
           !   zeta_o = vll/speed   ! Record the current pitch.
           !   zeta = zeta_o-zeta_o*(one-zeta_o*zeta_o)*dt*fact_kick*SQRT(ne_temp)/(modb*modb)
           !   vll = zeta*speed
           !END IF

           !------------------------------------------------------------
           !  Kick Model Scattering (new Energy, vll constant)
           !------------------------------------------------------------
           IF (modb>=B_kick_min .and. modb<=B_kick_max) THEN
              zeta_o = vll/speed   ! Record the current pitch.
              speed = speed*SQRT(one + fact_kick*modb*(1-zeta_o*zeta_o)*dt/SQRT(ne_temp))
           END IF

           !------------------------------------------------------------
           !  Final Moment and vll update (return q(4))
           !------------------------------------------------------------
           moment = half*mymass*(speed*speed - vll*vll)/modb
           q(4) = vll

         END IF

         RETURN

      END SUBROUTINE beams3d_physics

      !-----------------------------------------------------------------
      !     Function:      beams3d_follow_neut
      !     Authors:       S. Lazerson (lazerson@pppl.gov)
      !                    M. McMillan (matthew.mcmillan@my.wheaton.edu)
      !     Date:          12/01/2018
      !     Description:   Follows a particle into the plasma,
      !                    calculates deposition, and follows to wall
      !                    or plasma domain.
      !-----------------------------------------------------------------
      SUBROUTINE beams3d_follow_neut(t, q)
         !--------------------------------------------------------------
         !     Modules
         !--------------------------------------------------------------
         USE beams3d_grid
         USE beams3d_lines, ONLY: myline,xlast,ylast,zlast
         USE beams3d_runtime, ONLY: lvessel, to3, lplasma_only, &
                                    lvessel_beam, lsuzuki
         USE wall_mod, ONLY: collide, uncount_wall_hit

         !--------------------------------------------------------------
         !     Input Parameters
         !          t          Location along fieldline in t
         !          q            (q(1),q(2),q(3),q(4)) = (R,phi,Z,vll)
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(inout) :: t
         DOUBLE PRECISION, INTENT(inout) :: q(4)

         !--------------------------------------------------------------
         !     Local Parameters
         !--------------------------------------------------------------
         INTEGER, PARAMETER :: num_depo = 256
         DOUBLE PRECISION, PARAMETER :: dl = 5D-3
         DOUBLE PRECISION, PARAMETER :: stepsize(3)=(/0.25,0.05,0.01/)

         !--------------------------------------------------------------
         !     Local variables
         !--------------------------------------------------------------
         LOGICAL          :: ltest
         INTEGER          :: ier, l, m
         DOUBLE PRECISION :: rinv, phi_temp, dt_local, ti_temp, ne_temp,&
                             s_temp, x0, y0, z0, xw, yw, zw, te_temp, Zeff_temp
         DOUBLE PRECISION :: qf(3),qs(3),qe(3)
         DOUBLE PRECISION :: rlocal(num_depo), plocal(num_depo), zlocal(num_depo)
         DOUBLE PRECISION :: tilocal(num_depo), telocal(num_depo), nelocal(num_depo)
         DOUBLE PRECISION :: zefflocal(num_depo)
         DOUBLE PRECISION :: nilocal(NION,num_depo)
         DOUBLE PRECISION :: tau_inv(num_depo), energy(num_depo)
         DOUBLE PRECISION :: sigvii(num_depo), sigvcx(num_depo), sigvei(num_depo)
         ! For splines
         INTEGER :: i,j,k
         REAL*8 :: xparam, yparam, zparam
         INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
         REAL*8 :: fval(1)
         ! For Suzuki
         INTEGER :: A_IN(NION), Z_IN(NION)
         DOUBLE PRECISION :: ni_in(NION)

         !--------------------------------------------------------------
         !     Begin Subroutine
         !--------------------------------------------------------------
         ! Energy is needed in keV so 0.5*m*v*v/(ec*1000)

         ! This is the one that works for ADAS [kJ] E=0.5*m*v^2/1000
         ! Vll = V_neut (doesn't change durring neutral integration)
         ! energy in kJ
         energy = half*mymass*q(4)*q(4)*1D-3
         ! energy in keV (correct for Suzuki)
         !energy = half*mymass*q(4)*q(4)*1D-3/e_charge 
         
         qf(1) = q(1)*cos(q(2))
         qf(2) = q(1)*sin(q(2))
         qf(3) = q(3)

         !--------------------------------------------------------------
         !     Follow neutral into plasma using subgrid
         !--------------------------------------------------------------
         end_state(myline) = 3 ! It's a neutral for now
         xlast = qf(1)
         ylast = qf(2)
         zlast = qf(3)
         x0 = qf(1); y0 = qf(2); z0 = qf(3)
         DO l = 1, 3
            dt_local = stepsize(l)/q(4)
            DO
               qf = qf + myv_neut*dt_local
               q(1) = sqrt(qf(1)*qf(1)+qf(2)*qf(2))
               q(2) = ATAN2(qf(2),qf(1))
               q(3) = qf(3)
               t = t + dt_local
               phi_temp = MODULO(q(2), phimax)
               IF (phi_temp < 0) phi_temp = phi_temp + phimax
               !CALL EZspline_isInDomain(S_spl,q(1),phi_temp,q(3),ier)
               !IF (ier==0) THEN
               IF ((q(1) >= rmin-eps1) .and. (q(1) <= rmax+eps1) .and. &
                   (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
                   (q(3) >= zmin-eps3) .and. (q(3) <= zmax+eps3)) THEN
                  i = MIN(MAX(COUNT(raxis < q(1)),1),nr-1)
                  j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
                  k = MIN(MAX(COUNT(zaxis < q(3)),1),nz-1)
                  xparam = (q(1) - raxis(i)) * hri(i)
                  yparam = (phi_temp - phiaxis(j)) * hpi(j)
                  zparam = (q(3) - zaxis(k)) * hzi(k)
                  s_temp =1.5
                  !CALL EZspline_interp(S_spl,q(1),phi_temp,q(3),s_temp,ier)
                  CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                                  hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                                  S4D(1,1,1,1),nr,nphi,nz)
                  s_temp = fval(1)
                  IF (s_temp < one) EXIT
               END IF
               IF ((q(1) > 5*rmax)  .or. (q(1) < rmin)) THEN
                  t = my_end+dt_local
                  RETURN
               END IF  ! We're outside the grid
            END DO
            ! Take a step back
            qf = qf - myv_neut*dt_local
            t  =  t - dt_local
         END DO
         qs=qf

         !--------------------------------------------------------------
         !     Check to see if we hit the wall
         !--------------------------------------------------------------
         IF (lvessel_beam) THEN
            CALL collide(x0,y0,z0,qf(1),qf(2),qf(3),xw,yw,zw,ltest)
            IF (ltest) THEN
               q(1) = SQRT(qf(1)*qf(1)+qf(2)*qf(2))
               q(2) = ATAN2(qf(2),qf(1))
               q(3) = qf(3)
               end_state(myline) = 4
               CALL uncount_wall_hit
               RETURN
            END IF
         END IF

         xlast = qf(1)
         ylast = qf(2)
         zlast = qf(3)

         !--------------------------------------------------------------
         !     Follow particle track out of plasma
         !--------------------------------------------------------------
         OUTER: DO l = 1, 3
            dt_local = stepsize(l)/q(4)
            INNER: DO
               qf = qf + myv_neut*dt_local
               q(1) = sqrt(qf(1)*qf(1)+qf(2)*qf(2))
               q(2) = ATAN2(qf(2),qf(1))
               q(3) = qf(3)
               phi_temp = MODULO(q(2), phimax)
               IF (phi_temp < 0) phi_temp = phi_temp + phimax
               ! Assume we're in grid and only want to bug out if we're outside the grid
               IF ((q(1) >= rmin-eps1) .and. (q(1) <= rmax+eps1) .and. &
                   (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
                   (q(3) >= zmin-eps3) .and. (q(3) <= zmax+eps3)) THEN
                  i = MIN(MAX(COUNT(raxis < q(1)),1),nr-1)
                  j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
                  k = MIN(MAX(COUNT(zaxis < q(3)),1),nz-1)
                  xparam = (q(1) - raxis(i)) * hri(i)
                  yparam = (phi_temp - phiaxis(j)) * hpi(j)
                  zparam = (q(3) - zaxis(k)) * hzi(k)
                  CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                                  hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                                  S4D(1,1,1,1),nr,nphi,nz)
                  s_temp = fval(1)
                  IF (s_temp > one) EXIT INNER
               ELSE
                  EXIT INNER
               END IF
            END DO INNER
            ! Take a step back
            qf = qf - myv_neut*dt_local
         END DO OUTER
         qe=qf + myv_neut*dt_local


         !--------------------------------------------------------------
         !     Setup deposition arrays
         !--------------------------------------------------------------
         rlocal(1) = SQRT(qs(1)*qs(1)+qs(2)*qs(2))
         plocal(1) = ATAN2(qs(2),qs(1))
         zlocal(1) = qs(3)
         rlocal(num_depo) = SQRT(qe(1)*qe(1)+qe(2)*qe(2))
         plocal(num_depo) = ATAN2(qe(2),qe(1))
         zlocal(num_depo) = qe(3)
         DO i = 2, num_depo-1
            qf = (i-1)*(qe-qs)/(REAL(num_depo-1)) + qs
            rlocal(i) = sqrt(qf(1)*qf(1)+qf(2)*qf(2))
            plocal(i) = atan2(qf(2),qf(1))
            zlocal(i) = qf(3)
         END DO
         plocal = MODULO(plocal, phimax) ! Dont need to check for negative then
         ! Compute temp/density along path
         DO l = 1, num_depo
            i = MIN(MAX(COUNT(raxis < rlocal(l)),1),nr-1)
            j = MIN(MAX(COUNT(phiaxis < plocal(l)),1),nphi-1)
            k = MIN(MAX(COUNT(zaxis < zlocal(l)),1),nz-1)
            xparam = (rlocal(l) - raxis(i)) * hri(i)
            yparam = (plocal(l) - phiaxis(j)) * hpi(j)
            zparam = (zlocal(l) - zaxis(k)) * hzi(k)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            TI4D(1,1,1,1),nr,nphi,nz)
            tilocal(l) = MAX(fval(1),zero)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            TE4D(1,1,1,1),nr,nphi,nz)
            telocal(l) = MAX(fval(1),zero)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            NE4D(1,1,1,1),nr,nphi,nz)
            nelocal(l) = MAX(fval(1),zero)
            DO m = 1, NION
               CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            NI5D(1,1,1,1,m),nr,nphi,nz)
               nilocal(m,l) = MAX(fval(1),zero)
            END DO
         END DO
         tilocal = tilocal*1D-3
         telocal = telocal*1D-3
         tau_inv = zero
         
         IF (lsuzuki) THEN
            !--------------------------------------------------------------
            !     USE Suzuki to calcualte ionization rates
            !     Note: 10^18<ne<10^21
            !           E(keV/amu)/100 < Te < E(keV/amu)/2
            !--------------------------------------------------------------
            A_in = NINT(NI_AUX_M*inv_dalton)
            Z_in = NI_AUX_Z
            energy   = energy/(e_charge*mymass*inv_dalton) ! keV/amu
            DO l = 1, num_depo
               nelocal(l)  = MAX(MIN(nelocal(l),1E21),1E18)
               telocal(l)  = MAX(MIN(telocal(l),energy(l)*0.5),energy(l)*0.01)
               ni_in = nilocal(:,l)
               CALL suzuki_sigma(NION,energy(l),nelocal(l),telocal(l),ni_in,A_in,Z_in,tau_inv(l))
            END DO
            tau_inv = tau_inv*nelocal*ABS(q(4))*1E-4 !cm^2 to m^2 for sigma
         ELSE
            zefflocal = MATMUL(NI_AUX_Z*NI_AUX_Z,NILOCAL)/MATMUL(NI_AUX_Z,NILOCAL)
            zeff_temp = SUM(zefflocal)/DBLE(num_depo)
            !--------------------------------------------------------------
            !     USE ADAS to calcualte ionization rates
            !--------------------------------------------------------------
            ! Arguments to btsigv(irtype,beamchrg,eabeam,tatarg,n,izbeam,iztarg,btsigv, istat )hatom_btsigv.f90
            ! irtype 1:CX, 2:II
            ! beamchrg 1:neutral, 2:ion
            ! eabeam beam energy in keV/nucleon
            ! tatarg target energy in keV/nucleon
            ! n array size
            ! izbeam  Beam Z (always int)
            ! iztarg  Target Z (can be real)
            ! btsigv  cross section
            CALL adas_btsigv(2,1,energy,tilocal,num_depo,myZ,zeff_temp,sigvii,ier)  ! Ion Impact ionization cross-section term.
            CALL adas_btsigv(1,1,energy,tilocal,num_depo,myZ,zeff_temp,sigvcx,ier)  ! Charge Exchange ionization cross-section term.
            ! Arguments to sigvte(zneut,tevec,n1,sigv_adas,istat)
            ! zneut charge (=1)
            ! tevec electron temperature [keV]
            ! n1 array size
            ! This formula comes from the ADAS description of how to use the functions. (M. Gorelenkova)
            ! factor here is mp/me (assumes ion) Source NUBEAM: getsigs_adas.f, line 43
            telocal = telocal + to3*mpome*energy
            CALL adas_sigvte_ioniz(myZ,telocal,num_depo,sigvei,ier)        ! Electron Impact ionization cross-section term.
            ! Do this because Ztarg changes for each point.
            !DO l = 1, num_depo
            !   CALL adas_btsigv(2,1,energy,tilocal(l),1,myZ,zefflocal(l),sigvii(l),ier)  ! Ion Impact ionization cross-section term.
            !   CALL adas_btsigv(1,1,energy,tilocal(l),1,myZ,zefflocal(l),sigvcx(l),ier)  ! Charge Exchange ionization cross-section term.
            !END DO
            tau_inv = ((sigvii + sigvcx + sigvei)*nelocal) ! Delete a term if desired. (save a comment)
         END IF

         !--------------------------------------------------------------
         !     Calculate Ionization
         !--------------------------------------------------------------
         CALL RANDOM_NUMBER(rand_prob)
         cum_prob = one
         dt_local = SQRT(SUM((qe-qs)*(qe-qs)))/((num_depo-1)*q(4))
         tau_inv = EXP(-dt_local*tau_inv)
         DO l = 2, num_depo-1
            cum_prob = cum_prob*tau_inv(l)
            IF (cum_prob < rand_prob) EXIT
         END DO
         qf = qs + myv_neut*dt_local*(l-1)
         t  =  t + dt_local*(l-1)
         q(1) = SQRT(qf(1)*qf(1)+qf(2)*qf(2))
         q(2) = ATAN2(qf(2),qf(1))
         q(3) = qf(3)
         IF (l < num_depo-1) THEN
            IF ( (rlocal(l) <= rmin) .or. (rlocal(l) >= rmax) .or. &
                 (zlocal(l) <= zmin) .or. (zlocal(l) >= zmax) ) THEN 
               t = my_end + dt_local
               RETURN
            END IF
            i = MIN(MAX(COUNT(raxis < rlocal(l)),1),nr-1)
            j = MIN(MAX(COUNT(phiaxis < plocal(l)),1),nphi-1)
            k = MIN(MAX(COUNT(zaxis < zlocal(l)),1),nz-1)
            xparam = (rlocal(l) - raxis(i)) * hri(i)
            yparam = (plocal(l) - phiaxis(j)) * hpi(j)
            zparam = (zlocal(l) - zaxis(k)) * hzi(k)
            !CALL EZspline_interp(S_spl,rlocal(l),plocal(l),zlocal(l),s_temp,ier)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            S4D(1,1,1,1),nr,nphi,nz)
            s_temp = fval(1)
            lneut=.false.
            xlast = qf(1)
            ylast = qf(2)
            zlast = qf(3)
            RETURN
         END IF

         !--------------------------------------------------------------
         !     Follow neutral to wall or domain (big steps fine)
         !--------------------------------------------------------------
         x0 = qf(1); y0 = qf(2); z0 = qf(3)
         dt_local = 0.25/q(4)  
         ltest = .FALSE.
         xlast = qf(1)
         ylast = qf(2)
         zlast = qf(3)
         end_state(myline) = 3
         DO
            qf = qf + myv_neut*dt_local
            t = t + dt_local
            IF (lvessel_beam) CALL collide(x0,y0,z0,qf(1),qf(2),qf(3),xw,yw,zw,ltest)
            IF (ltest) THEN
               q(1) = SQRT(qf(1)*qf(1)+qf(2)*qf(2))
               q(2) = ATAN2(qf(2),qf(1))
               q(3) = qf(3)
               t = my_end+dt_local
               CALL uncount_wall_hit
               RETURN
            END IF
            !xlast = x0; ylast=y0; zlast=z0
            x0 = qf(1); y0 = qf(2); z0 = qf(3)
            q(1) = SQRT(qf(1)*qf(1)+qf(2)*qf(2))
            q(2) = ATAN2(qf(2),qf(1))
            q(3) = qf(3)
            IF ((q(1) > 2*rmax)  .or. (q(1) < rmin)) THEN; t = my_end+dt_local; RETURN; END IF  ! We're outside the grid
         END DO

         RETURN
      END SUBROUTINE beams3d_follow_neut

      !-----------------------------------------------------------------
      !     Function:      beams3d_ionize
      !     Authors:       S. Lazerson (lazerson@pppl.gov)
      !                    M. McMillan (matthew.mcmillan@my.wheaton.edu)
      !     Date:          12/01/2018
      !     Description:   Ionizes a particle relocating it to its
      !                    gyrocenter.
      !-----------------------------------------------------------------
      SUBROUTINE beams3d_ionize(t, q)
         USE beams3d_grid

         !--------------------------------------------------------------
         !     Input Parameters
         !          t          Location along fieldline in t
         !          q            (q(1),q(2),q(3),q(4)) = (R,phi,Z,vll)
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(inout) :: t
         DOUBLE PRECISION, INTENT(inout) :: q(4)

         !--------------------------------------------------------------
         !     Local parameters
         !--------------------------------------------------------------
         REAL*8, PARAMETER :: one = 1

         !--------------------------------------------------------------
         !     Local variables
         !--------------------------------------------------------------
         INTEGER          :: ier
         DOUBLE PRECISION :: r_temp, phi_temp, z_temp, x, y, &
                             br_temp, bp_temp, bz_temp, modb_temp, &
                             bx_temp, by_temp, binv, &
                             rho(3), rho2(3) 
         ! For splines
         INTEGER :: i,j,k
         REAL*8 :: xparam, yparam, zparam
         INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
         REAL*8 :: fval(1)

         !--------------------------------------------------------------
         !     Begin Subroutine
         !--------------------------------------------------------------
         lneut = .false.
         end_state(myline) = 0
         CALL RANDOM_NUMBER(rand_prob)

         ! Handle inputs
         ier = 0
         phi_temp = MOD(q(2),phimax)
         IF (phi_temp < 0) phi_temp = phi_temp + phimax
         r_temp = q(1)
         z_temp = q(3)

         ! Eval Spline
         i = MIN(MAX(COUNT(raxis < r_temp),1),nr-1)
         j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
         k = MIN(MAX(COUNT(zaxis < z_temp),1),nz-1)
         xparam = (r_temp - raxis(i)) * hri(i)
         yparam = (phi_temp - phiaxis(j)) * hpi(j)
         zparam = (z_temp - zaxis(k)) * hzi(k)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         BR4D(1,1,1,1),nr,nphi,nz)
         br_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         BPHI4D(1,1,1,1),nr,nphi,nz)
         bp_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         BZ4D(1,1,1,1),nr,nphi,nz)
         bz_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         MODB4D(1,1,1,1),nr,nphi,nz)
         modb_temp = fval(1)
         bx_temp = br_temp*cos(q(2))-bp_temp*sin(q(2))
         by_temp = br_temp*sin(q(2))+bp_temp*cos(q(2))
         binv = one/modb_temp

         ! rg = m*vperp/(q*B)
         ! Neutral position is average of gyrocenter |B|
         ! So step but then recalc at new position

         ! Calculate Gyroradius
         rho(1) = myv_neut(2)*bz_temp - myv_neut(3)*by_temp
         rho(2) = myv_neut(3)*bx_temp - myv_neut(1)*bz_temp
         rho(3) = myv_neut(1)*by_temp - myv_neut(2)*bx_temp
         !rho = rho*binv
         rho    = (mymass*binv*binv/mycharge)*rho 

         ! Calculate BxRg
         rho2(1) = by_temp*rho(3) - bz_temp*rho(2)
         rho2(2) = bz_temp*rho(1) - bx_temp*rho(3)
         rho2(3) = bx_temp*rho(2) - by_temp*rho(1)
         rho2 = rho2*binv

         ! Since rho==rg then rho2==rg but sqrt(rho+rho2)==rhog
         rho = inv_sqrt2 * rho * cos(pi2*rand_prob)
         rho2 = inv_sqrt2* rho2 * sin(pi2*rand_prob)

         ! Move to Gyrocenter
         x = q(1)*cos(q(2)) + rho(1) + rho2(1)
         y = q(1)*sin(q(2)) + rho(2) + rho2(2)
         z_temp= q(3) + rho(3) + rho2(3)
         r_temp = sqrt(x*x + y*y)
         phi_temp = ATAN2(y,x)
         IF (phi_temp<0) phi_temp = phi_temp + pi2

         ! Save on full toroidal grid
         q(1) = r_temp
         q(2) = phi_temp
         q(3) = z_temp

         ! Modify phi for splines
         phi_temp = MOD(phi_temp,phimax)

         ! Now recompute Splines
         i = MIN(MAX(COUNT(raxis < r_temp),1),nr-1)
         j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
         k = MIN(MAX(COUNT(zaxis < z_temp),1),nz-1)
         xparam = (r_temp - raxis(i)) * hri(i)
         yparam = (phi_temp - phiaxis(j)) * hpi(j)
         zparam = (z_temp - zaxis(k)) * hzi(k)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         BR4D(1,1,1,1),nr,nphi,nz)
         br_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         BPHI4D(1,1,1,1),nr,nphi,nz)
         bp_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         BZ4D(1,1,1,1),nr,nphi,nz)
         bz_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         MODB4D(1,1,1,1),nr,nphi,nz)
         modb_temp = fval(1)
         bx_temp = br_temp*cos(q(2))-bp_temp*sin(q(2))
         by_temp = br_temp*sin(q(2))+bp_temp*cos(q(2))
         binv = one/modb_temp


         ! Calculate the parallel velocity vll=v.B/B
         q(4) = binv*SUM(myv_neut*(/bx_temp,by_temp,bz_temp/))
         !q(4) = binv*( myv_neut(1)*bx_temp + myv_neut(2)*by_temp + myv_neut(3)*bz_temp )

         ! Calculate the magnetic moment mu = m*vperp^2/(2*B)=m*(v.v-vll.vll)/(2*B)
         moment = 0.5*binv*mymass*(SUM(myv_neut*myv_neut) - q(4)*q(4))
         !moment = 0.5*binv*mymass*(myv_neut(1)*myv_neut(1) + myv_neut(2)*myv_neut(2) + myv_neut(3)*myv_neut(3) - q(4)*q(4) )

         ! Check to see we didn't inject perfectly parallel (negative moment possible)
         IF (moment <= 0) THEN
            moment = 1000*TINY(moment)
            RETURN
         END IF

         RETURN

      END SUBROUTINE beams3d_ionize

      !-----------------------------------------------------------------
      !     Function:      beams3d_neutralize
      !     Authors:       M. McMillan (matthew.mcmillan@my.wheaton.edu)
      !     Date:          12/01/2018
      !     Description:   Calculate neutralization of an ion.  Note
      !                    that this never had the desired effect, 
      !                    is most deffinitely broken, and isn't called
      !                    anywhere in the code.
      !-----------------------------------------------------------------
      SUBROUTINE beams3d_neutralize(t, q)

         !--------------------------------------------------------------
         !     Input Parameters
         !          t          Location along fieldline in t
         !          q            (q(1),q(2),q(3),q(4)) = (R,phi,Z,vll)
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(inout) :: t
         DOUBLE PRECISION, INTENT(inout) :: q(4)

         !--------------------------------------------------------------
         !     Local parameters
         !--------------------------------------------------------------
         DOUBLE PRECISION, PARAMETER :: zero = 0
         DOUBLE PRECISION, PARAMETER :: one = 1

         !--------------------------------------------------------------
         !     Local Variables
         !--------------------------------------------------------------
         INTEGER          :: ier
         DOUBLE PRECISION :: r_temp, z_temp, phi_temp, modb_temp
         DOUBLE PRECISION :: bx_temp, by_temp, bz_temp, br_temp, bp_temp
         DOUBLE PRECISION :: binv, bnz, e1(3), e2(3), theta, rho(3), x, y, vperp

         !--------------------------------------------------------------
         !     Begin Subroutine
         !--------------------------------------------------------------
         ier = 0
         phi_temp = MOD(q(2),phimax)
         IF (phi_temp < zero) phi_temp = phi_temp + phimax
         r_temp = q(1)
         z_temp = q(3)

         ! Evaluate Splines
         CALL EZspline_interp(BR_spl,r_temp,phi_temp,z_temp,br_temp,ier)
         CALL EZspline_interp(BPHI_spl,r_temp,phi_temp,z_temp,bp_temp,ier)
         CALL EZspline_interp(BZ_spl,r_temp,phi_temp,z_temp,bz_temp,ier)
         CALL EZspline_interp(MODB_spl,r_temp,phi_temp,z_temp,modb_temp,ier)
         bx_temp = br_temp*cos(q(2))-bp_temp*sin(q(2))
         by_temp = br_temp*sin(q(2))+bp_temp*cos(q(2))
         binv = one/modb_temp
         bnz    = sqrt(modb_temp*modb_temp - bz_temp*bz_temp)
         vperp  = sqrt( 2*modb_temp*moment/mymass )

         ! Unit vectors
         IF (bnz == zero) THEN
            e1    = (/one,zero,zero/)
            e2    = (/zero,one,zero/)
         ELSE
            e1 = (/ -bz_temp*bx_temp, -bz_temp*by_temp, one /)
            e1 = e1*binv/bnz
            e2 = (/ by_temp*e1(3)-bz_temp*e1(2), bz_temp*e1(1)-bx_temp*e1(3), bx_temp*e1(2)-by_temp*e1(1) /)
            e2 = e2*binv
         END IF

         CALL RANDOM_NUMBER(theta)
         theta = theta*pi2
         rho         = ( sin(theta)*e1 + cos(theta)*e2 )

         myv_neut(1) = q(4)*binv*bx_temp - vperp*binv*( by_temp*rho(3)-bz_temp*rho(2) )
         myv_neut(2) = q(4)*binv*by_temp - vperp*binv*( bz_temp*rho(1)-bx_temp*rho(3) )
         myv_neut(3) = q(4)*binv*bz_temp - vperp*binv*( bx_temp*rho(2)-by_temp*rho(1) )

         rho         = ( mymass*vperp*binv/mycharge )*rho
         x           = q(1)*cos(q(2)) + rho(1)
         y           = q(1)*sin(q(2)) + rho(2)

         q(1)        = sqrt(x*x + y*y)
         q(2)        = ATAN2(y,x)
         IF (q(2) < 0) q(2) = q(2)+pi2
         q(3)        = q(3) + rho(3)

         moment = 0.5*binv*mymass*( myv_neut(1)*myv_neut(1) + myv_neut(2)*myv_neut(2) + myv_neut(3)*myv_neut(3) - q(4)*q(4) )

         lneut = .true.

         CALL RANDOM_NUMBER(rand_prob)

      END SUBROUTINE beams3d_neutralize

      !-----------------------------------------------------------------
      !     Function:      beams3d_DTRATE
      !     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
      !     Date:          09/30/2020
      !     Description:   Calculates the D-T Reaction rate, assumes
      !                    50/50 n_D/n_T based on n_e. See:
      !                    H.-S. Bosch and G. M. Hale 1992 Nucl. Fusion 32 611
      !                    https://doi.org/10.1088/0029-5515/32/4/I07
      !-----------------------------------------------------------------
      SUBROUTINE beams3d_DTRATE(q,reactrate)
         !--------------------------------------------------------------
         !     Input Parameters
         !          q            (q(1),q(2),q(3)) = (R,phi,Z)
         !          reactrate    Reaction rate (part/(m^3*s))
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(inout) :: q(3)
         DOUBLE PRECISION, INTENT(out) :: reactrate

         !--------------------------------------------------------------
         !     Local Variables
         !        r_temp     Helpers (r,phi,z, ne, ti, ze)
         !        zeta       Fusion helper
         !        theta      Fusion helper
         !        eta        Fusion Helper
         !        i,j,k      Spline Grid indicies
         !        xparam     Spline subgrid factor [0,1] (yparam,zparam)
         !        ict        Spline output control
         !        fval       Spline output array
         !--------------------------------------------------------------
         DOUBLE PRECISION :: r_temp, z_temp, phi_temp, &
                             ti_temp, zeta, theta, eta, &
                             nd_temp, nt_temp
         ! For splines
         INTEGER :: i,j,k
         REAL*8 :: xparam, yparam, zparam
         INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
         REAL*8 :: fval(1)

         !--------------------------------------------------------------
         !     Local Parameters
         !--------------------------------------------------------------
         INTEGER, PARAMETER :: mrc2 = 1124656
         DOUBLE PRECISION, PARAMETER :: BG = 34.3827
         DOUBLE PRECISION, DIMENSION(7), PARAMETER :: &
                CARR = (/ 1.17302E-09,  1.51361E-02,  7.51886E-02, &
                          4.60643E-03,  1.35000E-02, -1.06750E-04, &
                          1.36600E-05/)

         !--------------------------------------------------------------
         !     Begin Subroutine
         !--------------------------------------------------------------

         ! Setup position in a vll arrays
         r_temp   = q(1)
         phi_temp = MODULO(q(2), phimax)
         IF (phi_temp < 0) phi_temp = phi_temp + phimax
         z_temp   = q(3)

         ! Initialize values
         ti_temp = 0; nd_temp = 0; nt_temp = 0; reactrate = 0

         ! Check that we're inside the domain then proceed
         IF ((r_temp >= rmin-eps1) .and. (r_temp <= rmax+eps1) .and. &
             (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
             (z_temp >= zmin-eps3) .and. (z_temp <= zmax+eps3)) THEN
            i = MIN(MAX(COUNT(raxis < r_temp),1),nr-1)
            j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
            k = MIN(MAX(COUNT(zaxis < z_temp),1),nz-1)
            xparam = (r_temp - raxis(i)) * hri(i)
            yparam = (phi_temp - phiaxis(j)) * hpi(j)
            zparam = (z_temp - zaxis(k)) * hzi(k)
            ! Evaluate the Splines
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            TI4D(1,1,1,1),nr,nphi,nz)
            ti_temp = max(fval(1),zero)
            ! Assume 1 and 2 are D and T
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            NI5D(1,1,1,1,1),nr,nphi,nz)
            nd_temp = max(fval(1),zero)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            NI5D(1,1,1,1,2),nr,nphi,nz)
            nt_temp = max(fval(1),zero)
         ELSE
            RETURN
         END IF
         IF (ti_temp <= zero) RETURN ! Get out if ti small
         ti_temp = ti_temp*1E-3 ! to keV
         zeta =  one - ((((CARR(6)*ti_temp)+CARR(4))*ti_temp+CARR(2))*ti_temp)/ &
                       ((((CARR(7)*ti_temp)+CARR(5))*ti_temp+CARR(3))*ti_temp+one)
         theta = ti_temp/zeta
         eta   = (BG*BG/(4*theta))**(one/3.0)

         reactrate = 1E-6*CARR(1)*theta*SQRT(eta/(mrc2*ti_temp*ti_temp*ti_temp))*EXP(-3*eta)

         reactrate = reactrate*nd_temp*nt_temp
         RETURN

      END SUBROUTINE beams3d_DTRATE

      !-----------------------------------------------------------------
      !     Function:      beams3d_DDTRATE
      !     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
      !     Date:          09/30/2020
      !     Description:   Calculates the D-D->T Reaction rate, assumes
      !                    50/50 n_D/n_T based on n_e. See:
      !                    H.-S. Bosch and G. M. Hale 1992 Nucl. Fusion 32 611
      !                    https://doi.org/10.1088/0029-5515/32/4/I07
      !-----------------------------------------------------------------
      SUBROUTINE beams3d_DDTRATE(q,reactrate)
         !--------------------------------------------------------------
         !     Input Parameters
         !          q            (q(1),q(2),q(3)) = (R,phi,Z)
         !          reactrate    Reaction rate (part/(m^3*s))
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(inout) :: q(3)
         DOUBLE PRECISION, INTENT(out) :: reactrate

         !--------------------------------------------------------------
         !     Local Variables
         !        r_temp     Helpers (r,phi,z, ne, ti, ze)
         !        zeta       Fusion helper
         !        theta      Fusion helper
         !        eta        Fusion Helper
         !        i,j,k      Spline Grid indicies
         !        xparam     Spline subgrid factor [0,1] (yparam,zparam)
         !        ict        Spline output control
         !        fval       Spline output array
         !--------------------------------------------------------------
         DOUBLE PRECISION :: r_temp, z_temp, phi_temp, &
                             ti_temp, zeta, theta, eta, &
                             nd_temp
         ! For splines
         INTEGER :: i,j,k
         REAL*8 :: xparam, yparam, zparam
         INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
         REAL*8 :: fval(1)

         !--------------------------------------------------------------
         !     Local Parameters
         !--------------------------------------------------------------
         INTEGER, PARAMETER :: mrc2 = 937814
         DOUBLE PRECISION, PARAMETER :: BG = 31.3970
         DOUBLE PRECISION, DIMENSION(7), PARAMETER :: &
                CARR = (/ 5.65718E-12,  3.41267E-03,  1.99167E-03, &
                          0.00000E+00,  1.05060E-05,  0.00000E+00, &
                          0.00000E-00/)

         !--------------------------------------------------------------
         !     Begin Subroutine
         !--------------------------------------------------------------

         ! Setup position in a vll arrays
         r_temp   = q(1)
         phi_temp = MODULO(q(2), phimax)
         IF (phi_temp < 0) phi_temp = phi_temp + phimax
         z_temp   = q(3)

         ! Initialize values
         ti_temp = 0; nd_temp = 0; reactrate = 0

         ! Check that we're inside the domain then proceed
         IF ((r_temp >= rmin-eps1) .and. (r_temp <= rmax+eps1) .and. &
             (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
             (z_temp >= zmin-eps3) .and. (z_temp <= zmax+eps3)) THEN
            i = MIN(MAX(COUNT(raxis < r_temp),1),nr-1)
            j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
            k = MIN(MAX(COUNT(zaxis < z_temp),1),nz-1)
            xparam = (r_temp - raxis(i)) * hri(i)
            yparam = (phi_temp - phiaxis(j)) * hpi(j)
            zparam = (z_temp - zaxis(k)) * hzi(k)
            ! Evaluate the Splines
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            TI4D(1,1,1,1),nr,nphi,nz)
            ti_temp = max(fval(1),zero)
            ! Assume 1 and 2 are D and T
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            NI5D(1,1,1,1,1),nr,nphi,nz)
            nd_temp = max(fval(1),zero)
         ELSE
            RETURN
         END IF
         IF (ti_temp <= zero) RETURN ! Get out if ti small
         ti_temp = ti_temp*1E-3 ! to keV
         zeta =  one - ((((CARR(6)*ti_temp)+CARR(4))*ti_temp+CARR(2))*ti_temp)/ &
                       ((((CARR(7)*ti_temp)+CARR(5))*ti_temp+CARR(3))*ti_temp+one)
         theta = ti_temp/zeta
         eta   = (BG*BG/(4*theta))**(one/3.0)

         reactrate = 1E-6*CARR(1)*theta*SQRT(eta/(mrc2*ti_temp*ti_temp*ti_temp))*EXP(-3*eta)

         reactrate = reactrate*nd_temp*nd_temp*0.5

         RETURN

      END SUBROUTINE beams3d_DDTRATE

      !-----------------------------------------------------------------
      !     Function:      beams3d_DDHe3RATE
      !     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
      !     Date:          09/30/2020
      !     Description:   Calculates the D-D->He3 Reaction rate, assumes
      !                    50/50 n_D/n_T based on n_e. See:
      !                    H.-S. Bosch and G. M. Hale 1992 Nucl. Fusion 32 611
      !                    https://doi.org/10.1088/0029-5515/32/4/I07
      !-----------------------------------------------------------------
      SUBROUTINE beams3d_DDHe3RATE(q,reactrate)
         !--------------------------------------------------------------
         !     Input Parameters
         !          q            (q(1),q(2),q(3)) = (R,phi,Z)
         !          reactrate    Reaction rate (part/(m^3*s))
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(inout) :: q(3)
         DOUBLE PRECISION, INTENT(out) :: reactrate

         !--------------------------------------------------------------
         !     Local Variables
         !        r_temp     Helpers (r,phi,z, ne, ti, ze)
         !        zeta       Fusion helper
         !        theta      Fusion helper
         !        eta        Fusion Helper
         !        i,j,k      Spline Grid indicies
         !        xparam     Spline subgrid factor [0,1] (yparam,zparam)
         !        ict        Spline output control
         !        fval       Spline output array
         !--------------------------------------------------------------
         DOUBLE PRECISION :: r_temp, z_temp, phi_temp, nd_temp, &
                             ti_temp, zeta, theta, eta
         ! For splines
         INTEGER :: i,j,k
         REAL*8 :: xparam, yparam, zparam
         INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
         REAL*8 :: fval(1)

         !--------------------------------------------------------------
         !     Local Parameters
         !--------------------------------------------------------------
         INTEGER, PARAMETER :: mrc2 = 937814
         DOUBLE PRECISION, PARAMETER :: BG = 31.3970
         DOUBLE PRECISION, DIMENSION(7), PARAMETER :: &
                CARR = (/ 5.43360E-12,  5.85778E-03,  7.68222E-03, &
                          0.00000E+00, -2.96400E-06,  0.00000E+00, &
                          0.00000E+00/)

         !--------------------------------------------------------------
         !     Begin Subroutine
         !--------------------------------------------------------------

         ! Setup position in a vll arrays
         r_temp   = q(1)
         phi_temp = MODULO(q(2), phimax)
         IF (phi_temp < 0) phi_temp = phi_temp + phimax
         z_temp   = q(3)

         ! Initialize values
         ti_temp = 0; nd_temp = 0; reactrate = 0

         ! Check that we're inside the domain then proceed
         IF ((r_temp >= rmin-eps1) .and. (r_temp <= rmax+eps1) .and. &
             (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
             (z_temp >= zmin-eps3) .and. (z_temp <= zmax+eps3)) THEN
            i = MIN(MAX(COUNT(raxis < r_temp),1),nr-1)
            j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
            k = MIN(MAX(COUNT(zaxis < z_temp),1),nz-1)
            xparam = (r_temp - raxis(i)) * hri(i)
            yparam = (phi_temp - phiaxis(j)) * hpi(j)
            zparam = (z_temp - zaxis(k)) * hzi(k)
            ! Evaluate the Splines
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            TI4D(1,1,1,1),nr,nphi,nz)
            ti_temp = max(fval(1),zero)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            NI5D(1,1,1,1,1),nr,nphi,nz)
            nd_temp = max(fval(1),zero)
         ELSE
            RETURN
         END IF
         IF (ti_temp <= zero) RETURN ! Get out if ti small
         ti_temp = ti_temp*1E-3 ! to keV
         zeta =  one - ((((CARR(6)*ti_temp)+CARR(4))*ti_temp+CARR(2))*ti_temp)/ &
                       ((((CARR(7)*ti_temp)+CARR(5))*ti_temp+CARR(3))*ti_temp+one)
         theta = ti_temp/zeta
         eta   = (BG*BG/(4*theta))**(one/3.0)

         reactrate = 1E-6*CARR(1)*theta*SQRT(eta/(mrc2*ti_temp*ti_temp*ti_temp))*EXP(-3*eta)

         reactrate = reactrate*nd_temp*nd_temp*0.5
         RETURN

      END SUBROUTINE beams3d_DDHe3RATE

      !-----------------------------------------------------------------
      !     Function:      beams3d_MODB
      !     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
      !     Date:          09/30/2020
      !     Description:   Returns |B| at a point in space
      !-----------------------------------------------------------------
      SUBROUTINE beams3d_MODB(q,B)
         !--------------------------------------------------------------
         !     Input Parameters
         !          q            (q(1),q(2),q(3)) = (R,phi,Z)
         !          reactrate    Reaction rate (part/(m^3*s))
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(inout) :: q(3)
         DOUBLE PRECISION, INTENT(out) :: B

         !--------------------------------------------------------------
         !     Local Variables
         !        r_temp     Helpers (r,phi,z)
         !        i,j,k      Spline Grid indicies
         !        xparam     Spline subgrid factor [0,1] (yparam,zparam)
         !        ict        Spline output control
         !        fval       Spline output array
         !--------------------------------------------------------------
         DOUBLE PRECISION :: r_temp, z_temp, phi_temp
         ! For splines
         INTEGER :: i,j,k
         REAL*8 :: xparam, yparam, zparam
         INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
         REAL*8 :: fval(1)

         !--------------------------------------------------------------
         !     Begin Subroutine
         !--------------------------------------------------------------

         ! Setup position in a vll arrays
         r_temp   = q(1)
         phi_temp = MODULO(q(2), phimax)
         IF (phi_temp < 0) phi_temp = phi_temp + phimax
         z_temp   = q(3)

         ! Initialize values
         B = zero

         ! Check that we're inside the domain then proceed
         IF ((r_temp >= rmin-eps1) .and. (r_temp <= rmax+eps1) .and. &
             (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
             (z_temp >= zmin-eps3) .and. (z_temp <= zmax+eps3)) THEN
            i = MIN(MAX(COUNT(raxis < r_temp),1),nr-1)
            j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
            k = MIN(MAX(COUNT(zaxis < z_temp),1),nz-1)
            xparam = (r_temp - raxis(i)) * hri(i)
            yparam = (phi_temp - phiaxis(j)) * hpi(j)
            zparam = (z_temp - zaxis(k)) * hzi(k)
            ! Evaluate the Splines
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            MODB4D(1,1,1,1),nr,nphi,nz)
            B = max(fval(1),zero)
         ELSE
            RETURN
         END IF

         RETURN

      END SUBROUTINE beams3d_MODB

      !-----------------------------------------------------------------
      !     Function:      beams3d_SFLX
      !     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
      !     Date:          09/30/2020
      !     Description:   Returns normalized toroidal flux at a point
      !                    in space.
      !-----------------------------------------------------------------
      SUBROUTINE beams3d_SFLX(q,S)
         !--------------------------------------------------------------
         !     Input Parameters
         !          q    (q(1),q(2),q(3)) = (R,phi,Z)
         !          S    Backbround grid flux
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(inout) :: q(3)
         DOUBLE PRECISION, INTENT(out) :: S

         !--------------------------------------------------------------
         !     Local Variables
         !        r_temp     Helpers (r,phi,z)
         !        i,j,k      Spline Grid indicies
         !        xparam     Spline subgrid factor [0,1] (yparam,zparam)
         !        ict        Spline output control
         !        fval       Spline output array
         !--------------------------------------------------------------
         DOUBLE PRECISION :: r_temp, z_temp, phi_temp
         ! For splines
         INTEGER :: i,j,k
         REAL*8 :: xparam, yparam, zparam
         INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
         REAL*8 :: fval(1)

         !--------------------------------------------------------------
         !     Begin Subroutine
         !--------------------------------------------------------------

         ! Setup position in a vll arrays
         r_temp   = q(1)
         phi_temp = MODULO(q(2), phimax)
         IF (phi_temp < 0) phi_temp = phi_temp + phimax
         z_temp   = q(3)

         ! Initialize values
         S = 2

         ! Check that we're inside the domain then proceed
         IF ((r_temp >= rmin-eps1) .and. (r_temp <= rmax+eps1) .and. &
             (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
             (z_temp >= zmin-eps3) .and. (z_temp <= zmax+eps3)) THEN
            i = MIN(MAX(COUNT(raxis < r_temp),1),nr-1)
            j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
            k = MIN(MAX(COUNT(zaxis < z_temp),1),nz-1)
            xparam = (r_temp - raxis(i)) * hri(i)
            yparam = (phi_temp - phiaxis(j)) * hpi(j)
            zparam = (z_temp - zaxis(k)) * hzi(k)
            ! Evaluate the Splines
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            S4D(1,1,1,1),nr,nphi,nz)
            S = max(fval(1),zero)
         END IF

         RETURN

      END SUBROUTINE beams3d_SFLX

      !-----------------------------------------------------------------
      !     Function:      beams3d_suv2rzp
      !     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
      !     Date:          03/07/2021
      !     Description:   Returns R and Z given s,u,v coordinate
      !-----------------------------------------------------------------
      SUBROUTINE beams3d_suv2rzp(s,u,v,r_out,z_out,phi_out)
         !--------------------------------------------------------------
         !     Input Parameters
         !          s            Normalized Toroidal Flux
         !          u            Poloidal Angle
         !          v            Toroidal Angle
         !          r_out        Major Radius Distance R
         !          z_out        Vertical Distance Z
         !          phi_out      Toroidal Angle
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(inout) :: s
         DOUBLE PRECISION, INTENT(inout) :: u
         DOUBLE PRECISION, INTENT(inout) :: v
         DOUBLE PRECISION, INTENT(inout) :: r_out
         DOUBLE PRECISION, INTENT(inout) :: z_out
         DOUBLE PRECISION, INTENT(out) :: phi_out
         !REAL(rprec), POINTER, DIMENSION(:,:,:,:), INTENT(inout) :: X4D, Y4D

         !--------------------------------------------------------------
         !     Local Variables
         !        residual   Residual Error
         !--------------------------------------------------------------
         INTEGER          :: n
         DOUBLE PRECISION :: s0, u0, residual, detJ, delR, delZ, fnorm, &
                             factor, x, y, x0, y0, x_term, y_term, dxdR, dxdZ, dydR, dydZ

         ! For splines
         INTEGER :: i,j,k, ier
         REAL*8 :: xparam, yparam, zparam
         INTEGER, parameter :: ict(8)=(/1,1,1,1,0,0,0,0/)
         REAL*8 :: fvalx(1,4),fvaly(1,4) !(f,df/fR,df/dphi,dfdZ)


         !--------------------------------------------------------------
         !     Begin Subroutine
         !--------------------------------------------------------------

         !Begin Newton Method
         residual = 1.0
         factor = 1.0
         IF (r_out<0) r_out = raxis(1)+(raxis(nr)-raxis(1))*.75
         
         ! PHI does not change
         phi_out = MOD(v,MAXVAL(phiaxis))
         j = MIN(MAX(COUNT(phiaxis < phi_out),1),nphi-1)
         yparam = (phi_out - phiaxis(j)) * hpi(j)

         ! Adjust u
         u = MOD(u,pi2)

         x0 = s * COS(u)
         y0 = s * SIN(U)

         fnorm = x0*x0+y0*y0
         fnorm = MIN(1./fnorm,1E5)
         n = 1

         ! Loop Basically a NEWTON's METHOD
         !  F(R,Z) = (s0-s(R,Z))*(u0-u(R,Z))
         !  dF/dR  = -ds/dR*(u0-u(R,Z))-du/dR*(s0-s(R,Z))
         !  dF/dZ  = -ds/dZ*(u0-u(R,Z))-du/dZ*(s0-s(R,Z))
         DO WHILE (residual > 1.0E-9 .and. n<1000)
            i = MIN(MAX(COUNT(raxis < r_out),1),nr-1)
            k = MIN(MAX(COUNT(zaxis < z_out),1),nz-1)
            xparam = (r_out - raxis(i)) * hri(i)
            zparam = (z_out - zaxis(k)) * hzi(k)
            ! Evaluate the Splines
            CALL R8HERM3FCN(ict,1,1,fvalx,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            X4D(1,1,1,1),nr,nphi,nz)
            CALL R8HERM3FCN(ict,1,1,fvaly,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            Y4D(1,1,1,1),nr,nphi,nz)

            x_term   = x0 - fvalx(1,1)
            y_term   = y0 - fvaly(1,1)
            
            detJ = fvalx(1,2) * fvaly(1,4) - fvaly(1,2) * fvalx(1,4)
            detJ = MAX(detJ,0.0001) !Upper bound for step size as detJ enters in denominator
            delR = -(-x_term*fvaly(1,4) + y_term*fvalx(1,4))/detJ
            delZ = -( x_term*fvaly(1,2)  - y_term*fvalx(1,2))/detJ

            delR = MIN(MAX(delR,-hr(1)),hr(1))
            delZ = MIN(MAX(delZ,-hz(1)),hz(1))

            residual = (x_term*x_term+y_term*y_term)*fnorm
            !WRITE(6,*) '----- ',s,u,s0,u0,r_out,z_out,residual,tau,delR,delZ

            IF (residual < 0.01) THEN !"Damping" of oscillation
               delR = delR*0.5
               delZ = delZ*0.5
            END IF

            r_out = MAX(MIN(r_out + delR*factor,raxis(nr)),raxis(1))
            z_out = MAX(MIN(z_out + delZ*factor,zaxis(nz)),zaxis(1))
            !WRITE(6,*) '----- ',s,u,s0,u0,r_out,z_out,residual
            n=n+1
         END DO

 !        IF (n>=500) PRINT *,s,u,fnorm,residual


         RETURN

      END SUBROUTINE beams3d_suv2rzp

      !-----------------------------------------------------------------
      !     Function:      beams3d_fbounce
      !     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
      !     Date:          10/24/2020
      !     Description:   Calculates pseudo bounce frequency.
      !                    Uses large aspect ratio cicular cross
      !                    sction formulation.
      !-----------------------------------------------------------------
      SUBROUTINE beams3d_fbounce(q,mu,mass,fbounce)
         !--------------------------------------------------------------
         !     Input Parameters
         !          q            (q(1),q(2),q(3)) = (R,phi,Z)
         !          mu           Magnetic moment (J/T)
         !          mass         Particle mass (kg)
         !          fbounce      Reaction rate (part/(m^3*s))
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(inout) :: q(3)
         DOUBLE PRECISION, INTENT(inout) :: mu
         DOUBLE PRECISION, INTENT(inout) :: mass
         DOUBLE PRECISION, INTENT(out) :: fbounce

         !--------------------------------------------------------------
         !     Local Variables
         !        r_temp     Helpers (r,phi,z)
         !        i,j,k      Spline Grid indicies
         !        xparam     Spline subgrid factor [0,1] (yparam,zparam)
         !        ict        Spline output control
         !        fval       Spline output array
         !--------------------------------------------------------------
         DOUBLE PRECISION :: r_temp, z_temp, phi_temp, vperp, R0, Z0, r
         ! For splines
         INTEGER :: i,j,k
         REAL*8 :: xparam, yparam, zparam
         INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
         REAL*8 :: fval(1)

         !--------------------------------------------------------------
         !     Begin Subroutine
         !--------------------------------------------------------------

         ! Setup position in a vll arrays
         r_temp   = q(1)
         phi_temp = MODULO(q(2), phimax)
         IF (phi_temp < 0) phi_temp = phi_temp + phimax
         z_temp   = q(3)

         ! Initialize values
         fbounce = zero
         R0      = (raxis(nr)-raxis(1))*half+raxis(1)
         Z0      = (zaxis(nr)-zaxis(1))*half+zaxis(1)
         r       = SQRT((r_temp-R0)*(r_temp-R0)+(z_temp-Z0)*(z_temp-Z0))

         ! Check that we're inside the domain then proceed
         IF ((r_temp >= rmin-eps1) .and. (r_temp <= rmax+eps1) .and. &
             (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
             (z_temp >= zmin-eps3) .and. (z_temp <= zmax+eps3)) THEN
            i = MIN(MAX(COUNT(raxis < r_temp),1),nr-1)
            j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
            k = MIN(MAX(COUNT(zaxis < z_temp),1),nz-1)
            xparam = (r_temp - raxis(i)) * hri(i)
            yparam = (phi_temp - phiaxis(j)) * hpi(j)
            zparam = (z_temp - zaxis(k)) * hzi(k)
            ! Evaluate the Splines
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            MODB4D(1,1,1,1),nr,nphi,nz)
            !vperp = SQRT(2*mu*fval(1)/mass)
            !fbounce = vperp*SQRT(half*r/R0)/(R0*pi2) ! assume q=1
            fbounce = SQRT(r*mu*fval(1)/(mass*R0*R0*R0*pi2*pi2))
         ELSE
            RETURN
         END IF

         RETURN

      END SUBROUTINE beams3d_fbounce

      !-----------------------------------------------------------------
      !     Function:      beams3d_calc_dt
      !     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
      !     Date:          10/24/2020
      !     Description:   Calculates the timestep using the bounce
      !                    frequency from beasm3d_fbounce.
      !-----------------------------------------------------------------
      SUBROUTINE beams3d_calc_dt(q,mu,mass,dt)
         !--------------------------------------------------------------
         !     Input Parameters
         !          q            (q(1),q(2),q(3),q(4)) = (R,phi,Z,vll)
         !          mu           Magnetic moment (J/T)
         !          fbounce      Reaction rate (part/(m^3*s))
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(inout) :: q(4)
         DOUBLE PRECISION, INTENT(inout) :: mu
         DOUBLE PRECISION, INTENT(inout) :: mass
         DOUBLE PRECISION, INTENT(out) :: dt

         !--------------------------------------------------------------
         !     Local Variables
         !        freq_bounce Approx bounce frequency
         !--------------------------------------------------------------
         DOUBLE PRECISION :: freq_bounce, tf_max

         !--------------------------------------------------------------
         !     Begin Subroutine
         !--------------------------------------------------------------

         ! Define max time to follow particle
         tf_max = my_end

         ! Get bounce frequncy
         !CALL beams3d_fbounce(q(1:3),mu,mass,freq_bounce)

         ! Timestep is a fraction of bounce frequency
         !dt = one/(64*freq_bounce)
         !dt = SIGN(MAX(dt,1D-9),tf_max) ! Limiter and sign

         ! Use max distance
         dt = lendt_m/q(4)
         dt = SIGN(MAX(dt,1D-9),tf_max) ! Limiter and sign

         ! Make subtimestep fit (min 2 due to logic)
         ndt_max = MAX(CEILING(tf_max/(dt*NPOINC)),2)
         dt = tf_max/(ndt_max*NPOINC)

         RETURN

      END SUBROUTINE beams3d_calc_dt

END MODULE beams3d_physics_mod
