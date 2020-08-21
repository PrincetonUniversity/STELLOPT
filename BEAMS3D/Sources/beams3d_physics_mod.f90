!-----------------------------------------------------------------------
!     Module:        beams3d_physics_mod
!     Authors:       M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          07/20/2012
!     Description:   This module contains the physics subroutines for
!                    manipulating particles and finding when to do so.
!-----------------------------------------------------------------------
MODULE beams3d_physics_mod

      !-----------------------------------------------------------------------
      !     Libraries
      !-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime, ONLY: lneut, pi, pi2, dt, lverb, ADAS_ERR, &
                                 dt_save, lbbnbi, weight
      USE beams3d_lines, ONLY: R_lines, Z_lines, PHI_lines, &
                               myline, mytdex, moment, ltherm, &
                               nsteps, nparticles, vll_lines, &
                               moment_lines, mybeam, mycharge, myZ, &
                               mymass, myv_neut, B_temp, rand_prob, &
                               cum_prob, tau, &
                               epower_prof, ipower_prof, &
                               end_state, fact_crit, fact_pa, fact_vsound, &
                               ns_prof1, ns_prof2, ns_prof3, ns_prof4, &
                               ns_prof5
      USE beams3d_grid, ONLY: BR_spl, BZ_spl, delta_t, BPHI_spl, MODB_spl, MODB4D, &
                              phimax, S4D, TE4D, NE4D, TI4D, ZEFF4D, &
                              nr, nphi, nz, rmax, rmin, zmax, zmin, &
                              phimin, eps1, eps2, eps3, raxis, phiaxis, zaxis
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
      DOUBLE PRECISION, PRIVATE, PARAMETER :: inv_sqrt2     = 0.7071067812   !pi^(1/2)
      DOUBLE PRECISION, PRIVATE, PARAMETER :: mpome         = 5.44602984424355D-4 !e_c
      DOUBLE PRECISION, PRIVATE, PARAMETER :: inv_dalton    = 6.02214076208E+26 ! 1./AMU [1/kg]
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
                          inv_mymass, speed_cube, vcrit_cube, vfrac, modb, s_temp, vc3_tauinv
         DOUBLE PRECISION :: Ebench  ! for ASCOT Benchmark
         ! For splines
         INTEGER :: i,j,k, l
         REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
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
         te_temp  = 0; ne_temp  = 0; ti_temp  = 0; speed = 0; reduction = 0

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
            hx     = raxis(i+1) - raxis(i)
            hy     = phiaxis(j+1) - phiaxis(j)
            hz     = zaxis(k+1) - zaxis(k)
            hxi    = one / hx
            hyi    = one / hy
            hzi    = one / hz
            xparam = (r_temp - raxis(i)) * hxi
            yparam = (phi_temp - phiaxis(j)) * hyi
            zparam = (z_temp - zaxis(k)) * hzi
            !CALL R8HERM3xyz(r_temp,phi_temp,z_temp,&
            !                MODB_spl%x1(1),MODB_spl%n1,&
            !                MODB_spl%x2(1),MODB_spl%n2,&
            !                MODB_spl%x3(1),MODB_spl%n3,&
            !                MODB_spl%ilin1,MODB_spl%ilin2,MODB_spl%ilin3,&
            !                i,j,k,xparam,yparam,zparam,&
            !                hx,hxi,hy,hyi,hz,hzi,ier)
            ! Evaluate the Splines
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            MODB4D(1,1,1,1),nr,nphi,nz)
            modb = fval(1)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            TE4D(1,1,1,1),nr,nphi,nz)
            te_temp = max(fval(1),zero)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            NE4D(1,1,1,1),nr,nphi,nz)
            ne_temp = max(fval(1),zero)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            TI4D(1,1,1,1),nr,nphi,nz)
            ti_temp = max(fval(1),zero)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            S4D(1,1,1,1),nr,nphi,nz)
            s_temp = fval(1)

            ! Helpers
            te_cube = te_temp * te_temp * te_temp
            inv_mymass = 1/mymass

            !-----------------------------------------------------------
            !  Calculate Coulomb Logarithm (NRL pg. 35)
            !     te in eV and ne in cm^-3
            !-----------------------------------------------------------
            IF ((te_temp > 0).and.(ne_temp > 0)) THEN
               IF (te_temp < 10*myZ*myZ) THEN
                  coulomb_log = 23 - log( myZ*sqrt(ne_temp*1E-6/(te_cube) )   )
               ELSE
                  coulomb_log = 24 - log( sqrt(ne_temp*1E-6)/(te_temp) )
               END IF
               IF (coulomb_log .le. 1) coulomb_log = 1
               ! Callen Ch2 pg41 eq2.135 (fact*Vtherm; Vtherm = SQRT(2*E/mass) so E in J not eV)
               !v_crit = fact_crit*SQRT(2*te_temp*inv_mymass*e_charge)
               v_crit = fact_crit*SQRT(te_temp)
               !v_crit = (( 0.75*sqrt_pi*electron_mass*inv_mymass )**0.33333333333 )*sqrt(te_temp)*5.93096892024D5
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
            v_s = fact_vsound*sqrt(ti_temp)
            speed = sqrt( vll*vll + 2*moment*modb*inv_mymass )
            dve   = speed*tau_spit_inv
            dvi   = vc3_tauinv/(speed*speed)
            reduction = dve + dvi
            newspeed = speed - reduction*dt
            vfrac = newspeed/speed
            !Ebench = half*mymass*newspeed*newspeed/e_charge
            !IF ((Ebench <= 1000.) .or. (Ebench <= 1.5*te_temp)) THEN !Benchmark Thermalize

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
           !v_s = half*vc3_tauinv
           speed_cube = 2*vc3_tauinv*fact_pa*dt/(speed*speed*speed) ! redefine as inverse
           zeta_o = vll/speed   ! Record the current pitch.
           CALL gauss_rand(1,zeta)  ! A random from a standard normal (1,1)
           sigma = sqrt( ABS((1.0D0-zeta_o*zeta_o)*speed_cube) ) ! The standard deviation.
           zeta_mean = zeta_o *(1.0D0 - speed_cube )  ! The new mean in the distribution.
           zeta = zeta*sigma + zeta_mean  ! The new pitch angle.
           !!The pitch angle MUST NOT go outside [-1,1] nor be NaN; but could happen accidentally with the distribution.
           IF (ABS(zeta) >  0.999D+00) zeta =  SIGN(0.999D+00,zeta)
           vll = zeta*speed

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
         USE beams3d_runtime, ONLY: t_end, lvessel, to3, lplasma_only, &
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
         INTEGER          :: ier, l
         DOUBLE PRECISION :: rinv, phi_temp, dt_local, ti_temp, ne_temp,&
                             s_temp, x0, y0, z0, xw, yw, zw, te_temp, Zeff_temp
         DOUBLE PRECISION :: qf(3),qs(3),qe(3)
         DOUBLE PRECISION :: rlocal(num_depo), plocal(num_depo), zlocal(num_depo)
         DOUBLE PRECISION :: tilocal(num_depo), telocal(num_depo), nelocal(num_depo)
         DOUBLE PRECISION :: zefflocal(num_depo)
         DOUBLE PRECISION :: tau_inv(num_depo), energy(num_depo)
         DOUBLE PRECISION :: sigvii(num_depo), sigvcx(num_depo), sigvei(num_depo)
         ! For splines
         INTEGER :: i,j,k
         REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
         INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
         REAL*8 :: fval(1)
         ! For Suzuki
         INTEGER :: A_IN(1), Z_IN(1)
         DOUBLE PRECISION :: ni_in(1)

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
                  hx     = raxis(i+1) - raxis(i)
                  hy     = phiaxis(j+1) - phiaxis(j)
                  hz     = zaxis(k+1) - zaxis(k)
                  hxi    = one / hx
                  hyi    = one / hy
                  hzi    = one / hz
                  xparam = (q(1) - raxis(i)) * hxi
                  yparam = (phi_temp - phiaxis(j)) * hyi
                  zparam = (q(3) - zaxis(k)) * hzi
                  s_temp =1.5
                  !CALL EZspline_interp(S_spl,q(1),phi_temp,q(3),s_temp,ier)
                  CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                                  hx,hxi,hy,hyi,hz,hzi,&
                                  S4D(1,1,1,1),nr,nphi,nz)
                  s_temp = fval(1)
                  IF (s_temp < one) EXIT
               END IF
               IF ((q(1) > 5*rmax)  .or. (q(1) < rmin)) THEN
                  t = t_end(myline)+dt_local
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
         DO l = 1, 3
            dt_local = stepsize(l)/q(4)
            DO
               qf = qf + myv_neut*dt_local
               q(1) = sqrt(qf(1)*qf(1)+qf(2)*qf(2))
               q(2) = ATAN2(qf(2),qf(1))
               q(3) = qf(3)
               phi_temp = MODULO(q(2), phimax)
               IF (phi_temp < 0) phi_temp = phi_temp + phimax
               !CALL EZspline_isInDomain(S_spl,q(1),phi_temp,q(3),ier)
               !IF (ier==0) THEN
               ! Assume we're in grid and only want to bug out if we're outside the grid
               IF ((q(1) >= rmin-eps1) .and. (q(1) <= rmax+eps1) .and. &
                   (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
                   (q(3) >= zmin-eps3) .and. (q(3) <= zmax+eps3)) THEN
                  i = MIN(MAX(COUNT(raxis < q(1)),1),nr-1)
                  j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
                  k = MIN(MAX(COUNT(zaxis < q(3)),1),nz-1)
                  hx     = raxis(i+1) - raxis(i)
                  hy     = phiaxis(j+1) - phiaxis(j)
                  hz     = zaxis(k+1) - zaxis(k)
                  hxi    = one / hx
                  hyi    = one / hy
                  hzi    = one / hz
                  xparam = (q(1) - raxis(i)) * hxi
                  yparam = (phi_temp - phiaxis(j)) * hyi
                  zparam = (q(3) - zaxis(k)) * hzi
                  !CALL EZspline_interp(S_spl,q(1),phi_temp,q(3),s_temp,ier)
                  CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                                  hx,hxi,hy,hyi,hz,hzi,&
                                  S4D(1,1,1,1),nr,nphi,nz)
                  s_temp = fval(1)
                  IF (s_temp > one) EXIT
               ELSE
                  EXIT
               END IF
            END DO
            ! Take a step back
            qf = qf - myv_neut*dt_local
         END DO
         ! If we're outside the grid now then we hit the grid before the wall
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
            rlocal(i) = (i-1)*(rlocal(num_depo)-rlocal(1))/REAL(num_depo-1) + rlocal(1)
            plocal(i) = (i-1)*(plocal(num_depo)-plocal(1))/REAL(num_depo-1) + plocal(1)
            zlocal(i) = (i-1)*(zlocal(num_depo)-zlocal(1))/REAL(num_depo-1) + zlocal(1)
         END DO
         plocal = MODULO(plocal, phimax)
         ! Compute temp/density along path
         DO l = 1, num_depo
            !CALL R8HERM3xyz(rlocal(l),plocal(l),zlocal(l),&
            !                TI_spl%x1(1),TI_spl%n1,&
            !                TI_spl%x2(1),TI_spl%n2,&
            !                TI_spl%x3(1),TI_spl%n3,&
            !                TI_spl%ilin1,TI_spl%ilin2,TI_spl%ilin3,&
            !                i,j,k,xparam,yparam,zparam,&
            !                hx,hxi,hy,hyi,hz,hzi,ier)
            i = MIN(MAX(COUNT(raxis < rlocal(l)),1),nr-1)
            j = MIN(MAX(COUNT(phiaxis < plocal(l)),1),nphi-1)
            k = MIN(MAX(COUNT(zaxis < zlocal(l)),1),nz-1)
            hx     = raxis(i+1) - raxis(i)
            hy     = phiaxis(j+1) - phiaxis(j)
            hz     = zaxis(k+1) - zaxis(k)
            hxi    = one / hx
            hyi    = one / hy
            hzi    = one / hz
            xparam = (rlocal(l) - raxis(i)) * hxi
            yparam = (plocal(l) - phiaxis(j)) * hyi
            zparam = (zlocal(l) - zaxis(k)) * hzi
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            TI4D(1,1,1,1),nr,nphi,nz)
            tilocal(l) = MAX(fval(1),zero)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            TE4D(1,1,1,1),nr,nphi,nz)
            telocal(l) = MAX(fval(1),zero)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            NE4D(1,1,1,1),nr,nphi,nz)
            nelocal(l) = MAX(fval(1),zero)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            ZEFF4D(1,1,1,1),nr,nphi,nz)
            zefflocal(l) = MAX(fval(1),zero)
         END DO
         tilocal = tilocal*1D-3
         telocal = telocal*1D-3
         zeff_temp = SUM(zefflocal)/DBLE(num_depo)
         IF (lsuzuki) THEN
            !--------------------------------------------------------------
            !     USE Suzuki to calcualte ionization rates
            !     Note: 10^18<ne<10^21
            !           E(keV/amu)/100 < Te < E(keV/amu)/2
            !--------------------------------------------------------------
            !Z_in(1)  = NINT(mycharge/e_charge)
            A_in(1)  = NINT(plasma_mass*inv_dalton)
            energy   = energy/(e_charge*A_in(1)) ! keV/amu
            DO l = 1, num_depo
               nelocal(l)  = MAX(MIN(nelocal(l),1E21),1E18)
               telocal(l)  = MAX(MIN(telocal(l),energy(l)/2),1.0E-3)
               ni_in(1) = nelocal(l)/zefflocal(l)
               Z_in(1)  = NINT(zefflocal(l))
               CALL suzuki_sigma(1,energy(l),nelocal(l),telocal(l),ni_in,A_in,Z_in,tau_inv(l))
            END DO
            tau_inv = tau_inv*nelocal*ABS(q(4))*1E-4
         ELSE
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
            i = MIN(MAX(COUNT(raxis < rlocal(l)),1),nr-1)
            j = MIN(MAX(COUNT(phiaxis < plocal(l)),1),nphi-1)
            k = MIN(MAX(COUNT(zaxis < zlocal(l)),1),nz-1)
            hx     = raxis(i+1) - raxis(i)
            hy     = phiaxis(j+1) - phiaxis(j)
            hz     = zaxis(k+1) - zaxis(k)
            hxi    = one / hx
            hyi    = one / hy
            hzi    = one / hz
            xparam = (rlocal(l) - raxis(i)) * hxi
            yparam = (plocal(l) - phiaxis(j)) * hyi
            zparam = (zlocal(l) - zaxis(k)) * hzi
            !CALL EZspline_interp(S_spl,rlocal(l),plocal(l),zlocal(l),s_temp,ier)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
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
               CALL uncount_wall_hit
               RETURN
            END IF
            !xlast = x0; ylast=y0; zlast=z0
            x0 = qf(1); y0 = qf(2); z0 = qf(3)
            q(1) = SQRT(qf(1)*qf(1)+qf(2)*qf(2))
            q(2) = ATAN2(qf(2),qf(1))
            q(3) = qf(3)
            IF ((q(1) > 2*rmax)  .or. (q(1) < rmin)) THEN; t = t_end(myline)+dt_local; RETURN; END IF  ! We're outside the grid
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
         REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
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
         hx     = raxis(i+1) - raxis(i)
         hy     = phiaxis(j+1) - phiaxis(j)
         hz     = zaxis(k+1) - zaxis(k)
         hxi    = one / hx
         hyi    = one / hy
         hzi    = one / hz
         xparam = (r_temp - raxis(i)) * hxi
         yparam = (phi_temp - phiaxis(j)) * hyi
         zparam = (z_temp - zaxis(k)) * hzi
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BR4D(1,1,1,1),nr,nphi,nz)
         br_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BPHI4D(1,1,1,1),nr,nphi,nz)
         bp_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BZ4D(1,1,1,1),nr,nphi,nz)
         bz_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
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
         hx     = raxis(i+1) - raxis(i)
         hy     = phiaxis(j+1) - phiaxis(j)
         hz     = zaxis(k+1) - zaxis(k)
         hxi    = one / hx
         hyi    = one / hy
         hzi    = one / hz
         xparam = (r_temp - raxis(i)) * hxi
         yparam = (phi_temp - phiaxis(j)) * hyi
         zparam = (z_temp - zaxis(k)) * hzi
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BR4D(1,1,1,1),nr,nphi,nz)
         br_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BPHI4D(1,1,1,1),nr,nphi,nz)
         bp_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BZ4D(1,1,1,1),nr,nphi,nz)
         bz_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
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
         IF (phi_temp < 0) phi_temp = phi_temp + phimax
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
         IF (bnz == 0) THEN
            e1    = (/one,zero,zero/)
            e2    = (/zero,one,zero/)
         ELSE
            e1 = (/ -bz_temp*bx_temp, -bz_temp*by_temp, one /)
            e1 = e1*binv/bnz
            e2 = (/ by_temp*e1(3)-bz_temp*e1(2), bz_temp*e1(1)-bx_temp*e1(3), bx_temp*e1(2)-by_temp*e1(1) /)
            e2 = e2*binv
            !e1 = (/ -B_temp(3)*B_temp(1)*binv/bnz, -B_temp(3)*B_temp(2)*binv/bnz, binv*bnz /)
            !e2 = binv*(/ B_temp(2)*e1(3)-B_temp(3)*e1(2), B_temp(3)*e1(1)-B_temp(1)*e1(3), B_temp(1)*e1(2)-B_temp(2)*e1(1) /)
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

END MODULE beams3d_physics_mod
