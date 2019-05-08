!-----------------------------------------------------------------------
!     Subroutine:    chisq_gamma_c
!     Authors:       Aaron Bader
!     Date:          04/23/2019
!     Description:   Calculate gamma_c metric for energetic particle 
!                    transport
!                    
!-----------------------------------------------------------------------
      SUBROUTINE chisq_gamma_c(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals
      USE stel_tools
      USE stellopt_input_mod
      USE read_boozer_mod
      USE read_wout_mod, ONLY: ns
      USE math_utilities
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
      INTEGER :: dex, ik, i, j, k, ier
      REAL(rprec) :: s, rovera, theta, zeta, delzeta, u, v, coszeta, sinzeta
      REAL(rprec) :: iota, iotap, minB, maxB, B_refl, psi_a, B_zeta
      REAL(rprec) :: X, Y, dpsidr, dpsidz, Br, Bphi

      REAL(rprec), DIMENSION(3) :: sflCrd, Bxyz, crossnum, crossden
      REAL(rprec), DIMENSION(3) :: e_phi, e_r, e_z, grads
      INTEGER, PARAMETER :: nsteps = 10000
      INTEGER :: delzetadiv =200
      INTEGER, PARAMETER :: bpstep = 2 !division in b'


      REAL(rprec), DIMENSION(nsteps) :: R, Z, Bx, By, Bz, Bsupv, dBsupvdpsi, modB
      REAL(rprec), DIMENSION(nsteps) :: kappa_g, e_theta_norm, grad_psi_norm
      real(rprec), DIMENSION(nsteps,3) :: gradR, gradZ, gradB, grad_psi
      real(rprec), DIMENSION(nsteps,3) :: dxyzdu, dxyzdv, dxyzds, dBsupv 
      real(rprec), DIMENSION(nsteps) :: grad_psi_i, e_theta_i
      integer, parameter :: maxwells = 100

      integer, dimension(maxwells) :: well_start, well_stop
      integer :: in_well, cur_well, nwells

      REAL(rprec) :: deltabp, den
      REAL(rprec), DIMENSION(bpstep) :: bp

      real(rprec) :: grad_psi_min, e_theta_min, cur_Bmin
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      dex = COUNT(sigma < bigno) !count of surfaces to evaluate at
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I8))') 'GAMMA_C     ',dex,3
 
      IF (niter >= 0) THEN
        pi2 = 2*3.14159265358979
        delzeta = Rmajor/delzetadiv !step size
        !delzeta = -0.0034518 for comparison with ROSE
        rovera = sqrt(s) !rho = r/a
        psi_a = phiedge !Toroidal flux at the edge
        !nsteps = 100 !somewhat arbitrary number of steps

        DO ik = 1,nsd !go through each surface
          IF (sigma(ik) >= bigno) CYCLE
          mtargets = mtargets + 1

          !get normalized toroidal flux
          s = 1.0_rprec * (ik) / (ns-1) !vmec labels are in normalized flux 
          
          CALL EZspline_interp(iota_spl,s,iota,iflag) !get iota
          CALL EZspline_derivative(iota_spl,1,s,iotap,iflag) !get iota' necessary for later
          !write (*,*) 'iota', iota
          !write (*,*) 'iota', iotap

          !Initialize starting point
          theta = 0.0
          zeta = 0.0
          !zeta = -0.0034
          sflCrd(1) = s
          sflCrd(2) = theta
          sflCrd(3) = zeta

          !first time through calculate fields and some basic parameters
          DO j = 1,nsteps
            CALL pest2vmec(sflCrd) !convert to VMEC coordinates
            u = sflCrd(2)
            v = sflCrd(3)
            u = modulo(u,pi2)
            v = modulo(v,pi2)


            !Get rz values and gradients
            CALL get_equil_RZ(s, u, v, R(j), Z(j), ier, gradR(j,:), gradZ(j,:))

            !X and Y values from cyclindrical
            X=R(j)*cos(zeta)
            Y=R(j)*sin(zeta)
            !write (*,*) '--------------'
            !write (*,*) j,s,u,v
            !write (*,*) s,theta,zeta
            !write (*,*) R(j), X, Y, Z(j)
            !Note, the R and Z derivatives are with respect to rho and not s
            !drds = drdrho/2sqrt(s)  
            gradR(j,3) = gradR(j,3)/2/sqrt(s)
            gradZ(j,3) = gradZ(j,3)/2/sqrt(s)
            !now form dX/du, dX/dv, dX/ds where X=(x,y,z) 
            !start with the z terms, since we have them already
            !we note that in the PEST coordinates zeta = cylindrical angle phi
            !we'll use this relation a lot
            dxyzdu(j,3) = gradZ(j,1)
            dxyzdv(j,3) = gradZ(j,2)
            dxyzds(j,3) = gradZ(j,3)
            !now the x terms 
            coszeta = cos(zeta)
            sinzeta = sin(zeta)
            dxyzdu(j,1) = gradR(j,1)*coszeta
            dxyzdv(j,1) = gradR(j,2)*coszeta - R(j)/nfp*sinzeta
            dxyzds(j,1) = gradR(j,3)*coszeta
            !and the y terms
            dxyzdu(j,2) = gradR(j,1)*sinzeta
            dxyzdv(j,2) = gradR(j,2)*sinzeta + R(j)/nfp*coszeta
            dxyzds(j,2) = gradR(j,3)*sinzeta

            ! These derivatives have been checked and agree with ROSE
            ! however, the factor of 2*pi needs to be included because ROSE
            ! parametrizes u and v from 
            !write(*,*) 'dxyzdu',dxyzdu(j,1)*pi2,dxyzdu(j,2)*pi2,dxyzdu(j,3)*pi2
            !write(*,*) 'dxyzdv',dxyzdv(j,1)*pi2,dxyzdv(j,2)*pi2,dxyzdv(j,3)*pi2
            !write(*,*) 'dxyzds',dxyzds(j,1),dxyzds(j,2),dxyzds(j,3)

            !Now form gradpsi
            CALL cross_product(dxyzdu(j,:), dxyzdv(j,:), crossnum)
            !write(*,*) 'crossnum', crossnum
            CALL cross_product(dxyzdv(j,:), dxyzds(j,:), crossden)
            !write(*,*) 'crossden', crossden
            grads = crossnum/dot_product(crossden, dxyzdu(j,:))
            !write(*,*) 'denom',dot_product(crossden, dxyzdu(j,:))
            grad_psi(j,:) =grads * psi_a
            !write(*,*) 'grads',grads(1),grads(2),grads(3)
            !grad s agrees with ROSE, no normalizations needed
            grad_psi_norm(j) = sqrt(grad_psi(j,1)*grad_psi(j,1) + grad_psi(j,2)*grad_psi(j,2) & 
                                  & + grad_psi(j,3)*grad_psi(j,3))

            !Get B field
            !TODO there's some weird stuff going on here. get_equil_B sometimes
            !fails and gives 0.00 for the fields. So instead, I call the suv form
            !which is what get_equil_B calls, and then it doesn't fail
            !TODO figure out why exactly the code (copied from get_equil_B) divides
            !by nfp in the Bphi terms only. It seems weird
            !CALL get_equil_B(R(j), zeta, Z(j), Bx(j), By(j), Bz(j), ier, modB(j), gradB(j,:))
            CALL get_equil_Bcylsuv(s,u,v,Br,Bphi,Bz(j),ier,modB(j),gradB(j,:))
            Bx(j) = Br*cos(v) - Bphi*sin(v/nfp)
            By(j) = Br*sin(v) + Bphi*cos(v/nfp)
            Bxyz(1) = Bx(j)
            Bxyz(2) = By(j)
            Bxyz(3) = Bz(j)


            !write(*,*) 'modB',modB(j)
            !write(*,*) 'Bxyz',Bxyz(1),Bxyz(2),Bxyz(3),modB(j), ier

            !Get gaussian curvature
            CALL get_equil_kappa(s,u,v,kappa_g(j),ier)
            !write(*,*) 'kappa_g',kappa_g(j)
            !kappa_g disagrees with ROSE, need to look further
            !I think this kappa is the curvature of the surface, not the field line

            !Calculate |e_theta| 

            !e_r
            e_r(1) = X/R(j)
            e_r(2) = Y/R(j)
            e_r(3) = 0.0_rprec

            !e_phi
            e_phi(1) = -Y/R(j)
            e_phi(2) = X/R(j)
            e_phi(3) = 0.0_rprec

            !e_z (can move this outside the loop if desired for speed)
            e_z(1) = 0.0
            e_z(2) = 0.0
            e_z(3) = 1.0
            !write(*,*) 'e_r', e_r
            !write(*,*) 'e_z', e_z
            !write(*,*) 'e_phi', e_phi

            dpsidr = dot_product(grad_psi(j,:), e_r)
            dpsidz = dot_product(grad_psi(j,:), e_z)
            !B_zeta = dot_product(Bxyz, e_phi)

            !Note that e_theta_norm is negative iff Bphi is negative
            !this is the case in QH. However, all this does is multiply the
            !entire expression by -1 of gamma_c, and since the integrated 
            !quantity in Gamma_c is squared, this doesn't matter
            e_theta_norm(j) = sqrt(dpsidr*dpsidr + dpsidz*dpsidz)/Bphi
            !write(*,*) j,'e_theta_norm, B_zeta, Bphi',e_theta_norm(j),Bphi


            !Get B^v and derivatives
            CALL get_equil_Bsupv(s,u,v,Bsupv(j),ier,dBsupv(j,:))
            dBsupvdpsi(j) = dBsupv(j,3)*2*s/psi_a
            !write(*,*) 's',s,'psi_a',psi_a
            !write(*,*) 'Bsupv',Bsupv(j)
            !write(*,*) 'dBsupvdpsi', dBsupvdpsi(j)
            !DBsupvdpsi differs with ROSE, but I'm not sure which is correct

            !Rose version of BPHI, verified that these agree
            !Bsupv(j) = dot_product(Bxyz,e_phi)/R(j)
            !write(*,*) 'Bsupv rose',Bsupv(j)

            !advance the step
            zeta = zeta + delzeta
            theta = theta + (iota * delzeta)
            !zeta = modulo(zeta, pi2)
            !theta = modulo(theta, pi2)
            sflCrd(1) = s !I don't think this guy changes, but just in case
            sflCrd(2) = theta
            sflCrd(3) = zeta
          END DO
          minB = MINVAL(modB)
          maxB = MAXVAL(modB)
          
          !Make the bp array
          DO i=1,bpstep
            bp(i) = (1.0_rprec*(i-1))/(bpstep-1)
          END DO
          deltabp = (maxB - minB)/(bpstep+2) 

          !Go through each value of bp
          DO i=1,bpstep
            grad_psi_i = 1.0E10_rprec !This term appears in the denominator
            e_theta_i = 0

            !calculate the reflecting field
            B_refl = minB + i*deltabp
            !we follow along the line marking the well beginning and the well ends
            !when we have a well, we calculate the minimum of grad_psi (precomputed above)
            !then we set all values of grad_psi_i for that well range to the minimum
            !we also set all the well_mask indices to 1.
            in_well = 0

            !write (*,*) 'beginning bp',i
            !write (*,*) 'Bmin, Bmax, B_refl', minB, maxB, B_refl
            
            well_start = 0 !this is the array of all the indices where wells begin
            well_stop = 0 !this is the array of all the indices where wells end
            cur_well = 1 !this is the index of the current well

            !these values are calculated at the well minima
            grad_psi_min = 1.0E10_rprec
            e_theta_min = 0
            cur_Bmin = B_refl !this should be the maximum value in any well

            DO j=1,nsteps
              ! Not in a well and shouldn't be 
              IF ((in_well == 0) .and. (modB(j) > B_refl)) CYCLE

              ! We are in a well, but we just passed a boundary value
              IF ((in_well == 1) .and. (modB(j) > B_refl)) THEN
      
                ! reset the well status
                in_well = 0
            
                ! In this case we have exited a well that hasn't been started properly
                ! so we ignore it
                IF (well_start(cur_well) == 0) THEN
                  cur_Bmin = B_refl
                  grad_psi_min = 1.0E10_rprec
                  e_theta_min = 0
                  CYCLE
                END IF

                !mark the well end index
                well_stop(cur_well) = j-1

                !test some output
                !write(*,*) 'exiting well',cur_well
                !write(*,*) 'well start,stop',well_start(cur_well), well_stop(cur_well)
                !write(*,*) 'minima',grad_psi_min, e_theta_min

                !set gradpsi and etheta for all values
                grad_psi_i(well_start(cur_well):well_stop(cur_well)) = grad_psi_min
                e_theta_i(well_start(cur_well):well_stop(cur_well)) = e_theta_min

                !DO k = well_start(cur_well),well_stop(cur_well)
                  !write(*,*) k,grad_psi_i(k), e_theta_i(k)
                !END DO

                !reset minimum of modB, psi and e
                cur_Bmin = B_refl
                e_theta_min = 0
                grad_psi_min = 1.0E10_rprec

                cur_well = cur_well + 1
                
                !exit if we have gone through too many wells
                IF (cur_well >= maxwells) EXIT
                CYCLE

              END IF   
              ! We are not in a well, but are entering one
              IF (in_well == 0 .and. modB(j) < B_refl) THEN
                

                !add the beginning value, but only if we're not at the start
                IF (j > 1) THEN
                  well_start(cur_well) = j
                END IF 
                !mark that we're in a well
                in_well = 1
              END IF

              ! We are in a well and should be in one
              ! Note this gets executed if the previous IF statement is true
              IF (in_well == 1 .and. modB(j) < B_refl) THEN

                !check modB and update minima
                IF (modB(j) < cur_Bmin) THEN
                  cur_Bmin = modB(j)
                  e_theta_min = e_theta_norm(j)
                  grad_psi_min = grad_psi_norm(j)
                END IF
              END IF
            
            END DO
            !If we ended in the middle of a well, we'll have one more well_start
            !than well stop, let's correct that
            IF (well_start(cur_well) .NE. 0) THEN
              well_start(cur_well) = 0
            END IF 
            
            nwells = cur_well - 1

            !Write some test output
            !write(*,*) 'well number',nwells
            !DO j = 1,nwells
            !  write(*,*) 'well',j,well_start(j), well_stop(j)
            !END DO 
            !DO j = 1,nsteps
            !  write(*,*) j,B_refl,modB(j),grad_psi_norm(j),grad_psi_i(j),e_theta_norm(j),e_theta_i(j)
            !END DO

          END DO
        END DO
      ELSE !This is the initialization loop that just counts targets
        DO ik = 1, nsd
          IF (sigma(ik) < bigno) THEN
            mtargets = mtargets + 1
          END IF
        END DO
      END IF


      END SUBROUTINE chisq_gamma_c
