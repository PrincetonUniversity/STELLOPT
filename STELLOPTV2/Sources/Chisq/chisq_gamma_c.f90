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
      REAL(rprec) :: X, Y, Xp, Yp, dpsidr, dpsidz, Br, Bphi
      REAL(rprec) :: temp, sqrt_bbb, modbp, modbm, del 
      REAL(rprec) :: dIdb, dgdb, dbigGdb, dVdb, gamma_c, bigGamma_c
      REAL(rprec) :: dBds, dBdr, dBdphi, dBdz
      REAL(rprec) :: phip, phim, Bxt, Byt, Bzt, Brt, Bpt
      REAL(rprec) :: g,dsdx,dudx,dvdx,dsdy,dudy,dvdy,dsdz,dudz,dvdz

      REAL(rprec), DIMENSION(3) :: sflCrd, Bxyz, crossnum, crossden
      REAL(rprec), DIMENSION(3) :: e_phi, e_r, e_z, grads, gradbtest, gradB
      REAL(rprec), DIMENSION(3) :: grad_zeta, grad_psi_x_b, grad_psi_xyz
      INTEGER, PARAMETER :: nsteps = 10
      INTEGER :: delzetadiv =200
      INTEGER, PARAMETER :: bpstep = 2 !division in b'


      REAL(rprec), DIMENSION(nsteps) :: R, Z, Bx, By, Bz, Bsupv, dBsupvdpsi, modB
      REAL(rprec), DIMENSION(nsteps) :: kappa_g, e_theta_norm, grad_psi_norm, dBdpsi
      real(rprec), DIMENSION(nsteps,3) :: gradR, gradZ, grad_psi
      real(rprec), DIMENSION(nsteps,3) :: dxyzdu, dxyzdv, dxyzds
      real(rprec), DIMENSION(nsteps) :: grad_psi_i, e_theta_i, ds, dVdb_t1
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
        pi2 = 2.0_rprec*3.14159265358979_rprec
        !delzeta = Rmajor/delzetadiv !step size
        !delzeta = -0.00345175 !for comparison with ROSE
        delzeta = -0.0034433355_rprec
        rovera = sqrt(s) !rho = r/a
        psi_a = phiedge !Toroidal flux at the edge

        DO ik = 1,nsd !go through each surface
          IF (sigma(ik) >= bigno) CYCLE
          mtargets = mtargets + 1

          !get normalized toroidal flux
          s = 1.0_rprec * (ik) / (ns-1) !vmec labels are in normalized flux 
          
          CALL EZspline_interp(iota_spl,s,iota,iflag) !get iota
          CALL EZspline_derivative(iota_spl,1,s,iotap,iflag) !get iota' necessary for later
          write (*,*) 'iota', iota
          write (*,*) 'iota', iotap
          write (*,*) 'psi_a', psi_a

          !Initialize starting point
          theta = 0.0_rprec
          zeta = 0.0_rprec
          !zeta = -0.0034
          sflCrd(1) = s
          sflCrd(2) = theta
          sflCrd(3) = zeta

          
          !Do to some shitty behaviors in Stellopt, initializations may prevent
          !random crashes
          !Note that initializing to 0 instead of 0.0_rprec will definitely cause crashes
          gradR = 0.0_rprec
          gradZ = 0.0_rprec
          modB = 0.0_rprec 
          gradB = 0.0_rprec
          R = 0.0_rprec
          Z = 0.0_rprec

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

            !Calculate crude arclength
            IF (j > 1) THEN
              ds(j) = sqrt((X-Xp)*(X-Xp) + (Y-Yp)*(Y-Yp) + (Z(j) - Z(j-1))*(Z(j) - Z(j-1)))
              IF (j == 2) ds(1) = ds(j)
            END IF
            Xp = X
            Yp = Y


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
            !TODO figure out why exactly the code (copied from get_equil_B) divides
            !by nfp in the Bphi terms only. It seems weird but agrees with ROSE
            CALL get_equil_Bcyl(R(j),zeta,Z(j),Br,Bphi,Bz(j),ier)

            CALL get_equil_Bcylsuv(s,u,v,Br,Bphi,Bz(j),ier,modB(j))
            Bx(j) = Br*cos(v) - Bphi*sin(zeta)
            By(j) = Br*sin(y) + Bphi*cos(zeta)
            Bxyz(1) = Bx(j)
            Bxyz(2) = By(j)
            Bxyz(3) = Bz(j)
            !write(*,*) 'Bxyz',Bxyz(1),Bxyz(2),Bxyz(3),modB(j), ier
            
            
            
            !get_equil_xxx are supposed to return derivatives
            !but they seem to return nonsense, so let's just
            !create our own finite differences
            del = 0.005_rprec
            
            CALL get_equil_Bcylsuv(s-del,u,v,Brt,Bpt,Bzt,ier,modbm,gradbtest)
            CALL get_equil_Bcylsuv(s+del,u,v,Brt,Bpt,Bzt,ier,modbp,gradbtest)
            dBds = (modbp-modbm)/2/del
            dBdpsi(j) = dBds*pi2/psi_a !This has been verified with ROSE
            !write(*,*) 'dBdpsi',dBdpsi(j)
            
            !Handling u and v derivatives require dealing with boundaries
            !up = u+del
            !IF (up > pi2) up = up - pi2
            !CALL get_equil_Bcylsuv(s,up,v,Br,Bphi,Bz(j),ier,modbp,gradbtest)
            !um = u-del
            !IF (um < 0) um = um + pi2
            !CALL get_equil_Bcylsuv(s,um,v,Br,Bphi,Bz(j),ier,modbm,gradbtest)
            !dBdu = (modbp-modbm)/2/del
            !write(*,*) 'dBdu',dBdu,modbp,modbm
            ! 
            !vp = v+del
            !IF (vp > pi2) vp = vp - pi2
            !CALL get_equil_Bcylsuv(s,u,vp,Br,Bphi,Bz(j),ier,modbp,gradbtest)
            !vm = v-del
            !IF (vm < 0) vm = vm + pi2
            !CALL get_equil_Bcylsuv(s,u,vm,Br,Bphi,Bz(j),ier,modbm,gradbtest)
            !dBdv = (modbp-modbm)/2/del
            !write(*,*) 'dBdv', dBdv,modbp,modbm

            !Calculate dB/dr, dB/dth and dB/dz
            Brt = 0.0_rprec
            Bpt = 0.0_rprec
            Bzt = 0.0_rprec
            CALL get_equil_B(R(j)+del,MODULO(zeta, pi2), Z(j),Bxt,Byt,Bzt,ier)
            modbp = sqrt(Bxt*Bxt + Byt*Byt + Bzt*Bzt)
            CALL get_equil_B(R(j)-del,MODULO(zeta, pi2), Z(j),Bxt,Byt,Bzt,ier)
            modbm = sqrt(Bxt*Bxt + Byt*Byt + Bzt*Bzt)
            dBdr = (modbp-modbm)/2.0_rprec/del

            CALL get_equil_B(R(j),MODULO(zeta, pi2), Z(j)+del,Bxt,Byt,Bzt,ier)
            modbp = sqrt(Bxt*Bxt + Byt*Byt + Bzt*Bzt)
            CALL get_equil_B(R(j),MODULO(zeta, pi2), Z(j)-del,Bxt,Byt,Bzt,ier)
            modbm = sqrt(Bxt*Bxt + Byt*Byt + Bzt*Bzt)
            dBdz = (modbp-modbm)/2.0_rprec/del

            !TODO check to see if you need to restrict this to one field period
            phip = MODULO(zeta+del, pi2)
            CALL get_equil_B(R(j), phip, Z(j),Bxt,Byt,Bzt,ier)
            modbp = sqrt(Bxt*Bxt + Byt*Byt + Bzt*Bzt)
            phim = MODULO(zeta-del, pi2)
            CALL get_equil_B(R(j), phim, Z(j),Bxt,Byt,Bzt,ier)
            modbm = sqrt(Bxt*Bxt + Byt*Byt + Bzt*Bzt)
            dBdphi = (modbp-modbm)/2.0_rprec/del

            gradB(3) = dBdz
            gradB(1) = dBdr*X/R(j) - dBdphi*Y/R(j)
            gradB(2) = dBdr*Y/R(j) + dBdphi*X/R(j)
            
            !modbtest = modB(j)
            !CALL get_equil_Bcylsuv(s,ut+del,vt,Br,Bphi,Bz(j),ier,modB(j),gradbtest)
            
            !CALL get_equil_Bcylsuv(s,ut,vt-del,Br,Bphi,Bz(j),ier,modB(j),gradbtest)
            !modbtest = modB(j)
            !CALL get_equil_Bcylsuv(s,ut,vt+del,Br,Bphi,Bz(j),ier,modB(j),gradbtest)
            
            !write(*,*) 'modB',modB(j)
            !write(*,*) 'gradB', gradB(j,1), gradB(j,2), gradB(j,3)


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

            !geodesic curvature
            !TODO gradB is not in the correct coordinates
            kappa_g(j) = dot_product(Bxyz, gradB)/modB(j)/modB(j)
            !write (*,*) 'gradB', gradB
            !write (*,*) 'kappa_g', kappa_g(j)

            !The terms that go into gV
            grad_zeta(1) = -Y/R(j)/R(j)
            grad_zeta(2) = X/R(j)/R(j)
            grad_zeta(3) = 0.0_rprec

            !dvdB_t1 is the first term in the brackets of dVdb
            ! = iota' ( grad_psi cross b_hat) dot grad_zeta
            CALL cross_product(grad_psi(j,:), Bxyz, grad_psi_x_b)
            dVdb_t1(j) = iotap*dot_product(grad_psi_x_b, grad_zeta)/modB(j)

            !Get B^v and derivatives
            CALL get_equil_Bsupv(s,u,v,Bsupv(j),ier)
            !Use finite derivatives because ezspline derivs are broken
            CALL get_equil_Bsupv(s+del,u,v,modbp,ier)
            CALL get_equil_Bsupv(s-del,u,v,modbm,ier)
            dBsupvdpsi(j) = (modbp - modbm)/2.0_rprec/del
            dBsupvdpsi(j) = dBsupvdpsi(j)*pi2/psi_a
            !Both values verified with ROSE

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
            !DO k = 1,nwells
            !  write(*,*) 'well',k,well_start(k), well_stop(k)
            !END DO 
            !DO j = 1,nsteps
            !  write(*,*) j,B_refl,modB(j),grad_psi_norm(j),grad_psi_i(j),e_theta_norm(j),e_theta_i(j)
            !END DO

            !We've assembled all the info we need to compute the major quantities
            !for each well, we integrate the various quantities, dgdb, dGdb, dIdb, and dVdb
            DO k = 1,nwells

              dIdb = 0.0_rprec
              dgdb = 0.0_rprec
              dbigGdb = 0.0_rprec
              dVdb = 0.0_rprec

              DO j = well_start(k),well_stop(k)
                !double check that we're in a valid well
                if (grad_psi_i(j) > 1.0E8_rprec) CYCLE
                if (e_theta_i(j) == 0.0_rprec) CYCLE

                sqrt_bbb = sqrt(1 - modB(j)/B_refl)

                !dIdb
                temp = ds(j)*minB/B_refl/B_refl / sqrt_bbb
                dIdb = dIdb + temp
                
                !dgdb
                temp = ds(j) * grad_psi_norm(j) * kappa_g(j) * minB * minB
                temp = temp/B_refl/B_refl/2.0_rprec/modB(j)
                temp = temp*(sqrt_bbb + 1.0_rprec/sqrt_bbb)
                dgdb = dgdb + temp

                !dbigGdb TODO check the dBdpsi term as gradB(j,1)
                temp = dBdpsi(j) *ds(j) * minB / B_refl / B_refl / 2.0_rprec
                temp = temp*(sqrt_bbb + 1.0_rprec/sqrt_bbb)
                dbigGdb = dbigGdb + temp

                !dVdb
                temp = dVdb_t1(j) - (2.0_rprec * dBdpsi(j) - modB(j)/Bsupv(j)*dBsupvdpsi(j)) 
                temp = temp * 1.5_rprec * ds(j) / modB(j) / B_refl * sqrt_bbb
                dVdb = dVdb + temp
                

              END DO
            END DO 

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
