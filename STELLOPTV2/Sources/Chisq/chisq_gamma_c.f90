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
      USE equil_utils 
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
      INTEGER :: dex, ik, i, j, k, ier, jmin
      REAL(rprec) :: s, rovera, theta, zeta, delzeta, u, v, coszeta, sinzeta
      REAL(rprec) :: iota, iotap, minB, maxB, B_refl, psi_a, B_zeta
      REAL(rprec) :: X, Y, Xp, Yp, Z, Zp, R, dpsidr, dpsidz
      REAL(rprec) :: Bx, By, Bz, Br, Bphi
      REAL(rprec) :: Bxp, Byp, Bzp, Bx2p, By2p, Bz2p
      REAL(rprec) :: temp, sqrt_bbb, modbp, modbm, del 
      REAL(rprec) :: dIdb, dgdb, dbigGdb, dVdb, vrovervt, dloverb 
      REAL(rprec) :: gamma_c, dsoverb, wellGamma_c, bigGamma_c
      REAL(rprec) :: dBds, dBdr, dBdphi, dBdz
      REAL(rprec) :: phip, phim, Bxt, Byt, Bzt, Brt, Bpt
      REAL(rprec) :: Bxt2, Byt2, Bzt2
      REAL(rprec) :: g,dsdx,dudx,dvdx,dsdy,dudy,dvdy,dsdz,dudz,dvdz

      REAL(rprec), DIMENSION(3) :: sflCrd, Bxyz, crossnum, crossden
      REAL(rprec), DIMENSION(3) :: e_phi, e_r, e_z, grads, gradbtest, gradB
      REAL(rprec), DIMENSION(3) :: grad_zeta, grad_psi_x_b, grad_psi_xyz
      REAL(rprec), DIMENSION(3) :: bdotgradb
      INTEGER, PARAMETER :: ntransits = 400
      INTEGER, PARAMETER :: delzetadiv = 800
      INTEGER, PARAMETER :: nsteps = ntransits*delzetadiv
      INTEGER, PARAMETER :: bpstep = 80 !division in b'

      
      real(rprec), DIMENSION(3) :: gradR, gradZ, grad_psi
      real(rprec), DIMENSION(3) :: dxyzdu, dxyzdv, dxyzds

      REAL(rprec), DIMENSION(nsteps) :: ds, modB, dBdpsi, kappa_g
      REAL(rprec), DIMENSION(nsteps) :: grad_psi_norm, grad_psi_i
      REAL(rprec), DIMENSION(nsteps) :: e_theta_norm, e_theta_i
      REAL(rprec), DIMENSION(nsteps) :: Bsupv, dBsupvdpsi, dVdb_t1
      REAL(rprec), DIMENSION(nsteps, 3) :: binormal
    
      integer, parameter :: maxwells = 1000

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
      IF (iflag == 1) THEN 
        WRITE(iunit_out,'(A,2(2X,I8))') 'GAMMA_C     ',dex,4
        WRITE(iunit_out,'(A)') 'TARGET  SIGMA  GAMMA_C  K'
      END IF 
      IF (niter >= 0) THEN
        pi2 = 2.0_rprec*3.14159265358979_rprec
        delzeta = 1.0_rprec/delzetadiv !step size
        !delzeta = -0.0034433355_rprec !for comparison with ROSE
        rovera = sqrt(s) !rho = r/a
        psi_a = phiedge !Toroidal flux at the edge
        dloverb = 0.0_rprec

        DO ik = 1,nsd !go through each surface
          IF (sigma(ik) >= bigno) CYCLE
          mtargets = mtargets + 1

          !get normalized toroidal flux
          s = 1.0_rprec * (ik) / (ns-1) !vmec labels are in normalized flux 
          
          CALL EZspline_interp(iota_spl,s,iota,ier) !get iota
          CALL EZspline_derivative(iota_spl,1,s,iotap,ier) !get iota' necessary for later
          !write (*,*) 'iota', iota
          !write (*,*) 'iota', iotap
          !write (*,*) 'psi_a', psi_a

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
          Bxp = 0.0_rprec
          Byp = 0.0_rprec
          Bzp = 0.0_rprec
          Bx2p = 0.0_rprec
          By2p = 0.0_rprec
          Bz2p = 0.0_rprec

          !first time through calculate fields and some basic parameters
          DO j = 1,nsteps
            CALL pest2vmec(sflCrd) !convert to VMEC coordinates
            u = sflCrd(2)
            v = sflCrd(3)
            u = modulo(u,pi2)
            v = modulo(v,pi2)


            !write (*,*) '--------------------------------------------'
            !write (*,*) 'jsuv',j,s,u,v
            !Get rz values and gradients
            ier = 1
            DO i = 1,100
              CALL get_equil_RZ(s, u, v, R, Z, ier, gradR, gradZ)
              if (ier == 0) EXIT
            END DO
            if (ier .ne. 0) write (*,*) 'get_equil_RZ failed'
            if (i > 1) write(*,*) 'Attempts = ',i
            !X and Y values from cyclindrical
            X=R*cos(zeta)
            Y=R*sin(zeta)
            !write (*,*) 'rxyz',R(j), X, Y, Z(j)
            !write (*,*) 's,th,ze',s,theta,zeta

            !Calculate crude arclength
            !The integration over arclength is one of the key areas where stellopt/rose differ
            IF (j > 1) THEN
              ds(j) = sqrt((X-Xp)*(X-Xp) + (Y-Yp)*(Y-Yp) + (Z-Zp)*(Z-Zp))
              IF (j == 2) ds(1) = ds(j)
            END IF
            Xp = X
            Yp = Y
            Zp = Z


            !Note, the R and Z derivatives are with respect to rho and not s
            !drds = drdrho/2sqrt(s)  
            gradR(3) = gradR(3)/2/sqrt(s)
            gradZ(3) = gradZ(3)/2/sqrt(s)
            !now form dX/du, dX/dv, dX/ds where X=(x,y,z) 
            !start with the z terms, since we have them already
            !we note that in the PEST coordinates zeta = cylindrical angle phi
            !we'll use this relation a lot
            dxyzdu(3) = gradZ(1)
            dxyzdv(3) = gradZ(2)
            dxyzds(3) = gradZ(3)
            !now the x terms 
            coszeta = cos(zeta)
            sinzeta = sin(zeta)
            dxyzdu(1) = gradR(1)*coszeta
            dxyzdv(1) = gradR(2)*coszeta - R/nfp*sinzeta
            dxyzds(1) = gradR(3)*coszeta
            !and the y terms
            dxyzdu(2) = gradR(1)*sinzeta
            dxyzdv(2) = gradR(2)*sinzeta + R/nfp*coszeta
            dxyzds(2) = gradR(3)*sinzeta

            ! These derivatives have been checked and agree with ROSE
            ! however, the factor of 2*pi needs to be included because ROSE
            ! parametrizes u and v differently 
            !write(*,*) 'dxyzdu',dxyzdu(j,1)*pi2,dxyzdu(j,2)*pi2,dxyzdu(j,3)*pi2
            !write(*,*) 'dxyzdv',dxyzdv(j,1)*pi2,dxyzdv(j,2)*pi2,dxyzdv(j,3)*pi2
            !write(*,*) 'dxyzds',dxyzds(j,1),dxyzds(j,2),dxyzds(j,3)

            !Now form gradpsi
            CALL cross_product(dxyzdu, dxyzdv, crossnum)
            !write(*,*) 'crossnum', crossnum
            CALL cross_product(dxyzdv, dxyzds, crossden)
            !write(*,*) 'crossden', crossden
            grads = crossnum/dot_product(crossden, dxyzdu)
            !write(*,*) 'denom',dot_product(crossden, dxyzdu(j,:))
            grad_psi =grads * psi_a
            !write(*,*) 'grads',grads(1),grads(2),grads(3)
            !grad s agrees with ROSE, no normalizations needed
            grad_psi_norm(j) = sqrt(grad_psi(1)*grad_psi(1) + grad_psi(2)*grad_psi(2) & 
                                  & + grad_psi(3)*grad_psi(3))/pi2

            !Get B field
            !TODO figure out why exactly the code (copied from get_equil_B) divides
            !by nfp in the Bphi terms only. It seems weird but agrees with ROSE

            CALL get_equil_Bcylsuv(s,u,v,Br,Bphi,Bz,ier,modB(j))
            Bx = Br*cos(zeta) - Bphi*sin(zeta)
            By = Br*sin(zeta) + Bphi*cos(zeta)
            Bxyz(1) = Bx
            Bxyz(2) = By
            Bxyz(3) = Bz
            !write(*,*) 'Bxyz',Bxyz(1),Bxyz(2),Bxyz(3),modB(j), ier
            
            
            
            !get_equil_xxx are supposed to return derivatives
            !but they seem to return nonsense, so let's just
            !create our own finite differences
            del = 0.005_rprec

            
            CALL get_equil_Bcylsuv(s-del,u,v,Brt,Bpt,Bzt,ier,modbm,gradbtest)
            CALL get_equil_Bcylsuv(s+del,u,v,Brt,Bpt,Bzt,ier,modbp,gradbtest)
            dBds = (modbp-modbm)/2/del
            dBdpsi(j) = dBds*pi2/psi_a !This has been verified with ROSE
            !write(*,*) dBdpsi',dBdpsi(j)
            
            !Calculate dB/dr, dB/dth and dB/dz
            !Also calculate dB(xyz)d(rthetaz) which we'll need later
            Brt = 0.0_rprec
            Bpt = 0.0_rprec
            Bzt = 0.0_rprec
            CALL get_equil_B(R+del,MODULO(zeta, pi2), Z,Bxt,Byt,Bzt,ier)
            modbp = sqrt(Bxt*Bxt + Byt*Byt + Bzt*Bzt)
            CALL get_equil_B(R-del,MODULO(zeta, pi2), Z,Bxt2,Byt2,Bzt2,ier)
            modbm = sqrt(Bxt2*Bxt2 + Byt2*Byt2 + Bzt2*Bzt2)
            dBdr = (modbp-modbm)/2.0_rprec/del

            CALL get_equil_B(R,MODULO(zeta, pi2), Z+del,Bxt,Byt,Bzt,ier)
            modbp = sqrt(Bxt*Bxt + Byt*Byt + Bzt*Bzt)
            CALL get_equil_B(R,MODULO(zeta, pi2), Z-del,Bxt2,Byt2,Bzt2,ier)
            modbm = sqrt(Bxt2*Bxt2 + Byt2*Byt2 + Bzt2*Bzt2)
            dBdz = (modbp-modbm)/2.0_rprec/del
            

            phip = MODULO(zeta+del, pi2)
            CALL get_equil_B(R, phip, Z,Bxt,Byt,Bzt,ier)
            modbp = sqrt(Bxt*Bxt + Byt*Byt + Bzt*Bzt)
            phim = MODULO(zeta-del, pi2)
            CALL get_equil_B(R, phim, Z,Bxt2,Byt2,Bzt2,ier)
            modbm = sqrt(Bxt2*Bxt2 + Byt2*Byt2 + Bzt2*Bzt2)
            dBdphi = (modbp-modbm)/2.0_rprec/del


            gradB(3) = dBdz
            gradB(1) = dBdr*X/R - dBdphi*Y/R
            gradB(2) = dBdr*Y/R + dBdphi*X/R

            
            !Calculate |e_theta| 

            !e_r
            e_r(1) = X/R
            e_r(2) = Y/R
            e_r(3) = 0.0_rprec

            !e_phi
            e_phi(1) = -Y/R
            e_phi(2) = X/R
            e_phi(3) = 0.0_rprec

            !e_z (can move this outside the loop if desired for speed)
            e_z(1) = 0.0
            e_z(2) = 0.0
            e_z(3) = 1.0
            !write(*,*) 'e_r', e_r
            !write(*,*) 'e_z', e_z
            !write(*,*) 'e_phi', e_phi

            dpsidr = dot_product(grad_psi, e_r)
            dpsidz = dot_product(grad_psi, e_z)
            !B_zeta = dot_product(Bxyz, e_phi)

            !Note that e_theta_norm is negative iff Bphi is negative
            !this is the case in QH. However, all this does is multiply the
            !entire expression by -1 of gamma_c, and since the integrated 
            !quantity in Gamma_c is squared, this doesn't matter
            !the division by 2pi is needed to match with rose convention
            e_theta_norm(j) = sqrt(dpsidr*dpsidr + dpsidz*dpsidz)/Bphi/pi2
            !write(*,*) j,'e_theta_norm, B_zeta, Bphi',e_theta_norm(j),Bphi

            !geodesic curvature
            !calculate bdotgradb = partial b/ partial s
            bdotgradb = 0.0_rprec
            kappa_g(j) = 0.0_rprec
            binormal(j,:) = 0.0_rprec
            CALL cross_product(grad_psi/grad_psi_norm(j)/pi2, Bxyz/modB(j), binormal(j,:))
            !write (*,*) 'tangential', Bxyz/modB(j)
            !write (*,*) 'normal',grad_psi(j,:)/grad_psi_norm(j)
            !write (*,*) 'nxb',binormal(j,:)
            !write (*,*) 'bdotgradb',bdotgradb
            !write (*,*) 'kappa_g', kappa_g(j)

            IF (j > 2) THEN
              bdotgradb(1) = (Bx/modB(j) - Bx2p/modB(j-2))/(ds(j-1)+ds(j))
              bdotgradb(2) = (By/modB(j) - By2p/modB(j-2))/(ds(j-1)+ds(j))
              bdotgradb(3) = (Bz/modB(j) - Bz2p/modB(j-2))/(ds(j-1)+ds(j))
              !write (*,*) 'bdotgradb2 j=',j-1,bdotgradb
              kappa_g(j-1) = dot_product(bdotgradb, binormal(j-1, :))
              !write (*,*) 'kappa_g2 j=',j-1, kappa_g(j-1)
              IF (j == nsteps) THEN
                bdotgradb(1) = (Bx/modB(j) - Bxp/modB(j-1))/(ds(j))
                bdotgradb(2) = (By/modB(j) - Byp/modB(j-1))/(ds(j))
                bdotgradb(3) = (Bz/modB(j) - Bzp/modB(j-1))/(ds(j))
                kappa_g(j) = dot_product(bdotgradb, binormal(j,:))
                !write (*,*) 'bdotgrad b, j=',nsteps,bdotgradb
                !write (*,*) 'kappa_g, j=',nsteps, kappa_g(j)
              END IF  
               
            END IF
            IF (j == 2) THEN !handle the first index
              bdotgradb(1) = (Bx/modB(j) - Bxp/modB(j-1))/(ds(j-1))
              bdotgradb(2) = (By/modB(j) - Byp/modB(j-1))/(ds(j-1))
              bdotgradb(3) = (Bz/modB(j) - Bzp/modB(j-1))/(ds(j-1))

              kappa_g(j-1) = dot_product(bdotgradb, binormal(j-1,:))
              !write (*,*) 'bdotgrad b, j=1',bdotgradb
              !write (*,*) 'kappa_g, j=1', kappa_g(j-1)
            END IF 
            Bx2p = Bxp
            By2p = Byp
            Bz2p = Bzp
            Bxp = Bx
            Byp = By
            Bzp = Bz

            !kappa_g(j) = dot_product(Bxyz, gradB)/modB(j)/modB(j)
            !write (*,*) 'gradB', gradB
            !write (*,*) 'kappa_g', kappa_g(j)

            !The terms that go into gV
            grad_zeta(1) = -Y/R/R
            grad_zeta(2) = X/R/R
            grad_zeta(3) = 0.0_rprec

            !dvdB_t1 is the first term in the brackets of dVdb
            ! = iota' ( grad_psi cross b_hat) dot grad_zeta
            CALL cross_product(grad_psi, Bxyz, grad_psi_x_b)
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

            !integrate dloverb
            dloverb = dloverb + ds(j)/modB(j)

            !advance the step
            zeta = zeta + delzeta
            theta = theta + (iota * delzeta)
            !zeta = modulo(zeta, pi2)
            !theta = modulo(theta, pi2)
            sflCrd(1) = s !s should be constant, but in case of numerical precision errors
            sflCrd(2) = theta
            sflCrd(3) = zeta
          END DO


          bigGamma_c = 0.0_rprec
          minB = MINVAL(modB)
          maxB = MAXVAL(modB)
          
          !Make the bp array
          DO i=1,bpstep
            bp(i) = 1.0_rprec + (maxB-minB)/minB * (i-0.5_rprec) / bpstep 
          END DO
          deltabp = (maxB - minB)/minB/bpstep 

          !Go through each value of bp
          DO i=1,bpstep
            grad_psi_i = 1.0E10_rprec !This term appears in the denominator
            e_theta_i = 0

            !calculate the reflecting field
            B_refl = minB*bp(i)
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
                well_stop(cur_well) = j-2

                !test output
                !write(*,*) 'exiting well',cur_well
                !write(*,*) 'well start,stop',well_start(cur_well), well_stop(cur_well)
                !write(*,*) 'minima',grad_psi_min, e_theta_min

                !set gradpsi and etheta for all values
                grad_psi_i(well_start(cur_well):well_stop(cur_well)) = grad_psi_min
                e_theta_i(well_start(cur_well):well_stop(cur_well)) = e_theta_min
                !print out the Bmin
                !write (*,*) "Bmin, zeta for well ",cur_well, cur_Bmin, zeta*4.0_rprec/pi2
         

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
                  well_start(cur_well) = j+1
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
                  zeta = j*delzeta
                  jmin = j
                END IF
              END IF
            
            END DO
            !If we ended in the middle of a well, we'll have one more well_start
            !than well stop, let's correct that
            IF (well_start(cur_well) .NE. 0) THEN
              well_start(cur_well) = 0
            END IF 
            
            nwells = cur_well - 1
            !write(*,*) 'well number',nwells
            !DO k = 1,nwells
            !  write(*,*) 'well',k,well_start(k), well_stop(k)
            !END DO 
            !DO j = 1,nsteps
            !  write(*,*) j,B_refl,modB(j),grad_psi_norm(j),grad_psi_i(j),e_theta_norm(j),e_theta_i(j)
            !END DO

            !We've assembled all the info we need to compute the major quantities
            !for each well, we integrate the various quantities, dgdb, dGdb, dIdb, and dVdb
            gamma_c = 0.0_rprec
            vrovervt = 0.0_rprec
            wellgamma_c = 0.0_rprec
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
                temp = ds(j)/2.0_rprec/minB/bp(i)/bp(i) / sqrt_bbb
                dIdb = dIdb + temp
                
                !dgdb
                temp = ds(j) * grad_psi_norm(j) * kappa_g(j) 
                temp = temp/bp(i)/bp(i)/2.0_rprec/modB(j)
                temp = temp*(sqrt_bbb + 1.0_rprec/sqrt_bbb)
                dgdb = dgdb + temp

                !dbigGdb 
                temp = dBdpsi(j) *ds(j) /B_refl / bp(i) / modB(j) / 2.0_rprec
                temp = temp*(sqrt_bbb + 1.0_rprec/sqrt_bbb)
                dbigGdb = dbigGdb + temp

                !dVdb
                temp = dVdb_t1(j) - (2.0_rprec * dBdpsi(j) - modB(j)/Bsupv(j)*dBsupvdpsi(j)) 
                temp = temp * 1.5_rprec * ds(j) / modB(j) / B_refl * sqrt_bbb
                dVdb = dVdb + temp
              

              END DO !end integration over a single well


              !vrovervt ratio of radial to poloidal drifts
              !write (*,*) '---------------------------------------'
              !write (*,*) 'well k,b ', k, B_refl,well_start(k),well_stop(k)
              !write (*,*) 'dIdb ', dIdb
              !write (*,*) 'dgdb ' , dgdb
              !write (*,*) 'dbigGdb', dbigGdb
              !write (*,*) 'dVdb ', dVdb
              !write (*,*) 'etheta0 ', e_theta_i(j-1)
              !write (*,*) 'grad_psi ', grad_psi_i(j-1)
              IF (well_start(k) < well_stop(k)) THEN

                temp = dgdb/grad_psi_i(well_start(k))/dIdb / minB / e_theta_i(well_start(k))
                temp = temp / (dbigGdb/dIdb + 0.666666_rprec * dVdb/dIdb)
                vrovervt = temp
              ELSE
                vrovervt = 0.0
              END IF
              !write (*,*) 'vrovervt ', vrovervt

              gamma_c = 4.0_rprec/pi2 * atan(vrovervt)
              wellGamma_c = wellGamma_c + (gamma_c * gamma_c * dIdb)
              !write (*,*) 'wellGamma_c ', wellGamma_c
              !write (*,*) 'gamma_c, wellGamma_c', gamma_c, wellGamma_c
            END DO !end sum over all wells           
            bigGamma_c = bigGamma_c + wellGamma_c * pi2/4.0_rprec/sqrt(2.0_rprec)*deltabp
          END DO !end integration over bp
          bigGamma_c = bigGamma_c/dloverb

          vals(mtargets) = bigGamma_c
          sigmas(mtargets) = ABS(sigma(ik))
          targets(mtargets) = target(ik)
          IF (iflag ==1) THEN
            WRITE(iunit_out,'(3ES22.12E3,3(1X,I5))') targets(mtargets),sigmas(mtargets),vals(mtargets),ik
          END IF

          
        END DO !end loop over flux surfaces
      ELSE !This is the initialization loop that just counts targets
        DO ik = 1, nsd
          IF (sigma(ik) < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) THEN
              target_dex(mtargets)=jtarget_gamma_c
            END IF
          END IF
        END DO
      END IF


      END SUBROUTINE chisq_gamma_c
