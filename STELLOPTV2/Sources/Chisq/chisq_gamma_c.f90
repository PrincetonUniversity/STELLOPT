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
      INTEGER :: dex, ik, i, j, ier
      REAL(rprec) :: s, rovera, theta, zeta, delzeta, u, v, coszeta, sinzeta
      REAL(rprec) :: iota, iotap, minB, maxB, B_refl, psi_a

      REAL(rprec), DIMENSION(3) :: sflCrd, crossnum, crossden
      INTEGER, PARAMETER :: nsteps = 100
      INTEGER :: delzetadiv =200
      INTEGER, PARAMETER :: bpstep = 101 !division in b'


      REAL(rprec), DIMENSION(nsteps) :: R, Z, Bx, By, Bz, Bsupv, dBsupvdpsi, modB
      REAL(rprec), DIMENSION(nsteps) :: kappa_g
      real(rprec), DIMENSION(nsteps,3) :: gradR, gradZ, gradB, grads
      real(rprec), DIMENSION(nsteps,3) :: dxyzdu, dxyzdv, dxyzds, dBsupv 
      real(rprec), DIMENSION(nsteps) :: grad_psi_i
      integer, dimension(nsteps) :: well_mask

      REAL(rprec) :: deltabp, den
      REAL(rprec), DIMENSION(bpstep) :: bp
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      dex = COUNT(sigma < bigno) !count of surfaces to evaluate at
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I8))') 'GAMMA_C     ',dex,3
 
      IF (niter >= 0) THEN
        pi2 = 2*3.14159265358979
        delzeta = Rmajor/delzetadiv !step size
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

          !Initialize starting point
          theta = 0.0
          zeta = 0.0
          sflCrd(1) = s
          sflCrd(2) = theta
          sflCrd(3) = zeta

          !first time through calculate fields and some basic parameters
          DO j = 1,nsteps
            CALL pest2vmec(sflCrd) !convert to VMEC coordinates
            u = sflCrd(2)
            v = sflCrd(3)

            !Get rz values and gradients, store gradients in temporary arrays
            CALL get_equil_RZ(s, u, v, R(j), Z(j), ier, gradR(j,:), gradZ(j,:))
            !the R and Z derivatives are with respect to rho and not s
            !drds = drdrho/2sqrt(s)  
            gradR(j,3) = gradR(j,3)/2/sqrt(s)
            gradZ(j,3) = gradZ(j,3)/2/sqrt(s)
            !now form dX/du, dX/dv, dX/ds where X=(x,y,z) 
            !start with the z terms, since we have them already
            !we note that in the PEST coordinates zeta = cylindrical angle phi
            dxyzdu(j,3) = gradZ(j,1)
            dxyzdv(j,3) = gradZ(j,2)
            dxyzds(j,3) = gradZ(j,3)
            !now the x terms
            coszeta = cos(zeta)
            sinzeta = sin(zeta)
            dxyzdu(j,1) = gradR(j,1)*coszeta
            dxyzdv(j,1) = gradR(j,2)*coszeta - pi2*R(j)/nfp*sinzeta
            dxyzds(j,1) = gradR(j,3)*coszeta
            !and the y terms
            dxyzdu(j,2) = gradR(j,1)*sinzeta
            dxyzdv(j,2) = gradR(j,2)*sinzeta + pi2*R(j)/nfp*coszeta
            dxyzds(j,2) = gradR(j,3)*sinzeta

            CALL cross_product(dxyzdu, dxyzdv, crossnum)
            CALL cross_product(dxyzdv, dxyzds, crossden)
            grads(j,:) = crossnum/dot_product(crossden, dxyzdu(j,:))

            !Get B field
            CALL get_equil_B(R(j), zeta, Z(j), Bx(j), By(j), Bz(j), ier, modB(j), gradB(j,:))

            !Get B^v and derivatives
            CALL get_equil_Bsupv(s,u,v,Bsupv(j),ier,dBsupv(j,:))
            dBsupvdpsi(j) = dBsupv(j,1)*2*s/psi_a

            !Get gaussian curvature
            CALL get_equil_kappa(s,u,v,kappa_g(j),ier)


            !advance the step
            zeta = zeta + delzeta
            theta = theta + (iota * delzeta)
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

          !Go through each value of bp
          DO i=1,bpstep
            grad_psi_i = 1.0E10_rprec !This term appears in the denominator
            well_mask = 0

            !calculate the reflecting field
            B_refl = minB + (i-1)*deltabp
            !we follow along the line marking the well beginning and the well ends
            !when we have a well, we calculate the minimum of grad_psi (precomputed above)
            !then we set all values of grad_psi_i for that well range to the minimum
            !we also set all the well_mask indices to 1.

            !DO j=1,nsteps
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
