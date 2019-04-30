!-----------------------------------------------------------------------
!     Subroutine:    chisq_gamma_c
!     Authors:       Aaron Bader
!     Date:          04/23/2019
!     Description:   Calculate gamma_c metric for energetic particle 
!                    transport
!                    Temporarily, use helicity instead
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
      REAL(rprec) :: s, theta, zeta, delzeta, u, v
      REAL(rprec) :: iota, iotap, minB, maxB

      REAL(rprec), DIMENSION(3) :: sflCrd
      INTEGER, PARAMETER :: nsteps = 100
      INTEGER :: delzetadiv =200
      INTEGER, PARAMETER :: bpstep = 101 !division in b'


      REAL(rprec), DIMENSION(nsteps) :: R, Z, Bx, By, Bz, modB
      real(rprec), DIMENSION(nsteps,3) :: gradR, gradZ, gradB

      REAL(rprec) :: deltabp
      REAL(rprec), DIMENSION(bpstep) :: bp
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      dex = COUNT(sigma < bigno) !count of surfaces to evaluate at
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I8))') 'GAMMA_C     ',dex,3
 
      IF (niter >= 0) THEN
        delzeta = Rmajor/delzetadiv !step size
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

            !Get B field
            CALL get_equil_B(R(j), zeta, Z(j), Bx(j), By(j), Bz(j), ier, modB(j), gradB(j,:))
            !Max and min along field line


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
          DO j=1,bpstep
            bp(j) = (1.0_rprec*(j-1))/(bpstep-1)
          END DO

          !DO j=1,nsteps
          !  write(*,*) j, R(j), Z(j), modB(j)
          !END DO
          !write(*,*) 'MinB, maxB'
          !write(*,*) minB, maxB
        END DO
      ELSE !This is the initialization loop that just counts targets
        DO ik = 1, nsd
          IF (sigma(ik) < bigno) THEN
            mtargets = mtargets + 1
          END IF
        END DO
      END IF


      END SUBROUTINE chisq_gamma_c
