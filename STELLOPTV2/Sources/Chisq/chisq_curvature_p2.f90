!-----------------------------------------------------------------------
!     Subroutine:    chisq_curvature_P2
!     Authors:       A. Bader (lazerson@pppl.gov)
!     Date:          02/11/2013
!     Description:   This subroutine targets the second principal
!           curvature to reduce concave regions
!-----------------------------------------------------------------------
      SUBROUTINE chisq_curvature_P2(target,sigma,niter,iflag)
      !-----------------------------------------------------------------------
      !     Libraries
      !-----------------------------------------------------------------------

      USE stellopt_runtime
      USE stellopt_targets  
      USE read_wout_mod, ONLY: m => xm, n => xn, &
          rmnc_full => rmnc, zmns_full => zmns, nfp, &
          num_coeffs => mnmax, ns

      !-----------------------------------------------------------------------
      !     Input/Output Variables
      !
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag

      INTEGER, parameter :: np=128, nt=128 !hard code resolution for now
      REAL(rprec), DIMENSION(np, nt) :: x,y,z,r
      REAL(rprec), DIMENSION(np, nt) :: xu,xv,yu,yv,zu,zv !first derivs
      REAL(rprec), DIMENSION(np, nt) :: xuu, xuv, yuu, yuv, zuu, zuv !2nd d
      REAL(rprec), DIMENSION(np, nt) :: xvv, yvv, zvv !2nd d
      REAL(rprec), DIMENSION(np, nt) :: EE,FF,GG,Nx,Ny,Nz,magnorm
      REAL(rprec), DIMENSION(np, nt) :: LL,MM,NN
      REAL(rprec), DIMENSION(np, nt) ::HH,KK,P1,P2
      REAL(rprec), DIMENSION(num_coeffs) :: rmnc, zmns
      REAL(rprec) :: cosangle, sinangle, cosphi, sinphi
      INTEGER :: i,j,k
      REAL(rprec) :: u,v, angle, phi


      !----------------BEGIN SUBROUTINE --------------

      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'CURVATURE_P2 ',1,5
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  P1  P2'
      IF (niter >= 0) THEN
        !Only consider the boundary
        rmnc = rmnc_full(:, ns)
        zmns = zmns_full(:, ns)

        !---------------calculate x,y,z, and derivs--------------!
        !initialize
        r = 0
        x = 0
        y = 0
        z = 0
        xu = 0
        xv = 0
        xuu = 0
        xuv = 0
        xvv = 0
        yu = 0
        yv = 0
        yuu = 0
        yuv = 0
        yvv = 0
        zu = 0
        zv = 0
        zuu = 0
        zuv = 0
        zvv = 0

        DO i = 1,np
            u = REAL(i-1)/(np - 1)*2*pi
            DO j = 1,nt
                v = REAL(j-1)/(nt - 1)*2*pi
                cosphi = cos(v)
                sinphi = sin(v)
                DO k = 1,num_coeffs
                    angle = (m(k)*u - n(k)*v)
                    sinangle = sin(angle)
                    cosangle = cos(angle)
              
              
                    !X derivatives
                    xu(i,j) = xu(i,j) - rmnc(k)*sinangle*cosphi*m(k)
             

                    xv(i,j) = xv(i,j) + rmnc(k)*(sinangle*cosphi*n(k) - &
                        cosangle*sinphi)
                    xuu(i,j) = xuu(i,j) - rmnc(k)*cosangle*cosphi*m(k)*m(k)
                    xuv(i,j) = xuv(i,j) + rmnc(k)*(cosangle*cosphi*m(k)*n(k) + &
                        sinangle*sinphi*m(k))
                    xvv(i,j) = xvv(i,j) - rmnc(k)*(2*sinangle*sinphi*n(k) + &
                        cosangle*cosphi*(n(k)*n(k) + 1))
                    !Y derivatives
                    yu(i,j) = yu(i,j) - rmnc(k)*sinangle*sinphi*m(k)
                    yv(i,j) = yv(i,j) + rmnc(k)*(sinangle*sinphi*n(k) + &
                        cosangle*cosphi)
                    yuu(i,j) = yuu(i,j) - rmnc(k)*cosangle*sinphi*m(k)*m(k)
                    yuv(i,j) = yuv(i,j) + rmnc(k)*(cosangle*sinphi*m(k)*n(k) - &
                        sinangle*cosphi*m(k))
                    yvv(i,j) = yvv(i,j) + rmnc(k)*(2*sinangle*cosphi*n(k) - &
                        cosangle*sinphi*(n(k)*n(k) + 1))
                    !Z derivatives
                    zu(i,j) = zu(i,j) + zmns(k)*cosangle*m(k)
                    zv(i,j) = zv(i,j) - zmns(k)*cosangle*n(k)
                    zuu(i,j) = zuu(i,j) - zmns(k)*sinangle*m(k)*m(k)
                    zuv(i,j) = zuv(i,j) + zmns(k)*sinangle*m(k)*n(k)
                    zvv(i,j) = zvv(i,j) - zmns(k)*sinangle*n(k)*n(k)



                end do
            end do
        end do

     
     
     
        !---------------Curvature calculation-----------------------!
        EE = xu*xu + yu*yu + zu*zu
        FF = xu*xv + yu*yv + zu*zv
        GG = xv*xv + yv*yv + zv*zv

        Nx = yu*zv - zu*yv
        Ny = zu*xv - xu*zv
        Nz = xu*yv - yu*xv
        magnorm = sqrt(Nx*Nx + Ny*Ny + Nz*Nz)
        Nx = Nx/magnorm
        Ny = Ny/magnorm
        Nz = Nz/magnorm

        LL = xuu*Nx + yuu*Ny + zuu*Nz
        MM = xuv*Nx + yuv*Ny + zuv*Nz
        NN = xvv*Nx + yvv*Ny + zvv*Nz

        HH = (EE*NN +GG*LL - 2*FF*MM)/(2*(EE*GG - FF*FF))
        KK = (LL*NN - MM*MM)/(EE*GG - FF*FF)
        P1 = HH + sqrt(HH*HH - KK)
        P2 = HH - sqrt(HH*HH - KK)

        mtargets = mtargets + 1
        targets(mtargets) = target
        sigmas(mtargets) = sigma
        vals(mtargets) = -1*minval(P2)
        IF (iflag == 1) WRITE(iunit_out,'(5ES22.12E3)') target,sigma,vals(mtargets),MINVAL(P1),MINVAL(P2)

      ELSE
        IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_curvature_P2
        END IF
      END IF
      RETURN
      END SUBROUTINE chisq_curvature_P2
