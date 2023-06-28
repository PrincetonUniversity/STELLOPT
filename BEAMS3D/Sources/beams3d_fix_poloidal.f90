!-----------------------------------------------------------------------
!     Module:        beams3d_fix_poloidal
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          05/20/2021
!     Description:   This subroutine corrects U_lines using and X/Y
!                    formulation
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_fix_poloidal
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE beams3d_grid
      USE beams3d_lines
      USE beams3d_runtime
      USE mpi_inc
      USE mpi_sharmem
      USE EZspline_obj
      USE EZspline
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input Variables
!        NONE
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      INTEGER :: mystart, mystep,  ier !win_X4D, win_Y4D,
      INTEGER :: bcs1(2), bcs2(2), bcs3(2)
      !TYPE(EZspline3_r8) :: X_spl, Y_spl
      !REAL(rprec), POINTER, DIMENSION(:,:,:) :: X_ARR, Y_ARR
      !REAL(rprec), POINTER, DIMENSION(:,:,:,:) :: X4D, Y4D
      
      ! For splines
      INTEGER :: i,j,k,l
      REAL*8 :: xparam, yparam, zparam !, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(1), r1, p1, z1, xt, yt, s1
      INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      
      ! Recalculate U for all particles
      mystart = LBOUND(R_lines,2)
      DO myline = mystart,myend
         DO mystep = 0, npoinc
            s1 = S_lines(mystep,myline)
            r1 = R_lines(mystep,myline)
            p1 = MOD(PHI_lines(mystep,myline),phimax)
            IF (p1 < 0) p1 = p1 + phimax
            z1 = Z_lines(mystep,myline)
            i = MIN(MAX(COUNT(raxis < r1),1),nr-1)
            j = MIN(MAX(COUNT(phiaxis < p1),1),nphi-1)
            k = MIN(MAX(COUNT(zaxis < z1),1),nz-1)
            xparam = (r1 - raxis(i)) * hri(i)
            yparam = (p1 - phiaxis(j)) * hpi(j)
            zparam = (z1 - zaxis(k)) * hzi(k)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            X4D(1,1,1,1),nr,nphi,nz)
            xt = fval(1)/s1
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            Y4D(1,1,1,1),nr,nphi,nz)
            yt = fval(1)/s1
            U_lines(mystep,myline) = ATAN2(yt,xt)
         END DO
      END DO
      WHERE ( R_lines < 0 ) U_lines = 0

      ! Wait to deallocate till everyone is done with the memory
      CALL MPI_BARRIER(MPI_COMM_SHARMEM, ier)

      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_fix_poloidal