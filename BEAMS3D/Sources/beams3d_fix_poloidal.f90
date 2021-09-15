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
      INTEGER :: mystart, mystep, win_X4D, win_Y4D, ier
      INTEGER :: bcs1(2), bcs2(2), bcs3(2)
      TYPE(EZspline3_r8) :: X_spl, Y_spl
      REAL(rprec), POINTER, DIMENSION(:,:,:) :: X_ARR, Y_ARR
      REAL(rprec), POINTER, DIMENSION(:,:,:,:) :: X4D, Y4D
      ! For splines
      INTEGER :: i,j,k,l
      REAL*8 :: xparam, yparam, zparam !, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(1), r1, p1, z1, xt, yt
      INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Assume each processor will handle it's own Lines
      CALL mpialloc(X4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_X4D)
      CALL mpialloc(Y4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_Y4D)
      IF (myid_sharmem == 0) THEN
         ALLOCATE(X_ARR(nr,nphi,nz),Y_ARR(nr,nphi,nz))
         X_ARR = COS(U4D(1,:,:,:))
         Y_ARR = SIN(U4D(1,:,:,:))
         bcs1=(/ 0, 0/)
         bcs2=(/-1,-1/)
         bcs3=(/ 0, 0/)
         ier = 0
         CALL EZspline_init(X_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_fix_poloidal:X_spl',ier)
         CALL EZspline_init(Y_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_fix_poloidal:Y_spl',ier)
         X_spl%isHermite   = 1
         X_spl%x1   = raxis
         X_spl%x2   = phiaxis
         X_spl%x3   = zaxis
         Y_spl%isHermite   = 1
         Y_spl%x1   = raxis
         Y_spl%x2   = phiaxis
         Y_spl%x3   = zaxis
         CALL EZspline_setup(X_spl,X_ARR,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_fix_poloidal:X_spl',ier)
         CALL EZspline_setup(Y_spl,Y_ARR,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_fix_poloidal:Y_spl',ier)
         X4D = X_spl%fspl
         Y4D = Y_spl%fspl
         CALL EZspline_free(X_spl,ier)
         CALL EZspline_free(Y_spl,ier)
         DEALLOCATE(X_ARR,Y_ARR)
      END IF

      ! Don't do work till X4D and Y4D exist
      CALL MPI_BARRIER(MPI_COMM_SHARMEM, ier)
      
      ! Recalculate U for all particles
      mystart = LBOUND(R_lines,2)
      DO myline = mystart,myend
         DO mystep = 0, npoinc
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
            xt = fval(1)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                            Y4D(1,1,1,1),nr,nphi,nz)
            yt = fval(1)
            U_lines(mystep,myline) = ATAN2(yt,xt)
         END DO
      END DO
      WHERE ( R_lines < 0 ) U_lines = 0

      ! Wait to deallocate till everyone is done with the memory
      CALL MPI_BARRIER(MPI_COMM_SHARMEM, ier)

      ! DEALLOCATX X and Y
      CALL mpidealloc(X4D,win_X4D)
      CALL mpidealloc(Y4D,win_Y4D)

      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_fix_poloidal