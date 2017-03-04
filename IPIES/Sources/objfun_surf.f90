!-----------------------------------------------------------------------
!     Subroutine:    objfun_surf
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/12/2011
!     Description:   This subroutine calculates the distance between a
!                    surface and a set of points defining that surface.
!-----------------------------------------------------------------------
      SUBROUTINE objfun_surf(mode,n,x,objf,objgrd,nstate,iuser,ruser)
!      SUBROUTINE func_surf(n,xc,fc,iuser,ruser)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_runtime
      USE pies_background, ONLY: xm, xn, mnmax
      USE pies_fieldlines, ONLY: nintw, thetaln, philn, Rln, Zln, ik_help
!-----------------------------------------------------------------------
!     Input Parameters
!          n           Number of Fourier coefficients
!          xc          Fourier coefficients (1:n)=(1:mnmax)
!          fc          Total distance to points of field-line
!          iuser(1)    NAG work array - Number of points in field-line
!          iuser(2)    NAG work array - Signs
!          iuser(3)    NAG work array - Flag to perform pre-calc (0=No,1=YES)
!          ruser       NAG work array - xu,xv,Fln
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER             :: n, mode, nstate
      INTEGER             :: iuser(*)
      DOUBLE PRECISION    :: x(n), objf, objgrd(n), ruser(*)
!-----------------------------------------------------------------------
!     Local Variables
!           nsteps     Number of points in field line following
!           ier        Error flag
!           fmn        Fourier coefficient helper array
!           xu         Theta Array
!           xv         
!-----------------------------------------------------------------------
      INTEGER :: ier, j, mn, ni
      INTEGER, SAVE :: nsteps
      REAL(rprec) ::  fmn(1:n), dfmn(1:n), temp, pert
      REAL(rprec), ALLOCATABLE, SAVE :: cosmt(:,:), sinmt(:,:), cosnz(:,:), sinnz(:,:)
      REAL(rprec), ALLOCATABLE, SAVE :: xu(:),xv(:), Fln(:), freal(:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Quick note, to make stellarator non-symmetric we would need to
      ! also pass the NSS coefficients in xc so that when we compared
      ! in real space we could compare the full deal.  Will also need
      ! modify signs to indicate R or Z instead of cos or sin
      IF (nstate == 1) THEN
         ! Deallocate if allocated
         IF (ALLOCATED(xu)) DEALLOCATE(xu)
         IF (ALLOCATED(xv)) DEALLOCATE(xv)
         IF (ALLOCATED(Fln)) DEALLOCATE(Fln)
         IF (ALLOCATED(freal)) DEALLOCATE(freal)
         IF (ALLOCATED(cosmt)) DEALLOCATE(cosmt)
         IF (ALLOCATED(sinmt)) DEALLOCATE(sinmt)
         IF (ALLOCATED(cosnz)) DEALLOCATE(cosnz)
         IF (ALLOCATED(sinnz)) DEALLOCATE(sinnz)
         PRINT *,'got here 1'
         ! Allocate and initialize arrays
         nsteps = iuser(1)
         !nsteps = nintw
         ALLOCATE(cosmt(1:mnmax,1:nintw),sinmt(1:mnmax,1:nintw),&
                  cosnz(1:mnmax,1:nintw),sinnz(1:mnmax,1:nintw),STAT=ier)
         ALLOCATE(xu(1:nsteps),xv(1:nsteps),Fln(1:nsteps),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'xu xv Fln (func_surf)',ier)
         PRINT *,'got here 2'
         !xu(:) = thetaln(ik_help,:)
         !xv(:) = philn(ik_help,:)
         xu(:) = ruser(1:nsteps)
         xv(:) = ruser(nsteps+1:2*nsteps)
         Fln(:) = ruser(2*nsteps+1:3*nsteps)
         ALLOCATE(freal(1:nsteps),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'freal (func_surf)',ier)
         PRINT *,'got here 3'
         IF (iuser(2) == 0) THEN
            Fln(:) = Rln(ik_help,:)
         ELSE
            Fln(:) = Zln(ik_help,:)
         END IF
         DO mn = 1, mnmax
            cosmt(mn,:) = cos(xm(mn)*xu(:))
            sinmt(mn,:) = sin(xm(mn)*xu(:))
            cosnz(mn,:) = cos(xn(mn)*xv(:))
            sinnz(mn,:) = sin(xn(mn)*xv(:))
         END DO
         PRINT *,'got here 4'
         RETURN
      END IF
      fmn(:) = x(:)
      freal = 0.0
      IF (iuser(2) == 0) THEN
         DO mn = 1, mnmax
            freal(:) = freal(:) +fmn(mn)*(cosmt(mn,:)*cosnz(mn,:)-sinmt(mn,:)*sinnz(mn,:))
         END DO
         objf = ABS(SUM(freal(:) - Fln(:)))
         PRINT *,'objf',objf
         IF (mode == 2) THEN
            DO ni = 1, mnmax
               freal = 0
               temp = fmn(ni)
               pert = MAX(1e-12,0.01 * fmn(ni))
               fmn(ni) = fmn(ni)+pert
               DO mn = 1, mnmax
                  freal(:) = freal(:) +fmn(mn)*(cosmt(mn,:)*cosnz(mn,:)-sinmt(mn,:)*sinnz(mn,:))
               END DO
               objgrd(ni) = SUM(freal)/pert
               fmn(ni) = temp
            END DO
            PRINT *,'obgrd',objgrd
         END IF
      ELSE
         DO mn = 1, mnmax
            freal(:) = freal(:) + fmn(mn)*(sinmt(mn,:)*cosnz(mn,:)+cosmt(mn,:)*sinnz(mn,:))
         END DO
         objf = ABS(SUM(freal(:) - Fln(:)))
         IF (mode == 2) THEN
            DO ni = 1, mnmax
               freal = 0
               temp = fmn(ni)
               fmn(ni) = fmn(ni)+0.0001
               DO mn = 1, mnmax
                  freal(:) = freal(:) +fmn(mn)*(sinmt(mn,:)*cosnz(mn,:)+cosmt(mn,:)*sinnz(mn,:))
               END DO
               objgrd(ni) = ABS(SUM(freal))/0.0001
               fmn(ni) = temp
            END DO
         END IF
      ENDIF
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE objfun_surf
