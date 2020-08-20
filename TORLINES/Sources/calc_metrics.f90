!-----------------------------------------------------------------------
!     Subroutine:    calc_metrics
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/15/2011
!     Description:   This subroutine calculates the metric elements from
!                    the values currently stored in the background
!                    and realspace grids.
!
!                    We define our fuctions to be of the form
!                    f(rho,theta,phi) = sum{n=-ntor,ntor} sum{m=0,mpol) fmn(s) * cos(m*theta+n*N*phi)
!                    So we must chain rule to u and v where
!                    u = theta/2*pi  and  v = phi / 2*pi*N
!                    Thus
!                    df/du = df/dtheta * dtheta/du
!                    df/dv = df/dphi   * dphi/dv
!-----------------------------------------------------------------------
      SUBROUTINE calc_metrics
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY : rprec
      USE torlines_background, xm => xm_sav, xn => xn_sav
      USE torlines_realspace
      USE torlines_runtime
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!          ier         Error flag
!          u           Poloidal dummy index
!          v           Toroidal dummy index
!          mn          Fourier dummy index
!          factor      ratio dphi/dv
!          fmn_temp    Dummy array to hold derivative coefficients (ex. rmnc*m)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, u, v, mn, s, mnmax
      REAL(rprec) :: factor
      REAL(rprec), ALLOCATABLE :: fmn_temp(:,:)
      TYPE(EZspline3_r8) :: r1_spl, z1_spl
      INTEGER :: bcs1(2),bcs2(2)
      REAL(rprec) :: R_grad(3), Z_grad(3)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      bcs1=(/0,0/)
      bcs2=(/-1,-1/)
      !mnmax = SIZE(rmnc_sav,1)
      ! Calculate deriviatives
      !ALLOCATE(fmn_temp(1:mnmax,1:k),STAT=ier)
      !IF (ier /= 0) CALL handle_err(ALLOC_ERR,'FMN_TEMP (calc_metrics)',ier)
      rs = 0; zs = 0; ru=0; zu=0; rv=0; zv=0; 
      !FORALL(mn = 1:mnmax) fmn_temp(mn,1:k) = -rmnc_sav(mn,1:k)*xm(mn)
      !CALL mntouv(1,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,ru,1,1)
      !FORALL(mn = 1:mnmax) fmn_temp(mn,1:k) = -rmnc_sav(mn,1:k)*xn(mn)*nfp
      !CALL mntouv(1,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,rv,1,0)  
      !FORALL(mn = 1:mnmax) fmn_temp(mn,1:k) = zmns_sav(mn,1:k)*xm(mn)
      !CALL mntouv(1,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,zu,0,0)  
      !FORALL(mn = 1:mnmax) fmn_temp(mn,1:k) = zmns_sav(mn,1:k)*xn(mn)*nfp
      !CALL mntouv(1,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,zv,0,0) 
      !IF (lasym) THEN
      !   FORALL(mn = 1:mnmax) fmn_temp(mn,1:k) = rmns_sav(mn,1:k)*xm(mn)
      !   CALL mntouv(1,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,ru,0,0)
      !   FORALL(mn = 1:mnmax) fmn_temp(mn,1:k) = rmns_sav(mn,1:k)*xn(mn)*nfp
      !   CALL mntouv(1,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,rv,0,0)  
      !   FORALL(mn = 1:mnmax) fmn_temp(mn,1:k) = -zmnc_sav(mn,1:k)*xm(mn)
      !   CALL mntouv(1,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,zu,1,0)  
      !   FORALL(mn = 1:mnmax) fmn_temp(mn,1:k) = -zmnc_sav(mn,1:k)*xn(mn)*nfp
      !   CALL mntouv(1,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,zv,1,0)
      !END IF
      !DEALLOCATE(fmn_temp)
         
      ! Spline formulation
      ! xv*2*pi*nfp (didn't work)
      ! xv*2*pi     (worse than before)
      ! removing 2*pi dependancy is horrible
      ! full grid (nu-1 -> nu and nv-1 -> nv, better)
      ! Do full grid (no filling)
      CALL EZspline_init(r1_spl,k,nu,nv,bcs1,bcs2,bcs2,ier)
      CALL EZspline_init(z1_spl,k,nu,nv,bcs1,bcs2,bcs2,ier)
      r1_spl%x1 = rho
      z1_spl%x1 = rho
      r1_spl%x2(1:nu) = xu(1:nu)*pi2
      z1_spl%x2(1:nu) = xu(1:nu)*pi2
      r1_spl%x3(1:nv) = xv(1:nv)*pi2
      z1_spl%x3(1:nv) = xv(1:nv)*pi2
      r1_spl%isHermite    = 1
      z1_spl%isHermite    = 1
      CALL EZspline_setup(r1_spl,rreal(1:k,1:nu,1:nv),ier,exact_dim=.true.)
      CALL EZspline_setup(z1_spl,zreal(1:k,1:nu,1:nv),ier,exact_dim=.true.)
      
      DO v = 1, nv
         DO u = 1, nu
            DO s = 1, k
               CALL EZspline_gradient(r1_spl,rho(s),xu(u)*pi2,xv(v)*pi2,R_grad,ier)
               CALL EZspline_gradient(z1_spl,rho(s),xu(u)*pi2,xv(v)*pi2,Z_grad,ier)
               !IF (lverb) WRITE(327,'(3I5,6ES20.10)') s,u,v,rs(s,u,v),ru(s,u,v),rv(s,u,v),R_grad
               rs(s,u,v) = R_grad(1)
               ru(s,u,v) = R_grad(2)
               rv(s,u,v) = R_grad(3)*nfp
               zs(s,u,v) = Z_grad(1)
               zu(s,u,v) = Z_grad(2)
               zv(s,u,v) = Z_grad(3)*nfp
            END DO
         END DO
      END DO
      
      CALL EZspline_free(r1_spl,ier)
      CALL EZspline_free(z1_spl,ier)
      
      ! Fill out array
      !rs(:,nu,1:nv) = rs(:,1,1:nv); rs(:,1:nu,nv) = rs(:,1:nu,1); rs(:,nu,nv) = rs(:,1,1)
      !ru(:,nu,1:nv) = ru(:,1,1:nv); ru(:,1:nu,nv) = ru(:,1:nu,1); ru(:,nu,nv) = ru(:,1,1)
      !rv(:,nu,1:nv) = rv(:,1,1:nv); rv(:,1:nu,nv) = rv(:,1:nu,1); rv(:,nu,nv) = rv(:,1,1)
      !zs(:,nu,1:nv) = zs(:,1,1:nv); zs(:,1:nu,nv) = zs(:,1:nu,1); zs(:,nu,nv) = zs(:,1,1)
      !zu(:,nu,1:nv) = zu(:,1,1:nv); zu(:,1:nu,nv) = zu(:,1:nu,1); zu(:,nu,nv) = zu(:,1,1)
      !zv(:,nu,1:nv) = zv(:,1,1:nv); zv(:,1:nu,nv) = zv(:,1:nu,1); zv(:,nu,nv) = zv(:,1,1)
      
      ! Calculate metric elements
      g11 = rs*rs+zs*zs
      g12 = rs*ru+zs*zu
      g13 = rs*rv+zs*zv
      g22 = ru*ru+zu*zu
      g23 = ru*rv+zu*zv
      g33 = rv*rv+zv*zv
      detg= g11*(g22*g33-g23*g23) + g12*(g23*g13-g12*g33) + g13*(g12*g23-g22*g13)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE calc_metrics
