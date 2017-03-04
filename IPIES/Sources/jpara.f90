!-----------------------------------------------------------------------
!     Subroutine:    jpara
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/14/2011
!     Description:   This subroutine calculates the parallel component
!                    of the total current density from the perpendicular
!                    component calculated in jperp.
!                    This is achieved by calculating
!                     div(j) = div(jxB/B) + div(jdotB) = 0
!-----------------------------------------------------------------------
      SUBROUTINE jpara
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_background
      USE pies_realspace
      USE pies_runtime
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier,u,v,mn
      REAL(rprec) :: a11, a12, a13, a23, a31, a32, a33
      REAL(rprec), ALLOCATABLE :: fmn_temp(:,:)
      REAL(rprec), ALLOCATABLE :: gmnc(:,:),gmns(:,:)
      REAL(rprec), ALLOCATABLE :: gs(:,:,:),gu(:,:,:),gv(:,:,:)
      REAL(rprec), ALLOCATABLE :: br(:,:,:),bphi(:,:,:),bz(:,:,:)
      REAL(rprec), ALLOCATABLE :: muu(:,:,:)
      REAL(rprec), ALLOCATABLE :: muumnc(:,:),muumns(:,:)
      REAL(rprec), ALLOCATABLE :: jparmnc(:,:),jparmns(:,:)
      REAL(rprec), ALLOCATABLE :: jss(:,:,:),juu(:,:,:),jvv(:,:,:)
      REAL(rprec), ALLOCATABLE :: divj(:,:,:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! First we Fourier Transform sqrt(detg)
      ALLOCATE(fmn_temp(1:mnmax,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'FMN_TEMP (jpara)',ier)
      ALLOCATE(gmnc(1:mnmax,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'GMNC (jpara)',ier)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,gmnc,xm,xn,sqrt(detg),0,1)
      IF (lasym) THEN
         ALLOCATE(gmns(1:mnmax,0:k),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'GMNS (jpara)',ier)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,gmns,xm,xn,sqrt(detg),1,0)
      END IF
      ! Now calculate derivative quantities
      ALLOCATE(gs(0:k,1:nu,1:nv),gu(0:k,1:nu,1:nv),gv(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'GS GU GV (jpara)',ier)
      ALLOCATE(jss(0:k,1:nu,1:nv),juu(0:k,1:nu,1:nv),jvv(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'JSS JUU JVV (jpara)',ier)
      DO u=1,nu
         DO v=1,nv
            gs(1:k-1,u,v) = 0.5*(sqrt(detg(2:k,u,v))-sqrt(detg(0:k-2,u,v)))
            jss(1:k-1,u,v) = 0.5*(jsreal(2:k,u,v)-jsreal(0:k-2,u,v))
            gs(0,u,v)     = sqrt(detg(1,u,v)) - sqrt(detg(0,u,v))
            jss(0,u,v)     = jsreal(1,u,v) - jsreal(0,u,v)
            gs(k,u,v)     = sqrt(detg(k,u,v)) - sqrt(detg(k-1,u,v))
            jss(k,u,v)     = jsreal(k,u,v) - jsreal(k-1,u,v)
         END DO
      END DO
      FORALL(mn=1:mnmax) fmn_temp(mn,0:k) = -gmnc(mn,0:k)*xm(mn)
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,gu,1,1)
      FORALL(mn=1:mnmax) fmn_temp(mn,0:k) = -gmnc(mn,0:k)*xn(mn)*nfp
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,gv,1,0)
      FORALL(mn=1:mnmax) fmn_temp(mn,0:k) = -jumnc(mn,0:k)*xm(mn)
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,juu,1,0)
      FORALL(mn=1:mnmax) fmn_temp(mn,0:k) = -jvmnc(mn,0:k)*xn(mn)*nfp
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,jvv,1,0)
      IF (lasym) THEN
         FORALL(mn=1:mnmax) fmn_temp(mn,0:k) = gmns(mn,0:k)*xm(mn)
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,gu,0,0)
         FORALL(mn=1:mnmax) fmn_temp(mn,0:k) = gmns(mn,0:k)*xn(mn)*nfp
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,gv,0,0)
         FORALL(mn=1:mnmax) fmn_temp(mn,0:k) = jumns(mn,0:k)*xm(mn)
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,juu,0,0)
         FORALL(mn=1:mnmax) fmn_temp(mn,0:k) = jvmns(mn,0:k)*xn(mn)*nfp
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,jvv,0,0)
      END IF
      DEALLOCATE(fmn_temp)
      ! Now we Calculate the div(j_perp)
      ! note g=sqrt(g) in the follwing explanation
      !      js  = j^s
      !      gs  = d[sqrt(g)]/ds etc.
      !      jss = d[js]/ ds    etc.
      ! div(j_perp) = (gs*js+g*jss+gu*ju+g*juu+gv*jv+g*jvv)/g
      ALLOCATE(divj(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'DIVJ (jpara)',ier)
      divj = (gs*jsreal+gu*jureal+gv*jvreal+sqrt(detg)*(jss+juu+jvv))/sqrt(detg)
      PRINT *,' '
      PRINT *,' divj(:,1,1)'
      PRINT *,divj(:,1,1)
      ! Now calculated derivatives
      
      ! Now we calculate the parallel current (dmu/du)
      ALLOCATE(muu(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'MUU (jpara)',ier)
      ALLOCATE(br(0:k,1:nu,1:nv),bphi(0:k,1:nu,1:nv),bz(0:k,1:nu,1:nv), STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BR BPHI BZ (jperp)',ier)
      br   = bsreal*rs+bureal*ru+bvreal*rv
      bphi = rreal*bvreal
      bz   = bsreal*zs+bureal*zu+bvreal*zv
      muu = divj*( 1/bphi - (zs/br-rs/bz)/(rs*zu+zs*ru) )
      DEALLOCATE(br,bphi,bz)
      ! Now we Fourier transform dmu/du
      ALLOCATE(muumnc(1:mnmax,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'MUUMNC (jpara)',ier)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,muumnc,xm,xn,muu,1,1)
      jparmnc = muumnc
      DO mn = 1, mnmax
         IF (xm(mn) /=0) jparmnc(mn,:) = muumnc(mn,:)/xm(mn)
      END DO
      IF (lasym) THEN
         ALLOCATE(muumns(1:mnmax,0:k),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'MUUMNC (jpara)',ier)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,muumns,xm,xn,muu,0,0)
      END IF
      
      
      
      ! DEALLOCATIONS
      DEALLOCATE(gmnc)
      IF (lasym) THEN
         DEALLOCATE(gmns)
      END IF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE jpara
