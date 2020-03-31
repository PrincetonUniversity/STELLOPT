!-----------------------------------------------------------------------
!     Subroutine:    cyl2suv
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/15/2011
!     Description:   This subroutine takes a set of vectors in
!                    cylindrical (R,phi,Z) coordinates and calculates
!                    their curvilinear (s,u,v) components.
!-----------------------------------------------------------------------
      SUBROUTINE cyl2suv(k1,k2,nu,nv,r,ar,aphi,az,as,au,av,factor)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY : rprec
      USE pies_realspace, ONLY : rs, zs, ru, zu, rv, zv
!-----------------------------------------------------------------------
!     Input Parameters
!          k1           Starting index of surfaces
!          k2           Ending index of surfaces
!          nu           Number of Poloidal points
!          nv           Number of toroidal points
!          r            R-Vector
!          ar           R-hat component of vector
!          aphi         Phi-hat component of vector
!          az           Z-hat component of vector
!          as           s component of vector
!          au           u component of vector
!          av           v component of vector   
!          factor       dphi/dv scaling factor (optional)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in)     :: k1
      INTEGER, INTENT(in)     :: k2
      INTEGER, INTENT(in)     :: nu
      INTEGER, INTENT(in)     :: nv
      REAL(rprec), INTENT(in) :: r(k1:k2,1:nu,1:nv)
      REAL(rprec), INTENT(in) :: ar(k1:k2,1:nu,1:nv)
      REAL(rprec), INTENT(in) :: aphi(k1:k2,1:nu,1:nv)
      REAL(rprec), INTENT(in) :: az(k1:k2,1:nu,1:nv)
      REAL(rprec), INTENT(out) :: as(k1:k2,1:nu,1:nv)
      REAL(rprec), INTENT(out) :: au(k1:k2,1:nu,1:nv)
      REAL(rprec), INTENT(out) :: av(k1:k2,1:nu,1:nv)
      REAL(rprec), INTENT(in)  :: factor
!-----------------------------------------------------------------------
!     Local Variables
!          detinv       Inverse of the determinant
!-----------------------------------------------------------------------
      REAL(rprec) ::  detinv(k1:k2,1:nu,1:nv)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      av(k1:k2,:,:)     = aphi(k1:k2,:,:)/r(k1:k2,:,:)
      av = av / factor
      detinv(k1:k2,:,:) = 1./(rs(k1:k2,:,:)*zu(k1:k2,:,:)-zs(k1:k2,:,:)*ru(k1:k2,:,:))
      as(k1:k2,:,:)     = detinv(k1:k2,:,:)*( zu(k1:k2,:,:)*(ar(k1:k2,:,:)-av(k1:k2,:,:)*rv(k1:k2,:,:))&
                                             -ru(k1:k2,:,:)*(az(k1:k2,:,:)-av(k1:k2,:,:)*zv(k1:k2,:,:)))
      au(k1:k2,:,:)     = detinv(k1:k2,:,:)*(-zs(k1:k2,:,:)*(ar(k1:k2,:,:)-av(k1:k2,:,:)*rv(k1:k2,:,:))&
                                             +rs(k1:k2,:,:)*(az(k1:k2,:,:)-av(k1:k2,:,:)*zv(k1:k2,:,:)))

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE cyl2suv
