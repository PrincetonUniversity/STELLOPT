
      module  bnvariables
      use stel_kinds
      real(rprec) ::  pi,pi2,alp,alu,alv,alvp,fnuv,curpol
      real(rprec), allocatable :: sinu(:,:), conu(:,:), sinv(:,:),
     1    conv(:,:)
      real(rprec), allocatable :: tanu(:), tanv(:)
      integer, allocatable :: indu(:), indv(:)
      real(rprec), allocatable ::   x(:), y(:), z(:), djx(:),
     1    djy(:), djz(:)
      real(rprec), allocatable ::  xu(:), yu(:), zu(:), xv(:),
     1    yv(:), zv(:)
      real(rprec), allocatable :: guu(:), guv(:), gvv(:), dju(:),
     1    djv(:), sqf(:)
      real(rprec), allocatable ::  au(:), av(:)
      real(rprec), allocatable :: bsubu(:,:), bsubv(:,:), bu(:),
     1    bv(:)
      real(rprec), allocatable :: cr(:,:), cz(:,:), cl(:, :)
      real(rprec), allocatable :: bsubus(:,:), bsubvs(:,:)
      real(rprec), allocatable :: crs(:,:), czc(:,:), clc(:, :)

      end module bnvariables
