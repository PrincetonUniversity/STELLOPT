

      subroutine bn_alloc
      use meshes
      use bnvariables
      implicit none
      integer ierr

      IF (ALLOCATED(bu)) DEALLOCATE(bu)
      IF (ALLOCATED(bv)) DEALLOCATE(bv)
      IF (ALLOCATED(sinu)) DEALLOCATE(sinu)
      IF (ALLOCATED(conu)) DEALLOCATE(conu)
      IF (ALLOCATED(sinv)) DEALLOCATE(sinv)
      IF (ALLOCATED(conv)) DEALLOCATE(conv)
      IF (ALLOCATED(tanu)) DEALLOCATE(tanu)
      IF (ALLOCATED(tanv)) DEALLOCATE(tanv)
      IF (ALLOCATED(indu)) DEALLOCATE(indu)
      IF (ALLOCATED(indv)) DEALLOCATE(indv)
      IF (ALLOCATED(x)) DEALLOCATE(x)
      IF (ALLOCATED(y)) DEALLOCATE(y)
      IF (ALLOCATED(z)) DEALLOCATE(z)
      IF (ALLOCATED(djx)) DEALLOCATE(djx)
      IF (ALLOCATED(djy)) DEALLOCATE(djy)
      IF (ALLOCATED(djz)) DEALLOCATE(djz)
      IF (ALLOCATED(xu)) DEALLOCATE(xu)
      IF (ALLOCATED(yu)) DEALLOCATE(yu)
      IF (ALLOCATED(zu)) DEALLOCATE(zu)
      IF (ALLOCATED(xv)) DEALLOCATE(xv)
      IF (ALLOCATED(yv)) DEALLOCATE(yv)
      IF (ALLOCATED(zv)) DEALLOCATE(zv)
      IF (ALLOCATED(guu)) DEALLOCATE(guu)
      IF (ALLOCATED(guv)) DEALLOCATE(guv)
      IF (ALLOCATED(gvv)) DEALLOCATE(gvv)
      IF (ALLOCATED(dju)) DEALLOCATE(dju)
      IF (ALLOCATED(djv)) DEALLOCATE(djv)
      IF (ALLOCATED(sqf)) DEALLOCATE(sqf)
      IF (ALLOCATED(au)) DEALLOCATE(au)
      IF (ALLOCATED(av)) DEALLOCATE(av)
      allocate ( bu(nuv),                  bv(nuv)
     .          ,sinu(nuv,  0:md),conu(nuv,  0:md)
     .          ,sinv(nuv,-nd:nd),conv(nuv,-nd:nd)
     .          ,tanu(-nu:nu),      tanv(-nvp:nvp)
     .          ,indu(nuvp),            indv(nuvp)
     .          ,    x(nuvp),   y(nuvp),   z(nuvp)
     .          ,  djx(nuvp), djy(nuvp), djz(nuvp)
     .          ,    xu(nuv),   yu(nuv),   zu(nuv)

     .          ,    xv(nuv),   yv(nuv),   zv(nuv)
     .          , guu(nuv),guv(nuv),gvv(nuv)
     .          , dju(nuv),djv(nuv),sqf(nuv)
     .          , au(nuv) , av(nuv),               stat=ierr)

      if (ierr .ne. 0) stop 'Allocation in BNORM subroutine bn_alloc'

      end subroutine bn_alloc
