c---------------------------------------------------------
c  error check routine for 2d hybrid interpolation
c
      subroutine vecin3d_argchk(subname,jspline,
     >     icoeff,nx,ny,nz,ixdim,iydim,izdim,ier)

      implicit NONE

      character*(*), intent(in) :: subname   ! name of caller, for messages
      integer, intent(in) :: jspline(3)      ! interp. type, by dimension
      integer, intent(in) :: icoeff          ! coeffs per data point
      integer, intent(in) :: nx,ny,nz        ! grid sizes
      integer, intent(in) :: ixdim,iydim,izdim     ! coeff array dimensions

      integer, intent(out) :: ier            ! status coder (0=OK)

      !----------------------------------------
      !  local
      integer :: i,itest,jcoeff,ipcub
      !----------------------------------------

      ier=0

      jcoeff=1
      ipcub=-99

      do i=1,3
         if((jspline(i).lt.-1).or.(jspline(i).gt.2)) then
            write(6,*) ' ?'//trim(subname)//
     >           ': jspline(',i,')=',jspline(i)
            write(6,*) '  not in expected range -1:2.'
            ier=1

         else if(jspline(i).ge.1) then
            jcoeff=jcoeff*2
            if(ipcub.eq.-99) then
               ipcub=jspline(i)
            else
               if(ipcub.ne.jspline(i)) then
                  write(6,*) ' ?'//trim(subname)//
     >                 ': more than one cubic method: jspline=',jspline
               endif
            endif

         endif
      enddo

      if(ier.eq.0) then
         if(jcoeff.ne.icoeff) then
            write(6,*) ' ?'//trim(subname)//
     >           ': incorrect #coefficients per data point.'
            write(6,*) '  expected: ',jcoeff,'  received: ',icoeff
            write(6,*) '  jspline: ',jspline
            ier=1
         endif

         itest=nx
         if(jspline(1).eq.-1) itest=nx-1
         if(ixdim.ne.itest) then
            write(6,*) ' ?'//trim(subname)//': ixdim dimension error:'
            write(6,*) '  nx=',nx,' ixdim=',ixdim,' expected: ',itest
            write(6,*) '  jspline(1)=',jspline(1),'; -1 for zonal.'
            ier=1 
         endif

         itest=ny
         if(jspline(2).eq.-1) itest=ny-1
         if(iydim.ne.itest) then
            write(6,*) ' ?'//trim(subname)//': iydim dimension error:'
            write(6,*) '  ny=',ny,' iydim=',iydim,' expected: ',itest
            write(6,*) '  jspline(2)=',jspline(2),'; -1 for zonal.'
            ier=1
         endif

         itest=nz
         if(jspline(1).eq.-1) itest=nz-1
         if(izdim.ne.itest) then
            write(6,*) ' ?'//trim(subname)//': izdim dimension error:'
            write(6,*) '  nz=',nz,' izdim=',izdim,' expected: ',itest
            write(6,*) '  jspline(3)=',jspline(3),'; -1 for zonal.'
            ier=1 
         endif

      endif

      return
      end
