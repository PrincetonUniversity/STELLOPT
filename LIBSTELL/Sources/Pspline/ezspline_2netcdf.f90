!/////
! R8 !
!/////
   subroutine EZspline_2NetCDF_array3_r8(n1, n2, n3, x1, x2, x3, f, filename, ier)
     use ezspline_obj   
     use ezcdf
     implicit none
     integer, intent(in) :: n1, n2, n3
     real(ezspline_r8), intent(in) :: x1(:), x2(:), x3(:), f(:, :, :)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
     integer ifail

     integer dimlens(3), ncid
     
     ier = 0
     call cdfOpn(ncid, filename, 'w')
     if(ncid==0) then
        ier = 34
        return
     endif
     dimlens = (/ 1, 1, 1 /)
     call cdfDefVar(ncid, 'n1', dimlens, 'INT', ifail)
     call cdfDefVar(ncid, 'n2', dimlens, 'INT', ifail)
     call cdfDefVar(ncid, 'n3', dimlens, 'INT', ifail)
     dimlens = (/ n1, 1, 1 /)
     call cdfDefVar(ncid, 'x1', dimlens, 'R8', ifail)
     dimlens = (/ n2, 1, 1 /)
     call cdfDefVar(ncid, 'x2', dimlens, 'R8', ifail)
     dimlens = (/ n3, 1, 1 /)
     call cdfDefVar(ncid, 'x3', dimlens, 'R8', ifail)
     dimlens = (/ n1, n2, n3 /)
     call cdfDefVar(ncid, 'f', dimlens, 'R8', ifail)
     if(ifail/=0) then
        ier = 38 
        return
     endif
     call cdfPutVar(ncid, 'n1', n1, ifail)
     call cdfPutVar(ncid, 'n2', n2, ifail)
     call cdfPutVar(ncid, 'n3', n3, ifail)
     call cdfPutVar(ncid, 'x1', x1, ifail)
     call cdfPutVar(ncid, 'x2', x2, ifail)
     call cdfPutVar(ncid, 'x3', x3, ifail)
     call cdfPutVar(ncid, 'f', f, ifail)
     if(ifail/=0) then
        ier = 35 
        return
     endif

     call cdfCls(ncid)

   end subroutine EZspline_2NetCDF_array3_r8

   subroutine EZspline_2NetCDF_array2_r8(n1, n2, x1, x2, f, filename, ier)
     use ezspline_obj   
     use ezcdf
     implicit none
     integer, intent(in) :: n1, n2
     real(ezspline_r8), intent(in) :: x1(:), x2(:), f(:,:)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
     integer ifail

     integer dimlens(3), ncid
     
     ier = 0
     call cdfOpn(ncid, filename, 'w')
     if(ncid==0) then
        ier = 34
        return
     endif
     dimlens = (/ 1, 1, 1 /)
     call cdfDefVar(ncid, 'n1', dimlens, 'INT', ifail)
     call cdfDefVar(ncid, 'n2', dimlens, 'INT', ifail)
     dimlens = (/ n1, 1, 1 /)
     call cdfDefVar(ncid, 'x1', dimlens, 'R8', ifail)
     dimlens = (/ n2, 1, 1 /)
     call cdfDefVar(ncid, 'x2', dimlens, 'R8', ifail)
     dimlens = (/ n1, n2, 1 /)
     call cdfDefVar(ncid, 'f', dimlens, 'R8', ifail)
     if(ifail/=0) then
        ier = 38 
        return
     endif

     call cdfPutVar(ncid, 'n1', n1, ifail)
     call cdfPutVar(ncid, 'n2', n2, ifail)
     call cdfPutVar(ncid, 'x1', x1, ifail)
     call cdfPutVar(ncid, 'x2', x2, ifail)
     call cdfPutVar(ncid, 'f', f, ifail)
     if(ifail/=0) then
        ier = 35 
        return
     endif

     call cdfCls(ncid)

   end subroutine EZspline_2NetCDF_array2_r8

   subroutine EZspline_2NetCDF1_r8(n1, x1, f, filename, ier)
     use ezspline_obj   
     use ezcdf
     implicit none
     integer, intent(in) :: n1
     real(ezspline_r8), intent(in) :: x1(:), f(:)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
     integer ifail

     integer dimlens(3), ncid
     
     ier = 0
     call cdfOpn(ncid, filename, 'w')
     if(ncid==0) then
        ier = 34
        return
     endif
     dimlens = (/ 1, 1, 1 /)
     call cdfDefVar(ncid, 'n1', dimlens, 'INT', ifail)
     dimlens = (/ n1, 1, 1 /)
     call cdfDefVar(ncid, 'x1', dimlens, 'R8', ifail)
     dimlens = (/ n1, 1, 1 /)
     call cdfDefVar(ncid, 'f', dimlens, 'R8', ifail)
     if(ifail/=0) then
        ier = 38 
        return
     endif

     call cdfPutVar(ncid, 'n1', n1, ifail)
     call cdfPutVar(ncid, 'x1', x1, ifail)
     call cdfPutVar(ncid, 'f', f, ifail)
     if(ifail/=0) then
        ier = 35 
        return
     endif

     call cdfCls(ncid)

   end subroutine EZspline_2NetCDF1_r8

   subroutine EZspline_2NetCDF_cloud3_r8(n, x1, x2, x3, f, filename, ier)
     use ezspline_obj   
     use ezcdf
     implicit none
     integer, intent(in) :: n
     real(ezspline_r8), intent(in) :: x1(:), x2(:), x3(:), f(:)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
     integer ifail

     integer dimlens(3), ncid
     
     ier = 0
     call cdfOpn(ncid, filename, 'w')
     if(ncid==0) then
        ier = 34
        return
     endif
     dimlens = (/ 1, 1, 1 /)
     call cdfDefVar(ncid, 'n', dimlens, 'INT', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'x1', dimlens, 'R8', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'x2', dimlens, 'R8', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'x3', dimlens, 'R8', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'f', dimlens, 'R8', ifail)
     if(ifail/=0) then
        ier = 38 
        return
     endif

     call cdfPutVar(ncid, 'n', n, ifail)
     call cdfPutVar(ncid, 'x1', x1, ifail)
     call cdfPutVar(ncid, 'x2', x2, ifail)
     call cdfPutVar(ncid, 'x3', x3, ifail)
     call cdfPutVar(ncid, 'f', f, ifail)
     if(ifail/=0) then
        ier = 35 
        return
     endif

     call cdfCls(ncid)

   end subroutine EZspline_2NetCDF_cloud3_r8

   subroutine EZspline_2NetCDF_cloud2_r8(n, x1, x2, f, filename, ier)
     use ezspline_obj   
     use ezcdf
     implicit none
     integer, intent(in) :: n
     real(ezspline_r8), intent(in) :: x1(:), x2(:), f(:)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
     integer ifail

     integer dimlens(3), ncid
     
     ier = 0
     call cdfOpn(ncid, filename, 'w')
     if(ncid==0) then
        ier = 34
        return
     endif
     dimlens = (/ 1, 1, 1 /)
     call cdfDefVar(ncid, 'n', dimlens, 'INT', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'x1', dimlens, 'R8', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'x2', dimlens, 'R8', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'f', dimlens, 'R8', ifail)
     if(ifail/=0) then
        ier = 38 
        return
     endif

     call cdfPutVar(ncid, 'n', n, ifail)
     call cdfPutVar(ncid, 'x1', x1, ifail)
     call cdfPutVar(ncid, 'x2', x2, ifail)
     call cdfPutVar(ncid, 'f', f, ifail)
     if(ifail/=0) then
        ier = 35 
        return
     endif

     call cdfCls(ncid)

   end subroutine EZspline_2NetCDF_cloud2_r8
!/////
! R4 !
!/////
   subroutine EZspline_2NetCDF_array3_r4(n1, n2, n3, x1, x2, x3, f, filename, ier)
     use ezspline_obj   
     use ezcdf
     implicit none
     integer, intent(in) :: n1, n2, n3
     real(ezspline_r4), intent(in) :: x1(:), x2(:), x3(:), f(:, :, :)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
     integer ifail

     integer dimlens(3), ncid
     
     ier = 0
     call cdfOpn(ncid, filename, 'w')
     if(ncid==0) then
        ier = 34
        return
     endif
     dimlens = (/ 1, 1, 1 /)
     call cdfDefVar(ncid, 'n1', dimlens, 'INT', ifail)
     call cdfDefVar(ncid, 'n2', dimlens, 'INT', ifail)
     call cdfDefVar(ncid, 'n3', dimlens, 'INT', ifail)
     dimlens = (/ n1, 1, 1 /)
     call cdfDefVar(ncid, 'x1', dimlens, 'R4', ifail)
     dimlens = (/ n2, 1, 1 /)
     call cdfDefVar(ncid, 'x2', dimlens, 'R4', ifail)
     dimlens = (/ n3, 1, 1 /)
     call cdfDefVar(ncid, 'x3', dimlens, 'R4', ifail)
     dimlens = (/ n1, n2, n3 /)
     call cdfDefVar(ncid, 'f', dimlens, 'R4', ifail)
     if(ifail/=0) then
        ier = 38 
        return
     endif
     call cdfPutVar(ncid, 'n1', n1, ifail)
     call cdfPutVar(ncid, 'n2', n2, ifail)
     call cdfPutVar(ncid, 'n3', n3, ifail)
     call cdfPutVar(ncid, 'x1', x1, ifail)
     call cdfPutVar(ncid, 'x2', x2, ifail)
     call cdfPutVar(ncid, 'x3', x3, ifail)
     call cdfPutVar(ncid, 'f', f, ifail)
     if(ifail/=0) then
        ier = 35 
        return
     endif

     call cdfCls(ncid)

   end subroutine EZspline_2NetCDF_array3_r4

   subroutine EZspline_2NetCDF_array2_r4(n1, n2, x1, x2, f, filename, ier)
     use ezspline_obj   
     use ezcdf
     implicit none
     integer, intent(in) :: n1, n2
     real(ezspline_r4), intent(in) :: x1(:), x2(:), f(:,:)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
     integer ifail

     integer dimlens(3), ncid
     
     ier = 0
     call cdfOpn(ncid, filename, 'w')
     if(ncid==0) then
        ier = 34
        return
     endif
     dimlens = (/ 1, 1, 1 /)
     call cdfDefVar(ncid, 'n1', dimlens, 'INT', ifail)
     call cdfDefVar(ncid, 'n2', dimlens, 'INT', ifail)
     dimlens = (/ n1, 1, 1 /)
     call cdfDefVar(ncid, 'x1', dimlens, 'R4', ifail)
     dimlens = (/ n2, 1, 1 /)
     call cdfDefVar(ncid, 'x2', dimlens, 'R4', ifail)
     dimlens = (/ n1, n2, 1 /)
     call cdfDefVar(ncid, 'f', dimlens, 'R4', ifail)
     if(ifail/=0) then
        ier = 38 
        return
     endif

     call cdfPutVar(ncid, 'n1', n1, ifail)
     call cdfPutVar(ncid, 'n2', n2, ifail)
     call cdfPutVar(ncid, 'x1', x1, ifail)
     call cdfPutVar(ncid, 'x2', x2, ifail)
     call cdfPutVar(ncid, 'f', f, ifail)
     if(ifail/=0) then
        ier = 35 
        return
     endif

     call cdfCls(ncid)

   end subroutine EZspline_2NetCDF_array2_r4

   subroutine EZspline_2NetCDF1_r4(n1, x1, f, filename, ier)
     use ezspline_obj   
     use ezcdf
     implicit none
     integer, intent(in) :: n1
     real(ezspline_r4), intent(in) :: x1(:), f(:)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
     integer ifail

     integer dimlens(3), ncid
     
     ier = 0
     call cdfOpn(ncid, filename, 'w')
     if(ncid==0) then
        ier = 34
        return
     endif
     dimlens = (/ 1, 1, 1 /)
     call cdfDefVar(ncid, 'n1', dimlens, 'INT', ifail)
     dimlens = (/ n1, 1, 1 /)
     call cdfDefVar(ncid, 'x1', dimlens, 'R4', ifail)
     dimlens = (/ n1, 1, 1 /)
     call cdfDefVar(ncid, 'f', dimlens, 'R4', ifail)
     if(ifail/=0) then
        ier = 38 
        return
     endif

     call cdfPutVar(ncid, 'n1', n1, ifail)
     call cdfPutVar(ncid, 'x1', x1, ifail)
     call cdfPutVar(ncid, 'f', f, ifail)
     if(ifail/=0) then
        ier = 35 
        return
     endif

     call cdfCls(ncid)

   end subroutine EZspline_2NetCDF1_r4

   subroutine EZspline_2NetCDF_cloud3_r4(n, x1, x2, x3, f, filename, ier)
     use ezspline_obj   
     use ezcdf
     implicit none
     integer, intent(in) :: n
     real(ezspline_r4), intent(in) :: x1(:), x2(:), x3(:), f(:)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
     integer ifail

     integer dimlens(3), ncid
     
     ier = 0
     call cdfOpn(ncid, filename, 'w')
     if(ncid==0) then
        ier = 34
        return
     endif
     dimlens = (/ 1, 1, 1 /)
     call cdfDefVar(ncid, 'n', dimlens, 'INT', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'x1', dimlens, 'R4', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'x2', dimlens, 'R4', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'x3', dimlens, 'R4', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'f', dimlens, 'R4', ifail)
     if(ifail/=0) then
        ier = 38 
        return
     endif

     call cdfPutVar(ncid, 'n', n, ifail)
     call cdfPutVar(ncid, 'x1', x1, ifail)
     call cdfPutVar(ncid, 'x2', x2, ifail)
     call cdfPutVar(ncid, 'x3', x3, ifail)
     call cdfPutVar(ncid, 'f', f, ifail)
     if(ifail/=0) then
        ier = 35 
        return
     endif

     call cdfCls(ncid)

   end subroutine EZspline_2NetCDF_cloud3_r4

   subroutine EZspline_2NetCDF_cloud2_r4(n, x1, x2, f, filename, ier)
     use ezspline_obj   
     use ezcdf
     implicit none
     integer, intent(in) :: n
     real(ezspline_r4), intent(in) :: x1(:), x2(:), f(:)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
     integer ifail

     integer dimlens(3), ncid
     
     ier = 0
     call cdfOpn(ncid, filename, 'w')
     if(ncid==0) then
        ier = 34
        return
     endif
     dimlens = (/ 1, 1, 1 /)
     call cdfDefVar(ncid, 'n', dimlens, 'INT', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'x1', dimlens, 'R4', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'x2', dimlens, 'R4', ifail)
     dimlens = (/ n, 1, 1 /)
     call cdfDefVar(ncid, 'f', dimlens, 'R4', ifail)
     if(ifail/=0) then
        ier = 38 
        return
     endif

     call cdfPutVar(ncid, 'n', n, ifail)
     call cdfPutVar(ncid, 'x1', x1, ifail)
     call cdfPutVar(ncid, 'x2', x2, ifail)
     call cdfPutVar(ncid, 'f', f, ifail)
     if(ifail/=0) then
        ier = 35 
        return
     endif

     call cdfCls(ncid)

   end subroutine EZspline_2NetCDF_cloud2_r4
