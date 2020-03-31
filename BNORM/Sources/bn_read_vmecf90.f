
      subroutine bn_read_vmecf90(extension)
      use  meshes
      use  bnvariables
      use neswrite, only: coil_separation, mnmax_in => mnmax, ixm, ixn,
     1   raxis_in => raxis, zaxis_in => zaxis, nfp_in => nfp,
     2   iota_edge, phip_edge, raxis_s, zaxis_c, lasym_bn
      use read_wout_mod                  !Error-free reading of wout file
      implicit none
c-----------------------------------------------
c   local variables
c
      integer :: ierr, iopen, m, n, mn, mpol1,ia
      character*(*) :: extension

c-----------------------------------------------
      IF (ALLOCATED(bsubu)) DEALLOCATE(bsubu)
      IF (ALLOCATED(bsubv)) DEALLOCATE(bsubv)
      IF (ALLOCATED(cr)) DEALLOCATE(cr)
      IF (ALLOCATED(cz)) DEALLOCATE(cz)
      IF (ALLOCATED(cl)) DEALLOCATE(cl)
      allocate (bsubu(0:md,-nd:nd), bsubv(0:md,-nd:nd),
     1      cr(0:md,-nd:nd), cz(0:md,-nd:nd),cl(0:md,-nd:nd), stat=ierr)
      bsubu = 0;  bsubv = 0; cr = 0; cz = 0; cl = 0;
      if (ierr .ne. 0) stop 'Allocation error in bn_read_vmecf90'
c-----------------------------------------------
!
!     THIS MODULE SUBROUTINE LOADS UP ARRAYS, CONSTANTS READ IN FROM WOUT FILE
!
      ! SAL - PPPL modification to avoid bad allocations
      IF (.not. lwout_opened) THEN
         CALL read_wout_file(TRIM(extension),ierr,iopen)
      ELSE
         iopen = 0
         ierr  = 0
      END IF
!      ia = index(extension, '.nc')
!      IF (ia .gt. 0) THEN
!         call read_wout_file ('wout_' // TRIM(extension),ierr, iopen)
!         extension=extension(1:len_trim(extension)-3)
!      ELSE
!         call read_wout_file ('wout.' // TRIM(extension), ierr, iopen)
!      END IF
      if (iopen .ne. 0) stop 'error opening wout in bn_read_vmecf90'
      if (ierr .ne. 0) stop 'error reading wout in bn_read_vmecf90'

      ! IF LASYM is decteded then stop the code because we don't handle
      ! it yet.
      lasym_bn = .false.
      IF (lasym) THEN
         lasym_bn = .true.
      END IF
      
      mpol1  = mpol-1

      if(mf.lt.mpol1 .or. nf.lt.ntor) then
         print *, 'increase number of poloidal and/or toroidal modes:',
     1            ' mf,nf'
         print *, 'mpol1, ntor = ', mpol1, ntor, ' ; mf, nf = ', mf, nf
         stop
      endif
c---------------------------------------------------------------------
      IF (ALLOCATED(ixm)) DEALLOCATE(ixm)
      IF (ALLOCATED(ixn)) DEALLOCATE(ixn)
      IF (ALLOCATED(raxis_in)) DEALLOCATE(raxis_in)
      IF (ALLOCATED(zaxis_in)) DEALLOCATE(zaxis_in)
      IF (ALLOCATED(raxis_s)) DEALLOCATE(raxis_s)
      IF (ALLOCATED(zaxis_c)) DEALLOCATE(zaxis_c)
      IF (ALLOCATED(bsubus)) DEALLOCATE(bsubus)
      IF (ALLOCATED(bsubvs)) DEALLOCATE(bsubvs)
      IF (ALLOCATED(crs)) DEALLOCATE(crs)
      IF (ALLOCATED(czc)) DEALLOCATE(czc)
      IF (ALLOCATED(clc)) DEALLOCATE(clc)
      allocate (ixm(mnmax), ixn(mnmax), raxis_in(0:ntor),
     1          zaxis_in(0:ntor))
      allocate (raxis_s(0:ntor),zaxis_c(0:ntor))
      allocate (bsubus(0:md,-nd:nd), bsubvs(0:md,-nd:nd),
     1      crs(0:md,-nd:nd), czc(0:md,-nd:nd),clc(0:md,-nd:nd),
     1      stat=ierr)
      bsubus = 0;  bsubvs = 0; crs = 0; czc = 0; clc = 0;
      raxis_s = 0; zaxis_c = 0;

      raxis_in(0:ntor) = raxis(0:ntor,1)
      zaxis_in(0:ntor) = zaxis(0:ntor,1)
      mnmax_in = mnmax
      nfp_in = nfp

      do mn = 1, mnmax
         ixm(mn) = nint(xm(mn))
         ixn(mn) =-nint(xn(mn))/nfp   !!Flip sign: NESCOIL convention
         m = ixm(mn)
         n = ixn(mn)
         bsubu(m,n) = 1.5_dp*bsubumnc(mn,ns) - 0.5_dp*bsubumnc(mn,ns-1)
         bsubv(m,n) = 1.5_dp*bsubvmnc(mn,ns) - 0.5_dp*bsubvmnc(mn,ns-1)
         cl(m,n) = 1.5_dp*lmns(mn,ns) - 0.5_dp*lmns(mn,ns-1)
         cr(m,n) = rmnc(mn,ns)
         cz(m,n) = zmns(mn,ns)
      end do
      
      if (lasym_bn) then
         raxis_s(0:ntor)=raxis(0:ntor,2)
         zaxis_c(0:ntor)=zaxis(0:ntor,2)
         do mn = 1, mnmax
            m = ixm(mn)
            n = ixn(mn)
            bsubu(m,n) = 1.5_dp*bsubumns(mn,ns)-0.5_dp*bsubumns(mn,ns-1)
            bsubv(m,n) = 1.5_dp*bsubvmns(mn,ns)-0.5_dp*bsubvmns(mn,ns-1)
            clc(m,n) = 1.5_dp*lmnc(mn,ns) - 0.5_dp*lmnc(mn,ns-1)
            crs(m,n) = rmns(mn,ns)
            czc(m,n) = zmnc(mn,ns)
         end do
      end if

      iota_edge = 1.5_dp*iotas(ns) - 0.5_dp*iotas(ns-1)
      phip_edge = 1.5_dp*phip (ns) - 0.5_dp*phip (ns-1)

      if (coil_separation .le. 0._dp) then
         coil_separation = abs(cr(1,0))
         print *,' A default coil-plasma separation was chosen: ',
     1   coil_separation
         print *,' You may enter this value as the 2nd arg ',
     1           'on the command line'
      end if

      np  = nfp
      nvp = nv*np
      nuvp = nu*nv*np
      mb  = mpol1
      nb  = ntor

!
!     Deallocate memory in READ_WOUT module
!
!      call read_wout_deallocate

      end subroutine bn_read_vmecf90
