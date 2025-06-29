
      subroutine bn_read_boozer(extension)
      use  meshes
      use  bnvariables
      use neswrite, only: coil_separation, mnmax_in => mnmax, ixm, ixn,
     1   raxis_in => raxis, zaxis_in => zaxis, nfp_in => nfp,
     2   iota_edge, phip_edge, raxis_s, zaxis_c, lasym_bn
      !use read_wout_mod                  !Error-free reading of wout file
      use read_boozer_mod
      implicit none
c-----------------------------------------------
c   local variables
c
      integer :: ierr, iopen, m, n, mn, mpol1,ia
      character*(*) :: extension
      DOUBLE PRECISION, ALLOCATABLE :: mfact(:,:)

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
      !IF (.not. lwout_opened) THEN
      !   CALL read_wout_file(TRIM(extension),ierr,iopen)
      !ELSE
      !   iopen = 0
      !   ierr  = 0
      !END IF
      iopen = 0
      ierr  = 0
      IF (.not.ALLOCATED(rmnc_b))
     1   CALL read_boozer_file(TRIM(extension),ierr,iopen)

      if (iopen .ne. 0) stop 'error opening wout in bn_read_boozer'
      if (ierr .ne. 0) stop 'error reading wout in bn_read_boozer'

      ! IF LASYM is decteded then stop the code because we don't handle
      ! it yet.
      lasym_bn = .false.
      IF (lasym_b) THEN
         lasym_bn = .true.
      END IF
      
      mpol1  = mboz_b-1

      if(mf.lt.mpol1 .or. nf.lt.nboz_b) then
         print *, 'increase number of poloidal and/or toroidal modes:',
     1            ' mf,nf'
         print *, 'mpol1, ntor = ', mpol1, nboz_b,' ; mf, nf = ',mf, nf
         stop
      endif

      if(md.lt.MAXVAL(ixm_b) .or. nd.lt.MAXVAL(ixn_b)/nfp_b) then
         print *, 'increase number of poloidal and/or toroidal modes:',
     1            ' md,nd'
         print *, 'xm_nyq, xn_nyq = ', MAXVAL(ixm_b),  
     1            MAXVAL(ixn_b)/nfp_b, ' ; md, nd = ', md, nd
         stop
      endif

      if(md.lt.mf .or. nd.lt.nf) then
         print *, 'increase number of poloidal and/or toroidal modes:',
     1            ' md,nd'
         print *, 'mf, nf= ', mf, nf, ' ; md, nd = ', md, nd
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
      allocate (ixm(mnboz_b), ixn(mnboz_b), 
     1          raxis_in(0:nboz_b), zaxis_in(0:nboz_b))
      allocate (raxis_s(0:nboz_b),zaxis_c(0:nboz_b))
      allocate (bsubus(0:md,-nd:nd), bsubvs(0:md,-nd:nd),
     1      crs(0:md,-nd:nd), czc(0:md,-nd:nd),clc(0:md,-nd:nd),
     1      stat=ierr)
      bsubus = 0;  bsubvs = 0; crs = 0; czc = 0; clc = 0;
      raxis_s = 0; zaxis_c = 0;

      ! Not used in code
      !raxis_in(0:ntor) = raxis(0:ntor,1)
      !zaxis_in(0:ntor) = zaxis(0:ntor,1)
      mnmax_in = mnboz_b
      nfp_in = nfp_b


      do mn = 1, mnmax_in
         ixm(mn) = ixm_b(mn)
         ixn(mn) = ixn_b(mn)/nfp_b
         m = ixm(mn)
         n = ixn(mn)
      end do

      ALLOCATE(mfact(mnboz_b,2))
      WHERE (MOD(ixm_b(:),2) .eq. 0)
         mfact(:,1)= 1.5
         mfact(:,2)=-0.5
      ELSEWHERE
         mfact(:,1)= 1.5*SQRT((ns_b-1.0)/(ns_b-1.5))
         mfact(:,2)=-0.5*SQRT((ns_b-1.0)/(ns_b-2.5))
      ENDWHERE
      do mn = 1, mnmax_in
         m = ixm(mn)
         n = ixn(mn)
         cr(m,n) = mfact(mn,1)*rmnc_b(mn,ns_b)
     1           + mfact(mn,2)*rmnc_b(mn,ns_b-1)
         cz(m,n) = mfact(mn,1)*zmns_b(mn,ns_b)
     1           + mfact(mn,2)*zmns_b(mn,ns_b-1)
         cl(m,n) = 0.0 ! only used in bn_write_nescoil_input
         bsubu(m,n) = 1.5*buco_b(ns_b) - 0.5*buco_b(ns_b-1)
         bsubv(m,n) = 1.5*bvco_b(ns_b) - 0.5*bvco_b(ns_b-1)
      end do
      if (lasym_bn) then
         !raxis_s(0:ntor)=raxis(0:ntor,2)
         !zaxis_c(0:ntor)=zaxis(0:ntor,2)
         do mn = 1, mnmax_in
            m = ixm(mn)
            n = ixn(mn)
            crs(m,n) = mfact(mn,1)*rmns_b(mn,ns_b) 
     1              + mfact(mn,2)*rmns_b(mn,ns_b-1)
            czc(m,n) = mfact(mn,1)*zmnc_b(mn,ns_b) 
     1              + mfact(mn,2)*zmnc_b(mn,ns_b-1)
            clc(m,n) = 0.0 ! only used in bn_write_nescoil_input
         end do
      end if
      DEALLOCATE(mfact)

      iota_edge = 1.5_dp*iota_b(ns_b) - 0.5_dp*iota_b(ns_b-1)
      phip_edge = 1.5_dp*phip_b(ns_b) - 0.5_dp*phip_b(ns_b-1)

      if (coil_separation .le. 0._dp) then
         coil_separation = abs(cr(1,0))
         print *,' A default coil-plasma separation was chosen: ',
     1   coil_separation
         print *,' You may enter this value as the 2nd arg ',
     1           'on the command line'
      end if

      np  = nfp_b
      nvp = nv*np
      nuvp = nu*nv*np
      mb  = mpol1
      nb  = nboz_b

!
!     Deallocate memory in READ_WOUT module
!
!      call read_wout_deallocate

      end subroutine bn_read_boozer
