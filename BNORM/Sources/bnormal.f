
!----------------------------------------------------------------------
      subroutine bnormal(n_u,n_v,mfou,nfou,mdim,ndim,bnfou,
     1                   bnfou_c,extension,lbooz)
!----------------------------------------------------------------------
c                                                          15.08.98
c    purpose:
c

c    compuce
c                   /
c               1   |   [ j', x-x'] . n
c     B.n   =  ---- |  ----------------- df ,  x,x' on the boundary
c              4 pi |     |x-x'|^(3/2)
c                   /
c
c
c       j   = B_subu . xv - B_subv . xu
c
c the result ist normalized so that the toroidal integral
c
c              /
c              |
c     I_p  =   |  B. ds  per field period   is I_p  =1.
c              |
c              /
c ----------------------------------------------------------------------
      use meshes
      use stel_kinds
      implicit none

      integer :: n_u, n_v, mfou, nfou, mdim, ndim
      real(rprec), dimension(0:mfou,-nfou:nfou) :: bnfou, bnfou_c
      character*(*) :: extension
      logical :: lbooz
         nu     = n_u
         nv     = n_v
         nuv    = nu*nv
         nuvh1  = nu*nv/2+nu
         mf     = mfou
         nf     = nfou
         md     = mdim
         nd     = ndim
      IF (lbooz) THEN
         CALL bn_read_boozer(extension)
      ELSE
         call bn_read_vmecf90(extension)
      END IF
      ! Allocate arrays
      call bn_alloc
      ! Compute Fourier Transform helpers
      call bn_precal
      ! Compute B_u B_v
      call bn_bfield_parallel
      ! Compute X,Y,Z
      ! Compute dX/dU,dY/dU,dZ/dU
      ! Compute dX/dV,dY/dV,dZ/dV
      ! Compute djx,djy,djz
      call bn_surface
      !Comput vector potential
      call bn_vecpot
      call bn_fouri(bnfou,bnfou_c)

      end subroutine bnormal
