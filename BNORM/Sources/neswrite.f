
      module neswrite
      use stel_kinds
      logical :: lasym_bn
      integer :: ntheta, nzeta, nfp, mnmax, mnmax_ws
      integer :: nu, nv, mf, nf, md, nd
      integer, dimension(:), allocatable :: ixm, ixn, ixm_ws, ixn_ws
      real(rprec), allocatable, dimension(:) :: raxis, zaxis, rmnc, zmns
      real(rprec), allocatable, dimension(:) :: raxis_s, zaxis_c, 
     1                                          rmns, zmnc
      real(rprec), allocatable, dimension(:) :: rmnc_ws, zmns_ws
      real(rprec), allocatable, dimension(:) :: rmns_ws, zmnc_ws
      real(rprec) :: coil_separation
      real(rprec) :: iota_edge, phip_edge
      end module neswrite
