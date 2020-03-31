
      module Vsolver1
      use Vmeshes
      real(rprec), dimension(:), allocatable :: a, work,
     1   result, dumv, bpx, bpy, bpz, svdw, daccu, bmag
      real(rprec), dimension(:,:), allocatable :: ab,
     1      bfx, bfy, bfz, bfn, svdv
      end module Vsolver1
