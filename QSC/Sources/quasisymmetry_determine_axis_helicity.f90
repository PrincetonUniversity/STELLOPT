subroutine quasisymmetry_determine_axis_helicity

  use quasisymmetry_variables

  implicit none

  integer :: j, counter
  integer, dimension(:), allocatable :: quadrant

  allocate(quadrant(N_phi+1))

  do j = 1, N_phi
     if (normal_cylindrical(j,1) >= 0) then
        if (normal_cylindrical(j,3) >= 0) then 
           quadrant(j) = 1
        else
           quadrant(j) = 4
        end if
     else
        if (normal_cylindrical(j,3) >= 0) then 
           quadrant(j) = 2
        else
           quadrant(j) = 3
        end if
     end if
  end do
  quadrant(N_phi+1) = quadrant(1)

  counter = 0
  do j = 1, N_phi
     if (quadrant(j) == 4 .and. quadrant(j+1) == 1) then
        counter = counter + 1
     else if (quadrant(j) == 1 .and. quadrant(j+1) == 4) then
        counter = counter - 1
     else
        counter = counter + quadrant(j+1) - quadrant(j)
     end if
  end do

  axis_helicity = counter / 4
  if (verbose) print *,"Axis helicity counter:",counter, "  axis_helicity:",axis_helicity
  if (modulo(counter,4) .ne. 0) stop "Axis helicity counter was not a multiple of 4"

  deallocate(quadrant)

end subroutine quasisymmetry_determine_axis_helicity
