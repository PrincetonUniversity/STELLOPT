subroutine quasisymmetry_single_solve

  use quasisymmetry_variables
  use vmec_input, only: ntor   !1/23/19.

  implicit none

  integer :: iteration, new_N_phi

  iota_tolerance_achieved = .false.
  elongation_tolerance_achieved = .false.
  iteration = 0
  N_phi = N_phi_original
  do 
     iteration = iteration + 1
     if (verbose) then
        print "(a)"," -------------------------------------------------------"
        print "(a,i4)"," Solving system using N_phi=",N_phi
     end if
     
     call quasisymmetry_init_phi()

     call quasisymmetry_init_axis()

     call quasisymmetry_solve()

     if (trim(resolution_option) == resolution_option_fixed) exit

     if (iteration > 1) then
        if (verbose) print "(a,es10.3)","                                      abs(iota - last_iota) =",abs(iota - last_iota)
        if (verbose) print "(a,es10.3)"," abs(max_elongation - last_max_elongation) / max_elongation =",abs(max_elongation - last_max_elongation) / max_elongation
        if (abs(iota - last_iota) <= iota_tolerance) then
           if (verbose) print *,"iota_tolerance (absolute) achieved."
           iota_tolerance_achieved = .true.
        end if
        !if (abs(max_elongation - last_max_elongation) <= elongation_tolerance) then  ! Absolute tolerance
        if (abs(max_elongation - last_max_elongation) / max_elongation <= elongation_tolerance) then  ! Relative tolerance
           if (verbose) print *,"elongation_tolerance (relative) achieved."
           elongation_tolerance_achieved = .true.
        end if
        if (max_elongation > max_precise_elongation) then
           if (verbose) print "(a)"," Ignoring elongation_tolerance and iota_tolerance since max_elongation > max_precise_elongation."
           elongation_tolerance_achieved = .true.
           iota_tolerance_achieved = .true.
        end if
        if (iota_tolerance_achieved .and. elongation_tolerance_achieved) exit
     end if

     last_iota = iota
     last_max_elongation = max_elongation

     new_N_phi = N_phi * 2 + 1
     write(0,*)'qs_single. new_N_phi,N_phi,ntor=',new_N_phi,N_phi,ntor  !1/23/19.(8k11g)
     if (new_N_phi > max_N_phi) then
        if (verbose) print *,"Stopping N_phi refinement since max_N_phi exceeded."
        exit
     end if
     N_phi = new_N_phi
  end do

end subroutine quasisymmetry_single_solve
