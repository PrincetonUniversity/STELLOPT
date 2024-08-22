cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   -PENTA with Impurities-
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  To run,type: PENTA_imp [command line arguments]
c
c  where the command line arguments are:
c
c  1  unique identifier text on DKES database files
c      i.e., for (m,l,n)star_lijs_***, this would be
c      the text that replaces *** for the specific surface (ex: hsx_s10)
c  2  Er,min
c  3  Er,max  - args 2 and 3 specify the interval in Er over which
c      the search is done for the ambipolar electric field root.
c      Er,min and Er,max are dimensionless normalized electric
c      field variables, i.e.: e<a>Er,min/max/kT_e
c    - Note if the logical variable input_is_Er below is true then 
c      the search interval is set by args 2 and 3 as Er (V/cm)
c  4  flux surface number
c  5  parameter that determines whether output data should be
c      written into new files (=0) or appended to
c      existing files (=1). Typically,when a script is run for a
c      sequence of flux surfaces, this is set to 0 for the first
c      surface and to 1 for subsequent surfaces
c  6  extension of the profile_data_*** file (i.e., the *** text)
c      that comes from extracting the appropriate profiles from
c      the VMEC run using the profile_extractor.f code
c  7  unique identifier on plasma profiles file, i.e., the ***
c      text in the filename plasma_profiles***.dat. This allows
c      multiple profiles to be run without having to rename files
c  8  <B*E||> in T*V/m
c
c
c  Required input files:
c
c   [l,m,n]star_lijs_***_s# - where *** is arg1 above and # arg 4.
c       These files contain the values of efield and cmul and the
c       corresponding normalized monoenergetic coefficients.  Note
c       that L* should be L*/e**2 where e is the elementary charge.
c
c   ion_params - A file containing the namelist "ion_params" which
c       defines the variables num_ion_species, Z_ion_init and 
c       miomp_init which are the number of ion species (integer),
c       the corresponding ion charge numbers (real, array) and 
c       the ion to proton mass ratio (real, array) for each 
c       non-electron species.
c
c   profile_data_*** - where *** is arg1 above.  Contains geometry
c       variables, most likely from a VMEC run.
c
c   plasma_profilesXXX.dat - A file containing the plasma profile
c       data, where XXX is arg8 above.  This file has the format:
c       row 1: number of radial points for which data is specified.
c       all other rows: r/a, ne, Te, Ti1, ni1, Ti2, ni2 ...
c       where the densities are in units of 10**12/cc and the 
c       temperatures in eV and i1, i2... are the ion species.
c
c	Utilde2_profile - Contains the quantity <U**2>, where U is the
c		PS flow function as defined in Sugama and Nishimura.  The
c		first row is the number of points, then r/a points and the
c		corresponding <U**2> value.
c   
c  Output files:
c
c   "fluxes_vs_roa", "fluxes_vs_Er", "flows_vs_roa", "plasma_profiles_check"
c
c   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   NOTES
c    - Be sure to check parameters below
c    - based on Sugama, Nishimura PoP 9 4637 (2002) and 
c         Sugama, Nishimura PoP 15, 042502 (2008)
c    - Some code used from orignal PENTA code, written by Don Spong and 
c        modified by Jeremy Lore
c
c  7/2009 JL
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   PROGRAMMING NOTES
c   -speed optimization
c   -recursive ambipolar solver
c     -fix numargs
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   2024 Notes
c    - Created PENTA_IMP_MAIN program and changed PENTA_IMP to a
c      callable subroutine.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
      PROGRAM PENTA_MAIN
      use penta_kind_mod
      IMPLICIT NONE
      integer(iknd) :: numargs, js, i_append
      real(rknd) :: Er_min, Er_max, B_Eprl
      character(100) :: arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8
      character(100) :: coeff_ext, run_ident
      character(10)  :: pprof_char
      numargs = 8
      call getarg(1, arg1)
      call getarg(2, arg2)
      call getarg(3, arg3)
      call getarg(4, arg4)
      call getarg(5, arg5)
      call getarg(6, arg6)
      call getarg(7, arg7)
      if (numargs .eq. 8) call getarg(8, arg8)

      !store command line args
      coeff_ext = trim(adjustl(arg1))
      read(arg2,*) Er_min
      read(arg3,*) Er_max
      read(arg4,*) js
      read(arg5,*) i_append
      run_ident = trim(adjustl(arg6))
      pprof_char = trim(adjustl(arg7))
      if (numargs .eq. 8) read(arg8,*) B_Eprl

      ! Call subroutine
      CALL PENTA_IMP(coeff_ext, Er_min, Er_max, js, i_append, 
     1               run_ident, pprof_char, B_Eprl)
      END PROGRAM PENTA_MAIN

