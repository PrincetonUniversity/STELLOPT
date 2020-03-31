      SUBROUTINE GA_preset
      USE ga_mod
      IMPLICIT NONE

c
      kountmx=5
      irestrt=0
      itourny=0
      ielite=0
      iunifrm=0
      iniche=0
      iskip=0
      iend=0
      nchild=1
      unique_ind=0
      ibound=0
      nichflg(1:nparmax)=1
      parmax = 0
      parmin = 0
      nposibl = 0
      microga=0
      save_space = .false.

      END SUBROUTINE GA_preset
