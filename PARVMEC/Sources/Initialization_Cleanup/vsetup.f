      SUBROUTINE vsetup (iseq_count)
      USE vmec_main
      USE vacmod
      USE realspace
      USE vsvd
      USE mgrid_mod, ONLY: nbcoil_max, nlim_max, nextcur, mgrid_mode
      USE gmres_mod, ONLY: nfcn
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: iseq_count
C-----------------------------------------------
!
!     Reset default initial values
!
!     m=1 constraint (=t: apply correct, polar constraint; =f, apply approx. constraint)
      lconm1 = .true.
!      lconm1 = .false.
!	IF (lrfp) lconm1 = .false.    !SPH102109 converges better

!     2d preconditioner
      nfcn = 0

      lrecon = .true.
      loldout = .false.
      ledge_dump = .false.

      z00 = zero
      mgrid_mode = 'S'             !Assume scaled mode; read in from mgrid in free-bdy mode
      nextcur = 0
!
!     SOME DEFAULT VALUES WHICH MAY BE READ IN
!

!
!     Reconstruction stuff
!
      delbsq = one
      lpofr = .true.
      lpprof = .true.
      icurrout = 0
      total_chi_square_n = one
      total_chisq_n0 = total_chi_square_n
      nchistp = 0
      imse = -1
      itse = 0
      isnodes = 0
      ipnodes = 0
      iopt_raxis = 1
      imatch_phiedge = 1
      nflxs = 0
      nbfld = 0
      mseangle_offset = zero
      mseangle_offsetm = zero
      pres_offset = zero
      sigma_current = 1.e30_dp
      sigma_delphid = 1.e30_dp
      tensi = one
      tensp = one
      tensi2 = zero
      fpolyi = one
      presfac = one
      phidiam = 1.e30_dp

      mseprof = one
      indxflx = 0
      indxbfld = 0
      sigma_stark = 1.1*cbig
      sigma_thom = 1.1*cbig
      sigma_flux = 1.1*cbig
      sigma_b = 1.1*cbig

!
!     FREE-BOUNDARY STUFF, ONLY INITIALIZED FIRST TIME
!
      IF (iseq_count .eq. 0) THEN
        nbcoil_max = 0
        nlim_max = 0
      END IF

      END SUBROUTINE vsetup
