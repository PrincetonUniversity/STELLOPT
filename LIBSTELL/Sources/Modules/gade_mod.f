      MODULE gade_mod
      USE stel_kinds
      IMPLICIT NONE

!    shared
      INTEGER, PARAMETER :: indmax=500, nchrmax=5000, nparmax=500
      INTEGER, PARAMETER :: max_gen=1000
      INTEGER, PARAMETER :: gade_cleanup = -200
      INTEGER, PARAMETER :: pso_cleanup = -300

      INTEGER :: npopsiz, ngen, idum, ibound
      REAL(rprec)  :: pcross
      REAL(rprec), DIMENSION(nparmax) :: parmax, parmin
      LOGICAL :: save_space
!
!   GA Specific
      INTEGER :: nowrite,microga,unique_ind
      REAL(rprec)  :: pmutate,pcreep
      INTEGER :: iskip,iend,nchild,itourny,ielite,icreep,iunifrm,
     +           iniche
      INTEGER, DIMENSION(nparmax) :: nposibl, nichflg
!
!   DE specific
      INTEGER :: strategy, CR_strategy, out_iter
      REAL(rprec) :: f_cross

      NAMELIST /ga_de/ npopsiz,ngen,idum,ibound,pcross,
     +          strategy, CR_strategy, out_iter, f_cross,
     +          nowrite,microga,unique_ind,pmutate,pcreep,
     +          iskip,iend,nchild,itourny,ielite,icreep,iunifrm,
     +          iniche,save_space,
     +          parmin,parmax,nposibl,nichflg

      CONTAINS

      SUBROUTINE read_gade_namelist (iunit, istat)
      INTEGER :: iunit, istat

      READ (iunit, nml=ga_de, iostat=istat)

      END SUBROUTINE read_gade_namelist

      END MODULE gade_mod
