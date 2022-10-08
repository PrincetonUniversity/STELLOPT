!-----------------------------------------------------------------------
!     Module:        fieldlines_init_vars
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          11/06/2019
!     Description:   This subroutine just performs basic initalizations.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_vars
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------\
      USE fieldlines_runtime
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      pi = 4.0 * ATAN(1.0)
      pi2 = 8.0 * ATAN(1.0)
      mu0 = 16.0E-7 * ATAN(1.0)
      lverb    = .false.
      lvmec    = .false.
      leqdsk = .false.
      lpies    = .false.
      lspec    = .false.
      lcoil    = .false.
      lhint    = .false.
      lmgrid   = .false.
      lmu      = .false.
      lvessel  = .false.
      lvac     = .false.
      lpres    = .false.
      lrestart = .false.
      laxis_i  = .false.
      ladvanced = .false.
      lemc3 = .false.
      lerror_field = .false.
      lplasma_only = .false.
      lbfield_only = .false.
      lafield_only = .false.
      lreverse  = .false.
      lhitonly  = .false.
      lraw   = .false.
      lwall_trans = .false.
      ledge_start = .false.
      lnescoil    = .false.
      lmodb       = .false.
      lfield_start = .false.
      nruntype = runtype_old
      id_string     = ''
      coil_string   = ''
      mgrid_string  = ''
      vessel_string = ''
      restart_string = ''
      eqdsk_string = ''
      line_select = 96
      ldex_default = 0
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_vars
