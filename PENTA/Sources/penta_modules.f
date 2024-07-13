c---------------------------------------------------------------------------------
c   This module contains the kind specifications for PENTA and all subroutines 
c    and modules.
c   
      module penta_kind_mod
      implicit none
      integer, parameter :: rknd=SELECTED_REAL_KIND(15,300) 
      integer, parameter :: iknd=SELECTED_INT_KIND(8)       
      end module penta_kind_mod
c
c-----------------------------------------------------------------
c   Constants used in several routines
c-----------------------------------------------------------------
c
      module phys_const
      use penta_kind_mod
      implicit none
      real(rknd),parameter :: p_mass = 1.672621637e-27_rknd !proton mass
      real(rknd),parameter :: e_mass = 9.10938215e-31_rknd  !electron mass
      real(rknd),parameter :: elem_charge = 1.602176487e-19_rknd !elementary charge
      real(rknd),parameter :: eps0 = 8.854187817e-12_rknd 
      real(rknd),parameter :: pi = 3.14159265358979323846264338328_rknd 
      end module phys_const
c
c-----------------------------------------------------------------
c   File I/O unit numbers
c-----------------------------------------------------------------
c
      module io_unit_spec
      integer, parameter :: iu_nl=21, iu_vmec=22, iu_pprof=23, 
     1  iu_coeff=24, iu_Ufile=25, iu_flux_out=10, iu_pprof_out=11,
     2  iu_fvEr_out=12, iu_flows_out=13
      end module io_unit_spec
c
c------------------------------------------------------------------
c   Ambipolar roots passed from find_Er_roots to the main program
c------------------------------------------------------------------
c
      module Er_roots_pass
      use penta_kind_mod
      implicit none
      integer(iknd) :: num_roots
      real(rknd), allocatable :: er_roots(:)
      end module Er_roots_pass
c
c------------------------------------------------------------------
c   Parameters set in main program passed as variables
c------------------------------------------------------------------
c
      module parameter_pass
      use penta_kind_mod
      implicit none
      integer(iknd),parameter :: kcord=2_iknd, keord=2_iknd
      logical :: use_quanc8, log_interp
      real(rknd) :: Kmin, Kmax, epsabs, epsrel
      integer(iknd) :: num_K_steps
      end module parameter_pass
c
c------------------------------------------------------------------
c   Variables passed from energy_conv to intfun
c------------------------------------------------------------------
c
      module intfun_pass
      use penta_kind_mod
      implicit none
      integer(iknd) :: jval, num_c_pass, num_e_pass
      real(rknd) :: emax_pass, emin_pass, cmax_pass, cmin_pass
      real(rknd), dimension(:), allocatable :: xt_c_pass,xt_e_pass
      real(rknd), dimension(:,:), allocatable :: c_spl_pass
      logical :: log_coeff
      end module intfun_pass
c
c------------------------------------------------------------------
c   Variables passed from find_Er_roots to rad_flux
c------------------------------------------------------------------
c
      module Er_fit_pass
      use penta_kind_mod
      implicit none
      real(rknd), allocatable, dimension(:) ::  diff_qg, Er_fit
      integer(iknd),parameter :: num_pts=2_iknd
      end module Er_fit_pass
c
c------------------------------------------------------------------
c   Variables passed from define_thermal_mats to energy_conv
c      and intfun
c------------------------------------------------------------------
c
      module thermal_pass
      use penta_kind_mod
      implicit none
      real(rknd):: ma, qa, na, vta, abs_Er, loglam
      real(rknd), allocatable, dimension(:) :: vtb, nb, qb
      end module thermal_pass
c
c------------------------------------------------------------------
c   Variables read from profile_data_***
c------------------------------------------------------------------
c
      module vmec_var_pass
      use penta_kind_mod
      implicit none
      real(rknd) :: arad, Rmajor, r_surf, roa_surf, chip
     1 ,psip, btheta, bzeta, vp, bsq, iota, vol_p
      end module vmec_var_pass
c
c------------------------------------------------------------------
c   Variables read from plasma_profiles*.dat
c------------------------------------------------------------------
c
      module pprof_pass
      use penta_kind_mod
      implicit none
      real(rknd) :: ne, Te, dnedr, dTedr
      real(rknd), allocatable :: ni(:), Ti(:), dnidr(:), dTidr(:)
      real(rknd), dimension(:,:), allocatable :: ni_prof, Ti_prof
      end module pprof_pass
c
c------------------------------------------------------------------
c   Variables read from L,M,N * coefficient files
c------------------------------------------------------------------
c
      module coeff_var_pass
      use penta_kind_mod
      implicit none
      real(rknd), dimension(:), allocatable :: cmul_ls, efield_ls,
     >  cmul_ms, efield_ms, cmul_ns, efield_ns
      real(rknd), dimension(:,:), allocatable :: coef2d_ls,
     >  coef2d_ms, coef2d_ns
      integer(iknd) :: num_cl, num_el, num_cm, num_em, num_cn, num_en
      real(rknd), dimension(:), allocatable :: xt_cl, xt_el, xt_cm,
     1 xt_em, xt_cn, xt_en
      real(rknd), dimension(:,:), allocatable :: c_spll, c_splm, c_spln
      real(rknd) :: eminl, emaxl, eminm, emaxm, eminn, emaxn,
     1  cminl, cmaxl, cminm, cmaxm, cminn, cmaxn
      end module coeff_var_pass