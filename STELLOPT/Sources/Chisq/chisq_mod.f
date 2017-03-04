      MODULE chisq_mod
      USE stel_kinds
      IMPLICIT NONE
      INTEGER, PARAMETER :: ntargets=76
      INTEGER, PARAMETER :: unit_dkes=29, unit_cobra=28,
     1                      unit_outdata=29, unit_diagno=30,
     2                      unit_isote=31
      INTEGER, PARAMETER :: ivar_aspect=1, ivar_maxj=2, ivar_beta=3,
     1  ivar_pedge=4,    ivar_pgrad=5,   ivar_iota=6,   ivar_iota_p=7,
     2  ivar_bmn=8,      ivar_bmin=9,    ivar_bmax=10,  ivar_jstar=11,
     3  ivar_jinvar=12,  ivar_ripple=13, ivar_dkes=14,  ivar_neo=15,
     4  ivar_pseudo=16,  ivar_orbit=17,  ivar_well=18,  ivar_mercier=19,
     5  ivar_balloon=20, ivar_kink=21,   ivar_jac=22,   ivar_dsubr=23,
     6  ivar_bootsj=24,  ivar_fluxp=25,  ivar_rbtor=26, ivar_curve=27,
     7  ivar_center=28,  ivar_ellipticity=29,           ivar_kappa=30,
     8  ivar_zmax=31,    ivar_vv=32,     ivar_vv_rms=33,
     9  ivar_coil_complexity=34,         ivar_coil_jmax=35,
     A  ivar_berr_ave=36,  ivar_coil_len=37,   ivar_coil_sep=38,
     B  ivar_coil_curv=39, ivar_coil_access=40,   ivar_coil_icurv=41,
     C  ivar_coil_ymin=42, ivar_coil_vfreg=43, ivar_extcur=44,
     D  ivar_oh=45, ivar_vac_island=46, ivar_coil_polcur=47,
     E  ivar_coil_pldist=48, ivar_coil_plcross=49, ivar_coil_reg=50,
     F  ivar_coil_curdens=51, ivar_iota_bounds=52, ivar_coil_rmax=53,
     G  ivar_coil_vfrmax=54, ivar_vacuum_field=55, ivar_coil_sep_mat=56,
     H  ivar_diagnostics=57, ivar_vacislandres=58, ivar_coil_dis_bds=59,
     I  ivar_vv_max=60, ivar_coil_aBerr=61,
     J  ivar_coil_mBerr=62, ivar_Jconf=63, ivar_p_prof=64,
     K  ivar_emis=65, ivar_pdamp=66, ivar_emisdamp=67, ivar_diagno=68,
     L  ivar_curtor=69, ivar_jedge=70, ivar_mse=71, ivar_isote= 72,
     M  ivar_iota_extrap = 73, ivar_ne = 74, ivar_te = 75, ivar_ti = 76

      INTEGER, DIMENSION(:), ALLOCATABLE :: index_array
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: wegt,
     1   chisq_match, chisq_target
      CHARACTER(len=2) :: comm1
      CHARACTER(len=18), DIMENSION(:), ALLOCATABLE :: chisq_descript

      CHARACTER(len=*), DIMENSION(ntargets), PARAMETER :: descript =
     1   ( /
     1   'Aspect Ratio      ','Max. Current      ','<Beta>            ',
     2   'Edge Pressure     ','Pressure Gradient ',
     3   'Iota (1/q)        ','d-Iota/ds         ',
     4   'Bmn Spectrum      ','Bmin Spectrum     ','Bmax Spectrum     ',
     5   'J*  Spectrum      ','J-Invariant       ','<Magnetic Ripple> ',
     6   'DKES L11 Transport','Ripple Diffusion  ','Pseudo Symmetry   ',
     7   'Orbit Confinement ','Magnetic Well     ','Mercier           ',
     8   'Ballooning Stab.  ','External Kink     ','Resonant Jacobian ',
     9   'Resistive D_R     ',
     A   'Bootstrap Current ','Poloidal Flux     ','R-Btor            ',
     B   'Surface Curvature ','Radial Centering  ','Ellipticity       ',
     C   '<Elongation>      ','Zmax Boundary     ',
     D   '3D Boundary       ','3D Shape-RMS      ',
     E   'Nescoil Complexity','Nescoil Sheet Jmax','Nescoil Berr-ave  ',
     F   'Coil Length       ','Coil Separation   ','Coil Curvature    ',
     G   'Coil-Access Reqts ','Coil Integ. Curv. ','Coil ymin         ',
     H   'Coil VF Regulariz.','Coil Current Reg. ','OH Constraint     ',
     I   'Vac Island Width  ','Coil I-Poloidal   ',
     J   'Coil-Plasma Dist. ',
     J   'Coil-Pl X CoilCoil','Coil MC Regulariz.','Coil Current Dens.',
     K   'Iota Bounds       ','Coil MC max radius','Coil VF max radius',
     L   'Vac Fld Bnorm Res.','C-C Separ Matrix  ','Diagnostics Match ',
     M   'Vac Island Residue','Coil Displ Bounds ',
     N   '3D Shape-Max dev  ','Coil Bnorm avg err','Coil Bnorm max err',
     O   'J-contour confine.','Pressure Profile  ','SXR Emissivity    ',
     P   'Press Prof Damping','Emiss Prof Damping','DIAGNO mag. diag. ',
     Q   'Total tor. current','Edge current dens.','MSE Diagnostic    ',
     R   'Iso-Te Constraint ','Iota Extrapolated ','Electron Density  ',
     S   'Electron Temp.    ','Ion Temperature   '/ )

      END MODULE chisq_mod
