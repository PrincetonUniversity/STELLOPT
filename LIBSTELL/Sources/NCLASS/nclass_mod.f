      SUBROUTINE NCLASS(k_config,k_order,k_potato,m_i,m_z,c_den,c_potb,
     &                  c_potl,p_b2,p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,
     &                  p_grphi,p_gr2phi,p_ngrth,amu_i,grt_i,temp_i,
     &                  den_iz,fex_iz,grp_iz,m_s,jm_s,jz_s,p_bsjb,
     &                  p_etap,p_exjb,calm_i,caln_ii,capm_ii,capn_ii,
     &                  bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,sqz_s,upar_s,
     &                  utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,chip_ss,
     &                  chit_ss,dp_ss,dt_ss,iflag)
!***********************************************************************
!NCLASS calculates the neoclassical transport properties of a multiple
!  species axisymmetric plasma using k_order parallel and radial force
!  balance equations for each species
!References:
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079
!  W.A.Houlberg 6/99
!Input:
!  k_config -configuration
!       =1 stellarator
!       =else tokamak       
!  k_order-order of v moments to be solved (-)
!       =2 u and q
!       =3 u, q, and u2
!       =else error
!  k_potato-option to include potato orbits (-)
!          =0 off
!          =else on
!  m_i-number of isotopes (1<m_i<mx_mi+1)
!  m_z-highest charge state of all species (0<m_z<mx_mz+1)
!  c_den-density cutoff below which species is ignored (/m**3)
!  c_potb-kappa(0)*Bt(0)/[2*q(0)**2] (T)
!  c_potl-q(0)*R(0) (m)
!  p_b2-<B**2> (T**2)
!  p_bm2-<1/B**2> (/T**2)
!  p_eb-<E.B> (V*T/m)
!  p_fhat-mu_0*F/(dPsi/dr) (rho/m)
!  p_fm(3)-poloidal moments of geometric factor for PS viscosity (-)
!  p_ft-trapped fraction (-)
!  p_grbm2-<grad(rho)**2/B**2> (rho**2/m**2/T**2)
!  p_grphi-radial electric field Phi' (V/rho)
!  p_gr2phi-radial electric field gradient Psi'(Phi'/Psi')' (V/rho**2)
!  p_ngrth-<n.grad(Theta)> (1/m)
!  amu_i(i)-atomic mass number of i (-)
!  grt_i(i)-temperature gradient of i (keV/rho)
!  temp_i(i)-temperature of i (keV)
!  den_iz(i,z)-density of i,z (/m**3)
!  fex_iz(3,i,z)-moments of external parallel force on i,z (T*j/m**3)
!  grp_iz(i,z)-pressure gradient of i,z (keV/m**3/rho)
!Output:
!  m_s-number of species (1<ms<mx_ms+1)
!  jm_s(s)-isotope number of s (-)
!  jz_s(s)-charge state of s (-)
!  p_bsjb-<J_bs.B> (A*T/m**2)
!  p_etap-parallel electrical resistivity (Ohm*m)
!  p_exjb-<J_ex.B> current response to fex_iz (A*T/m**2)
!  calm_i(3,3,i)-tp eff friction matrix for i (kg/m**3/s)
!  caln_ii(3,3,i1,i2)-fp eff friction matrix for i1 on i2 (kg/m**3/s)
!  capm_ii(3,3,i1,i2)-test part (tp) friction matrix for i1 on i2 (-)
!  capn_ii(3,3,i1,i2)-field part (fp) friction matrix for i1 on i2 (-)
!  bsjbp_s(s)-<J_bs.B> driven by unit p'/p of s (A*T*rho/m**3)
!  bsjbt_s(s)-<J_bs.B> driven by unit T'/T of s (A*T*rho/m**3)
!  dn_s(s)-diffusion coefficient (diag comp) of s (rho**2/s)
!  gfl_s(m,s)-radial particle flux comps of s (rho/m**3/s)
!             m=1, banana-plateau, p' and T'
!             m=2, Pfirsch-Schluter
!             m=3, classical
!             m=4, banana-plateau, <E.B>
!             m=5, banana-plateau, external parallel force fex_iz
!  qfl_s(m,s)-radial heat conduction flux comps of s (W*rho/m**3)
!             m=1, banana-plateau, p' and T'
!             m=2, Pfirsch-Schluter
!             m=3, classical
!             m=4, banana-plateau, <E.B>
!             m=5, banana-plateau, external parallel force fex_iz
!  sqz_s(s)-orbit squeezing factor for s (-)
!  upar_s(3,m,s)-parallel flow of s from force m (T*m/s)
!                m=1, p', T', Phi'
!                m=2, <E.B>
!                m=3, fex_iz
!  utheta_s(3,m,s)-poloidal flow of s from force m (m/s/T)
!                  m=1, p', T'
!                  m=2, <E.B>
!                  m=3, fex_iz
!  vn_s(s)-convection velocity (off diag comps-p', T') of s (rho/s)
!  veb_s(s)-<E.B> particle convection velocity of s (rho/s)
!  qeb_s(s)-<E.B> heat convection velocity of s (rho/s)
!  xi_s(s)-charge weighted density factor of s (-)
!  ymu_s(s)-normalized viscosity for s (kg/m**3/s)
!  chip_ss(s1,s2)-heat cond coefficient of s2 on p'/p of s1 (rho**2/s)
!  chit_ss(s1,s2)-heat cond coefficient of s2 on T'/T of s1 (rho**2/s)
!  dp_ss(s1,s2)-diffusion coefficient of s2 on p'/p of s1 (rho**2/s)
!  dt_ss(s1,s2)-diffusion coefficient of s2 on T'/T of s1 (rho**2/s)
!  iflag-warning and error flag
!       =-4 warning: no viscosity
!       =-3 warning: no banana viscosity
!       =-2 warning: no Pfirsch-Schluter viscosity
!       =-1 warning: no potato orbit viscosity
!       =0 no warnings or errors
!       =1 error: order of v moments to be solved must be 2 or 3
!       =2 error: number of species must be 1<m_i<mx_mi+1
!       =3 error: number of species must be 0<m_z<mx_mz+1
!       =4 error: number of species must be 1<m_s<mx_ms+1
!       =5 error: inversion of flow matrix failed
!       =6 error: trapped fraction must be 0.0.le.p_ft.le.1.0
!***********************************************************************
      IMPLICIT NONE
      INCLUDE 'pamx_mi.inc'
      INCLUDE 'pamx_ms.inc'
      INCLUDE 'pamx_mz.inc'
!Declaration of input variables
      INTEGER        k_config
      INTEGER        k_order,                 k_potato
      INTEGER        m_i,                     m_z
      REAL           c_den,                   c_potb,
     &               c_potl
      REAL           p_b2,                    p_bm2,
     &               p_eb,                    p_fhat,
     &               p_fm(3),                 p_ft,
     &               p_grbm2,                 p_grphi,
     &               p_gr2phi,                p_ngrth
      REAL           amu_i(mx_mi),            grt_i(mx_mi),
     &               temp_i(mx_mi)
      REAL           den_iz(mx_mi,mx_mz),     fex_iz(3,mx_mi,mx_mz),
     &               grp_iz(mx_mi,mx_mz)
!Declaration of output variables
      INTEGER        iflag,                   m_s
      INTEGER        jm_s(mx_ms),             jz_s(mx_ms)
      REAL           p_bsjb,                  p_etap,
     &               p_exjb
      REAL           calm_i(3,3,mx_mi)
      REAL           caln_ii(3,3,mx_mi,mx_mi),capm_ii(3,3,mx_mi,mx_mi),
     &               capn_ii(3,3,mx_mi,mx_mi)
      REAL           bsjbp_s(mx_ms),          bsjbt_s(mx_ms),
     &               dn_s(mx_ms),             gfl_s(5,mx_ms),
     &               qfl_s(5,mx_ms),          sqz_s(mx_ms),
     &               upar_s(3,3,mx_ms),       utheta_s(3,3,mx_ms),
     &               vn_s(mx_ms),             veb_s(mx_ms),
     &               qeb_s(mx_ms),            xi_s(mx_ms),
     &               ymu_s(3,3,mx_ms)
      REAL           chip_ss(mx_ms,mx_ms),    chit_ss(mx_ms,mx_ms),
     &               dp_ss(mx_ms,mx_ms),      dt_ss(mx_ms,mx_ms)
!Declaration of local variables
      INTEGER        k_banana,                k_pfirsch
      INTEGER        i,                       im,
     &               iz,                      iza,
     &               jflag,                   jm,
     &               k,                       l
      REAL           dent
      REAL           z_coulomb,               z_j7kv,
     &               z_protonmass
      REAL           denz2(mx_mi),            vt_i(mx_mi)
      REAL           pgrp_iz(mx_mi,mx_mz)
      REAL           amnt_ii(mx_mi,mx_mi)
      REAL           tau_ss(mx_ms,mx_ms)
!Initialization
!  Error flag
      iflag=0
!  Consistency checks
!     Order must be two or three moments
      IF(k_order.lt.2.or.k_order.gt.3) THEN
        iflag=1
        GOTO 1000
      ENDIF
!     At least two but not greater than mx_mi species
      IF(m_i.lt.2.or.m_i.gt.mx_mi) THEN
        iflag=2
        GOTO 1000
      ENDIF
!     Highest charge state at least 1 but not greater than mx_mz
      IF(m_z.lt.1.or.m_i.gt.mx_mz) THEN
        iflag=3
        GOTO 1000
      ENDIF
!     Trapped fraction between 0 and 1 inclusive
      IF((p_ft.lt.0.0).or.(p_ft.gt.1.0)) THEN
        iflag=6
        GOTO 1000
      ENDIF
      IF(k_config .eq. 1) THEN
!     Workaround for the stellarator case
        k_potato = 0
        k_pfirsch = 0
        k_banana = 1
      ELSE  
!       Potato orbit contribution to viscosity
        
        IF((ABS(c_potb).gt.0.0).and.(ABS(c_potl).gt.0.0)
     &     .and.(k_potato.ne.0)) THEN
          k_potato=1
        ELSE
          k_potato=0
          iflag=-1
        ENDIF
!       Pfirsch-Schluter contribution to viscosity
        IF(ABS(p_fm(1)+p_fm(2)+p_fm(3)).gt.0.0) THEN
          k_pfirsch=1
        ELSE
          k_pfirsch=0
          iflag=-2
        ENDIF
!       Banana contribution to viscosity
        IF(ABS(p_ft).gt.0.0) THEN
          k_banana=1
        ELSE
          k_banana=0
          iflag=-3
        ENDIF
!       No viscsoity
        IF((k_banana.eq.0).and.(k_pfirsch.eq.0)) THEN
          k_potato=0
          iflag=-4
        ENDIF
      END IF  
!  Physical and conversion constants
      z_coulomb=1.6022e-19
      z_j7kv=1.6022e-16
      z_protonmass=1.6726e-27
!Find significant charge states and mapping
      m_s=0
      DO im=1,m_i
        DO iza=1,m_z
          IF(den_iz(im,iza).gt.c_den) THEN
            m_s=m_s+1
!           Set isotope number and charge state for this species
            jm_s(m_s)=im
            IF(amu_i(im).lt.0.5) THEN
              jz_s(m_s)=-iza
            ELSE
              jz_s(m_s)=iza
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      IF(m_s.lt.2.or.m_s.gt.mx_ms) THEN
        iflag=4
        GOTO 1000
      ENDIF
!Get friction coefficients
      CALL NCLASS_MN(k_order,m_i,amu_i,temp_i,capm_ii,capn_ii)
!Calculate thermal velocity
      DO im=1,m_i
        vt_i(im)=SQRT(2.0*z_j7kv*temp_i(im)/amu_i(im)/z_protonmass)
      ENDDO
!Get collision times
      CALL NCLASS_TAU(m_i,m_s,jm_s,jz_s,amu_i,temp_i,vt_i,den_iz,
     &                amnt_ii,tau_ss)
!Calculate reduced friction coefficients
      CALL RARRAY_ZERO(9*m_i,calm_i)
      DO im=1,m_i
        DO jm=1,m_i
          DO k=1,k_order
            DO l=1,k_order
!             Sum over isotopes b for test particle component
              calm_i(k,l,im)=calm_i(k,l,im)
     &                      +amnt_ii(im,jm)*capm_ii(k,l,im,jm)
!             Field particle component.
              caln_ii(k,l,im,jm)=amnt_ii(im,jm)*capn_ii(k,l,im,jm)            
            ENDDO   
          ENDDO    
        ENDDO   
      ENDDO   
!Calculate species charge state density factor, total nT, squeezing
      dent=0.0
      DO im=1,m_i
        denz2(im)=0.0
        DO iza=1,m_z
          IF(den_iz(im,iza).gt.c_den) THEN
            denz2(im)=denz2(im)+den_iz(im,iza)*iza**2
            dent=dent+den_iz(im,iza)*temp_i(im)
          ENDIF
        ENDDO   
      ENDDO   
      DO i=1,m_s                                                            
        im=jm_s(i)
        iz=jz_s(i)
        iza=IABS(iz)
        xi_s(i)=den_iz(im,iza)*iz**2/denz2(im)
        sqz_s(i)=1.0+p_fhat**2/p_b2*amu_i(im)*z_protonmass
     &               *ABS(p_gr2phi/(z_coulomb*iz))
      ENDDO
!Get normalized viscosities
      CALL NCLASS_MU(k_config,k_order,k_banana,k_pfirsch,k_potato,m_s,
     &            jm_s,jz_s,c_potb,c_potl,p_fm,p_ft,p_ngrth,amu_i,
     &            temp_i,vt_i,den_iz,sqz_s,ymu_s,tau_ss)
!Add potential gradient to pressure gradient
      DO i=1,m_s
        im=jm_s(i)
        iz=jz_s(i)
        iza=IABS(iz)
        pgrp_iz(im,iza)=grp_iz(im,iza)
     &                  +p_grphi*den_iz(im,iza)*iz*z_coulomb/z_j7kv
      ENDDO
!Get normalized parallel flows within a surface
      jflag=0
      CALL NCLASS_FLOW(k_order,m_i,m_s,jm_s,jz_s,p_b2,p_bm2,p_eb,p_fhat,
     &                 p_grbm2,grt_i,temp_i,calm_i,caln_ii,den_iz,
     &                 fex_iz,pgrp_iz,xi_s,ymu_s,p_bsjb,p_etap,p_exjb,
     &                 bsjbp_s,bsjbt_s,gfl_s,qfl_s,upar_s,chip_ss,
     &                 chit_ss,dp_ss,dt_ss,jflag)
      IF(jflag.gt.0) THEN
        iflag=5
        GOTO 1000
      ENDIF
!Calculate poloidal velocity from parallel velocity
      DO i=1,m_s
        im=jm_s(i)
        iz=jz_s(i)
        iza=IABS(iz)
        DO k=1,k_order
          DO l=1,3
            utheta_s(k,l,i)=upar_s(k,l,i)/p_b2
            IF(l.eq.1) THEN
              IF(k.eq.1) THEN
                utheta_s(k,l,i)=utheta_s(k,l,i)
     &                          +p_fhat*pgrp_iz(im,iza)*z_j7kv
     &                          /(z_coulomb*iz*den_iz(im,iza))/p_b2
              ELSEIF(k.eq.2) THEN
                utheta_s(k,l,i)=utheta_s(k,l,i)
     &                          +p_fhat*z_j7kv*grt_i(im)
     &                          /(iz*z_coulomb*p_b2)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!Convert to diffusivities and conductivities
!  Full coefficient matrices 
      DO i=1,m_s
        im=jm_s(i)
        iza=IABS(jz_s(i))
!  Diagonal diffusivity
        dn_s(i)=dp_ss(i,i)
!  Off-diagonal expressed as radial velocity
        vn_s(i)=(gfl_s(1,i)+gfl_s(2,i)+gfl_s(3,i)
     &          +dn_s(i)*(grp_iz(im,iza)-den_iz(im,iza)*grt_i(im))
     &          /temp_i(im))/den_iz(im,iza)
!  <E.B> particle and heat convection velocities
        veb_s(i)=gfl_s(4,i)/den_iz(im,iza)
        qeb_s(i)=qfl_s(4,i)/den_iz(im,iza)/temp_i(im)/z_j7kv
      ENDDO
 1000 RETURN
      END
      SUBROUTINE NCLASS_FLOW(k_order,m_i,m_s,jm_s,jz_s,p_b2,p_bm2,p_eb,
     &                       p_fhat,p_grbm2,grt_i,temp_i,calm_i,caln_ii,
     &                       den_iz,fex_iz,grp_iz,xi_s,ymu_s,p_bsjb,
     &                       p_etap,p_exjb,bsjbp_s,bsjbt_s,gfl_s,qfl_s,
     &                       upar_s,chip_ss,chit_ss,dp_ss,dt_ss,iflag)
!***********************************************************************
!NCLASS_FLOW calculates the k_order neoclassical flows u and q/p (and
!  u2) plus other transport properties
!References:                                                     
!  Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  W.A.Houlberg 6/99
!Input:
!  k_order-order of v moments to be solved (-)
!       =2 u and q
!       =3 u, q, and u2
!       =else error
!  m_i-number of isotopes (1<m_i<mx_mi+1)
!  m_s-number of species (1<ms<mx_ms+1)
!  jm_s(s)-isotope number of s (-)
!  jz_s(s)-charge state of s (-)
!  p_b2-<B**2> (T**2)
!  p_bm2-<1/B**2> (/T**2)
!  p_eb-<E.B> (V*T/m)
!  p_fhat-mu_0*F/(dPsi/dr) (rho/m)
!  p_grbm2-<grad(rho)**2/B**2> (rho**2/m**2/T**2)
!  grt_i(i)-temperature gradient of i (keV/rho)
!  temp_i(i)-temperature of i (keV)
!  calm_i(3,3,i)-tp eff friction matrix for i (kg/m**3/s)
!  caln_ii(3,3,i1,i2)-fp eff friction matrix for i1 on i2 (kg/m**3/s)
!  den_iz(i,z)-density of i,z (/m**3)
!  fex_iz(3,i,z)-moments of external parallel force on i,z (T*j/m**3)
!  grp_iz(i,z)-pressure gradient of i,z (keV/m**3/rho)
!  xi_s(s)-charge weighted density factor of s (-)
!  ymu_s(s)-normalized viscosity for s (kg/m**3/s)
!Output:
!  p_bsjb-<J_bs.B> (A*T/m**2)
!  p_etap-parallel electrical resistivity (Ohm*m)
!  p_exjb-<J_ex.B> current response to fex_iz (A*T/m**2)
!  bsjbp_s(s)-<J_bs.B> driven by unit p'/p of s (A*T*rho/m**3)
!  bsjbt_s(s)-<J_bs.B> driven by unit T'/T of s (A*T*rho/m**3)
!  gfl_s(m,s)-radial particle flux comps of s (rho/m**3/s)
!             m=1, banana-plateau, p' and T'
!             m=2, Pfirsch-Schluter
!             m=3, classical
!             m=4, banana-plateau, <E.B>
!             m=5, banana-plateau, external parallel force fex_iz
!  qfl_s(m,s)-radial heat conduction flux comps of s (W*rho/m**3)
!             m=1, banana-plateau, p' and T'
!             m=2, Pfirsch-Schluter
!             m=3, classical
!             m=4, banana-plateau, <E.B>
!             m=5, banana-plateau, external parallel force fex_iz
!  upar_s(3,m,s)-parallel flow of s from force m (T*m/s)
!                m=1, p' and T'
!                m=2, <E.B>
!                m=3, fex_iz
!  chip_ss(s1,s2)-heat cond coefficient of s2 on p'/p of s1 (rho**2/s)
!  chit_ss(s1,s2)-heat cond coefficient of s2 on T'/T of s1 (rho**2/s)
!  dp_ss(s1,s2)-diffusion coefficient of s2 on p'/p of s1 (rho**2/s)
!  dt_ss(s1,s2)-diffusion coefficient of s2 on T'/T of s1 (rho**2/s)
!  iflag-error flag
!       =0 no errors
!       =1 inversion of flow matrix failed
!***********************************************************************
      IMPLICIT NONE
      INCLUDE 'pamx_mi.inc'
      INCLUDE 'pamx_ms.inc'
      INCLUDE 'pamx_mz.inc'
!Declaration of input variables
      INTEGER        k_order,                 m_i,
     &               m_s
      INTEGER        jm_s(mx_ms),             jz_s(mx_ms)       
      REAL           p_b2,                    p_bm2,
     &               p_eb,                    p_fhat,                    
     &               p_grbm2
      REAL           temp_i(mx_mi),           grt_i(mx_mi),
     &               calm_i(3,3,mx_mi)
      REAL           caln_ii(3,3,mx_mi,mx_mi)
      REAL           den_iz(mx_mi,mx_mz),     grp_iz(mx_mi,mx_mz),
     &               fex_iz(3,mx_mi,mx_mz)
      REAL           xi_s(mx_ms),             ymu_s(3,3,mx_ms)
!Declaration of output variables
      INTEGER        iflag
      REAL           p_bsjb,                  p_etap,
     &               p_exjb
      REAL           bsjbp_s(mx_ms),          bsjbt_s(mx_ms),
     &               gfl_s(5,mx_ms),          qfl_s(5,mx_ms),
     &               upar_s(3,3,mx_ms)
      REAL           chip_ss(mx_ms,mx_ms),    chit_ss(mx_ms,mx_ms),
     &               dp_ss(mx_ms,mx_ms),      dt_ss(mx_ms,mx_ms)
!Declaration of local variables
      INTEGER        i,                       im,
     &               iz,                      iza,
     &               j,                       jm,
     &               k,                       l,                     
     &               l1,                      m,
     &               m1
      INTEGER        indx(3),                 indxb(3*mx_mi)
      REAL           cbp,                     cbpa,
     &               cbpaq,                   cc,
     &               ccl,                     ccla,
     &               cclaq,                   cclb,
     &               cclbq,                   cps,
     &               cpsa,                    cpsaq,
     &               cpsb,                    cpsbq,
     &               d,                       denzc,
     &               p_ohjb,                  z_coulomb,
     &               z_j7kv
      REAL           aa(3,3),                 xl(3)
      REAL           crhat(3,6,mx_mi),        rhat(3,6,mx_ms)
      REAL           srcth(3,mx_ms),          srcthp(mx_ms),            
     &               srctht(mx_ms)
      REAL           ab(3*mx_mi,3*mx_mi),     xab(3*mx_mi,3)
      REAL           crhatp(3,mx_ms,mx_mi),   crhatt(3,mx_ms,mx_mi)
      REAL           rhatp(3,mx_ms,mx_ms),    rhatt(3,mx_ms,mx_ms),
     &               uaip(3,mx_ms,mx_ms),     uait(3,mx_ms,mx_ms)
      REAL           xabp(3*mx_mi,mx_ms),     xabt(3*mx_mi,mx_ms)
!Initialization
!  Error flag
      iflag=0
!  Physical and conversion constants
      z_coulomb=1.6022e-19
      z_j7kv=1.6022e-16
!  Zero out arrays
      CALL RARRAY_ZERO(3*6*m_i,crhat)
      CALL RARRAY_ZERO(3*mx_ms*m_i,crhatp)
      CALL RARRAY_ZERO(3*mx_ms*m_i,crhatt)
      CALL RARRAY_ZERO(3*mx_mi*3*m_i,ab)
      CALL RARRAY_ZERO(5*m_s,gfl_s)
      CALL RARRAY_ZERO(5*m_s,qfl_s)
      p_etap=0.0
      p_bsjb=0.0
      p_exjb=0.0
      p_ohjb=0.0
      cc=(p_fhat/z_coulomb)*z_j7kv
      cbp=p_fhat/p_b2/z_coulomb
      cps=(p_fhat/z_coulomb)*(1.0/p_b2-p_bm2)
      ccl=(p_grbm2/z_coulomb)/p_fhat
!Calculate responses for each species
      DO i=1,m_s
        im=jm_s(i)
        iz=jz_s(i)
        iza=IABS(iz)
!  Set up response matrix for each charge state 
        DO k=1,k_order
          DO l=1,k_order
            aa(k,l)=xi_s(i)*calm_i(k,l,im)-ymu_s(k,l,i)
          ENDDO
        ENDDO
!  Get lu decomposition of response matrix
        CALL U_LU_DECOMP(aa,k_order,3,indx,d,iflag)
        IF(iflag.ne.0) GOTO 1000
!  Get sources and evaluate responses from back substitution 
!       Lambda terms involving isotopic flows 
        DO l=1,k_order
          DO k=1,k_order
            IF(k.eq.l) THEN
              rhat(k,l,i)=xi_s(i)
            ELSE
              rhat(k,l,i)=0.0
            ENDIF
          ENDDO
          CALL U_LU_BACKSUB(aa,k_order,3,indx,rhat(1,l,i))
          DO k=1,k_order
            crhat(k,l,im)=crhat(k,l,im)+xi_s(i)*rhat(k,l,i)
          ENDDO
        ENDDO
!       Poloidal source (p' and T') terms 
        srcth(1,i)=(cc/iz)*grp_iz(im,iza)/den_iz(im,iza)
        srcth(2,i)=(cc/iz)*grt_i(im)
        srcth(3,i)=0.0
        DO k=1,k_order
          rhat(k,4,i)=0.0
          DO l=1,k_order
            rhat(k,4,i)=rhat(k,4,i)+srcth(l,i)*ymu_s(k,l,i)
          ENDDO
        ENDDO
        CALL U_LU_BACKSUB(aa,k_order,3,indx,rhat(1,4,i))
        DO k=1,k_order
          crhat(k,4,im)=crhat(k,4,im)+xi_s(i)*rhat(k,4,i)
        ENDDO
!       Unit p'/p and T'/T terms for decomposition of fluxes 
        srcthp(i)=-(cc/iz)*temp_i(im)
        srctht(i)=-(cc/iz)*temp_i(im)
        DO k=1,k_order
          rhatp(k,i,i)=srcthp(i)*ymu_s(k,1,i)
          rhatt(k,i,i)=srctht(i)*ymu_s(k,2,i)
        ENDDO
        CALL U_LU_BACKSUB(aa,k_order,3,indx,rhatp(1,i,i))
        CALL U_LU_BACKSUB(aa,k_order,3,indx,rhatt(1,i,i))
        DO k=1,k_order
          crhatp(k,i,im)=crhatp(k,i,im)+xi_s(i)*rhatp(k,i,i)
          crhatt(k,i,im)=crhatt(k,i,im)+xi_s(i)*rhatt(k,i,i)
        ENDDO
!       Parallel electric field terms for resistivity 
        rhat(1,5,i)=-iz*z_coulomb*den_iz(im,iza)
        rhat(2,5,i)=0.0
        rhat(3,5,i)=0.0
        CALL U_LU_BACKSUB(aa,k_order,3,indx,rhat(1,5,i))
        DO k=1,k_order
          crhat(k,5,im)=crhat(k,5,im)+xi_s(i)*p_eb*rhat(k,5,i)
        ENDDO
!       External force terms 
        rhat(1,6,i)=-fex_iz(1,im,iza)
        rhat(2,6,i)=-fex_iz(2,im,iza)
        IF(k_order.eq.3) THEN
          rhat(3,6,i)=-fex_iz(3,im,iza)
        ELSE
          rhat(3,6,i)=0.0
        ENDIF
        CALL U_LU_BACKSUB(aa,k_order,3,indx,rhat(1,6,i))
        DO k=1,k_order
          crhat(k,6,im)=crhat(k,6,im)+xi_s(i)*rhat(k,6,i)
        ENDDO
      ENDDO   
!Load coefficient matrix and source terms for isotopic flows
      DO im=1,m_i
        DO m=1,k_order
          m1=im+(m-1)*m_i
!  Diagonal coefficients       
          ab(m1,m1)=1.0
!  Source terms 
!         p' and T'
          xab(m1,1)=crhat(m,4,im)
!         Unit p'/p and T'/T
          DO j=1,m_s
            xabp(m1,j)=crhatp(m,j,im)
            xabt(m1,j)=crhatt(m,j,im)
          ENDDO
!         <E.B>
          xab(m1,2)=crhat(m,5,im)
!         External source
          xab(m1,3)=crhat(m,6,im)
!  Field particle friction       
          DO jm=1,m_i
            DO l=1,k_order
              l1=jm+(l-1)*m_i
              DO k=1,k_order
                ab(m1,l1)=ab(m1,l1)+caln_ii(k,l,im,jm)*crhat(m,k,im)
              ENDDO  
            ENDDO
          ENDDO
        ENDDO        
      ENDDO   
!Get lu decomposition of coefficient matrix 
      CALL U_LU_DECOMP(ab,k_order*m_i,3*mx_mi,indxb,d,iflag)
      IF(iflag.ne.0) GOTO 1000
!Evaluate isotopic flows from back substitution for each source 
!  xab(1,k) to xab(m_i,k) are the isotopic velocities 
!  xab(m_i+1,k) to xab(2*m_i,k) are the isotopic heat flows 
!  xab(2*m_i+1,k) to xab(3*m_i,k) are the u2 flows 
!  Evaluate species flows 
      DO k=1,3
        CALL U_LU_BACKSUB(ab,k_order*m_i,3*mx_mi,indxb,xab(1,k))
      ENDDO
      DO i=1,m_s
        CALL U_LU_BACKSUB(ab,k_order*m_i,3*mx_mi,indxb,xabp(1,i))
        CALL U_LU_BACKSUB(ab,k_order*m_i,3*mx_mi,indxb,xabt(1,i))
      ENDDO
      DO i=1,m_s
        im=jm_s(i)
!  Source contributions 
        DO m=1,3        
          DO k=1,k_order
            IF(m.eq.2) THEN
              upar_s(k,m,i)=p_eb*rhat(k,5,i)
            ELSE
              upar_s(k,m,i)=rhat(k,m+3,i)
            ENDIF
          ENDDO
!         Response contributions           
          DO jm=1,m_i
            CALL RARRAY_ZERO(k_order,xl)
            DO l=1,k_order
              l1=jm+(l-1)*m_i
              DO k=1,k_order
                xl(k)=xl(k)-caln_ii(k,l,im,jm)*xab(l1,m)
              ENDDO
            ENDDO
            DO l=1,k_order
              DO k=1,k_order
                upar_s(k,m,i)=upar_s(k,m,i)+xl(l)*rhat(k,l,i)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!  Unit p'/p and T'/T
        DO j=1,m_s
          DO k=1,k_order
            uaip(k,j,i)=rhatp(k,j,i)
            uait(k,j,i)=rhatt(k,j,i)
          ENDDO
!         Response contributions           
          DO jm=1,m_i
            CALL RARRAY_ZERO(k_order,xl)
            DO l=1,k_order
              l1=jm+(l-1)*m_i
              DO k=1,k_order
                xl(k)=xl(k)-caln_ii(k,l,im,jm)*xabp(l1,j)
              ENDDO
            ENDDO
            DO l=1,k_order
              DO k=1,k_order
                uaip(k,j,i)=uaip(k,j,i)+xl(l)*rhat(k,l,i)
              ENDDO
            ENDDO
            CALL RARRAY_ZERO(k_order,xl)
            DO l=1,k_order
              l1=jm+(l-1)*m_i
              DO k=1,k_order
                xl(k)=xl(k)-caln_ii(k,l,im,jm)*xabt(l1,j)
              ENDDO
            ENDDO
            DO l=1,k_order
              DO k=1,k_order
                uait(k,j,i)=uait(k,j,i)+xl(l)*rhat(k,l,i)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!Currents and fluxes
      CALL RARRAY_ZERO(m_s,bsjbp_s)
      CALL RARRAY_ZERO(m_s,bsjbt_s)
      CALL RARRAY_ZERO(mx_ms*m_s,dp_ss)
      CALL RARRAY_ZERO(mx_ms*m_s,dt_ss)
      CALL RARRAY_ZERO(mx_ms*m_s,chip_ss)
      CALL RARRAY_ZERO(mx_ms*m_s,chit_ss)
      DO i=1,m_s
        im=jm_s(i)
        iz=jz_s(i)
        iza=IABS(iz)
!<J_bs.B> 
        denzc=den_iz(im,iza)*iz*z_coulomb
        p_bsjb=p_bsjb+denzc*upar_s(1,1,i)
!<J_OH.B>
        p_ohjb=p_ohjb+denzc*upar_s(1,2,i)
!<J_ex.B>
        p_exjb=p_exjb+denzc*upar_s(1,3,i)
!  Unit p'/p and T'/T
        DO j=1,m_s
          bsjbp_s(j)=bsjbp_s(j)+denzc*uaip(1,j,i)
          bsjbt_s(j)=bsjbt_s(j)+denzc*uait(1,j,i)
        ENDDO
!Fluxes
!  Banana-Plateau
        cbpa=cbp/iz
        cbpaq=cbpa*(z_j7kv*temp_i(im))
!       Unit p'/p and T'/T
        dp_ss(i,i)=dp_ss(i,i)-cbpa*ymu_s(1,1,i)*srcthp(i)
        dt_ss(i,i)=dt_ss(i,i)-cbpa*ymu_s(1,2,i)*srctht(i)
        chip_ss(i,i)=chip_ss(i,i)-cbpaq*ymu_s(2,1,i)*srcthp(i)
        chit_ss(i,i)=chit_ss(i,i)-cbpaq*ymu_s(2,2,i)*srctht(i)
        DO k=1,k_order
!         p' and T'
          gfl_s(1,i)=gfl_s(1,i)-cbpa*ymu_s(1,k,i)
     &                          *(upar_s(k,1,i)+srcth(k,i))
          qfl_s(1,i)=qfl_s(1,i)-cbpaq*ymu_s(2,k,i)
     &                          *(upar_s(k,1,i)+srcth(k,i))
!         Unit p'/p and T'/T            
          DO j=1,m_s
            dp_ss(j,i)=dp_ss(j,i)-cbpa*ymu_s(1,k,i)*uaip(k,j,i)
            dt_ss(j,i)=dt_ss(j,i)-cbpa*ymu_s(1,k,i)*uait(k,j,i)
            chip_ss(j,i)=chip_ss(j,i)-cbpaq*ymu_s(2,k,i)*uaip(k,j,i)
            chit_ss(j,i)=chit_ss(j,i)-cbpaq*ymu_s(2,k,i)*uait(k,j,i)
          ENDDO
!         <E.B> 
          gfl_s(4,i)=gfl_s(4,i)-cbpa*ymu_s(1,k,i)*upar_s(k,2,i)
          qfl_s(4,i)=qfl_s(4,i)-cbpaq*ymu_s(2,k,i)*upar_s(k,2,i)
!         External force 
          gfl_s(5,i)=gfl_s(5,i)-cbpa*ymu_s(1,k,i)*upar_s(k,3,i)
          qfl_s(5,i)=qfl_s(5,i)-cbpaq*ymu_s(2,k,i)*upar_s(k,3,i)
        ENDDO
!  Pfirsch-Schluter and classical 
!       Test particle               
        cpsa=cps*(xi_s(i)/iz)
        cpsaq=cpsa*(z_j7kv*temp_i(im))
        ccla=ccl*(xi_s(i)/iz)
        cclaq=ccla*(z_j7kv*temp_i(im))
        DO k=1,k_order
!         Pfirsch-Schluter 
          gfl_s(2,i)=gfl_s(2,i)-cpsa*calm_i(1,k,im)*srcth(k,i)
          qfl_s(2,i)=qfl_s(2,i)-cpsaq*calm_i(2,k,im)*srcth(k,i)
!         Classical 
          gfl_s(3,i)=gfl_s(3,i)+ccla*calm_i(1,k,im)*srcth(k,i)
          qfl_s(3,i)=qfl_s(3,i)+cclaq*calm_i(2,k,im)*srcth(k,i)
        ENDDO
!       Unit p'/p and T'/T 
        dp_ss(i,i)=dp_ss(i,i)-(cpsa-ccla)*calm_i(1,1,im)*srcthp(i)
        dt_ss(i,i)=dt_ss(i,i)-(cpsa-ccla)*calm_i(1,2,im)*srctht(i)
        chip_ss(i,i)=chip_ss(i,i)-(cpsaq-cclaq)*calm_i(2,1,im)*srcthp(i)
        chit_ss(i,i)=chit_ss(i,i)-(cpsaq-cclaq)*calm_i(2,2,im)*srctht(i)
!       Field particle 
        DO j=1,m_s
          jm=jm_s(j)
          cpsb=cpsa*xi_s(j)
          cpsbq=cpsb*(z_j7kv*temp_i(im))
          cclb=ccla*xi_s(j)
          cclbq=cclb*(z_j7kv*temp_i(im))
          DO k=1,k_order
!           Pfirsch-Schluter 
            gfl_s(2,i)=gfl_s(2,i)-cpsb*caln_ii(1,k,im,jm)*srcth(k,j)
            qfl_s(2,i)=qfl_s(2,i)-cpsbq*caln_ii(2,k,im,jm)*srcth(k,j)
!           Classical 
            gfl_s(3,i)=gfl_s(3,i)+cclb*caln_ii(1,k,im,jm)*srcth(k,j)
            qfl_s(3,i)=qfl_s(3,i)+cclbq*caln_ii(2,k,im,jm)*srcth(k,j)
          ENDDO
!         Unit p'/p and T'/T 
          dp_ss(j,i)=dp_ss(j,i)-(cpsb-cclb)*caln_ii(1,1,im,jm)*srcthp(j)
          dt_ss(j,i)=dt_ss(j,i)-(cpsb-cclb)*caln_ii(1,2,im,jm)*srctht(j)
          chip_ss(j,i)=chip_ss(j,i)
     &                 -(cpsbq-cclbq)*caln_ii(2,1,im,jm)*srcthp(j)
          chit_ss(j,i)=chit_ss(j,i)
     &                 -(cpsbq-cclbq)*caln_ii(2,2,im,jm)*srctht(j)
        ENDDO
      ENDDO
!Electrical resistivity
      p_etap=p_eb/p_ohjb
!Convert to diffusivities and conductivities
!  Full coefficient matrices 
      DO i=1,m_s
        im=jm_s(i)
        iza=IABS(jz_s(i))
        DO j=1,m_s
          dp_ss(j,i)=dp_ss(j,i)/den_iz(im,iza)
          dt_ss(j,i)=dt_ss(j,i)/den_iz(im,iza)
          chip_ss(j,i)=chip_ss(j,i)/den_iz(im,iza)/temp_i(im)/z_j7kv
          chit_ss(j,i)=chit_ss(j,i)/den_iz(im,iza)/temp_i(im)/z_j7kv
        ENDDO
      ENDDO
 1000 RETURN
      END
      SUBROUTINE NCLASS_K(k_config,k_banana,k_pfirsch,k_potato,m_s,jm_s,
     &                jz_s,c_potb,c_potl,p_fm,p_ft,p_ngrth,x,amu_i,
     &                temp_i,vt_i,sqz_s,ykb_s,ykp_s,ykpo_s,ykpop_s,
     &                tau_ss)
!***********************************************************************
!NCLK calculates the velocity-dependent neoclassical viscosity
!  coefficients, K
!References:                                                      
!  Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079                     
!  Kessel, Nucl Fusion 34 (1994) 1221                              
!  Shaing, Yokoyama, Wakatani, Hsu, Phys Plasmas 3 (1996) 965      
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230                     
!  W.A.Houlberg 6/99                          
!Input:
!  k_config -configuration
!       =1 stellarator
!       =else tokamak                   
!  k_banana-option to include banana viscosity (-)
!          =0 off
!          =else on
!  k_pfirsch-option to include Pfirsch-Schluter viscosity (-)
!          =0 off
!          =else on
!  k_potato-option to include potato orbits (-)
!          =0 off
!          =else on
!  m_s-number of species (1<m_s<mx_ms+1)
!  jm_s(s)-isotope number of s (-)
!  jz_s(s)-charge state of s (-)
!  c_potb-kappa(0)*Bt(0)/[2*q(0)**2] (T)
!  c_potl-q(0)*R(0) (m)
!  p_fm(3)-poloidal moments of geometric factor for PS viscosity (-)
!  p_ft-trapped fraction (-)
!  p_ngrth-<n.grad(Theta)> (1/m)
!  x-velocity normalized to thermal velocity v/(2kT/m)**0.5 (-)
!  amu_i(i)-atomic mass number of i (-)
!  temp_i(i)-temperature of i (keV)
!  vt_i(i)-thermal velocity of i (m/s)
!  sqz_s(s)-orbit squeezing factor for s (-)
!Output:
!  ykb_s(s)-banana viscosity for s (kg/m**3/s)
!  ykp_s(s)-Pfirsch-Schluter viscosity for s (kg/m**3/s)
!  ykpo_s(s)-potato viscosity for s (kg/m**3/s)
!  ykpop_s(s)-potato-plateau viscosity for s (kg/m**3/s)
!  tau_ss(s1,s2)-90 degree scattering time of s1 on s2 (s)
!***********************************************************************
      IMPLICIT NONE
      INCLUDE 'pamx_mi.inc'
      INCLUDE 'pamx_ms.inc'
!Declaration of input variables
      INTEGER        k_banana,                k_pfirsch,
     &               k_potato,                k_config,
     &               m_s
      INTEGER        jm_s(mx_ms),             jz_s(mx_ms)
      REAL           c_potb,                  c_potl,
     &               p_fm(3),                 p_ft,
     &               p_ngrth,                 x
      REAL           amu_i(mx_mi),            temp_i(mx_mi),
     &               vt_i(mx_mi)
      REAL           sqz_s(mx_ms)
!Declaration of output variables
      REAL           ykb_s(mx_ms),            ykp_s(mx_ms),
     &               ykpo_s(mx_ms),           ykpop_s(mx_ms)
      REAL           tau_ss(mx_ms,mx_ms)
!Declaration of local variables
      INTEGER        i,                       im,
     &               iz,                      m
      REAL           c1,                      c2,
     &               c3,                      c4,
     &               c5,                      c6,
     &               z_coulomb,               z_pi,
     &               z_protonmass
      REAL           ynud_s(mx_ms),           ynut_s(mx_ms),
     &               ynutis(3,mx_ms)
!Initialization
!  Physical and conversion constants
      z_coulomb=1.6022e-19
      z_pi=ACOS(-1.0)
      z_protonmass=1.6726e-27
!  Zero out arrays
      CALL RARRAY_ZERO(m_s,ykb_s)
      CALL RARRAY_ZERO(m_s,ykp_s)
      CALL RARRAY_ZERO(m_s,ykpo_s)
      CALL RARRAY_ZERO(m_s,ykpop_s)
!  Get collisional frequencies
      CALL NCLASS_NU(m_s,jm_s,p_ngrth,x,temp_i,vt_i,tau_ss,ynud_s,
     &               ynut_s,ynutis)
!  Set velocity dependent viscosities (K's)
      c1=1.5*x**2
      c2=3.0*2.19/(2.0**1.5)*x**(1.0/3.0)
      IF(k_potato.ne.0) THEN
        c3=3.0*z_pi/(64.0*2.0**0.33333)/ABS(c_potl)
      ENDIF
      DO i=1,m_s
        im=jm_s(i)
        iz=jz_s(i)
        IF(k_banana.ne.0) THEN
!         Provide cutoff to eliminate failure at unity trapped fraction
!         At A=>1 viscosity will go over to Pfirsch-Schluter value
          c4=1.0-p_ft
          IF(c4.lt.1.0e-3) c4=1.0e-3
          IF (k_config .eq.1) THEN
            c5=2.48*p_fm(1)/vt_i(im)/x
            c6=1.96*p_fm(1)*p_fm(2)/vt_i(im)/x
            ykb_s(i)=p_ft/c4*ynud_s(i)/(1.0+c5*ynud_s(i))/
     &                                 (1.0+c6*ynud_s(i))
          ELSE  
            ykb_s(i)=p_ft/c4/sqz_s(i)**1.5*ynud_s(i)
          ENDIF  
        ENDIF
        IF(k_pfirsch.ne.0) THEN
          DO m=1,3
            ykp_s(i)=ykp_s(i)+c1*vt_i(im)**2*p_fm(m)
     &                        *(ynutis(m,i)/ynut_s(i))
          ENDDO
        ENDIF
        IF(k_potato.ne.0) THEN
          c4=ABS(amu_i(im)*z_protonmass*vt_i(im)
     &       /(iz*z_coulomb*c_potb*c_potl))
          ykpo_s(i)=c2*c4**(1.0/3.0)*ynud_s(i)/sqz_s(i)**(5.0/3.0)
          ykpop_s(i)=c3*vt_i(im)*c4**(4.0/3.0)
        ENDIF
      ENDDO         
      RETURN
      END
      SUBROUTINE NCLASS_MN(k_order,m_i,amu_i,temp_i,capm_ii,capn_ii)
!***********************************************************************
!NCLASS_MN calculates the k_order*korder matrix of test particle (M) and
!  field particle (N) coefficients of the collision operator using the
!  Laguerre polynomials of order 3/2 as basis functions for each
!  isotopic species combination
!References:                                         
!  Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079 (HS81)
!  Hirshman, Phys Fluids 20 (1977) 589
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  W.A.Houlberg 6/99
!Input:
!  k_order-order of v moments to be solved (-)
!       =2 u and q
!       =3 u, q, and u2
!       =else error
!  m_i-number of isotopes (1<m_i<mx_mi+1)
!  amu_i(a)-atomic mass of a (-)
!  temp_i(a)-temperature of a (keV)
!Output:
!  capm_ii(3,3,i1,i2)-test part (tp) friction matrix for i1 on i2 (-)
!  capn_ii(3,3,i1,i2)-field part (fp) friction matrix for i1 on i2 (-)
!Comments:                                                     
!The indices on the M and N matrices are one greater than the notation
!  in the review article so as to avoid 0 as an index
!***********************************************************************
      IMPLICIT NONE
      INCLUDE 'pamx_mi.inc'
!Declaration of input variables
      INTEGER        k_order,                 m_i
      REAL           amu_i(mx_mi),            temp_i(mx_mi)
!Declaration of output variables
      REAL           capm_ii(3,3,mx_mi,mx_mi), capn_ii(3,3,mx_mi,mx_mi)
!Declaration of local variables
      INTEGER        im,                      jm
      REAL           xab,                     xab2,                    
     &               xmab,                    xtab,
     &               yab32,                   yab52,
     &               yab72,                   yab92
!Loop over isotope a
      DO im=1,m_i
!Loop over isotope b
        DO jm=1,m_i
!         Ratio of masses
          xmab=amu_i(im)/amu_i(jm)
!         Ratio of temperatures
          xtab=temp_i(im)/temp_i(jm)
!         Ratio of thermal velocities, vtb/vta
          xab=SQRT(xmab/xtab)
!  Elements of M
          xab2=xab**2
          yab32=(1.0+xab2)*SQRT(1.0+xab2)
          yab52=(1.0+xab2)*yab32
          IF(k_order.eq.3) THEN
            yab72=(1.0+xab2)*yab52
            yab92=(1.0+xab2)*yab72
          ENDIF
!         Eqn 4.11 for M00 (HS81)
          capm_ii(1,1,im,jm)=-(1.0+xmab)/yab32
!         Eqn 4.12 for M01 (HS81)
          capm_ii(1,2,im,jm)=3.0/2.0*(1.0+xmab)/yab52
!         Eqn 4.8 for M10 (HS81)
          capm_ii(2,1,im,jm)=capm_ii(1,2,im,jm)
!         Eqn 4.13 for M11 (HS81)
          capm_ii(2,2,im,jm)=-(13.0/4.0+xab2*(4.0+xab2*15.0/2.0))/yab52
          IF(k_order.eq.3) THEN
!           Eqn 4.15 for M02 (HS81)
            capm_ii(1,3,im,jm)=-15.0/8.0*(1.0+xmab)/yab72          
!           Eqn 4.16 for M12 (HS81)
            capm_ii(2,3,im,jm)=(69.0/16.0+xab2*(6.0+xab2*63.0/4.0))
     &                         /yab72
!           Eqn 4.8 for M20 (HS81)
            capm_ii(3,1,im,jm)=capm_ii(1,3,im,jm)
!           Eqn 4.8 for M21 (HS81)
            capm_ii(3,2,im,jm)=capm_ii(2,3,im,jm)
!           Eqn 5.21 for M22 (HS81)
            capm_ii(3,3,im,jm)=-(433.0/64.0+xab2*(17.0+xab2*(459.0/8.0
     &                         +xab2*(28.0+xab2*175.0/8.0))))/yab92
          ENDIF     
!  Elements of N
!         Momentum conservation, Eqn 4.11 for N00 (HS81)
          capn_ii(1,1,im,jm)=-capm_ii(1,1,im,jm)
!         Eqn 4.9 and 4.12 for N01 (HS81)
          capn_ii(1,2,im,jm)=-xab2*capm_ii(1,2,im,jm)
!         Momentum conservation, Eqn 4.12 for N10 (HS81)
          capn_ii(2,1,im,jm)=-capm_ii(2,1,im,jm)
!         Eqn 4.14 for N11 (HS81)	- corrected rhs
          capn_ii(2,2,im,jm)=(27.0/4.0)*SQRT(xtab)*xab2/yab52
          IF(k_order.eq.3) THEN
!           Eqn 4.15 for N02 (HS81) - corrected rhs by Ta/Tb
            capn_ii(1,3,im,jm)=-xab2**2*capm_ii(1,3,im,jm)
!           Eqn 4.17 for N12 (HS81)
            capn_ii(2,3,im,jm)=-225.0/16.0*xtab*xab2**2/yab72
!           Momentum conservation for N20 (HS81)
            capn_ii(3,1,im,jm)=-capm_ii(3,1,im,jm)
!           Eqn 4.9 and 4.17 for N21 (HS81)
            capn_ii(3,2,im,jm)=-225.0/16.0*xab2**2/yab72
!           Eqn 5.22 for N22 (HS81) 
            capn_ii(3,3,im,jm)=2625.0/64.0*xtab*xab2**2/yab92
          ENDIF
        ENDDO   
      ENDDO   
      RETURN
      END
      SUBROUTINE NCLASS_MU(k_config,k_order,k_banana,k_pfirsch,k_potato,
     &                  m_s,jm_s,jz_s,c_potb,c_potl,p_fm,p_ft,p_ngrth,
     &                  amu_i,temp_i,vt_i,den_iz,sqz_s,ymu_s,tau_ss)
!***********************************************************************
!NCLASS_MU calculates the k_order*k_order matrix of neoclassical
!  viscosities by integrating the velocity-dependent banana and Pfirsch-
!  Schluter contributions
!References:                                             
!  Shaing, Yokoyama, Wakatani, Hsu, Phys Plasmas 3 (1996) 965
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  W.A.Houlberg 6/99
!Input:
!  k_config -configuration
!       =1 stellarator
!       =else tokamak             
!  k_order-order of v moments to be solved (-)
!       =2 u and q
!       =3 u, q, and u2
!       =else error
!  k_banana-option to include banana viscosity (-)
!          =0 off
!          =else on
!  k_pfirsch-option to include Pfirsch-Schluter viscosity (-)
!          =0 off
!          =else on
!  k_potato-option to include potato orbits (-)
!          =0 off
!          =else on
!  m_s-number of species (1<m_s<mx_ms+1)
!  jm_s(s)-isotope number of species s (-)
!  jz_s(s)-charge state of s (-)
!  c_potb-kappa(0)*Bt(0)/[2*q(0)**2] (T)
!  c_potl-q(0)*R(0) (m)
!  p_fm(3)-poloidal moments of geometric factor for PS viscosity (-)
!  p_ft-trapped fraction (-)
!  p_ngrth-<n.grad(Theta)> (1/m)
!  amu_i(i)-atomic mass number of i (-)
!  temp_i(i)-temperature of i (keV)
!  vt_i(i)-thermal velocity of i (m/s)
!  den_iz(i,z)-density of i,z (/m**3)
!  sqz_s(s)-orbit squeezing factor for s (-)
!Output:
!  ymu_s(s)-normalized viscosity for s (kg/m**3/s)
!  tau_ss(s1,s2)-90 degree scattering time of s1 on s2 (s)
!***********************************************************************
      IMPLICIT NONE
      INCLUDE 'pamx_mi.inc'
      INCLUDE 'pamx_ms.inc'
      INCLUDE 'pamx_mz.inc'
!Declaration of input variables
      INTEGER        k_banana,                k_order,
     &               k_pfirsch,               k_potato,
     &               k_config,                m_s
      INTEGER        jm_s(mx_ms),             jz_s(mx_ms)
      REAL           c_potb,                  c_potl,
     &               p_ft,                    p_ngrth,
     &               p_fm(3)
      REAL           amu_i(mx_mi),            temp_i(mx_mi),
     &               vt_i(mx_mi)
      REAL           den_iz(mx_mi,mx_mz)
      REAL           sqz_s(mx_ms)
!Declaration of output variables
      REAL           ymu_s(3,3,mx_ms)
      REAL           tau_ss(mx_ms,mx_ms)
!Declaration of local variables
      INTEGER        mpnts
      PARAMETER     (mpnts=13)
      INTEGER        i,                       init,
     &               im,                      iza,
     &               k,                       l,
     &               m
      REAL           bmax,                    c1,
     &               c2,                      ewt,
     &               expmx2,                  h,
     &               x2,                      x4,
     &               x6,                      x8,
     &               x10,                     x12,
     &               xk,                      xx,
     &               z_pi,                    z_protonmass
      REAL           x(mpnts),                w(mpnts,5)
      REAL           ykb_s(mx_ms),            ykp_s(mx_ms),
     &               ykpo_s(mx_ms),           ykpop_s(mx_ms),
     &               ymubs(3,3,mx_ms),        ymubps(3,3,mx_ms),
     &               ymupps(3,3,mx_ms),       ymupos(3,3,mx_ms)
      REAL           dum(3)
      SAVE           init
      SAVE           c1,                      h,
     &               x,                       w,
     &               z_protonmass
      DATA bmax/     3.2/
      DATA init/     0/
!Initialization
      IF(init.eq.0) THEN
!  Physical and conversion constants
        z_pi=ACOS(-1.0)
        z_protonmass=1.6726e-27
!  Set integration points and weights
        h=bmax/(mpnts-1)
        x(1)=0.0
        w(1,1)=0.0
        w(1,2)=0.0
        w(1,3)=0.0
        w(1,4)=0.0
        w(1,5)=0.0
        DO m=2,mpnts
          x(m)=h*(m-1)
          x2=x(m)*x(m)
          expmx2=EXP(-x2)
          x4=x2*x2
          w(m,1)=x4*expmx2
          x6=x4*x2
          w(m,2)=x6*expmx2
          x8=x4*x4
          w(m,3)=x8*expmx2
          x10=x4*x6
          w(m,4)=x10*expmx2
          x12=x6*x6
          w(m,5)=x12*expmx2
        ENDDO   
        c1=8.0/3.0/SQRT(z_pi)*h
        init=1
      ENDIF
      CALL RARRAY_ZERO(9*m_s,ymu_s)
      CALL RARRAY_ZERO(9*m_s,ymubs)
      CALL RARRAY_ZERO(9*m_s,ymubps)
      CALL RARRAY_ZERO(9*m_s,ymupos)
      CALL RARRAY_ZERO(9*m_s,ymupps)
      IF((k_banana.ne.0).or.(k_pfirsch.ne.0)) THEN
!Loop over grid (first node has null value)
        DO m=2,mpnts
          IF(m.eq.mpnts) THEN
!           Use half weight for end point
            ewt=0.5
          ELSE
!           Use full weight
            ewt=1.0
          ENDIF
          xx=x(m)
!  Get velocity-dependent k values        
          CALL NCLASS_K(k_config,k_banana,k_pfirsch,k_potato,m_s,jm_s,
     &              jz_s,c_potb,c_potl,p_fm,p_ft,p_ngrth,xx,amu_i,
     &              temp_i,vt_i,sqz_s,ykb_s,ykp_s,ykpo_s,ykpop_s,tau_ss)
!  Loop over species
          DO i=1,m_s
            im=jm_s(i)
            iza=IABS(jz_s(i))
            c2=c1*ewt*den_iz(im,iza)*amu_i(im)*z_protonmass
            dum(1)=c2*w(m,1)
            dum(2)=c2*(w(m,2)-5.0/2.0*w(m,1))
            dum(3)=c2*(w(m,3)-5.0*w(m,2)+(25.0/4.0)*w(m,1))
            IF(k_banana.eq.0) THEN
              xk=ykp_s(i)
            ELSEIF(k_pfirsch.eq.0) THEN
              xk=ykb_s(i)
            ELSE
              xk=ykb_s(i)*ykp_s(i)/(ykb_s(i)+ykp_s(i))
            ENDIF
            ymubs(1,1,i)=ymubs(1,1,i)+ykb_s(i)*dum(1)
            ymubs(1,2,i)=ymubs(1,2,i)+ykb_s(i)*dum(2)
            ymubs(2,2,i)=ymubs(2,2,i)+ykb_s(i)*dum(3)
            ymubps(1,1,i)=ymubps(1,1,i)+xk*dum(1)
            ymubps(1,2,i)=ymubps(1,2,i)+xk*dum(2)
            ymubps(2,2,i)=ymubps(2,2,i)+xk*dum(3)
            IF(k_potato.ne.0) THEN
              ymupos(1,1,i)=ymupos(1,1,i)+ykpo_s(i)*dum(1)
              ymupos(1,2,i)=ymupos(1,2,i)+ykpo_s(i)*dum(2)
              ymupos(2,2,i)=ymupos(2,2,i)+ykpo_s(i)*dum(3)
              xk=ykpo_s(i)*ykpop_s(i)/(ykpo_s(i)+ykpop_s(i))
              ymupps(1,1,i)=ymupps(1,1,i)+xk*dum(1)
              ymupps(1,2,i)=ymupps(1,2,i)+xk*dum(2)
              ymupps(2,2,i)=ymupps(2,2,i)+xk*dum(3)
            ENDIF
            IF(k_order.eq.3) THEN
              dum(1)=c2*((1.0/2.0)*w(m,3)-(7.0/2.0)*w(m,2)
     &               +(35.0/8.0)*w(m,1))
              dum(2)=c2*((1.0/2.0)*w(m,4)-(19.0/4.0)*w(m,3)
     &               +(105.0/8.0)*w(m,2)-(175.0/16.0)*w(m,1))
              dum(3)=c2*((1.0/4.0)*w(m,5)-(7.0/2.0)*w(m,4)
     &               +(133.0/8.0)*w(m,3)-(245.0/8.0)*w(m,2)
     &               +(1225.0/64.0)*w(m,1))
              ymubs(1,3,i)=ymubs(1,3,i)+ykb_s(i)*dum(1)
              ymubs(2,3,i)=ymubs(2,3,i)+ykb_s(i)*dum(2)
              ymubs(3,3,i)=ymubs(3,3,i)+ykb_s(i)*dum(3)
              xk=ykb_s(i)*ykp_s(i)/(ykb_s(i)+ykp_s(i))
              ymubps(1,3,i)=ymubps(1,3,i)+xk*dum(1)
              ymubps(2,3,i)=ymubps(2,3,i)+xk*dum(2)
              ymubps(3,3,i)=ymubps(3,3,i)+xk*dum(3)
              IF(k_potato.ne.0) THEN
                ymupos(1,3,i)=ymupos(1,3,i)+ykpo_s(i)*dum(1)
                ymupos(2,3,i)=ymupos(2,3,i)+ykpo_s(i)*dum(2)
                ymupos(3,3,i)=ymupos(3,3,i)+ykpo_s(i)*dum(3)
                xk=ykpo_s(i)*ykpop_s(i)/(ykpo_s(i)+ykpop_s(i))
                ymupps(1,3,i)=ymupps(1,3,i)+xk*dum(1)
                ymupps(2,3,i)=ymupps(2,3,i)+xk*dum(2)
                ymupps(3,3,i)=ymupps(3,3,i)+xk*dum(3)
              ENDIF
            ENDIF
          ENDDO    
        ENDDO
!Load net viscosity
        DO i=1,m_s
          DO l=1,k_order
            DO k=1,l
              IF(k_potato.eq.0) THEN
!               Banana Pfirsch-Schluter
                ymu_s(k,l,i)=ymubps(k,l,i)
              ELSE
!               Banana Pfirsch-Schluter plus potato potato-plateau
                ymu_s(k,l,i)=(ymupos(k,l,i)**3*ymupps(k,l,i)
     &                       +ymubs(k,l,i)**3*ymubps(k,l,i))
     &                       /(ymupos(k,l,i)**3+ymubs(k,l,i)**3)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!Fill viscosity matrix using symmetry
        DO i=1,m_s
          DO l=1,k_order-1
            DO k=l+1,k_order
                ymu_s(k,l,i)=ymu_s(l,k,i)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END
      SUBROUTINE NCLASS_NU(m_s,jm_s,p_ngrth,x,temp_i,vt_i,tau_ss,ynud_s,
     &                     ynut_s,ynutis)
!***********************************************************************
!NCLASS_NU calculates the velocity dependent pitch angle diffusion and
!  anisotropy relaxation rates, nu_D, nu_T, and nu_T*I_Rm
!References:                                                     
!  Hirshman, Sigmar, Phys Fluids 19 (1976) 1532
!  Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079
!  Shaing, Yokoyama, Wakatani, Hsu, Phys Plasmas 3 (1996) 965
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  W.A.Houlberg 6/99
!Input:
!  m_s-number of species (1<m_s<mx_ms+1)
!  jm_s(s)-isotope number of s (-)
!  p_ngrth-<n.grad(Theta)> (1/m)
!  x-velocity normalized to thermal velocity v/(2kT/m)**0.5 (-)
!  temp_i(i)-temperature of i (keV)
!  vt_i(i)-thermal velocity of i (m/s)
!  tau_ss(s1,s2)-90 degree scattering time of s1 on s2 (s)
!Output:
!  ynud_s(s)-pitch angle diffusion rate for s (/s)
!  ynut_s(s)-anisotropy relaxation rate for s (/s)
!  ynutis(3,s)-PS anisotropy relaxation rates for s (/s)
!***********************************************************************
      IMPLICIT NONE
      INCLUDE 'pamx_mi.inc'
      INCLUDE 'pamx_ms.inc'
!Declaration of input variables
      INTEGER        m_s
      INTEGER        jm_s(mx_ms)       
      REAL           p_ngrth,                  x
      REAl           temp_i(mx_mi),            vt_i(mx_mi)
      REAL           tau_ss(mx_ms,mx_ms)
!Declaration of output variables
      REAL           ynud_s(mx_ms),            ynut_s(mx_ms),
     &               ynutis(3,m_s)
!Declaration of local variables
      INTEGER        i,                       im,
     &               j,                       jm,
     &               m
      REAL           c1,                      c2,
     &               c3,                      g,
     &               phi,                     z_pi
!Declaration of external functions
      REAL           U_ERF
!Initializaton
!  Physical and conversion constants
      z_pi=ACOS(-1.0)
!  Zero out arrays
      CALL RARRAY_ZERO(m_s,ynud_s)
      CALL RARRAY_ZERO(m_s,ynut_s)
      CALL RARRAY_ZERO(3*m_s,ynutis)
      DO i=1,m_s
        im=jm_s(i)
!Calculate nu_D and nu_T
        DO j=1,m_s
          jm=jm_s(j)
          c1=vt_i(jm)/vt_i(im)
          c2=x/c1
          phi=U_ERF(c2)
          g=(phi-c2*(2.0/SQRT(z_pi))*EXP(-c2**2))/(2.0*c2**2)
          ynud_s(i)=ynud_s(i)+(3.0*SQRT(z_pi)/4.0)*(phi-g)/x**3
     &                        /tau_ss(i,j)
          ynut_s(i)=ynut_s(i)+((3.0*SQRT(z_pi)/4.0)*((phi-3.0*g)/x**3
     &                       +4.0*(temp_i(im)/temp_i(jm)
     &                       +1.0/c1**2)*g/x))/tau_ss(i,j)
        ENDDO   
! Calculate nu_T*I_m
        DO m=1,3
          IF(ABS(p_ngrth).gt.0.0) THEN
            c1=x*vt_i(im)*m*p_ngrth
            c2=(ynut_s(i)/c1)**2
            IF(c2.gt.9.0) THEN
              ynutis(m,i)=0.4
            ELSE
              c3=ynut_s(i)/c1*ATAN(c1/ynut_s(i))
              ynutis(m,i)=0.5*c3+c2*(3.0*(c3-0.5)+c2*4.5*(c3-1.0))
            ENDIF
          ELSE
            ynutis(m,i)=0.4
          ENDIF
        ENDDO   
      ENDDO
      RETURN
      END
      SUBROUTINE NCLASS_TAU(m_i,m_s,jm_s,jz_s,amu_i,temp_i,vt_i,den_iz,
     &                      amnt_ii,tau_ss)
!***********************************************************************
!NCLASS_TAU calculates the collision time for 90 degree scattering
!  assuming a common value for the Coulomb logarithm for each isotope
!References:                                      
!  Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  W.A.Houlberg 6/99
!Input:
!  m_i-number of isotopes (1<m_i<mx_mi+1)
!  m_s-number of (1<m_s<mx_ms+1)
!  jm_s(s)-isotope number of s (-)
!  jz_s(s)-charge state of s (-)
!  amu_i(i)-atomic mass number of i (-)
!  temp_i(i)-temperature of i (keV)
!  vt_i(i)-thermal velocity of i (m/s)
!  den_iz(i,z)-density of i,z (/m**3)
!Output:
!  amnt_ii(s1,s2)-eff relaxation rate for s1 on s2 (kg/m**3/s)
!  tau_ss(s1,s2)-90 degree scattering time of s1 on s2 (s)
!***********************************************************************
      IMPLICIT NONE
      INCLUDE 'pamx_mi.inc'
      INCLUDE 'pamx_ms.inc'
      INCLUDE 'pamx_mz.inc'
!Declaration of input variables
      INTEGER        m_i,                      m_s
      INTEGER        jm_s(mx_ms),              jz_s(mx_ms)
      REAL           amu_i(mx_mi),             temp_i(mx_mi),
     &               vt_i(mx_mi)
      REAL           den_iz(mx_mi,mx_mz)
!Declaration of output variables
      REAL           amnt_ii(mx_mi,mx_mi)
      REAL           tau_ss(mx_ms,mx_ms)
!Declaration of local variables
      INTEGER        i,                       im,
     &               iz,                      iza,
     &               j,                       jm,
     &               jz,                      jza
      REAL           c1,                      c2,
     &               clnab,                   z_coulomb,
     &               z_epsilon0,              z_pi,
     &               z_protonmass
      REAL           xlnab(mx_mi,mx_mi),      xn(mx_mi),
     &               xnz(mx_mi),              xz(mx_mi),
     &               xnz2(mx_mi)
!Initialization
!  Physical and conversion constants
      z_coulomb=1.6022e-19
      z_epsilon0=8.8542e-12
      z_pi=ACOS(-1.0)
      z_protonmass=1.6726e-27
!  Zero out arrays
      CALL RARRAY_ZERO(mx_mi*m_i,amnt_ii)
      CALL RARRAY_ZERO(mx_ms*m_s,tau_ss)
      CALL RARRAY_ZERO(m_i,xn)
      CALL RARRAY_ZERO(m_i,xnz)
      CALL RARRAY_ZERO(m_i,xnz2)
      c1=4.0/3.0/SQRT(z_pi)*4.0*z_pi*(z_coulomb
     &   /(4.0*z_pi*z_epsilon0))**2*(z_coulomb/z_protonmass)**2
!Coulomb logarithm
      DO i=1,m_s
        im=jm_s(i)
        iz=jz_s(i)
        iza=IABS(iz)
        IF(iz.lt.0) THEN
!         Electrons
          clnab=37.8-ALOG(SQRT(den_iz(im,iza))/temp_i(im))
        ENDIF
        xn(im)=xn(im)+den_iz(im,iza)
        xnz(im)=xnz(im)+iz*den_iz(im,iza)
        xnz2(im)=xnz2(im)+iz**2*den_iz(im,iza)
      ENDDO
      DO i=1,m_i
        xz(i)=xnz(i)/xn(i)
      ENDDO
      DO im=1,m_i
        iz=jz_s(im)
        DO jm=1,m_i
          jz=jz_s(jm)
          IF(iz.lt.0.or.jz.lt.0) THEN
!  Electrons
            xlnab(im,jm)=clnab
          ELSE
!  Ions
            xlnab(im,jm)=40.3
     &                   -ALOG(xz(im)*xz(jm)*(amu_i(im)+amu_i(jm))
     &                   /(amu_i(im)*temp_i(jm)+amu_i(jm)*temp_i(im))
     &                   *SQRT(xnz2(im)/temp_i(im)+xnz2(jm)/temp_i(jm)))
          ENDIF
        ENDDO
      ENDDO
!Collision times and mass density weighted collision rates
      DO i=1,m_s
        im=jm_s(i)
        iz=jz_s(i)
        iza=IABS(iz)
        c2=(vt_i(im)**3)*amu_i(im)**2/c1
        DO j=1,m_s
          jm=jm_s(j)
          jz=jz_s(j)
          jza=IABS(jz)
          tau_ss(i,j)=c2/xlnab(im,jm)/iz**2
     &                /(den_iz(jm,jza)*jz**2)
          amnt_ii(im,jm)=amnt_ii(im,jm)+amu_i(im)*z_protonmass
     &                   *den_iz(im,iza)/tau_ss(i,j)
        ENDDO          
      ENDDO   
      RETURN
      END
