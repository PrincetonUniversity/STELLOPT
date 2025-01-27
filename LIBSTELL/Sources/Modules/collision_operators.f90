!-----------------------------------------------------------------------
!     Module:        collision_operators
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          12/18/2024
!     Description:   This module contains routines for calculating
!                    various collision based quantities.
!-----------------------------------------------------------------------
      MODULE collision_operators
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, PRIVATE, PARAMETER :: &
                              electron_mass   = 9.1093837139D-31 !m_e
      DOUBLE PRECISION, PRIVATE, PARAMETER :: &
                              electron_charge = 1.602176634D-19 !e_c
      DOUBLE PRECISION, PRIVATE, PARAMETER :: &
                              sqrt_pi         = SQRT(4.0 * ATAN(1.0))   !pi^(1/2)
      DOUBLE PRECISION, PRIVATE, PARAMETER :: &
                              inv_dalton      = 6.02214076208E+26 ! 1./AMU [1/kg]
      DOUBLE PRECISION, PRIVATE, PARAMETER :: &
                              hbar            = 1.054571817E+34 ! hbar [J.s]
      DOUBLE PRECISION, PRIVATE :: coulomb_factor, fact_crit, &
                                   fact_crit_weiland
      
!-----------------------------------------------------------------------
!     Function Interfaces
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!     Subroutines and Functions
!-----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE SET_COULOMB_FACTOR(mass,z,mass_plasma)
         !--------------------------------------------------------------
         !   MASS         Particle mass [kg]
         !   Z            Particle charge number [e]
         !   MASS_PLASMA  Plasma effective mass [kg]
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(in) :: mass,z,mass_plasma
         coulomb_factor = z * ( mass + mass_plasma ) / &
                              ( mass * mass_plasma * 6.02214076208D+26 )
         RETURN
      END SUBROUTINE SET_COULOMB_FACTOR

      SUBROUTINE SET_CRIT_FACTOR(z_plasma,mass_plasma)
         !--------------------------------------------------------------
         !   Z_PLASMA     Particle charge number [e]
         !   MASS_PLASMA  Plasma effective mass [kg]
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(in) :: z_plasma,mass_plasma
         fact_crit = SQRT( 2.0 * electron_charge /mass_plasma) * ( 0.75 * sqrt_pi * &
                            SQRT( mass_plasma / electron_mass ) )**(1.0/3.0)
         fact_crit_weiland=fact_crit*z_plasma** (1.0/3.0) ! WESSON PG 226 5.4.9
         RETURN
      END SUBROUTINE SET_CRIT_FACTOR

      FUNCTION COULOMB_LOG_NRL_EE(ne,te) RESULT(result)
         !--------------------------------------------------------------
         !   NE      Electron Density [m^-3]
         !   TE      Electron Temperature [eV]
         !   RESULT  Coulomb Logarithm Thermal electron-electron
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(in) :: ne,te
         DOUBLE PRECISION :: result
         result = 23.5 - LOG( SQRT( te * 1D-6 ) * te**-1.25 ) &
                - SQRT( 1D-5 + 0.0625 * ( LOG(te) - 2.0 )**2.0 )
         RETURN
      END FUNCTION COULOMB_LOG_NRL_EE

      FUNCTION COULOMB_LOG_NRL_IEA(Z,ne,te) RESULT(result)
         !--------------------------------------------------------------
         !           For Ti*me/mi < Te < 10*Z*Z [eV]
         !   Z       Ion charge number [e]
         !   NE      Electron Density [m^-3]
         !   TE      Electron Temperature [eV]
         !   RESULT  Coulomb Logarithm Thermal ion-electron
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(in) :: Z,ne,te
         DOUBLE PRECISION :: result
         result = 23.0 - LOG( SQRT( ne * 1.0E-6 ) * Z * te**(-1.50))
         RETURN
      END FUNCTION COULOMB_LOG_NRL_IEA

      FUNCTION COULOMB_LOG_NRL_IEB(ne,te) RESULT(result)
         !--------------------------------------------------------------
         !           For Ti*me/mi < 10*Z*Z [eV] < Te
         !   NE      Electron Density [m^-3]
         !   TE      Electron Temperature [eV]
         !   RESULT  Coulomb Logarithm Thermal ion-electron
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(in) :: ne,te
         DOUBLE PRECISION :: result
         result = 24.0 - LOG( SQRT( ne * 1.0E-6 ) / te)
         RETURN
      END FUNCTION COULOMB_LOG_NRL_IEB

      FUNCTION COULOMB_LOG_NRL_IEC(mass,Z,ni,ti) RESULT(result)
         !--------------------------------------------------------------
         !           For Te < Ti*me/mi
         !   Z       Ion mass [kg]
         !   Z       Ion charge number [e]
         !   NI      Ion Density [m^-3]
         !   TI      Ion Temperature [eV]
         !   RESULT  Coulomb Logarithm Thermal ion-electron
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(in) :: mass,Z,ni,ti
         DOUBLE PRECISION :: result
         result = 16.0 - LOG( SQRT( ni * 1.0E-6 ) * Ti**(-1.5) * Z * Z * mass * inv_dalton)
         RETURN
      END FUNCTION COULOMB_LOG_NRL_IEC

      FUNCTION COULOMB_LOG_NRL_II(mass1,z1,ni1,ti1,mass2,z2,ni2,ti2) RESULT(result)
         !--------------------------------------------------------------
         !   MASS1   First Ion mass [kg]
         !   Z1      First Ion charge number [e]
         !   NI1     First Ion Density [m^-3]
         !   TI1     First Ion Temperature [eV]
         !   MASS2   Second Ion mass [kg]
         !   Z2      Second Ion charge number [e]
         !   NI2     Second Ion Density [m^-3]
         !   TI2     Second Ion Temperature [eV]
         !   RESULT  Coulomb Logarithm Thermal ion-ion
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(in) :: mass1,z1,ni1,ti1
         DOUBLE PRECISION, INTENT(in) :: mass2,z2,ni2,ti2
         DOUBLE PRECISION :: result
         result = 23.0 - LOG( ( z1 * z2 * ( mass1 + mass2 ) ) * &
                    ( ni1 * z1 * z1 / ti1 + ni2 * z2 * z2 / ti2) * &
                    1.0E-6 / ( mass1 * ti2 + mass2 * ti1 ) )
         RETURN
      END FUNCTION COULOMB_LOG_NRL_II

      FUNCTION COULOMB_LOG_NRL_GENERAL(mass1,Z1,mass2,Z2,ld,u) RESULT(result)
         !--------------------------------------------------------------
         !   MASS1   First Species mass [kg]
         !   Z1      First Species charge number [e]
         !   MASS2   Second Species mass [kg]
         !   Z2      Second Species charge number [e]
         !   LD      Plasma Debye Length [m]
         !   U       Relative velocities (|v1-v2|)
         !   RESULT  Coulomb Logarithm species-species
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(in) :: mass1,z1,mass2,z2,ld,u
         DOUBLE PRECISION :: result
         DOUBLE PRECISION :: rmina,rminb
         rmina   =  1.0 / mass2 + 1.0 / mass1
         rminb   =  0.5 * hbar * rmina / u
         rmina   =  rmina*electron_charge*electron_charge*Z1*Z2/( u * u )
         result  = LOG( LD / MAX(rmina,rminb) )
      END FUNCTION COULOMB_LOG_NRL_GENERAL

      FUNCTION COULOMB_LOG_NRL_COUNTERSTREAM(ne,te,vbeta,Zeff) RESULT(result)
         !--------------------------------------------------------------
         !   NE      Electron Density [m^-3]
         !   TE      Electron Temperature [eV]
         !   VBETA   Normalized Particle Velocity (v/c)
         !   ZEFF    Plasma Effective Charge [arb]
         !   RESULT  Coulomb Logarithm Counterstreaming ions
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(in) :: ne,te,vbeta,Zeff
         DOUBLE PRECISION :: result
         result = 43.0 - log(Zeff*coulomb_factor*sqrt(ne*1D-6/te)/(vbeta*vbeta))
         RETURN
      END FUNCTION COULOMB_LOG_NRL_COUNTERSTREAM

      FUNCTION V_CRITICAL(te) RESULT(result)
         !--------------------------------------------------------------
         !   TE      Electron Temperature [eV]
         !   RESULT  Critical Velocity [m/s]
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(in) :: te
         DOUBLE PRECISION :: result
         result = fact_crit*SQRT(te)
         RETURN
      END FUNCTION V_CRITICAL

      FUNCTION V_CRITICAL_WEILAND(te,ion_coulomb_log,electron_coulomb_log) RESULT(result)
         !--------------------------------------------------------------
         !   TE                    Electron Temperature [eV]
         !   ION_COULOMB_LOG       Ion Coulomb Logarithm [arb]
         !   ELECTRON_COULOMB_LOG  Eelctron Coulomb Logarithm [arb]
         !   RESULT                Critical Velocity [m/s]
         !--------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(in) :: te,ion_coulomb_log,electron_coulomb_log
         DOUBLE PRECISION :: result
         !The coulomb ratio is from Weiland (2018) eq.11
         result = fact_crit_weiland*SQRT(te)*(ion_coulomb_log/electron_coulomb_log)**(1.0/3.0)
         RETURN
      END FUNCTION V_CRITICAL_WEILAND

      FUNCTION TAU_SPITZER(mass,ne,te,Z,coulomb_log) RESULT(result)
         !--------------------------------------------------------------
         !   MASS         Particle mass [kg]
         !   Z            Particle charge number [e]
         !   NE           Electron Density [m^-3]
         !   TE           Electron Temperature [eV]
         !   COULOMB_LOG  Coulomb Logarithm
         !   RESULT       Spitzer Collisional Timescale [s]
         !--------------------------------------------------------------
         DOUBLE PRECISION, INTENT(in) :: mass,ne,te,Z,coulomb_log
         DOUBLE PRECISION :: result
         result = 3.777183D41 * mass * SQRT( te * te * te ) / &
                    (ne * Z * Z * coulomb_log)
         RETURN
      END FUNCTION TAU_SPITZER
      
!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE collision_operators