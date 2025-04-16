!-----------------------------------------------------------------------
!     Module:        fusion_mod
!     Authors:       A. J. Coelho
!     Date:          04/16/2025
!     Description:   This module contains routines for working with
!                    quantities realated to nuclear fusion   
!-----------------------------------------------------------------------
      MODULE fusion_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, PRIVATE, PARAMETER :: &
                              electron_charge = 1.602176634D-19 !e_c
!-----------------------------------------------------------------------
      CONTAINS

        FUNCTION BREMSSTRAHLUNG_POWER(Z,ni,ne,Te) RESULT(result)
            !--------------------------------------------------------------
            !   Z            Ion charge number [-]
            !   NI           Ion Density [m^-3]
            !   NE           Electron Density [m^-3]
            !   TE           Electron Temperature [eV]
            !   RESULT       Bremstrahlung power density [W/m^3]
            !
            ! This func computes the Bremsstrahlung irradiated power density
            ! according to Freidberg's Plasma Physics and Fusion Energy, 
            ! page 56, Eq. (3.42)
            !--------------------------------------------------------------
            DOUBLE PRECISION, INTENT(in) :: ni,ne,Te
            INTEGER, INTENT(in) :: Z
            DOUBLE PRECISION :: result, CB, n20, Tk, Zeff
            CB = 5.35D3
            n20 = ne / 1.0D20
            Tk = Te / 1.0D3 
            Zeff = Z*Z * ni/ne
            result = CB * Zeff * n20**2 * SQRT(Tk)  !W/m^3
            RETURN
        END FUNCTION BREMSSTRAHLUNG_POWER

        FUNCTION DT_CROSS_SECTION(Ti) RESULT(result)
            !--------------------------------------------------------------
            !   TI           Ion Temperature [eV]
            !   RESULT       Cross section <sigma.v> [m^3/s]
            !
            ! This function calculates the fusion cross section <sigma.v> 
            ! given the ion temperature using the Bosch Hale model.
            ! H.-S. Bosch and G. M. Hale 1992 Nucl. Fusion 32 611
            ! https://doi.org/10.1088/0029-5515/32/4/I07
            !--------------------------------------------------------------
            DOUBLE PRECISION, INTENT(in) :: Ti
            INTEGER, PARAMETER :: mrc2 = 1124656
            DOUBLE PRECISION, PARAMETER :: BG = 34.3827
            DOUBLE PRECISION, DIMENSION(7), PARAMETER :: &
                CARR = (/ 1.17302D-09,  1.51361D-02,  7.51886D-02, &
                          4.60643D-03,  1.35000D-02, -1.06750D-04, &
                          1.36600D-05/)
            DOUBLE PRECISION :: Tk, zeta, theta, eta, result

            Tk = Ti*1.0D-3 ! to keV
            zeta =  1.0D0 - ((((CARR(6)*Tk)+CARR(4))*Tk+CARR(2))*Tk)/ &
                        ((((CARR(7)*Tk)+CARR(5))*Tk+CARR(3))*Tk+1.0D0)
            theta = Tk/zeta
            eta   = (BG*BG/(4*theta))**(1.0D0/3.0D0)
    
            result = 1.0D-6*CARR(1)*theta*SQRT(eta/(mrc2*Tk*Tk*Tk))*EXP(-3*eta)
            RETURN
        END FUNCTION DT_CROSS_SECTION

        FUNCTION DD_CROSS_SECTION(Ti) RESULT(result)
            !--------------------------------------------------------------
            !   TI           Ion Temperature [eV]
            !   RESULT       Cross section <sigma.v> [m^3/s]
            !
            ! This function calculates the fusion cross section <sigma.v> 
            ! given the ion temperature using the Bosch Hale model.
            ! H.-S. Bosch and G. M. Hale 1992 Nucl. Fusion 32 611
            ! https://doi.org/10.1088/0029-5515/32/4/I07
            !--------------------------------------------------------------
            DOUBLE PRECISION, INTENT(in) :: Ti
            INTEGER, PARAMETER :: mrc2 = 937814
            DOUBLE PRECISION, PARAMETER :: BG = 31.3970
            DOUBLE PRECISION, DIMENSION(7), PARAMETER :: &
                CARR = (/ 5.65718D-12,  3.41267D-03,  1.99167D-03, &
                          0.00000D+00,  1.05060D-05,  0.00000D+00, &
                          0.00000D-00/)
            DOUBLE PRECISION :: Tk, zeta, theta, eta, result

            Tk = Ti*1.0D-3 ! to keV
            zeta =  1.0D0 - ((((CARR(6)*Tk)+CARR(4))*Tk+CARR(2))*Tk)/ &
                        ((((CARR(7)*Tk)+CARR(5))*Tk+CARR(3))*Tk+1.0D0)
            theta = Tk/zeta
            eta   = (BG*BG/(4*theta))**(1.0D0/3.0D0)
    
            result = 1.0D-6*CARR(1)*theta*SQRT(eta/(mrc2*Tk*Tk*Tk))*EXP(-3*eta)
            RETURN
        END FUNCTION DD_CROSS_SECTION

        FUNCTION DDHe3_CROSS_SECTION(Ti) RESULT(result)
            !--------------------------------------------------------------
            !   TI           Ion Temperature [eV]
            !   RESULT       Cross section <sigma.v> [m^3/s]
            !
            ! This function calculates the fusion cross section <sigma.v> 
            ! given the ion temperature using the Bosch Hale model.
            ! H.-S. Bosch and G. M. Hale 1992 Nucl. Fusion 32 611
            ! https://doi.org/10.1088/0029-5515/32/4/I07
            !--------------------------------------------------------------
            DOUBLE PRECISION, INTENT(in) :: Ti
            INTEGER, PARAMETER :: mrc2 = 937814
            DOUBLE PRECISION, PARAMETER :: BG = 31.3970
            DOUBLE PRECISION, DIMENSION(7), PARAMETER :: &
                CARR = (/ 5.43360D-12,  5.85778D-03,  7.68222D-03, &
                          0.00000D+00, -2.96400D-06,  0.00000D+00, &
                          0.00000D+00/)
            DOUBLE PRECISION :: Tk, zeta, theta, eta, result

            Tk = Ti*1.0D-3 ! to keV
            zeta =  1.0D0 - ((((CARR(6)*Tk)+CARR(4))*Tk+CARR(2))*Tk)/ &
                        ((((CARR(7)*Tk)+CARR(5))*Tk+CARR(3))*Tk+1.0D0)
            theta = Tk/zeta
            eta   = (BG*BG/(4*theta))**(1.0D0/3.0D0)
    
            result = 1.0D-6*CARR(1)*theta*SQRT(eta/(mrc2*Tk*Tk*Tk))*EXP(-3*eta)
            RETURN
        END FUNCTION DDHe3_CROSS_SECTION

        FUNCTION DHe3_CROSS_SECTION(Ti) RESULT(result)
            !--------------------------------------------------------------
            !   TI           Ion Temperature [eV]
            !   RESULT       Cross section <sigma.v> [m^3/s]
            !
            ! This function calculates the fusion cross section <sigma.v> 
            ! given the ion temperature using the Bosch Hale model.
            ! H.-S. Bosch and G. M. Hale 1992 Nucl. Fusion 32 611
            ! https://doi.org/10.1088/0029-5515/32/4/I07
            !--------------------------------------------------------------
            DOUBLE PRECISION, INTENT(in) :: Ti
            INTEGER, PARAMETER :: mrc2 = 1124572
            DOUBLE PRECISION, PARAMETER :: BG = 68.7508
            DOUBLE PRECISION, DIMENSION(7), PARAMETER :: &
                CARR = (/ 5.51036D-10,  6.41918D-03, -2.02896D-03, &
                         -1.91080D-05,  1.35776D-04,  0.00000D+00, &
                          0.00000D+00/)
            DOUBLE PRECISION :: Tk, zeta, theta, eta, result

            Tk = Ti*1.0D-3 ! to keV
            zeta =  1.0D0 - ((((CARR(6)*Tk)+CARR(4))*Tk+CARR(2))*Tk)/ &
                        ((((CARR(7)*Tk)+CARR(5))*Tk+CARR(3))*Tk+1.0D0)
            theta = Tk/zeta
            eta   = (BG*BG/(4*theta))**(1.0D0/3.0D0)
    
            result = 1.0D-6*CARR(1)*theta*SQRT(eta/(mrc2*Tk*Tk*Tk))*EXP(-3*eta)
            RETURN
        END FUNCTION DHe3_CROSS_SECTION

        FUNCTION ALPHA_POWER(nD,nT,TD,TT) RESULT(result)
            !--------------------------------------------------------------
            !   ND           Deuterium Density [m^-3]
            !   NT           Tritium Density [m^-3]
            !   TD           Deuterium Temperature [eV]
            !   TT           Tritium Temperature [eV]
            !   result       Alpha power density [W/m^3]
            !
            ! This func computes the alpha power density = 
            ! nD*nT*<sigma.v>*E_alpha
            !--------------------------------------------------------------
            DOUBLE PRECISION, INTENT(in) :: nD,nT,TD,TT
            DOUBLE PRECISION, PARAMETER :: E_alpha = 3.52D6*electron_charge
            DOUBLE PRECISION :: Ti, sigmav, result
            Ti = 0.5D+00 * (TD+TT)
            sigmav = DT_CROSS_SECTION(Ti) ! m^3/s
            result = nD * nT * sigmav * E_alpha ! W/m^3
            RETURN
        END FUNCTION ALPHA_POWER

    END MODULE fusion_mod

        