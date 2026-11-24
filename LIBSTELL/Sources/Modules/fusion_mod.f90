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

        FUNCTION MAVRIN_COOLING_RATE(Z,Te) RESULT(result)
            !--------------------------------------------------------------
            !   Z            Atomic number [-]
            !   TE           Electron temperature [eV]
            !   RESULT       Total impurity cooling coefficient Lz [W m^3]
            !
            ! Mavrin (2018), doi:10.1080/10420150.2018.1462361.
            ! Coronal fits include bremsstrahlung, line and recombination
            ! radiation. Power density is ne*ni*Lz, where ni includes all
            ! charge states. Fits were derived for ne = 1e20 m^-3.
            ! As in pySTEL/libstell/radiation.py, clip Te to 0.1--100 keV
            ! and use the hotter interval at each internal boundary.
            ! Hydrogen isotopes (Z=1) must use BREMSSTRAHLUNG_POWER.
            !--------------------------------------------------------------
            USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
            USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
            INTEGER, INTENT(in) :: Z
            DOUBLE PRECISION, INTENT(in) :: Te
            DOUBLE PRECISION :: result, Tk, x, log_rate
            DOUBLE PRECISION :: A(0:4)
            INTEGER :: i

            IF (.NOT. IEEE_IS_FINITE(Te)) THEN
                ERROR STOP 'MAVRIN_COOLING_RATE: electron temperature must be finite'
            END IF
            Tk = MIN(MAX(Te,100.0D0),100000.0D0)*1.0D-3

            ! Coefficients A0,...,A4 in ascending powers of log10(Te [keV]).
            SELECT CASE (Z)
            CASE (2) ! Helium
                A = (/ -3.5551D1, 3.1469D-1, 1.0156D-1, -9.3730D-2, 2.5020D-2 /)
            CASE (3) ! Lithium
                A = (/ -3.5115D1, 1.9475D-1, 2.5082D-1, -1.6070D-1, 3.5190D-2 /)
            CASE (4) ! Beryllium
                A = (/ -3.4765D1, 3.7270D-2, 3.8363D-1, -2.1384D-1, 4.1690D-2 /)
            CASE (6) ! Carbon
                IF (Tk < 0.5D0) THEN
                    A = (/ -3.4738D1, -5.0085D0, -1.2788D1, -1.6637D1, -7.2904D0 /)
                ELSE
                    A = (/ -3.4174D1, -3.6687D-1, 6.8856D-1, -2.9191D-1, 4.4470D-2 /)
                END IF
            CASE (7) ! Nitrogen
                IF (Tk < 0.5D0) THEN
                    A = (/ -3.4065D1, -2.3614D0, -6.0605D0, -1.1570D1, -6.9621D0 /)
                ELSE IF (Tk < 2.0D0) THEN
                    A = (/ -3.3899D1, -5.9668D-1, 7.6272D-1, -1.7160D-1, 5.8770D-2 /)
                ELSE
                    A = (/ -3.3913D1, -5.2628D-1, 7.0047D-1, -2.2790D-1, 2.8350D-2 /)
                END IF
            CASE (8) ! Oxygen
                IF (Tk < 0.3D0) THEN
                    A = (/ -3.7257D1, -1.5635D1, -1.7141D1, -5.3765D0, 0.0D0 /)
                ELSE
                    A = (/ -3.3640D1, -7.6211D-1, 7.9655D-1, -2.0850D-1, 1.4360D-2 /)
                END IF
            CASE (10) ! Neon
                IF (Tk < 0.7D0) THEN
                    A = (/ -3.3132D1, 1.7309D0, 1.5230D1, 2.8939D1, 1.5648D1 /)
                ELSE IF (Tk < 5.0D0) THEN
                    A = (/ -3.3290D1, -8.7750D-1, 8.6842D-1, -3.9544D-1, 1.7244D-1 /)
                ELSE
                    A = (/ -3.3410D1, -4.5345D-1, 2.9731D-1, 4.3960D-2, -2.6930D-2 /)
                END IF
            CASE (18) ! Argon
                IF (Tk < 0.6D0) THEN
                    A = (/ -3.2155D1, 6.5221D0, 3.0769D1, 3.9161D1, 1.5353D1 /)
                ELSE IF (Tk < 3.0D0) THEN
                    A = (/ -3.2530D1, 5.4490D-1, 1.5389D0, -7.6887D0, 4.9806D0 /)
                ELSE
                    A = (/ -3.1853D1, -1.6674D0, 6.1339D-1, 1.7480D-1, -8.2260D-2 /)
                END IF
            CASE (36) ! Krypton
                IF (Tk < 0.447D0) THEN
                    A = (/ -3.4512D1, -2.1484D1, -4.4723D1, -4.0133D1, -1.3564D1 /)
                ELSE IF (Tk < 2.364D0) THEN
                    A = (/ -3.1399D1, -5.0091D-1, 1.9148D0, -2.5865D0, -5.2704D0 /)
                ELSE
                    A = (/ -2.9954D1, -6.3683D0, 6.6831D0, -2.9674D0, 4.8356D-1 /)
                END IF
            CASE (54) ! Xenon
                IF (Tk < 0.5D0) THEN
                    A = (/ -2.9303D1, 1.4351D1, 4.7081D1, 5.9580D1, 2.5615D1 /)
                ELSE IF (Tk < 2.5D0) THEN
                    A = (/ -3.1113D1, 5.9339D-1, 1.2808D0, -1.1628D1, 1.0748D1 /)
                ELSE IF (Tk < 10.0D0) THEN
                    A = (/ -2.5813D1, -2.7526D1, 4.8614D1, -3.6885D1, 1.0069D1 /)
                ELSE
                    A = (/ -2.2138D1, -2.2592D1, 1.9619D1, -7.5181D0, 1.0858D0 /)
                END IF
            CASE (74) ! Tungsten
                IF (Tk < 1.5D0) THEN
                    A = (/ -3.0374D1, 3.8304D-1, -9.5126D-1, -1.0311D0, -1.0103D-1 /)
                ELSE IF (Tk < 4.0D0) THEN
                    A = (/ -3.0238D1, -2.9208D0, 2.2824D1, -6.3303D1, 5.1849D1 /)
                ELSE
                    A = (/ -3.2153D1, 5.2499D0, -6.2740D0, 2.6627D0, -3.6759D-1 /)
                END IF
            CASE DEFAULT
                WRITE(ERROR_UNIT,*) 'No Mavrin radiation fit for atomic number ',Z
                ERROR STOP 'MAVRIN_COOLING_RATE: unsupported atomic number'
            END SELECT

            x = LOG10(Tk)
            log_rate = A(4)
            DO i=3,0,-1
                log_rate = log_rate*x + A(i)
            END DO
            result = 10.0D0**log_rate
            RETURN
        END FUNCTION MAVRIN_COOLING_RATE

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
