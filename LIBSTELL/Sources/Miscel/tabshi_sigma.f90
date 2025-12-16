!-----------------------------------------------------------------------
!     Module:        tabshi_sigma
!     Authors:       L van Ham (lucas.van.ham@ipp.mpg.de)
!     Date:          01/12/2025
!     Description:   This module includes analytical functions used to
!                    calculate the cross-sections for different reactions
!                    for different energetic hydrogenic species upon
!                    colliding with a "cold" H2 gas. Used in modelling
!                    NBI neutralizer physics. These functions were
!                    lifted from the paper
!
!                       T. Tabata, T. Shirai,       
!                  ANALYTIC CROSS SECTIONS FOR COLLISIONS OF H+, H2+,
!                  H3+, H, H2, AND H- WITH HYDROGEN MOLECULES
!                  At. Data Nucl. Data Tables 76 (2000), Issue 1,p1-25
!         
!-----------------------------------------------------------------------
MODULE tabshi_sigma
!-----------------------------------------------------------------------
!   MODULE PARAMETERS
!-----------------------------------------------------------------------
CONTAINS
    FUNCTION get_functional_1(x,c1,c2)                     result(f1)
        !-------------------------------------------------------------------
        !     Definition (i) in Tabata (2000) 
        !       
        !       Input parameters
        !           x, c1, c2
        !       Output parameters
        !           f1
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: f1
        DOUBLE PRECISION, INTENT(in) :: x, c1, c2
        DOUBLE PRECISION :: sigma0 = 1.0E-20 ! In m^-2
        DOUBLE PRECISION :: ERyd =  1.361E-2 ! keV Rydberg constant

        f1 = sigma0*c1*(x/ERyd)**c2
        RETURN

    END FUNCTION get_functional_1

    FUNCTION get_functional_2(x,c1,c2,c3,c4)               result(f2)
        !-------------------------------------------------------------------
        !     Definition (ii) in Tabata (2000) 
        !       
        !       Input parameters
        !           x, c1, c2, c3, c4
        !       Output parameters
        !           f2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: f2
        DOUBLE PRECISION, INTENT(in) :: x, c1, c2, c3, c4
        DOUBLE PRECISION :: f1, denom

        f1 = get_functional_1(x,c1,c2)
        denom = 1.0+(x/c3)**(c2+c4)
        f2 = f1/denom 
        RETURN

    END FUNCTION get_functional_2

    FUNCTION get_functional_3(x,c1,c2,c3,c4,c5,c6)         result(f3)
        !-------------------------------------------------------------------
        !     Definition (iii) in Tabata (2000) 
        !       
        !       Input parameters
        !           x, c1, c2, c3, c4, c5, c6
        !       Output parameters
        !           f3
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: f3
        DOUBLE PRECISION, INTENT(in) :: x, c1, c2, c3, c4, c5, c6
        DOUBLE PRECISION :: f1, denom

        f1 = get_functional_1(x,c1,c2)
        denom = 1.0+(x/c3)**(c2+c4)+(x/c5)**(c2+c6)
        f3 = f1/denom 
        RETURN

    END FUNCTION get_functional_3

    FUNCTION get_functional_4(x,c1,c2,c3,c4,c5,c6,c7,c8)   result(f4)
        !-------------------------------------------------------------------
        !     Definition (iv) in Tabata (2000) 
        !       
        !       Input parameters
        !           x, c1, c2, c3, c4, c5, c6, c7, c8
        !       Output parameters
        !           f4
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: f4
        DOUBLE PRECISION, INTENT(in) :: x, c1, c2, c3, c4, c5, c6, c7, c8
        DOUBLE PRECISION :: f1, num, denom

        f1 = get_functional_1(x,c1,c2)
        num = 1+(x/c3)**(c4-c2)
        denom = 1.0+(x/c5)**(c4+c6)+(x/c7)**(c4+c8)
        f4 = f1*num/denom 
        RETURN

    END FUNCTION get_functional_4

    FUNCTION get_sigma_eq2(E1,a1,a2,a3,a4,a5,a6,a7,a8,a9)   result(sigma2)
        !-------------------------------------------------------------------
        !     Equation (2) in Tabata (2000) 
        !       
        !       Input parameters
        !           E1, a2, a3, a4, a5, a6
        !       Output parameters
        !           sigma2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma2
        DOUBLE PRECISION, INTENT(in) :: E1, a1, a2, a3, a4, a5, a6
        DOUBLE PRECISION :: f2_a, f2_b

        f2_a = get_functional_2(E1,   a1,a2,a3,a4)
        f2_b = get_functional_2(E1/a6,a1,a2,a3,a4)
        sigma2 = f2_a+a5*f2_b
        RETURN

    END FUNCTION get_sigma_eq2    

    FUNCTION get_sigma_eq6(E1,a1,a2,a3,a4,a5,a6)   result(sigma6)
        !-------------------------------------------------------------------
        !     Equation (6) in Tabata (2000) 
        !       
        !       Input parameters
        !           E1, a1, a2, a3, a4, a5, a6
        !       Output parameters
        !           sigma6
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma6
        DOUBLE PRECISION, INTENT(in) :: E1, a1, a2, a3, a4, a5, a6
        DOUBLE PRECISION :: f3

        f3 = get_functional_3(E1,a1,a2,a3,a4,a5,a6)
        sigma6 = f3 
        RETURN

    END FUNCTION get_sigma_eq6

    FUNCTION get_sigma_eq8(E1,a1,a2,a3,a4,a5,a6,a7,a8,a9)   result(sigma8)
        !-------------------------------------------------------------------
        !     Equation (8) in Tabata (2000) 
        !       
        !       Input parameters
        !           E1, a1, a2, a3, a4, a5, a6, a7, a8, a9
        !       Output parameters
        !           sigma8
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma8
        DOUBLE PRECISION, INTENT(in) :: E1, a1, a2, a3, a4, a5, a6, a7, a8, a9
        DOUBLE PRECISION :: f2, f3

        f2 = get_functional_2(E1,a1,a2,a3,a4)
        f3 = get_functional_3(E1,a5,a2,a6,a7,a8,a9)
        sigma8 = f2+f3
        RETURN

    END FUNCTION get_sigma_eq8        

    FUNCTION get_sigma_eq11(E1,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)   result(sigma11)
        !-------------------------------------------------------------------
        !     Equation (11) in Tabata (2000) 
        !       
        !       Input parameters
        !           E1, a2, a3, a4, a5, a6, a7, a8, a9, a10
        !       Output parameters
        !           sigma11
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma11
        DOUBLE PRECISION, INTENT(in) :: E1, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10
        DOUBLE PRECISION :: f2, f3

        f3 = get_functional_3(E1,a1,a2,a3,a4,a5,a6)
        f2 = get_functional_3(E1,a7,a8,a9,a10)
        sigma11 = f3 + f2
        RETURN

    END FUNCTION get_sigma_eq11

    FUNCTION get_sigma_eq13(E1,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)   result(sigma13)
        !-------------------------------------------------------------------
        !     Equation (13) in Tabata (2000) 
        !       
        !       Input parameters
        !           E1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12
        !       Output parameters
        !           sigma13
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma13
        DOUBLE PRECISION, INTENT(in) :: E1, a1, a2, a3, a4,  a5,  a6
        DOUBLE PRECISION, INTENT(in) ::     a7, a8, a9, a10, a11, a12
        DOUBLE PRECISION :: f3_a, f3_b

        f3_a = get_functional_3(E1,a1,a2,a3,a4,a5,a6)
        f3_b = get_functional_3(E1,a7,a8,a9,a10,a11,a12)
        sigma13 = f3_a + f3_b
        RETURN

    END FUNCTION get_sigma_eq13

    FUNCTION get_sigma_neut_Hplus(E) result(sigma)
        !-------------------------------------------------------------------
        !     Reaction (6) in Tabata (2000) [H+ + H2 -> fast H]
        !       
        !       Input parameters
        !           E       energy in keV
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E 
        DOUBLE PRECISION :: Eth, E1
        DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6, a7, a8, a9

        Eth = 2.5E-3 ! Threshold energy in keV
        E1   = E - Eth
        a1 = 2.12E2;  a2 = 1.721; a3 = 6.7E-4; a4 = 3.239E-1; a5 = 4.34E-3; a6 = 1.296; 
        a7 = 1.42E-1; a8 = 9.34;  a9 = 2.997;
        sigma = get_sigma_eq8(E1,a1,a2,a3,a4,a5,a6,a8,a9)*1e-4
        RETURN

    END FUNCTION get_sigma_neut_Hplus

    FUNCTION get_sigma_ionp_Hneut(E) result(sigma)
        !-------------------------------------------------------------------
        !     Reaction (31) in Tabata (2000) [H + H2 -> fast H+]
        !       
        !       Input parameters
        !           E       energy in keV
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E 
        DOUBLE PRECISION :: Eth, E1
        DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6

        Eth = 2.0E-2 ! Threshold energy in keV
        E1 = E - Eth
        a1 = 2.53E-4; a2 = 1.728; a3 = 2.164; a4 = 7.74e-1; a5 = 1.639; a6 = 1.43E1;
        sigma = get_sigma_eq2(E1,a1,a2,a3,a4,a5,a6)*1e-4
        RETURN

    END FUNCTION get_sigma_ionp_Hneut

    FUNCTION get_sigma_ionn_Hneut(E) result(sigma)
        !-------------------------------------------------------------------
        !     Reaction (29) in Tabata (2000) [H + H2 -> fast H-]
        !       
        !       Input parameters
        !           E       energy in keV
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E 
        DOUBLE PRECISION :: Eth, E1
        DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12

        Eth = 2.1E-2 ! Threshold energy in keV
        E1 = E - Eth
        a1 = 9.73E-3; a2 = 2.38;  a3 = 1.39E-2; a4 = -5.51E-1; a5  = 7.7E-2; a6  = 2.12;
        a7 = 1.97E-6; a8 = 2.051; a9 = 5.5;     a10 = 6.62E-1; a11 = 2.02E1; a12 = 3.62;
        sigma = get_sigma_eq13(E1,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)*1e-4
        RETURN

    END FUNCTION get_sigma_ionn_Hneut   

    FUNCTION get_sigma_neut_Hmin(E) result(sigma)
        !-------------------------------------------------------------------
        !     Reaction (47) in Tabata (2000) [H- + H2 -> fast H]
        !       
        !       Input parameters
        !           E       energy in keV
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E 
        DOUBLE PRECISION :: Eth, E1
        DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12

        Eth = 2.25E-3 ! Threshold energy in keV
        E1 = E - Eth
        a1 = 4.19E-2; a2 = 1.89;  a3 = 1.78E-1; a4 = -2.3E-1; a5 = 1.04; a6 = 8.7E-1;
        a7 = 1.65E1;  a8 = 1.088; a9 = 5.33E-3; a10 = 1.66E-1;
        sigma = get_sigma_eq11(E1,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*1e-4; 
        RETURN

    END FUNCTION get_sigma_neut_Hmin
    
    FUNCTION get_sigma_ionp_Hmin(E) result(sigma)
        !-------------------------------------------------------------------
        !     Reaction (49) in Tabata (2000) [H- + H2 -> fast H+]
        !       
        !       Input parameters
        !           E       energy in keV
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E 
        DOUBLE PRECISION :: Eth, E1
        DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6

        Eth = 0 ! Threshold energy in keV
        E1 = E - Eth
        a1 = 1.75E-8; a2 = 3.88; a3 = 9.06E-1; a4 = -2.74E-1; a5 = 3.19; a6 = 1.19;
        sigma = get_sigma_eq6(E1,a1,a2,a3,a4,a5,a6)*1e-4; 
        RETURN

    END FUNCTION get_sigma_ionp_Hmin 

END MODULE 