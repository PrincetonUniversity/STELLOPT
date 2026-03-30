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
    PUBLIC :: get_sigma_Hp_H, get_sigma_Hn_H,  &
              get_sigma_H_Hp, get_sigma_Hm_Hp, &
              get_sigma_H_Hm
    PUBLIC :: get_sigma_H2p_H2, get_sigma_H2_H2p, &
              get_sigma_H2p_HHp, get_sigma_H2_HHp
    PUBLIC :: get_sigma_H3p_H2H, get_sigma_H3p_3H, &
              get_sigma_H3p_H2pH, get_sigma_H3p_H2Hp,&
              get_sigma_H3p_Hp2H, get_sigma_H3p_H2pHp, &
              get_sigma_H3p_2HpH, get_sigma_H3p_3Hp
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

    FUNCTION get_sigma_eq2(E1,a1,a2,a3,a4,a5,a6)   result(sigma2)
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
    
    FUNCTION get_sigma_eq3(E1,a1,a2,a3,a4,a5,a6,a7,a8)   result(sigma3)
        !-------------------------------------------------------------------
        !     Equation (3) in Tabata (2000) 
        !       
        !       Input parameters
        !           E1, a2, a3, a4, a5, a6, a7, a8
        !       Output parameters
        !           sigma3
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma3
        DOUBLE PRECISION, INTENT(in) :: E1, a1, a2, a3, a4, a5, a6, a7, a8
        DOUBLE PRECISION :: f2_a, f2_b

        f2_a = get_functional_2(E1,a1,a2,a3,a4)
        f2_b = get_functional_2(E1,a5,a6,a7,a8)
        sigma3 = f2_a+f2_b

        RETURN

    END FUNCTION get_sigma_eq3    

    FUNCTION get_sigma_eq4(E1,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)   result(sigma4)
        !-------------------------------------------------------------------
        !     Equation (4) in Tabata (2000) 
        !       
        !       Input parameters
        !           E1, a2, a3, a4, a5, a6, a7, a8, a9, a10
        !       Output parameters
        !           sigma4
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma4
        DOUBLE PRECISION, INTENT(in) :: E1, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10
        DOUBLE PRECISION :: f2_a, f2_b, f2_c

        f2_a = get_functional_2(E1,    a1,a2,a3,a4)
        f2_b = get_functional_2(E1,    a5,a6,a7,a8)
        f2_c = get_functional_2(E1/a10,a5,a6,a7,a8)
        sigma4 = f2_a+f2_b+a9*f2_c
        RETURN

    END FUNCTION get_sigma_eq4   

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
    
    FUNCTION get_sigma_eq10(E1,a1,a2,a3,a4,a5,a6,a7,a8)   result(sigma10)
        !-------------------------------------------------------------------
        !     Equation (10) in Tabata (2000) 
        !       
        !       Input parameters
        !           E1, a1, a2, a3, a4, a5, a6, a7, a8
        !       Output parameters
        !           sigma10
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma10
        DOUBLE PRECISION, INTENT(in) :: E1, a1, a2, a3, a4, a5, a6, a7, a8
        DOUBLE PRECISION :: f3_a, f3_b

        f3_a = get_functional_3(E1,a1,a2,a3,a4,a5,a6)
        f3_b = get_functional_3(E1/a8,a1,a2,a3,a4,a5,a6)
        sigma10 = f3_a+a7*f3_b
        RETURN

    END FUNCTION get_sigma_eq10

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
        f2 = get_functional_2(E1,a7,a8,a9,a10)
        sigma11 = f3 + f2
        RETURN

    END FUNCTION get_sigma_eq11

    FUNCTION get_sigma_eq12(E1,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)   result(sigma12)
        !-------------------------------------------------------------------
        !     Equation (12) in Tabata (2000) 
        !       
        !       Input parameters
        !           E1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12
        !       Output parameters
        !           sigma12
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma12
        DOUBLE PRECISION, INTENT(in) :: E1, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12
        DOUBLE PRECISION :: f2_a, f2_b, f3

        f3 = get_functional_3(E1,a1,a2,a3,a4,a5,a6)
        f2_a = get_functional_2(E1,a7,a8,a9,a10)
        f2_b = get_functional_2(E1/a12,a7,a8,a9,a10)
        sigma12 = f3 + f2_a + a11*f2_b
        RETURN

    END FUNCTION get_sigma_eq12

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

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!       SINGLE-PROTON REACTIONS (H+, H, H-)
!           Reaction      Reactant      Fast product(s)           
!               6            H+           H
!               29           H            H-
!               31           H            H+
!               47           H-           H
!               49           H-           H+
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

    FUNCTION get_sigma_Hp_H(E) result(sigma)
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
        sigma = get_sigma_eq8(E1,a1,a2,a3,a4,a5,a6,a7,a8,a9)
        RETURN

    END FUNCTION get_sigma_Hp_H

    FUNCTION get_sigma_H_Hp(E) result(sigma)
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
        sigma = get_sigma_eq2(E1,a1,a2,a3,a4,a5,a6)
        RETURN

    END FUNCTION get_sigma_H_Hp

    FUNCTION get_sigma_H_Hm(E) result(sigma)
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
        sigma = get_sigma_eq13(E1,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)
        RETURN

    END FUNCTION get_sigma_H_Hm   

    FUNCTION get_sigma_Hn_H(E) result(sigma)
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
        sigma = get_sigma_eq11(E1,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
        RETURN

    END FUNCTION get_sigma_Hn_H
    
    FUNCTION get_sigma_Hm_Hp(E) result(sigma)
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
        sigma = get_sigma_eq6(E1,a1,a2,a3,a4,a5,a6)
        RETURN

    END FUNCTION get_sigma_Hm_Hp 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!       DOUBLE-PROTON REACTIONS (H2+, H2)
!           Reaction      Reactant      Fast product(s)           
!               14           H2           H+, H 
!               38           H2           H2+
!               43           H2+          H+, H
!               50           H2+          H2 
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

    FUNCTION get_sigma_H2p_H2(E) result(sigma)
        !-------------------------------------------------------------------
        !     Reaction (50) in Tabata (2000) [H2+ + H2 -> fast H2]
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

        Eth = 0.0E-3 ! Threshold energy in keV
        E1 = E - Eth
        a1 = 2.29E+2; a2 = 2.78;    a3 = 4.75E-3; a4  = 1.248E-1; a5 = 2.14E-1;  a6 = 2.33;
        a7 = 7.96;    a8 = 6.82E-1; a9 = 6.59E-3; a10 = 4.51;    a11 = 1.67E-1; a12 = 1.164E+4;
        sigma = get_sigma_eq12(E1,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)
        RETURN

    END FUNCTION get_sigma_H2p_H2

    FUNCTION get_sigma_H2_H2p(E) result(sigma)
        !-------------------------------------------------------------------
        !     Reaction (38) in Tabata (2000) [H2 + H2 -> fast H2+]
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
        DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6, a7, a8

        Eth = 3.2E-2 ! Threshold energy in keV
        E1 = E - Eth
        a1 = 1.879E-3; a2 = 2.497; a3 = 6.62E-2; a4 = -4.67E-1; a5 = 3.58E-1; a6 = 5.0E-1;
        a7 = 7.67;     a8 = 2.01E2;
        sigma = get_sigma_eq10(E1,a1,a2,a3,a4,a5,a6,a7,a8)
        RETURN

    END FUNCTION get_sigma_H2_H2p

    FUNCTION get_sigma_H2p_HHp(E) result(sigma)
        !-------------------------------------------------------------------
        !     Reaction (14) in Tabata (2000) [H2+ + H2 -> fast H+ + fast H]
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
        DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10

        Eth = 5.0E-3 ! Threshold energy in keV
        E1 = E - Eth
        a1 = 6.34E+1; a2 = 1.78; a3 = 1.38E-3; a4 = 4.06E-1; a5 = 1.63E-1; a6 = 3.27E-1;
        a7 = 1.554E+1; a8 = 3.903; a9 = 1.735; a10 = 1.02E+1
        sigma = get_sigma_eq4(E1,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
        RETURN

    END FUNCTION get_sigma_H2p_HHp

    FUNCTION get_sigma_H2_HHp(E) result(sigma)
        !-------------------------------------------------------------------
        !     Reaction (43) in Tabata (2000) [H2 + H2 -> fast H+ + fast H]
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
        a1 = 1.307E-5; a2 = 1.586; a3 = 1.066E+1; a4 = 2.03; a5 = 2.73; a6 = 4.71
        sigma = get_sigma_eq2(E1,a1,a2,a3,a4,a5,a6)
        RETURN

    END FUNCTION get_sigma_H2_HHp

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!       TRIPLE-PROTON REACTIONS (H3+)
!           Reaction      Reactant      Product(s)     Note
!               18           H3+         at least 1 H+
!               19           H3+         at least 1 H2+
!               20           H3+         at least 1 H
!               21           H3+         at least 1 H2
!         Most of these reactions do not translate nicely to a simple MC
!         cross-section (e.g. reaction 18 can imply (H+,H2) dissociation or
!         (H+,H,H) dissociation) so have to do some linear algebra
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

    FUNCTION get_sigma_18(E) result(sigma)
        !-------------------------------------------------------------------
        !     Reaction (18) in Tabata (2000)
        !       Input parameters
        !           E       energy in keV
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E 
        DOUBLE PRECISION :: Eth, E1
        DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6, a7, a8

        Eth = 1.1E-2 ! Threshold energy in keV
        E1 = E - Eth
        a1 = 6.67E-1; a2 = 1.35; a3 = 4.42E-2; a4 = 7.1E-1; a5 = 6.7E-5; a6 = 1.54;
        a7 = 1.1E1; a8 = -1.0E-1;
        sigma = get_sigma_eq3(E1,a1,a2,a3,a4,a5,a6,a7,a8)
        RETURN

    END FUNCTION get_sigma_18

    FUNCTION get_sigma_19(E) result(sigma)
        !-------------------------------------------------------------------
        !     Reaction (19) in Tabata (2000)
        !       Input parameters
        !           E       energy in keV
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E 
        DOUBLE PRECISION :: Eth, E1
        DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6, a7, a8

        Eth = 1.55E-2 ! Threshold energy in keV
        E1 = E - Eth
        a1 = 5.03E-1; a2 = 1.0; a3 = 2.5E-2; a4 = 2.0; a5 = 1.17E-1; a6 = 3.18E-1;
        a7 = 9.4E+1; a8 = 1.35
        sigma = get_sigma_eq3(E1,a1,a2,a3,a4,a5,a6,a7,a8)

        RETURN
    END FUNCTION get_sigma_19

    FUNCTION get_sigma_20(E) result(sigma)
        !-------------------------------------------------------------------
        !     Reaction (20) in Tabata (2000)
        !       Input parameters
        !           E       energy in keV
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E 
        DOUBLE PRECISION :: Eth, E1
        DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6, a7, a8

        Eth = 1.55E-2 ! Threshold energy in keV
        E1 = E - Eth
        a1 = 5.89E-1; a2 = 1.0; a3 = 2.5E-2; a4 = 1.5; a5 = 4.05E-2; a6 = 7.59E-1;
        a7 = 4.64E+1; a8 = 1.1
        
        sigma = get_sigma_eq3(E1,a1,a2,a3,a4,a5,a6,a7,a8)

        RETURN
    END FUNCTION get_sigma_20

    FUNCTION get_sigma_21(E) result(sigma)
        !-------------------------------------------------------------------
        !     Reaction (21) in Tabata (2000) 
        !       Input parameters
        !           E       energy in keV
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E 
        DOUBLE PRECISION :: Eth, E1
        DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6, a7, a8

        Eth = 0.0 ! Threshold energy in keV
        E1 = E - Eth
        a1 = 3.78E+1; a2 = 1.0; a3 = 2.0E-3; a4 = 2.5E-1; a5 = 4.14E-2; a6 = 6.25E-1;
        a7 = 4.89E+1; a8 = 1.69
        sigma = get_sigma_eq3(E1,a1,a2,a3,a4,a5,a6,a7,a8)

        RETURN
    END FUNCTION get_sigma_21
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!       ACTUAL TRIPLE-PROTON REACTIONS (H3+)
!           1. H3+ + H2 -> H2  + H        + (H2+) 
!           2. H3+ + H2 -> H   + H  + H   + (H2+)
!           3. H3+ + H2 -> H2+ + H        + (H2)
!           4. H3+ + H2 -> H2  + H+       + (H2) 
!           5. H3+ + H2 -> H+  + H  + H   + (H2) 
!           6. H3+ + H2 -> H2+ + H+       + (H2 +  e-)
!           7. H3+ + H2 -> H+  + H+ + H   + (H2 +  e-)
!           8. H3+ + H2 -> H+  + H+ + H+  + (H2 + 2e-)
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

!-------------------------------------------------------------------
!    GET_SIGMAS_H3P : Calculates cross-sections of H3+ reactions.
!-------------------------------------------------------------------
    FUNCTION get_sigmas_H3p(E) result(sigmas_H3p)

        IMPLICIT NONE
        DOUBLE PRECISION :: E ! Energy in keV
        DOUBLE PRECISION :: sigmas_H3p(8)
        DOUBLE PRECISION :: s(4), mat(4,4)
        DOUBLE PRECISION, PARAMETER :: BR_cd = 0.4d0,  &     ! "Free" branching ratios
                                    BR_fd = 0.01d0, &
                                    BR_ge = 0.1d0,  &
                                    BR_he = 0.05d0
        DOUBLE PRECISION, PARAMETER :: BR_dd = 1.0d0-BR_cd-BR_fd, & ! Fixed branching ratios
                                    BR_ee = 1.0d0-BR_ge-BR_he
        INTEGER :: IPIV(4), INFO

        ! Get cross-sections
        s(1)=get_sigma_18(E)
        s(2)=get_sigma_19(E)
        s(3)=get_sigma_20(E)
        s(4)=get_sigma_21(E)

        ! build matrix
        mat(1,:) = (/ 0d0, 0d0, 1*BR_dd+1*BR_fd, 1*BR_ee+2*BR_ge+3*BR_he/)
        mat(2,:) = (/ 0d0, 0d0, 1*BR_cd+1*BR_fd, 0d0/)
        mat(3,:) = (/ 1d0, 3d0, 1*BR_cd, 2*BR_ee+1*BR_ge/)
        mat(4,:) = (/ 1d0, 0d0, 1*BR_dd, 0d0/)
        
        ! solve
        CALL DGESV(4, 1, mat, 4, IPIV, s, 4, INFO)
        IF (INFO /= 0) THEN
            PRINT *, 'DGESV failed, info = ', info
        END IF

        sigmas_H3p(1) = s(1)        
        sigmas_H3p(2) = s(2)
        sigmas_H3p(3) = BR_cd*s(3)
        sigmas_H3p(4) = BR_dd*s(3)
        sigmas_H3p(5) = BR_ee*s(4)
        sigmas_H3p(6) = BR_fd*s(3)
        sigmas_H3p(7) = BR_ge*s(4)
        sigmas_H3p(8) = BR_he*s(4)

        RETURN
    END FUNCTION get_sigmas_H3p


    FUNCTION get_sigma_H3p_H2H(E) result(sigma)
        !-------------------------------------------------------------------
        !    H3+ + H2 -> H2 + H + (H2+) 
        !       Input parameters
        !           E       energy in keV 
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E
        DOUBLE PRECISION :: sigmas_H3p(8)
        sigmas_H3p = get_sigmas_H3p(E)
        sigma = sigmas_H3p(1)
        RETURN
    END FUNCTION get_sigma_H3p_H2H

    FUNCTION get_sigma_H3p_3H(E) result(sigma)
        !-------------------------------------------------------------------
        !    H3+ + H2 -> H + H + H + (H2+) 
        !       Input parameters
        !           E       energy in keV 
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E
        DOUBLE PRECISION :: sigmas_H3p(8)
        sigmas_H3p = get_sigmas_H3p(E)
        sigma = sigmas_H3p(2)
        RETURN
    END FUNCTION get_sigma_H3p_3H

    FUNCTION get_sigma_H3p_H2pH(E) result(sigma)
        !-------------------------------------------------------------------
        !    H3+ + H2 -> H2+ + H + (H2) 
        !       Input parameters
        !           E       energy in keV 
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E
        DOUBLE PRECISION :: sigmas_H3p(8)
        sigmas_H3p = get_sigmas_H3p(E)
        sigma = sigmas_H3p(3)
        RETURN
    END FUNCTION get_sigma_H3p_H2pH

    FUNCTION get_sigma_H3p_H2Hp(E) result(sigma)
        !-------------------------------------------------------------------
        !    H3+ + H2 -> H2 + H+ + (H2) 
        !       Input parameters
        !           E       energy in keV 
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E
        DOUBLE PRECISION :: sigmas_H3p(8)
        sigmas_H3p = get_sigmas_H3p(E)
        sigma = sigmas_H3p(4)
        RETURN
    END FUNCTION get_sigma_H3p_H2Hp

    FUNCTION get_sigma_H3p_Hp2H(E) result(sigma)
        !-------------------------------------------------------------------
        !    H3+ + H2 -> H+ + H + H + (H2) 
        !       Input parameters
        !           E       energy in keV 
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E
        DOUBLE PRECISION :: sigmas_H3p(8)
        sigmas_H3p = get_sigmas_H3p(E)
        sigma = sigmas_H3p(5)
        RETURN
    END FUNCTION get_sigma_H3p_Hp2H

    FUNCTION get_sigma_H3p_H2pHp(E) result(sigma)
        !-------------------------------------------------------------------
        !    H3+ + H2 -> H2+ + H+ + (H2 + e-)
        !       Input parameters
        !           E       energy in keV 
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E
        DOUBLE PRECISION :: sigmas_H3p(8)
        sigmas_H3p = get_sigmas_H3p(E)
        sigma = sigmas_H3p(6)
        RETURN
    END FUNCTION get_sigma_H3p_H2pHp

    FUNCTION get_sigma_H3p_2HpH(E) result(sigma)
        !-------------------------------------------------------------------
        !    H3+ + H2 -> H+ + H+ + H + (H2 + e-)
        !       Input parameters
        !           E       energy in keV 
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E
        DOUBLE PRECISION :: sigmas_H3p(8)
        sigmas_H3p = get_sigmas_H3p(E)
        sigma = sigmas_H3p(7)
        RETURN
    END FUNCTION get_sigma_H3p_2HpH

    FUNCTION get_sigma_H3p_3Hp(E) result(sigma)
        !-------------------------------------------------------------------
        !    H3+ + H2 -> H+ + H+ + H+ + (H2 + 2e-)
        !       Input parameters
        !           E       energy in keV 
        !       Output parameters
        !           sigma   cross-section in m^-2
        !-------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION :: sigma
        DOUBLE PRECISION, INTENT(in) :: E
        DOUBLE PRECISION :: sigmas_H3p(8)
        sigmas_H3p = get_sigmas_H3p(E)
        sigma = sigmas_H3p(8)
        RETURN
    END FUNCTION get_sigma_H3p_3Hp


END MODULE 
