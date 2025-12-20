!-----------------------------------------------------------------------
!     Module:        tabshi_db
!     Authors:       L van Ham (lucas.van.ham@ipp.mpg.de)
!     Date:          02/12/2025
!     Description:   This module initializes the database for reactions_db 
!                    used in tabshi_sigma
!-----------------------------------------------------------------------

MODULE tabshi_db
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE tabshi_sigma
    !-------------------------------------------------------------------
    !     Interface for calculating cross-sections
    !       Input parameter: E
    !       Output parameter: sigma
    !-------------------------------------------------------------------
    ABSTRACT INTERFACE 
        FUNCTION sigma_interface(E)   result(sigma)
            DOUBLE PRECISION, INTENT(in) :: E ! energy
            DOUBLE PRECISION :: sigma ! cross-section
        END FUNCTION sigma_interface
    END INTERFACE 
    !-------------------------------------------------------------------
    !     box_reaction: Reaction template
    !-------------------------------------------------------------------
    TYPE :: box_reaction
        CHARACTER(len=32) :: name 
        INTEGER :: input_Z, input_A ! charge and mass in 
        INTEGER :: output_Z_1, output_A_1 ! charge and mass out
        LOGICAL :: enabled ! for testing purposes
        PROCEDURE(sigma_interface), POINTER, NOPASS :: calc_sigma => NULL()
    END TYPE box_reaction

    TYPE(box_reaction), ALLOCATABLE :: reactions_db(:) 
    INTEGER :: n_reactions = 0                       
    PUBLIC :: tabshi_init_reactions, reactions_db, n_reactions

CONTAINS
    !-------------------------------------------------------------------
    !    TABSHI_INIT_REACTIONS : Initializes reactions_db database.
    !-------------------------------------------------------------------
        SUBROUTINE tabshi_init_reactions()

        IMPLICIT NONE

        ! Initialize database
        ALLOCATE(reactions_db(12))

    !-------------------------------------------------------------------
    !   SINGLE PROTON
    !-------------------------------------------------------------------

        ! CX of fast H+ with H2 to produce fast H
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H+ + H2 -> fast H", &
            input_Z = 1, input_A = 1, &
            output_Z_1 = 0, output_A_1 = 1, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_neut_Hplus

        ! Interaction of H with H2 to produce fast H+
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H + H2 -> fast H+", &
            input_Z = 0, input_A = 1, &
            output_Z_1 = 1, output_A_1 = 1, &            
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_ionp_Hneut

        ! Interaction of H with H2 to produce fast H-
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H + H2 -> fast H-", &
            input_Z = 0, input_A = 1, &
            output_Z_1 = -1, output_A_1 = 1, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_ionn_Hneut

        ! Detachment of H- electron with H2 to produce fast H
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H- + H2 -> fast H", &
            input_Z = -1, input_A = 1, &
            output_Z_1 = 0, output_A_1 = 1, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_neut_Hmin

        ! Double electron loss of H- with H2 to produce fast H+
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H- + H2 -> fast H+", &
            input_Z = -1, input_A = 1, &
            output_Z_1 = 1, output_A_1 = 1, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_ionp_Hmin

    !-------------------------------------------------------------------
    !   DOUBLE PROTON
    !-------------------------------------------------------------------
        ! Neutralization of H2+ to produce fast H2
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H2+ + H2 -> fast H2", &
            input_Z = 1, input_A = 2, &
            output_Z_1 = 0, output_A_1 = 2, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_neut_H2plus

        ! Ionization of H2 to produce fast H2+
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H2 + H2 -> fast H2+", &
            input_Z = 0, input_A = 2, &
            output_Z_1 = 1, output_A_1 = 2, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_ionp_H2neut

        ! Dissociation of H2+ into H+ and H
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H2 + H2 -> fast H+, fast H", &
            input_Z = 1, input_A = 2, &
            output_Z_1 = 1, output_A_1 = 1, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_diss_H2plus

    !-------------------------------------------------------------------
    !   TRIPLE PROTON
    !-------------------------------------------------------------------
        ! Kinetic dissociation of H3+ forming H+ and H2
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H3+ + H2 -> H+ + H2 + H2", &
            input_Z = 1, input_A = 3, &
            output_Z_1 = 1, output_A_1 = 1, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_kindiss_H3neut_Hplus

        ! Kinetic dissociation of H3+ forming H2+ and H
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H3+ + H2 -> H2+ + H + H2", &
            input_Z = 1, input_A = 3, &
            output_Z_1 = 1, output_A_1 = 2, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_kindiss_H3neut_H2plus

        ! Dissociation of H3+ after charge exchange forming H and H2+
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H3+ + (H2) -> H + H2+", &
            input_Z = 1, input_A = 3, &
            output_Z_1 = 0, output_A_1 = 1, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_cxdiss_H3neut_Hneut

        ! Dissociation of H3+ after charge exchange forming H2 and H+
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H3+ + (H2) -> H + H2+", &
            input_Z = 1, input_A = 3, &
            output_Z_1 = 0, output_A_1 = 2, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_cxdiss_H3neut_H2neut

        END SUBROUTINE tabshi_init_reactions
        
END MODULE tabshi_db