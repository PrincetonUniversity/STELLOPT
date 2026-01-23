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
        INTEGER :: output_Z(3), output_A(3) ! charge and mass out
        INTEGER :: nproducts
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
        INTEGER :: unused = -10

        ! Initialize database
        ALLOCATE(reactions_db(12))

    !-------------------------------------------------------------------
    !   SINGLE PROTON
    !-------------------------------------------------------------------

        ! CX of fast H+ with H2 to produce fast H
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H+ + H2 -> fast H", &
            nproducts = 1, &
            input_Z = 1, input_A = 1, &
            output_Z = [0, unused, unused], output_A = [1, unused, unused], &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_neut_Hplus

        ! Interaction of H with H2 to produce fast H+
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H + H2 -> fast H+", &
            nproducts = 1, &
            input_Z = 0, input_A = 1, &
            output_Z = [1, unused,unused], output_A = [1,unused,unused], &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_ionp_Hneut

        ! Interaction of H with H2 to produce fast H-
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H + H2 -> fast H-", &
            nproducts = 1, &
            input_Z = 0, input_A = 1, &
            output_Z = [-1, unused,unused], output_A = [1,unused,unused], &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_ionn_Hneut

        ! Detachment of H- electron with H2 to produce fast H
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H- + H2 -> fast H", &
            nproducts = 1, &
            input_Z = -1, input_A = 1, &
            output_Z = [0, unused,unused], output_A = [1,unused,unused], &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_neut_Hmin

        ! Double electron loss of H- with H2 to produce fast H+
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H- + H2 -> fast H+", &
            nproducts = 1, &
            input_Z = -1, input_A = 1, &
            output_Z = [1, unused,unused], output_A = [1,unused,unused], &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_ionp_Hmin

    !-------------------------------------------------------------------
    !   DOUBLE PROTON
    !-------------------------------------------------------------------
        ! Neutralization of H2+ to produce fast H2
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H2+ + H2 -> fast H2", &
            nproducts = 1, &
            input_Z = 1, input_A = 2, &
            output_Z = [0, unused,unused], output_A = [2,unused,unused], &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_neut_H2plus

        ! Ionization of H2 to produce fast H2+
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H2 + (H2) -> fast H2+", &
            nproducts = 1, &
            input_Z = 0, input_A = 2, &
            output_Z = [1, unused,unused], output_A = [2,unused,unused], &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_ionp_H2neut

        ! Dissociation of H2+ into H+ and H
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H2 + H2 -> fast H+, fast H", &
            input_Z = 1, input_A = 2, &
            nproducts = 2, &
            output_Z = [1, 0, unused], output_A = [1, 1, unused], &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_diss_H2plus

    !-------------------------------------------------------------------
    !   TRIPLE PROTON
    !-------------------------------------------------------------------
        ! Dissociation of H3+ forming H+ and H2
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H3+ + H2 -> H+ + H2 + H2", &
            nproducts = 2, &
            input_Z = 1, input_A = 3, &
            output_Z = [1, 2, unused], output_A = [1, 0, unused], &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_diss_H3plus_Hplus

        ! Dissociation of H3+ forming H2+ and H
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H3+ + H2 -> H2+ + H + H2", &
            nproducts = 2, &
            input_Z = 1, input_A = 3, &
            output_Z = [2, 1, unused], output_A = [1, 0, unused], &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_diss_H3plus_H2plus

        ! Dissociation of H3+ with charge exchange to gas
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H3+ + H2 -> H + H2 + H2+", &
            nproducts = 2, &
            input_Z = 1, input_A = 3, &
            output_Z = [1, 2, unused], output_A = [0, 0, unused], &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_diss_H3plus_neut

        ! Full breakup of H3+ forming H+, H, H
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H3+ + (H2) -> H+ + H + H + H2", &
            nproducts = 3, &
            input_Z = 1, input_A = 3, &
            output_Z = [1, 1, 1], output_A = [1, 0, 0], &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_diss_H3plus_triple

        END SUBROUTINE tabshi_init_reactions
        
END MODULE tabshi_db