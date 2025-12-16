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
        INTEGER :: output_Z, output_A ! charge and mass out
        LOGICAL :: enabled ! for testing purposes
        PROCEDURE(sigma_interface), POINTER :: calc_sigma => NULL()
    END TYPE box_reaction

    TYPE(box_reaction), ALLOCATABLE :: reactions_db(:) 
    INTEGER :: n_reactions = 0                       
    PUBLIC :: init_reactions, reactions_db, n_reactions

CONTAINS
    !-------------------------------------------------------------------
    !    TABSHI_INIT_REACTIONS : Initializes reactions_db database.
    !-------------------------------------------------------------------
        SUBROUTINE tabshi_init_reactions()

        IMPLICIT NONE

        ! Initialize database
        ALLOCATE(reactions_db(5))

        ! CX of fast H+ with H2 to produce fast H
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H+ + H2 -> fast H", &
            input_Z = 1, input_A = 1, &
            output_Z = 0, output_A = 1, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_neut_Hplus

        ! Interaction of H with H2 to produce fast H+
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H + H2 -> fast H+", &
            input_Z = 0, input_A = 1, &
            output_Z = 1, output_A = 1, &            
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_ionp_Hneut

        ! Interaction of H with H2 to produce fast H-
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H + H2 -> fast H-", &
            input_Z = 0, input_A = 1, &
            output_Z = -1, output_A = 1, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_ionn_Hneut

        ! Detachment of H- electron with H2 to produce fast H
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H- + H2 -> fast H", &
            input_Z = -1, input_A = 1, &
            output_Z = 0, output_A = 1, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_neut_Hmin

        ! Double electron loss of H- with H2 to produce fast H+
        n_reactions = n_reactions + 1
        reactions_db(n_reactions) = box_reaction( & 
            name = "H- + H2 -> fast H+", &
            input_Z = -1, input_A = 1, &
            output_Z = 1, output_A = 1, &
            enabled = .TRUE.)
        reactions_db(n_reactions)%calc_sigma => get_sigma_ionp_Hmin

        END SUBROUTINE tabshi_init_reactions
        
END MODULE tabshi_db