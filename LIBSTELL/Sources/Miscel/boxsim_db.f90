!-----------------------------------------------------------------------
!     Module:        boxsim_db
!     Authors:       L van Ham (lucas.van.ham@ipp.mpg.de)
!     Date:          02/12/2025
!     Description:   This module initializes the database for reactions_db 
!                    used in boxsim_sigma
!-----------------------------------------------------------------------
MODULE boxsim_db
  !-----------------------------------------------------------------------
  !     Libraries
  !-----------------------------------------------------------------------
  USE boxsim_sigma
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
  INTEGER, PARAMETER :: boxsim_max_products = 3
  INTEGER, PARAMETER :: boxsim_unused = -10
  TYPE :: box_reaction
      CHARACTER(len=64) :: name 
      INTEGER :: input_Q ! charge in
      INTEGER :: input_CAT ! family
      INTEGER :: input_A ! atomicity
      INTEGER :: nproducts ! Number of products out
      INTEGER :: output_Q(boxsim_max_products) = boxsim_unused! charge out
      INTEGER :: output_A(boxsim_max_products) = boxsim_unused! mass out
      LOGICAL :: enabled ! for testing purposes
      PROCEDURE(sigma_interface), POINTER, NOPASS :: calc_sigma => NULL()
  END TYPE box_reaction

  TYPE(box_reaction), ALLOCATABLE :: reactions_db(:) 
  INTEGER :: n_reactions = 0                    
  PUBLIC :: boxsim_init_reactions, reactions_db, n_reactions


  ! Particle kinds 
  INTEGER, PARAMETER :: boxsim_kind_H   = 1
  INTEGER, PARAMETER :: boxsim_kind_D   = 2
  INTEGER, PARAMETER :: boxsim_kind_T   = 3
  INTEGER, PARAMETER :: boxsim_kind_He3 = 4
  INTEGER, PARAMETER :: boxsim_kind_He4 = 5

  INTEGER, PARAMETER :: boxsim_nkinds   = 5 ! one for each of below
  INTEGER, PARAMETER :: boxsim_ncats    = 2

  INTEGER, PARAMETER :: boxsim_cat_H  = 1
  INTEGER, PARAMETER :: boxsim_cat_He = 2
  INTEGER, PARAMETER :: boxsim_cat_Z(boxsim_ncats) = [1, 2]

  INTEGER, PARAMETER :: boxsim_kind_cat(boxsim_nkinds) = &
                      [boxsim_cat_H, boxsim_cat_H, boxsim_cat_H, &
                       boxsim_cat_He, boxsim_cat_He]
  DOUBLE PRECISION, PARAMETER :: boxsim_kind_mass(boxsim_nkinds) = &
                      [ 1.67262192D-27, &
                        3.34358378D-27, &
                        5.00735675D-27, &
                        5.00641279D-27, &
                        6.64465735D-27 &
                        ]
CONTAINS
  !-------------------------------------------------------------------
  !    boxsim_INIT_REACTIONS : Initializes reactions_db database.
  !-------------------------------------------------------------------
  SUBROUTINE boxsim_init_reactions()

  IMPLICIT NONE

  ! Initialize database
  ALLOCATE(reactions_db(17)) ! H: 5, H2: 4, H3: 8(!)

  !-------------------------------------------------------------------
  !   SINGLE PROTON
  !-------------------------------------------------------------------

  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H+ + H2 -> fast H", &
      nproducts = 1, &
      input_Q   = 1, &
      input_A   = 1, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_Hp_H
  reactions_db(n_reactions)%output_A(1) = 1 ! H
  reactions_db(n_reactions)%output_Q(1) = 0 ! H

  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H + H2 -> fast H+", &
      nproducts = 1, &
      input_Q   = 0, &
      input_A   = 1, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H_Hp
  reactions_db(n_reactions)%output_A(1) = 1 ! H+
  reactions_db(n_reactions)%output_Q(1) = 1 ! H+

  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H + H2 -> fast H-", &
      nproducts = 1, &
      input_Q   = 0, &
      input_A   = 1, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H_Hm
  reactions_db(n_reactions)%output_A(1) = 1  ! H-
  reactions_db(n_reactions)%output_Q(1) = -1 ! H-

  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H- + H2 -> fast H", &
      nproducts = 1, &
      input_Q   = -1,  &
      input_A   = 1, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_Hm_H
  reactions_db(n_reactions)%output_A(1) = 1 ! H
  reactions_db(n_reactions)%output_Q(1) = 0 ! H

  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H- + H2 -> fast H+", &
      nproducts = 1, &
      input_Q   = -1,&
      input_A   = 1, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_Hm_Hp
  reactions_db(n_reactions)%output_A(1) = 1 ! H+
  reactions_db(n_reactions)%output_Q(1) = 1 ! H+

!-------------------------------------------------------------------
!   DOUBLE PROTON
!-------------------------------------------------------------------
  ! Neutralization of H2+ to produce fast H2
  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H2+ + H2 -> fast H2", &
      nproducts = 1, &
      input_Q   = 1, &
      input_A   = 2, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H2p_H2
  reactions_db(n_reactions)%output_A(1) = 2 ! H2
  reactions_db(n_reactions)%output_Q(1) = 0 ! H2

  ! Ionization of H2 to produce fast H2+
  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H2 + (H2) -> fast H2+", &
      nproducts = 1, &
      input_Q   = 0, &
      input_A   = 2, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H2_H2p
  reactions_db(n_reactions)%output_A(1) = 2 ! H2+
  reactions_db(n_reactions)%output_Q(1) = 1 ! H2+

  ! Dissociation of H2+ into H+ and H
  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H2+ + H2 -> fast H+, fast H", &
      nproducts = 2, &
      input_Q   = 1, &
      input_A   = 2, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H2p_HHp
  reactions_db(n_reactions)%output_A(1) = 1 ! H+
  reactions_db(n_reactions)%output_Q(1) = 1 ! H+
  reactions_db(n_reactions)%output_A(2) = 1 ! H
  reactions_db(n_reactions)%output_Q(2) = 0 ! H

  ! Dissociation of H2 into H+ and H
  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H2 + H2 -> fast H+, fast H", &
      nproducts = 2, &
      input_Q = 0, &
      input_A = 2, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H2_HHp
  reactions_db(n_reactions)%output_A(1) = 1 ! H+
  reactions_db(n_reactions)%output_Q(1) = 1 ! H+
  reactions_db(n_reactions)%output_A(2) = 1 ! H
  reactions_db(n_reactions)%output_Q(2) = 0 ! H

!-------------------------------------------------------------------
!   TRIPLE PROTON
!-------------------------------------------------------------------

  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H3+ + H2 -> H2 + H + H2+", &
      nproducts = 2, &
      input_Q   = 1, &
      input_A   = 3, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H3p_H2H
  reactions_db(n_reactions)%output_A(1) = 2 ! H2
  reactions_db(n_reactions)%output_Q(1) = 0 ! H2
  reactions_db(n_reactions)%output_A(2) = 1 ! H
  reactions_db(n_reactions)%output_Q(2) = 0 ! H

  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H3+ + H2 -> H + H + H + H2+", &
      nproducts = 3, &
      input_Q   = 1, &
      input_A   = 3, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H3p_3H
  reactions_db(n_reactions)%output_A(1) = 1 ! H
  reactions_db(n_reactions)%output_Q(1) = 0 ! H
  reactions_db(n_reactions)%output_A(2) = 1 ! H
  reactions_db(n_reactions)%output_Q(2) = 0 ! H
  reactions_db(n_reactions)%output_A(3) = 1 ! H
  reactions_db(n_reactions)%output_Q(3) = 0 ! H

  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H3+ + H2 -> H2+ + H + H2", &
      nproducts = 2, &
      input_Q   = 1, &
      input_A   = 3, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H3p_H2pH
  reactions_db(n_reactions)%output_A(1) = 2 ! H2+
  reactions_db(n_reactions)%output_Q(1) = 1 ! H2+
  reactions_db(n_reactions)%output_A(2) = 1 ! H
  reactions_db(n_reactions)%output_Q(2) = 0 ! H

  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H3+ + H2 -> H2 + H+ + H2", &
      nproducts = 2, &
      input_Q   = 1, &
      input_A   = 3, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H3p_H2Hp
  reactions_db(n_reactions)%output_A(1) = 2 ! H2
  reactions_db(n_reactions)%output_Q(1) = 0 ! H2
  reactions_db(n_reactions)%output_A(2) = 1 ! H+
  reactions_db(n_reactions)%output_Q(2) = 1 ! H+

  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H3+ + H2 -> H+ + H + H + H2", &
      nproducts = 3, &
      input_Q   = 1, &
      input_A   = 3, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H3p_Hp2H
  reactions_db(n_reactions)%output_A(1) = 1 ! H+
  reactions_db(n_reactions)%output_Q(1) = 1 ! H+
  reactions_db(n_reactions)%output_A(2) = 1 ! H
  reactions_db(n_reactions)%output_Q(2) = 0 ! H
  reactions_db(n_reactions)%output_A(3) = 1 ! H
  reactions_db(n_reactions)%output_Q(3) = 0 ! H

  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H3+ + H2 -> H2+ + H+ + H2 + e-", &
      nproducts = 2, &
      input_Q   = 1, &
      input_A   = 3, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H3p_H2pHp
  reactions_db(n_reactions)%output_A(1) = 2 ! H2+
  reactions_db(n_reactions)%output_Q(1) = 1 ! H2+
  reactions_db(n_reactions)%output_A(2) = 1 ! H+
  reactions_db(n_reactions)%output_Q(2) = 1 ! H+

  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H3+ + H2 -> H+ + H+ + H + H2 + e-", &
      nproducts = 3, &
      input_Q   = 1, &
      input_A   = 3, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H3p_2HpH
  reactions_db(n_reactions)%output_A(1) = 1 ! H+
  reactions_db(n_reactions)%output_Q(1) = 1 ! H+
  reactions_db(n_reactions)%output_A(2) = 1 ! H+
  reactions_db(n_reactions)%output_Q(2) = 1 ! H+
  reactions_db(n_reactions)%output_A(3) = 1 ! H
  reactions_db(n_reactions)%output_Q(3) = 0 ! H

  n_reactions = n_reactions + 1
  reactions_db(n_reactions) = box_reaction( & 
      name = "H3+ + H2 -> H+ + H+ + H+ + H2 + 2e-", &
      nproducts = 3, &
      input_Q   = 1, &
      input_A   = 3, &
      input_CAT = boxsim_cat_H, &
      enabled = .TRUE.)
  reactions_db(n_reactions)%calc_sigma => get_sigma_H3p_3Hp
  reactions_db(n_reactions)%output_A(1) = 1 ! H+
  reactions_db(n_reactions)%output_Q(1) = 1 ! H+
  reactions_db(n_reactions)%output_A(2) = 1 ! H+
  reactions_db(n_reactions)%output_Q(2) = 1 ! H+
  reactions_db(n_reactions)%output_A(3) = 1 ! H+
  reactions_db(n_reactions)%output_Q(3) = 0 ! H+

  END SUBROUTINE boxsim_init_reactions
      
  SUBROUTINE boxsim_vector_from_string(str, C, ierr)
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(in) :: str !! token string
  INTEGER, INTENT(OUT) :: C(boxsim_nkinds), ierr
  CHARACTER(LEN=8) :: tokens(boxsim_nkinds)
  CHARACTER(LEN=2) :: letters
  INTEGER :: ntokens, itoken, boxsim_count, mass_number, boxsim_kind
  LOGICAL :: ldash

  C = 0
  CALL tokenize_input(str)
  IF (ierr/=0) RETURN

  DO itoken = 1, ntokens
    CALL parse_token(tokens(itoken),boxsim_count)
    IF (ierr/=0) RETURN
    CALL token_to_kind()
    IF (ierr/=0) RETURN
    C(boxsim_kind) = C(boxsim_kind) + boxsim_count
  END DO

  RETURN

  CONTAINS 

    SUBROUTINE tokenize_input(str_in)

    CHARACTER(LEN=*), INTENT(in) :: str_in
    INTEGER :: n, tok_start, i
    CHARACTER(LEN=8) :: str_mod

    ierr = 0
    str_mod = TRIM(str_in)
    n = LEN_TRIM(str_mod)
    IF (n==0 .OR. .NOT. is_upper(str_mod(1:1))) THEN
      ierr = 1
      RETURN
    END IF

    ! Split
    ntokens = 0
    tok_start = 1
    DO i = 2, n
      IF (is_upper(str_mod(i:i))) THEN 
        ntokens = ntokens + 1
        tokens(ntokens) = str_mod(tok_start:i-1)
        tok_start = i
      END IF
    END DO
    ntokens = ntokens+1
    tokens(ntokens) = str_mod(tok_start:n)

    END SUBROUTINE tokenize_input

    SUBROUTINE parse_token(t,num)
    CHARACTER(LEN=*), INTENT(in) :: t
    INTEGER, INTENT(out) :: num
    INTEGER :: letter_end, dash_pos, n

    ierr = 0
    n = LEN_TRIM(t)
    letter_end = 1
    IF (n>=2) THEN
      IF (is_lower(t(2:2))) letter_end = 2
    END IF
    letters = t(1:letter_end)
    IF (letter_end<2) letters(2:2) = ' '

    num = 1
    mass_number = 0
    ldash = .FALSE.
    IF (n > letter_end) THEN
      dash_pos = INDEX(t(letter_end+1:n), '-')
      ldash = dash_pos > 0
      IF (ldash) THEN
        READ(t(letter_end+dash_pos+1:n),*,IOSTAT=ierr) mass_number
      ELSE
        READ(t(letter_end+1:n),*,IOSTAT=ierr) num
      END IF
      IF (ierr/=0 .OR. num<1) ierr = 1
    END IF

    END SUBROUTINE parse_token

    SUBROUTINE token_to_kind()
    ierr = 0
    boxsim_kind = 0

    SELECT CASE (letters)
      CASE ('H ') ! Protium
        IF (ldash) THEN
          ierr = 1
        ELSE
          boxsim_kind = boxsim_kind_H
        END IF
      CASE ('D ') ! Deuterium
        IF (ldash) THEN
          ierr = 1
        ELSE
          boxsim_kind = boxsim_kind_D
        END IF
      CASE ('T ') ! Tritium
        IF (ldash) THEN
          ierr = 1
        ELSE
          boxsim_kind = boxsim_kind_T
        END IF          
      CASE ('He') ! Helium 
        IF (.NOT.ldash .OR. mass_number==4) THEN
          boxsim_kind = boxsim_kind_He4
        ELSEIF (mass_number==3) THEN
          boxsim_kind = boxsim_kind_He3
        ELSE
          ierr = 1
        END IF
      CASE DEFAULT
        ierr = 1
    END SELECT

    END SUBROUTINE token_to_kind

    FUNCTION is_upper(car)
    CHARACTER(LEN=1), INTENT(in) :: car
    LOGICAL :: is_upper
    is_upper = (IACHAR(car)>=IACHAR('A')) .AND. (IACHAR(car)<=IACHAR('Z'))
    END FUNCTION is_upper

    FUNCTION is_lower(car)
    CHARACTER(LEN=1), INTENT(in) :: car
    LOGICAL :: is_lower
    is_lower = (IACHAR(car)>=IACHAR('a')) .AND. (IACHAR(car)<=IACHAR('z'))
    END FUNCTION is_lower

    
  END SUBROUTINE boxsim_vector_from_string

  SUBROUTINE boxsim_string_from_counts(part_counts, str, ierr)
!!-----------------------------------------------------------
!! Converts a part_counts array to a human-readable string
!!-----------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(in)            :: part_counts(boxsim_nkinds)
  CHARACTER(LEN=8), INTENT(out)  :: str
  INTEGER, INTENT(out)           :: ierr

  INTEGER :: ikind, icount
  CHARACTER(LEN=4)  :: piece, count_str
  CHARACTER(LEN=24) :: buffer

  ierr = 0
  buffer = ''
  DO ikind = 1, boxsim_nkinds ! loop over kinds
      IF (part_counts(ikind) <= 0) CYCLE

      SELECT CASE (ikind)
        CASE (boxsim_kind_H);   piece = 'H'
        CASE (boxsim_kind_D);   piece = 'D'
        CASE (boxsim_kind_T);   piece = 'T'
        CASE (boxsim_kind_He4); piece = 'He'
        CASE (boxsim_kind_He3); piece = 'He-3'
      END SELECT

      IF (ikind == boxsim_kind_He3) THEN
          DO icount = 1, part_counts(ikind)
              buffer = TRIM(buffer) // TRIM(piece)
          END DO
      ELSE IF (part_counts(ikind) > 1) THEN
          WRITE(count_str,'(I0)') part_counts(ikind)
          buffer = TRIM(buffer) // TRIM(piece) // TRIM(count_str)
      ELSE
          buffer = TRIM(buffer) // TRIM(piece)
      END IF
  END DO

  IF (LEN_TRIM(buffer) == 0) THEN
      ierr = 1   ! nothing present
      RETURN
  END IF
  IF (LEN_TRIM(buffer) > LEN(str)) THEN
      ierr = 2   ! does not fit 
      RETURN
  END IF
  str = TRIM(buffer)

  END SUBROUTINE boxsim_string_from_counts

  SUBROUTINE boxsim_strip_charge_suffix(str_in, comp_str, charge_int, ierr)
!!-------------------------------------------------------------------------
!! Parses particle charge from a string.
!!----------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(in)  :: str_in
  CHARACTER(LEN=8), INTENT(out) :: comp_str
  INTEGER, INTENT(out) :: charge_int, ierr
  CHARACTER(LEN=8) :: charge_tok
  CHARACTER(LEN=1) :: sign_char
  INTEGER :: n_total, ipos, n_tok, mag, ios

  ierr = 0
  charge_int = 0
  n_total = LEN_TRIM(str_in)
  IF (n_total == 0) THEN
      ierr = 1
      RETURN
  END IF

  ! Find space
  ipos = INDEX(str_in(1:n_total), ' ')

  IF (ipos == 0) THEN ! No space so neutral
      comp_str = str_in(1:n_total)
      charge_int = 0
      RETURN
  END IF

  comp_str  = str_in(1:ipos-1)
  charge_tok = str_in(ipos+1:n_total)

  IF (LEN_TRIM(comp_str) == 0 .OR. LEN_TRIM(charge_tok) == 0) THEN
      ierr = 1
      RETURN
  END IF

  n_tok = LEN_TRIM(charge_tok)
  sign_char = charge_tok(n_tok:n_tok)
  IF (sign_char /= '+' .AND. sign_char /= '-') THEN
      ierr = 1   ! bad symbol 
      RETURN
  END IF

  IF (n_tok == 1) THEN
      mag = 1   ! No number means magnitude 1
  ELSE
      READ(charge_tok(1:n_tok-1), *, IOSTAT=ios) mag
      IF (ios /= 0 .OR. mag < 1) THEN ! This better be a number or I'll break
          ierr = 1
          RETURN
      END IF
  END IF

  IF (sign_char == '+') THEN
      charge_int = mag
  ELSE
      charge_int = -mag
  END IF

  END SUBROUTINE boxsim_strip_charge_suffix


  SUBROUTINE boxsim_parse_species(str, part_counts, charge_int, Z_int, ierr)
!!-------------------------------------------------------------------
!! Fully parse a species (kind, charge) from a string
!!-------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(in) :: str
  INTEGER, INTENT(out) :: part_counts(boxsim_nkinds), charge_int,  Z_int, ierr
  CHARACTER(LEN=8) :: comp_str
  INTEGER :: cat

  CALL boxsim_strip_charge_suffix(str, comp_str, charge_int, ierr)
  IF (ierr /= 0) RETURN
  CALL boxsim_vector_from_string(comp_str, part_counts, ierr)
  IF (ierr /= 0) RETURN
  CALL boxsim_cat_from_counts(part_counts, cat, ierr)
  IF (ierr /= 0) RETURN
  Z_int = boxsim_cat_Z(cat)

  END SUBROUTINE boxsim_parse_species

  SUBROUTINE boxsim_cat_from_counts(part_counts, cat, ierr)
  !!-----------------------------------------------------
  !! Determines the category from a counts array (H or He)
  !!-----------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(in)  :: part_counts(boxsim_nkinds) !! counts array
  INTEGER, INTENT(out) :: cat !! category out
  INTEGER, INTENT(out) :: ierr !! error flag
  INTEGER :: ikind, this_cat
  INTEGER, PARAMETER :: no_cat_yet = 0

  ierr = 0
  cat = no_cat_yet
  DO ikind = 1, boxsim_nkinds ! loop over all possible kinds
      IF (part_counts(ikind) > 0) THEN
          this_cat = boxsim_kind_cat(ikind)
          IF (cat==no_cat_yet) THEN  ! first category set
            cat = this_cat
          ELSE IF (this_cat /= cat) THEN ! does not match
              ierr = 1   
              RETURN
          END IF
      END IF
  END DO
  IF (cat==no_cat_yet) ierr = 1   ! no cat???
  
  END SUBROUTINE boxsim_cat_from_counts

  SUBROUTINE boxsim_charge_suffix(charge_int, suffix, ierr)
!!-----------------------------------------------------------
!! Get charge suffix ('+','-','X+','X-', or blank for neutral).
!!-----------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(in) :: charge_int
  CHARACTER(LEN=4), INTENT(out) :: suffix
  INTEGER, INTENT(out) :: ierr
  INTEGER :: mag
  CHARACTER(LEN=3) :: mag_str
  
  ierr = 0
  mag = ABS(charge_int)
  
  suffix = '    ' ! Default neutral
  IF (charge_int == 0) RETURN  
  mag_str = '' ! Default 1
  IF (mag > 1) WRITE(mag_str,'(I0)') mag
  
  IF (charge_int > 0) THEN
      suffix = TRIM(mag_str)//'+'
  ELSE
      suffix = TRIM(mag_str)//'-'
  END IF
  
  END SUBROUTINE boxsim_charge_suffix


  SUBROUTINE boxsim_species_from_counts(part_counts, charge_int, str, ierr)
!!------------------------------------------------------------------------
!! Converts part_counts and charge into a string
!!------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(in)           :: part_counts(boxsim_nkinds), charge_int
  CHARACTER(LEN=8), INTENT(out) :: str
  INTEGER, INTENT(out) :: ierr
  CHARACTER(LEN=8) :: comp_str
  CHARACTER(LEN=4) :: suffix 
  INTEGER :: ierr2, needed_len

  CALL boxsim_string_from_counts(part_counts, comp_str, ierr)
  IF (ierr /= 0) RETURN
  CALL boxsim_charge_suffix(charge_int, suffix, ierr2)
  IF (ierr2 /= 0) THEN
      ierr = ierr2
      RETURN
  END IF

  IF (charge_int == 0) THEN
      str = TRIM(comp_str)
      RETURN
  END IF

  needed_len = LEN_TRIM(comp_str) + 1 + LEN_TRIM(suffix) 
  IF (needed_len > LEN(str)) THEN
      ierr = 2
      RETURN
  END IF
  str = TRIM(comp_str)//' '//TRIM(suffix)

  END SUBROUTINE boxsim_species_from_counts

  SUBROUTINE boxsim_split_counts(parent_counts, split_pattern, nproducts, & 
                                 product_counts, ierr)
!!------------------------------------------------------------------------
!! Returns the boxsim counts of all reaction products
!!------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(in) :: parent_counts(boxsim_nkinds)
  INTEGER, INTENT(in) :: split_pattern(*)
  INTEGER, INTENT(in) :: nproducts
  INTEGER, INTENT(out):: product_counts(boxsim_nkinds, nproducts)
  INTEGER, INTENT(out):: ierr

  INTEGER :: a, ikind, iatom, ipos, jpos, tmp, iprod, csr
  INTEGER, ALLOCATABLE :: bag(:)
  DOUBLE PRECISION :: r

  ierr = 0
  product_counts = 0
  a = SUM(parent_counts)
  IF ((a/=SUM(split_pattern(1:nproducts))).OR.(a==0)) THEN
    ierr = 1
    RETURN 
  END IF
  ALLOCATE(bag(a))
  ipos = 0
  DO ikind = 1, boxsim_nkinds
    DO iatom = 1, parent_counts(ikind)
      ipos = ipos + 1
      bag(ipos) = ikind
    END DO
  END DO

  ! Shuffle bag
  DO ipos = a, 2, -1
    CALL RANDOM_NUMBER(r)
    jpos = INT(r*ipos)+1
    tmp = bag(ipos)
    bag(ipos) = bag(jpos)
    bag(jpos) = tmp
  END DO
  ! Deal
  csr = 0
  DO iprod = 1, nproducts
    DO iatom = 1, split_pattern(iprod)
      csr = csr + 1
      ikind = bag(csr)
      product_counts(ikind,iprod) = product_counts(ikind,iprod)+1
    END DO
  END DO
  DEALLOCATE(bag)

  END SUBROUTINE boxsim_split_counts

  SUBROUTINE boxsim_ref_from_cat(cat, ref_M, ierr)
!!------------------------------------------------------------------------
!! Returns the reference mass of a family
!!------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(in) :: cat
  DOUBLE PRECISION, INTENT(out) :: ref_M
  INTEGER, INTENT(out) :: ierr

  ierr = 0
  ref_M = 0.0d0
  SELECT CASE (cat)
    CASE (boxsim_cat_H)
      ref_M = boxsim_kind_mass(boxsim_kind_H)
    CASE (boxsim_cat_He)
      ref_M = boxsim_kind_mass(boxsim_kind_He4)
    CASE DEFAULT
      ierr = 1
  END SELECT

  END SUBROUTINE boxsim_ref_from_cat
  
END MODULE boxsim_db