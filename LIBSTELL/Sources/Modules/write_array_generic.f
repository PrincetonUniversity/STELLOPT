      MODULE write_array_generic
      USE stel_kinds
C-----------------------------------------------
C   L o c a l  V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
      INTEGER :: start_count, end_count, length, line
      CHARACTER(LEN=*), PARAMETER :: adyes = "yes", adno = "no"
      CHARACTER(LEN=3) :: adv
      
      INTERFACE write_array
         MODULE PROCEDURE write_array_section_real,
     1                    write_array_section_logical,
     2                    write_array_section_integer
      END INTERFACE
 
      PRIVATE :: parse_name

      CONTAINS

      SUBROUTINE write_array_section_real(iunit, name, array, n1,  
     1                                    ndim2, low_index)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: iunit, n1
      INTEGER, INTENT(in), OPTIONAL :: ndim2, low_index
      CHARACTER(len=*), INTENT(in) :: name
      REAL(rprec), DIMENSION(n1), INTENT(in) :: array
      REAL(rprec) :: temp
      CHARACTER(len=LEN_TRIM(name)+20) :: parse
      CHARACTER(len=40+LEN_TRIM(name)) :: line1,line2
C-----------------------------------------------
!
!     this subroutine writes repeated namelist entries in n*value format to save space
!     in namelist file connected to "iunit"
!
!     n1 = extent of 1st dimension of array
!     ndim2 = row (2nd dimension) of array, if present
!
      parse = name
      CALL parse_name (parse, n1, ndim2, low_index)

      IF (n1 <= 0) THEN
         RETURN
      ELSE IF (n1 == 1) THEN
         WRITE (line2, 990) array(1)
         line1 = TRIM(parse) // ADJUSTL(line2) 
         WRITE (iunit,980) TRIM(line1)

      ELSE IF (ALL(array(1) == array(2:n1))) THEN
         WRITE (line1, *) n1
         WRITE (line2, 990) array(1)
         line1 = TRIM(parse) // TRIM(ADJUSTL(line1)) 
     1         // '*' // ADJUSTL(line2)
         WRITE (iunit,980) TRIM(line1)

      ELSE
!       Look for repeats within the array
         start_count = 1
         line = 1
         WRITE(iunit, 980, advance="no") TRIM(parse)
         DO WHILE (start_count .le. n1)
            temp = array(start_count)
            end_count = start_count
            DO WHILE (end_count .lt. n1)
               IF (array(end_count+1) .ne. temp) EXIT
               end_count = end_count+1
            END DO
!              Limit no. records/line and start new record for each array
!              If maximum packing desired, eliminate the end_count==n in following test
            IF (line==3 .or. end_count.eq.n1) THEN
               adv = "yes"
               line = 0
            ELSE
               adv = "no"
               line = line+1
            END IF
            length = end_count - start_count + 1
            IF (length > 1) THEN
               WRITE (line1, *) length
               WRITE (line2, 990) temp
               line1 = TRIM(ADJUSTL(line1)) // '*' // ADJUSTL(line2)
               WRITE (iunit,980, advance=TRIM(adv)) TRIM(line1)
            ELSE
               WRITE (line1, 990) temp
               WRITE (iunit, 980, advance=TRIM(adv))TRIM(ADJUSTL(line1))
            END IF
            start_count = end_count+1
         END DO

  980 FORMAT(2x,a)
  990 FORMAT(1p,e21.14)

      END IF

      END SUBROUTINE write_array_section_real


      SUBROUTINE write_array_section_logical(iunit, name, larray, n1, 
     1                                       ndim2, low_index)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: iunit, n1
      INTEGER, INTENT(in), OPTIONAL :: ndim2, low_index
      CHARACTER(LEN=*), INTENT(in) :: name
      LOGICAL, DIMENSION(n1), INTENT(in) :: larray
      LOGICAL :: ltemp
      CHARACTER(len=LEN_TRIM(name)+20) :: parse
C-----------------------------------------------
!
!     this subroutine writes repeated namelist entries in n*value format to save space
!     in namelist file connected to "iunit"
!
!     n1 = extent of 1st dimension of array
!     ndim2 = row (2nd dimension) of array, if present
!
      parse = name
      CALL parse_name (parse, n1, ndim2, low_index)

      IF (n1 <= 0) THEN
         RETURN
      ELSE IF (n1 == 1) THEN
         WRITE(iunit,990) TRIM(parse), larray(1)

      ELSE IF (ALL(larray(1) .eqv. larray(2:n1))) THEN
          WRITE(iunit,993) TRIM(parse), n1, larray(1)

      ELSE
!       Look for repeats within the array
         start_count = 1
         line = 1
         WRITE(iunit, 980, advance="no") TRIM(parse)
         DO WHILE (start_count .le. n1)
            ltemp = larray(start_count)
            end_count = start_count
            DO WHILE (end_count .lt. n1)
               IF (larray(end_count+1) .neqv. ltemp) EXIT
               end_count = end_count+1
            END DO
!              Limit no. records/line and start new record for each array
!              If maximum packing desired, eliminate the end_count==n in following test
            IF (line==20 .or. end_count.eq.n1) THEN
               adv = "yes"
               line = 0
            ELSE
               adv = "no"
               line = line+1
            END IF
            length = end_count - start_count + 1
            IF (length > 1) THEN
               WRITE (iunit, 982, advance=TRIM(adv)) length, ltemp
            ELSE
               WRITE (iunit, 984, advance=TRIM(adv)) ltemp
            END IF
            start_count = end_count+1
         END DO

  980 FORMAT(2x,a)
  982 FORMAT(2x,i4,'*',l1)
  984 FORMAT(2x,l1)
  990 FORMAT(2x,a,l1)
  993 FORMAT(2x,a,i4,'*',l1)

      END IF

      END SUBROUTINE write_array_section_logical

      SUBROUTINE write_array_section_integer(iunit, name, iarray, n1, 
     1                                       ndim2, low_index)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: iunit, n1
      INTEGER, INTENT(in), OPTIONAL :: ndim2, low_index
      CHARACTER(len=*), INTENT(in) :: name
      INTEGER, DIMENSION(n1), INTENT(in) :: iarray
!      INTEGER, DIMENSION(:), INTENT(in) :: iarray
      INTEGER :: temp
      CHARACTER(len=LEN_TRIM(name)+20)   :: parse
      CHARACTER(len=40+LEN_TRIM(name))   :: line1,line2
C-----------------------------------------------
!
!     this subroutine writes repeated namelist entries in n*value format to save space
!     in namelist file connected to "iunit"
!
!     n1 = extent of 1st dimension of array
!     ndim2 = row (2nd dimension) of array, if present
!
      parse = name
      CALL parse_name (parse, n1, ndim2, low_index)

      IF (n1 <= 0) THEN
         RETURN
      ELSE IF (n1 == 1) THEN
         WRITE (line1, *) iarray(1)
         line2 = TRIM(parse) // ADJUSTL(line1)
         WRITE (iunit, 980) TRIM(line2)

      ELSE IF (ALL(iarray(1) == iarray(2:n1))) THEN
        WRITE (line1, *) n1
        WRITE (line2, *) iarray(1)
        line1 = TRIM(parse) // TRIM(ADJUSTL(line1)) 
     1        // '*' // ADJUSTL(line2)
        WRITE (iunit,980) TRIM(line1)

      ELSE
!       Look for repeats within the array
         start_count = 1
         line = 1
         WRITE(iunit, 980, advance="no") TRIM(parse)
         DO WHILE (start_count .le. n1)
            temp = iarray(start_count)
            end_count = start_count
            DO WHILE (end_count .lt. n1)
               IF (iarray(end_count+1) .ne. temp) EXIT
               end_count = end_count+1
            END DO
!              Limit no. records/line and start new record for each array
!              If maximum packing desired, eliminate the end_count==n in following test
            IF (line==8 .or. end_count.eq.n1) THEN
               adv = "yes"
               line = 0
            ELSE
               adv = "no"
               line = line+1
            END IF
            length = end_count - start_count + 1
            IF (length > 1) THEN
               WRITE (line1, *) length
               WRITE (line2, *) temp
               line1 = TRIM(ADJUSTL(line1)) // '*' // ADJUSTL(line2)
               WRITE (iunit,980, advance=TRIM(adv)) TRIM(line1)
            ELSE
               WRITE (line1, *) temp
               WRITE (iunit, 980, advance=TRIM(adv)) 
     1               TRIM(ADJUSTL(line1))
            END IF
            start_count = end_count+1
         END DO

  980 FORMAT(2x,a)

      END IF

      END SUBROUTINE write_array_section_integer


      SUBROUTINE parse_name (name, n1, ndim2, low_index)
      CHARACTER(LEN=*), INTENT(inout) :: name
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in), OPTIONAL :: ndim2, low_index
      CHARACTER(LEN=*), PARAMETER :: fmt(5) = (/ '(a,i1)',
     1       '(a,i2)', '(a,i3)', '(a,i4)', '(a,i5)' /)
      CHARACTER(len=LEN(name)+100) :: parse
      CHARACTER(len=100)           :: fmt_string
      INTEGER :: i1, istart, iend

      istart = 1
      IF (PRESENT(low_index)) istart = low_index 
      i1 = IndexFormat(istart)
      WRITE (parse, fmt(i1)) TRIM(name) // "(", istart 
      
      iend = istart + n1 - 1
      i1 = IndexFormat(iend)
      WRITE (parse, fmt(i1)) TRIM(parse) // ":", iend

      IF (PRESENT(ndim2)) THEN
         i1 = IndexFormat(ndim2)
         WRITE (parse, fmt(i1)) TRIM(parse) // ",", ndim2
      END IF
      
      name = TRIM(parse) // ")="

      END SUBROUTINE parse_name

      INTEGER FUNCTION IndexFormat(number)
      INTEGER, INTENT(in) :: number
      INTEGER :: i1

      SELECT CASE (ABS(number))
      CASE (1000:9999)
          i1 = 4
      CASE (100:999)
          i1 = 3
      CASE (10:99)
          i1 = 2
      CASE DEFAULT
          i1 = 1
      END SELECT

      IF (number < 0) i1 = i1 + 1
      IndexFormat = i1

      END FUNCTION IndexFormat

      END MODULE write_array_generic
