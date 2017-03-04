      SUBROUTINE read_coilgeom_elem(iunit, num, nwant, weight, lwant,
     1                              ivars)
      USE stel_kinds
      USE chisq_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: iunit, nwant, ivars
      INTEGER :: num
      REAL(rprec), DIMENSION(nwant), INTENT(in) :: weight
      LOGICAL, INTENT(in) :: lwant
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: icount, istat, nums
      REAL(rprec), PARAMETER :: one = 1, epsm = EPSILON(one)
      REAL(rprec), DIMENSION (nwant) :: values
      CHARACTER(LEN=200) :: parse_line
C-----------------------------------------------
!
!     NOTE: THE VALUES ALREADY HAVE THE WEIGHTS ON THEM
!     TO BE CONSISTENT WITH REST OF STELLOPT, WE DIVIDE THEM OUT HERE
!     AND MULTIPLY THEM BACK LATER
!
!
!     MAY BE READING IN MATRIX OF VALUES, LINE BY LINE
!     SO ICOUNT VALUE MAY ONLY EXIST FIRST TIME THROUGH
!
      READ (iunit,'(a)') parse_line

      READ (parse_line, *, iostat=istat) icount
      IF (istat .eq. 0) THEN
         READ (iunit,*, iostat=istat) values
      ELSE
         READ (parse_line, *, iostat=istat) values
      END IF

      IF (istat .ne. 0) THEN
         WRITE (6, *) 'ISTAT = ', istat,' in read_coilgeom_elem'
         STOP 
      END IF

      IF (lwant) THEN
         nums = num+nwant
         wegt(num+1:nums) = one/(weight(:nwant) + epsm)
         chisq_match(num+1:nums) = values(:nwant) * wegt(num+1:nums)
         chisq_target(num+1:nums) = 0
         index_array(num+1:nums) = ivars
         num = nums
      ENDIF

      END SUBROUTINE read_coilgeom_elem
