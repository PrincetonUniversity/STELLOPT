      SUBROUTINE chisq_vacopt (num, nopt, iflag, extension, lscreen)
      USE stel_kinds
      USE bsc_T
      USE chisq_mod
      USE coilsnamin
      USE safe_open_mod
      USE optim, ONLY: home_dir, exe_suffix, remove
      USE vacfield_mod, ONLY: numres, nstell_coils, sigma_angles_lo,
     1    sigma_angles_hi, sigma_shifts_lo, sigma_shifts_hi, vac_extcur
      USE system_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nopt
      INTEGER :: num, iflag
      CHARACTER(LEN=*) :: extension
      LOGICAL, INTENT(in) :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iunit=33, istat, no, mres, i, idum
      REAL(rprec) :: resid, dum, Target_b, sigma_b
      CHARACTER(LEN=200) :: cmd
      CHARACTER(LEN=2) :: screen
!***************************************
!
!  PHYSICS: Computes residues of (targeted) periodic orbits in vacuum magnetic field 
!           configuration using the tangent map method of Cary/Hanson
!
!  **** NOTE OF IMPORTANCE ****
!
!  This routine must be maintained in conjunction with VACOPT code
!
!****************************************************************************
!
      IF (nopt > 0) THEN
!
!     1. WRITE OUT NLST_VACOPT NAMELIST WITH ROTATION ANGLES AND SHIFTS
!     BASED ON xc_opt VALUES
!     xvacopt WILL USE THESE VALUES TO COMPUTE A ROTATED COILS FILE BEFORE
!     DOING THE ISLAND CALCULATION. SET A FLAG TO TELL XVACOPT TO DO THIS,
!     EITHER IN COMMAND LINE OR IN NAMELIST
!
!     LOAD angles, shifts FROM xc_opt ARRAY. SET ldisplace = .true IN NAMELIST
!     PASS NAME OF ORIGINAL coils_dot FILE, SINCE ALL ROTATIONS/DISPLACEMENTS
!     ARE RELATIVE TO THE ORIGINAL CONFIGURATION!!!
!      
!
!     1. CALL xvacopt TO WRITE OUT CHI-SQ FILE
!
         cmd = TRIM(home_dir) // 'xvacopt' // TRIM(exe_suffix)
         WRITE (screen,'(L2)') lscreen
         CALL load_physics_codes (cmd, 'input', screen, 'chisq_vacopt',
     1        extension, iunit, istat)

!
!     2. READ chisq_vacopt FILE AND FORM WEIGHTED CHISQ
!     

 !     a. island residues
         DO no = 1, numres
            READ (iunit,'(i3,1p,4e15.6)') mres, dum, dum, resid, sigma_b
            num = num+1
            chisq_match(num) = 0
            chisq_target(num)= resid
            wegt(num) = sigma_b
            index_array(num) = ivar_vacislandres
         END DO


!     b. angle, shift bounds
         DO no = 1, nstell_coils
            DO i = 1, 12             !(3 angles, 3 shifts, lo and hi bounds)
               READ (iunit,'(i3,1p,2e15.6)') idum, sigma_b, Target_b
               IF (ABS(sigma_b) .ne. zero) THEN
                  num = num+1
                  chisq_match(num)  = 0
                  chisq_target(num) = Target_b
                  wegt(num)         = sigma_b
                  index_array(num)  = ivar_coil_dis_bds
               END IF
            END DO
         END DO

         CLOSE (iunit)

      ELSE

!     COMPUTE NO. STELLARATOR SYMMETRIC CHI-SQ's
         DO no = 1, numres
            num = num+1
            IF (nopt .eq. -2) chisq_descript(num) = 
     1                        descript(ivar_vacislandres)
         END DO
         DO no = 1, nstell_coils
            i = 0
            sigma_b = sigma_angles_lo(no)%as_array(1)     !theta
            IF (ABS(sigma_b).gt.zero) i = i+1
            sigma_b = sigma_angles_lo(no)%as_array(2)     !phi
            IF (ABS(sigma_b).gt.zero) i = i+1
            sigma_b = sigma_angles_lo(no)%as_array(3)     !rot-angle
            IF (ABS(sigma_b).gt.zero) i = i+1
            sigma_b = sigma_angles_hi(no)%as_array(1)
            IF (ABS(sigma_b).gt.zero) i = i+1
            sigma_b = sigma_angles_hi(no)%as_array(2)
            IF (ABS(sigma_b).gt.zero) i = i+1
            sigma_b = sigma_angles_hi(no)%as_array(3)
            IF (ABS(sigma_b).gt.zero) i = i+1
            sigma_b = sigma_shifts_lo(no)%as_array(1)     !x-shift
            IF (ABS(sigma_b).gt.zero) i = i+1
            sigma_b = sigma_shifts_lo(no)%as_array(2)     !y-shift
            IF (ABS(sigma_b).gt.zero) i = i+1
            sigma_b = sigma_shifts_lo(no)%as_array(3)     !z-shift
            IF (ABS(sigma_b).gt.zero) i = i+1
            sigma_b = sigma_shifts_hi(no)%as_array(1)
            IF (ABS(sigma_b).gt.zero) i = i+1
            sigma_b = sigma_shifts_hi(no)%as_array(2)
            IF (ABS(sigma_b).gt.zero) i = i+1
            sigma_b = sigma_shifts_hi(no)%as_array(3)
            IF (ABS(sigma_b).gt.zero) i = i+1
            IF (i .gt. 0) THEN
               DO idum = 1, i
                  num = num+1
                  IF (nopt .eq. -2) chisq_descript(num) = 
     1                              descript(ivar_coil_dis_bds)
               ENDDO
            ENDIF
         ENDDO

         iflag = 0

      ENDIF

      END SUBROUTINE chisq_vacopt
