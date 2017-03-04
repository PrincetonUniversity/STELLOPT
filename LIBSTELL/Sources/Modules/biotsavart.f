      MODULE biotsavart
      USE stel_kinds
      USE bsc_T
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nfp_bs
      TYPE (bsc_coil),     POINTER            :: single_coil => null()
      TYPE (bsc_coilcoll), DIMENSION(:), ALLOCATABLE, TARGET ::
     &   coil_group

!
!     coil_group            Collection of coils, used to store coils by group id    
!     single coil           A single coil
!
      CONTAINS

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------

      SUBROUTINE initialize_biotsavart (extcur_in, extension, 
     1                                  xpt, scaled)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), INTENT(in), DIMENSION(:) :: extcur_in
      REAL(rprec), DIMENSION(:,:), OPTIONAL :: xpt
      LOGICAL, OPTIONAL                     :: scaled
      CHARACTER(LEN=*), OPTIONAL            :: extension
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      TYPE(bsc_coil)     :: coil_temp
      INTEGER            :: istat
      INTEGER     :: nextcur, ig, nc
      REAL(rprec)        :: current, current_first
      CHARACTER(len=200) :: coil_file
      LOGICAL            :: scaled_or_raw
!-----------------------------------------------
      scaled_or_raw = .true.
      IF (PRESENT(scaled)) scaled_or_raw = scaled

      IF (PRESENT(extension)) THEN
!  Parse coils.extension file and initialize bsc routines
         coil_file = 'coils.' // TRIM(extension)      
         CALL parse_coils_file (TRIM(coil_file))
!  Set currents in coils in each group
!  Note: current_first is read in from "coil_file"
         nextcur = SIZE(coil_group)
         IF (scaled_or_raw) THEN
            DO ig = 1, nextcur
               DO nc = 1, coil_group(ig) % ncoil
                  current = coil_group(ig) % coils(nc) % current
                  IF (nc .eq. 1) current_first = current
                  IF (current_first .ne. zero) 
     1            coil_group(ig) % coils(nc) % current =
     1            (current/current_first)*extcur_in(ig)
               END DO
            END DO
         END IF

      ELSE IF (PRESENT(xpt)) THEN
!  Initialize biotsavart MODULE
         nextcur = 1
         CALL cleanup_biotsavart
!  Create the coil (DO NOT CLOSE IT!)
         nc = SIZE(xpt,2)
         ALLOCATE (single_coil)
         CALL bsc_construct(single_coil,'fil_loop','','',                      &
     &         extcur_in(1),xpt(1:3,1:nc))
       
      ELSE
         STOP 'Fatal: initialize_bs: xpt or extension must be specified'
      END IF

      END SUBROUTINE initialize_biotsavart

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------

      SUBROUTINE parse_coils_file (coil_file, lgrps)
      USE safe_open_mod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*)      :: coil_file
      LOGICAL, OPTIONAL     :: lgrps           
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: iou_coil0=22
      INTEGER, PARAMETER :: maxgroups = 1000
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER            :: istat, i_line
      INTEGER, DIMENSION(:) ::  ngroup(maxgroups)
      INTEGER     :: nmaxnodes, iou_coil, nextcur
      INTEGER     :: n_line_start_string, n_line_skip
      CHARACTER(LEN=200) :: line, group, line_lc
      CHARACTER(LEN=28) :: start_string = '** coils_dot_starts_below **'
      LOGICAL            :: local_lgrps
!-----------------------------------------------
!
!  OPEN COILS.EXT (COIL_FILE)
!
      iou_coil = iou_coil0
      CALL safe_open(iou_coil, istat, TRIM(coil_file), 'old',                  &
     &   'formatted')
      IF (istat .ne. 0) STOP 'Error opening input coil file'

!  start_string added JDH, April 2011.
!  Find the line number where the start_string occurs (if at all)
!  Leave the file positioned at the line after the start_string
!  (or at the beginning of the file, if the start_string does not occur)
      i_line = 0
      n_line_start_string = 0
      read_loop : DO
         i_line = i_line + 1
         READ(iou_coil,'(a)',END = 100, IOSTAT=istat) line
         IF (istat .ne. 0) THEN
            WRITE(6,*) 'Problem in parse_coils_file. istat =',istat
            WRITE(6,*) ' Line number is ', i_line
            WRITE(6,*) line
            STOP
         END IF
         line_lc = line
         CALL tolower(line_lc)
         istat = INDEX(line, start_string)
         IF (istat .eq. 0) THEN  ! did not find start_string
            CYCLE
         ELSE   !  did find start_string
            n_line_start_string = i_line
            WRITE(6,*) 'Found start_string: ', start_string
            WRITE(6,*) ' in line ', n_line_start_string
            EXIT
         ENDIF
      END DO read_loop
      
100   IF (n_line_start_string .eq. 0) THEN
         REWIND (iou_coil, IOSTAT=istat)
         IF (istat .ne. 0) THEN
            WRITE(6,*) 'Problem 2 in parse_coils_file. istat =',istat
            STOP
         END IF
      END IF
!    Define number of lines to skip on subsequent re-reading
      n_line_skip = 3 + n_line_start_string

!  READ IN NUMBER OF FIELD PERIODS
      READ (iou_coil, '(a)' , iostat=istat) line
      istat = INDEX(line, 'periods')
      IF (istat .eq. 0)                                                        &
     &   STOP 'First line of coils file must contain # periods'
      READ (line, *, iostat=istat) group, nfp_bs

      local_lgrps = .false.
      IF (PRESENT(lgrps)) local_lgrps = lgrps
      
!   First pass through the input file, to find out how many coil groups, 
!   and the maximum number of nodes in ANY one coil.
      CALL read_coils_pass1(iou_coil, nextcur, nmaxnodes, ngroup,              &
     &                      local_lgrps, n_line_skip)

!   Now that we know how many coil groups there are, ALLOCATE
!   the bsc_coilcoll array.
      IF (nextcur .gt. 0) THEN
         CALL cleanup_biotsavart
         ALLOCATE (coil_group(nextcur), stat=istat)
         IF (istat .ne. 0) STOP 'ERROR ALLOCATION COIL COLLECTION'
      ELSE
         WRITE(*,*) 'number coilgroups = ',nextcur,' <= 0 '
         STOP
      END IF

!   Second pass through the file. Read ALL the coils, and store
!   them in the appropriate groups.
      CALL read_coils_pass2(iou_coil, nmaxnodes, coil_group, ngroup,           &
     &                      local_lgrps, n_line_skip)
 
      END SUBROUTINE parse_coils_file

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------

      SUBROUTINE read_coils_pass1 (iou, n_coilgroups, nmaxnodes, ngroup,       & 
     &                             lgrps, n_skip)
!     Subroutine to do a first pass read of the "coils" file. The purpose is
!     to find the number of coil groups, and the maximum number of
!     nodes in ANY one coil.

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)  :: iou, n_skip
      INTEGER, INTENT(out) :: n_coilgroups, nmaxnodes
      INTEGER, DIMENSION(:), INTENT(out) :: ngroup
      LOGICAL, INTENT(in)  :: lgrps

!     iou           Default integer - Fortran I/O unit number
!     n_coilgroups  Number of external current coil groups.
!                   (In MAKEGRID this is nextcur. Read as n_extcur, not next_cur)
!     n_skip        Number of line to skip at the beginning of the file
!     nmaxnodes     Maximum number of nodes in ANY one coil
!     ngroup        Array of coil group numbers
!                      igroup, when lgrps is False
!                      Sequential, when lgrps is True
!     lgrps         Logical.
!                   True - Each separate coil is stored in its own coil group
!                   False - Placement into coil groups determined by integer igroup
!                      read from the coils. file
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(200) :: line
      CHARACTER(100) :: group_id
      REAL(rprec)    :: xw, yw, zw
      REAL(rprec)    :: currin
      INTEGER        :: i, igroup, inodes
      INTEGER        :: istat
      LOGICAL        :: lparsed

!     line          Used for reading (and rereading) a line from the input file
!     group_id      Identifier of an external current coil group
!                   In MAKEGRID, this is variable group
!     xw, yw, zw    Scalars used for parsing Cartisian coordinates of wires
!     currin        current
!     igroup        Integer, number of an external current coil group
!     inodes        Integer, counter for the number of nodes in a coil
!     istat         Default integer, for I/O status
!     lparsed       Indicates previous line was parsed (used in case "end" statement missing)
!-----------------------------------------------

!  Start of Executable Code
!     Read in first n_skip lines
      REWIND (iou, IOSTAT=istat)
      DO i = 1,n_skip
         READ(iou,'(a)') line
      END DO
    
!     Initialize counters
      inodes = 0
      nmaxnodes = 0
      n_coilgroups = 0
      ngroup(:) = -1

!     Loop to read in the rest of the lines from the file
!     Last line of the file should be just 'end'.
      read_loop : DO
         READ(iou,'(a)',END = 100, IOSTAT=istat) line
         IF (istat .ne. 0) THEN
            WRITE(6,*) 'Problem in read_coils_pass1. istat =',istat
            WRITE(6,*) line
            STOP
         END IF
         
         IF (line(1:3) .eq. 'end') EXIT
         inodes = inodes + 1

!     Reread the line, assuming igroup and group_id are there.
!     For most of the file, they won't be there, and istat will be nonzero
!     But, when igroup and group_id are there, istat will be zero.
!     This is the indication of the END of a coil 
!
!     Note: IGROUP numbering may NOT be contiguous in file NOR does it
!           have to be sequential (1,2,3). Thus, igroup = 55,11,101,... 
!           MIGHT occur in file and should be accounted for in this logic.

         READ(line,*,iostat=istat) xw, yw, zw, currin, igroup, group_id
         lparsed = (istat .eq. 0)
         IF (lparsed) THEN 
            i = MINVAL(ABS(igroup - ngroup(:)))
            IF (i.ne.0 .or. lgrps) THEN       !ADD NEW GROUP AND STORE IN NGROUP
               n_coilgroups = n_coilgroups + 1
               IF (n_coilgroups .gt. SIZE(ngroup)) THEN
                  STOP ' read_coils_pass1: coil groups > SIZE(ngroup)'
               ENDIF
               IF (lgrps) THEN 
                  ngroup(n_coilgroups) = n_coilgroups
               ELSE ! lgrps is False
                  ngroup(n_coilgroups) = igroup
               END IF
            END IF
            nmaxnodes = MAX(nmaxnodes,inodes)
            inodes = 0
         END IF
      END DO read_loop

!
!     CATCH OLD STYLE (ONE BIG COIL) FILE
!
      IF (nmaxnodes .eq. 0) THEN
         nmaxnodes = inodes
         n_coilgroups = 1
         ngroup(1) = 1
      END IF
      
      RETURN

100   IF (.not. lparsed) THEN
         WRITE(6,*) 'Problems in read_coils_pass1'
         WRITE(6,*) 'EOF reached before END'
         WRITE(6,*) 'Make sure last line of file is "end"' 
      END IF
      
      END SUBROUTINE read_coils_pass1
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------

      SUBROUTINE read_coils_pass2 (iou, nmaxnodes, coil_group, ngroup,         &
     &                             lgrps, n_skip)
!     Subroutine to do a second pass read of the coils file. The data is then 
!     used tocreate a filamentary loop coil (fil_loop). The fil_loop is then appended 
!     to the correct coil group.

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: iou, n_skip
      INTEGER, INTENT(in) :: nmaxnodes
      INTEGER, DIMENSION(:), INTENT(in) :: ngroup
      TYPE (bsc_coilcoll), DIMENSION(:), INTENT(inout) :: coil_group
      LOGICAL, INTENT(in) :: lgrps

!     iou           Default Integer - Fortran I/O unit number
!     nmaxnodes     Maximum number of nodes in ANY one coil
!     ngroup        Array of integer coil group numbers
!                      igroup, when lgrps is False
!                      Sequential, when lgrps is True
!     n_skip        Number of line to skip at the beginning of the file
!     coil_group    Array of special type from module bsc. Each element in the array
!                      stores a coil collection
!     lgrps         Logical.
!                   True - Each separate coil is stored in its own coil group
!                   False - Placement into coil groups determined by integer igroup
!                      read from the coils. file
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(:)  :: index1(1)
      INTEGER            :: istat
      INTEGER     :: n_coilgroups, id_group
      INTEGER     :: i, igroup, inodes, nnod, icoil
      TYPE(bsc_coil)     :: coil_temp
      REAL(rprec), DIMENSION(3,nmaxnodes) :: xnod_in
      REAL(rprec)        :: currin, currin_first
      REAL(rprec)               :: this_rcirc
      REAL(rprec), DIMENSION(3) :: this_xcent, this_enhat
      CHARACTER          :: line*200, group_id*100
      CHARACTER(len=LEN(coil_group(1)%s_name)) :: s_name
      CHARACTER(len=LEN(coil_group(1)%l_name)) :: l_name
      LOGICAL            :: lparsed

!     index1         Default integer array - result from MINLOC
!     istat         Default integer, for I/O status
!     n_coilgroups   Number of external current coil groups.
!                    (In MAKEGRID this is nextcur. Read as n_extcur, not next_cur)
!     id_group       sequential id (from ngroup array) for coil group igroup
!     igroup         number of an externalcurrent coil group
!     inodes         counter for the number of nodes in a coil
!     nnod           number of nodes in a coil.
!     coil_temp      Special type from module bsc. Stores a single coil.
!     xnod_in        Real array to store coil nodes as they're being read in.
!     currin         current
!     currin_first   current of first node in a coil.
!     line           CHARACTER, for reading (and rereading) a line from the input file
!     group_id       CHARACTER, identifier of an external current coil group
!                    In MAKEGRID, this is variable group
!-----------------------------------------------
!     Find out how many coil groups there are
      n_coilgroups = SIZE(coil_group)

!     Create the coil collections
      DO i = 1,n_coilgroups
         WRITE(l_name,*) 'i = ',i
         CALL bsc_construct(coil_group(i),'boring id',l_name)
      END DO

!     Rewind the file
      REWIND (iou, iostat=istat)

!     Read in first n_skip lines
      DO i = 1,n_skip
         READ(iou,'(a)') line
      END DO
    
!     Initialize counters
      inodes = 0
      id_group = 0

!     Loop to read in the rest of the lines from the file
!     Last line of the file should be just 'end'.
      read_loop : DO
         READ(iou,'(a)',END = 100) line
         IF (line(1:3) .eq. 'end') EXIT
         inodes = inodes + 1

!     Reread the line
         READ(line,*,iostat=istat) xnod_in(1:3,inodes), currin
        
!     Save the current from the first node
         IF (inodes .eq. 1) currin_first = currin

!     Reread the line, assuming igroup and group_id are there.
!     For most of the file, they won't be there, and istat will be nonzero
!     But, when igroup and group_id are there, istat will be zero.
!     This is the indication of the END of a coil 
         READ(line,*,iostat=istat) xnod_in(1:3,inodes), currin,                &              
     &     igroup, group_id
         lparsed = (istat .eq. 0)
         IF (lparsed) THEN

!  JDH 2007-09-03. Modified below to make consistent with module coils_dot
!    behavior. In particular, interpretation of a single node coils is as
!    a circular coil. short name for the coil is unchanged from previous makegrid.
!     Find sequential group id no. for this igroup value
            IF (lgrps) THEN
               id_group = id_group + 1
            ELSE 
               index1 = MINLOC(ABS(igroup - ngroup(:)))
               id_group = index1(1)
               IF (igroup .ne. ngroup(id_group))                               & 
     &            STOP 'ID_GROUP != IGROUP in coils_dot_pass2'
            END IF

!     Define a short name for the coil
            icoil = coil_group(id_group) % ncoil + 1
            WRITE(s_name, '(a4,i5.5)') 'ID #', icoil

!  Various cases, depending on how many lines have been read in
            SELECT CASE (inodes)
            
            CASE (1) ! Circular coil
               this_xcent(1:3) = (/ zero, zero, xnod_in(3,1) /)
               this_enhat(1:3) = (/ zero, zero, 1.0_rprec /)
               this_rcirc = xnod_in(1,1)
               CALL bsc_construct(coil_temp,'fil_circ',s_name,'',              &                   
     &            currin_first, rcirc = this_rcirc,                            &
     &            xcent = this_xcent(1:3),enhat = this_enhat(1:3))

            CASE (2) ! Two point filament - treat like infinite straight coil
               nnod = inodes
               CALL bsc_construct(coil_temp,'fil_loop',s_name,'',              &                   
     &            currin_first,xnod_in(1:3,1:nnod))

            CASE DEFAULT ! Filamentary loop
!     Last point is supposed to be identical with the
!     first point of the coil. bsc_construct assumes that the
!     coil is not yet closed, so don't include the last point
               nnod = inodes - 1
               CALL bsc_construct(coil_temp,'fil_loop',s_name,'',              &                   
     &            currin_first,xnod_in(1:3,1:nnod))

            END SELECT
            
!     Append the coil to the appropriate coil group
            CALL bsc_append(coil_group(id_group),coil_temp)

!     Save the EXTERNAL coil group identifier group_id
            coil_group(id_group) % s_name = TRIM(group_id)
            WRITE (coil_group(id_group) % l_name,'(a7,i6.6)')
     1         ' IGROUP',igroup

!   Reset the number of nodes to zero
            inodes = 0

         END IF  ! of istat .eq. 0 IF 
      END DO read_loop
      
      RETURN

100   IF (.not. lparsed) THEN
         WRITE(6,*) 'Problems in read_coils_pass2'
         WRITE(6,*) 'EOF reached before END'
         WRITE(6,*) 'Make sure last line of file is "end"' 
      END IF
      
      END SUBROUTINE  read_coils_pass2

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------

      SUBROUTINE bfield (rp, phi, zp, br, bp, bz, ig)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, OPTIONAL :: ig
      REAL(rprec), INTENT(in)  :: rp, phi, zp
      REAL(rprec), INTENT(out) :: br, bp, bz
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: igroup
      REAL :: cosp, sinp
      REAL(rprec), DIMENSION(3) :: xpt, bvec
!-----------------------------------------------

!     Convert to cartesian coordinates
      cosp = COS(phi);    sinp = SIN(phi)
      xpt(1) = rp*cosp
      xpt(2) = rp*sinp
      xpt(3) = zp

      igroup = 1
      IF (PRESENT(ig)) igroup = ig

      CALL bsc_b (coil_group(igroup), xpt, bvec)

!     Convert back to cylindrical coordinates from cartesian vector components of B
      br = bvec(1)*cosp + bvec(2)*sinp
      bp =-bvec(1)*sinp + bvec(2)*cosp
      bz = bvec(3)

      END SUBROUTINE bfield

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------

      SUBROUTINE write_coils_file (extension)
      USE safe_open_mod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*) :: extension
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: cunit=30, ierr, ig, ncoils, n, nwire, iwire, nextcur
      REAL(rprec) :: current
      CHARACTER(len=LEN(coil_group(1)%s_name)) :: g_name
!-----------------------------------------------
!
!     WRITES (OR OVERWRITE) A COILS.EXT FILE BASED ON PRESENT COIL_GROUP DATA
!
      CALL safe_open(cunit, ierr, 'coils.' // TRIM(extension),
     1    'replace', 'formatted')
      IF (ierr .ne. 0)STOP 'Error opening coils-dot file in write_coils'

!
!     WRITE CANONICAL HEADER
!
      WRITE (cunit,100) nfp_bs
  100 FORMAT("periods ",i2,/,"begin filament",/,"mirror NUL")

!
!     WRITE (x, y, z, cur [,label]) INFO
!
      nextcur = SIZE(coil_group)
      DO ig = 1, nextcur
         ncoils = coil_group(ig) % ncoil
         g_name = coil_group(ig) % s_name
         DO n = 1, ncoils
            current = coil_group(ig) % coils(n) % current
            nwire = SIZE(coil_group(ig) % coils(n) % xnod, 2)
            IF (ANY(coil_group(ig) % coils(n) % xnod(:,1) .ne.
     1              coil_group(ig) % coils(n) % xnod(:,nwire))) 
     1      PRINT *, 'Coil did not close in WRITE_COILS_DOT for group ',
     2      ig,' COIL ',n
            DO iwire = 1, nwire-1
               WRITE(cunit,'(1p,4e22.14)') 
     1            coil_group(ig) % coils(n) % xnod(:,iwire), current
            END DO
            WRITE(cunit,'(1p,4e22.14,i4,1x,a)')
     1         coil_group(ig) % coils(n) % xnod(:,nwire), zero, ig,
     2         TRIM(g_name)
         END DO
      END DO

      WRITE (cunit, '(a3)') "end"

      CLOSE (cunit)

      END SUBROUTINE write_coils_file

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------

      SUBROUTINE afield (rp, phi, zp, ar, ap, az, ig)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, OPTIONAL :: ig
      REAL(rprec), INTENT(in)  :: rp, phi, zp
      REAL(rprec), INTENT(out) :: ar, ap, az
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: igroup
      REAL :: cosp, sinp
      REAL(rprec), DIMENSION(3) :: xpt, avec
!-----------------------------------------------

!     Convert to cartesian coordinates
      cosp = COS(phi);    sinp = SIN(phi)
      xpt(1) = rp*cosp
      xpt(2) = rp*sinp
      xpt(3) = zp

      igroup = 1
      IF (PRESENT(ig)) igroup = ig

      CALL bsc_a (coil_group(igroup), xpt, avec)

!     Convert back to cylindrical coordinates from cartesian vector components of B
      ar = avec(1)*cosp + avec(2)*sinp
      ap =-avec(1)*sinp + avec(2)*cosp
      az = avec(3)

      END SUBROUTINE afield

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------

      SUBROUTINE cleanup_biotsavart
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i
!-----------------------------------------------
      
      IF (ALLOCATED(coil_group)) THEN
         DO i = 1, SIZE(coil_group)
            CALL bsc_destroy(coil_group(i))
         END DO
         DEALLOCATE(coil_group)
      END IF

      IF (ASSOCIATED(single_coil)) THEN
         CALL bsc_destroy(single_coil)
         DEALLOCATE(single_coil)
      END IF

      END SUBROUTINE cleanup_biotsavart

      END MODULE biotsavart

!MODIFICATION HISTORY
!03.16.04 (SPH) Added l_name to coil_group (in read_pass2), based on original igroup in coils-dot file
!10.09.08 (SPH) REPLACE INTEGER(iprec) with INTEGER
!11.24.10 (SPH) added lparsed logical to avoid eof error message if coils file does not terminate with "end"
! 2011-04-01 (JDH) Added start_string to parse_coils_file. Also, miscellaneous cleanup.
! 2014-08-01 (SAL) Added afield routine (copy of bfield)
