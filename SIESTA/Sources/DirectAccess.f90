      MODULE DirectAccess
      USE safe_open_mod
      USE stel_kinds, ONLY: dp

!     DIRECT ACCESS FILE HANDLING ROUTINES
      INTEGER, PARAMETER :: CreateNew=0, OpenExisting=1, Scratch=2
      INTEGER :: rec_length, data_size, block_size, num_rows
      INTEGER :: iunit_da, blocks_per_row, recs_per_block
      INTEGER :: irec_pos, byte_size_rec, byte_size_dp
      CHARACTER(LEN=256) :: filename

      CONTAINS 

      SUBROUTINE OpenDAFile(datasize, blksize, blocksperrow,            &
     &                      filename_in, iunit, iflag)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: datasize, blksize, blocksperrow,           &
     &                       iflag
      INTEGER, INTENT(inout) :: iunit 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER*(*), INTENT(in) :: filename_in
      INTEGER :: ierr
      REAL(dp) :: dummy
      CHARACTER(LEN=10) :: Status
!-----------------------------------------------
   
!GET EFFECTIVE "byte_size" OF UNFORMATTED dp VARIABLE (machine dependent!)
!ON Compaq, rec length size is in units of 4 byte chunks, so byte_size_rec = 4
!ON lf95, byte_size_rec = 1
      INQUIRE(iolength=byte_size_rec) dummy
      byte_size_dp = KIND(dummy)
   
      data_size = datasize
      rec_length = byte_size_rec*datasize
      block_size = byte_size_dp*blksize
      blocks_per_row = blocksperrow
      filename = filename_in
!skip this distance to next block, if datasize different from blksize
      recs_per_block  = MAX(1,blksize/datasize)
      irec_pos = 0

!  create disk file for doing direct access i/o.

      IF (iflag .eq. CreateNew) THEN
         Status = "replace"
      ELSE IF (iflag .eq. OpenExisting) THEN
         Status = "old"
      ELSE
         Status = "scratch"
      END IF


      CALL safe_open(iunit, ierr, filename, Status, 'unformatted',   &
     &     rec_length, 'DIRECT')

      iunit_da = iunit

      IF (ierr .ne. 0) THEN
         WRITE (6, '(a7,i4)') 'Status code: ', Status, ' Error stat: ', ierr
         STOP 'Error creating Direct Access file!'
      END IF

      END SUBROUTINE OpenDAFile

      
      SUBROUTINE ChangeDAFileParams(datasize, blksize,                  &
     &                              blocksperrow, new_filename, nrows)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: datasize, blksize, blocksperrow, nrows
      INTEGER :: ierr, new_data_size, new_block_size, new_rec_length,   &
     &           new_blocks_per_row, new_recs_per_block, inew_da
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, j, k, boffset, recloc, nsplit, isplit, new_row_offset
      CHARACTER*(*)  :: new_filename
      REAL(dp), ALLOCATABLE :: DataItem(:)
!-----------------------------------------------
  
!store new parameters
      new_data_size = datasize
      new_rec_length = byte_size_rec*datasize
      new_block_size = byte_size_dp*blksize
      new_blocks_per_row = blocksperrow
!skip this distance to next block, if rec_length different from block_size
      new_recs_per_block  = MAX(1,blksize/datasize)
      inew_da      = iunit_da+1

      IF ((rec_length.eq.new_rec_length) .and.                          &
     &    (block_size.eq.new_block_size)) RETURN


!  create disk file for doing direct access i/o.
      ierr = INDEX(filename,'.',back=.true.)

      filename = new_filename

      CALL safe_open(inew_da, ierr, filename, 'replace', 'unformatted', &
     &     new_rec_length, 'DIRECT')
      IF (ierr .ne. 0) STOP 'Error opening existing Direct Access file!'

      ALLOCATE (DataItem(new_data_size), stat=ierr)
      IF (ierr .ne. 0) STOP 'Allocation error in ChangeDAFileParams'

!Load up block by reading a column at a time
      IF (irec_pos .eq. 0) THEN
         DO i = 1, nrows
            DO j = 1, blocks_per_row
               boffset = 1
               DO k = 1, recs_per_block
                  CALL ReadDAItem1(DataItem(boffset), i, j, k)
                  boffset = boffset + recs_per_block
               END DO
               !write full block in new da file
               recloc = 1 + new_recs_per_block*((j-1) + new_blocks_per_row*(i-1))
               WRITE (inew_da, rec=recloc, iostat=ierr) DataItem
            END DO
         END DO

      ELSE 
         num_rows = nrows

!        BLOCK POSITION CORRELATED WITH :: blpls=1 (U), bldia=2 (D), blmin=3 (L)
!        Data stored as follows for each COLUMN (m,n,ntype):
!        U0,D1,L2 U3,D4,L5  .... U(3*n1)  ,D(3*n1+1),L(3*n1+2)
!        U1,D2,L3 U4,D5,L6       U(3*n2+1),D(3*n2+2),L(3*n2+3)
!        U2,D3,L4 U5,D6,L7       U(3*n3+2),D(3*n3+3),L(3*n3+4)
!
!        Therefore we must split the i=1,nrows loop into 3
         isplit = 0
         DO nsplit = 1, 3
         DO i = nsplit, nrows, 3
            isplit = isplit+1         !sequential index of records in original file
            DO j = 1, blocks_per_row    
!        j=1=>U  2=>D  3=>L block
               boffset = 1
               new_row_offset = 2-j+i
               IF (new_row_offset.lt.1 .or. new_row_offset.gt.nrows) CYCLE
               DO k = 1, recs_per_block
!        Read one column (k) at a time
                  CALL ReadDAItem_SEQ(DataItem(boffset), isplit, j, k)
                  boffset = boffset + recs_per_block
               END DO
               !write full block to new da file
               recloc = 1 + new_recs_per_block*((j-1) + new_blocks_per_row*(new_row_offset-1))
               WRITE (inew_da, rec=recloc, iostat=ierr) DataItem
            END DO
         END DO
         END DO

         IF (isplit .ne. nrows) STOP 'isplit != nrows'

      END IF

!     Close old scratch file when finished writing out new file
      CLOSE (iunit_da, status='delete', iostat=ierr)
 
      iunit_da = inew_da
      data_size = new_data_size
      rec_length = new_rec_length
      block_size = new_block_size
      blocks_per_row = new_blocks_per_row
      recs_per_block = new_recs_per_block
      irec_pos = 0

      DEALLOCATE (DataItem)

      END SUBROUTINE ChangeDAFileParams

      
      SUBROUTINE CloseDAFile
      IMPLICIT NONE
      INTEGER :: istat

      IF (iunit_da .gt. 6) CLOSE (iunit_da, IOSTAT=istat)
      iunit_da = 0

      END SUBROUTINE CloseDAFile

      
      SUBROUTINE RemoveDAFile
      IMPLICIT NONE
      INTEGER :: istat

      CALL CloseDAFile

      iunit_da=100
      OPEN(FILE=filename, UNIT=iunit_da, IOSTAT=istat)
      IF (istat .eq. 0) CLOSE(iunit_da,STATUS='delete',IOSTAT=istat)

      END SUBROUTINE RemoveDAFile

      
      SUBROUTINE WriteDAItem_RA(DataItem, BlockRowIndex, ColIndex, IndexInBlock)
      IMPLICIT NONE
      REAL(dp), INTENT(in) :: DataItem(data_size)
      INTEGER, INTENT(in)  :: BlockRowIndex, ColIndex, IndexInBlock
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: recloc, ierr
      INTEGER :: StartIndex
!-----------------------------------------------

      IF (ColIndex > blocks_per_row) STOP 'ColIndex > Block_Per_Row in WriteDAItem'
      IF (IndexInBlock > recs_per_block)  STOP 'IndexInBloc > skip_size in WriteDAItem'
         

      StartIndex = IndexInBlock
      IF (recs_per_block .eq. 1) StartIndex = 1
      recloc = StartIndex + recs_per_block*((ColIndex-1) + blocks_per_row*(BlockRowIndex-1))

      WRITE (iunit_da, rec=recloc, iostat=ierr) DataItem
      IF (ierr .ne. 0) THEN
         WRITE (6,*) 'Ierr = ', ierr, ' in WriteDAItem'
         STOP
      END IF

      END SUBROUTINE WriteDAItem_RA

      
      SUBROUTINE WriteDAItem_SEQ(DataItem)
      IMPLICIT NONE
      REAL(dp), INTENT(in) :: DataItem(data_size)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ierr
!-----------------------------------------------

      !Perform "sequential" write of direct access elements (buffering helps)

      irec_pos = irec_pos+1

      WRITE (iunit_da, rec=irec_pos, iostat=ierr) DataItem
      IF (ierr .ne. 0) THEN
         WRITE (6,*) 'Ierr = ', ierr, ' in WriteDAItem'
         STOP
      END IF

      END SUBROUTINE WriteDAItem_SEQ


      SUBROUTINE ReadDAItem1(DataItem, BlockRowIndex, ColIndex, StartIndex)
      IMPLICIT NONE
      REAL(dp), INTENT(out) :: DataItem(data_size)
      INTEGER, INTENT(in)  :: BlockRowIndex, ColIndex, StartIndex
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: recloc, ierr
!      INTEGER :: StartIndex=1
!-----------------------------------------------

      recloc = StartIndex + recs_per_block*((ColIndex-1) + blocks_per_row*(BlockRowIndex-1))
      READ (iunit_da, rec=recloc, iostat=ierr) DataItem
      IF (ierr .ne. 0) THEN
         WRITE (6,*) 'Ierr = ', ierr, ' in ReadDAItem'
         STOP
      END IF

      END SUBROUTINE ReadDAItem1


      SUBROUTINE ReadDAItem2(DataItem, BlockRowIndex, ColIndex)
      IMPLICIT NONE
      REAL(dp), INTENT(out) :: DataItem(data_size)
      INTEGER, INTENT(in)  :: BlockRowIndex, ColIndex
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: recloc, ierr
!      INTEGER :: StartIndex=1
!-----------------------------------------------

      recloc = 1 + recs_per_block*((ColIndex-1) + blocks_per_row*(BlockRowIndex-1))
      READ (iunit_da, rec=recloc, iostat=ierr) DataItem
      IF (ierr .ne. 0) THEN
         WRITE (6,*) 'Ierr = ', ierr, ' in ReadDAItem'
         STOP
      END IF

      END SUBROUTINE ReadDAItem2


      SUBROUTINE ReadDAItem_SEQ(DataItem, BlockRowIndex, BlockTypeIndex, ColIndex)
      IMPLICIT NONE
      REAL(dp), INTENT(out) :: DataItem(data_size)
      INTEGER, INTENT(in)  :: BlockRowIndex, ColIndex, BlockTypeIndex
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: recloc, ierr
!-----------------------------------------------

      recloc = BlockTypeIndex + blocks_per_row*(BlockRowIndex-1+num_rows*(ColIndex-1))
      READ (iunit_da, rec=recloc, iostat=ierr) DataItem
      IF (ierr .ne. 0) THEN
         WRITE (6,*) 'Ierr = ', ierr, ' in ReadDAItem'
         STOP
      END IF

      END SUBROUTINE ReadDAItem_SEQ


      END MODULE DirectAccess
