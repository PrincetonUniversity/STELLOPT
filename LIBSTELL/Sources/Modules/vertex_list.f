      MODULE vertex_list
      USE stel_kinds
      IMPLICIT NONE
      INTEGER, PARAMETER           :: nvertex = 4

      PRIVATE                       :: FIND_VERTEX

      TYPE POINT
         INTEGER                   :: xpt, ypt, phi_plane
      END TYPE POINT

      TYPE MAGFIELD
         REAL(rprec)               :: br, bp, bz
      END TYPE MAGFIELD
      
      
      TYPE PBFIELD                                     
         TYPE (MAGFIELD), POINTER  :: pbfld            !POINTER to BFIELD
      END TYPE PBFIELD
      
      
      TYPE VERTEX
         TYPE (POINT)            :: p1                 !integer-INDEXed position of vertex (lower-left point)
         TYPE (VERTEX), POINTER  :: pprev              !POINTER to previously ALLOCATED  vertex (next in list)
         TYPE (PBFIELD)          :: bvector(nvertex)   !pointers to the bfield vector at nvertex nodes of square
         LOGICAL                 :: lalloc(nvertex)    !true IF bvector(k) is ALLOCATED
      END TYPE VERTEX

      TYPE PVERTEX                                     !POINTER to VERTEX
         TYPE (VERTEX), POINTER  :: pverx
      END TYPE PVERTEX
      
      INTEGER                             :: nlist        !number of lists
      REAL(rprec)                         :: hR, hZ, hP   !scale factors for converting xpt to R, etc.
      INTEGER, ALLOCATABLE                :: ibs_cnt(:)   !counts number of biot-savart calls made
      TYPE (PVERTEX), ALLOCATABLE, TARGET :: tail(:)      !array of root vertices for separate lists      
      TYPE (VERTEX), POINTER              :: root         !last vertex most recently added to list
      TYPE (VERTEX), POINTER              :: new_root, prev_vertex    !stores previously computed nodes
      LOGICAL                             :: lsearch(nvertex)   !true IF bfield at ith node already computed

      CONTAINS
      
      SUBROUTINE initialize_list (number_lists, hR_in, hZ_in, hP_in)
      IMPLICIT NONE
      INTEGER, INTENT(in)        :: number_lists
      REAL(rprec), INTENT(in)    :: hR_in, hZ_in, hP_in
      INTEGER                    :: ilist

      nlist = number_lists
      IF (nlist .le. 0) STOP ' Must specify nlist > 0'

      ALLOCATE (tail(nlist), ibs_cnt(nlist))
            
      DO ilist = 1, nlist 
         NULLIFY (tail(ilist)%pverx)
      END DO

      ibs_cnt(1:nlist) = 0

      hR = hR_in;   hZ = hZ_in;    hP = hP_in

      END SUBROUTINE initialize_list


      SUBROUTINE add_vertex (rin, zin, phi_pt, list_index, bfield)
      IMPLICIT NONE
      TYPE (VERTEX), POINTER     :: new_vertex, swap
      TYPE (MAGFIELD), POINTER   :: MagFld
      INTEGER, INTENT(in)        :: phi_pt, list_index
      REAL(rprec), INTENT(in)    :: rin, zin
      REAL(rprec)                :: Rreal, Zreal, PhiReal
      INTEGER                    :: rpt, zpt, xpt, ypt
      INTEGER                    :: icount
      LOGICAL                    :: lexists
      EXTERNAL bfield                              !Could be magnetic field or something ELSE...
!
!     adds new_vertex to tail of list (if it does not exist already), assigning it to root
!     If it exists already, moves vertex to tail of list, assigning it to root
!     thus, this vertex can be identified with "root" for use by external programs

      IF (list_index .gt. nlist) STOP 'list_index > nlist in add_vertex'
      IF (list_index .le. 0) STOP 'list_index <= 0 in add_vertex'

!     Get lower left INDEX of nearest vertex. Z can be positive or negative so USE FLOOR    
      rpt    = floor(rin/hR)
      zpt    = floor(zin/hZ)

      root => tail(list_index)%pverx
      
      ALLOCATE (new_vertex)
      IF (ASSOCIATED(root)) THEN
         new_vertex%pprev => root                  !point to previous node unless this is the first in list
      ELSE
         NULLIFY(new_vertex%pprev)
      END IF
      
      new_vertex%p1 = POINT(rpt, zpt, phi_pt)      !store INTEGER INDEX of lower-left vertex on grid
      
!     Assign bfield values to this vertex from exist vertex values
      CALL find_vertex (new_vertex, lexists)

      IF (lexists) THEN
!        new_node already existed!
!        move it to root position. and assign root to new_root 
!        prev_vertex is the preceding vertex in list
         DEALLOCATE (new_vertex)
         IF (ASSOCIATED(prev_vertex)) THEN
            swap => root
            root => new_root
            root%pprev => swap
         END IF
!        PRINT *,' NODE ALREADY EXISTED IN LIST #', list_index

      ELSE
         root => new_vertex
         
         DO icount = 1, nvertex 
            IF (.not.lsearch(icount)) THEN
!              Node is one of four corners of square, depending on value of icount (1-4)
               ALLOCATE (new_vertex%bvector(icount)%pbfld)
               new_vertex%lalloc(icount) = .true.
               MagFld => new_vertex%bvector(icount)%pbfld
               !Simulate bio-savart CALL; must add toroidal plane in arg list ...
               xpt = rpt;   ypt = zpt
               IF (icount.eq.2 .or. icount.eq.4) xpt = xpt+1
               IF (icount.ge.3) ypt = ypt+1
               Rreal = xpt*hR;    Zreal = ypt*hZ
               PhiReal = (phi_pt-1)*hP
               CALL Bfield (Rreal, PhiReal, Zreal, 
     1                      MagFld%br, MagFld%bp, MagFld%bz)
               ibs_cnt(list_index) = ibs_cnt(list_index) + 1
            ELSE
               new_vertex%lalloc(icount) = .false.
            END IF
         END DO
      END IF
      
      tail(list_index)%pverx => root

      END SUBROUTINE add_vertex


      SUBROUTINE find_vertex (new_node, lexists)
      IMPLICIT NONE
      TYPE (VERTEX), POINTER :: new_node
      LOGICAL, INTENT(out)   :: lexists
      TYPE (VERTEX), POINTER :: stepit, previous
      INTEGER                :: xpt(nvertex), ypt(nvertex), xvert, 
     1                          yvert, xstep, ystep, icount, ivertex
      
!     Scans vertex (assumed in a given phi plane) ASSOCIATED with new_node and assigns bfield values
!     IF they already exist in list. If it finds the bfield ASSOCIATED with the new_node
!     or ANY of the other 3 corners, the global pointers new_root, prev_vertex become ASSOCIATED.

!     Find nodes associate with the vertex new_node WHERE Bfield is required
      xpt(1) = new_node%p1%xpt;    ypt(1) = new_node%p1%ypt
      xpt(2) = xpt(1)+1;           ypt(2) = ypt(1)
      xpt(3) = xpt(1);             ypt(3) = ypt(1)+1
      xpt(4) = xpt(2);             ypt(4) = ypt(3)

      stepit => root;      
      NULLIFY(previous); NULLIFY(new_root); NULLIFY(prev_vertex)
      lexists = .false.
      
      lsearch(1:nvertex) = .false.

      icount = 0
      DO WHILE (icount.lt.nvertex .and. ASSOCIATED(stepit))
         xstep = stepit%p1%xpt;     ystep = stepit%p1%ypt
!
!          Logic to catch case whether vertex already exists
!          This may miss the case if all nodes are filled first by adjacent vertices, but does not matter
!
         IF (xstep.eq.xpt(1) .and.  ystep.eq.ypt(1)) THEN
            new_root => stepit
            IF (ASSOCIATED(previous)) previous%pprev => stepit%pprev
            prev_vertex => previous 
            lexists = .true.                  
            EXIT
         END IF
         DO ivertex = 1, nvertex 
            xvert = xpt(ivertex); yvert = ypt(ivertex)
   
            IF ((xstep.gt.xvert) .or. (ystep.gt.yvert) .or.
     1           (xstep.lt.xvert-1) .or. (ystep.lt.yvert-1) .or. 
     2           lsearch(ivertex) ) CYCLE

            SELECT CASE (ivertex)
            CASE (1)                           !!Lower left vertex
               IF (xstep.eq.xvert) THEN
                  new_node%bvector(1) = stepit%bvector(3)
               ELSE 
                  IF (ystep .eq. yvert) THEN
                     new_node%bvector(1) = stepit%bvector(2)
                  ELSE                         !!ystep < yvert
                     new_node%bvector(1) = stepit%bvector(4)
                  END IF
               END IF
                  
            CASE (2)                           !!Lower right vertex
               IF (xstep.lt.xvert) THEN
                  new_node%bvector(2) = stepit%bvector(4)
               ELSE                            !!xstep == xvert
                  IF (ystep .eq. yvert) THEN
                     new_node%bvector(2) = stepit%bvector(1)
                  ELSE                         !!ystep < yvert
                     new_node%bvector(2) = stepit%bvector(3)
                  END IF
               END IF
            
            CASE (3)                           !!Upper left vertex
               IF (xstep .eq. xvert) THEN
                  new_node%bvector(3) = stepit%bvector(1)
               ELSE                            !!xstep < xvert
                  IF (ystep .eq. yvert) THEN
                     new_node%bvector(3) = stepit%bvector(2)
                  ELSE                         !!ystep < yvert
                     new_node%bvector(3) = stepit%bvector(4)
                  END IF
               END IF

            CASE (4)                           !!Upper right vertex
               IF (xstep .lt. xvert) THEN
                  new_node%bvector(4) = stepit%bvector(2)
               ELSE                            !!xstep == xvert
                  IF (ystep .eq. yvert) THEN
                     new_node%bvector(4) = stepit%bvector(1)
                  ELSE                         !!ystep < yvert
                     new_node%bvector(4) = stepit%bvector(3)
                  END IF
               END IF
            
            END SELECT

            icount = icount + 1
            IF (icount .ge. nvertex) EXIT
            lsearch(ivertex) = .true.

         END DO

         previous => stepit
         stepit => stepit%pprev

      END DO

      END SUBROUTINE find_vertex


      SUBROUTINE b_at_point (rin, zin, br, bp, bz)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(in)           :: rin, zin
      REAL(rprec), INTENT(out)          :: br, bp, bz
      REAL(rprec), DIMENSION(nvertex)   :: brvac, bpvac, bzvac
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      TYPE (MAGFIELD), POINTER   :: bvac
      REAL(rprec), PARAMETER     :: zero = 0, one = 1
      INTEGER                    :: ir, jz, i
      REAL(rprec)                :: ri, zj, dri, dzj, fint(nvertex)
C-----------------------------------------------
!
!      THIS ROUTINE FINDS THE CYLINDRICAL COMPONENTS OF THE EXTERNAL
!      MAGNETIC FIELD AT A FIXED PHI PLANE AT THE POINT (rin, zin) BY 2-D INTERPOLATION
!      
!      THIS CALL MUST BE MADE JUST AFTER "ADD_VERTEX" CALL, SO THAT
!      "ROOT" CONTAINS THE CORRECT LOWER-LEFT VERTEX IN THE PROPER PHI-PLANE       
!

!
!       DETERMINE R, Z VERTEX CORNER INDICES (IR,IZ)
!
         ir = root%p1%xpt
         jz = root%p1%ypt

!
!       COMPUTE R , Z AND DR , DZ WITH RESPECT TO LOWER-LEFT VERTEX
!

         ri = ir*hR
         zj = jz*hZ
         
         dri = (rin - ri)/hR
         dzj = (zin - zj)/hZ
         
!
!      RETRIEVE B FIELD AT 4 VERTICES OF ROOT GRID ELEMENT
!

         DO i = 1, nvertex
            bvac => root%bvector(i)%pbfld
            brvac(i) = bvac%br
            bpvac(i) = bvac%bp
            bzvac(i) = bvac%bz
         END DO
         
!
!       COMPUTE INTERPOLATED B FIELD
!
         fint(4) = dri*dzj
         fint(3) = dzj - fint(4)
         fint(2) = dri - fint(4)
         fint(1) = 1 + fint(4) - (dri + dzj)

         IF (dri.lt.zero .or. dri.gt.one .or.
     1        dzj.lt.zero .or. dzj.gt.one) THEN
            WRITE (6, '(2(a,1pe12.4))')
     1       'Possible interp. error in b_from_coil: dri = ', 
     2       dri,' dzj = ', dzj
         END IF
         
         br = SUM(fint*brvac)
         bp = SUM(fint*bpvac)
         bz = SUM(fint*bzvac)
         
      END SUBROUTINE b_at_point

      SUBROUTINE cleanup_list
      IMPLICIT NONE
      TYPE (VERTEX), POINTER  :: stepit
      INTEGER                :: icount, ilist

      DO ilist = 1, nlist
         root => tail(ilist)%pverx
         DO WHILE (ASSOCIATED(root))
            stepit => root%pprev
            DO icount = 1, nvertex
               IF (root%lalloc(icount)) 
     1             DEALLOCATE (root%bvector(icount)%pbfld)
            END DO
            DEALLOCATE (root)
            root => stepit
         END DO
      END DO

      DEALLOCATE (tail, ibs_cnt)
      
      END SUBROUTINE cleanup_list

      END MODULE vertex_list
