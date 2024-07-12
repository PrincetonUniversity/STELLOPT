!-----------------------------------------------------------------------
!     Subroutine:    write_pies_coil
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/20/2012
!     Description:   This subroutine reads the coils file and creates a
!                    PIES coil_data file.  The coil currents are
!                    replaced by the values found in the EXTCUR array.
!-----------------------------------------------------------------------
      SUBROUTINE write_pies_coil
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_runtime
!-----------------------------------------------------------------------
!     Local Variables
!          iunit       File Unit Number
!          ierr        Error flag
!          mn          Dummy index over modes
!          im          Dummy index over modes
!          in          Dummy index over modes
!          uv          Dummy index over real-space
!          ik          Dummy index over radial surfaces
!          vsurf       VMEC surface in PIES background coordinates
!          rho_in      VMEC Rho
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: maxseg, maxcoil, maxgroup
      INTEGER :: iunit, ncoil, nseg, ngroup,c,i,j,ier, cg
      INTEGER, ALLOCATABLE :: dex1(:), dex2(:), cgroup(:),&
                              ncoilgrp(:)
      REAL(rprec) :: dumx, dumy, dumz, dumc
      REAL(rprec), ALLOCATABLE :: x(:), y(:), z(:), curre(:), &
                                  coilcur(:), groupcur(:)
      CHARACTER*32, ALLOCATABLE :: group(:)
      CHARACTER*256 :: top(3)
      
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      iunit = 321
      ncoil = 0
      nseg  = 0
      maxgroup = nextcur
      WHERE (ABS(extcur) < 1) extcur = 0.0
      OPEN(UNIT=iunit,FILE=TRIM(coils_file))
      READ(iunit,'(A)') (top(i), i=1,3)                               ! Read Header
      ! Count up the file
      DO WHILE(.TRUE.)
         dumc = 1
         DO WHILE(dumc .ne. 0)
            READ(iunit,*,ERR=22,END=22) dumx, dumy, dumz, dumc
            nseg = nseg + 1
         END DO
         ncoil = ncoil + 1
      END DO
 22   REWIND(iunit)
      maxseg = nseg
      maxcoil = ncoil
      ALLOCATE(x(maxseg),y(maxseg),z(maxseg),curre(maxseg),STAT=ier)
      ALLOCATE(dex1(maxcoil),dex2(maxcoil),coilcur(maxcoil),&
               cgroup(maxcoil),group(maxcoil), STAT = ier)
      ALLOCATE(ncoilgrp(maxgroup),groupcur(maxgroup), STAT = ier)
      ncoil = 0
      nseg = 0
      READ(iunit,'(a32)') (top(i), i=1,3)                               ! Read Header 
      DO WHILE(.TRUE.)
         ncoil = ncoil + 1
         dex1(ncoil) = nseg + 1
         i=0
         c=1
         DO WHILE (c .ne. 0)
            i = i + 1
            nseg = nseg + 1
            READ(iunit,*,ERR=98,END=98) x(nseg),y(nseg),z(nseg),curre(nseg)
            c=curre(nseg)
            IF (i .eq. 1) coilcur(ncoil) = c
         END DO
         ! Now read the group info
         BACKSPACE(iunit)
         READ(iunit,*) dumx, dumy, dumz, dumc, cgroup(ncoil), group(ncoil)
         ncoilgrp(cgroup(ncoil))=ncoilgrp(cgroup(ncoil)) + 1
         ! Normalize the current to the first segment in the group
         IF (ncoilgrp(cgroup(ncoil)) .eq. 1) THEN
            groupcur(cgroup(ncoil)) = coilcur(ncoil)
            coilcur(ncoil) = 1.
         ELSE
            coilcur(ncoil) = coilcur(ncoil)/groupcur(cgroup(ncoil))
         END IF
         ! End each coil
         dex2(ncoil) = nseg
         !WRITE(6,'(A,I7,A,I7,A,I7,A)')'  Coil Number: ',ncoil,&
         !          ' dex=[',dex1(ncoil),',',dex2(ncoil),']'
      END DO
 98   ncoil = ncoil - 1
      nseg = nseg - 1
      ngroup = MAXVAL(cgroup)
      CLOSE(iunit)
      ! Output some stuff
      write(6,*) '-----Coil File Parameters-----'
      write(6,'(A,A)')       '        file: ',TRIM(coils_file)
      write(6,'(A,I8)') '  datapoints: ',nseg
      write(6,'(A,I8)') '       coils: ',ncoil
      write(6,'(A,I8)') '      groups: ',ngroup
      DO i = 1, ngroup
         WRITE(6,'(A,I3,A,E21.14)')'  EXTCUR(',i,') = ',extcur(i)
      END DO
      ! Now write out the coil_data file
      OPEN(UNIT=iunit,FILE=TRIM('coil_data'))
      DO i = 1, ncoil
         cg = cgroup(i)
         DO j = dex1(i), dex2(i)-1
            WRITE(iunit,"(1X,4(1PE21.14,1X))") x(j),y(j),z(j),extcur(cg)
         END DO
         j = dex2(i)
         WRITE(iunit,"(1X,4(1PE21.14,1X))") x(j),y(j),z(j),0.0
!         WRITE(iunit,"(1X,4(1PE21.14,1X),i3,1X,A32)") x(j),y(j),z(j),&
!                     0.0,cgroup(i),group(i)
      END DO
      CLOSE(iunit)    
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE write_pies_coil
