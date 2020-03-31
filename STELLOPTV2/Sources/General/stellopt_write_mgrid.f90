!-----------------------------------------------------------------------
!     Subroutine:    stellopt_write_mgrid(MPI_COMM)
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/15/2016
!     Description:   This subroutine reads a coils file and generates
!                    the appropriate vacuum grid file.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_write_mgrid(COMM_LOCAL,proc_string,lscreen)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE vmec_input,  ONLY: nigroup, extcur, lasym, nfp, nzeta,&
                             mgrid_file_vmec => mgrid_file
      USE biotsavart
      USE write_mgrid, ONLY: nr=>ir, nz=>jz, nphi => kp, rmin, rmax, &
                            zmin, zmax, lstell_sym, br,bp,bz, mgrid_mode, &
                            mgrid_file, nextcur
      USE mgrid_mod, ONLY: vn_nextcur, vn_mgmode, vn_ir, &
          vn_jz, vn_kp, vn_nfp, vn_rmin, vn_rmax, &
          vn_zmin, vn_zmax, vn_coilgrp, vn_coilcur, & 
          vn_br0, vn_bz0, vn_bp0
      USE ezcdf
      USE mpi_params
      USE mpi_inc
      USE stellopt_runtime, ONLY: pi2,NETCDF_OPEN_ERR,handle_err
      IMPLICIT NONE
!-----------------------------------------------------------------------
!     Input Variables
!        COMM_LOCAL   Local MPI communicator
!----------------------------------------------------------------------
      INTEGER :: COMM_LOCAL
      CHARACTER(*) :: proc_string
      LOGICAL, INTENT(inout)        :: lscreen
!-----------------------------------------------------------------------
!     Local Parameters
!----------------------------------------------------------------------
      CHARACTER(LEN=*), PARAMETER ::&
         coildim(2) = (/'stringsize          ',&
                        'external_coil_groups'/),&
         groupdim(1)= (/'external_coils'/),&
         cylcoord(3)= (/'rad','zee','phi'/)
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER(KIND=BYTE_8) :: chunk
      INTEGER :: numprocs_local, mystart, myend
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
      LOGICAL :: lfile_found
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: lfile_found_global
      INTEGER :: ig, i, j, k, s, ngrid
      REAL(rprec) :: r, phi, z, dr, dphi, dz, current, current_first
      CHARACTER(LEN=100), ALLOCATABLE, DIMENSION(:) :: vn_br, vn_bp, vn_bz
      CHARACTER(256) :: coil_string
      CHARACTER(LEN=100) :: temp
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      numprocs_local=1
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_COMM_SIZE( COMM_LOCAL, numprocs_local, ierr_mpi )
      CALL MPI_COMM_RANK( COMM_LOCAL, myworkid, ierr_mpi )
      ALLOCATE(lfile_found_global(numprocs_local))
!DEC$ ENDIF

      IF (lscreen) WRITE(6,*) '---------------------------  WRITING MGRID  ------------------------'
      ! Read Coil
      coil_string = 'coils.'//TRIM(proc_string)
!DEC$ IF DEFINED (MPI_OPT)
      lfile_found = .false.
      DO 
         INQUIRE(FILE=TRIM(coil_string),EXIST=lfile_found)
         CALL MPI_ALLGATHER(lfile_found,1,MPI_LOGICAL,lfile_found_global,1,MPI_LOGICAL,COMM_LOCAL,ierr_mpi)
         IF (ALL(lfile_found_global)) EXIT
      END DO
      DEALLOCATE(lfile_found_global)
!DEC$ ENDIF
      CALL parse_coils_file(TRIM(coil_string))

      ! Parse EXTCUR Array
      nextcur = SIZE(coil_group) !SAL
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BCAST(extcur(1:nextcur),nextcur,MPI_REAL8, master, COMM_LOCAL, ierr_mpi)
!DEC$ ENDIF
      DO i = 1, nextcur
         DO j = 1, coil_group(i) % ncoil
            current = coil_group(i) % coils(j) % current
            IF (j .eq. 1) current_first = current
            SELECT CASE(mgrid_mode)
               CASE('S')
                  IF (current_first .ne. zero) coil_group(i) % coils(j) % current = (current/current_first)
!               CASE DEFAULT
!                  coil_group(i) % coils(j) % current = current*extcur(i)
            END SELECT
         END DO
      END DO

      ! Print out some info
      IF (lscreen) THEN
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   R   = [',rmin,',',rmax,'];  NR:   ',nr
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   PHI = [',0.0,',',pi2/nfp,'];  NPHI: ',nphi
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   Z   = [',zmin,',',zmax,'];  NZ:   ',nz
         CALL FLUSH(6)
      END IF

      ! Setup Grid Helpers
      dr = (rmax - rmin)/REAL(nr-1)
      dz = (zmax - zmin)/REAL(nz-1)
      dphi = pi2/REAL(nphi)/REAL(nfp)
      ALLOCATE(br(nr,nz,nphi),bz(nr,nz,nphi),bp(nr,nz,nphi))

      ! Break up the Work
      chunk = FLOOR(REAL(nr*nphi*nz) / REAL(numprocs_local))
      mystart = myworkid*chunk + 1
      myend = mystart + chunk - 1

      ! Setup Chunking
!DEC$ IF DEFINED (MPI_OPT)
      chunk = FLOOR(REAL(nr*nphi*nz) / REAL(numprocs_local))
      mystart = myworkid*chunk + 1
      myend = mystart + chunk - 1
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
      ALLOCATE(mnum(numprocs_local), moffsets(numprocs_local))
      CALL MPI_ALLGATHER(chunk,1,MPI_INTEGER,mnum,1,MPI_INTEGER,COMM_LOCAL,ierr_mpi)
      CALL MPI_ALLGATHER(mystart,1,MPI_INTEGER,moffsets,1,MPI_INTEGER,COMM_LOCAL,ierr_mpi)
      i = 1
      DO
         IF ((moffsets(numprocs_local)+mnum(numprocs_local)-1) == nr*nphi*nz) EXIT
         IF (i == numprocs_local) i = 1
         mnum(i) = mnum(i) + 1
         moffsets(i+1:numprocs_local) = moffsets(i+1:numprocs_local) + 1
         i=i+1
      END DO
      mystart = moffsets(myworkid+1)
      chunk  = mnum(myworkid+1)
      myend   = mystart + chunk - 1
!DEC$ ENDIF

      ! Master handles netCDF file
      IF (myworkid == master) THEN
         ! Open file
         i = 0
         mgrid_file = 'mgrid_'//TRIM(proc_string)// '.nc'
         CALL cdf_open(ngrid,mgrid_file,'w',i)
         if (i .ne. 0) CALL handle_err(NETCDF_OPEN_ERR,'stellopt_write_mgrid FILE:'//TRIM(mgrid_file),i)

         ! Define Coordinates
         ALLOCATE (vn_br(nextcur), vn_bz(nextcur), vn_bp(nextcur))
         CALL cdf_define(ngrid, vn_ir, nr)
         CALL cdf_define(ngrid, vn_jz, nz)
         CALL cdf_define(ngrid, vn_kp, nphi)
         CALL cdf_define(ngrid, vn_nfp, nfp)
         CALL cdf_define(ngrid, vn_nextcur, nextcur)
         CALL cdf_define(ngrid, vn_rmin, rmin)
         CALL cdf_define(ngrid, vn_zmin, zmin)
         CALL cdf_define(ngrid, vn_rmax, rmax)
         CALL cdf_define(ngrid, vn_zmax, zmax)
         CALL cdf_define(ngrid, vn_mgmode, mgrid_mode)
         CALL cdf_define(ngrid, vn_coilcur, extcur(1:nextcur),dimname=groupdim)
         IF (nextcur .eq. 1) THEN
            CALL cdf_define(ngrid, vn_coilgrp,coil_group(1)%s_name) 
         ELSE
            CALL cdf_define(ngrid, vn_coilgrp,coil_group(1:nextcur)%s_name, &
                            dimname=coildim)
         END IF
         DO ig = 1, nextcur
            write (temp, '(a,i3.3)') "_",ig
            vn_br(ig) = vn_br0 // temp
            vn_bp(ig) = vn_bp0 // temp
            vn_bz(ig) = vn_bz0 // temp
            CALL cdf_define(ngrid, vn_br(ig), br, dimname=cylcoord)
            CALL cdf_define(ngrid, vn_bp(ig), bp, dimname=cylcoord)
            CALL cdf_define(ngrid, vn_bz(ig), bz, dimname=cylcoord)
         END DO

         ! Begin writing data
         CALL cdf_write(ngrid, vn_ir, nr)
         CALL cdf_write(ngrid, vn_jz, nz)
         CALL cdf_write(ngrid, vn_kp, nphi)
         CALL cdf_write(ngrid, vn_nfp, nfp)
         CALL cdf_write(ngrid, vn_nextcur, nextcur)
         CALL cdf_write(ngrid, vn_rmin, rmin)
         CALL cdf_write(ngrid, vn_zmin, zmin)
         CALL cdf_write(ngrid, vn_rmax, rmax)
         CALL cdf_write(ngrid, vn_zmax, zmax)
         IF (nextcur .eq. 1) THEN
            CALL cdf_write(ngrid, vn_coilgrp, coil_group(1)%s_name)
         ELSE
            CALL cdf_write(ngrid, vn_coilgrp, coil_group(1:nextcur)%s_name)
         END IF
         CALL cdf_write(ngrid, vn_mgmode, mgrid_mode)
         CALL cdf_write(ngrid, vn_coilcur, extcur(1:nextcur))

      END IF

      ! Now loop over coil Groups
      DO ig = 1, nextcur
         DO s = mystart, myend
            !i = MOD(s-1,nr)+1
            !j = MOD(s-1,nr*nphi)
            !j = FLOOR(REAL(j) / REAL(nr))+1
            !k = CEILING(REAL(s) / REAL(nr*nphi))
            ! Changed from nr,nphi,nz to nr,nz,nphi ordering
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nz)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nz))
            r = rmin + dr*(i-1)
            z = zmin + dz*(j-1)
            phi = dphi*(k-1)
            CALL bfield(r,phi,z, br(i,j,k), bp(i,j,k), bz(i,j,k), IG = ig)
         END DO

         ! Gather the Results
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BARRIER(COMM_LOCAL,ierr_mpi)
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        br,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        COMM_LOCAL,ierr_mpi)
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        bp,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        COMM_LOCAL,ierr_mpi)
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        bz,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        COMM_LOCAL,ierr_mpi)
!DEC$ ENDIF
  

         ! Write to netCDF file
         IF (myworkid == master) THEN
            IF (.not. lasym) THEN
               !s = nz/2
               !DO j = 1,nz
               !   br(:,nz+1-j,:) = -br(:,j,:)
               !   bz(:,nz+1-j,:) =  bz(:,j,:)
               !   bp(:,nz+1-j,:) =  bp(:,j,:)
               !END DO
               
            END IF
            CALL cdf_write(ngrid, vn_br(ig), br)
            CALL cdf_write(ngrid, vn_bp(ig), bp)
            CALL cdf_write(ngrid, vn_bz(ig), bz)
         END IF

         IF (lscreen) THEN 
            IF (ABS(extcur(ig)).ge. 1.0E9) THEN
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(ig)%ncoil,'  EXTCUR = ',extcur(ig)/1.0E9,' [GA]'
            ELSE IF (ABS(extcur(ig)).ge. 1.0E6) THEN
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(ig)%ncoil,'  EXTCUR = ',extcur(ig)/1.0E6,' [MA]'
            ELSE IF (ABS(extcur(ig)).ge. 1.0E3) THEN
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(ig)%ncoil,'  EXTCUR = ',extcur(ig)/1.0E3,' [kA]'
            ELSE
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(ig)%ncoil,'  EXTCUR = ',extcur(ig),' [A]'
            END IF
            CALL FLUSH(6)
         END IF

      END DO

      ! Clsoe file
      IF (myworkid == master) CALL cdf_close(ngrid)
      CALL cleanup_biotsavart

      ! Free memory
      IF (ALLOCATED(br)) DEALLOCATE(br)
      IF (ALLOCATED(bp)) DEALLOCATE(bp)
      IF (ALLOCATED(bz)) DEALLOCATE(bz)
      IF (ALLOCATED(vn_br)) DEALLOCATE(vn_br)
      IF (ALLOCATED(vn_bp)) DEALLOCATE(vn_bp)
      IF (ALLOCATED(vn_bz)) DEALLOCATE(vn_bz)
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)

      ! Repoint MGRID
      mgrid_file_vmec = mgrid_file
      nzeta = nphi

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_write_mgrid
