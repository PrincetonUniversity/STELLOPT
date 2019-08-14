!-----------------------------------------------------------------------
!     Module:        fieldlines_init_restart
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/27/2012
!     Description:   This subroutine loads the fields from a restart file.
!                    This is still under development.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_restart
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_runtime
      USE fieldlines_grid
      USE fieldlines_lines, ONLY: nlines
      USE ez_hdf5
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, i
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (lverb) THEN
         WRITE(6,'(A)')  '----- Reading Restart File -----'
         WRITE(6,'(A)')  '   FILE: '//TRIM(restart_string)
      END IF
      !PRINT *,'myid:',myid,TRIM(restart_string)
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
!DEC$ ENDIF
!DEC$ IF DEFINED (LHDF5)
      IF (myid == master) THEN
      CALL open_hdf5(TRIM(restart_string),fid,ier)
      PRINT *,'ier',ier
      CALL read_scalar_hdf5(fid,'nr',ier,INTVAR=nr)
      CALL read_scalar_hdf5(fid,'nphi',ier,INTVAR=nphi)
      CALL read_scalar_hdf5(fid,'nz',ier,INTVAR=nz)
      ALLOCATE(raxis(nr),phiaxis(nphi),zaxis(nz))
      CALL read_var_hdf5(fid,'raxis',nr,ier,DBLVAR=raxis)
      CALL read_var_hdf5(fid,'phiaxis',nphi,ier,DBLVAR=phiaxis)
      CALL read_var_hdf5(fid,'zaxis',nz,ier,DBLVAR=zaxis)
      ALLOCATE(B_R(nr,nphi,nz),B_PHI(nr,nphi,nz),B_Z(nr,nphi,nz))
      CALL read_var_hdf5(fid,'B_R',nr,nphi,nz,ier,DBLVAR=B_R)
      CALL read_var_hdf5(fid,'B_PHI',nr,nphi,nz,ier,DBLVAR=B_PHI)
      CALL read_var_hdf5(fid,'B_Z',nr,nphi,nz,ier,DBLVAR=B_Z)
      CALL close_hdf5(fid,ier)
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
!DEC$ ENDIF
      rmin = MINVAL(raxis)
      rmax = MAXVAL(raxis)
      phimin = MINVAL(phiaxis)
      phimax = MAXVAL(phiaxis)
      zmin   = MINVAL(zaxis)
      zmax   = MAXVAL(zaxis)
      IF (lverb) THEN
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   R   = [',rmin,',',rmax,'];  NR:   ',nr
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   PHI = [',phimin,',',phimax,'];  NPHI: ',nphi
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   Z   = [',zmin,',',zmax,'];  NZ:   ',nz
         !IF (lauto) WRITE(6,'(A)') '   AUTO CALCULATED STARTING POINTS!'
         WRITE(6,'(A,I4)')               '   # of Fieldlines: ',nlines
         CALL FLUSH(6)
      END IF
      RETURN
      !IF (lauto) THEN
      !   nlines = nr
      !   IF (nlines > MAXLINES) nlines = MAXLINES-2
      !   DO i = 1, nlines
      !      r_start(i)   = i*(rmax-rmin)/nlines + rmin
      !      z_start(i)   = 0.0
      !      phi_start(i) = 0.0
      !      phi_end(i)   = 2*pi*1000.
      !   END DO
      !END IF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_restart
