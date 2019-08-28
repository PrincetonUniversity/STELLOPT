!-----------------------------------------------------------------------
!     Module:        beams3d_init_restart
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/27/2012
!     Description:   This subroutine loads the fields from a restart file.
!                    This is still under development.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_init_restart
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime
      USE beams3d_grid
!DEC$ IF DEFINED (LHDF5)
      USE ez_hdf5
!DEC$ ENDIF  
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (lverb) THEN
         WRITE(6,'(A)')  '----- Reading Restart File -----'
         WRITE(6,'(A)')  '   FILE: '//TRIM(restart_string)
      END IF
!      CALL open_hdf5(TRIM(restart_string),fid,ier)
!      CALL read_var_hdf5(fid,'nr',ier,INTVAR=nr)
!      CALL read_var_hdf5(fid,'nphi',ier,INTVAR=nphi)
!      CALL read_var_hdf5(fid,'nz',ier,INTVAR=nz)
!      ALLOCATE(raxis(nr),phiaxis(nphi),zaxis(nz))
!      CALL read_var_hdf5(fid,'raxis',nr,ier,DBLVAR=raxis)
!      CALL read_var_hdf5(fid,'phiaxis',nphi,ier,DBLVAR=phiaxis)
!      CALL read_var_hdf5(fid,'zaxis',nz,ier,DBLVAR=zaxis)
!      ALLOCATE(B_R(nr,nphi,nz),B_Z(nr,nphi,nz))
!      CALL read_var_hdf5(fid,'B_R',nr,nphi,nz,ier,DBLVAR=B_R)
!      CALL read_var_hdf5(fid,'B_Z',nr,nphi,nz,ier,DBLVAR=B_Z)
!      CALL close_hdf5(fid,ier)
      ! No redefine id_string
      id_string = TRIM(id_string) // '_new'
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_init_restart
