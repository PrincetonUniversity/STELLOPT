!-----------------------------------------------------------------------
!     Subroutine:    follow
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/8/2011
!     Description:   This subroutine preforms field line following in
!                    the TORLINES background coordinates system.
!-----------------------------------------------------------------------
      SUBROUTINE follow
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE torlines_fieldlines
      USE torlines_background
      USE torlines_realspace
      USE torlines_runtime
      !USE pies_magco
      USE EZspline_obj
      USE EZspline
      USE ez_hdf5
      USE mpi_params
!-----------------------------------------------------------------------
!     Local Variables
!          u           Poloidal dummy index
!          v           Toroidal dummy index
!          ierr        Error flag
!          bcs1/2/3    Boundary Condition Arrays for EZspline
!          factor      Factor used to convert B^v B^phi
!          brho_real   bsreal/bvreal
!          btheta_real bureal/bvreal
!-----------------------------------------------------------------------
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'
!DEC$ ENDIF  
      INTEGER :: ier,u,v,ik,form,mn
      REAL(rprec) :: r_good, r_bad, dr
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (bound_separation > 1) THEN
         IF (lverb) WRITE(6,'(a)')'----     Rough GRID     ----'
         ! First follow fieldlines from axis to edge
         nlines  = 201
         r_start = -1; phi_start = 0; z_start = 0; phi_end = -1;
         DO ik = 1, nlines
            r_start(ik)   = REAL(ik)/REAL(nlines+1)  ! We know the axis is good and the edge must be bad.
            z_start(ik)   = 0.0
            phi_start(ik) = phimn
            phi_end(ik)   = 100*phmx
         END DO
         IF (lverb) WRITE(6,'(A,F8.5,A,F8.5,A)') '   rho   = [',MINVAL(r_start,MASK=(r_start >=0)),',',MAXVAL(r_start,MASK=(r_start >=0)),'];'
         CALL FLUSH(6)
         CALL torlines_follow
         
         ! From last good line to first bad
         IF (myid == master) THEN
            ik     = COUNT(goodline)
            IF (ik == nlines) THEN
               r_good = 1.0
               r_bad  = 1.0
            ELSE
               r_good = r_start(ik)
               r_bad  = r_start(ik+1)
               dr     = r_bad-r_good
            END IF
         END IF
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'follow',ierr_mpi)
         CALL MPI_BCAST(r_good,1,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'follow',ierr_mpi)
         CALL MPI_BCAST(r_bad,1,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'follow',ierr_mpi)
         CALL MPI_BCAST(dr,1,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'follow',ierr_mpi)
!DEC$ ENDIF
         nlines  = 201
         r_start = -1; phi_start = 0; z_start = 0; phi_end = -1;
         IF (r_good /= r_bad) THEN
            IF (lverb) WRITE(6,*) ''
            IF (lverb) WRITE(6,'(a)')'----     Edge GRID     ----'
            DO ik = 1, nlines
               r_start(ik)   = r_good+dr*REAL(ik)/REAL(nlines+1)
               z_start(ik)   = 0.0
               phi_start(ik) = phimn
               phi_end(ik)   = 100*phmx
            END DO
            IF (lverb) WRITE(6,'(A,F8.5,A,F8.5,A)') '   rho   = [',MINVAL(r_start,MASK=(r_start >=0)),',',MAXVAL(r_start,MASK=(r_start >=0)),'];'
            CALL FLUSH(6)
            CALL torlines_follow
         END IF
         
         ! Now from axis to edge
         IF (lverb) WRITE(6,*) ''
         IF (lverb) WRITE(6,'(a)')'----     Final GRID     ----'
         IF (myid == master) THEN
            ik     = COUNT(goodline)
            IF (ik == nlines) THEN
               r_good = 1.0
            ELSE
               r_good = r_start(ik+2)
            END IF
         END IF
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'follow',ierr_mpi)
         CALL MPI_BCAST(r_good,1,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'follow',ierr_mpi)
!DEC$ ENDIF
         nlines  = 201
         r_start = -1; phi_start = 0; z_start = 0; phi_end = -1;
         DO ik = 1, nlines
            r_start(ik)   = r_good*(sqrt(REAL(ik-1)/REAL(nlines-1)))
            z_start(ik)   = 0.0
            phi_start(ik) = phimn
            phi_end(ik)   = 1000*phmx
            !WRITE(6,'(4(ES20.10))') r_start(ik), z_start(ik),phi_start(ik),phi_end(ik)
         END DO
         IF (lverb) WRITE(6,'(A,F8.5,A,F8.5,A)') '   rho   = [',MINVAL(r_start,MASK=(r_start >=0)),',',MAXVAL(r_start,MASK=(r_start >=0)),'];'
         CALL FLUSH(6)
         CALL torlines_follow
      ELSE
         IF (lverb) WRITE(6,'(a)')'----     BASIC GRID     ----'
         ! First follow fieldlines from axis to edge
         nlines  = 201
         r_start = -1; phi_start = 0; z_start = 0; phi_end = -1;
         DO ik = 1, nlines
            r_start(ik)   = REAL(ik-1)/REAL(nlines-1)
            z_start(ik)   = 0.0
            phi_start(ik) = phimn
            phi_end(ik)   = 1000*phmx
         END DO
         IF (lverb) WRITE(6,'(A,F8.5,A,F8.5,A)') '   rho   = [',MINVAL(r_start,MASK=(r_start >=0)),',',MAXVAL(r_start,MASK=(r_start >=0)),'];'
         CALL FLUSH(6)
         CALL torlines_follow
      END IF
      
      CALL torlines_write('FIELDLINES')

      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE follow
