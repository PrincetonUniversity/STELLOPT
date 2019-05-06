!-----------------------------------------------------------------------
!     Module:        mpi_sharmem
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/05/2019
!     Description:   This module helps simplify creation of shared
!                    memory arrays using MPI-3.
!-----------------------------------------------------------------------
MODULE MPI_SHARMEM
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
   IMPLICIT NONE
!-----------------------------------------------------------------------
!     Subroutines
!         mpialloc_1d:    Allocate a 1D variable (int/double)
!         mpialloc_2d:    Allocate a 2D variable (int/double)
!         mpialloc_3d:    Allocate a 3D variable (int/double)
!-----------------------------------------------------------------------
   INTERFACE mpialloc_1d
      MODULE PROCEDURE mpialloc_1d_int, mpialloc_1d_dbl
   END INTERFACE
   INTERFACE mpialloc_2d
      MODULE PROCEDURE mpialloc_2d_int, mpialloc_2d_dbl
   END INTERFACE
   INTERFACE mpialloc_3d
      MODULE PROCEDURE mpialloc_3d_int, mpialloc_3d_dbl
   END INTERFACE
   CONTAINS

      SUBROUTINE mpialloc_1d_int(array,n1,subid,mymaster,share_comm,win)
      ! Libraries
!DEC$ IF DEFINED (MPI_OPT)
      USE MPI
!DEC$ ENDIF
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      INTEGER, POINTER, INTENT(inout) :: array(:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(inout) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(1)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      disp_unit = 1
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
      RETURN
      END SUBROUTINE mpialloc_1d_int

      SUBROUTINE mpialloc_1d_dbl(array,n1,subid,mymaster,share_comm,win)
      ! Libraries
!DEC$ IF DEFINED (MPI_OPT)
      USE MPI
!DEC$ ENDIF
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      DOUBLE PRECISION, POINTER, INTENT(inout) :: array(:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(inout) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(1)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      disp_unit = 1
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
      RETURN
      END SUBROUTINE mpialloc_1d_dbl

      SUBROUTINE mpialloc_2d_int(array,n1,n2,subid,mymaster,share_comm,win)
      ! Libraries
!DEC$ IF DEFINED (MPI_OPT)
      USE MPI
!DEC$ ENDIF
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      INTEGER, POINTER, INTENT(inout) :: array(:,:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: n2
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(inout) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(2)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      array_shape(2) = n2
      disp_unit = 1
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1*n2,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
      RETURN
      END SUBROUTINE mpialloc_2d_int

      SUBROUTINE mpialloc_2d_dbl(array,n1,n2,subid,mymaster,share_comm,win)
      ! Libraries
!DEC$ IF DEFINED (MPI_OPT)
      USE MPI
!DEC$ ENDIF
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      DOUBLE PRECISION, POINTER, INTENT(inout) :: array(:,:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: n2
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(inout) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(2)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      array_shape(2) = n2
      disp_unit = 1
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1*n2,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
      RETURN
      END SUBROUTINE mpialloc_2d_dbl

      SUBROUTINE mpialloc_3d_int(array,n1,n2,n3,subid,mymaster,share_comm,win)
      ! Libraries
!DEC$ IF DEFINED (MPI_OPT)
      USE MPI
!DEC$ ENDIF
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      INTEGER, POINTER, INTENT(inout) :: array(:,:,:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: n2
      INTEGER, INTENT(in) :: n3
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(inout) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(3)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      array_shape(2) = n2
      array_shape(3) = n3
      disp_unit = 1
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1*n2*n3,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
      RETURN
      END SUBROUTINE mpialloc_3d_int

      SUBROUTINE mpialloc_3d_dbl(array,n1,n2,n3,subid,mymaster,share_comm,win)
      ! Libraries
!DEC$ IF DEFINED (MPI_OPT)
      USE MPI
!DEC$ ENDIF
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      DOUBLE PRECISION, POINTER, INTENT(inout) :: array(:,:,:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: n2
      INTEGER, INTENT(in) :: n3
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(inout) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(3)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      array_shape(2) = n2
      array_shape(3) = n3
      disp_unit = 1
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1*n2*n3,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
      RETURN
      END SUBROUTINE mpialloc_3d_dbl

      SUBROUTINE mpidealloc_1d_dbl(array,win)
      ! Libraries
!DEC$ IF DEFINED (MPI_OPT)
      USE MPI
!DEC$ ENDIF
      IMPLICIT NONE
      ! Arguments
      DOUBLE PRECISION, POINTER, INTENT(inout) :: array(:)
      INTEGER, INTENT(inout) :: win
      INTEGER :: ier
      CALL MPI_WIN_FREE(win,ier)
      IF (ASSOCIATED(array)) NULLIFY(array)
      RETURN
      END SUBROUTINE mpidealloc_1d_dbl
   
END MODULE MPI_SHARMEM
