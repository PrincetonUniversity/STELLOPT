SUBROUTINE mpialloc_2d_int(array,n1,n2,subid,mymaster,share_comm,win)
    ! Libraries
    USE MPI
    USE ISO_C_BINDING

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(inout), POINTER :: array(:,:)
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

    ier = 0
    array_shape(1) = n1; array_shape(2) = n2;
    disp_unit = 8_MPI_ADDRESS_KIND
       window_size = 0_MPI_ADDRESS_KIND
    IF (subid == mymaster) window_size = INT(PRODUCT(ARRAY_SHAPE),MPI_ADDRESS_KIND)
    CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
    IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
    CALL C_F_POINTER(baseptr, array, array_shape)
    CALL MPI_WIN_FENCE(0, win, ier)

    RETURN

END SUBROUTINE mpialloc_2d_int
