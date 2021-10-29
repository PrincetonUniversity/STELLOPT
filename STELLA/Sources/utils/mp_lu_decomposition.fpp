module mp_lu_decomposition

#if defined MPI && ISO_C_BINDING

  implicit none

  public :: lu_decomposition_local
  public :: lu_inverse_local
  public :: lu_matrix_multiply_local

  interface lu_decomposition_local
!    module procedure lu_decomposition_local_real
     module procedure lu_decomposition_local_complex
  end interface

  interface lu_inverse_local
!    module procedure lu_inverse_local_real
     module procedure lu_inverse_local_complex
  end interface

  interface lu_matrix_multiply_local
!    module procedure lu_decomposition_local_real
     module procedure lu_matrix_multiply_local_complex
  end interface

contains 

  subroutine lu_decomposition_local_complex (mp_comm, root, win, lu, idx, d)

    use mpi

    implicit none
    
    integer, intent (in) :: win, mp_comm, root
    complex, dimension (:,:), intent (in out) :: lu
    integer, dimension (:), intent (out) :: idx
    real, intent (out) :: d

    real, parameter :: zero = 1.0e-20
    real, dimension (size(lu,1)) :: vv
    complex, dimension (size(lu,2)) :: dum
    integer, dimension(:), allocatable :: row_limits

    integer :: i, j, k, n, imax, rdiv, rmod
    integer :: iproc, nproc, ierr
    real :: dmax, tmp

    n = size(lu,1)

    call mpi_comm_size(mp_comm, nproc, ierr)
    call mpi_comm_rank(mp_comm, iproc, ierr)

    allocate (row_limits(0:nproc))

    d = 1.0
    vv = maxval(cabs(lu),dim=2)
    if (any(vv==0.0)) &
      write (*,*) 'singular matrix in lu_decomposition on process ', iproc
    vv = 1.0/vv
    do j = 1, n
        !divide up the work using row_limits
        rdiv = (n-j)/nproc
        rmod = mod(n-j,nproc)
        row_limits(0) = j+1
        if(rdiv.eq.0) then
          row_limits(rmod+1:) = -1
          do k=1,rmod
            row_limits(k)  = row_limits(k-1) + 1
          enddo
        else
          do k=1,nproc
            row_limits(k) = row_limits(k-1) + rdiv
            if(k.le.rmod) row_limits(k) = row_limits(k) + 1
          enddo
        endif

        !pivot if needed
        dmax = -1.0
        do k = j, n
          tmp = vv(k)*abs(lu(k,j))
          if(tmp.gt.dmax) then 
            dmax = tmp
            imax = k
          endif
        enddo

        if(iproc.eq.root) then
          idx(j) = imax
          if (j /= imax) then
            dum = lu(imax,:)
            lu(imax,:) = lu(j,:)
            lu(j,:) = dum
            vv(imax) = vv(j)
            d = -d
          end if
          if (lu(j,j)==0.0) lu(j,j) = zero
        else
          if (j /= imax) vv(imax) = vv(j)
        endif

        call mpi_win_fence(0,win,ierr)

        !get the lead multiplier
        do i = row_limits(iproc), row_limits(iproc+1)-1
          lu(i,j) = lu(i,j)/lu(j,j)
        enddo
         
        call mpi_win_fence(0,win,ierr)

        do k=row_limits(iproc), row_limits(iproc+1)-1
          do i = j+1,n
            lu(i,k) = lu(i,k) - lu(i,j)*lu(j,k)
          enddo
        enddo

        call mpi_win_fence(0,win,ierr)
    enddo

    deallocate (row_limits)

  end subroutine lu_decomposition_local_complex

  subroutine lu_inverse_local_complex (mp_comm, root, win, lu, idx, inverse)

    use linear_solve, only: lu_back_substitution

    implicit none

    integer, intent (in) :: win, mp_comm, root
    complex, dimension (:,:), intent (in) :: lu
    integer, dimension (:), intent (in) :: idx
    complex, dimension (:,:), intent (out) :: inverse

    integer :: i, n, nproc, iproc, rdiv, rmod, ierr
    integer :: llim, ulim

    n = size(lu,1)

    call mpi_comm_size(mp_comm, nproc, ierr)
    call mpi_comm_rank(mp_comm, iproc, ierr)

    rdiv = n/nproc
    rmod = mod(n,nproc)

    llim = 1 + iproc*rdiv + min(rmod,iproc)
    ulim = (iproc+1)*rdiv + min(rmod,iproc+1)

    do i = llim, ulim
      inverse(i,:) = 0
      inverse(i,i) = 1.0
    enddo

    call mpi_win_fence(0,win,ierr)

    do i = llim, ulim
      call lu_back_substitution(lu,idx,inverse(:,i))
    enddo

    call mpi_win_fence(0,win,ierr)

  end subroutine lu_inverse_local_complex

  subroutine lu_matrix_multiply_local_complex (mp_comm, root, win, mat, b)

    implicit none

    integer, intent (in) :: win, mp_comm, root
    complex, dimension (:,:), intent (in) :: mat
    complex, dimension (:), intent (out) :: b
    complex, dimension (size(b)) :: a

    integer :: i, n, nproc, iproc, rdiv, rmod
    integer :: llim, ulim, ierr

    n = size(mat,1)

    call mpi_comm_size(mp_comm, nproc, ierr)
    call mpi_comm_rank(mp_comm, iproc, ierr)

    rdiv = n/nproc
    rmod = mod(n,nproc)

    llim = 1 + iproc*rdiv + min(rmod,iproc)
    ulim = (iproc+1)*rdiv + min(rmod,iproc+1)

    do i = llim, ulim
      a(i) = sum(mat(i,:)*b(:))
    enddo

    call mpi_win_fence(0,win,ierr)

    do i = llim, ulim
      b(i) = a(i)
    enddo

    call mpi_win_fence(0,win,ierr)

  end subroutine lu_matrix_multiply_local_complex

#endif

end module mp_lu_decomposition
