! Modifications for optimised local copy in c_redist_22 and c_redist_32 
! (and their inverse routines):
! (c) The Numerical Algorithms Group (NAG) Ltd, 2012
! on behalf of EPSRC for the HECToR project
module redistribute
!
! Redistribute distributed (integer, real, complex or logical) 
! (1, 2, 3, or 4) dimensional arrays into two dimensional arrays with 
! first index on local processor, and vice versa.
!
! The first operation is called 'gather' and the second is called 'scatter.'
!  
! One can also do a 'fill' operation.  This consists of copying 
! values from a (2, 3, or 4) dimensional array of 
! (integer, real, complex, or logical ) values into 
! another array with the same number of dimensions.
!
! One can also do a three index to four index redistribution for complex numbers.
!
  implicit none
  private

  public :: index_list_type, delete_list
  public :: redist_type, delete_redist
! TT>
  public :: report_map_property, measure_gather, measure_scatter
  public :: gather_count, scatter_count, time_redist
! <TT

  public :: init_redist, gather, scatter
  public :: init_fill, fill
  public :: set_redist_character_type

  interface gather
     module procedure c_redist_22, r_redist_22, i_redist_22, l_redist_22
     module procedure c_redist_32, r_redist_32, i_redist_32, l_redist_32
     module procedure c_redist_42, r_redist_42, i_redist_42, l_redist_42
     module procedure c_redist_23
     module procedure c_redist_34, r_redist_34
     module procedure c_redist_33
     module procedure c_redist_35, r_redist_35
  end interface

  interface scatter
     module procedure c_redist_12,     r_redist_12,     i_redist_12,     l_redist_12
     module procedure c_redist_22_inv, r_redist_22_inv, i_redist_22_inv, l_redist_22_inv
     module procedure c_redist_32_inv, r_redist_32_inv, i_redist_32_inv, l_redist_32_inv
     module procedure c_redist_42_inv, r_redist_42_inv, i_redist_42_inv, l_redist_42_inv
     module procedure c_redist_33_inv
     module procedure c_redist_34_inv, r_redist_34_inv
     module procedure c_redist_35_inv, r_redist_35_inv
  end interface

! TT>
  interface measure_gather
     module procedure measure_gather_32, measure_gather_33
     module procedure measure_gather_22
  end interface

  interface measure_scatter
     module procedure measure_scatter_23, measure_scatter_33
     module procedure measure_scatter_22
  end interface

  integer :: gather_count=0, scatter_count=0
  real, save :: time_redist(2)=0.
! <TT

  interface fill
     module procedure c_fill_2, c_fill_3, c_fill_4
     module procedure r_fill_2, r_fill_3, r_fill_4
     module procedure i_fill_2, i_fill_3, i_fill_4
     module procedure l_fill_2, l_fill_3, l_fill_4
  end interface

  type :: index_map
     integer :: nn
     integer, dimension (:), pointer :: k => null()
     integer, dimension (:), pointer :: l => null()
     integer, dimension (:), pointer :: m => null()
     integer, dimension (:), pointer :: n => null()
     integer, dimension (:), pointer :: o => null()
  end type index_map

  ! TT: want to add map name, from_layout and to_layout
  type :: redist_type
     private
     integer, dimension (5) :: to_low, from_low, to_high, from_high
     type (index_map), dimension (:), pointer :: to   => null()
     type (index_map), dimension (:), pointer :: from => null()
     complex, dimension (:), pointer :: complex_buff  => null()
     real,    dimension (:), pointer :: real_buff     => null()
     integer, dimension (:), pointer :: integer_buff  => null()
     logical, dimension (:), pointer :: logical_buff  => null()
     character (len=3) :: redistname = ""
  end type redist_type

  type :: index_list_type
     integer, dimension (:), pointer :: first  => null()
     integer, dimension (:), pointer :: second => null()
     integer, dimension (:), pointer :: third  => null()
     integer, dimension (:), pointer :: fourth => null()
     integer, dimension (:), pointer :: fifth => null()
  end type index_list_type


contains

  subroutine set_redist_character_type(r, chartype)

    type (redist_type), intent (inout) :: r
    character (3), intent (in) :: chartype

    r%redistname = chartype

  end subroutine set_redist_character_type

  subroutine init_redist (r, char, to_low, to_high, to_list, &
       from_low, from_high, from_list, ierr)

    use mp, only: iproc, nproc, proc0
    type (redist_type), intent (inout) :: r
    character(1), intent (in) :: char
!    integer, intent (in) :: to_low
! TT> caused a problem on PGI compiler
!    type (index_list_type), dimension (0:) :: to_list, from_list
    type (index_list_type), dimension (0:nproc-1), intent (in) :: to_list, from_list
! <TT
    integer, dimension(:), intent (in) :: from_low, to_high, from_high, to_low
    ! TT: to_high, from_high seem never used

    integer :: j, ip, n_to, n_from, buff_size
    integer, optional, intent (out) :: ierr

    allocate (r%to(0:nproc-1), r%from(0:nproc-1))

    if (present(ierr)) ierr = 0
    buff_size = 0

! TT> added possibility of higher dimension
!    r%to_low(1) = 1
!    r%to_low(2) = to_low ! TT: assumes two-dimensionality
    ! Actually, did above work for redist_23 or 34?
    ! It is used to determine the lower bound of the index
    ! when redist receives array
    ! The new description below hangs up if dimension of to_high is
    ! set incorrectly
!    r%to_low(:) = 1
!    r%to_low(size(to_high)) = to_low
!    print *, 'size(to_high) is', size(to_high), iproc
!    print *, 'r%to_low', r%to_low, iproc
! <TT

    do ip = 0, nproc - 1
       if (associated(to_list(ip)%first)) then
          n_to = size(to_list(ip)%first)
          r%to(ip)%nn = n_to

          allocate (r%to(ip)%k(n_to))
          allocate (r%to(ip)%l(n_to))

          r%to(ip)%k = to_list(ip)%first
          r%to(ip)%l = to_list(ip)%second

          if (associated(to_list(ip)%third)) then
             allocate (r%to(ip)%m(n_to))
             r%to(ip)%m = to_list(ip)%third
          end if

          if (associated(to_list(ip)%fourth)) then
             allocate (r%to(ip)%n(n_to))
             r%to(ip)%n = to_list(ip)%fourth
          end if

          if (associated(to_list(ip)%fifth)) then
             allocate (r%to(ip)%o(n_to))
             r%to(ip)%o = to_list(ip)%fifth
          end if

          if (ip /= iproc) buff_size = max(buff_size, n_to)
       else
          r%to(ip)%nn = 0
       endif
    enddo

    do j = 1, size(from_low)
       r%from_low(j) = from_low(j)
    enddo

    do j = 1, size(from_high)
       r%from_high(j) = from_high(j)
    enddo

    do j = 1, size(to_high)
       r%to_high(j) = to_high(j)
    enddo

    do j = 1, size(to_low)
       r%to_low(j) = to_low(j)
    end do

    do ip = 0, nproc - 1
       if (associated(from_list(ip)%first)) then

          n_from = size(from_list(ip)%first)
          r%from(ip)%nn = n_from

          allocate (r%from(ip)%k(n_from))
          allocate (r%from(ip)%l(n_from))

          r%from(ip)%k = from_list(ip)%first
          r%from(ip)%l = from_list(ip)%second

          if (associated (from_list(ip)%third)) then
             allocate (r%from(ip)%m(n_from))
             r%from(ip)%m = from_list(ip)%third
          endif

          if (associated (from_list(ip)%fourth)) then
             allocate (r%from(ip)%n(n_from))
             r%from(ip)%n = from_list(ip)%fourth
          endif

          if (associated (from_list(ip)%fifth)) then
             allocate (r%from(ip)%o(n_from))
             r%from(ip)%o = from_list(ip)%fifth
          endif

          if (ip /= iproc) buff_size = max(buff_size, n_from)
       else
          r%from(ip)%nn = 0
       endif
    enddo

    select case (char)
       case ('c') 
          if (buff_size > 0) allocate (r%complex_buff(buff_size))
       case ('r') 
          if (buff_size > 0) allocate (r%real_buff(buff_size))
       case ('i') 
          if (buff_size > 0) allocate (r%integer_buff(buff_size))
       case ('l') 
          if (buff_size > 0) allocate (r%logical_buff(buff_size))
       case default
          if (proc0) then
             write(*,*) 'Type to be redistributed invalid.  Must stop.'
             write(*,*) char
          endif
          stop
    end select

  end subroutine init_redist

  !  subroutine init_fill (f, char, to_low, to_high, to_list, &
  !       from_low, from_high, from_list, ierr)
  subroutine init_fill (f, char, to_low, to_list, &
       from_low, from_list, ierr)

    use mp, only: nproc, proc0, iproc
    type (redist_type), intent (out) :: f
    character(1), intent (in) :: char
! TT> caused a problem on PGI compiler
!    type (index_list_type), dimension (0:) :: to_list, from_list
    type (index_list_type), dimension (0:nproc-1), intent (in) :: to_list, from_list
! <TT
    !    integer, dimension(:), intent (in) :: to_low, from_low, to_high, from_high
        integer, dimension(:), intent (in) :: to_low, from_low
    integer, optional, intent (out) :: ierr
    
    integer :: j, ip, n_to, n_from, buff_size

    if (present(ierr)) ierr = 0
    do j = 1, size(to_low)
       f%to_low(j) = to_low(j)
    enddo

    do j = 1, size(from_low)
       f%from_low(j) = from_low(j)
    enddo

    allocate (f%to(0:nproc-1), f%from(0:nproc-1))

    buff_size = 0
    do ip = 0, nproc - 1
       if (associated(to_list(ip)%first)) then
          n_to = size(to_list(ip)%first)
          f%to(ip)%nn = n_to
          allocate (f%to(ip)%k(n_to))
          f%to(ip)%k = to_list(ip)%first

          if (associated (to_list(ip)%second)) then
             allocate (f%to(ip)%l(n_to))
             f%to(ip)%l = to_list(ip)%second
          endif

          if (associated (to_list(ip)%third)) then
             allocate (f%to(ip)%m(n_to))
             f%to(ip)%m = to_list(ip)%third
          endif

          if (associated (to_list(ip)%fourth)) then
             allocate (f%to(ip)%n(n_to))
             f%to(ip)%n = to_list(ip)%fourth
          endif

          if (associated (to_list(ip)%fifth)) then
             allocate (f%to(ip)%o(n_to))
             f%to(ip)%o = to_list(ip)%fifth
          endif

          if (ip /= iproc) buff_size = max(buff_size, n_to)
       else
          f%to(ip)%nn = 0
       endif
    enddo

    do ip = 0, nproc - 1
       if (associated(from_list(ip)%first)) then
          n_from = size(from_list(ip)%first)
          f%from(ip)%nn = n_from
          allocate (f%from(ip)%k(n_from))
          f%from(ip)%k = from_list(ip)%first

          if (associated (from_list(ip)%second)) then
             allocate (f%from(ip)%l(n_from))
             f%from(ip)%l = from_list(ip)%second
          endif

          if (associated (from_list(ip)%third)) then
             allocate (f%from(ip)%m(n_from))
             f%from(ip)%m = from_list(ip)%third
          endif

          if (associated (from_list(ip)%fourth)) then
             allocate (f%from(ip)%n(n_from))
             f%from(ip)%n = from_list(ip)%fourth
          endif

          if (associated (from_list(ip)%fifth)) then
             allocate (f%from(ip)%o(n_from))
             f%from(ip)%o = from_list(ip)%fifth
          endif

          if (ip /= iproc) buff_size = max(buff_size, n_from)
       else
          f%from(ip)%nn = 0
       endif
    enddo

    select case (char)
       case ('c') 
          if (buff_size > 0) allocate (f%complex_buff(buff_size))          
       case ('r') 
          if (buff_size > 0) allocate (f%real_buff(buff_size))
       case ('i') 
          if (buff_size > 0) allocate (f%integer_buff(buff_size))
       case ('l') 
          if (buff_size > 0) allocate (f%logical_buff(buff_size))
       case default
          if (proc0) then
             write(*,*) 'Type to be redistributed invalid.  Must stop.'
             write(*,*) char
          endif
          stop
    end select
   
  end subroutine init_fill

  subroutine delete_redist(r)

    use mp, only: nproc
    type (redist_type), intent (in out) :: r
    
    integer :: i

    if (associated(r%to)) then
       do i = 0, nproc-1
          if (associated (r%to(i)%k)) deallocate (r%to(i)%k)
          if (associated (r%to(i)%l)) deallocate (r%to(i)%l)
          if (associated (r%to(i)%m)) deallocate (r%to(i)%m)
          if (associated (r%to(i)%n)) deallocate (r%to(i)%n)
       enddo
       deallocate (r%to)
    endif

    if (associated(r%from)) then
       do i = 0, nproc-1
          if (associated (r%from(i)%k)) deallocate (r%from(i)%k)
          if (associated (r%from(i)%l)) deallocate (r%from(i)%l)
          if (associated (r%from(i)%m)) deallocate (r%from(i)%m)
          if (associated (r%from(i)%n)) deallocate (r%from(i)%n)
       enddo
       deallocate (r%from)
    endif

    if (associated (r%complex_buff)) deallocate (r%complex_buff)
    if (associated (r%real_buff))    deallocate (r%real_buff)
    if (associated (r%integer_buff)) deallocate (r%integer_buff)
    if (associated (r%logical_buff)) deallocate (r%logical_buff)

  end subroutine delete_redist

  subroutine delete_list(list)
    use mp, only: nproc
! TT> caused a problem on PGI compiler
!    type (index_list_type), dimension(0:) :: list
    type (index_list_type), dimension(0:nproc-1), intent (inout) :: list
! <TT

    integer :: ip

    do ip = 0, nproc-1
       if(associated(list(ip)%first)) deallocate(list(ip)%first)
       if(associated(list(ip)%second)) deallocate(list(ip)%second)
       if(associated(list(ip)%third)) deallocate(list(ip)%third)
       if(associated(list(ip)%fourth)) deallocate(list(ip)%fourth)
       if(associated(list(ip)%fifth)) deallocate(list(ip)%fifth)
    enddo

  end subroutine delete_list

  subroutine c_redist_12 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):), intent (in) :: from_here

    complex, dimension (r%to_low(1):,r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i), r%to(iproc)%l(i)) &
            = from_here(r%from(iproc)%k(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i)) 
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_redist_12

  subroutine c_redist_22 (r, from_here, to_here)

    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    call c_redist_22_old_copy(r, from_here, to_here)
    call c_redist_22_mpi_copy(r, from_here, to_here)

  end subroutine c_redist_22

  subroutine c_redist_22_old_copy (r, from_here, to_here)

    use mp, only: iproc
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i

!CMR 
! In the stella standard FFT situation this routine maps 
!         xxf(it,ixxf) to yxf(ik,iyxf) data type 
!         where it is kx (or x) index, ik is ky (or y) index, 
!         ixxf is (y,ig,isgn,"les") and iyxf is (x,ig,isgn,"les") 
!
    do i = 1, r%from(iproc)%nn
!
! redistribute from local processor to local processor
! NB r%from(iproc)%nn is #elements sent by THIS processor to THIS processor
!    In this situation the data at (r%to(iproc)%k(i),r%to(iproc)%l(i)) 
!    should come from (r%from(iproc)%k(i), r%from(iproc)%l(i)). 
!
! This do loop, in stella standard FFT situation, corresponds to:
!    to_here(ik,iyxf)=from_here(it,ixxf)
!
       to_here(r%to(iproc)%k(i), &
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i))
    end do

  end subroutine c_redist_22_old_copy

  subroutine c_redist_22_mpi_copy (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_redist_22_mpi_copy


  subroutine c_redist_22_inv (r, from_here, to_here)

    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here

    call c_redist_22_inv_old_copy(r, from_here, to_here)
    call c_redist_22_inv_mpi_copy(r, from_here, to_here)

  end subroutine c_redist_22_inv

  subroutine c_redist_22_inv_old_copy (r, from_here, to_here)

    use mp, only: iproc
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here

    integer :: i

!CMR 
! In the stella standard FFT situation this routine maps 
!         yxf(ik,iyxf) to xxf(it,ixxf) data type 
!         where it is kx (or x) index, ik is ky (or y) index, 
!         ixxf is (y,ig,isgn,"les") and iyxf is (x,ig,isgn,"les") 
!
    do i = 1, r%to(iproc)%nn
!
! redistribute from local processor to local processor
! NB r%from(iproc)%nn is #elements sent by THIS processor to THIS processor
!    In this situation the data at (r%from(iproc)%k(i),r%from(iproc)%l(i)) 
!    should come from (r%to(iproc)%k(i), r%to(iproc)%l(i)). 
!
! This do loop, in stella standard FFT situation, corresponds to:
!    to_here(it,ixxf)=from_here(ik,iyxf)
!
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

  end subroutine c_redist_22_inv_old_copy


  subroutine c_redist_22_inv_mpi_copy (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here


    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine c_redist_22_inv_mpi_copy


  subroutine c_redist_32 (r, from_here, to_here)

    use job_manage, only: time_message

    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here


    call c_redist_32_old_copy(r, from_here, to_here)
    call c_redist_32_mpi_copy(r, from_here, to_here)

  end subroutine c_redist_32


  subroutine c_redist_32_old_copy(r, from_here, to_here)

    use mp, only: iproc
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i

!CMR 
! In the stella standard FFT situation this routine maps 
!         g(ig, isgn, iglo) to xxf(it,ixxf) data type 
!         where it is kx (or x) index, ixxf is (y,ig,isgn,"les") 
!         and iglo is ("xyles") 
    do i = 1, r%from(iproc)%nn
!
! redistribute from local processor to local processor
! NB r%from(iproc)%nn is #elements sent by THIS processor to THIS processor
!    In this situation the data at (r%to(iproc)%k(i),r%to(iproc)%l(i)) 
!    should come from (r%from(iproc)%k(i),r%from(iproc)%l(i),r%from(iproc)%m(i)). 
!
! This do loop, in stella standard FFT situation, corresponds to:
!    to_here(it,ixxf)=from_here(ig,isgn,iglo)
!
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

  end subroutine c_redist_32_old_copy

  subroutine c_redist_32_mpi_copy(r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_redist_32_mpi_copy


  subroutine c_redist_32_inv (r, from_here, to_here)

    use job_manage, only: time_message

    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    call c_redist_32_inv_old_copy(r, from_here, to_here)
    call c_redist_32_inv_mpi_copy(r, from_here, to_here)

  end subroutine c_redist_32_inv

  subroutine c_redist_32_inv_old_copy(r, from_here, to_here)

    use mp, only: iproc

    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i

!CMR 
! In the stella standard FFT situation this routine maps 
!         xxf(it,ixxf) to g(ig, isgn, iglo) data type 
!         where it is kx (or x) index, ixxf is (y,ig,isgn,"les") 
!         and iglo is ("xyles") 

    do i = 1, r%to(iproc)%nn
!
! redistribute from local processor to local processor
! NB r%to(iproc)%nn is #elements sent by THIS processor to THIS processor
!    In this situation the data at (r%from(iproc)%k(i),r%from(iproc)%l(i),r%from(iproc)%m(i)) 
!    should come from (r%to(iproc)%k(i),r%to(iproc)%l(i)). 
!
! This do loop, in stella standard FFT situation, corresponds to:
!    to_here(ig,isgn,iglo)=from_here(it,ixxf)
!
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

  end subroutine c_redist_32_inv_old_copy

  subroutine c_redist_32_inv_mpi_copy(r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive

    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here


    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine c_redist_32_inv_mpi_copy



  subroutine c_redist_42 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i), &
                           r%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_redist_42

  subroutine c_redist_42_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i), &
               r%from(iproc)%n(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine c_redist_42_inv

  subroutine c_redist_23 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i),&
               r%to(iproc)%m(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_redist_23

  subroutine c_redist_34 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i),&
               r%to(iproc)%m(i),&
               r%to(iproc)%n(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i), &
                        r%to(ipfrom)%n(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i), &
                        r%to(ipfrom)%n(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_redist_34

  subroutine c_redist_34_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i), &
                           r%to(iproc)%m(i), &
                           r%to(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i), &
                                              r%to(ipto)%m(i), &
                                              r%to(ipto)%n(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i), &
                                              r%to(ipto)%m(i), &
                                              r%to(ipto)%n(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine c_redist_34_inv

  subroutine r_redist_34 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    real, dimension (r%from_low(1):, &
         r%from_low(2):, &
         r%from_low(3):), intent (in) :: from_here

    real, dimension (r%to_low(1):, &
         r%to_low(2):, &
         r%to_low(3):, &
         r%to_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i),&
               r%to(iproc)%m(i),&
               r%to(iproc)%n(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i), &
                        r%to(ipfrom)%n(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i), &
                        r%to(ipfrom)%n(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine r_redist_34

  subroutine r_redist_34_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    real, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):), intent (in) :: from_here

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i), &
                           r%to(iproc)%m(i), &
                           r%to(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i), &
                                              r%to(ipto)%m(i), &
                                              r%to(ipto)%n(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i), &
                                              r%to(ipto)%m(i), &
                                              r%to(ipto)%n(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine r_redist_34_inv

  subroutine c_redist_35 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):, &
                        r%to_low(5):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i),&
               r%to(iproc)%m(i),&
               r%to(iproc)%n(i),&
               r%to(iproc)%o(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i), &
                        r%to(ipfrom)%n(i), &
                        r%to(ipfrom)%o(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i), &
                        r%to(ipfrom)%n(i), &
                        r%to(ipfrom)%o(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_redist_35

  subroutine r_redist_35 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    real, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):, &
                        r%to_low(5):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i),&
               r%to(iproc)%m(i),&
               r%to(iproc)%n(i),&
               r%to(iproc)%o(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i), &
                        r%to(ipfrom)%n(i), &
                        r%to(ipfrom)%o(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i), &
                        r%to(ipfrom)%n(i), &
                        r%to(ipfrom)%o(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine r_redist_35

  subroutine c_redist_35_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):, &
                        r%to_low(5):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i), &
                           r%to(iproc)%m(i), &
                           r%to(iproc)%n(i), &
                           r%to(iproc)%o(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i), &
                                              r%to(ipto)%m(i), &
                                              r%to(ipto)%n(i), &
                                              r%to(ipto)%o(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i), &
                                              r%to(ipto)%m(i), &
                                              r%to(ipto)%n(i), &
                                              r%to(ipto)%o(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine c_redist_35_inv

  subroutine r_redist_35_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    real, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):, &
                        r%to_low(5):), intent (in) :: from_here

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i), &
                           r%to(iproc)%m(i), &
                           r%to(iproc)%n(i), &
                           r%to(iproc)%o(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i), &
                                              r%to(ipto)%m(i), &
                                              r%to(ipto)%n(i), &
                                              r%to(ipto)%o(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i), &
                                              r%to(ipto)%m(i), &
                                              r%to(ipto)%n(i), &
                                              r%to(ipto)%o(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine r_redist_35_inv

  subroutine r_redist_12 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    real, dimension (r%from_low(1):), intent (in) :: from_here

    real, dimension (r%to_low(1):,r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i), r%to(iproc)%l(i)) &
            = from_here(r%from(iproc)%k(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i)) 
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine r_redist_12

  subroutine r_redist_22 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    real, dimension (r%from_low(1):, &
                     r%from_low(2):), intent (in) :: from_here

    real, dimension (r%to_low(1):, &
                     r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine r_redist_22

  subroutine r_redist_22_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    real, dimension (r%to_low(1):, &
                     r%to_low(2):), intent (in) :: from_here

    real, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine r_redist_22_inv

  subroutine r_redist_32 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    real, dimension (r%to_low(1):, &
                     r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine r_redist_32

  subroutine r_redist_32_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    real, dimension (r%to_low(1):, &
                     r%to_low(2):), intent (in) :: from_here

    real, dimension (r%from_low(1):, &
                     r%from_low(2):, &
                     r%from_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine r_redist_32_inv

  subroutine r_redist_42 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in) :: from_here

    real, dimension (r%to_low(1):, &
                     r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i), &
                           r%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine r_redist_42

  subroutine r_redist_42_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    real, dimension (r%to_low(1):, &
                     r%to_low(2):), intent (in) :: from_here

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i), &
               r%from(iproc)%n(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do
 
  end subroutine r_redist_42_inv

  subroutine i_redist_12 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    integer, dimension (r%from_low(1):), intent (in) :: from_here

    integer, dimension (r%to_low(1):,r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i), r%to(iproc)%l(i)) &
            = from_here(r%from(iproc)%k(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i)) 
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine i_redist_12

  subroutine i_redist_22 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine i_redist_22

  subroutine i_redist_22_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%integer_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%integer_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine i_redist_22_inv

  subroutine i_redist_32 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine i_redist_32

  subroutine i_redist_32_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%integer_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%integer_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine i_redist_32_inv

  subroutine i_redist_42 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in) :: from_here

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i), &
                           r%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine i_redist_42

  subroutine i_redist_42_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i), &
               r%from(iproc)%n(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%integer_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%integer_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine i_redist_42_inv

  subroutine l_redist_12 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    logical, dimension (r%from_low(1):), intent (in) :: from_here

    logical, dimension (r%to_low(1):,r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i), r%to(iproc)%l(i)) &
            = from_here(r%from(iproc)%k(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i)) 
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine l_redist_12

  subroutine l_redist_22 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine l_redist_22

  subroutine l_redist_22_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%logical_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%logical_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine l_redist_22_inv

  subroutine l_redist_32 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine l_redist_32

  subroutine l_redist_32_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%logical_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%logical_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine l_redist_32_inv

  subroutine l_redist_42 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in) :: from_here

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i), &
                           r%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine l_redist_42

  subroutine l_redist_42_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i), &
               r%from(iproc)%n(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%logical_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%logical_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine l_redist_42_inv

! TT>
  subroutine c_redist_33 (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i),&
               r%to(iproc)%m(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_redist_33

  subroutine c_redist_33_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i), &
                           r%to(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i), &
                                              r%to(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i), &
                                              r%to(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine c_redist_33_inv
! <TT

  subroutine c_fill_2 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    complex, dimension (f%from_low(1):, &
                        f%from_low(2):), intent (in) :: from_here

    complex, dimension (f%to_low(1):, &
                        f%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%complex_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%complex_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%complex_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%complex_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%complex_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%complex_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_fill_2

  subroutine c_fill_3 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    complex, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here

    complex, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%complex_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%complex_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%complex_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%complex_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%complex_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%complex_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_fill_3

  subroutine c_fill_4 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    complex, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):, &
                        f%from_low(4):), intent (in) :: from_here

    complex, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):, &
                        f%to_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i), &
               f%to(iproc)%n(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i), &
                           f%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%complex_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%complex_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%complex_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%complex_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%complex_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%complex_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_fill_4

  subroutine r_fill_2 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    real, dimension (f%from_low(1):, &
                        f%from_low(2):), intent (in) :: from_here

    real, dimension (f%to_low(1):, &
                        f%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%real_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%real_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%real_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%real_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%real_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%real_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine r_fill_2

  subroutine r_fill_3 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    real, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here

    real, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%real_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%real_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%real_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%real_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%real_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%real_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine r_fill_3

  subroutine r_fill_4 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    real, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):, &
                        f%from_low(4):), intent (in) :: from_here

    real, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):, &
                        f%to_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i), &
               f%to(iproc)%n(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i), &
                           f%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%real_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%real_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%real_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%real_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%real_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%real_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine r_fill_4

  subroutine i_fill_2 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    integer, dimension (f%from_low(1):, &
                        f%from_low(2):), intent (in) :: from_here

    integer, dimension (f%to_low(1):, &
                        f%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%integer_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%integer_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%integer_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%integer_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%integer_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%integer_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine i_fill_2

  subroutine i_fill_3 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    integer, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here

    integer, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%integer_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%integer_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%integer_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%integer_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%integer_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%integer_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine i_fill_3

  subroutine i_fill_4 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    integer, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):, &
                        f%from_low(4):), intent (in) :: from_here

    integer, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):, &
                        f%to_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i), &
               f%to(iproc)%n(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i), &
                           f%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%integer_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%integer_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%integer_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%integer_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%integer_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%integer_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine i_fill_4

  subroutine l_fill_2 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    logical, dimension (f%from_low(1):, &
                        f%from_low(2):), intent (in) :: from_here

    logical, dimension (f%to_low(1):, &
                        f%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%logical_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%logical_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%logical_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%logical_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%logical_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%logical_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine l_fill_2

  subroutine l_fill_3 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    logical, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here

    logical, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%logical_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%logical_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%logical_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%logical_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%logical_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%logical_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine l_fill_3

  subroutine l_fill_4 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    logical, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):, &
                        f%from_low(4):), intent (in) :: from_here

    logical, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):, &
                        f%to_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i), &
               f%to(iproc)%n(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i), &
                           f%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%logical_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%logical_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%logical_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%logical_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%logical_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%logical_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine l_fill_4

  subroutine report_map_property (r)

    use mp, only: iproc, nproc, proc0, sum_reduce, max_reduce
    type (redist_type), intent (in) :: r
    type :: redist_prp
       integer :: local_max, local_total
       integer :: comm_max, comm_total
       integer :: elm_max, elm_total
    end type redist_prp
    type (redist_prp) :: prp
    integer :: ip, rank_from, rank_to
    integer, dimension(:), allocatable :: lbd_from, lbd_to

    prp%comm_max = 0
    prp%comm_total = 0
    prp%elm_total = 0

    do ip=0, nproc-1
       if (ip == iproc) then
          prp%local_total = r%to(ip)%nn
       else
          if (r%to(ip)%nn > 0) then
             prp%comm_total = prp%comm_total + 1
             prp%elm_total = prp%elm_total + r%to(ip)%nn
          end if
       end if
    end do
    prp%local_max = prp%local_total
    prp%comm_max = prp%comm_total
    prp%elm_max = prp%elm_total

    call max_reduce (prp%local_max, 0)
    call sum_reduce (prp%local_total, 0)
    call max_reduce (prp%comm_max, 0)
    call sum_reduce (prp%comm_total, 0)
    call max_reduce (prp%elm_max, 0)
    call sum_reduce (prp%elm_total, 0)

    if (proc0) then
       rank_from = 1
       if (associated(r%from(0)%l)) rank_from = 2
       if (associated(r%from(0)%m)) rank_from = 3
       if (associated(r%from(0)%n)) rank_from = 4
       rank_to = 1
       if (associated(r%to(0)%l)) rank_to = 2
       if (associated(r%to(0)%m)) rank_to = 3
       if (associated(r%to(0)%n)) rank_to = 4
       allocate (lbd_from(rank_from), lbd_to(rank_to))
       lbd_from = r%from_low(1:rank_from)
       lbd_to = r%to_low(1:rank_to)
       print '(a,i2,a,i2)', 'From rank', rank_from, ' to rank', rank_to
       print '(a,t20,4i10)', 'From lbound (proc0)', r%from_low(1:rank_from)
       print '(a,t20,4i10)', 'To lbound (proc0)',   r%to_low(1:rank_to)
       print '(a,t49,a,t64,a)', '--- Redistribution statistics ---', &
            'max', 'avg'
       print '(a,t40,i12,t55,f15.2)', 'Number of local move elements', &
            prp%local_max, real(prp%local_total)/real(nproc)
       print '(a,t40,i12,t60,f10.2)', &
            'Number of inter-processor communications', &
            prp%comm_max, real(prp%comm_total)/real(nproc)
       print '(a,t40,i12,t55,f15.2)', &
            'Number of inter-processor move elements', &
            prp%elm_max, real(prp%elm_total)/real(nproc)
       print *
    end if

  end subroutine report_map_property

  subroutine measure_gather_32 (map, gin, gout)

    use job_manage, only: time_message
    use mp, only: proc0

    type (redist_type), intent (in out) :: map
    complex, dimension (:,:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    if (proc0) call time_message(.false.,time_redist,' Redistribution')
    call c_redist_32 (map, gin, gout)
    if (proc0) call time_message(.false.,time_redist,' Redistribution')
    gather_count = gather_count + 1

  end subroutine measure_gather_32

  subroutine measure_scatter_23 (map, gin, gout)

    use job_manage, only: time_message
    use mp, only: proc0

    type (redist_type), intent (in out) :: map
    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:,:), intent (out) :: gout

    if (proc0) call time_message(.false.,time_redist,' Redistribution')
    call c_redist_32_inv (map, gin, gout)
    if (proc0) call time_message(.false.,time_redist,' Redistribution')
    scatter_count = scatter_count + 1

  end subroutine measure_scatter_23

  subroutine measure_gather_33 (map, gin, gout)

    use job_manage, only: time_message
    use mp, only: proc0

    type (redist_type), intent (in out) :: map
    complex, dimension (:,:,:), intent (in) :: gin
    complex, dimension (:,:,:), intent (out) :: gout

    if (proc0) call time_message(.false.,time_redist,' Redistribution')
    call c_redist_33 (map, gin, gout)
    if (proc0) call time_message(.false.,time_redist,' Redistribution')
    gather_count = gather_count + 1

  end subroutine measure_gather_33

  subroutine measure_scatter_33 (map, gin, gout)

    use job_manage, only: time_message
    use mp, only: proc0

    type (redist_type), intent (in out) :: map
    complex, dimension (:,:,:), intent (in) :: gin
    complex, dimension (:,:,:), intent (out) :: gout

    if (proc0) call time_message(.false.,time_redist,' Redistribution')
    call c_redist_33_inv (map, gin, gout)
    if (proc0) call time_message(.false.,time_redist,' Redistribution')
    scatter_count = scatter_count + 1

  end subroutine measure_scatter_33

  subroutine measure_gather_22 (map, gin, gout)

    use job_manage, only: time_message
    use mp, only: proc0

    type (redist_type), intent (in out) :: map
    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    if (proc0) call time_message(.false.,time_redist,' Redistribution')
    call c_redist_22 (map, gin, gout)
    if (proc0) call time_message(.false.,time_redist,' Redistribution')
    gather_count = gather_count + 1

  end subroutine measure_gather_22

  subroutine measure_scatter_22 (map, gin, gout)

    use job_manage, only: time_message
    use mp, only: proc0

    type (redist_type), intent (in out) :: map
    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    if (proc0) call time_message(.false.,time_redist,' Redistribution')
    call c_redist_22_inv (map, gin, gout)
    if (proc0) call time_message(.false.,time_redist,' Redistribution')
    scatter_count = scatter_count + 1

  end subroutine measure_scatter_22

end module redistribute
