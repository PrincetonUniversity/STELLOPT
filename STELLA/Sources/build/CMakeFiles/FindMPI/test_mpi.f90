      program hello
      implicit none
      include 'mpif.h'
      integer ierror
      call MPI_INIT(ierror)
      call MPI_FINALIZE(ierror)
      end program
