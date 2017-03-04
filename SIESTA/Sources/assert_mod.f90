      MODULE assert_mod
        USE stel_kinds
        interface assert
        module  procedure assert_int, assert_real, assert_string
        end interface assert

        CONTAINS

        subroutine assert_int(lcond,msg,ival)
        implicit none
        logical,intent(in) :: lcond
        character*(*),intent(in) :: msg
        integer,intent(in) :: ival

        if (.not.lcond) then
          write(*,*) msg,ival
          stop '** assertion error ** '
        endif

        return
        end subroutine assert_int


        subroutine assert_real(lcond,msg,ival)
        implicit none
        logical,intent(in) :: lcond
        character*(*),intent(in) :: msg
        real(rprec),intent(in) :: ival

        if (.not.lcond) then
          write(*,*) msg,ival
          stop '** assertion error ** '
        endif

        return
        end subroutine assert_real


        subroutine assert_string(lcond,msg,ival)
        implicit none
        logical,intent(in) :: lcond
        character*(*),intent(in) :: msg
        character*(*),intent(in) :: ival

        if (.not.lcond) then
          write(*,*) msg,ival
          stop '** assertion error ** '
        endif

        return
        end subroutine assert_string

      END MODULE assert_mod
