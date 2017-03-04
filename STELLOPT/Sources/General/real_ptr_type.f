      MODULE real_ptr_type
      USE stel_kinds

      TYPE real_ptr
            REAL(rprec), DIMENSION(:), POINTER :: x
      END TYPE real_ptr

      END MODULE real_ptr_type
