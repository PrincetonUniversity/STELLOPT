      MODULE coils
      USE coiltypes

! This module declares an allocatable array of coils for each type of coil.
! The actual allocation is done in the initialization routine.

         TYPE(modularcoil), DIMENSION(:), ALLOCATABLE :: modular
         TYPE(saddlecoil),  DIMENSION(:), ALLOCATABLE :: saddle
         TYPE(vfcoil),      DIMENSION(:), ALLOCATABLE :: vertical

      END MODULE coils
