!-----------------------------------------------------------------------
!     Module:        thrift_vars
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This module contains various variables used in the
!                    thrift code.
!-----------------------------------------------------------------------
MODULE thrift_vars
    !-------------------------------------------------------------------
    !     Libraries
    !-------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    !-------------------------------------------------------------------
    !     Module Variables
    !          lverb         Logical to control screen output
    !-------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: ntimesteps, nrho, win_thrift_j
    REAL(rprec) :: tend
    REAL(rprec), DIMENSION(:,:), POINTER :: THRIFT_J, THRIFT_S11, &
                 THRIFT_S12, THRIFT_S22, THRIFT_JBOOT


END MODULE thrift_vars
