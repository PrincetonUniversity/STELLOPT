!-----------------------------------------------------------------------
!     Module:        thrift_vars
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This module contains various variables used in the
!                    thrift code.
!-----------------------------------------------------------------------
MODULE thrift_vars
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    !-----------------------------------------------------------------------
    !     Module Variables
    !          lverb         Logical to control screen output
    !----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: ntimesteps, nrho
    REAL(rprec) :: tend


END MODULE thrift_vars
