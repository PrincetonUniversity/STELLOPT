!-----------------------------------------------------------------------
!     Module:        thrift_equil
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This module contains various variables which
!                    directly relate to equilibrium values.
!-----------------------------------------------------------------------
MODULE thrift_equil
    !-------------------------------------------------------------------
    !     Libraries
    !-------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE EZspline_obj
    !-------------------------------------------------------------------
    !     Module Variables
    !          bt0     B-Field on axis
    !          el0     Elongation
    !        iota0     Rotational Transform on axis
    !     iota_spl     Spline of iota profile
    !     phip_spl     Spline of dphi/drho
    !       vp_spl     Spline of dV/dphi
    !-------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: ns_eq
    REAL(rprec) :: bt0, el0, iota0, eq_beta, eq_Aminor, eq_Rmajor, &
                   eq_phiedge, eq_volume

    ! Spline helpers
    INTEGER :: bcs1(2)
    TYPE(EZspline1_r8) :: iota_spl, phip_spl, vp_spl, bu_spl, bv_spl, &
                            bsq_spl


END MODULE thrift_equil
