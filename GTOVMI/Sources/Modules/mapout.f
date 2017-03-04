
       MODULE mapout
        USE precision
        TYPE :: RSZS
!	All MKS units !!!!!!!!!!!!
         INTEGER :: npsi, nthet
         REAL(rprec) :: raxis, zaxis, eaxis, qaxis, fraction_bndry
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: rcentr
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: zcentr
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: aminor
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: elong
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: tnedni	! indentation
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: triang
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: square	! squareness
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: psival
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: vplas
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: area
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: curint
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: curavg	! <Jtor>
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: curintp	! d(Itor/ds)
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: phi
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: qsi
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: fpol
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: pressure
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: tflx
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: vprime
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: ffprime
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: pprime
         REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: rs
         REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: zs
         REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: arcsur
         REAL(rprec),  DIMENSION(0:19,0:0) :: rbc, rbs, zbc, zbs
        END TYPE RSZS
      END MODULE mapout
