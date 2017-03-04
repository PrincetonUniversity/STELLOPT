        MODULE gfile
         USE precision
         TYPE :: GEQDSK
         INTEGER :: shot
         INTEGER :: time
         INTEGER :: error
         CHARACTER(len=63) :: source
         CHARACTER(len=10), DIMENSION(6) :: ecase
         INTEGER :: mw
         INTEGER :: mh
         REAL(rprec)  :: xdim
         REAL(rprec) :: zdim
         REAL(rprec) :: rzero
         REAL(rprec) :: rgrid1
         REAL(rprec)  :: zmid
         REAL(rprec) :: rmaxis
         REAL(rprec) :: zmaxis
         REAL(rprec) :: ssimag
         REAL(rprec)  :: ssibry
         REAL(rprec) :: bcentr
         REAL(rprec)  :: cpasma
         REAL(rprec) ,DIMENSION(:), ALLOCATABLE :: fpol
         REAL(rprec) ,DIMENSION(:), ALLOCATABLE :: pres
         REAL(rprec) ,DIMENSION(:), ALLOCATABLE  :: ffprim
         REAL(rprec) ,DIMENSION(:), ALLOCATABLE :: pprime
         REAL(rprec) ,DIMENSION(:,:), ALLOCATABLE :: psirz
         REAL(rprec) ,DIMENSION(:), ALLOCATABLE  :: qpsi
         INTeger(rprec) :: nbdry
         INTeger(rprec) :: limitr
         REAL(rprec) ,DIMENSION(:), ALLOCATABLE  :: rbdry
         REAL(rprec) ,DIMENSION(:), ALLOCATABLE  :: zbdry
         REAL(rprec) ,DIMENSION(:), ALLOCATABLE :: xlim
         REAL(rprec) ,DIMENSION(:), ALLOCATABLE :: ylim
         REAL(rprec) ,DIMENSION(:), ALLOCATABLE :: R
         REAL(rprec) ,DIMENSION(:), ALLOCATABLE :: Z
         REAL(rprec) ,DIMENSION(:), ALLOCATABLE :: rhovn
         REAL(rprec) ,DIMENSION(:), ALLOCATABLE :: epoten
        END TYPE GEQDSK
      END MODULE gfile
