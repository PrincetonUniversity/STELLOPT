!-------------------------------------------------------------------------------
! J. Breslau 8/18/17
MODULE windingsurface
  USE stel_kinds, ONLY : rprec
  IMPLICIT NONE

  TYPE vsurf
     REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: rctab, zstab
     INTEGER                                  :: mmax, nmax, nfp
  END TYPE vsurf

  TYPE(vsurf)                                 :: windsurf

CONTAINS
  SUBROUTINE read_winding_surface(filename, ierr)
    !USE stellopt_vars
    USE mpi_params
    USE mpi_inc
    IMPLICIT NONE
    INTRINSIC ALL, MAXVAL, ABS

    ! Arguments
    CHARACTER(256), INTENT(IN) :: filename
    INTEGER, INTENT(OUT)       :: ierr

    ! Local variables
    REAL(rprec), DIMENSION(:), ALLOCATABLE :: rmnc, zmns
    INTEGER, DIMENSION(:), ALLOCATABLE     :: mnum, nnum
    INTEGER  :: iunit, istat, count, mode
    LOGICAL  :: lexist

    ! Reset winding surface
    IF (ALLOCATED(windsurf%rctab)) DEALLOCATE(windsurf%rctab)
    IF (ALLOCATED(windsurf%zstab)) DEALLOCATE(windsurf%zstab)
    windsurf%mmax = -1;  windsurf%nmax = -1

    IF (myid.eq.master) THEN   ! Master process does I/O
       ! Open coil winding surface file
       lexist = .false.
       INQUIRE(FILE=TRIM(filename),EXIST=lexist)
       IF (lexist) THEN
          iunit=12;  istat=1
          OPEN(iunit, file=TRIM(filename), status='old', action='read', iostat=istat)
          IF (istat.eq.0) THEN  !Read the file
             WRITE(6,*)'Reading winding surface from file ',TRIM(filename),'.'

             count = 0
             DO ! Count lines
                READ(iunit,*,END=100)
                count = count + 1
             ENDDO

100          IF (count.gt.0) THEN
                WRITE(6,'(I6,A)')count,' modes in winding surface file.'

                ! Allocate temporary (sparse) space for values
                ALLOCATE(mnum(count), nnum(count), rmnc(count), zmns(count))

                ! Read values
                REWIND(iunit)
                DO mode=1,count
                   READ(iunit,*) mnum(mode), nnum(mode), rmnc(mode), zmns(mode)
                ENDDO !mode

                IF (ALL(mnum.ge.0)) THEN  ! Avoid error condition

                   ! Find max poloidal & toroidal mode numbers
                   windsurf%mmax = MAXVAL(mnum)
                   windsurf%nmax = MAXVAL(ABS(nnum))
110                FORMAT(I6,A,I5)
                   WRITE(6,110) -windsurf%nmax,' <= n <=',windsurf%nmax
                   WRITE(6,110) 0,' <= m <=',windsurf%mmax

                   ! Initialize the (dense) data structure
                   ALLOCATE(windsurf%rctab(-windsurf%nmax:windsurf%nmax, 0:windsurf%mmax))
                   ALLOCATE(windsurf%zstab(-windsurf%nmax:windsurf%nmax, 0:windsurf%mmax))
                   windsurf%rctab = 0d0;  windsurf%zstab = 0d0

                   !Load temporary data into permanent structure
                   DO mode=1,count
                      windsurf%rctab(nnum(mode),mnum(mode)) = rmnc(mode)
                      windsurf%zstab(nnum(mode),mnum(mode)) = zmns(mode)
                   ENDDO !mode
                ENDIF ! No negative m values

                ! Free up temporary space
                DEALLOCATE(mnum, nnum, rmnc, zmns)
             ENDIF !count > 0

             CLOSE(iunit)
          ENDIF !File opened successfully
       ENDIF !File exists
    ENDIF !Master process

!DEC$ IF DEFINED (MPI_OPT)
    ! Broadcast mode counts to all processors
    CALL MPI_BCAST(windsurf%nmax, 1, MPI_INTEGER, master, MPI_COMM_STEL, ierr)
    IF (ierr /= MPI_SUCCESS) RETURN
!DEC$ ELSE
    ierr = 0
!DEC$ ENDIF
    IF (windsurf%nmax.lt.0) THEN
       ierr = 1;  RETURN
    ENDIF
!DEC$ IF DEFINED (MPI_OPT)
    CALL MPI_BCAST(windsurf%mmax, 1, MPI_INTEGER, master, MPI_COMM_STEL, ierr)
    IF (ierr /= MPI_SUCCESS) RETURN

    ! Broadcast data to all processors
    IF (myid.ne.master) THEN
       ALLOCATE(windsurf%rctab(-windsurf%nmax:windsurf%nmax, 0:windsurf%mmax))
       ALLOCATE(windsurf%zstab(-windsurf%nmax:windsurf%nmax, 0:windsurf%mmax))
    ENDIF
    count = (2*windsurf%nmax + 1)*(windsurf%mmax + 1)
    CALL MPI_BCAST(windsurf%rctab, count, MPI_REAL8, master, MPI_COMM_STEL, ierr)
    IF (ierr /= MPI_SUCCESS) RETURN
    CALL MPI_BCAST(windsurf%zstab, count, MPI_REAL8, master, MPI_COMM_STEL, ierr)
    IF (ierr == MPI_SUCCESS) ierr = 0
!DEC$ ENDIF
  END SUBROUTINE read_winding_surface

!-----------------------------------------------------------------------
!     Subroutine:    stellopt_uv_to_xyz
!     Author:        J. Breslau (jbreslau@pppl.gov)
!     Date:          8/21/2017
!     Description:   This subroutine converts a u,v coordinate pair
!                    on a winding surface to Cartesian coordinates.
!     Inputs:
!                    Poloidal coordinate 0 < u < 1.
!                    Toroidal coordinate 0 < v < nfp.
!     Outputs:
!                    Cartesian values (x, y, z).
!-----------------------------------------------------------------------
  SUBROUTINE stellopt_uv_to_xyz(u, v, x, y, z)
    IMPLICIT NONE
    INTRINSIC COS, SIN

    ! Arguments
    REAL(rprec), INTENT(IN)  :: u, v
    REAL(rprec), INTENT(OUT) :: x, y, z

    ! Constants
    REAL(rprec), PARAMETER :: twopi = 6.283185307179586476925286766559D0

    ! Local variables
    REAL(rprec) :: theta, Nphi         ! Winding surface angular coordinates, in radians
    REAL(rprec) :: cmt, smt            ! cos(m theta), sin(m theta)
    REAL(rprec) :: cap, cam, sap, sam  ! cos(m theta +/- n N phi), sin(m theta +/- n N phi)
    REAL(rprec) :: R                   ! Cylindrical R coordinate
    REAL(rprec) :: Cu, Su, Cv, Sv, tmp
    INTEGER     :: m, n                ! Poloidal, toroidal mode numbers

    ! Initialize N phi trig. functions
    Nphi = twopi*v
    Cv = COS(Nphi);  Sv = SIN(Nphi)

    ! Initialize m theta trig. fcns, recurrence relation
    theta = twopi*u
    Cu = COS(theta);  Su = SIN(theta)
    cmt = 1.0;  smt = 0.0

    ! Outer loop over poloidal mode numbers
    R = 0.0;  z = 0.0
    DO m=0,windsurf%mmax
       ! Initialize +/- N phi recurrence relations
       cam = cmt;  cap = cmt;  sam = smt;  sap = smt

       ! n=0 mode contributions
       R = R + windsurf%rctab(0,m)*cap
       z = z + windsurf%zstab(0,m)*sap

       ! Inner loop over nonzero toroidal mode numbers
       DO n=1,windsurf%nmax
          ! Increment/decrement sin(arg), cos(arg)
          tmp = cap*Cv - sap*Sv;  sap = sap*Cv + cap*Sv;  cap = tmp
          tmp = cam*Cv + sam*Sv;  sam = sam*Cv - cam*Sv;  cam = tmp

          ! Add terms to series
          R = R + windsurf%rctab(n,m)*cap + windsurf%rctab(-n,m)*cam
          z = z + windsurf%zstab(n,m)*sap + windsurf%zstab(-n,m)*sam
       END DO !n

       ! Increment sin(m theta), cos(m theta)
       tmp = cmt*Cu - smt*Su;  smt = smt*Cu + cmt*Su;  cmt = tmp
    END DO !m

    ! Convert from cylindrical to Cartesian coordinates
    x = R*COS(Nphi/windsurf%nfp);  y = R*SIN(Nphi/windsurf%nfp)
  END SUBROUTINE stellopt_uv_to_xyz

!-------------------------------------------------------------------------------
!     Subroutine:    stellopt_uv_to_xyz_prime
!     Author:        J. Breslau (jbreslau@pppl.gov)
!     Date:          8/24/2017
!     Description:   This subroutine takes a u,v coordinate pair
!                    on a winding surface and returns the first and second
!                    derivatives of Cartesian coordinates (x,y,z) with respect
!                    to u and v at that location.
!     Inputs:
!                    Poloidal coordinate 0 < u < 1.
!                    Toroidal coordinate 0 < v < nfp.
!
!     Output:        xyzprime(3,5), where
!                      xyzprime(:,1) = d(xyz)/du
!                      xyzprime(:,2) = d(xyz)/dv
!                      xyzprime(:,3) = d2(xyz)/du2
!                      xyzprime(:,4) = d2(xyz)/dudv
!                      xyzprime(:,5) = d2(xyz)/dv2
!-------------------------------------------------------------------------------
  SUBROUTINE stellopt_uv_to_xyz_prime(u, v, xyzprime)
    IMPLICIT NONE
    INTRINSIC COS, SIN

    ! Arguments
    REAL(rprec),                 INTENT(IN)  :: u, v
    REAL(rprec), DIMENSION(3,5), INTENT(OUT) :: xyzprime

    ! Constants
    REAL(rprec), PARAMETER :: twopi = 6.283185307179586476925286766559D0

    ! Local variables
    REAL(rprec) :: theta, Nphi         ! Winding surface angular coordinates, in radians
    REAL(rprec) :: cmt, smt            ! cos(m theta), sin(m theta)
    REAL(rprec) :: cap, cam, sap, sam  ! cos(m theta +/- n N phi), sin(m theta +/- n N phi)
    REAL(rprec) :: R                   ! Cylindrical R coordinate
    REAL(rprec) :: Rprime(5)           ! 1st, 2nd Derivatives of R w/r/t u,v
    REAL(rprec) :: Cu, Su, Cv, Sv, tmp, dphidv
    INTEGER     :: m, n                ! Poloidal, toroidal mode numbers

    ! First get derivatives of cylindrical coordinates R,z
    ! Initialize N phi trig. functions
    Nphi = twopi*v
    Cv = COS(Nphi);  Sv = SIN(Nphi)

    ! Initialize m theta trig. fcns, recurrence relation
    theta = twopi*u
    Cu = COS(theta);  Su = SIN(theta)
    cmt = 1.0;  smt = 0.0

    ! Outer loop over poloidal mode numbers
    R = 0.0;  Rprime = 0.0;  xyzprime(3,:) = 0.0
    DO m=0,windsurf%mmax
       ! Initialize +/- N phi recurrence relations
       cam = cmt;  cap = cmt;  sam = smt;  sap = smt

       ! n=0 mode contributions
       R         = R         +     windsurf%rctab(0,m)*cap         ! R
       Rprime(1) = Rprime(1) +   m*windsurf%rctab(0,m)*sap         ! dR/du
       Rprime(3) = Rprime(3) + m*m*windsurf%rctab(0,m)*cap         ! d2R/du^2
       xyzprime(3,1) = xyzprime(3,1) +   m*windsurf%zstab(0,m)*cap ! dz/du
       xyzprime(3,3) = xyzprime(3,3) + m*m*windsurf%zstab(0,m)*sap ! d2z/du2

       ! Inner loop over nonzero toroidal mode numbers
       DO n=1,windsurf%nmax
          ! Increment/decrement sin(arg), cos(arg)
          tmp = cap*Cv - sap*Sv;  sap = sap*Cv + cap*Sv;  cap = tmp
          tmp = cam*Cv + sam*Sv;  sam = sam*Cv - cam*Sv;  cam = tmp

          ! Add terms to series
          R         = R         +      windsurf%rctab(n,m)*cap + windsurf%rctab(-n,m)*cam   ! R
          Rprime(1) = Rprime(1) +   m*(windsurf%rctab(n,m)*sap + windsurf%rctab(-n,m)*sam)  ! dR/du
          Rprime(2) = Rprime(2) +   n*(windsurf%rctab(n,m)*sap - windsurf%rctab(-n,m)*sam)  ! dR/dv
          Rprime(3) = Rprime(3) + m*m*(windsurf%rctab(n,m)*cap + windsurf%rctab(-n,m)*cam)  ! d2R/du2
          Rprime(4) = Rprime(4) + m*n*(windsurf%rctab(n,m)*cap - windsurf%rctab(-n,m)*cam)  ! d2R/dudv
          Rprime(5) = Rprime(5) + n*n*(windsurf%rctab(n,m)*cap + windsurf%rctab(-n,m)*cam)  ! d2R/dv2

          xyzprime(3,1) = xyzprime(3,1) +   m*(windsurf%zstab(n,m)*cap + windsurf%zstab(-n,m)*cam) ! dz/du
          xyzprime(3,2) = xyzprime(3,2) +   n*(windsurf%zstab(n,m)*cap - windsurf%zstab(-n,m)*cam) ! dz/dv
          xyzprime(3,3) = xyzprime(3,3) + m*m*(windsurf%zstab(n,m)*sap + windsurf%zstab(-n,m)*sam) ! d2z/du2
          xyzprime(3,4) = xyzprime(3,4) + m*n*(windsurf%zstab(n,m)*sap - windsurf%zstab(-n,m)*sam) ! d2z/dudv
          xyzprime(3,5) = xyzprime(3,5) + n*n*(windsurf%zstab(n,m)*sap + windsurf%zstab(-n,m)*sam) ! d2z/dv2
       END DO !n

       ! Increment sin(m theta), cos(m theta)
       tmp = cmt*Cu - smt*Su;  smt = smt*Cu + cmt*Su;  cmt = tmp
    END DO !m

    ! Scale derivatives by appropriate multiples of two*pi
    Rprime(1) =      -twopi * Rprime(1)
    Rprime(2) =      -twopi * Rprime(2)
    Rprime(3) = -(twopi**2) * Rprime(3)
    Rprime(4) = -(twopi**2) * Rprime(4)
    Rprime(5) = -(twopi**2) * Rprime(5)

    xyzprime(3,1) =       twopi * xyzprime(3,1)
    xyzprime(3,2) =       twopi * xyzprime(3,2)
    xyzprime(3,3) = -(twopi**2) * xyzprime(3,3)
    xyzprime(3,4) = -(twopi**2) * xyzprime(3,4)
    xyzprime(3,5) = -(twopi**2) * xyzprime(3,5)

    ! Now convert from cylindrical R to Cartesian (x,y)
    dphidv = twopi/REAL(windsurf%nfp)
    Cv = COS(dphidv*v);  Sv = SIN(dphidv*v)

    xyzprime(1,1) = Rprime(1)*Cv                                            !dx/du
    xyzprime(2,1) = Rprime(1)*Sv                                            !dy/du

    xyzprime(1,2) = Rprime(2)*Cv - R*Sv*dphidv                              !dx/dv
    xyzprime(2,2) = Rprime(2)*Sv + R*Cv*dphidv                              !dy/dv

    xyzprime(1,3) = Rprime(3)*Cv                                            !d2x/du2
    xyzprime(2,3) = Rprime(3)*Sv                                            !d2y/du2

    xyzprime(1,4) = Rprime(4)*Cv - Rprime(1)*Sv*dphidv                      !d2x/dudv
    xyzprime(2,4) = Rprime(4)*Sv + Rprime(1)*Cv*dphidv                      !d2y/dudv

    xyzprime(1,5) = Rprime(5)*Cv - (2.0*Rprime(2)*Sv + R*Cv*dphidv)*dphidv  !d2x/dv2
    xyzprime(2,5) = Rprime(5)*Sv + (2.0*Rprime(2)*Cv - R*Sv*dphidv)*dphidv  !d2y/dv2
  END SUBROUTINE stellopt_uv_to_xyz_prime
END MODULE windingsurface
