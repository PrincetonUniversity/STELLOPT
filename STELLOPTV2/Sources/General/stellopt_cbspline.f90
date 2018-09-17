! Evaluate 1D basis splines and their first two derivatives over their domain
! given lists of knots and control points.
MODULE stellopt_cbspline
  IMPLICIT NONE
  INTEGER, PARAMETER :: CBR8 = selected_real_kind(15,100)

  TYPE cbspline
     REAL(CBR8), DIMENSION(:), ALLOCATABLE :: knot, coef
     INTEGER :: ncoefs, mu_max
  END type cbspline
CONTAINS
!---------------------------------------------------------------
  SUBROUTINE cbspline_init(cs, n, ierr)
    IMPLICIT NONE
    TYPE(cbspline), INTENT(OUT) :: cs
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: ierr

    ierr = 0
    IF (n.LT.1) THEN
       ierr = 1
       RETURN
    ENDIF

    cs%ncoefs = n;  cs%mu_max = n

    IF (ALLOCATED(cs%knot)) DEALLOCATE(cs%knot)
    ALLOCATE(cs%knot(0:n+3))
    IF (ALLOCATED(cs%coef)) DEALLOCATE(cs%coef)
    ALLOCATE(cs%coef(0:n-1))

    IF (.NOT.ALLOCATED(cs%knot) .OR. .NOT.ALLOCATED(cs%coef)) ierr = 2
  END SUBROUTINE cbspline_init
!---------------------------------------------------------------
  SUBROUTINE cbspline_delete(cs)
    IMPLICIT NONE
    TYPE(cbspline), INTENT(INOUT) :: cs
    cs%ncoefs = 0;  cs%mu_max = 0
    IF (ALLOCATED(cs%knot)) DEALLOCATE(cs%knot)
    IF (ALLOCATED(cs%coef)) DEALLOCATE(cs%coef)
  END SUBROUTINE cbspline_delete
!---------------------------------------------------------------
  SUBROUTINE cbspline_setup(cs, n, kval, cval, ierr)
    IMPLICIT NONE
    TYPE(cbspline), INTENT(INOUT) :: cs
    INTEGER, INTENT(IN) :: n
    REAL(CBR8), DIMENSION(n+4), INTENT(IN) :: kval
    REAL(CBR8), DIMENSION(n), INTENT(IN) :: cval
    INTEGER, INTENT(OUT) :: ierr

    ierr = 0
    IF (n.ne.cs%ncoefs) THEN
       ierr = 1
       RETURN
    ENDIF

    cs%knot(0:n+3) = kval(1:n+4)
    cs%coef(0:n-1) = cval(1:n)
  END SUBROUTINE cbspline_setup
!---------------------------------------------------------------
  SUBROUTINE cbspline_eval(cs, s, f, ierr)
    IMPLICIT NONE
    TYPE(cbspline), INTENT(IN) :: cs
    REAL(CBR8), INTENT(IN) :: s
    REAL(CBR8), INTENT(OUT) :: f
    INTEGER, INTENT(OUT) :: ierr

    REAL(CBR8) :: kms, smk, cm4, cm1, c2_0, c2_1, c2_2, c1_0, c1_1
    INTEGER mu

    ierr = 0
    IF (cs%ncoefs.LT.1) THEN
       ierr = 1
       RETURN
    ENDIF

    ! Test whether in-bounds
    IF ((s.LT.cs%knot(2)).OR.(s.GT.cs%knot(cs%mu_max))) THEN
       ierr = 2
       RETURN
    ENDIF

    ! Find the knot interval containing the point
    DO mu=3,cs%mu_max-1
       IF (cs%knot(mu).GT.s) EXIT
    ENDDO

    ! Construct the spline at s using three matrix multiplies.
    kms = cs%knot(mu) - s;  smk = s - cs%knot(mu - 1)
    IF (mu.GT.3) THEN
       cm4 = cs%coef(mu-4)
    ELSE
       cm4 = cs%coef(cs%ncoefs - 2)
    ENDIF
    IF (mu.LT.cs%ncoefs+1) THEN
       cm1 = cs%coef(mu-1)
    ELSE
       cm1 = cs%coef(1)
    ENDIF

    c2_0 = (kms*cm4 + (s - cs%knot(mu-3))*cs%coef(mu-3))/&
         (cs%knot(mu) - cs%knot(mu-3))
    c2_1 = ((cs%knot(mu+1) - s)*cs%coef(mu-3) + &
         (s - cs%knot(mu-2))*cs%coef(mu-2))/(cs%knot(mu+1) - cs%knot(mu-2))
    c2_2 = ((cs%knot(mu+2) - s)*cs%coef(mu-2) + &
         smk*cm1)/(cs%knot(mu+2) - cs%knot(mu-1))

    c1_0 = (kms*c2_0 + &
         (s - cs%knot(mu-2))*c2_1)/(cs%knot(mu) - cs%knot(mu-2))
    c1_1 = ((cs%knot(mu+1) - s)*c2_1 + &
         smk*c2_2)/(cs%knot(mu+1) - cs%knot(mu-1))

    f = (kms*c1_0 + smk*c1_1)/(cs%knot(mu) - cs%knot(mu-1))
  END SUBROUTINE cbspline_eval
!---------------------------------------------------------------
  SUBROUTINE cbspline_derivs(cs, s, f, fp, fpp, ierr)
    IMPLICIT NONE
    TYPE(cbspline), INTENT(IN) :: cs
    REAL(CBR8), INTENT(IN) :: s
    REAL(CBR8), INTENT(OUT) :: f, fp, fpp
    INTEGER, INTENT(OUT) :: ierr

    REAL(CBR8) :: kms, smk, cm4, cm1, c2_0, c2_1, c2_2, c1_0, c1_1
    REAL(CBR8) :: d2_0, d2_1, d2_2, d1_0, d1_1, e1_0, e1_1
    INTEGER mu

    ierr = 0
    IF (cs%ncoefs.LT.1) THEN
       ierr = 1
       RETURN
    ENDIF

    ! Test whether in-bounds
    IF ((s.LT.cs%knot(2)).OR.(s.GT.cs%knot(cs%mu_max))) THEN
       ierr = 2
       RETURN
    ENDIF

    ! Find the knot interval containing the point
    DO mu=3,cs%mu_max-1
       IF (cs%knot(mu).GT.s) EXIT
    ENDDO

    ! Construct the spline & derivs at s using three matrix multiplies.
    kms = cs%knot(mu) - s;  smk = s - cs%knot(mu - 1)
    IF (mu.GT.3) THEN
       cm4 = cs%coef(mu-4)
    ELSE
       cm4 = cs%coef(cs%ncoefs - 2)
    ENDIF
    IF (mu.LT.cs%ncoefs+1) THEN
       cm1 = cs%coef(mu-1)
    ELSE
       cm1 = cs%coef(1)
    ENDIF

    d2_0 = (cs%coef(mu-3) - cm4)/(cs%knot(mu) - cs%knot(mu-3))
    d2_1 = (cs%coef(mu-2) - cs%coef(mu-3))/(cs%knot(mu+1) - cs%knot(mu-2))
    d2_2 = (cm1 - cs%coef(mu-2))/(cs%knot(mu+2) - cs%knot(mu-1))

    c2_0 = (kms*cm4 + (s - cs%knot(mu-3))*cs%coef(mu-3))/&
         (cs%knot(mu) - cs%knot(mu-3))
    c2_1 = ((cs%knot(mu+1) - s)*cs%coef(mu-3) + &
         (s - cs%knot(mu-2))*cs%coef(mu-2))/(cs%knot(mu+1) - cs%knot(mu-2))
    c2_2 = ((cs%knot(mu+2) - s)*cs%coef(mu-2) + &
         smk*cm1)/(cs%knot(mu+2) - cs%knot(mu-1))

    e1_0 = 2.0*(d2_1 - d2_0)/(cs%knot(mu) - cs%knot(mu-2))
    e1_1 = 2.0*(d2_2 - d2_1)/(cs%knot(mu+1) - cs%knot(mu-1))

    d1_0 = ((kms*d2_0 - c2_0) + &
         ((s - cs%knot(mu-2))*d2_1 + c2_1))/(cs%knot(mu) - cs%knot(mu-2))
    d1_1 = (((cs%knot(mu+1) - s)*d2_1 - c2_1) + &
         (smk*d2_2 + c2_2))/(cs%knot(mu+1) - cs%knot(mu-1))

    c1_0 = (kms*c2_0 + &
         (s - cs%knot(mu-2))*c2_1)/(cs%knot(mu) - cs%knot(mu-2))
    c1_1 = ((cs%knot(mu+1) - s)*c2_1 + &
         smk*c2_2)/(cs%knot(mu+1) - cs%knot(mu-1))

    f = (kms*c1_0 + smk*c1_1)/(cs%knot(mu) - cs%knot(mu-1))
    fp = ((kms*d1_0 - c1_0) + &
         (smk*d1_1 + c1_1))/(cs%knot(mu) - cs%knot(mu-1))
    fpp = ((kms*e1_0 - 2.0*d1_0) + &
         (smk*e1_1 + 2.0*d1_1))/(cs%knot(mu) - cs%knot(mu-1));
  END SUBROUTINE cbspline_derivs
END MODULE stellopt_cbspline
