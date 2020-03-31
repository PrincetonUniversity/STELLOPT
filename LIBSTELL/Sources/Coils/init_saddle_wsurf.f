      SUBROUTINE init_saddle_wsurf (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE saddle_surface
      USE mpi_params
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nb, mb, ik
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: rbc_in, zbs_in
      REAL(rprec) :: delta
      INTEGER, INTENT(out) :: nvariables
      REAL(rprec) :: xvariables(*)
      EXTERNAL unique_boundary, convert_boundary,
     1   unique_boundary_PG, convert_boundary_PG
!-----------------------------------------------
      IF (nopt_wsurf .eq. -1) THEN
         xvariables(1:numsurf_sad) = rmn_sad(:numsurf_sad)
         xvariables(numsurf_sad+1:2*numsurf_sad) = zmn_sad(:numsurf_sad)
         nvariables = 2*numsurf_sad
         RETURN
      END IF

      irm0_bdy = 0;      izm0_bdy = 0;      irho_bdy = 0;

      ntor_opt = MAX(MAXVAL(n_sad(1:numsurf_sad)), 
     1           ABS(MINVAL(n_sad(1:numsurf_sad))))
      mpol_opt = MAX(MAXVAL(m_sad(1:numsurf_sad)), 
     1           ABS(MINVAL(m_sad(1:numsurf_sad))))
      IF (nopt_wsurf .eq. 1) mpol_opt = mpol_opt + 1
      ik = (2*ntor_opt+1)*(mpol_opt+1)

      ALLOCATE (nbrho_opt(ik), mbrho_opt(ik),
     1  rbc(-ntor_opt:ntor_opt,0:mpol_opt),
     2  zbs(-ntor_opt:ntor_opt,0:mpol_opt),
     3  rbc_in(-ntor_opt:ntor_opt,0:mpol_opt),
     4  zbs_in(-ntor_opt:ntor_opt,0:mpol_opt),
     5  rhobc(-ntor_opt:ntor_opt,0:mpol_opt),
     6  nrz0_opt(2*(ntor_opt+1)),
     7  delta_mn(-ntor_opt:ntor_opt,-mpol_opt:mpol_opt))

      rbc = 0
      zbs = 0
      nvariables = 0

      DO ik = 1, numsurf_sad
         mb = m_sad(ik)

         IF (mb .eq. 0) THEN
            nb = n_sad(ik)
            IF (nb .ge. 0) THEN
               rbc(nb, mb) = rmn_sad(ik)
               zbs(nb, mb) = -zmn_sad(ik)
            ELSE
               nb = -nb
               rbc(nb, mb) = rmn_sad(ik)
               zbs(nb, mb) = zmn_sad(ik)
            END IF                              ! mb=0

         ELSE
            nb = -n_sad(ik)
            rbc(nb, mb) = rmn_sad(ik)
            zbs(nb, mb) = zmn_sad(ik)
         END IF

      END DO

      IF (nopt_wsurf .eq. 0) THEN

         rbc_in = rbc
         zbs_in = zbs                                                    !copy and store m=0 components

!
!        check IF conversion was made or IF original bdy already in proper form
!        IF NOPT_BOUNDARY=0, USE HIRSHMAN/BRESLAU REPRESENTATION
!                        =1, USE PG (PAUL GARABEDIAN) DELTA_MN REPRESENTATION
!
         CALL convert_boundary(rbc, zbs, rhobc, mpol_opt, ntor_opt)

         CALL unique_boundary(rbc_in, zbs_in, rhobc, mpol_opt, 
     1                        ntor_opt, mpol_opt, ntor_opt, mpol_opt)
         delta = SUM((rbc - rbc_in)**2)/rbc(0,1)**2
     1         + SUM((zbs - zbs_in)**2)/zbs(0,1)**2

         IF (delta.gt.1.e-8_dp .and. myid.eq.master) 
     1      WRITE(*,10) 100*(one-delta)

 10   FORMAT(' Input boundary representation was converted!',/,
     1       ' Reliability of conversion = ',f7.2,'%')

         DO nb = -ntor_opt, ntor_opt
            IF (rbc(nb,0).ne.zero .or. zbs(nb,0).ne.zero) THEN
               nvariables = nvariables + 1
               n_sad(nvariables) = nb
               m_sad(nvariables) = 0
               irm0_bdy = irm0_bdy + 1
               nrz0_opt(irm0_bdy) = nb
               xvariables(nvariables) = rbc(nb,0)
            END IF
         END DO

         DO ik = 1, irm0_bdy
            nb = n_sad(ik)
            IF (nb .ne. 0) THEN
               izm0_bdy = izm0_bdy + 1
               nvariables = nvariables + 1
               nrz0_opt(nvariables) = nb
               xvariables(nvariables) = zbs(nb,0)
            END IF
         END DO

         ik = irm0_bdy

         DO mb = 0, mpol_opt
            DO nb = -ntor_opt, ntor_opt
               IF (mb .ne. 0) THEN
                  ik = ik+1
                  n_sad(ik) = nb
                  m_sad(ik) = mb
               END IF
               IF (rhobc(nb,mb) .ne. zero
     1             .and. (mb .ne. 0 .or. nb .ge. 0)) THEN
                  nvariables = nvariables + 1
                  irho_bdy = irho_bdy + 1
                  nbrho_opt(irho_bdy) = nb
                  mbrho_opt(irho_bdy) = mb
                  xvariables(nvariables) = rhobc(nb,mb)
               END IF
            END DO
         END DO

         numsurf_sad = ik

         IF (numsurf_sad .gt. nsurf) STOP ' NUMSURF_SAD > NSURF '

         END IF

         IF (nopt_wsurf.eq.1) THEN
           CALL convert_boundary_PG(rbc,zbs,delta_mn,mpol_opt,ntor_opt)

           DO nb = -ntor_opt, ntor_opt
             DO mb = -mpol_opt, mpol_opt
               IF (delta_mn(nb,mb) .ne. zero
     1            .and. .not.(nb .eq.0 .and. mb .eq. 0)) THEN
                  irho_bdy = irho_bdy + 1
                  nvariables = nvariables + 1
                  nbrho_opt(irho_bdy) = nb
                  mbrho_opt(irho_bdy) = mb
                  xvariables(nvariables) = delta_mn(nb,mb)
               END IF
             END DO
           END DO

         END IF

      DEALLOCATE (rbc_in, zbs_in)

      END SUBROUTINE init_saddle_wsurf
