      SUBROUTINE load_saddle_wsurf (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE saddle_surface
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nvariables, mb, nb, ik, j1
      REAL(rprec) :: xvariables(*)
      EXTERNAL unique_boundary, unique_boundary_PG
!-----------------------------------------------

      IF (nopt_wsurf .eq. -1) THEN
         rmn_sad(:numsurf_sad) = xvariables(1:numsurf_sad)
         zmn_sad(:numsurf_sad) = xvariables(numsurf_sad+1:2*numsurf_sad)
         nvariables = 2*numsurf_sad

      ELSE IF (nopt_wsurf .eq. 0) THEN

         rbc = 0;   zbs = 0

         DO ik = 1, irm0_bdy
            rbc(nrz0_opt(ik),0) = xvariables(ik)
         END DO
         DO ik = 1, izm0_bdy
            zbs(nrz0_opt(ik+irm0_bdy),0) = xvariables(ik+irm0_bdy)
         END DO
         DO ik = 1, irho_bdy
            j1 = ik + irm0_bdy + izm0_bdy 
            rhobc(nbrho_opt(ik),mbrho_opt(ik)) = xvariables(j1)
         END DO

         nvariables = irm0_bdy + izm0_bdy + irho_bdy

         CALL unique_boundary(rbc, zbs, rhobc, mpol_opt, ntor_opt, 
     1                               mpol_opt, ntor_opt, mpol_opt)

!        CONVERT FROM VMEC TO NESCOIL FORM FOR WINDING SURFACE

         rmn_sad = 0
         zmn_sad = 0

         DO ik = 1, numsurf_sad
            mb = m_sad(ik)

            IF (mb .eq. 0) THEN
               nb = n_sad(ik)
               IF (nb .ge. 0) THEN
                  rmn_sad(ik) =  rbc(nb, mb)
                  zmn_sad(ik) = -zbs(nb, mb)
               ELSE
                  nb = -nb
                  rmn_sad(ik) = rbc(nb, mb)
                  zmn_sad(ik) = zbs(nb, mb)
               END IF                              ! mb=0

            ELSE
               nb = -n_sad(ik)
               rmn_sad(ik) = rbc(nb, mb)
               zmn_sad(ik) = zbs(nb, mb)
            END IF
         END DO


      ELSE IF (nopt_wsurf .eq. 1) THEN
       
         DO ik = 1, irho_bdy
            j1 = ik
            delta_mn(nbrho_opt(ik),mbrho_opt(ik)) = xvariables(j1)
         END DO

         CALL unique_boundary_PG(rbc, zbs, delta_mn, ntor_opt, mpol_opt,
     1                           mpol_opt, ntor_opt)

         nvariables = 0
         DO mb = 0, mpol_opt
            DO nb = -ntor_opt, ntor_opt
               IF( ABS(rbc(nb,mb)) .gt. zero .or.
     1            ABS(zbs(nb,mb)) .gt. zero ) THEN
                  nvariables = nvariables + 1
                  m_sad(nvariables) = mb
                  IF (mb .eq. 0) THEN
                     n_sad(nvariables) = nb
                     rmn_sad(nvariables) = rbc(nb,mb)
                     zmn_sad(nvariables) = -zbs(nb,mb)
                  ELSE
                     n_sad(nvariables) = -nb
                     rmn_sad(nvariables) = rbc(nb,mb)
                     zmn_sad(nvariables) = zbs(nb,mb)
                  END IF
               END IF
            END DO
         END DO

         IF (nvariables .ne. numsurf_sad) THEN
            WRITE(6,*) "nvariables != numsurf_sad in load_saddle_surf"
            STOP
         END IF

         nvariables = irho_bdy

      END IF

      END SUBROUTINE load_saddle_wsurf
