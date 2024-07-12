      SUBROUTINE coil_group_report(cg,igroup,rR,zZ,nphi,phi_array)
!  Subroutine to print out information about a coil group
!  (Stored in a bsc_coilcoll)

!  JDH 2011-08-15 - First Version

      USE stel_kinds
      USE stel_constants
      
      USE bsc_T
       !  Access to derived types, magnetic field computation 
       
      USE sym_check, ONLY: check_symmetry

      IMPLICIT NONE
!-----------------------------------------------
!   A R G U M E N T   V a r i a b l e s
!-----------------------------------------------
      
      TYPE (bsc_coilcoll),  INTENT(inout) :: cg
      INTEGER, INTENT(in) :: igroup, nphi
      REAL(rprec), INTENT(in) :: rR, zZ
      REAL(rprec), DIMENSION(1:nphi), INTENT(in) :: phi_array

!  cg          A bsc_coilcoll (hodling a coil group) to report on
!  igroup      Coil group number (merely reported here)
!  rR          Cylindrical R at which to report B
!  zZ          Cylindrical Z at which to report B
!  nphi        Number of phi value at which to report B (length of phi_array)
!  phi_array   Array of cylindrical phi values at which to report B
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: big = 1.E40_rprec
      INTEGER, PARAMETER :: int_big = 987654321

      INTEGER :: n_coil_fl, n_coil_fc, n_coil_other, nfil_max_fl,              &
     &   nfil_min_fl, nfil_sum_fl, nfil_sumsq_fl

      INTEGER :: ic, ncoil, nfil, iphi, i

      REAL(rprec) :: cur_max_fl, cur_min_fl, cur_sum_fl,                       &
     &   cur_sumsq_fl, cur_max_fc, cur_min_fc, cur_sum_fc,                     &
     &   cur_sumsq_fc, ave_nfil_fl, sd_nfil_fl, cur_ave_fl,                    &
     &   cur_sd_fl, cur_ave_fc, cur_sd_fc
     
      TYPE(bsc_coil), POINTER :: coil
      
      REAL(rprec), DIMENSION(3) :: x, b
      
      REAL(rprec) :: c, s, br, bphi, bz, cur, phi

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------  
      
!  Find out how many coils
      ncoil = cg % ncoil
      IF (ncoil .le. 0) THEN
         WRITE(*,4000)  igroup, cg % s_name
4000  FORMAT(/80("*")/"No coils in coil group # ",i4,", group ID ",a30)
         RETURN
      ELSE
         WRITE(*,4100) igroup, cg % s_name
4100  FORMAT(/80("*")/"Coil group # ",i4,", with group ID ",a30)
      ENDIF

!  Initialize variables before loop over coils
      n_coil_fl = 0
      n_coil_fc = 0
      n_coil_other = 0
      nfil_max_fl = 0
      nfil_min_fl = int_big
      nfil_sum_fl = 0
      nfil_sumsq_fl = 0
      cur_max_fl = - big
      cur_min_fl = big
      cur_sum_fl = zero
      cur_sumsq_fl = zero
      cur_max_fc = - big
      cur_min_fc = big
      cur_sum_fc = zero
      cur_sumsq_fc = zero

!  Loop over coils
      DO ic = 1,ncoil
         coil => cg % coils(ic)
         cur = coil % current

!  Different coding, depending on c_type
         SELECT CASE (coil % c_type)
         CASE ('fil_loop','floop')
            n_coil_fl = n_coil_fl + 1
            
            nfil = SIZE(coil % xnod,2) - 1
            nfil_max_fl = MAX(nfil_max_fl,nfil)
            nfil_min_fl = MIN(nfil_min_fl,nfil)
            nfil_sum_fl = nfil_sum_fl + nfil
            nfil_sumsq_fl = nfil_sumsq_fl + nfil * nfil
            
            cur_max_fl = MAX(cur_max_fl,cur)
            cur_min_fl = MIN(cur_min_fl,cur)
            cur_sum_fl = cur_sum_fl + cur
            cur_sumsq_fl = cur_sumsq_fl + cur * cur

         CASE ('fil_circ','fcirc')
            n_coil_fc = n_coil_fc + 1
            cur_max_fc = MAX(cur_max_fc,cur)
            cur_min_fc = MIN(cur_min_fc,cur)
            cur_sum_fc = cur_sum_fc + cur
            cur_sumsq_fc = cur_sumsq_fc + cur * cur
         
         CASE DEFAULT
            n_coil_other = n_coil_other + 1
         
         END SELECT
      END DO
      
      ave_nfil_fl = nfil_sum_fl * one / MAX(n_coil_fl,1)
      sd_nfil_fl = SQRT(nfil_sumsq_fl * one / MAX(n_coil_fl,1) -               &
     &   ave_nfil_fl ** 2)
      cur_ave_fl = cur_sum_fl / MAX(n_coil_fl,1)
      cur_sd_fl =  SQRT(cur_sumsq_fl / MAX(n_coil_fl,1) -                      &
     &   cur_ave_fl ** 2)
      cur_ave_fc = cur_sum_fc / MAX(n_coil_fc,1)
      cur_sd_fc =  SQRT(cur_sumsq_fc / MAX(n_coil_fc,1) -                      &
     &   cur_ave_fc ** 2)

!  Write out results for all coils
      WRITE(*,5000)  
      WRITE(*,5100) n_coil_fl, n_coil_fc, n_coil_other, ncoil
      
      IF (n_coil_fl .gt. 0) THEN
         WRITE(*,5200) nfil_min_fl, nfil_max_fl, ave_nfil_fl,                  &
     &      sd_nfil_fl
         WRITE(*,5300) cur_min_fl, cur_max_fl, cur_ave_fl, cur_sd_fl
      ENDIF
      
      IF (n_coil_fc .gt. 0) THEN
         WRITE(*,5400) cur_min_fc, cur_max_fc, cur_ave_fc, cur_sd_fc
      ENDIF

5000  FORMAT(/"               Loops     Circles   Other   Total")
5100  FORMAT(" # coils ",4(3x,i7))
5200  FORMAT(/" Filamentary loops, number of filaments:"/                      &
     &     "      Min     Max      Average          SD "/                      &
     &   i8,3x,i8,2(2x,f11.1))
5300  FORMAT(/" Filamentary loops, current:"/                                  &
     &     "      Min          Max         Average          SD "/              &
     &   3(2x,es12.5),2x,es9.2)
5400  FORMAT(/" Filamentary circles, current:"/                                &
     &     "      Min          Max         Average          SD "/              &
     &   3(2x,es12.5),2x,es9.2)
     
!  Loop over phi values
      WRITE(*,6000) rR, zZ
      DO iphi = 1,nphi
         phi = phi_array(iphi)
         c = COS(phi)
         s = SIN(phi)
         x(1) = rR * c
         x(2) = rR * s
         x(3) = zZ
         CALL bsc_b(cg,x,b)
         br = b(1) * c + b(2) * s
         bphi = - b(1) * s + b(2) * c
         bz = b(3)
         WRITE(*,6100) phi, br, bphi, bz
      END DO
      
      CALL check_symmetry(igroup)

6000  FORMAT(/"Magnetic field with unit current multiplier, at R,Z=",          &
     &   2(2x,es12.5)/t7,"phi",t24,"B . rhat",t38,"B . phihat",                &
     &   t52,"B . zhat")
6100  FORMAT(4(3x,es12.5))
      
      END SUBROUTINE coil_group_report