!============================================================================
      MODULE v3read_wout
!--------------------------------------------------------------------------
! 
! FUNCTION: Module containing routines necessary to transfer data from a VMEC
!           output wout file into AJAX memory for flux coordinate conversions, etc
!
! Module created by J. Shields. Last modified 12/13/05
!
! Note:  The original code that these routines are based on can be found in 
!         files "setup_ajax.f90" & "track_ajax_dr.f90" that were written
!         by Wayne Houlberg and included in the 
!         downloadable .tar file containing the TRACK package.
!        (see e.g. ~shields/track/track/src/setup_ajax.f90
!                  ~shields/track/track/src/track_ajax_dr.f90)
!
!------------------------------------------------------------------------------

!      USE stel_kinds
!      USE stel_constants
!      USE track_mod
      USE AJAX_MOD
!      USE eq_T         ! module containing VMEC output info
!      USE math_utilities   ! module containing vector manipulation routines
!      USE intpol_cdf       ! module containing intpol netCDF i/o routines
!      USE ip_beamline  ! module containing int_pol TYPE definition
      IMPLICIT NONE


      CONTAINS

!=====================================================================
      SUBROUTINE READEQ_FROM_WOUT(in_woutfile)
!-----------------------------------------------------------------------------
!
! FUNCTION:  read in equilbrium info from an input VMEC wout file
!
! created 4/26/05 by J. Shields
! updated 1/6/06 (J. Shields) to take in_woutfile as a an input
!
! ***IMPORTANT CAVEAT:  note that VMEC and AJAX define arrays rmn, zmn, lmn, etc
!  as *TRANSPOSES* of one another.  For example VMEC defines rmn(mn_mode, radius),
!  while AJAX defines rmn  as  rmn(radius, mn_mode) !
!
!---------------------------------------------------------------------------- 
!      USE stel_kinds
      USE SPEC_KIND_MOD
      USE stel_constants      
      USE read_wout_mod, only : read_wout_file, rmnc, zmns, gmnc,              &
     &   currumnc, currvmnc, bsubumnc, bsubvmnc, lasym, rmns, zmnc,            &
     &   lmns, gmns, currumns, currvmns, bsubumns, bsubvmns, xm, xn,           &
     &   xm_nyq, xn_nyq, ns, mnmax, mnmax_nyq, iotaf,  phi
     
!      USE mgrid_mod, only : mgrid_mode
     
      IMPLICIT NONE

!  Declare Arguments 
      CHARACTER(len=*), INTENT(in) :: in_woutfile

!  Declare local variables

      LOGICAL ::                                                               &
     & l_axisymmetric, l_mfilter_ajax

      INTEGER:: ierr, iflag
      INTEGER, DIMENSION(3)   :: dims
      CHARACTER(len=*), PARAMETER ::  sub_name = 'READEQ_FROM_WOUT: '
      CHARACTER(len=60) ::     message
 
      INTEGER             ::                                                   &
     & n_eq,nrho_ajax,ntheta_ajax,nzeta_ajax


      INTEGER ::                                                               &
     & i,j,k, k_grid, k_pflx,                                                  &
     & nr_pflx, nr_rz, nk_rz, nk_lam, nr_lam, size_m, size_n, mn_size,         &
     & size_iotaf, nr_rzl, nk_rzl,                                             &
     & size1_rmnc, size2_rmnc, size1_zmns, size2_zmns,                         &
     & size1_lmns, size2_lmns, size1_gmnc, size2_gmnc, size_phi


  
      INTEGER, ALLOCATABLE ::                                                  &
     & my_xm(:), my_xn(:)

       REAL(KIND=rspec), ALLOCATABLE ::                                        &
     & rho_iotaf(:), rho_phi(:),                                               &
     & rho_rz(:),rho_lam(:),  my_iotaf(:), my_phi(:),                          &
     & my_rmnc(:,:),my_zmns(:,:), my_lmns(:,:), my_gmnc(:,:)

       REAL(KIND=rspec), ALLOCATABLE ::                                        &
     & rmnc_ajax(:,:),zmns_ajax(:,:), lmns_ajax(:,:)    ! AJAX compatible rmnc, etc.

       REAL(KIND=rspec) ::  phitot, phi_sum, my_temp, my_junk(23,56)
       REAL(KIND=rspec) ::  rho_rand(5), myphitot(5)
       REAL(KIND=rspec) ::  phi_ratio(5)
     
!  Start of executable code
!       write(*,*) 'Now entering', sub_name

!  Read the wout file
      ierr = 0
!      in_woutfile = 'wout_.nc'

      CALL read_wout_file(in_woutfile, ierr)
      if (ierr .ne. 0 ) write(*,*) sub_name // 'read_wout error'


!!  Non-array variables
!      this % mgrid_mode = mgrid_mode

!......apparently variable "radius" is not a PUBLIC variable.  JS 4/27/05      
!      my_radius = radius

!............array variables............................!
      dims(1:2) = shape(rmnc)
!      write(*,*) 'dims 1 and 2 for rmnc = ', dims(1), dims(2)
!      write(*,*) 'shape(rmnc) = ', shape(rmnc)
      ALLOCATE(my_rmnc(dims(1),dims(2)))
      my_rmnc = rmnc
      size1_rmnc = dims(1)
      size2_rmnc = dims(2)

      dims(1:2) = shape(zmns)
      ALLOCATE(my_zmns(dims(1),dims(2)))
      my_zmns = zmns
      size1_zmns = dims(1)
      size2_zmns = dims(2)

      dims(1:2) = shape(lmns)
      ALLOCATE(my_lmns(dims(1),dims(2)))
      my_lmns = lmns
      size1_lmns = dims(1)
      size2_lmns = dims(2)

      dims(1:2) = shape(gmnc)
      ALLOCATE(my_gmnc(dims(1),dims(2)))
      my_gmnc = gmnc
      size1_gmnc = dims(1)
      size2_gmnc = dims(2)

!      write(*,*) 'size of rmnc = ', size1_rmnc, size2_rmnc
!      write(*,*) 'size of zmns = ', size1_zmns, size2_zmns
!      write(*,*) 'size of lmnc = ', size1_lmns, size2_lmns
!      write(*,*) 'size of gmnc = ', size1_gmnc, size2_gmnc

!       do i = 1,12
!         do j = 1, 8
!           write(*,*) ' ORIGINAL rmnc(',i,j, ')',  rmnc(i,j)
!           write(*,*) ' COPIED  rmnc(',i,j, ')',  my_rmnc(i,j)
!         end do
!       end do

      dims(1:1) = shape(xm)
      ALLOCATE(my_xm(dims(1)))
      my_xm = INT(xm)
      size_m = dims(1)


      dims(1:1) = shape(xn)
      ALLOCATE(my_xn(dims(1)))
      my_xn = INT(xn)
      size_n = dims(1)

      mn_size = MAX(size_m, size_n)
      if ( size_m .ne. size_n) then
        write(*,*) '&&& WARNING: m and n array sizes are NOT EQUAL &&&'
      end if

!       write(*,*) ' ORIGINAL xm array = ', xm
!       write(*,*) ' COPIED xm array = ', my_xm

!...........now fill AJAX-compatible arrays as the *TRANSPOSE* of VMEC arrays....!
      
      ALLOCATE( rmnc_ajax(size2_rmnc,size1_rmnc) )
      ALLOCATE( zmns_ajax(size2_zmns,size1_zmns) )
      ALLOCATE( lmns_ajax(size2_lmns,size1_lmns) )

      rmnc_ajax = TRANSPOSE(my_rmnc)
      zmns_ajax = TRANSPOSE(my_zmns)
      lmns_ajax = TRANSPOSE(my_lmns)

!!......test to verify that transpose was done correctly.  JS 4/28/05....!
!      do i=1,size2_rmnc
!        do j=1,size1_rmnc
!          write(*,*) 'rmnc AJAX, reg', rmnc_ajax(i,j), my_rmnc(j,i)
!        end do
!      end do


!..........read in "iota on a full mesh" variable.............!
      dims(1:1) = shape(iotaf)
      ALLOCATE(my_iotaf(dims(1)))
      my_iotaf = iotaf
      size_iotaf = dims(1)

!       write(*,*) ' ORIGINAL iotaf array = ', iotaf
!       write(*,*) ' COPIED iotaf array = ', my_iotaf


!..........read in "Toroidal flux on full mesh" variable.............!
!...note that phi is defined(1:n) and my_phi is defined (O,(n-1) )
      dims(1:1) = shape(phi)
      ALLOCATE( my_phi(1:dims(1)) )
      my_phi = phi
      size_phi = dims(1)

!      write(*,*) 'size_phi = ', size_phi
!      write(*,*) 'phi array = ', phi
!      write(*,*) 'my_phi array = ', my_phi


!      IF (lasym) THEN
!         dims(1:2) = shape(rmns)
!         ALLOCATE(my_rmns(dims(1),dims(2)))
!         my_rmns = rmns
!         dims(1:2) = shape(zmnc)
!         ALLOCATE(my_zmnc(dims(1),dims(2)))
!         my_zmnc = zmnc
!         dims(1:2) = shape(gmns)
!         ALLOCATE(my_gmns(dims(1),dims(2)))
!         my_gmns = gmns
!      END IF

!..........at this point, all needed wout variables should be read in.......!





!!.......define nr_rzl, nk_rzl  based on size of rmnc,zmns, gmnc arrays..............!
      nr_rzl = MAX(size1_rmnc, size1_zmns, size1_gmnc) 
      nk_rzl = MAX(size2_rmnc, size2_zmns, size2_gmnc,size_m, size_n)

!......set nr_rz, nk_rz to rmnc size and nr_lam and nk_lam to lmns size.  JS 4/18/05
!.....recall that VMEC and AJAX use *transposed* array conventions! ............!
!.......( i.e VMEC uses rmnc(nk_rz,nr_rz), while AJAX uses rmn(nr_rz,nk_rz)......!
      nr_rz  = size2_rmnc
      nk_rz  = size1_rmnc
      nr_lam = size2_lmns
      nk_lam = size1_lmns

!......verify order of array filling.......................!

!      write(*,*) 'nr_rz SIZE = ', nr_rz
!      write(*,*) 'nk_rz SIZE = ', nk_rz
!      write(*,*) 'AJAX rmnc shape =  ', nr_rz, nk_rz, SHAPE(rmnc_ajax)
!      write(*,*) 'AJAX lmns shape =  ', nr_lam, nk_lam, SHAPE(lmns_ajax)

!.........allocate rho_rz, rho_lam arrays.......................................!
!      ALLOCATE( rho_rz(1:nr_rzl), rho_lam(1:nr_rzl) )
      ALLOCATE( rho_rz(1:nr_rz), rho_lam(1:nr_lam) )

!......fill rho_rz as an evenly spaced "full-grid", following Houlberg.........!
      do i = 1, nr_rz
        rho_rz(i) = SQRT( REAL(i-1) / REAL(nr_rz - 1) )
      end do

!.......fill rho_phi as an evenly spaced FULL grid as well..............!
      ALLOCATE( rho_phi(1:nr_rz)  )
      rho_phi = rho_rz

!......fill rho_lam as an evenly spaced "half-grid", following Houlberg......!
!........note that rho_lam(1) is automatically set to 0.0...............
      rho_lam(1) = 0.0
      do i = 2,nr_lam
        rho_lam(i) = SQRT( (REAL(i) - 1.5) / REAL(nr_lam - 1) )
      end do

!      do i = 1,nr_rz
!         my_temp = phi(i) / (rho_rz(i))**2.0
!         write(*,*) 'phi/rho^2 = ', my_temp 
!      end do

!.......allocate the array for radial polodial flux (associated w. iotaf)........!      
      ALLOCATE( rho_iotaf(size_iotaf) )

!......fill rho_iotaf as an evenly spaced "half-grid", following Houlberg. JS 4/5/05...!
!........note that rho_pflx(1) is automatically set to 0.0...............
      rho_iotaf(1) = 0.0
      do i = 2,size_iotaf
        rho_iotaf(i) = SQRT( (REAL(i) - 1.5) / REAL(size_iotaf - 1) )
      end do


!.........Verify that phi/(rho**2) is a constant equal to phitot.  Note that.....!
!...........rho for phi arrays must be a FULL grid to work correctly.............!
      do i = 2,4
        phi_ratio(i) = my_phi(i) / ( rho_phi(i)**2.0)
!        write(*,*) '(i, rho_phi, my_phi) = ', '(', i, rho_phi(i),              &
!     &               my_phi(i), ')'
!        write(*,*) 'phi ratio', i, ' = ', phi_ratio(i)
      end do

!.....temporarily "hardwire" phitot value to  0.507375002
!      phitot = 0.507375002
!.....instead, set phitot to the phi/rho^2 ratio................!
      phitot = phi_ratio(4)
 

!........set grid to be proportional to flux (*NOT* to SQRT(flux)!! ), (ie k_grid =1)
!.....note that Houlberg test data uses k_grid = 0 ! [ grid ~ SQRT(flux) ]..
      k_grid = 1

!......set AJAX poloidal flux representation = iotabar (i.e. "rotational transform")...!
      k_pflx = 1

!.......give l_mfilter_ajax a default value, even though we are not currently using
!.........this optional AJAX variable for the non-axisymmetric case.  (JS 3/23/05)
      l_mfilter_ajax = .TRUE.

!.......explicitly assign a value to nrho_ajax, etc.  
!.........for reference: AJAX default is (21,21,21) and WH uses:  (31,33,21)
!      nrho_ajax = 21
!      ntheta_ajax = 21
!      nzeta_ajax = 21
      nrho_ajax = nr_rz
      ntheta_ajax = 99
      nzeta_ajax = 21

      call AJAX_LOAD_INTERFACE(phitot, size_iotaf, rho_iotaf, my_iotaf,        &
     &           nr_rz, nk_rz, rho_rz, nr_lam, nk_lam, rho_lam,                &
     &           mn_size, my_xm, my_xn, rmnc_ajax, zmns_ajax, lmns_ajax,       &
     &           k_pflx, k_grid, l_mfilter_ajax, nrho_ajax, ntheta_ajax,       &
     &            nzeta_ajax,iflag, message)


      RETURN
      
      END SUBROUTINE READEQ_FROM_WOUT


!================================================================================
      SUBROUTINE JS_READEQ()
!-------------------------------------------------------------------------------
!
! FUNCTION: reads a simplified form of a 3-D VMEC equilbrium solution &
!           acts as an interface btwn the V3FIT "model" TYPE & AJAX ! (JS 3/17/05)
!
! **Based on the Wayne H. TRACK driver subroutine READEQ
!
!
!....modified to create rz & lambda arrays of UNIFORM size.     JS 4/6/05
!
!....modified to take lambda and iotaf from "model".  JS 7/28/05
!
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      USE track_mod
      USE stel_constants
      USE v3f_global   ! module containing VMEC output info (i.e. it SAVEs model_a)

      IMPLICIT NONE

!Declaration of input variables
      INTEGER             ::                                                   &
     & n_eq,nrho_ajax,ntheta_ajax,nzeta_ajax

!Declaration of output variables
      INTEGER              ::     iflag

      CHARACTER(len=60)     ::     message

!Declaration of local variables
      LOGICAL ::                                                               &
     & l_axisymmetric, l_mfilter_ajax

      INTEGER ::                                                               &
     & i,j,k,                                                                  &
     & k_grid,k_pflx,                                                          &
     & nr_pflx, nr_rz, nk_rz, nk_lam, nr_lam, size_m, size_n, mn_size,         &
     & nr_rz_temp, nk_rz_temp, nr_lam_temp, nk_lam_temp, ai_size,              &
     & nr_rzl, nk_rzl,                                                         &
     & size1_rmnc, size2_rmnc, size1_zmns, size2_zmns,                         &
     & size1_lmns, size2_lmns,                                                 &
     & nr_pflx_iterations, max_fill,                                           &
     & n_degree         ! degree of coeff() polynomial                         


  
      INTEGER, ALLOCATABLE ::                                                  &
     & m(:),n(:)

      REAL(KIND=rspec) ::                                                      &
     & phitot, stepsize, stepsum,                                              &
     & rho_pflx_stepsum, pflx_stepsum,                                         &
     & rho_stepsize, pflx_stepsize

      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rho_pflx(:),pflx(:),                                                    &
     & rho_rz(:),rho_lam(:),                                                   &
     & rmn(:,:),zmn(:,:),lmn(:,:), coeff(:)

      INTEGER :: ii, jj


!      write(*,*) 'now in JS_READEQ'

!..........verify that we are able to access the "model_a" data...........!
!      write(*,*) 'original model_a curtor = ',                                 &
!     &            model_a%eqstate%varp%curtor
!
!      write(*,*) 'MODEL_A aux1 rmnc(2,3) = ',                                  &
!     &            model_a%eqstate%aux1%rmnc(2,3)


 
!&&&
!........set grid to be proportional to flux (*NOT* to SQRT(flux)!! ), (ie k_grid =1)
!.....note that Houlberg test data uses k_grid = 0 ! [ grid ~ SQRT(flux) ]..
      k_grid = 1

      size_m = SIZE(model_a%eqstate%aux1%xm)
      size_n = SIZE(model_a%eqstate%aux1%xn)

      size1_rmnc = SIZE(model_a%eqstate%aux1%rmnc, 1)
      size2_rmnc = SIZE(model_a%eqstate%aux1%rmnc, 2)

      size1_zmns = SIZE(model_a%eqstate%aux1%zmns, 1)
      size2_zmns = SIZE(model_a%eqstate%aux1%zmns, 2)

      size1_lmns = SIZE(model_a%eqstate%aux1%lmns, 1)
      size2_lmns = SIZE(model_a%eqstate%aux1%lmns, 2)



!      write(*,*) 'rmnc dim 1 SIZE = ', size1_rmnc
!      write(*,*) 'rmnc dim 2 SIZE = ', size2_rmnc

!      write(*,*) 'zmns dim 1 SIZE = ', size1_zmns
!      write(*,*) 'zmns dim 2 SIZE = ', size2_zmns

!      write(*,*) 'lmns dim 1 SIZE = ', size1_lmns
!      write(*,*) 'lmns dim 2 SIZE = ', size2_lmns

      
!      write(*,*) ' m size = ', size_m
!      write(*,*) ' n size = ', size_n
!      write(*,*) ' m components = ', model_a%eqstate%aux1%xm


!!......define nr_rzl, nk_rzl  based on size of rmnc,zmns, gmnc arrays...........!
!......modified to be based on lmns instead of gmnc.  JS 7/15/05
      nr_rzl = MAX(size1_rmnc, size1_zmns, size1_lmns) 
      nk_rzl = MAX(size2_rmnc, size2_zmns, size2_lmns,size_m, size_n)


!......set nr_rz, nk_rz, nr_lam, nk_lam to rmnc size. Also note that the AJAX
!..... arrays will be the *TRANSPOSES* of the VMEC arrays!   JS 5/2/05
      nr_rz  = size2_rmnc
      nk_rz  = size1_rmnc
      nr_lam = size2_rmnc
      nk_lam = size1_rmnc

!      write(*,*) 'nr_rz SIZE = ', nr_rz
!      write(*,*) 'nk_rz SIZE = ', nk_rz




!......verify that nk_rz, nk_lam match with m,n array sizing...............!
      if( (nk_rz > size_m)  .OR. (nk_rz > size_n) .OR.                         &
     &    (nk_lam > size_m) .OR. (nk_lam > size_n)      ) then
        write(*,*) 'WARNING-- mismatch in m,n array sizing.'
      end if


!.......allocate and fill integer m,n arrays..................................!
      mn_size = MAX(nk_rz, nk_lam)
      ALLOCATE( m(1:mn_size), n(1:mn_size) )

      IF( mn_size <= size_m) then
        FORALL(i =1:mn_size ) m(i) = INT( model_a%eqstate%aux1%xm(i) )
      ELSE
        FORALL(i =1:size_m ) m(i) = INT( model_a%eqstate%aux1%xm(i) )
        FORALL(i =(size_m+1):mn_size ) m(i) = 0.0
      END IF

      IF( mn_size <= size_n) then
        FORALL(i =1:mn_size ) n(i) = INT( model_a%eqstate%aux1%xn(i) )
      ELSE
        FORALL(i =1:size_n ) n(i) = INT( model_a%eqstate%aux1%xn(i) )
        FORALL(i =(size_n+1):mn_size ) n(i) = 0.0
      END IF
      

!.......allocate integer rmn,zmn, lmn arrays..................................!
      ALLOCATE( rmn( 1:nr_rz, 1:nk_rz ) )
      ALLOCATE( zmn( 1:nr_rz, 1:nk_rz ) )
      ALLOCATE( lmn( 1:nr_lam, 1:nk_lam ) )

!.........fill rmn, zmn  (note: AJAX rmn, zmn arrays are TRANSPOSES of VMEC arrays!   
      rmn = TRANSPOSE(model_a%eqstate%aux1%rmnc)
      zmn = TRANSPOSE(model_a%eqstate%aux1%zmns)





!.........allocate (and initialize) the poloidal flux array......................!. 
!......(note that when k=1, pflx is the rotational transform (iotabar)..........!
!......also note: rho~sqrt(toroidal flux) & normalized to unity at boundary.....!
      nr_pflx = SIZE(model_a%eqstate%aux1%iotaf)
      ALLOCATE( pflx(1:nr_pflx) )
      pflx = model_a%eqstate%aux1%iotaf

!      write(*,*) 'pflx in JS_READEQ = ', pflx

!.........fill lambda array...........................................................!
      lmn = TRANSPOSE(model_a%eqstate%aux1%lmns)



!.........allocate rho_rz, rho_lam arrays.......................................!
      ALLOCATE( rho_rz(1:nr_rz), rho_lam(1:nr_lam) )


!&&&&.........ALL array allocations should now be done!..............&&&&&&........!


!......fill rho_rz as an evenly spaced "full-grid", following Houlberg. JS 4/5/05...!
      do i = 1, nr_rz
        rho_rz(i) = SQRT( REAL(i-1) / REAL(nr_rz - 1) )
      end do


!......fill rho_lam as an evenly spaced "half-grid", following Houlberg. JS 4/5/05...!
!........note that rho_lam(1) is automatically set to 0.0...............
      rho_lam(1) = 0.0
      do i = 2,nr_lam
        rho_lam(i) = SQRT( (REAL(i) - 1.5) / REAL(nr_lam - 1) )
      end do


!&&&
!.....now do the interface for  AJAX_LOAD_MAGFLUX..............!

!.........poloidal flux representation = iotabar  (i.e. "rotational transform")...!
      k_pflx = 1


!.............equate phitot with "phiedge" in subTYPE eq_param_var..............!
      phitot = model_a%eqstate%varp%phiedge
!      write(*,*) 'phitot in JS_READEQ = ', phitot


!.........allocate the array for the rho gridding of the poloidal flux.........!. 
      ALLOCATE( rho_pflx(1:nr_pflx) )
      rho_pflx = 0.0

!......fill rho_pflx as an evenly spaced "half-grid", following Houlberg. JS 4/5/05...!
!........note that rho_pflx(1) is automatically set to 0.0...............
      rho_pflx(1) = 0.0
      do i = 2,nr_pflx
        rho_pflx(i) = SQRT( (REAL(i) - 1.5) / REAL(nr_pflx - 1) )
      end do


!.......explicitly assign a value to nrho_ajax, etc.  
!.........for reference: AJAX default is (21,21,21) and WH uses:  (31,33,21)
       nrho_ajax = 21
       ntheta_ajax = 21
       nzeta_ajax = 21


!.......give l_mfilter_ajax a default value, even though we are not currently using
!.........this optional AJAX variable for the non-axisymmetric case.  (JS 3/23/05)
      l_mfilter_ajax = .TRUE.


!......new version 4/18/05 capable of calling nr_rz, nr_lam separately!
      call AJAX_LOAD_INTERFACE(phitot, nr_pflx, rho_pflx, pflx,                &
     &           nr_rz, nk_rz, rho_rz, nr_lam, nk_lam, rho_lam,                &
     &           mn_size, m, n, rmn, zmn, lmn,                                 &
     &           k_pflx, k_grid, l_mfilter_ajax,nrho_ajax, ntheta_ajax,        &
     &            nzeta_ajax, iflag, message)



 9999 CONTINUE


      DEALLOCATE(rho_pflx,pflx)
      DEALLOCATE(rho_rz,rho_lam,m,n,rmn,zmn,lmn)

!      write(*,*) 'now leaving JS_READEQ'

      RETURN
      END SUBROUTINE JS_READEQ





!===================================================================================
      SUBROUTINE TEST_AJAX_COORD_TRANS
!--------------------------------------------------------------------------------
!
! Function:  test pointers, etc in Fortran 95
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      USE track_mod
      USE eq_T         ! module containing VMEC output info
      IMPLICIT NONE

      INTEGER :: mysize = 99
      INTEGER :: counter = 0
      INTEGER :: myint, istat =1, i, iflag=0

      CHARACTER(len=80) :: mychar1
      CHARACTER(len=80) :: mychar2, message

      REAL(KIND=rspec), DIMENSION(3) ::                                        &
     &             r_car_test, r_cyl_test, r_flx_out,                          &
     &             r_flx_in, r_cyl_out, r_flx_diff,                            &
     &             r_cyl, r_cyl0, r_cyl1   ! cylindrical coordinate vectors



      REAL(KIND=rspec)             :: rtemp, x

    
      write(*,*) 'Hi.  Now in TEST_AJAX_COORD_TRANS'

!........assign DEFAULT values to arrays that will be output by AJAX subroutines....!
      r_cyl_out = (/ 0.700001, 0.7000001, 0.700001 /)
      r_flx_out  = (/ 0.080001, 0.750001, 0.750001 /)


!......assign cylinderical coord vectors  (rho, phi, z)   [note: phi in radians]...!
      r_cyl0 = (/ 1.0, 0.0, 2.0 /)
      r_cyl1 = (/ 0.0085052, 0.9000000 , 0.0189542 /)
!      write(*,*) 'r_cyl0 = ', r_cyl0
!      write(*,*) 'r_cyl1 = ', r_cyl1


!..........enter in default flux coord values.  We'll then do a transform and an
!...........inverse transform to see if these values are returned.  JS 3/16/05....!
!      r_flx_in  = (/ 0.12, 0.91, 0.91 /)
!      r_flx_in  = (/ 0.44, 0.56, 0.78 /)
!      r_flx_in  = (/ 0.69999, 0.56, 0.78 /)
      r_flx_in  = (/ 0.5555599, 3.002, 2.11 /)
!      r_flx_in  = (/ 0.6999, 1.41111, 0.95555 /)






!!.....now verify that rcar and rcyl arrays refer to same point..............!
!      call CYLIN_TO_CARTESIAN(r_cyl, r_car_test, istat)
!      write(*,*) 'output array r_car_test = ', r_car_test

      write(*,*) '--------------------------------------------------'

      write(*,*) 'DEFAULT r_cyl_out prior to AJAX_FLX2CYL = ', r_cyl_out
      write(*,*) 'r_flx_in into AJAX_FLX2CYL = ', r_flx_in


!.......debug flx to cyl trans FIRST before the reverse trans.......!
      call AJAX_FLX2CYL(r_flx_in, r_cyl_out, iflag, message)

      if ( iflag == 0) then
         write(*,*) 'output array r_cyl_out = ', r_cyl_out
      else if (iflag == -1 ) then
         write(*,*) 'AJAX_FLX2CYL warning.  Message is : ', message
      else if (iflag == 1 ) then
         write(*,*) 'AJAX_FLX2CYL error.  Message is : ', message
      else
         write(*,*) ' WARNING.  unknown AJAX_FLX2CYL iflag value.'
      end if


      write(*,*) '   '

!.........now try to convert r_cyl to flux
      write(*,*) '--------------------------------------------'
!      r_cyl_test = r_cyl1
      r_cyl_test = r_cyl_out
      write(*,*) 'DEFAULT r_flx_out prior to AJAX_CYL2FLX =', r_flx_out
      write(*,*) 'r_cyl_test into AJAX_CYL2FLX = ', r_cyl_test


      call AJAX_CYL2FLX(r_cyl_test, r_flx_out, iflag, message)


      if ( iflag == 0) then
         write(*,*) 'output array r_flx_out = ', r_flx_out
         r_flx_diff = r_flx_out - r_flx_in
         write(*,*) '------------------------------------------'
         write(*,*) 'NET RESULT : '
         write(*,*) 'ORIGINAL ( input) r_flx_in = ', r_flx_in
         write(*,*) 'FINAL (returned) r_flx_out = ', r_flx_out

      else if (iflag == -1 ) then
         write(*,*) 'AJAX_CYL2FLX warning.  Message is : ', message
      else if (iflag == 1 ) then
         write(*,*) 'AJAX_CYL2FLX error.  Message is : ', message
      else
         write(*,*) ' WARNING.  unknown AJAX_CYL2FLX iflag value.'
      end if


      write(*,*) ' '
      write(*,*) 'Hi.  Now leaving TEST_AJAX_COORD_TRANS'

      RETURN
      END SUBROUTINE TEST_AJAX_COORD_TRANS




!=============================================================================
      SUBROUTINE JS_TRACK_AJAX_DR
!-------------------------------------------------------------------------------
! FUNCTION: routine to test the module TRACK with AJAX
!   THIS version created 1/31/05 as a *simplified* version that assumes k_equil ==1
!   and thereby is able to combine the original TRACK_AJAX_DR and SETUP_AJAX
!   subroutines into a shorter and simpler program.
!
!References:
!  W.A.Houlberg, D.McCune 11/2001
!Comments:
!  Input options for the equilibrium geometry for the test cases include:
!    AJAX   - a simplified approximation to a tokamak (2D)
!           - a condensed set of VMEC output data (2D or 3D)
!
!  **This version is very similar to Houlberg's original TRACK_AJAX_DR driver subroutine,
!    except that the call to TRACK has been commented out and SETUP_AJAX has been
!  replaced by a direct call to WH_READEQ  JS  3/18/05
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      USE TRACK_MOD
      USE AJAX_MOD
      USE WRITE_MOD
      IMPLICIT NONE

!Declaration of namelist input variables
      CHARACTER(len=25) ::                                                     &
     & cn_eq

      INTEGER ::                                                               &
     & k_equil,k_seg,                                                          &
     & n_rho

      REAL(KIND=rspec) ::                                                      &
     & r0,a0,s0,e0,e1,d1,                                                      &
     & bt0,q0,q1,                                                              &
     & r_seg(3,2)

!Declaration of local variables
      CHARACTER(len=25) ::                                                     &
     & cn_msg,cn_sum,cn_tmp

      CHARACTER(len=120) ::                                                    &
     & label,message

      INTEGER ::                                                               &
     & n_msg,n_sum,n_tmp

      INTEGER ::                                                               &
     & i,                                                                      & 
     & iflag,                                                                  &
     & npro,n_int,nrho_ajax,ntheta_ajax,nzeta_ajax,                            &
     & n_seg

      INTEGER, ALLOCATABLE ::                                                  &
     & irho_int(:),izone_int(:)

      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rho(:),                                                                 &
     & s_int(:),rcar_int(:,:),rcyl_int(:,:),rflx_int(:,:),sdotb_int(:),        &
     & sdotphi_int(:)

!Output arrays
      INTEGER ::                                                               &
     & mx_npro

      PARAMETER(mx_npro=100)

      CHARACTER(len=15) ::                                                     &
     & namepro(mx_npro),unitpro(mx_npro)

      CHARACTER(len=79) ::                                                     &
     & descpro(mx_npro)

      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & valpro(:,:)

      NAMELIST/indata/k_equil,cn_eq,                                           &
     &                r0,a0,s0,e0,e1,d1,bt0,q0,q1,                             &
     &                k_seg,r_seg,                                             &
     &                n_rho

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Warning and error flag and message
      iflag=0
      message=''

!Points defining the test segment
      n_seg=2

!Points for internal AJAX grids
      nrho_ajax=31
      ntheta_ajax=33
      nzeta_ajax=21

!-------------------------------------------------------------------------------
!Description of namelist variables
!-------------------------------------------------------------------------------
!Units are enclosed in []

!Equilibrium specification:

           !Required:

!k_equil   -MHD equilibrium/geometry option [-]
!          =0 simple approximation to 2D plasma using AJAX
!            (specify r0,a0,s0,e0,e1,d1,bt0,q0,q1)
!          =1 read 2D or 3D inverse coordinate data file using AJAX
!            (specify file name in cn_eq)
!          =2 read 2D treq data file using XPLASMA
!            (specify file name in cn_eq)
!          =else not allowed

           !Depending on above choices:

!cn_eq     -name of input file for inverse coordinate data [character]
!r0        -major radius of geometric center [m]
!a0        -minor radius in midplane [m]
!s0        -axis shift normalized to a0 [-]
!e0        -axis elongation normalized to a0 [-]
!e1        -edge elongation normalized to a0 [-]
!d1        -edge triangularity normalized to a0 [-]
!bt0       -toroidal field at R0 [T]
!q0        -safety factor on axis [-]
!q1        -safety factor at edge [-]

!TRACK segment specification:

           !Required:

!k_seg     -option for coordinate system specifying the segments [-]
!          =1 flux coordinates (rho,theta,zeta)
!          =2 Cartesian coordinates (x,y,z)
!          =else default cylindrical coordinates (R,phi,Z)
!r_seg(3,2)-coordinates at beginning and end of test segment
!           (see k_seg for allowable choices)

!Plasma specification:

           !Required:

!n_rho     -no. radial plasma nodes [-]

!-------------------------------------------------------------------------------
!Read namelist file
!-------------------------------------------------------------------------------
      n_tmp=20
      cn_tmp='nml_track.dat'
      OPEN(UNIT=n_tmp,                                                         &
     &     FILE=cn_tmp,                                                        &
     &     STATUS='old',                                                       &
     &     ACCESS='sequential')

      READ(n_tmp,indata)

      CLOSE(UNIT=n_tmp)
!
!-------------------------------------------------------------------------------
!Open standard output files
!-------------------------------------------------------------------------------
!Summary
      n_sum=10
      cn_sum='sum_track.dat'
      OPEN(UNIT=n_sum,                                                         &
     &     FILE=cn_sum,                                                        &
     &     STATUS='unknown',                                                   &
     &     FORM='formatted')

!Messages
      n_msg=11
      cn_msg='msg_track.dat'
      OPEN(UNIT=n_msg,                                                         &
     &     FILE=cn_msg,                                                        &
     &     STATUS='unknown',                                                   &
     &     FORM='formatted')

!-------------------------------------------------------------------------------
!Set radial grid
!-------------------------------------------------------------------------------
      ALLOCATE(rho(1:n_rho))
      rho(1:n_rho)=REAL((/ (i-1,i=1,n_rho) /),rspec)/REAL(n_rho-1,rspec)

!-------------------------------------------------------------------------------
!Set up MHD equilibrium information
!-------------------------------------------------------------------------------
      IF(k_equil ==0 .OR. k_equil == 1) THEN

        !AJAX option for MHD equilibrium interface
        IF(k_equil == 1) THEN

          !Open data file for reading
          OPEN(UNIT=n_tmp,                                                     &
     &         FILE=cn_eq,                                                     &
     &         STATUS='old',                                                   &
     &         ACCESS='sequential')


        ENDIF


!.............determine values of info entering SETUP_AJAX.........!
!........(note that r0 thru q1 values are ALL 0.0.  n_rho, rho are NON-ZERO).

        write(*,*) 'r0 = ', r0, ' a0 = ', a0,  '  s0 = ', s0
        write(*,*) 'e0 = ', e0, ' e1 = ', e1, ' d1 = ', d1
        write(*,*) 'bt0 = ', bt0, ' q0 = ', q0, ' q1 = ', q1
        write(*,*) 'n_rho = ', n_rho
        write(*,*) 'nrho_ajax = ',nrho_ajax,'ntheta_ajax = ',ntheta_ajax 
!        write(*,*) 'Values in ARRAY rho = ', rho 

!        !Set up the internal data in AJAX
!        CALL SETUP_AJAX(k_equil,n_tmp,r0,a0,s0,e0,e1,d1,bt0,q0,q1,n_rho,       &
!     &                  rho,nrho_ajax,ntheta_ajax,nzeta_ajax,iflag,            &
!     &                  message)

!.......replace SETUP_AJAX call with a direct call to READEQ().  JS 1/31/05
        iflag=0
        message=''

!.........replaced READEQ call with WH_READEQ call (only difference btwn the two
!..........READEQ routines is that WH_READEQ bundles the AJAX calls into a....!
!..........separate subroutine AJAX_LOAD_INTERFACE (JS 3/22/05).......................!
!        CALL READEQ(n_tmp,nrho_ajax,ntheta_ajax,nzeta_ajax,                     &
!     &              iflag,message)

        CALL WH_READEQ(n_tmp,nrho_ajax,ntheta_ajax,nzeta_ajax,                  &
     &              iflag,message)


        !Check messages
        IF(iflag /= 0) THEN

          CALL WRITE_LINE(n_msg,message,1,1)
          IF(iflag > 0) GOTO 9999
          iflag=0
          message=''

        ENDIF

        IF(k_equil == 1) CLOSE(UNIT=n_tmp)


      ELSE

       !Illegal choice of k_equil
       iflag=1
        message='TRACK_DR/ERROR:illegal k_equil (use 0,1,2)'
        CALL WRITE_LINE(n_msg,message,1,1)
        GOTO 9999

      ENDIF

!!-------------------------------------------------------------------------------
!!Allocate arrays and call TRACK
!!-------------------------------------------------------------------------------
!!Dimensioned to allow up to 6 times as many intersections as surfaces
!      ALLOCATE(irho_int(1:6*n_rho),                                            &
!     &         izone_int(1:6*n_rho),                                           &
!     &         s_int(1:6*n_rho),                                               &
!     &         rcar_int(1:3,1:6*n_rho),                                        &
!     &         rcyl_int(1:3,1:6*n_rho),                                        &
!     &         rflx_int(1:3,1:6*n_rho),                                        &
!     &         sdotb_int(1:6*n_rho),                                           &
!     &         sdotphi_int(1:6*n_rho))
!
!      CALL TRACK(n_rho,rho,n_seg,r_seg,n_int,irho_int,s_int,iflag,             &
!     &           message,                                                      &
!     &           K_SEG=k_seg,                                                  &
!     &           IZONE_INT=izone_int,                                          &
!     &           RCAR_INT=rcar_int,                                            &
!     &           RCYL_INT=rcyl_int,                                            &
!     &           RFLX_INT=rflx_int,                                            &
!     &           SDOTB_INT=sdotb_int,                                          &
!     &           SDOTPHI_INT=sdotphi_int)
!
!      !Check messages
!      IF(iflag /= 0) THEN
!
!!...........replaced WRITE_LINE w. standard write(*,*) stmnt.  JS 1/28/05...!
!        CALL WRITE_LINE(n_msg,message,1,1)
!!        write(*,*) n_msg, message
! 
!        IF(iflag > 0) GOTO 9999
!        iflag=0
!        message=''
!
!      ENDIF
!
!!-------------------------------------------------------------------------------
!!Record information along test segment for output files
!!-------------------------------------------------------------------------------
!!Allocate and initialize output radial arrays
!      ALLOCATE(valpro(1:n_int,1:mx_npro))
!      namepro(:)=''
!      unitpro(:)=''
!      descpro(:)=''
!      valpro(:,:)=0.0_rspec
!      npro=0
!
!!Data
!      npro=npro+1
!      namepro(npro)='s'
!      unitpro(npro)='m'
!      descpro(npro)='Distance along path)'
!      valpro(1:n_int,npro)=s_int(1:n_int)
!
!      npro=npro+1
!      namepro(npro)='surf_no'
!      unitpro(npro)='-'
!      descpro(npro)='Surface intersected'
!      valpro(1:n_int,npro)=REAL(irho_int(1:n_int),rspec)
!
!      npro=npro+1
!      namepro(npro)='zone_no'
!      unitpro(npro)='-'
!      descpro(npro)='Zone entered'
!      valpro(1:n_int,npro)=REAL(izone_int(1:n_int),rspec)
!
!      npro=npro+1
!      namepro(npro)='x'
!      unitpro(npro)='m'
!      descpro(npro)='x of intersections'
!      valpro(1:n_int,npro)=rcar_int(1,1:n_int)
!
!      npro=npro+1
!      namepro(npro)='y'
!      unitpro(npro)='m'
!      descpro(npro)='y of intersections'
!      valpro(1:n_int,npro)=rcar_int(2,1:n_int)
!
!      npro=npro+1
!      namepro(npro)='z'
!      unitpro(npro)='m'
!      descpro(npro)='z of intersections'
!      valpro(1:n_int,npro)=rcar_int(3,1:n_int)
!
!      npro=npro+1
!      namepro(npro)='R'
!      unitpro(npro)='m'
!      descpro(npro)='R of intersections'
!      valpro(1:n_int,npro)=rcyl_int(1,1:n_int)
!
!      npro=npro+1
!      namepro(npro)='phi'
!      unitpro(npro)='-'
!      descpro(npro)='Toroidal angle of intersections'
!      valpro(1:n_int,npro)=rcyl_int(2,1:n_int)
!
!      npro=npro+1
!      namepro(npro)='Z'
!      unitpro(npro)='m'
!      descpro(npro)='Z of intersections'
!      valpro(1:n_int,npro)=rcyl_int(3,1:n_int)
!
!      npro=npro+1
!      namepro(npro)='rho'
!      unitpro(npro)='-'
!      descpro(npro)='Normalized minor radius of intersections'
!      valpro(1:n_int,npro)=rflx_int(1,1:n_int)
!
!      npro=npro+1
!      namepro(npro)='theta'
!      unitpro(npro)='-'
!      descpro(npro)='Poloidal angle of intersections'
!      valpro(1:n_int,npro)=rflx_int(2,1:n_int)
!
!      npro=npro+1
!      namepro(npro)='zeta'
!      unitpro(npro)='-'
!      descpro(npro)='Toroidal angle of intersections'
!      valpro(1:n_int,npro)=rflx_int(3,1:n_int)
!
!      npro=npro+1
!      namepro(npro)='sdotb'
!      unitpro(npro)='-'
!      descpro(npro)='Cos of angle between segment and B'
!      valpro(1:n_int,npro)=sdotb_int(1:n_int)
!
!      npro=npro+1
!      namepro(npro)='sdotphi'
!      unitpro(npro)='-'
!      descpro(npro)='Cos of angle with toroidal direction'
!      valpro(1:n_int,npro)=sdotphi_int(1:n_int)
!
!!Print path profiles to summary data file
!      label='*** Intersection Data ***'
!      CALL WRITE_OUT1(n_sum,'sum',n_int,npro,valpro,namepro,unitpro,1,         &
!     &                LABEL=label)
!
!!Print path profiles to 1D data file
!      n_tmp=20
!      cn_tmp='1d_s_track.dat'
!      OPEN(UNIT=n_tmp,                                                         &
!     &     FILE=cn_tmp,                                                        &
!     &     STATUS='unknown',                                                   &
!     &     FORM='formatted')
!
!      CALL WRITE_OUT1(n_tmp,'1d',n_int,npro,valpro,namepro,unitpro,-1)
!
!      CLOSE(UNIT=n_tmp)



 9999 CONTINUE

      return
      END SUBROUTINE JS_TRACK_AJAX_DR



!=============================================================================
      SUBROUTINE TRACK_AJAX_DR
!-------------------------------------------------------------------------------
!TRACK_AJAX_DR is a program to test the module TRACK with AJAX
!References:
!  W.A.Houlberg, D.McCune 11/2001
!Comments:
!  Input options for the equilibrium geometry for the test cases include:
!    AJAX   - a simplified approximation to a tokamak (2D)
!           - a condensed set of VMEC output data (2D or 3D)
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      USE TRACK_MOD
      USE AJAX_MOD
      USE WRITE_MOD
      IMPLICIT NONE

!Declaration of namelist input variables
      CHARACTER(len=25) ::                                                     &
     & cn_eq

      INTEGER ::                                                               &
     & k_equil,k_seg,                                                          &
     & n_rho

      REAL(KIND=rspec) ::                                                      &
     & r0,a0,s0,e0,e1,d1,                                                      &
     & bt0,q0,q1,                                                              &
     & r_seg(3,2)

!Declaration of local variables
      CHARACTER(len=25) ::                                                     &
     & cn_msg,cn_sum,cn_tmp

      CHARACTER(len=120) ::                                                    &
     & label,message

      INTEGER ::                                                               &
     & n_msg,n_sum,n_tmp

      INTEGER ::                                                               &
     & i,                                                                      & 
     & iflag,                                                                  &
     & npro,n_int,nrho_ajax,ntheta_ajax,nzeta_ajax,                            &
     & n_seg

      INTEGER, ALLOCATABLE ::                                                  &
     & irho_int(:),izone_int(:)

      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rho(:),                                                                 &
     & s_int(:),rcar_int(:,:),rcyl_int(:,:),rflx_int(:,:),sdotb_int(:),        &
     & sdotphi_int(:)

!Output arrays
      INTEGER ::                                                               &
     & mx_npro

      PARAMETER(mx_npro=100)

      CHARACTER(len=15) ::                                                     &
     & namepro(mx_npro),unitpro(mx_npro)

      CHARACTER(len=79) ::                                                     &
     & descpro(mx_npro)

      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & valpro(:,:)

      NAMELIST/indata/k_equil,cn_eq,                                           &
     &                r0,a0,s0,e0,e1,d1,bt0,q0,q1,                             &
     &                k_seg,r_seg,                                             &
     &                n_rho

       write(*,*) 'Now entering TRACK_AJAX_DR'

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Warning and error flag and message
      iflag=0
      message=''

!Points defining the test segment
      n_seg=2

!Points for internal AJAX grids
      nrho_ajax=31
      ntheta_ajax=33
      nzeta_ajax=21

!-------------------------------------------------------------------------------
!Description of namelist variables
!-------------------------------------------------------------------------------
!Units are enclosed in []

!Equilibrium specification:

           !Required:

!k_equil   -MHD equilibrium/geometry option [-]
!          =0 simple approximation to 2D plasma using AJAX
!            (specify r0,a0,s0,e0,e1,d1,bt0,q0,q1)
!          =1 read 2D or 3D inverse coordinate data file using AJAX
!            (specify file name in cn_eq)
!          =2 read 2D treq data file using XPLASMA
!            (specify file name in cn_eq)
!          =else not allowed

           !Depending on above choices:

!cn_eq     -name of input file for inverse coordinate data [character]
!r0        -major radius of geometric center [m]
!a0        -minor radius in midplane [m]
!s0        -axis shift normalized to a0 [-]
!e0        -axis elongation normalized to a0 [-]
!e1        -edge elongation normalized to a0 [-]
!d1        -edge triangularity normalized to a0 [-]
!bt0       -toroidal field at R0 [T]
!q0        -safety factor on axis [-]
!q1        -safety factor at edge [-]

!TRACK segment specification:

           !Required:

!k_seg     -option for coordinate system specifying the segments [-]
!          =1 flux coordinates (rho,theta,zeta)
!          =2 Cartesian coordinates (x,y,z)
!          =else default cylindrical coordinates (R,phi,Z)
!r_seg(3,2)-coordinates at beginning and end of test segment
!           (see k_seg for allowable choices)

!Plasma specification:

           !Required:

!n_rho     -no. radial plasma nodes [-]

!-------------------------------------------------------------------------------
!Read namelist file
!-------------------------------------------------------------------------------
      n_tmp=20
      cn_tmp='nml_track.dat'
      OPEN(UNIT=n_tmp,                                                         &
     &     FILE=cn_tmp,                                                        &
     &     STATUS='old',                                                       &
     &     ACCESS='sequential')

      READ(n_tmp,indata)

      CLOSE(UNIT=n_tmp)
!
!-------------------------------------------------------------------------------
!Open standard output files
!-------------------------------------------------------------------------------
!Summary
      n_sum=10
      cn_sum='sum_track.dat'
      OPEN(UNIT=n_sum,                                                         &
     &     FILE=cn_sum,                                                        &
     &     STATUS='unknown',                                                   &
     &     FORM='formatted')

!Messages
      n_msg=11
      cn_msg='msg_track.dat'
      OPEN(UNIT=n_msg,                                                         &
     &     FILE=cn_msg,                                                        &
     &     STATUS='unknown',                                                   &
     &     FORM='formatted')

!-------------------------------------------------------------------------------
!Set radial grid
!-------------------------------------------------------------------------------
      ALLOCATE(rho(1:n_rho))
      rho(1:n_rho)=REAL((/ (i-1,i=1,n_rho) /),rspec)/REAL(n_rho-1,rspec)

!-------------------------------------------------------------------------------
!Set up MHD equilibrium information
!-------------------------------------------------------------------------------
      IF(k_equil ==0 .OR. k_equil == 1) THEN

        !AJAX option for MHD equilibrium interface
        IF(k_equil == 1) THEN

          !Open data file for reading
          OPEN(UNIT=n_tmp,                                                     &
     &         FILE=cn_eq,                                                     &
     &         STATUS='old',                                                   &
     &         ACCESS='sequential')


        ENDIF


!.............determine values of info entering SETUP_AJAX.........!
!........(note that r0 thru q1 values are ALL 0.0.  n_rho, rho are NON-ZERO).
!        write(*,*) 'r0 = ', r0, ' a0 = ', a0,  '  s0 = ', s0
!        write(*,*) 'e0 = ', e0, ' e1 = ', e1, ' d1 = ', d1
!        write(*,*) 'bt0 = ', bt0, ' q0 = ', q0, ' q1 = ', q1
!        write(*,*) 'n_rho = ', n_rho
!        write(*,*) 'nrho_ajax = ',nrho_ajax,'ntheta_ajax = ',ntheta_ajax 
!        write(*,*) 'Values in ARRAY rho = ', rho 

        !Set up the internal data in AJAX
        CALL SETUP_AJAX(k_equil,n_tmp,r0,a0,s0,e0,e1,d1,bt0,q0,q1,n_rho,       &
     &                  rho,nrho_ajax,ntheta_ajax,nzeta_ajax,iflag,            &
     &                  message)

        !Check messages
        IF(iflag /= 0) THEN

!...........replaced WRITE_LINE w. standard write(*,*) stmnt.  JS 1/28/05...!
          CALL WRITE_LINE(n_msg,message,1,1)
!          write(*,*) n_msg, message

          IF(iflag > 0) GOTO 9999
          iflag=0
          message=''

        ENDIF

        IF(k_equil == 1) CLOSE(UNIT=n_tmp)

!      ELSEIF(k_equil == 2) THEN
!
!        !XPLASMA option for MHD equilibrium interface
!
!        !Open data file for reading
!        OPEN(UNIT=n_tmp,                                                       &
!     &       FILE=cn_eq,                                                       &
!     &       STATUS='old',                                                     &
!     &       ACCESS='sequential')
!
!        !Set up the internal data in XPLASMA
!         CALL SETUP_XPLASMA(n_tmp,cn_eq,iflag,message)
!
!        !Check messages
!        IF(iflag /= 0) THEN
!
!          CALL WRITE_LINE(n_msg,message,1,1)
!          IF(iflag > 0) GOTO 9999
!          iflag=0
!          message=''
!
!        ENDIF
!
!        CLOSE(UNIT=n_tmp)

      ELSE

       !Illegal choice of k_equil
       iflag=1
        message='TRACK_DR/ERROR:illegal k_equil (use 0,1,2)'
        CALL WRITE_LINE(n_msg,message,1,1)
        GOTO 9999

      ENDIF

       write(*,*) 'Now calling TRACK subroutine from TRACK_AJAX_DR'

!-------------------------------------------------------------------------------
!Allocate arrays and call TRACK
!-------------------------------------------------------------------------------
!Dimensioned to allow up to 6 times as many intersections as surfaces
      ALLOCATE(irho_int(1:6*n_rho),                                            &
     &         izone_int(1:6*n_rho),                                           &
     &         s_int(1:6*n_rho),                                               &
     &         rcar_int(1:3,1:6*n_rho),                                        &
     &         rcyl_int(1:3,1:6*n_rho),                                        &
     &         rflx_int(1:3,1:6*n_rho),                                        &
     &         sdotb_int(1:6*n_rho),                                           &
     &         sdotphi_int(1:6*n_rho))

      CALL TRACK(n_rho,rho,n_seg,r_seg,n_int,irho_int,s_int,iflag,             &
     &           message,                                                      &
     &           K_SEG=k_seg,                                                  &
     &           IZONE_INT=izone_int,                                          &
     &           RCAR_INT=rcar_int,                                            &
     &           RCYL_INT=rcyl_int,                                            &
     &           RFLX_INT=rflx_int,                                            &
     &           SDOTB_INT=sdotb_int,                                          &
     &           SDOTPHI_INT=sdotphi_int)

      !Check messages

      write(*,*) 'iflag value after TRACK call = ', iflag
      IF(iflag /= 0) THEN

!...........replaced WRITE_LINE w. standard write(*,*) stmnt.  JS 1/28/05...!
        CALL WRITE_LINE(n_msg,message,1,1)
!        write(*,*) n_msg, message
 
        IF(iflag > 0) GOTO 9999
        iflag=0
        message=''

      ENDIF

!-------------------------------------------------------------------------------
!Record information along test segment for output files
!-------------------------------------------------------------------------------
!Allocate and initialize output radial arrays
      ALLOCATE(valpro(1:n_int,1:mx_npro))
      namepro(:)=''
      unitpro(:)=''
      descpro(:)=''
      valpro(:,:)=0.0_rspec
      npro=0

!Data
      npro=npro+1
      namepro(npro)='s'
      unitpro(npro)='m'
      descpro(npro)='Distance along path)'
      valpro(1:n_int,npro)=s_int(1:n_int)

      npro=npro+1
      namepro(npro)='surf_no'
      unitpro(npro)='-'
      descpro(npro)='Surface intersected'
      valpro(1:n_int,npro)=REAL(irho_int(1:n_int),rspec)

      npro=npro+1
      namepro(npro)='zone_no'
      unitpro(npro)='-'
      descpro(npro)='Zone entered'
      valpro(1:n_int,npro)=REAL(izone_int(1:n_int),rspec)

      npro=npro+1
      namepro(npro)='x'
      unitpro(npro)='m'
      descpro(npro)='x of intersections'
      valpro(1:n_int,npro)=rcar_int(1,1:n_int)

      npro=npro+1
      namepro(npro)='y'
      unitpro(npro)='m'
      descpro(npro)='y of intersections'
      valpro(1:n_int,npro)=rcar_int(2,1:n_int)

      npro=npro+1
      namepro(npro)='z'
      unitpro(npro)='m'
      descpro(npro)='z of intersections'
      valpro(1:n_int,npro)=rcar_int(3,1:n_int)

      npro=npro+1
      namepro(npro)='R'
      unitpro(npro)='m'
      descpro(npro)='R of intersections'
      valpro(1:n_int,npro)=rcyl_int(1,1:n_int)

      npro=npro+1
      namepro(npro)='phi'
      unitpro(npro)='-'
      descpro(npro)='Toroidal angle of intersections'
      valpro(1:n_int,npro)=rcyl_int(2,1:n_int)

      npro=npro+1
      namepro(npro)='Z'
      unitpro(npro)='m'
      descpro(npro)='Z of intersections'
      valpro(1:n_int,npro)=rcyl_int(3,1:n_int)

      npro=npro+1
      namepro(npro)='rho'
      unitpro(npro)='-'
      descpro(npro)='Normalized minor radius of intersections'
      valpro(1:n_int,npro)=rflx_int(1,1:n_int)

      npro=npro+1
      namepro(npro)='theta'
      unitpro(npro)='-'
      descpro(npro)='Poloidal angle of intersections'
      valpro(1:n_int,npro)=rflx_int(2,1:n_int)

      npro=npro+1
      namepro(npro)='zeta'
      unitpro(npro)='-'
      descpro(npro)='Toroidal angle of intersections'
      valpro(1:n_int,npro)=rflx_int(3,1:n_int)

      npro=npro+1
      namepro(npro)='sdotb'
      unitpro(npro)='-'
      descpro(npro)='Cos of angle between segment and B'
      valpro(1:n_int,npro)=sdotb_int(1:n_int)

      npro=npro+1
      namepro(npro)='sdotphi'
      unitpro(npro)='-'
      descpro(npro)='Cos of angle with toroidal direction'
      valpro(1:n_int,npro)=sdotphi_int(1:n_int)

!Print path profiles to summary data file
      label='*** Intersection Data ***'
      CALL WRITE_OUT1(n_sum,'sum',n_int,npro,valpro,namepro,unitpro,1,         &
     &                LABEL=label)

!Print path profiles to 1D data file
      n_tmp=20
      cn_tmp='1d_s_track.dat'
      OPEN(UNIT=n_tmp,                                                         &
     &     FILE=cn_tmp,                                                        &
     &     STATUS='unknown',                                                   &
     &     FORM='formatted')

      CALL WRITE_OUT1(n_tmp,'1d',n_int,npro,valpro,namepro,unitpro,-1)

      CLOSE(UNIT=n_tmp)

 9999 CONTINUE

      return
      END SUBROUTINE TRACK_AJAX_DR



!================================================================================
      SUBROUTINE SETUP_AJAX(k_equil,n_eq,r0,a0,s0,e0,e1,d1,bt0,q0,q1,          &
     &                      n_rho,rho,nrho_ajax,ntheta_ajax,nzeta_ajax,        &
     &                      iflag,message)
!-------------------------------------------------------------------------------
!SETUP_AJAX sets up the equilibrium information for AJAX
!References:
!  W.A.Houlberg 11/2001
!Input:
!  k_equil             -MHD equilibrium/geometry option [-]
!                      =1 read inverse coordinate data file
!                      =else use simple approximation to 2-D plasma
!  n_eq                -input unit number [-]  (ie input FILE unit number) 
!  r0                  -major radius of geometric center [m]
!  a0                  -minor radius in midplane [m]
!  s0                  -axis shift normalized to a0 [-]
!  e0                  -axis elongation normalized to a0 [-]
!  e1                  -edge elongation normalized to a0 [-]
!  d1                  -edge triangularity normalized to a0 [-]
!  bt0                 -toroidal field at r0 [T]
!  q0                  -safety factor on axis [-]
!  q1                  -safety factor at edge [-]
!  n_rho               -no. radial nodes for output [-]
!  rho(n_rho)          -radial grid [rho]
!  nrho_ajax           -radial nodes to use in internal AJAX data [-]
!  ntheta_ajax         -poloidal nodes to use in internal AJAX data [-]
!  nzeta_ajax          -toroidal nodes to use in 3-D internal AJAX data [-]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      USE AJAX_MOD
      IMPLICIT NONE

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & k_equil,n_eq,n_rho,nrho_ajax,ntheta_ajax,nzeta_ajax

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & r0,a0,s0,e0,e1,d1,bt0,q0,q1,                                            &
     & rho(1:n_rho)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of local variables
      INTEGER ::                                                               &
     & k_pflx

      REAL(kind=rspec) ::                                                      &
     & phitot,z_pi

      REAL(KIND=rspec) ::                                                      &
     & q(1:n_rho),v1(1:n_rho),v2(1:n_rho)

!Machine and physical constants
      z_pi=ACOS(-1.0)

      IF(k_equil == 1) THEN
!-------------------------------------------------------------------------------
!Read inverse coordinate expansion
!-------------------------------------------------------------------------------
        iflag=0
        message=''
        CALL READEQ(n_eq,nrho_ajax,ntheta_ajax,nzeta_ajax,iflag,message)

        !Check messages
        IF(iflag > 0) THEN

          message='SETUP_AJAX(1)/'//message
          GOTO 9999

        ENDIF

      ELSE
!-------------------------------------------------------------------------------
!Load geometry using approximations to axisymmetric plasma
!-------------------------------------------------------------------------------
!Load boundary values of R,Z and fill in with approximations
        iflag=0
        message=''
        CALL AJAX_LOAD_RZBDY(r0,a0,s0,e0,e1,d1,iflag,message,                  &
     &                       NRHO_AJAX=nrho_ajax,                              &
     &                       NTHETA_AJAX=ntheta_ajax)

        !Check messages
        IF(iflag > 0) THEN

          message='SETUP_AJAX(2)/'//message
          GOTO 9999

        ENDIF

!Get metrics to evaluate total toroidal magnetic flux
        iflag=0
        message=''
        CALL AJAX_FLUXAV_G(n_rho,rho,iflag,message,RM2_R=v1,VP_R=v2)
        phitot=r0*bt0/4.0_rspec/z_pi*rho(n_rho)*v1(n_rho)*v2(n_rho)

        !Check messages
        IF(iflag > 0) THEN

          message='SETUP_AJAX(3)/'//message
          GOTO 9999

        ENDIF

        k_pflx=0
        q(1:n_rho)=q0+(q1-q0)*rho(1:n_rho)**2
        iflag=0
        message=''
        CALL AJAX_LOAD_MAGFLUX(phitot,k_pflx,n_rho,rho,q,iflag,message)

        !Check messages
        IF(iflag > 0) THEN

          message='SETUP_AJAX(4)/'//message
          GOTO 9999

        ENDIF

      ENDIF

 9999 CONTINUE

      END SUBROUTINE SETUP_AJAX

!================================================================================
      SUBROUTINE READEQ(n_eq,nrho_ajax,ntheta_ajax,nzeta_ajax,iflag,           &
     &                  message)
!-------------------------------------------------------------------------------
!READEQ reads a simplified form of a 3-D VMEC equilbrium solution
!
!     ***ie it is an interface btwn VMEC and TRACK!  (JS 1/31/05)
!
!References:
!  W.A.Houlberg 11/2001
!Input:
!  n_eq                -input unit number [-]    (i.e input FILE unit number!!)
!  nrho_ajax           -radial nodes to use in internal AJAX data [-]
!  ntheta_ajax         -poloidal nodes to use in internal AJAX data [-]
!  nzeta_ajax          -toroidal nodes to use in 3-D internal AJAX data [-]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      USE AJAX_MOD
      IMPLICIT NONE

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & n_eq,nrho_ajax,ntheta_ajax,nzeta_ajax

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of local variables
      LOGICAL ::                                                               &
     & l_axisymmetric

      INTEGER ::                                                               &
     & i,k,                                                                    &
     & k_grid,k_pflx,                                                          &
     & nr_iota,nr_rzl,nk_rzl

      INTEGER, ALLOCATABLE ::                                                  &
     & m(:),n(:)

      REAL(KIND=rspec) ::                                                      &
     & phitot

      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rho_iota(:),iotabar(:),                                                 &
     & rho_rz(:),rho_lam(:),                                                   &
     & rmn(:,:),zmn(:,:),lmn(:,:)

      INTEGER :: ii, jj

!phitot [Wb]
      READ(n_eq,*) phitot

!iotabar [-], rho~sqrt(toroidal flux) and normalized to unity at boundary
      READ(n_eq,*) nr_iota
      ALLOCATE(rho_iota(1:nr_iota),iotabar(1:nr_iota))
      READ(n_eq,*) (rho_iota(i),i=1,nr_iota)
      READ(n_eq,*) (iotabar(i),i=1,nr_iota)


      write(*,*) 'phitot in READEQ = ', phitot
      write(*,*) 'nr_iota in READEQ = ', nr_iota
!      write(*,*) 'rho_iota in READEQ = ', rho_iota

!R,Z and lambda grids, rho~sqrt(toroidal flux) and arbitrary normalization
      k_grid=0
      READ(n_eq,*) nr_rzl
      ALLOCATE(rho_rz(1:nr_rzl),rho_lam(1:nr_rzl))
      READ(n_eq,*) (rho_rz(i),i=1,nr_rzl)
      READ(n_eq,*) (rho_lam(i),i=1,nr_rzl)

      write(*,*)  'nr_rzl in READEQ = ', nr_rzl
!      write(*,*)  'rho_rz in READEQ = ', rho_rz

!R,Z and lambda modes
      READ(n_eq,*) nk_rzl

!Modes and amplitudes
      ALLOCATE(m(nk_rzl),n(nk_rzl),rmn(1:nr_rzl,1:nk_rzl),                     &
     &           zmn(1:nr_rzl,1:nk_rzl),lmn(1:nr_rzl,1:nk_rzl))

      !Check to see whether this is a tokamak or stellarator
      l_axisymmetric=.TRUE.

      DO k=1,nk_rzl

        READ(n_eq,*) m(k),n(k)
!        write(*,*) 'k, m(k), n(k) in READEQ = ', k, m(k), n(k)

        IF(n(k) /= 0 .AND. l_axisymmetric) l_axisymmetric=.FALSE.
        READ(n_eq,*) (rmn(i,k),i=1,nr_rzl)
        READ(n_eq,*) (zmn(i,k),i=1,nr_rzl)
        READ(n_eq,*) (lmn(i,k),i=1,nr_rzl)

      ENDDO

!      write(*,*) 'm array in READEQ = ', m
!      write(*,*) 'n array in READEQ = ', n

!!..........see which r_ij  and z_ij coeffs are zero.  JS 2/3/05
!      do ii = 1, nr_rzl
!        do jj = 1, nk_rzl
! 
!           if( rmn(ii,jj) == 0.0 ) then
!              write(*,*) 'zero Rij pair in READEQ = ', ii, jj
!           else if( zmn(ii,jj) == 0.0 ) then
!              write(*,*) 'zero Zij pair in READEQ = ', ii, jj
!           end if
!        end do
!      end do
!
!      do ii = 1, nr_rzl
!        do jj = 1, nk_rzl
! 
!           if( rmn(ii,jj) == 0.0  .AND. zmn(ii,jj) == 0.0 ) then
!              write(*,*) 'zero Rij AND Zij pair in READEQ = ', ii, jj
!           end if
!        end do
!      end do

      IF(l_axisymmetric) THEN

        !Axisymmetric plasma, calculate lambdas in AJAX on the internal grid
        iflag=0
        message=''
        CALL AJAX_LOAD_RZLAM(nr_rzl,nk_rzl,rho_rz,m,n,rmn,zmn,iflag,           &
     &                       message,                                          &
     &                       K_GRID=k_grid,                                    &
     &                       L_MFILTER_AJAX=.TRUE.,                            &
     &                       NRHO_AJAX=nrho_ajax,                              &
     &                       NTHETA_AJAX=ntheta_ajax)

        !Check messages
        IF(iflag > 0) THEN

          message='READEQ(1)/'//message
          GOTO 9999

        ENDIF

      ELSE

        !Non-axisymmetric plasma, use the lambdas from the data file
        iflag=0
        message=''
        CALL AJAX_LOAD_RZLAM(nr_rzl,nk_rzl,rho_rz,m,n,rmn,zmn,iflag,           &
     &                       message,                                          &
     &                       K_GRID=k_grid,                                    &
     &                       NRHO_AJAX=nrho_ajax,                              &
     &                       NTHETA_AJAX=ntheta_ajax,                          &
     &                       NZETA_AJAX=nzeta_ajax,                            &
     &                       NR_LAM=nr_rzl,                                    &
     &                       NK_LAM=nk_rzl,                                    &
     &                       RHO_LAM=rho_lam,                                  &
     &                       LAM=lmn)

        !Check messages
        IF(iflag > 0) THEN

          message='READEQ(2)/'//message
          GOTO 9999

        ENDIF

      ENDIF


!      write(*,*) 'iotabar in READEQ = ', iotabar
!      write(*,*) 'rho_iota in READEQ = ', rho_iota

      iflag=0
      message=''
      k_pflx=1
      CALL AJAX_LOAD_MAGFLUX(phitot,k_pflx,nr_iota,rho_iota,iotabar,           &
     &                       iflag,message)

      !Check messages
      IF(iflag > 0) message='READEQ(3)/'//message

 9999 CONTINUE

      DEALLOCATE(rho_iota,iotabar)

      DEALLOCATE(rho_rz,rho_lam,m,n,rmn,zmn,lmn)

      END SUBROUTINE READEQ

!================================================================================
      SUBROUTINE WH_READEQ(n_eq,nrho_ajax,ntheta_ajax,nzeta_ajax,iflag,           &
     &                  message)
!-------------------------------------------------------------------------------
!READEQ reads a simplified form of a 3-D VMEC equilbrium solution
!
!     ***ie it is an interface btwn VMEC and TRACK!  (JS 1/31/05)
!        
!   *this subroutine is essentially the same as Wayne Houlberg's original READEQ()
!    subroutine, except that the calls to AJAX_LOAD_RZLAM & AJAX_LOAD_MAGFLUX have
!    been moved to a separate subroutine called AJAX_LOAD_INTERFACE
!
!References:
!  W.A.Houlberg 11/2001
!Input:
!  n_eq                -input unit number [-]    (i.e input FILE unit number!!)
!  nrho_ajax           -radial nodes to use in internal AJAX data [-]
!  ntheta_ajax         -poloidal nodes to use in internal AJAX data [-]
!  nzeta_ajax          -toroidal nodes to use in 3-D internal AJAX data [-]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      USE AJAX_MOD
      IMPLICIT NONE

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & n_eq,nrho_ajax,ntheta_ajax,nzeta_ajax

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of local variables
      LOGICAL ::                                                               &
     & l_axisymmetric, l_mfilter_ajax

      INTEGER ::                                                               &
     & i,k,                                                                    &
     & k_grid,k_pflx,                                                          &
     & nr_iota,nr_rzl,nk_rzl

      INTEGER, ALLOCATABLE ::                                                  &
     & m(:),n(:)

      REAL(KIND=rspec) ::                                                      &
     & phitot

      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rho_iota(:),iotabar(:),                                                 &
     & rho_rz(:),rho_lam(:),                                                   &
     & rmn(:,:),zmn(:,:),lmn(:,:)

      INTEGER :: ii, jj

!phitot [Wb]
      READ(n_eq,*) phitot

!iotabar [-], rho~sqrt(toroidal flux) and normalized to unity at boundary
      READ(n_eq,*) nr_iota
      ALLOCATE(rho_iota(1:nr_iota),iotabar(1:nr_iota))
      READ(n_eq,*) (rho_iota(i),i=1,nr_iota)
      READ(n_eq,*) (iotabar(i),i=1,nr_iota)


      write(*,*) 'phitot in WH_READEQ = ', phitot
      write(*,*) 'nr_iota in WH_READEQ = ', nr_iota
!      write(*,*) 'rho_iota in WH_READEQ = ', rho_iota

!R,Z and lambda grids, rho~sqrt(toroidal flux) and arbitrary normalization
      k_grid=0
      READ(n_eq,*) nr_rzl
      ALLOCATE(rho_rz(1:nr_rzl),rho_lam(1:nr_rzl))
      READ(n_eq,*) (rho_rz(i),i=1,nr_rzl)
      READ(n_eq,*) (rho_lam(i),i=1,nr_rzl)

      write(*,*)  'nr_rzl in WH_READEQ = ', nr_rzl
!      write(*,*)  'rho_rz in WH_READEQ = ', rho_rz

!R,Z and lambda modes
      READ(n_eq,*) nk_rzl

!Modes and amplitudes
      ALLOCATE(m(nk_rzl),n(nk_rzl),rmn(1:nr_rzl,1:nk_rzl),                     &
     &           zmn(1:nr_rzl,1:nk_rzl),lmn(1:nr_rzl,1:nk_rzl))

      !Check to see whether this is a tokamak or stellarator
      l_axisymmetric=.TRUE.

      DO k=1,nk_rzl

        READ(n_eq,*) m(k),n(k)
!        write(*,*) 'k, m(k), n(k) in WH_READEQ = ', k, m(k), n(k)

        IF(n(k) /= 0 .AND. l_axisymmetric) l_axisymmetric=.FALSE.
        READ(n_eq,*) (rmn(i,k),i=1,nr_rzl)
        READ(n_eq,*) (zmn(i,k),i=1,nr_rzl)
        READ(n_eq,*) (lmn(i,k),i=1,nr_rzl)

      ENDDO

!      write(*,*) 'm array in WH_READEQ = ', m
!      write(*,*) 'n array in WH_READEQ = ', n

!!..........see which r_ij  and z_ij coeffs are zero.  JS 2/3/05
!      do ii = 1, nr_rzl
!        do jj = 1, nk_rzl
! 
!           if( rmn(ii,jj) == 0.0 ) then
!              write(*,*) 'zero Rij pair in WH_READEQ = ', ii, jj
!           else if( zmn(ii,jj) == 0.0 ) then
!              write(*,*) 'zero Zij pair in WH_READEQ = ', ii, jj
!           end if
!        end do
!      end do
!
!      do ii = 1, nr_rzl
!        do jj = 1, nk_rzl
! 
!           if( rmn(ii,jj) == 0.0  .AND. zmn(ii,jj) == 0.0 ) then
!              write(*,*) 'zero Rij AND Zij pair in READEQ = ', ii, jj
!           end if
!        end do
!      end do

      k_pflx = 1
      l_mfilter_ajax = .TRUE.

      call AJAX_LOAD_INTERFACE(phitot, nr_iota, rho_iota, iotabar,             &
     &           nr_rzl, nk_rzl, rho_rz, nr_rzl, nk_rzl, rho_lam,              &
     &           nk_rzl, m, n, rmn, zmn, lmn,                                  &
     &           k_pflx, k_grid, l_mfilter_ajax,nrho_ajax, ntheta_ajax,        &
     &            nzeta_ajax, iflag, message)

 9999 CONTINUE

      DEALLOCATE(rho_iota,iotabar)

      DEALLOCATE(rho_rz,rho_lam,m,n,rmn,zmn,lmn)

      END SUBROUTINE WH_READEQ


!================================================================================
      SUBROUTINE AJAX_LOAD_INTERFACE(phitot, nr_iota, rho_iota, iotabar,       &
     &           nr_rz, nk_rz, rho_rz, nr_lam, nk_lam, rho_lam,                &
     &           mn_size, m,n, rmn, zmn, lmn,                                  &
     &           k_pflx, k_grid, l_mfilter_ajax, nrho_ajax, ntheta_ajax,       &
     &            nzeta_ajax,iflag, message)

!-------------------------------------------------------------------------------
! FUNCTION: This is an Interface subroutine that takes input values from e.g.
!           WH_READEQ() and feeds them into the AJAX loading subroutines
!        AJAX_LOAD_RZLAM & AJAX_LOAD_MAGFLUX have
!
!Input:
!  n_eq                -input unit number [-]    (i.e input FILE unit number!!)
!  nrho_ajax           -radial nodes to use in internal AJAX data [-]
!  ntheta_ajax         -poloidal nodes to use in internal AJAX data [-]
!  nzeta_ajax          -toroidal nodes to use in 3-D internal AJAX data [-]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      USE AJAX_MOD
      IMPLICIT NONE

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & nrho_ajax,ntheta_ajax,nzeta_ajax

      INTEGER :: n_eq, i, k

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of local variables
      LOGICAL ::                                                               &
     & l_axisymmetric, l_mfilter_ajax


      INTEGER, INTENT(IN) ::                                                   &
     & k_grid,k_pflx,                                                          &
     & nr_iota,nr_rz, nk_rz, nr_lam, nk_lam, mn_size

      INTEGER, INTENT(IN), DIMENSION(mn_size) :: m, n

      REAL(KIND=rspec) ::                                                      &
     & phitot

      REAL(KIND=rspec), INTENT(IN)::                                           &
     & rho_iota(nr_iota),iotabar(nr_iota),                                     &
     & rho_rz(nr_rz),rho_lam(nr_lam)
 

      REAL(KIND=rspec), DIMENSION(nr_rz, nk_rz), INTENT(IN) ::  rmn, zmn

      REAL(KIND=rspec), DIMENSION(nr_lam, nk_lam), INTENT(IN) ::   lmn


      INTEGER :: ii, jj


!      write(*,*) 'phitot in INTERFACE = ', phitot
!      write(*,*) 'nr_iota in INTERFACE = ', nr_iota
!      write(*,*) 'nr_rz, nk_rz in INTERFACE = ', nr_rz, nk_rz
!      write(*,*) 'nr_lam, nk_lam in INTERFACE = ', nr_lam, nk_lam
!      write(*,*) 'k_pflx in INTERFACE = ', k_pflx
!      write(*,*) 'k_grid in INTERFACE = ', k_grid


!!R,Z and lambda modes
!      READ(n_eq,*) nk_rzl
!
!!Modes and amplitudes
!      ALLOCATE(m(nk_rzl),n(nk_rzl),rmn(1:nr_rzl,1:nk_rzl),                     &
!     &           zmn(1:nr_rzl,1:nk_rzl),lmn(1:nr_rzl,1:nk_rzl))
!
!      !Check to see whether this is a tokamak or stellarator
!      l_axisymmetric=.TRUE.
!
!      DO k=1,nk_rzl
!
!        READ(n_eq,*) m(k),n(k)
!        write(*,*) 'k, m(k), n(k) in READEQ = ', k, m(k), n(k)
!        IF(n(k) /= 0 .AND. l_axisymmetric) l_axisymmetric=.FALSE.
!        READ(n_eq,*) (rmn(i,k),i=1,nr_rzl)
!        READ(n_eq,*) (zmn(i,k),i=1,nr_rzl)
!        READ(n_eq,*) (lmn(i,k),i=1,nr_rzl)
!      ENDDO
!
!      write(*,*) 'l_axisymmetric in INTERFACE = ', l_axisymmetric
!
!      IF(l_axisymmetric) THEN
!
!        !Axisymmetric plasma, calculate lambdas in AJAX on the internal grid
!        iflag=0
!        message=''
!        CALL AJAX_LOAD_RZLAM(nr_rzl,nk_rzl,rho_rz,m,n,rmn,zmn,iflag,           &
!     &                       message,                                          &
!     &                       K_GRID=k_grid,                                    &
!     &                       L_MFILTER_AJAX=l_mfilter_ajax,                    &
!     &                       NRHO_AJAX=nrho_ajax,                              &
!     &                       NTHETA_AJAX=ntheta_ajax)
!
!        !Check messages
!        IF(iflag > 0) THEN
!
!          message='READEQ(1)/'//message
!          GOTO 9999
!
!        ENDIF
!
!      ELSE

        !Non-axisymmetric plasma, use the lambdas from the data file
        iflag=0
        message=''
        CALL AJAX_LOAD_RZLAM(nr_rz,nk_rz,rho_rz,m,n,rmn,zmn,iflag,           &
     &                       message,                                          &
     &                       K_GRID=k_grid,                                    &
     &                       NRHO_AJAX=nrho_ajax,                              &
     &                       NTHETA_AJAX=ntheta_ajax,                          &
     &                       NZETA_AJAX=nzeta_ajax,                            &
     &                       NR_LAM=nr_lam,                                    &
     &                       NK_LAM=nk_lam,                                    &
     &                       RHO_LAM=rho_lam,                                  &
     &                       LAM=lmn)

        !Check messages
        IF(iflag > 0) THEN

          message='READEQ(2)/'//message
          GOTO 9999

        ENDIF

!      ENDIF


      iflag=0
      message=''
      CALL AJAX_LOAD_MAGFLUX(phitot,k_pflx,nr_iota,rho_iota,iotabar,           &
     &                       iflag,message)

      !Check messages
      IF(iflag > 0) message='READEQ(3)/'//message


 9999 CONTINUE

!      DEALLOCATE(rho_iota,iotabar)
!      DEALLOCATE(rho_rz,rho_lam,m,n,rmn,zmn,lmn)

      END SUBROUTINE AJAX_LOAD_INTERFACE

      END MODULE v3read_wout

