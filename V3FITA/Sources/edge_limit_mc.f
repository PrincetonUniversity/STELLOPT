!*******************************************************************************
!  File edge_limit_mc.f
!  Contains module edge_limit_mc

!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine edge_limit_mc_model_compute
!  Both an edge_limit_desc and a model are needed for this computation.

!*******************************************************************************
!  MODULE edge_limit_mc
!    (edge_limit - Model Compute Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
      MODULE edge_limit_mc

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations, constants, utilities
!-------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      USE v3_utilities

!-------------------------------------------------------------------------------
!  MDDC Derived Types
!-------------------------------------------------------------------------------
      USE edge_limit_T
      
!-------------------------------------------------------------------------------
!  Model Derived Types
!-------------------------------------------------------------------------------
      USE model_T

!-------------------------------------------------------------------------------
!  Module for limiter_iso function
!-------------------------------------------------------------------------------
      USE limiter_iso_T

      IMPLICIT NONE
      
!*******************************************************************************
! SECTION II.  INTERFACE BLOCKS
!*******************************************************************************

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      CONTAINS
          
!*******************************************************************************
! SECTION III.  MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Compute an edge_limit signal
!
!    Information comes from the edge_limit_desc and the model
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    s_type = geometric
!      g_type = edge_limit
!        edge limit
!        signal_model_compute_edge_limit
!           % data(1) - limiter function - zero on limiter, negative inside plasma
!           % data(2) - R position of maximum limiter function
!           % data(3) - phi position of maximum limiter function
!           % data(4) - Z position of maximum limiter function
!
!-------------------------------------------------------------------------------

      SUBROUTINE edge_limit_mc_model_compute(a_el_desc,a_model,                &
     &   arr_data,arr_sigma)
!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  Declare Arguments 
      TYPE (edge_limit_desc), INTENT (inout), TARGET :: a_el_desc
      TYPE (model), INTENT (inout), TARGET :: a_model
      REAL(rprec), POINTER, DIMENSION(:) :: arr_data, arr_sigma
!  Note - convention is that arr_data and arr_sigma
!    1) are deallocated in this subroutine, and then
!    2) are allocated in this subroutine

!  Declare local variables
      TYPE (eq_state), POINTER             :: state => null()
      TYPE (edge_limit_desc), POINTER      :: my_el_desc => null()
      TYPE (limiter_iso), POINTER          :: my_iso => null()
      REAL(rprec), DIMENSION(:,:), POINTER :: r_one_uv => null()
      REAL(rprec), DIMENSION(:,:), POINTER :: z_one_uv => null()
      
      REAL(rprec) :: fval_max, fval
      REAL(rprec), DIMENSION(3) :: rpz_at_max, rpz_edge
      INTEGER :: n_data
      INTEGER :: ist, jumax, kvmax, j, k, ier1, istat1

      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'edge_limit_mc_model_compute: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

!  Point state to the model equilibrium state
      state => a_model % eqstate
! 
!  Point to edge_limit
      my_el_desc => a_el_desc

      SELECT CASE (TRIM(ADJUSTL(my_el_desc % edge_limit_type)))
      CASE ('iso_fun')
         my_iso => my_el_desc % lim_iso
      
!  Get the s=1 surface information
         CALL edge_limit_mc_aux_rz_one(state,my_iso % numin,                   &
     &      my_iso % vgrid,r_one_uv,z_one_uv)

!  Allocate data and sigma
         n_data = 4
         IF (ASSOCIATED(arr_data)) THEN
            DEALLOCATE(arr_data,arr_sigma, stat=istat1)
            CALL assert_eq(0,istat1,sub_name // 'arr_data, s dealloc')
         ENDIF
         ALLOCATE(arr_data(n_data),arr_sigma(n_data),stat=istat1)
         CALL assert_eq(0,istat1,sub_name // 'arr_data, sigma alloc')

!----------------------------------------------------------------------
!-- Loop over toroidal planes                                        --
!----------------------------------------------------------------------
!  loop around the s = 1 surface
!  Evaluate limiter function for each point
!  take maximum, store R, Phi, and Z at location of maximum

         fval_max = -1.e10
         rpz_at_max = (/ zero, zero, zero /)
      
         jumax = SIZE(r_one_uv,1)
         kvmax = SIZE(my_iso % vgrid,1)
         CALL assert_eq(kvmax,SIZE(r_one_uv,2),
     &       sub_name // 'Array size mismatch, r_one_uv, vgrid')

         TOR_PLANE: DO k = 1, kvmax

            DO j = 1, jumax   ! loop over poloidal positions
               rpz_edge(1) = r_one_uv(j,k)
               rpz_edge(2) = my_iso % vgrid(k)
               rpz_edge(3) = z_one_uv(j,k)
               CALL limiter_iso_f_cyl(my_iso,rpz_edge,fval)
               IF (fval .gt. fval_max) THEN
                  fval_max = fval
                  rpz_at_max = rpz_edge
               ENDIF
            END DO

         END DO TOR_PLANE
         
         IF (my_el_desc % l_on_edge) THEN
            arr_data(1) = fval_max
         ELSE
            arr_data(1) = MAX(fval_max,zero)
         ENDIF
         arr_data(2:4) = rpz_at_max
         arr_sigma(1:4) = 1.E-10

      CASE ('polygon')
         CALL err_fatal(sub_name // 'This type NYI: ',                         &
     &      char=my_el_desc % edge_limit_type)

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized edge_limit_type: ',          &
     &      char=my_el_desc % edge_limit_type)

      END SELECT ! Different coding depending on edge_limit_type

!  Deallocate the _one_uv pointers
      IF (ASSOCIATED(r_one_uv)) THEN
         DEALLOCATE(r_one_uv,z_one_uv,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'deallocate error')
      END IF
       
      RETURN
      END SUBROUTINE edge_limit_mc_model_compute

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE edge_limit_mc_aux_rz_one(state,numin,vgrid,r_one_uv,          &
     &   z_one_uv)
!
!  Subroutine to compute R and Z on the s=1 surface, for various values 
!  of u and v.
!  Called by the subroutine edge_limit_mc_model_compute
!  First version 30 January 2009
!  Revised 2009-03-03
!  Based on similar coding in eq_aux_2_construct

!  Arguments
!    state     An eq_state. (in)
!    numin     Minimum number of evenly spaced poloidal angles. (in)
!    vgrid     array of toroidal angles. (in)
!    r_one_uv  Array of R values. (out) Allocated here. Deallocation must be 
!      done elsewhere, presumably in the calling subroutine
!    z_one_uv  Pointer to Z values. (out) Allocated here. Deallocation must be 
!      done elsewhere, presumably in the calling subroutine

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_state), INTENT (inout), TARGET             :: state
      INTEGER, INTENT(in)                                 :: numin
      REAL(rprec), DIMENSION(:), INTENT(in)               :: vgrid
      REAL(rprec), DIMENSION(:,:), POINTER                :: r_one_uv,         &
     &   z_one_uv

!  Pointers to components
      TYPE (eq_param_fix), POINTER  :: fixp => null()
      TYPE (eq_aux_1), POINTER :: aux1 => null()

!  Declare local variables
      INTEGER :: ns, ju, kv, n_field_periods
      INTEGER :: mnmax
      LOGICAL :: lasym
      REAL(rprec)    :: delu
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: cosz, sinz
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosu, sinu,                  & 
     &                                            cosv, sinv
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: ugrid

      INTEGER :: i, j, k, js             
      INTEGER, DIMENSION(3)   :: dims           

      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'edge_limit_mc_aux_rz_one: '

!----------------------------------------------------------------------
!-- Start of executable code                                         --
!----------------------------------------------------------------------
!  Set Pointers, DEFINE the aux1 if needed
      fixp => state % fixp
      
      IF (.not. state % l_def_aux1) THEN
         CALL eq_aux_1_define(state)
      END IF
       
      aux1 => state % aux1

!  local variable definitions
!  Consistency checks need to be done
      lasym = fixp % lasym
!  Check to see that aux1 is in correct state.
      ns = SIZE(aux1 % rmnc,2)
      ju = MAX(numin,ns,2 * fixp % mpol + 6)
      kv = SIZE(vgrid,1)
      n_field_periods = fixp % nfp
      mnmax = size(aux1 % xm,1)
      delu = twopi / ju
!----------------------------------------------------------------------
!-- allocate arrays                                                  --
!----------------------------------------------------------------------
!  u VMEC flux grid
      ALLOCATE (ugrid(ju),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc ugrid')

!  non-nyq variables
      ALLOCATE (r_one_uv(ju,kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc rsuv')

      ALLOCATE (z_one_uv(ju,kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc zsuv')

!  Sin and cosine arrays
      ALLOCATE (cosz(mnmax), sinz(mnmax),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc cosz')

      ALLOCATE (cosu(mnmax, ju), sinu(mnmax, ju),
     &          cosv(mnmax, kv), sinv(mnmax, kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc cosu')

!----------------------------------------------------------------------
!-- set up (s,u,v) VMEC flux grids                                   --
!----------------------------------------------------------------------
!  index ns corresponds to s = 1
      DO i = 1, ju
         ugrid(i) = (i - 1) * delu
         cosu(:,i) = COS(aux1 % xm * ugrid(i))
         sinu(:,i) = SIN(aux1 % xm * ugrid(i))
      END DO
      DO i = 1, kv
!  vgrid comes n as an argument - values of toroidal angle
         cosv(:,i) = COS(aux1 % xn * vgrid(i))
         sinv(:,i) = SIN(aux1 % xn * vgrid(i))
      END DO
!
!       DO FOURIER MODE SUM 
!
      DO j = 1, ju
      DO k = 1, kv
!          zetamn = xm*ugrid(j) - xn*vgrid(k)
!          cosz = COS(zetamn)
!          sinz = SIN(zetamn)
         cosz = cosu(:,j) * cosv(:,k) + sinu(:,j) * sinv(:,k)
         sinz = sinu(:,j) * cosv(:,k) - cosu(:,j) * sinv(:,k)
         r_one_uv(j,k) = SUM(aux1 % rmnc(:,ns) * cosz)
         z_one_uv(j,k) = SUM(aux1 % zmns(:,ns) * sinz)

!----------------------------------------------------------------------
!-- stellarator asymmetric terms                                     --
!----------------------------------------------------------------------
         IF (lasym) THEN
            r_one_uv(j,k) = r_one_uv(j,k) + SUM(aux1 % rmns(:,ns)*sinz)
            z_one_uv(j,k) = z_one_uv(j,k) + SUM(aux1 % zmnc(:,ns)*cosz)
         ENDIF
      END DO
      END DO

!  Deallocate space and nullify pointers
      DEALLOCATE (cosz, sinz, cosu, cosv, sinu, sinv)
      aux1 => null()
      fixp => null()

      RETURN
      
      END SUBROUTINE edge_limit_mc_aux_rz_one

!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JDH 2009-03-03. First version of module. Based on coding from previous
!    module signal_model
!
!  2009-03-14 JDH
!     Fixed so that aux1 is DEFINED, if need be.
!
      END MODULE edge_limit_mc
