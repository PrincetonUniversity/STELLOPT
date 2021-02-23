!-----------------------------------------------------------------------
!     Subroutine:    chisq_gamma_c
!     Authors:       Aaron Bader
!     Date:          04/23/2019
!     Description:   Calculate gamma_c metric for energetic particle 
!                    transport
!                    
!-----------------------------------------------------------------------
      SUBROUTINE chisq_gamma_c(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals
!      USE stel_tools
      USE stellopt_input_mod
!      USE read_boozer_mod
!      USE read_wout_mod, ONLY: ns
      USE equil_utils 
      USE safe_open_mod, ONLY: safe_open
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      REAL(rprec), INTENT(in)    ::  target(nsd)
      REAL(rprec), INTENT(in)    ::  sigma(nsd)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!        lreset_s    Gets set to true if using R,PHI,Z Specification
!        ik          Dummy index
!        ier         Error Flag
!        te_val      Holds profile evaulation
!-----------------------------------------------------------------------
      LOGICAL :: lsym, lgammac_debug
      INTEGER :: dex, ik, i, j, k, ier, jmin, igc2, igc3, igc4, igc5, igc6, igc7, ierrgc2
      INTEGER :: igc8
      REAL(rprec) :: s, rovera, theta, zeta_p, delzeta_p, u, v, coszeta_p, sinzeta_p
      REAL(rprec) :: coszeta_fp, sinzeta_fp

      REAL(rprec) :: s_initA, u_initA, v_initA, s_initB, u_initB, v_initB
      REAL(rprec) :: iota, iotap, minB, maxB, B_refl, psi_a, B_zeta_p
      REAL(rprec) :: X, Y, Xp, Yp, Z, Zp, R, dpsidr, dpsidz, X_fp, Y_fp
      REAL(rprec) :: Bx, By, Bz, Br, Bphi, Bx_fp, By_fp
      REAL(rprec) :: bxn, byn, bzn, bxnp, bynp, bznp, bxn2p, byn2p, bzn2p
      REAL(rprec) :: Bxp, Byp, Bzp, Bx2p, By2p, Bz2p
      REAL(rprec) :: temp, sqrt_bbb, modbp, modbm, del 
      REAL(rprec) :: dIdb, dgdb, dbigGdb, dVdb, vrovervt, dloverb 
      REAL(rprec) :: gamma_c, dsoverb, wellGamma_c, bigGamma_c
      REAL(rprec) :: dBds, dBdr, dBdphi, dBdz
      REAL(rprec) :: phip, phim, Bxt, Byt, Bzt, Brt, Bpt
      REAL(rprec) :: Bxt2, Byt2, Bzt2
      REAL(rprec) :: g,dsdx,dudx,dvdx,dsdy,dudy,dvdy,dsdz,dudz,dvdz
      REAL(rprec) :: sqrtg, jacobian, jac_suvxyz
      REAL(rprec) :: norm_grad_psi_xyz
      REAL(rprec), DIMENSION(3) :: sflCrd, Bxyz, bnxyz, crossnum, crossden
      REAL(rprec), DIMENSION(3) :: Bxyz_fp
      REAL(rprec), DIMENSION(3) :: e_phi, e_r, e_z, grads, gradbtest
      REAL(rprec), DIMENSION(3) :: grads_xyz, gradu_xyz, gradv_xyz
      REAL(rprec), DIMENSION(3) :: grad_zeta_p, grad_psi_x_b, grad_psi_xyz
      REAL(rprec), DIMENSION(3) :: bdotgradb, grad_psi_x_b_xyz, grad_zeta_p_xyz,bdotgradb_sp
      INTEGER, PARAMETER :: gammac_ntransits = 400
      INTEGER, PARAMETER :: gammac_delzetadiv = 400
      INTEGER, PARAMETER :: nsteps = gammac_ntransits*gammac_delzetadiv
      INTEGER, PARAMETER :: gammac_bpstep = 80 !division in b'

      
      real(rprec), DIMENSION(3) :: gradR, gradZ, gradB, gradB_xyz
      real(rprec), DIMENSION(3) :: gradR_init, gradZ_init, gradB_init
      real(rprec), DIMENSION(3) :: grad_psi
      REAL(rprec), DIMENSION(3) :: esubs, esubu, esubv
      REAL(rprec), DIMENSION(3) :: es, eu, ev, es_init, eu_init, ev_init
      real(rprec), DIMENSION(3) :: dxyzdu, dxyzdv, dxyzds

      REAL(rprec), DIMENSION(nsteps) :: ds, modB, dBdpsi, kappa_g
      REAL(rprec), DIMENSION(nsteps) :: grad_psi_norm, grad_psi_i
      REAL(rprec), DIMENSION(nsteps) :: e_theta_norm, e_theta_i
      REAL(rprec), DIMENSION(nsteps) :: dBsupvdpsi, dVdb_t1, dBsupphidpsi
      REAL(rprec), DIMENSION(nsteps) :: dBsupv
      REAL(rprec), DIMENSION(nsteps) :: Bsups, Bsupu, Bsupv
      REAL(rprec), DIMENSION(nsteps, 3) :: binormal
    
      integer, parameter :: maxwells = 5000

      integer, dimension(maxwells) :: well_start, well_stop
      integer :: in_well, cur_well, nwells

      REAL(rprec) :: deltabp, den
      REAL(rprec), DIMENSION(gammac_bpstep) :: bp

      real(rprec) :: grad_psi_min, e_theta_min, cur_Bmin

      REAL(rprec), PARAMETER :: zero   = 0.0_rprec
      REAL(rprec), PARAMETER :: one    = 1.0_rprec
      REAL(rprec), PARAMETER :: two    = 2.0_rprec
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      dex   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'GAMMA_C ',dex,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  GAMMA_C  K'
      IF (niter >= 0) THEN
        pi2 = 2.0_rprec*3.14159265358979_rprec
        delzeta_p = -1.0_rprec/gammac_delzetadiv !step size
        !delzeta = -0.0034433355_rprec !for comparison with ROSE
        psi_a = phiedge !Toroidal flux at the edge

!--------------------------------- DO ik = 1,nsd over each surface
        DO ik = 1,nsd !go through each surface
          IF (sigma(ik) >= bigno) CYCLE

          dloverb = 0.0_rprec
          mtargets = mtargets + 1

          !get normalized toroidal flux
          s = shat(ik)
          !s = 1.0_rprec * (ik) / (ns-1) !vmec labels are in normalized flux 
          rovera = sqrt(s) !rho = r/a   
          
          CALL get_equil_iota(s,iota,ier)  ! get iota
          CALL get_equil_iota(s,iota,ier,iotap) ! get d(iota)/ds
          !CALL EZspline_interp(iota_spl,s,iota,ier)
          !CALL EZspline_derivative(iota_spl,1,s,iotap,ier)
          IF (lgammac_debug == 1) then
            write (iunit_out,*) 's', s
            write (iunit_out,*) 'nfp', nfp
            write (iunit_out,*) 'iota', iota
            write (iunit_out,*) 'iotap', iotap
            write (iunit_out,*) 'psi_a', psi_a
          end if

          !Initialize starting point
          theta = zero; zeta_p = zero
          sflCrd(1) = s
          sflCrd(2) = zero
          sflCrd(3) = zero

          ! Initailziations
          gradR = zero; gradZ = zero; modB = zero
          R = zero; Z = zero
          Bxp = zero; Byp = zero; Bzp = zero;
          Bx2p = zero; By2p = zero; Bz2p = zero;

          if (lgammac_debug == 1) write (iunit_out,*) '--------------------------------------------'
!          CALL safe_open(igc2, ierrgc2, 'gc2.'//trim(proc_string),  &
!                       'replace','formatted')
!          CALL safe_open(igc3, ierrgc2, 'gc3.'//trim(proc_string),  &
!                       'replace','formatted')
!          CALL safe_open(igc4, ierrgc2, 'gc4.'//trim(proc_string),  &
!                       'replace','formatted')
!          CALL safe_open(igc5, ierrgc2, 'gc5.'//trim(proc_string),  &
!                       'replace','formatted')
!          CALL safe_open(igc6, ierrgc2, 'gc6.'//trim(proc_string),  &
!                       'replace','formatted')
!          CALL safe_open(igc7, ierrgc2, 'gc7.'//trim(proc_string),  &
!                       'replace','formatted')
!          CALL safe_open(igc8, ierrgc2, 'gc8.'//trim(proc_string),  &
!                       'replace','formatted')
!          write(igc2,*), 'j rho suv_pest suv_vmec uv_mod_vmec_fix R Z X Y'
!          write(igc3,*), 'j modB Br Bphi Bz Bx By Bsups(j) Bsupu(j) Bsupv(j)'
!          write(igc4,*), 'j gradR_init3 gradZ_init3 gradB_init3 gradR3 gradZ3 gradB3'
!          write(igc5,*), 'j esubu3 esubv3 es_init3 eu_init3 ev_init3 gradS3 grad_psi_norm(j) dBdpsi(j) sqrtg '
!          write(igc6,*), 'j X Y Z bnxyz(1:3) Bxyz(1:3) modB dxyzds(1:3) dxyzdu(1:3) dxyzdv(1:3) grads_xyz(1:3) gradu_xyz(1:3) gradv_xyz(1:3) e_theta_norm binormal(1:3) grad_psi_xyz(1:3) dVdb_t1 dBsupphidpsi dBdpsi sqrtg ds jac_suvxyz'
!          write(igc7,*), 'j bdotgradb(1:3) kappa_g'
!          write(igc8,*), 'j PEST(Rad,Pol,Tor) X Y Z Bx By Bz gradPsi_x gradPsi_y gradPsi_z e_theta_norm dBsupphidpsi'
          !first time through calculate fields and some basic parameters
!--------------------------------- DO j = 1,nsteps over each step along line
          DO j = 1,nsteps
            ! coming into the do-loop, sflCrd contains the PEST coordinate to be considered
            ! these the s, theta and zeta values for the pest coorinate are advanced at the 
            ! end of this do loop. 
           s_initA = sflCrd(1)  ! sflCrd(1) = s
            u_initA = sflCrd(2)  ! sflCrd(2) = theta (pest)
            v_initA = sflCrd(3)  ! sflCrd(3) = zeta (pest) (see above and at end of this do-loop)
            
            ! pest2vmec in: s, theta, phi ('laboratory' 0->2pi, full torus)
            !  on input, phi is converted to 'vmec' phi which goes
            !  from 0->2pi over a field period.
            !           out: s, theta, phi('vmec'
            ! For compatibility with stell_tools/pest2vmec, invert the poloidal and toroidal angles
            sflCrd(2) = -sflCrd(2)
            sflCrd(3) = -sflCrd(3)
            CALL pest2vmec(sflCrd) !convert to VMEC coordinates
            s_initB = sflCrd(1)  ! sflCrd now has VMEC coordinates!!! s=s
            u_initB = sflCrd(2)  ! u (vmec) <- changes to VMEC coords
            v_initB = sflCrd(3)  ! v (vmec) <- goes from 0-> 2pi over a field period 
            u = modulo(u_initB,pi2) ! does nothing in this position (jcs)
            v = modulo(v_initB,pi2) ! does nothing in this position (jcs)
            ! Be careful mixing the two sets of ccordinates! 
            ! the variables 's, theta and zeta' contain the original pest coordinate, as
            ! do s_initA, u_initA, and v_initA

           ! write (*,*) 'jsuv',j,s,u,v
            ! First get all the values you need
            ier = 0
           ! the following function is expecting (s,u,v)_vmec, where
           !  v_vmec goes from 0->2pi over a single field period. 
           !  The equivalent 'Laboratory' Phi = v_vmec/NFP 
            CALL get_equil_RZ(s, u, v, R, Z, ier, gradR_init, gradZ_init)
           ! at this point, gradR, gradZ are calculated for a single field period
           ! derivatives are w.r.t  (U,V, RHO)!!!!
            CALL get_equil_Bflx(s, u, v, Bsups(j), Bsupu(j), Bsupv(j), ier, B_GRAD = gradB_init)
           ! at this point, Bsups, Bsupu, Bsupv, gradB are calculated for a single field period
           ! derivatives are w.r.t rho=sqrt(s),u,v

            !Get B field
            !TODO figure out why exactly the code (copied from get_equil_B) divides
            !by nfp in the Bphi terms only. It seems weird but agrees with ROSE
! comment above may be out of date.
            CALL get_equil_Bcylsuv(s, u, v, Br, Bphi, Bz, ier, modB(j))
            CALL get_equil_Bcylsuv2(s, u, v, ier, dBsupv)
           ! at this point, Br, Bphi, Bz, modB, are calculated for a single field period
            ! Warning - Br and Bphi come from s,u,v conversion.

            ! Now Calc Geometric values
            ! zeta was the same as zeta (pest). Now it is being mixed with R
            ! to make X, Y.
            ! X, Y are full torus, X_fp, Y_fp are single field period
            X=R*cos(-zeta_p)
            !Y=R*sin(zeta_p)  ! JCS / AB Check Sign
            Y=R*sin(-zeta_p)  ! JCS / AB Check Sign
            !X_fp=R*cos(v)
            !Y_fp=R*sin(v)


            ! Convert radial gradients from d/drho to d/ds
            ! s = rho^2; ds/drho = 2*rho; drho/ds = 1/(2*rho) = 0.5/rovera
            ! - use rovera from above
            ! gradX(1) = dX/du, gradX(2) = dX/dv, gradX(3) = dX/dsqrt(s)
            gradB(1) = gradB_init(1)
            gradR(1) = gradR_init(1)
            gradZ(1) = gradZ_init(1)
            gradB(2) = gradB_init(2)
            gradR(2) = gradR_init(2)
            gradZ(2) = gradZ_init(2)
            gradB(3) = gradB_init(3) / (two * rovera)
            gradR(3) = gradR_init(3) / (two * rovera)
            gradZ(3) = gradZ_init(3) / (two * rovera)
           
!new
           ! dBsupphidpsi calculation begins here
           ! B^phi = Toroidal component of B vevtor   psi=toroidal flux / 2pi
           ! see VMEC notes
           dBsupphidpsi(j) = dBsupv(3)*(pi2/psi_a) / (two * rovera)
           ! Following AB's description (Eqs 39-45)
           !


            ! Calc grad(s)
            ! grad(s) = ( (d (X,Y,Z) / du) x (d (X,Y,Z) / dv) ) /  
            ! (d (X,Y,Z) / du dot ( d (X,Y,Z) / dv  x d (X,Y,Z) / ds) 
            ! grad(s) = ((d (R,Phi,Z) / du) x (d (R,Phi,Z) / dv)) /  
            ! (d (R,Phi,Z) / du dot ( d (R,Phi,Z) / dv  x d (R,Phi,Z) / ds) 
            ! grad(s) = ...
            !    (dR/du + dPhi/du +dZ/du) x (dR/dv + dPhi/dv + dZ/dv) / jacobian
            ! jacobian =  (dR/du + dPhi/du +dZ/du) dot 
            !    (dR/dv + dPhi/dv + dZ/dv) x (dR/ds + dPhi/ds + dZ/ds)
            ! esubu = dR/du +dPhi/du + dZ/du
            ! esubv = dR/dv +dPhi/dv + dZ/dv
            ! esubs = dR/ds +dPhi/ds + dZ/ds
            ! Note: dPhi/dv = 1/NFP 
            ! jacobian = e_u . e_v x e_s

           ! ??? JCS check- derivatives are w.r.t s,u,v

            ! Calc grad(s) (copied from txport)
            esubs(1) = gradR(3)    ! dR/ds
            esubs(2) = zero        ! dPhi/ds
            esubs(3) = gradZ(3)    ! dZ/ds

            esubu(1) = gradR(1)    ! dR/du
            esubu(2) = zero        ! dPhi/du
            esubu(3) = gradZ(1)    ! dZ/du

            esubv(1) = gradR(2)    ! dR/dv
            !esubv(2) = one
            esubv(2) = one/nfp         ! dPhi/dv:  v= nfp *mod(phi,2*pi/nfp) -> dv ~ nfp * dphi
            esubv(3) = gradZ(2)    ! dZ/dv

            !esubv(1) = esubv(1)*nfp
            !esubv(3) = esubv(3)*nfp
            !esubs(2) = esubs(2)*R
            !esubu(2) = esubu(2)*R
            !esubv(2) = esubv(2)*R
!  esubv x esubs = (R) esubv(2) * esubs(3) - esubv(3) * esubs(2) +
!                 (Phi) esubv(3) * esubs(1) - esubv(1) * esubs(3) +
!                  (Z) esubv(1) * esubs(2) - esubv(2) * esubs(1)
!              jacobian = esubu(1) * (esubv(2) * esubs(3)) +  &
!                         esubu(2) * (esubv(3) * esubs(1)) +  &
!                         esubu(3) * (-esubv(2) * esubs(1))
              jacobian = esubu(1) * (esubv(2) * esubs(3)) +  &
                         esubu(2) * (esubv(3) * esubs(1)) +  &
                         esubu(3) * (-esubv(2) * esubs(1))

            !sqrtg = R*(gradR(1)*gradZ(3)-gradR(3)*gradZ(1))
            sqrtg = jacobian

            CALL cross_product(esubu,esubv,es_init)
            CALL cross_product(esubv,esubs,eu_init)
            CALL cross_product(esubs,esubu,ev_init)

! JCS   es, eu, ev = e^s,e^u, e^v =  grad(s), grad(u), and grad(v) in cylindrical coordinates
!  (R, phi, Z) for the points on the single field period.
            es = es_init/sqrtg
            eu = eu_init/sqrtg
            ev = ev_init/sqrtg
! gradS is grad psi
            gradS = es * phiedge
            ! grad_psi is for the points on a single field period (not full torus)
            grad_psi = gradS ! in cylindrical coordinates
            ! grad_psi norma is the norm of grad_psi / (2*pi)  JCS - why 2pi?
            grad_psi_norm(j) = sqrt(grad_psi(1)*grad_psi(1) + grad_psi(2)*grad_psi(2) & 
                                  & + grad_psi(3)*grad_psi(3))/pi2

!old`            ! dB/dpsi
!            dBdpsi(j) = gradB(3)*pi2/psi_a !This has been verified with ROSE
!new            ! dB/dpsi
            dBdpsi(j) = gradB(3)*pi2/psi_a !This has been verified with nothing (JCS)

            !Calculate crude arclength
            ! ds gives the distance between the current point and the previous point.
            !   doesn't have a meaning for the first index - will leave it at zero 
            !The integration over arclength is one of the key areas where stellopt/rose differ
            IF (j > 1) THEN
              ! X, and Xp are alright to use in this context (JCS) X_fp, and Y_fp would be wrong 
              ds(j) = sqrt((X-Xp)*(X-Xp) + (Y-Yp)*(Y-Yp) + (Z-Zp)*(Z-Zp))
!              IF (j == 2) ds(1) = ds(j)
              IF (j == 2) ds(1) = 0
            END IF
            ! update Xp, Yp, Zp for next iteration
            Xp = X
            Yp = Y
            Zp = Z
            !integrate dloverb
            dloverb = dloverb + ds(j)/modB(j)

!      dB/ds, dB/dpsi, dB/dr, dB/phi, dB/dz, gradB=dB/d(x,y,z)
! I have gradB = dB/d(s,u,v)
! Later , will want dB/dx, dB/dy, and dB/dz for B dot grad(|B|) = d (bx (x) + by(y) + bz(z)) / d(l)


! JCS The following block is probably needed ? and possibly out-of-date
!           !now form dX/du, dX/dv, dX/ds where X=(x,y,z) 
!           !we note that in the PEST coordinates zeta = cylindrical angle phi
!           !we'll use this relation a lot
!           ! Also note, that the terms above (inc. derivatives)
!           ! were evaluated on the VMEC grid for a single field period
!           ! and may need to be projected/rotataed to account for
!           ! full-torus coordinates

!          coszeta, sinzeta are the cos and sin values of the vmec coordinate
!          v, which goes from 0->2pi on a field period

           ! the above comment is out of date

           ! note the '-' sign
           coszeta_p = cos(-zeta_p) ! dX/dR
           sinzeta_p = sin(-zeta_p) ! dY/dR 
!  go from dR/(rho,u,v) and dZ(rho,u,v) and dPhi(rho,u,v) to
!  dx/d(rho,u,v), dy/d(rho,u,v),  dz/d(rho,u,v),
! on the full torus
! (x) - component of d(xyz)/d(sqrt(s))
           dxyzds(1) = gradR(3)*coszeta_p ! dR/dsqrt(s) * dX/dR=dX/dsqrt(s)
! (y) - component of d(xyz)/d(sqrt(s))
           dxyzds(2) = gradR(3)*sinzeta_p  ! dR/dsqrt(s) * dy/dR = dy/dsqrt(s)
! (z) - component of d(xyz)/d(sqrt(s))
           dxyzds(3) = gradZ(3)

! (x) - components of d(xyz)/d(u)
           dxyzdu(1) = gradR(1)*coszeta_p ! dR/du * dX/dR = dX/du
! (y) - components of d(xyz)/d(u)
           dxyzdu(2) = gradR(1)*sinzeta_p  ! dR/du * dY/dR = dY/du
! (z) - components of d(xyz)/d(u)
           dxyzdu(3) = gradZ(1) 

! (x) - components of d(xyz)/d(v)
! dX/dv = dR/dv * dX/dR + dX/dPhi * dPhi/dv
           dxyzdv(1) = gradR(2)*coszeta_p - R*sinzeta_p/nfp
! (y) - components of d(xyz)/d(v)
! dY/dv = dR/dv * dY/dR + dY/dPhi * dPhi/dv
           dxyzdv(2) = gradR(2)*sinzeta_p + R*coszeta_p /nfp!
! (z) - components of d(xyz)/d(v)
           dxyzdv(3) = gradZ(2)
! Now have d(xyz)/dsqrt(s), d(xyz)/du, and d(xyz)/dv projected in cartestian coordinates (x,y,z)
! for points on the full torus (as a consequence of using sinzeta, coszeta)

! from d(x,y,z)/d(sqrt(s),u,v), generate d(sqrt(s),u,v) / d(x,y,z), aka grad(s), grad(u) and grad(v) in (x,y,z) coords
! JCS the following may not be true anymore. Code has changed
! . No checking has been performed
!           ! These derivatives have been checked and agree with ROSE
!           ! however, the factor of 2*pi needs to be included because ROSE
!           ! parametrizes u and v differently 
!           !write(*,*) 'dxyzdu',dxyzdu(j,1)*pi2,dxyzdu(j,2)*pi2,dxyzdu(j,3)*pi2
!           !write(*,*) 'dxyzdv',dxyzdv(j,1)*pi2,dxyzdv(j,2)*pi2,dxyzdv(j,3)*pi2
!           !write(*,*) 'dxyzds',dxyzds(j,1),dxyzds(j,2),dxyzds(j,3)
!
! need jacobian in _xyz
           ! cross den is first part of denominator
           CALL cross_product(dxyzdv, dxyzds, crossden)
!           !write(*,*) 'crossden', crossden
           jac_suvxyz = -1.0_rprec * dot_product(crossden, dxyzdu) 

           !Now form gradpsi: crossnum is numerator part
           CALL cross_product(dxyzdu, dxyzdv, crossnum)
           !write(*,*) 'crossnum', crossnum
           ! and grad(sqrt(s)) - including
           grads_xyz = crossnum/jac_suvxyz
!           !write(*,*) 'denom',dot_product(crossden, dxyzdu(j,:))
           grad_psi_xyz = grads_xyz * psi_a / pi2  ! JCS added pi2 
           norm_grad_psi_xyz = sqrt(grad_psi_xyz(1)**2 + grad_psi_xyz(2)**2 &
                                    + grad_psi_xyz(3)**2)
!           !write(*,*) 'grads',grads(1),grads(2),grads(3)
!           !grad s agrees with ROSE, no normalizations needed
!           !grad_psi_norm(j) = sqrt(grad_psi(1)*grad_psi(1) + grad_psi(2)*grad_psi(2) & 
!           !                      & + grad_psi(3)*grad_psi(3))/pi2
!
           !Now form gradu: crossnum is numerator part
           CALL cross_product(dxyzdv, dxyzds, crossnum)
           gradu_xyz = crossnum/jac_suvxyz

           !Now form gradv: crossnum is numerator part
           CALL cross_product(dxyzds, dxyzdu, crossnum)
           gradv_xyz = crossnum/jac_suvxyz
      
           

            ! Warning - Br and Bphi come from s,u,v conversion.
            ! Here, the pest zeta is mixed in.
            ! Bxyz contains the B vector on the full torus, Bxyz_fp contain the
            ! B vectors on a sigle field period
            ! CHECK THIS USAGE OF ZETA_P!!!! JCS
            ! zeta_p = -zeta_vmec
            Bx = Br*cos(-zeta_p) - Bphi*sin(-zeta_p)
            By = Br*sin(-zeta_p) + Bphi*cos(-zeta_p)
            !Bx_fp = Br*cos(v) - Bphi*sin(v)
            !By_fp = Br*sin(v) + Bphi*cos(v)
            Bxyz(1) = Bx
            Bxyz(2) = By
            Bxyz(3) = Bz
            bnxyz(1) = Bx / modB(j)
            bnxyz(2) = By / modB(j)
            bnxyz(3) = Bz / modB(j)
            !Bxyz_fp(1) = Bx_fp
            !Bxyz_fp(2) = By_fp
            !Bxyz_fp(3) = Bz
            !write(*,*) 'Bxyz',Bxyz(1),Bxyz(2),Bxyz(3),modB(j), ier
           
            
            
! JCS the 'old' gradB was w.r.t (x,y,z) coordinates
! need a grad_psi.  check for normalization grad_psi ~ grad_s * psi_lcfs
            !Calculate |e_theta|  (cylindrical ?? JCS)
! e_r, e_phi and e_z are expressed in the (x,y,z)-components 
            ! This is the unit vector pointing from the origin to the point on the
            ! surface.  Note: the original version used the full torus coords,
            ! the new version uses the signle-fp coords
            !e_r
!old - is only used for the dot product with grad_psi_xyz
            e_r(1) = X/R
            e_r(2) = Y/R
            e_r(3) = 0.0_rprec
! new
!            e_r(1) = X_fp/R
!            e_r(2) = Y_fp/R
!            e_r(3) = 0.0_rprec

            !e_phi
! old - e_phi is not used
            e_phi(1) = -Y/R
           e_phi(2) = X/R
            e_phi(3) = 0.0_rprec
!new
!            e_phi(1) = -Y_fp/R
!            e_phi(2) = X_fp/R
!            e_phi(3) = 0.0_rprec
! e_z is only used for the dot products with grad_psi_xyz
            !e_z (can move this outside the loop if desired for speed)
            e_z(1) = 0.0
            e_z(2) = 0.0
            e_z(3) = 1.0
            !write(*,*) 'e_r', e_r
            !write(*,*) 'e_z', e_z
            !write(*,*) 'e_phi', e_phi
! new: this is the dot product between vectors and gradients on the 
! field period (not full torus)
! update JCS - grad_psi is in cyl coords, e_r is in (x,y,z)? 
!            dpsidr = dot_product(grad_psi, e_r)
!            dpsidz = dot_product(grad_psi, e_z)
! update JCS - grad_psi_xyz is in (x,y,z) coords, e_r is in (x,y,z)? 
            dpsidr = dot_product(grad_psi_xyz, e_r)
            dpsidz = dot_product(grad_psi_xyz, e_z)
            !B_zeta = dot_product(Bxyz, e_phi)
!JCS I should check dpsidr and dpsidz here

            !Note that e_theta_norm is negative iff Bphi is negative
            !this is the case in QH. However, all this does is multiply the
            !entire expression by -1 of gamma_c, and since the integrated 
            !quantity in Gamma_c is squared, this doesn't matter
            !the division by 2pi is needed to match with rose convention
!old            e_theta_norm(j) = sqrt(dpsidr*dpsidr + dpsidz*dpsidz)/Bphi/pi2
!new  JCS added Abs part. Eqn 41 of nemov has || operator
            e_theta_norm(j) = abs(sqrt(dpsidr*dpsidr + dpsidz*dpsidz)/Bphi)
            !write(*,*) j,'e_theta_norm, B_zeta, Bphi',e_theta_norm(j),Bphi


! JCS e_theta_norm - I did not check the 2pi normalization (or any NFP mod),
!   but I did take out the 2pi

! JCS - I think this can be simplified if the ezspline stuff returned a tensor grad(b) object
! for b dot grad b, but maybe that would require extra over head - but, it could be built 
! into the stel_tools library, as a benchmark.

! At the end, I need binormal(j,:), bdotgradb (intermediate) and kappa_g(j)
            !geodesic curvature
            !calculate bdotgradb = partial b/ partial s
!            bdotgradb = 0.0_rprec
!            kappa_g(j) = 0.0_rprec
!            binormal(j,:) = 0.0_rprec
! step one: binormal vectr on full torux
!  -construct from B vector (normed, in (x,y,z)-coordinates) 
! and grad psi (normed, also in (x,y,z)-coordinates)
! the old method used x,y,z coordinates. for grad_psi and Bxyz
            binormal(j,:) = 0.0_rprec
!old            !CALL cross_product(grad_psi/grad_psi_norm(j)/pi2, Bxyz/modB(j), binormal(j,:))
! new
            CALL cross_product(grad_psi_xyz/norm_grad_psi_xyz, bnxyz, binormal(j,:))


! step two: bdotgradb on full torus
! I have gradB = d|B|/d(sqrt(s),u,v) and Br, BPhi, Bz
! The loop uses x,y,z coordinates - will convert to x,y,z?
! bdotgradb = d(b-vector)/dl - will calculate using central differences for most points
! except for the first and last which will use forward- or backward- differences

            bdotgradb = 0.0_rprec
! step three: kappa_g on full torus
!  -- cross product of bdotgradb and binormal
            !geodesic curvature
            kappa_g(j) = 0.0_rprec

!new
! Use central differences to calculated bdotgradb - 
! Since this is done in a loop, three (3) values need to be calculated.
! The first and last points will use forward/backward differences
           bxn = bnxyz(1)
           byn = bnxyz(2)
           bzn = bnxyz(3)

            IF (j == 2) THEN !handle the first index
              bdotgradb(1) = (bxn - bxnp)/ds(j)
              bdotgradb(2) = (byn - bynp)/ds(j)
              bdotgradb(3) = (bzn - bznp)/ds(j)
              ! set the value of kappa_g for the 1st index
              kappa_g(j-1) = dot_product(bdotgradb, binormal(j-1,:))
!             write(igc7,'(1X,I8,4(2X,E16.8))') 1,&
!                       bdotgradb(1), bdotgradb(2), bdotgradb(3), kappa_g(1)
            END IF 

            IF (j >= 3) THEN ! handles the j-1 index
              bdotgradb(1) = (bxn - bxn2p)/(ds(j-1)+ds(j))
              bdotgradb(2) = (byn - byn2p)/(ds(j-1)+ds(j))
              bdotgradb(3) = (bzn - bzn2p)/(ds(j-1)+ds(j))
              ! set the value of kappa_g for the 2nd thru 2nd-to-last index
              kappa_g(j-1) = dot_product(bdotgradb, binormal(j-1, :))
!             write(igc7,'(1X,I8,4(2X,E16.8))') (j-1),&
!                       bdotgradb(1), bdotgradb(2), bdotgradb(3), kappa_g(j-1)
              IF (j == nsteps) THEN ! handle the last point
                bdotgradb(1) = (bxn - bxnp)/(ds(j))
                bdotgradb(2) = (byn - bynp)/(ds(j))
                bdotgradb(3) = (bzn - bznp)/(ds(j))
              ! set the value of kappa_g for the last index

                kappa_g(j) = dot_product(bdotgradb, binormal(j,:))
!             write(igc7,'(1X,I8,4(2X,E16.8))') j, &
!                       bdotgradb(1), bdotgradb(2), bdotgradb(3), kappa_g(j)
              END IF  
               
            END IF

            ! Keep track off the previous 2 sets of Bx, By, Bz  (B?p, and B?2p)
            Bx2p = Bxp
            By2p = Byp
            Bz2p = Bzp
            Bxp = Bx
            Byp = By
            Bzp = Bz
            ! Keep track off the previous 2 sets of bx, by, bz  (b?p, and b?2p)
            bxn2p = bxnp
            byn2p = bynp
            bzn2p = bznp
            bxnp = bxn
            bynp = byn
            bznp = bzn





            !kappa_g(j) = dot_product(Bxyz, gradB)/modB(j)/modB(j)
            !write (*,*) 'gradB', gradB
            !write (*,*) 'kappa_g', kappa_g(j)







            !The terms that go into gV
!old
!            grad_zeta(1) = -Y/R/R
!            grad_zeta(2) = X/R/R
!new - JCS does grad_zeta = grad_phi?
            grad_zeta_p_xyz(1) = -Y_fp
            grad_zeta_p_xyz(2) = X_fp
            grad_zeta_p_xyz(3) = 0.0_rprec




!new
            !dvdB_t1 is the first term in the brackets of dVdb
            ! = iota' ( grad_psi cross b_hat) dot grad_zeta
            CALL cross_product(grad_psi_xyz, bnxyz, grad_psi_x_b_xyz)
            dVdb_t1(j) = iotap*dot_product(grad_psi_x_b_xyz, grad_zeta_p_xyz)

!NOTE TO SELF; I am here



!  JCS  - check this.  dBsubvdpsi - check which index is appropriate and
! the normalization


!           write(igc2,'(1X,I8,13(2X,E16.8))') j,rovera,s_initA, u_initA, v_initA, s_initB, u_initB, v_initB, u, v,R,Z,X,Y
!            write(igc3,'(1X,I8,9(2X,E16.8))') j,modB(j),Br,Bphi,Bz,Bx,By,Bsups(j),Bsupu(j),Bsupv(j)
!            write(igc4,'(1X,I8,18(2X,E16.8))') j,gradR_init(1),gradR_init(2),gradR_init(3), & 
!                                                 gradZ_init(1),gradZ_init(2),gradZ_init(3), &
!                                                 gradB_init(1),gradB_init(2),gradB_init(3), &
!                                                 gradR(1),gradR(2),gradR(3), & 
!                                                 gradZ(1),gradZ(2),gradZ(3), &
!                                                 gradB(1),gradB(2),gradB(3)
!!
!            write(igc5,'(1X,I8,21(2X,E16.8))') j,esubu(1), esubu(2), esubu(3), &
!                                                 esubv(1), esubv(2), esubv(3), &
!                                                 es_init(1), es_init(2), es_init(3), &
!                                                 eu_init(1), eu_init(2), eu_init(3), &
!                                                 ev_init(1), ev_init(2), ev_init(3), &
!                                                 gradS(1), gradS(2), gradS(3), &
!                                                 grad_psi_norm(j), dBdpsi(j), &
!                                                 sqrtg
! igc6: has full torus X,Y,Z coords, Bvector in (x,y,z) coords). d(xyz)/ds, d(xyz)/du and d(xyz)/dv in (x,y,z) values
! for easy plotting.  Also has gradS, gradU and gradV in (x,y,z) values
!  Also has d(b-unit vector)/dl = (b dot grad)B
!   which is kappa = curvature vectors in (x,y,z) coordintes
!             write(igc6,'(1X,I8,41(2X,E16.8))') j,X,Y,Z,bnxyz(1), bnxyz(2), bnxyz(3), &
!                                                 Bxyz(1), Bxyz(2), Bxyz(3), modB(j), &
!                                                 dxyzds(1),dxyzds(2),dxyzds(3), &
!!                                                 dxyzdu(1),dxyzdu(2),dxyzdu(3), &
!                                                 dxyzdv(1),dxyzdv(2),dxyzdv(3), &
!                                                 grads_xyz(1), grads_xyz(2), grads_xyz(3), &
!                                                 gradu_xyz(1), gradu_xyz(2), gradu_xyz(3), &
!                                                 gradv_xyz(1), gradv_xyz(2), gradv_xyz(3), &
!                                                 e_theta_norm(j), binormal(j,1), binormal(j,2), &
!                                                 binormal(j,3), &
!                                                 grad_psi_xyz(1), grad_psi_xyz(2), grad_psi_xyz(3), &
!                                                 dVdb_t1(j), dBsupphidpsi(j), dBdpsi(j), &
!                                                 sqrtg, ds(j), jac_suvxyz
!            write(igc8,'(1X,I8,14(2X,E16.8))') j, s_initA, u_initA, &
!                       v_initA, X, Y, Z, Bx, By, Bz, grad_psi_xyz(1), &
!                       grad_psi_xyz(2), grad_psi_xyz(3), &
!                       e_theta_norm(j), dBsupphidpsi(j)
!! JCS Moded Y, grad_psi_xyz & e_theta_norm
!

             !advance the step- note, this zeta and theta are the PEST coordinates
            zeta_p = zeta_p + delzeta_p
            theta = theta + (iota * delzeta_p)
            !zeta = modulo(zeta, pi2)
            !theta = modulo(theta, pi2)
            ! reset sflCrd to contain the next PEST coordinate
            sflCrd(1) = s !s should be constant, but in case of numerical precision errors
            sflCrd(2) = theta
            sflCrd(3) = zeta_p
            ! on exit, the sflCrd variable now contains the PEST coordinate of the next point
           END DO
          !------------------------------END DO j = 1,nsteps over each step along line
          ! Print things out to a file for post-analysis
          !CALL safe_open(igc2, ierrgc2, 'gc2.'//trim(proc_string), 'replace','formatted')
!          close(igc2)
!          close(igc3)
!          close(igc4)
!          close(igc5)
!          close(igc6)
!          close(igc7)
!          close(igc8)
          !write(igc2,*), "u,v,R,Z,X,Y,Br,Bphi,Bz,Bsups,Bsupu,Bsupv"
          !do j = 1,nsteps
            !  WRITE(6,'(2X,I3,8(2X,E11.4))')
          !end do

          bigGamma_c = 0.0_rprec
          minB = MINVAL(modB)
          maxB = MAXVAL(modB)

          !Make the bp array
          DO i=1,gammac_bpstep
            bp(i) = 1.0_rprec + (maxB-minB)/minB * (i-0.5_rprec) / gammac_bpstep 
          END DO
          deltabp = (maxB - minB)/minB/gammac_bpstep 
          if (lgammac_debug == 1) then
            write(*,*) 'minB',minB
            write(*,*) 'maxB',maxB
            write(*,*) 'deltabp',deltabp
          end if

          !Go through each value of bp
!          write (*,*) '<---Starting loop i=1,bpstep'
!--------------------------------- DO i = 1,bpstep over each bp step
          DO i=1,gammac_bpstep
            grad_psi_i = 1.0E10_rprec !This term appears in the denominator
            e_theta_i = 0

            !calculate the reflecting field
            B_refl = minB*bp(i)
            !we follow along the line marking the well beginning and the well ends
            !when we have a well, we calculate the minimum of grad_psi (precomputed above)
            !then we set all values of grad_psi_i for that well range to the minimum
            !we also set all the well_mask indices to 1.
            in_well = 0

!            write (*,*) 'beginning bp',i
!            write (*,*) 'Bmin, Bmax, B_refl', minB, maxB, B_refl
            
            well_start = 0 !this is the array of all the indices where wells begin
            well_stop = 0 !this is the array of all the indices where wells end
            cur_well = 1 !this is the index of the current well

            !these values are calculated at the well minima
            grad_psi_min = 1.0E10_rprec
            e_theta_min = 0
            cur_Bmin = B_refl !this should be the maximum value in any well

!            write (*,*) '<---Starting loop j=1,nsteps'
!--------------------------------- DO j = 1,nsteps 
            DO j=1,nsteps
              ! Not in a well and shouldn't be 
              IF ((in_well == 0) .and. (modB(j) > B_refl)) CYCLE

              ! We are in a well, but we just passed a boundary value
              IF ((in_well == 1) .and. (modB(j) > B_refl)) THEN
      
                ! reset the well status
                in_well = 0
            
                ! In this case we have exited a well that hasn't been started properly
                ! so we ignore it
                IF (well_start(cur_well) == 0) THEN
                  cur_Bmin = B_refl
                  grad_psi_min = 1.0E10_rprec
                  e_theta_min = 0
                  CYCLE
                END IF

                !mark the well end index
                well_stop(cur_well) = j-2

                !set gradpsi and etheta for all values
                grad_psi_i(well_start(cur_well):well_stop(cur_well)) = grad_psi_min
                e_theta_i(well_start(cur_well):well_stop(cur_well)) = e_theta_min

                !test output
!                write(*,*) 'exiting well',cur_well
!                write(*,*) 'well start,stop',well_start(cur_well), well_stop(cur_well)
!                write(*,*) 'minima',grad_psi_min, e_theta_min
                !print out the Bmin
!                write (*,*) "Bmin, zeta for well ",cur_well, cur_Bmin, zeta_p*4.0_rprec/pi2
         

                !reset minimum of modB, psi and e
                cur_Bmin = B_refl
                e_theta_min = 0
                grad_psi_min = 1.0E10_rprec

                cur_well = cur_well + 1
                
                !exit if we have gone through too many wells
                IF (cur_well >= maxwells) EXIT
                CYCLE

              END IF   
              ! We are not in a well, but are entering one
              IF (in_well == 0 .and. modB(j) < B_refl) THEN
                

                !add the beginning value, but only if we're not at the start
                IF (j > 1) THEN
                  well_start(cur_well) = j+1
                END IF 
                !mark that we're in a well
                in_well = 1
              END IF

              ! We are in a well and should be in one
              ! Note this gets executed if the previous IF statement is true
              IF (in_well == 1 .and. modB(j) < B_refl) THEN

                !check modB and update minima
                IF (modB(j) < cur_Bmin) THEN
                  cur_Bmin = modB(j)
                  e_theta_min = e_theta_norm(j)
                  grad_psi_min = grad_psi_norm(j)
                  zeta_p = j*delzeta_p
                  jmin = j
                END IF
              END IF
            
            END DO
!------------------------------END DO j = 1,nsteps 
!            write (*,*) '<---Finished loop j=1,nsteps'
            !If we ended in the middle of a well, we'll have one more well_start
            !than well stop, let's correct that
            IF (well_start(cur_well) .NE. 0) THEN
              well_start(cur_well) = 0
            END IF 
            
            nwells = cur_well - 1
!            write(*,*) 'well number',nwells
!            DO k = 1,nwells
!              write(*,*) 'well',k,well_start(k), well_stop(k)
!            END DO 
!            DO j = 1,nsteps
!              write(*,*) j,B_refl,modB(j),grad_psi_norm(j),grad_psi_i(j),e_theta_norm(j),e_theta_i(j)
!            END DO

            !We've assembled all the info we need to compute the major quantities
            !for each well, we integrate the various quantities, dgdb, dGdb, dIdb, and dVdb
            gamma_c = 0.0_rprec
            vrovervt = 0.0_rprec
            wellgamma_c = 0.0_rprec
!            write (*,*) '<---Starting loop k=1,nwells'
!----------------------------- DO k = 1,nwells
            DO k = 1,nwells

              dIdb = 0.0_rprec
              dgdb = 0.0_rprec
              dbigGdb = 0.0_rprec
              dVdb = 0.0_rprec
!---------------------------- DO j = well_start(k),well_stop(k)
              DO j = well_start(k),well_stop(k)
                !double check that we're in a valid well
                if (grad_psi_i(j) > 1.0E8_rprec) CYCLE
                if (e_theta_i(j) == 0.0_rprec) CYCLE

                sqrt_bbb = sqrt(1 - modB(j)/B_refl)

                !dIdb
                temp = ds(j)/2.0_rprec/minB/bp(i)/bp(i) / sqrt_bbb
                dIdb = dIdb + temp
                
                !dgdb
                temp = ds(j) * grad_psi_norm(j) * kappa_g(j) 
                temp = temp/bp(i)/bp(i)/2.0_rprec/modB(j)
                temp = temp*(sqrt_bbb + 1.0_rprec/sqrt_bbb)
                dgdb = dgdb + temp

                !dbigGdb 
                temp = dBdpsi(j) *ds(j) /B_refl / bp(i) / modB(j) / 2.0_rprec
                temp = temp*(sqrt_bbb + 1.0_rprec/sqrt_bbb)
                dbigGdb = dbigGdb + temp

                !dVdb
                !temp = dVdb_t1(j) - (2.0_rprec * dBdpsi(j) - modB(j)/Bsupv(j)*dBsupvdpsi(j)) 
                temp = dVdb_t1(j) - (2.0_rprec * dBdpsi(j) - modB(j)/Bsupv(j)*dBsupphidpsi(j)) 
                temp = temp * 1.5_rprec * ds(j) / modB(j) / B_refl * sqrt_bbb
                dVdb = dVdb + temp
              

              END DO !end integration over a single well
!-------------------------END DO j = well_start(k),well_stop(k)


              !vrovervt ratio of radial to poloidal drifts
!              write (*,*) '---------------------------------------'
!              write (*,*) 'well k,b ', k, B_refl,well_start(k),well_stop(k)
!              write (*,*) 'dIdb ', dIdb
!              write (*,*) 'dgdb ' , dgdb
!              write (*,*) 'dbigGdb', dbigGdb
!              write (*,*) 'dVdb ', dVdb
!              write (*,*) 'etheta0 ', e_theta_i(j-1)
!              write (*,*) 'grad_psi ', grad_psi_i(j-1)
              IF (well_start(k) < well_stop(k)) THEN

                temp = dgdb/grad_psi_i(well_start(k))/dIdb / minB / e_theta_i(well_start(k))
                temp = temp / (dbigGdb/dIdb + 0.666666_rprec * dVdb/dIdb)
                vrovervt = temp
              ELSE
                vrovervt = 0.0
              END IF
!              write (*,*) 'vrovervt ', vrovervt

              gamma_c = 4.0_rprec/pi2 * atan(vrovervt)
              wellGamma_c = wellGamma_c + (gamma_c * gamma_c * dIdb)
!              write (*,*) 'wellGamma_c ', wellGamma_c
!              write (*,*) 'gamma_c, wellGamma_c', gamma_c, wellGamma_c
            END DO !end sum over all wells           
!--------------------------END DO k = 1,nwells
!            write (*,*) '<---finished loop k=1,nwells'
            bigGamma_c = bigGamma_c + wellGamma_c * pi2/4.0_rprec/sqrt(2.0_rprec)*deltabp
          END DO !end integration over bp
!------------------------------END DO i = 1,bpstep over each bp step
!          write (*,*) '<---Finished loop i=1,bpstep'
          bigGamma_c = bigGamma_c/dloverb
          if (lgammac_debug == 1) then
            write(*,*) 'dloverb',dloverb
            write(*,*) 'bigGamma_c',bigGamma_c
          end if

          vals(mtargets) = bigGamma_c
          sigmas(mtargets) = ABS(sigma(ik))
          targets(mtargets) = target(ik)
          IF (iflag ==1) THEN
            !WRITE(iunit_out,'(3ES22.12E3,3(1X,I5))') targets(mtargets),sigmas(mtargets),vals(mtargets),ik
            WRITE(iunit_out,'(3ES22.12E3,3(1X,I5))') targets(mtargets),sigmas(mtargets),vals(mtargets),ik
          END IF

          
        END DO !end loop over flux surfaces
!------------------------------END DO ik = 1,nsd over each surface
      ELSE !This is the initialization loop that just counts targets
        DO ik = 1, nsd
          IF (sigma(ik) < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) THEN
              target_dex(mtargets)=jtarget_gamma_c
            END IF
          END IF
        END DO
      END IF


      END SUBROUTINE chisq_gamma_c
