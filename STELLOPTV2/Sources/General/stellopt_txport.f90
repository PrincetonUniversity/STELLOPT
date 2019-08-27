!-----------------------------------------------------------------------
!     Subroutine:    stellopt_txport
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/18/2013
!     Description:   This subroutine calculated turbulent transport.
!                    Calculation of geometric coefficients is provided by
!  =================================================================================
!  =========        Geometry Interface for Stellarators and Tokamaks       =========
!  =========          (P.Xanthopoulos, W.A.Cooper, and Yu.Turkin)          =========
!  =========          pax@ipp.mpg.de  http://www.ipp.mpg.de/~pax/          =========
!  =========    Xanthopoulos et al. Physics of Plasmas 16, 082303 (2009)   =========
!  =================================================================================
!
!                    Calculation of turbulent transport proxy functions by
!  =================================================================================
!  =========              Turbulent Transport Proxy Functions              =========
!  =========                  (H.Mynick, P.Xanthopoulos)                   =========
!  =========           mynick@pppl.gov  http://w3.pppl.gov/~mynick/        =========
!  =========            Mynick Physics of Plasmas 13, 058102 (2006)        =========
!  =================================================================================
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_txport(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE stellopt_targets
      USE equil_vals
      USE equil_utils
      USE booz_persistent
      USE read_boozer_mod
      USE EZspline_obj
      USE EZspline
!DEC$ IF DEFINED (GENE)
      USE gene_subroutine, ONLY: rungene
      USE par_in, ONLY: diagdir,file_extension
      USE parameters_IO, ONLY: read_parameters, write_parameters, PAR_OUT_FILE
      USE coordinates, ONLY: kymin
      USE geometry, ONLY: geomdir, magn_geometry, geomfile
      USE file_io, ONLY: erase_stop_file
      USE mpi_params
      USE discretization, ONLY: mype_gl, n_procs_sim
      USE communications, ONLY: omp_level
      USE mpi_inc
!DEC$ ENDIF
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      LOGICAL :: file_test_eig,file_test_omega, lscreen_gene
      INTEGER :: maxPnt, nalpha0_, ialpha, i, iunit, ik, pol_turns_gene, &
                 iunit_curv, irtmx, irt, iunit_eigen, nu, nv, spawn_procs, &
                 spawn_root, spawn_info, spawn_comm, spawn_intercomm, ier
      REAL(rprec) :: a, s, Ba, Fa, iot,iotp,qprim,shat,&
                     pval, pprime, dtheta, dalpha, alpha0_start_, phi0, &
                     th, jac1, c, &
                     gaa, gat, gst, gsa, gss, alpha, thetastar, &
                     dBds, dBda, q, dpdx, sqrtg, modb, ztemp, &
                     Rp, Rm, Pp, Pm, Zp, Zm, Bp, Bm, u, v, R, Z, &
                     rky, rkalf, sgm2, rkpar0sq, rkparsq, dln_jac_rbb2, &
                     wstare, wde, bs, trm1a, trm1b, trm2a, trm3a, trm4a, &
                     trm6a, zt, dlzt, cj, gene_eig_r, gene_eig_i, gene_kymin,&
                     temp1,temp2,temp3, abserr,alpha0_end,alpha0_start
      COMPLEX :: ctrm1b,ccof0,cbb2a,cwa, cimag
      REAL(rprec), DIMENSION(3) :: sflCrd0,sflCrd, sflCrd_sav, gradS,gradThetaStar,&
                                   gradPhi,mag,gradAlpha, wrk, gradB,R_grad,Z_grad,&
                                   esubs, esubu, esubv, es, eu, ev, gradlam, ea, et
      REAL(rprec),ALLOCATABLE :: phi0_gl(:), cosmt(:),cosnt(:),sinmt(:),sinnt(:), temp_arr(:), tout_arr(:)
      REAL(rprec),ALLOCATABLE :: g11(:,:),g12(:,:),g22(:,:),Bhat(:,:),abs_jac(:,:),&
                                 L1(:,:),L2(:,:),dBdt(:,:),rll(:,:),rllm(:,:),vrbl(:,:), &
                                 kp1(:,:), sloc(:,:), rkp1av(:,:), slocav(:,:), rkpn(:,:), &
                                 rkpg(:,:),rkk(:,:), rllmp(:,:), rllmpp(:,:), vv(:,:), qqfac(:,:),&
                                 dkp(:,:), dkpfac(:,:), vkp1fac(:,:), vqqprox(:,:)
      DOUBLE PRECISION        :: cof(0:3),wa(3,3),wk(8)
      CHARACTER(256) :: temp_str,num_str,gene_str,gene_dir_str,ch_in,gist_filename
      CHARACTER(256) :: command, argv(3)

      ! For J. Proll Proxy
      REAL(rprec) :: bmax, bmin, b0, bamp
      DOUBLE PRECISION :: result,abserr2,A2,B2
      INTEGER :: neval,last
      INTEGER, PARAMETER :: KEY = 3, NW=1024
      DOUBLE PRECISION, PARAMETER :: EPSABS=1.0D-9, EPSREL=1.0D-3
      DOUBLE PRECISION, DIMENSION(nw) :: work
      DOUBLE PRECISION, DIMENSION(nw) :: iwork

      
      REAL(rprec), PARAMETER :: qpthk  = 0.0_rprec
      REAL(rprec), PARAMETER :: dl_z   = 0.2073_rprec
      REAL(rprec), PARAMETER :: rhsrkp = 0.002_rprec
      REAL(rprec), PARAMETER :: tau_s  = 1.1267_rprec
      REAL(rprec), PARAMETER :: rkp_p  = 3.0_rprec
      REAL(rprec), PARAMETER :: rkp_cr = 5.335D-02
      REAL(rprec), PARAMETER :: dk     = 0.9596_rprec
      REAL(rprec), PARAMETER :: tau    = 1.0000_rprec
      REAL(rprec), PARAMETER :: fft    = 0.0_rprec
      REAL(rprec), PARAMETER :: mm     = 5.0_rprec
      REAL(rprec), PARAMETER :: rkymin = 0.05_rprec
      REAL(rprec), PARAMETER :: omn    = 5.0_rprec
      REAL(rprec), PARAMETER :: omti   = 3.0_rprec
      REAL(rprec), PARAMETER :: omte   = 3.0_rprec
      REAL(rprec), PARAMETER :: rr0    = 4.5_rprec
      REAL(rprec), PARAMETER :: rkparfrac = 1.0_rprec
      REAL(rprec), PARAMETER :: rhstar = 0.01_rprec
      REAL(rprec), PARAMETER :: zt_k   = 0.00_rprec
      REAL(rprec), PARAMETER :: dl_Ror = 0.01_rprec
      REAL(rprec), PARAMETER :: zero   = 0.0_rprec
      REAL(rprec), PARAMETER :: one    = 1.0_rprec
      REAL(rprec), PARAMETER :: two    = 2.0_rprec
      
      EXTERNAL C02AFG
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      abserr = 1.0E-9
!DEC$ IF DEFINED (TXPORT_OPT)
      IF (lscreen) WRITE(6,'(a)') ' --------------------  TURBULENT TRANSPORT CALC  ----------------------'
      CALL FLUSH(6)
      IF (ALLOCATED(txport_q)) DEALLOCATE(txport_q)
      IF (ALLOCATED(txport_q_all)) DEALLOCATE(txport_q_all)
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','parvmec','paravmec','vboot','vmec2000_oneeq')
            !IF (lscreen) WRITE(6,'(A)') '     s    alpha0   iota    shear    Bref    Lref  Q_turb'
            maxPnt = nz_txport
            dtheta = pi2/maxPnt
            nalpha0_=nalpha_txport
            IF (lglobal_txport) THEN
               alpha0_start = +pi2/(2*nfp)
               alpha0_end   = -pi2/(2*nfp)
            ELSE
               ! Old
               !alpha0_start  = 0 
               !alpha0_end    = +pi2/(2*nfp)
               ! New
               alpha0_start = alpha_start_txport
               alpha0_end = alpha_end_txport
               IF (lasym)  alpha0_start = -pi2/(2*nfp) ! Get full flux tube
            END IF
            alpha0_start_ = alpha0_start
            dalpha = (alpha0_end - alpha0_start) / REAL(nalpha0_-1,rprec)
            IF (nalpha_txport == 1) dalpha = zero
            ALLOCATE(phi0_gl(1:nalpha0_),g11(nalpha0_,maxPnt),g12(nalpha0_,maxPnt),g22(nalpha0_,maxPnt),&
                     Bhat(nalpha0_,maxPnt),abs_jac(nalpha0_,maxPnt),L1(nalpha0_,maxPnt),L2(nalpha0_,maxPnt),&
                     dBdt(nalpha0_,maxPnt),kp1(nalpha0_,maxPnt))
            ALLOCATE(temp_arr(maxPnt),tout_arr(maxPnt))
            ALLOCATE(rll(nalpha0_,maxPnt),rllm(nalpha0_,maxPnt),vrbl(nalpha0_,maxPnt),sloc(nalpha0_,maxPnt),&
                     rkp1av(nalpha0_,maxPnt),slocav(nalpha0_,maxPnt),rkpn(nalpha0_,maxPnt),rkpg(nalpha0_,maxPnt), &
                     rkk(nalpha0_,maxPnt),rllmp(nalpha0_,maxPnt),rllmpp(nalpha0_,maxPnt),vv(nalpha0_,maxPnt),&
                     dkp(nalpha0_,maxPnt),dkpfac(nalpha0_,maxPnt),vkp1fac(nalpha0_,maxPnt),qqfac(nalpha0_,maxPnt),&
                        vqqprox(nalpha0_,maxPnt))
            IF (ALLOCATED(txport_q)) DEALLOCATE(txport_q)
            IF (ALLOCATED(txport_q_all)) DEALLOCATE(txport_q_all)
            ALLOCATE(txport_q(nrad,nalpha0_,maxPnt))
            IF (TRIM(txport_proxy(1:3)) == 'all')  ALLOCATE(txport_q_all(6,nrad,nalpha0_,maxPnt))
            DO ik = 1, nrad
               IF (sigma_txport(ik) >= bigno) CYCLE
               ! Set Some defaults
               s = s_txport(ik)
               Fa = phiedge
               a = Aminor
               Ba = ABS(Fa/(pi*a**2))
               CALL EZspline_interp(iota_spl,s,iot,iflag)
               CALL EZspline_derivative(iota_spl,1,s,iotp,iflag)
               !iot = iota_b(ik)
               !iotp = 0.5*(iota_b(ik+1)-iota_b(ik-1))/drho
               q = one/iot
               qprim = -iotp*q**2
               shat = two*s/q*qprim
               CALL get_equil_p(s,pval,iflag,pprime)
               dpdx = -4.0_rprec*sqrt(s)/Ba**2 * pprime*mu0
               sflCrd0(1) = s
               sflCrd0(2) = zero
               !dtheta = pi2*pol_turns_gene/maxPnt
               phi0 = alpha0_start_
               ! Open the output file
               IF (lglobal_txport) THEN
                  temp_str = 'gist_geney_'
               ELSE
                  temp_str = 'gist_genet_'
               END IF
               !ALLOCATE(cosmt(mnboz_b),cosnt(mnboz_b),sinmt(mnboz_b),sinnt(mnboz_b))
               IF (lglobal_txport) THEN
                  WRITE(num_str,'(1(A,I3.3))') '_',ik
                  iunit = 400
                  gist_filename = TRIM(TRIM(temp_str)//TRIM(proc_string)//TRIM(num_str))
                  CALL safe_open(iunit,iflag,TRIM(gist_filename),'unknown','formatted')
                  iunit_curv = iunit+1
                  CALL safe_open(iunit_curv,iflag,TRIM('curv_stellopt_'//TRIM(proc_string)//TRIM(num_str)),'unknown','formatted')
               END IF
               ! main Loop
               DO ialpha = 1, nalpha0_
                  IF (.not. lglobal_txport) THEN
                     WRITE(num_str,'(2(A,I3.3))') '_',ik,'_',ialpha
                     iunit = 400
                     gist_filename = TRIM(TRIM(temp_str)//TRIM(proc_string)//TRIM(num_str))
                     CALL safe_open(iunit,iflag,TRIM(gist_filename),'unknown','formatted')
                     iunit_curv = iunit+1
                     CALL safe_open(iunit_curv,iflag,TRIM('curv_stellopt_'//TRIM(proc_string)//TRIM(num_str)),'unknown','formatted')
                  END IF
                  phi0_gl(ialpha) = alpha0_start_ - (ialpha-1)*dalpha
                  phi0 = phi0_gl(ialpha)
                  IF (lglobal_txport .and. ialpha==1) THEN
                     WRITE(iunit,'(A)') '&PARAMETERS'
                     WRITE(iunit,"(A,3F12.7,I5)") "!s0, alpha0_start, alpha0_end, nalpha0 = ", &
                                        s,alpha0_start, alpha0_end,nalpha0_
                     WRITE(iunit,"(A,I3)") "n0_global = ",nfp
                     WRITE(iunit,"(A,F6.3)") "s0 =",s
                     WRITE(iunit,"(A,I4)") "nalpha0 = ",nalpha0_
                     WRITE(iunit,"(A,2F12.7)") "!major, minor radius[m]= ",Rmajor, a
                     WRITE(iunit,"(A,F12.7)") "my_dpdx = ",dpdx
                     WRITE(iunit,"(A,F12.7)") "q0 = ",ABS(q)
                     IF (ABS(shat)<0.15) THEN
                        WRITE(iunit,"(A,F12.7)") "!shat = ",shat 
                        WRITE(iunit,"(A,F12.7)") "shat = ",0.0 
                     ELSE
                        WRITE(iunit,"(A,F12.7)") "shat = ",shat 
                     END IF
                     WRITE(iunit,"(A,I5)") "gridpoints = ",maxPnt
                     !WRITE(iunit,"(A,I5)") "n_pol = ",pol_turns_gene
                     WRITE(iunit,"(A,I5)") "n_pol = ",1
                     WRITE(iunit,"(A)") "/"
                  ELSE IF (.not. lglobal_txport) THEN
                     WRITE(iunit,'(A)') '&PARAMETERS'
                     WRITE(iunit,"(A,2F12.7)") "!s0, alpha0 = ",s,phi0
                     WRITE(iunit,"(A,2F12.7)") "!major, minor radius[m]= ",Rmajor, a
                     WRITE(iunit,"(A,F12.7)") "my_dpdx = ",dpdx
                     WRITE(iunit,"(A,F12.7)") "q0 = ",ABS(q)
                     IF (ABS(shat)<0.15) THEN
                        WRITE(iunit,"(A,F12.7)") "!shat = ",shat 
                        WRITE(iunit,"(A,F12.7)") "shat = ",0.0 
                     ELSE
                        WRITE(iunit,"(A,F12.7)") "shat = ",shat 
                     END IF
                     WRITE(iunit,"(A,I5)") "gridpoints = ",maxPnt
                     !WRITE(iunit,"(A,I5)") "n_pol = ",pol_turns_gene
                     WRITE(iunit,"(A,I5)") "n_pol = ",1
                     WRITE(iunit,"(A)") "/"
                  END IF
                  CALL FLUSH(iunit)
                  sflCrd0(3) = phi0
                  ! Follow Line
                  DO i = 1, maxPnt
                     th = -pi + (i-1)*dtheta
                     sflCrd(1) = sflCrd0(1)
                     sflCrd(2) = th
                     sflCrd(3) = sflCrd0(3) + q*(th-sflCrd0(2))
                     ! Get Metric Elements
                     !IF (pest) THEN
                        sflCrd_sav = sflCrd
                        CALL pest2vmec(sflCrd) ! Returns on 2pi grid
                     !ELSE
                     !END IF
                     ! Calc the Jacobian
                     !  Vector3d F (ct().Jmatr[0]); // (R,fi,Z)
                     !  Vector3d Fs(ct().Jmatr[1]);// d(R,fi,Z)/ds   -- partial derivative on s
                     !  Vector3d Ft(ct().Jmatr[2]);// d(R,fi,Z)/dtheta
                     !  Vector3d Fp(ct().Jmatr[3]);// d(R,fi,Z)/dphi
                     !   double maxFlux = mFlux*mBnorm;
                     !  ct().jac = ct().jacVmec/maxFlux;      // Jacobian in coordinates (Flux,theta,phi)
                     !   // Multiplying by R we take into account the cylindrical coordinate system
                     !  Fs[1] *=R; // R*dfi/ds
                     !  Ft[1] *=R; // R*dfi/dtheta
                     !  Fp[1] *=R; // R*dfi/dphi
                     !  Vector3d Ff(Fs); // create a copy of Fs
                     !  Ff /= maxFlux;
                     !  ct().gradtorFlux = (Ft^Fp)/ct().jac; // grad(Psi_tor) in cyl. coord, Psi_tor is the toroidal flux
                     !  ct().gradS = ct().gradtorFlux;
                     !  ct().gradS/= maxFlux; // grad(s) in cyl. coord
                     !  ct().gradTheta = (Fp^Ff)/ct().jac; // grad(theta) in cyl. coord; theta is the poloidal angle
                     !  ct().gradPhi = (Ff^Ft)/ct().jac; // grad(phi)   in cyl. coord; phi is the toroidal angle
                     ! Now let's try some Brute force
                     u = sflCrd(2)
                     v = sflCrd(3)
                     IF (u < 0) THEN
                        u = -MOD(ABS(u),pi2)
                        u = u + pi2
                     END IF
                     IF (v < 0) THEN
                        v = -MOD(ABS(v),pi2)
                        v = v + pi2
                     END IF
                     IF (u > pi2) u = MOD(u,pi2)
                     IF (v > pi2) v = MOD(v,pi2)
                     CALL get_equil_RZ(sflCrd(1),u,v,R,Z,iflag,R_GRAD=R_grad,Z_GRAD=Z_grad)
                     ! Get equil_RZ returns dR/drho and dZ/drho
                     !   dR/ds=(0.5/rho)*dR/drho=(0.5/sqrt(s))*dR/drho
                     R_grad(3) = 0.5*R_grad(3)/SQRT(sflCrd(1))
                     Z_grad(3) = 0.5*Z_grad(3)/SQRT(sflCrd(1))
                     ! e_s
                     esubs(1) = R_grad(3)
                     esubs(2) = zero
                     esubs(3) = Z_grad(3)
                     ! e_u
                     esubu(1) = R_grad(1)
                     esubu(2) = zero
                     esubu(3) = Z_grad(1)
                     ! e_v
                     esubv(1) = R_grad(2)
                     esubv(2) = one
                     esubv(3) = Z_grad(2) 
                     esubv(1) = esubv(1)*nfp
                     esubv(3) = esubv(3)*nfp
                     ! Cylindrical Coordianates
                     !CALL EZspline_interp(R_spl,u,v,sflCrd(1),R,iflag)
                     esubs(2) = esubs(2)*R
                     esubu(2) = esubu(2)*R
                     esubv(2) = esubv(2)*R
                     ! sqrt(g) = R*(Ru*Zs-Rs*Zu)
                     sqrtg = R*(R_grad(1)*Z_grad(3)-R_grad(3)*Z_grad(1))
                     ! e^i = (e_j x e_k)/sqrt(g)
                     !CALL EZspline_interp(G_spl,u,v,sflCrd(1),sqrtg,iflag)
                     CALL cross_product(esubu,esubv,es)
                     CALL cross_product(esubv,esubs,eu)
                     CALL cross_product(esubs,esubu,ev)
                     es = es/sqrtg
                     eu = eu/sqrtg
                     ev = ev/sqrtg
                     ! Get Field (before we adjust eu so gradB is correct)
                     CALL get_equil_Bflx(sflCrd(1),u,v,temp1,temp2,temp3,iflag,modb,gradB)
                     gradB = gradB(3)*es + gradB(1)*eu + gradB(2)*ev*nfp
                     ! Now Adjust e^u for lambda
                     !IF (pest) THEN
                        CALL get_equil_L(sflCrd(1),u,v,temp1,iflag,gradlam)
                        eu = eu + gradlam(3)*es + gradlam(1)*eu + gradlam(2)*ev*nfp
                     !ELSE
                     !END IF
                     ! Now do some calculations
                     sflCrd = sflCrd_sav
                     thetastar = sflCrd(2)
                     gradS = es
                     gradThetaStar = eu
                     gradPhi = ev
                     gradAlpha = qprim * thetaStar*gradS+q*gradThetaStar-gradPhi
                     alpha     = q*thetaStar - sflCrd(3)
                     ! Metrice and Jacobian
                     gss = DOT_PRODUCT(gradS,gradS)
                     gsa = DOT_PRODUCT(gradS,gradAlpha)
                     gst = DOT_PRODUCT(gradS,gradThetaStar)
                     gaa = DOT_PRODUCT(gradAlpha,gradAlpha)
                     gat = DOT_PRODUCT(gradAlpha,gradThetaStar)
                     CALL cross_product(gradS,gradAlpha,wrk)
                     jac1 = one/DOT_PRODUCT(wrk,gradThetaStar)
                     
                     Bhat(ialpha,i) = modb/Ba
                     
                     g11(ialpha,i) = gss*a**2/(4*s)
                     !g12(ialpha,i) = gsa*a**2*iot/2*sloc_fac
                     g12(ialpha,i) = gsa*a**2*iot/2
                     g22(ialpha,i) = (Bhat(ialpha,i)**2+g12(ialpha,i)**2)/g11(ialpha,i)
                     abs_jac(ialpha,i) = ABS(jac1*2*q/a**3)
                     
                     
                     CALL cross_product(gradAlpha,gradThetaStar, es) 
                     CALL cross_product(gradThetaStar,gradS,     ea) 
                     CALL cross_product(gradS,gradAlpha,         et)  

                     es(:) = es(:)*jac1
                     ea(:) = ea(:)*jac1
                     et(:) = et(:)*jac1

                     gradB = gradB/Ba
                     dBds = DOT_PRODUCT(gradB,es)
                     dBda = DOT_PRODUCT(gradB,ea)
                     dBdt(ialpha,i) = DOT_PRODUCT(gradB,et)
                     
                     c = iot*iot*a**4
                     L1(ialpha,i) = q/sqrt(s)*(dBda + c*(gss*gat-gsa*gst)*dBdt(ialpha,i)/(4*Bhat(ialpha,i)**2))
                     L2(ialpha,i) = two*sqrt(s)*(dBds + c*(gaa*gst-gsa*gat)*dBdt(ialpha,i)/(4*Bhat(ialpha,i)**2))
                     !IF (nodrift) THEN
                     !   L2(ialpha,i) = zero
                     !   dpdx = zero
                     !END IF
                     kp1(ialpha,i) = L2(ialpha,i) - dpdx/two/Bhat(ialpha,i)
                     WRITE(iunit,"(8ES20.10)") g11(ialpha,i),g12(ialpha,i),g22(ialpha,i),&
                                               Bhat(ialpha,i),abs_jac(ialpha,i),L2(ialpha,i),&
                                               L1(ialpha,i),dBdt(ialpha,i); CALL FLUSH(iunit)
                     WRITE(iunit_curv,"(2ES20.10)") kp1(ialpha,i), L1(ialpha,i); CALL FLUSH(iunit_curv)
                  END DO
                  CALL FLUSH(iunit)
                  CALL FLUSH(iunit_curv)
                  IF (.not. lglobal_txport) THEN
                     CLOSE(iunit)
                     CLOSE(iunit_curv)
                  END IF
               END DO
               IF (lglobal_txport) THEN
                  CLOSE(iunit)
                  CLOSE(iunit_curv)
               END IF
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !!         TXPORT Section
               !!             Here we calculate the transport
               !!             Note GENE modifies L2 by subtracting 0.5*dpdx, if you wish to use
               !!             this curvature as GENE does you may simply reference the kp1
               !!             variable.
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               temp_str = 'txport_out.'
               iunit = 25
               WRITE(num_str,'(A,I3.3)') '_',ik
               CALL safe_open(iunit,iflag,TRIM(TRIM(temp_str)//TRIM(proc_string)//TRIM(num_str)),'unknown','formatted')
               rll  = g12/g11 + qpthk
               rllm = Bhat*Bhat/g11+(one+g11*rll)*(one+g11*rll)
               vrbl = g12/g11
               dtheta = pi2/(maxPnt-1)
               sloc(:,2:maxPnt) = (vrbl(:,2:maxPnt)-vrbl(:,1:maxPnt-1))/dtheta
               sloc(:,1) = sloc(:,2)
               ! Disabled for test 08192016
               DO ialpha = 1, nalpha0_
                   !rkp1av(ialpha,:) = kp1(ialpha,:)
                   !slocav(ialpha,:) = sloc(ialpha,:)
                  temp_arr=kp1(ialpha,:); tout_arr(:) = 0
                  CALL smoothg(temp_arr,maxPnt,dl_z,tout_arr)
                  rkp1av(ialpha,:) = tout_arr(:)
                  temp_arr=sloc(ialpha,:); tout_arr(:) = 0
                  CALL smoothg(temp_arr,maxPnt,dl_z,tout_arr)
                  slocav(ialpha,:) = tout_arr(:)
               END DO
               rkpn = kp1
               rkpg = L1
               rkk  = 2*dpdx*abs_jac*(shat+iot*ABS(q))*(rkpn+rkpg*rll)
               ! Calculate numerical derivatives
               rllmp(:,2:maxPnt-1)  = (rllm(:,3:maxPnt)-rllm(:,1:maxPnt-2))/(2*dtheta)
               rllmp(:,1)           = two*rllmp(:,2)-rllmp(:,3)
               rllmpp(:,2:maxPnt-1) = (rllm(:,3:maxPnt)-two*rllm(:,2:maxPnt-1)+rllm(:,1:maxPnt-2))/(dtheta*dtheta)
               rllmpp(:,1)          = two*rllmpp(:,2)-rllmpp(:,3)
               rllmp                = iot*rllmp
               rllmpp               = iot*iot*rllmpp
               vv                   = rkk / rllm
               ! Calculated q_proxy     
               CALL tolower(txport_proxy)
               SELECT CASE(TRIM(txport_proxy))
                  CASE('prox1')
                     qqfac=(one+one/(rhsrkp*(one+(tau_s*slocav)**2)))
                     dkp  = rkp_p-rkp_cr
                     dkpfac = zero
                     WHERE(dkp > zero) dkpfac = dkp
                     vkp1fac = zero
                     WHERE(rkp1av < zero) vkp1fac = - rkp1av
                     vqqprox = (dk/rkp_p)*sqrt(tau*vkp1fac*dkpfac)*qqfac
                  CASE('prox1_sam')
                     qqfac=(one+one/(rhsrkp*(one+(tau_s*sloc)**2)))
                     dkp  = rkp_p-rkp_cr
                     dkpfac = zero
                     WHERE(dkp > zero) dkpfac = dkp
                     vkp1fac = kp1
                     WHERE(kp1 < zero) vkp1fac = - kp1
                     vqqprox = (dk/rkp_p)*sqrt(tau*vkp1fac*dkpfac)*qqfac
                     WHERE(kp1 > zero) vqqprox = -0.01*vqqprox
                     !WHERE(rkp1av < zero) vqqprox = vqqprox
                  CASE('prox1b')
                     qqfac=(one+one/(rhsrkp*(one+(zero*slocav)**2)))
                     dkp  = rkp_p-rkp_cr
                     dkpfac = zero
                     WHERE(dkp > zero) dkpfac = dkp
                     vkp1fac = zero
                     WHERE(rkp1av < zero) vkp1fac = - rkp1av
                     vqqprox = (dk/rkp_p)*sqrt(tau*vkp1fac*dkpfac)*qqfac
                  CASE('prox1d')
                     qqfac=g11
                     dkp  = rkp_p-rkp_cr
                     dkpfac = zero
                     WHERE(dkp > zero) dkpfac = dkp
                     vkp1fac = zero
                     WHERE(rkp1av < zero) vkp1fac = - rkp1av
                     vqqprox = (dk/rkp_p)*sqrt(tau*vkp1fac*dkpfac)*qqfac
                  CASE('prox_l2')
                     vqqprox = L2
                     WHERE(vqqprox > zero) vqqprox=vqqprox*0.01
                  CASE('prox_tem_proll','temproxy','tem_overlap')
                     DO ialpha = 1, nalpha0_
                        bmax = MAXVAL(Bhat(ialpha,:))
                        bmin = MINVAL(Bhat(ialpha,:))
                        b0   = 0.5_rprec*(bmin+bmax)
                        bamp = 0.5_rprec*(bmax-bmin)
                        vqqprox(ialpha,:) = kp1(ialpha,:)*(Bhat(ialpha,:)-b0)/bamp  ! no minus sign necessary
                     END DO
                  CASE('temproxyadvanced','tem_bounce')
                     DO ialpha = 1, nalpha0_
                        bmax = MAXVAL(Bhat(ialpha,:))
                        bmin = MINVAL(Bhat(ialpha,:))
                        ier = 0
                        IF (EZspline_allocated(Bhat_spl)) CALL EZspline_free(Bhat_spl,ier)
                        IF (EZspline_allocated(L2_spl)) CALL EZspline_free(L2_spl,ier)
                        CALL EZspline_init(Bhat_spl,maxPnt,bcs1,ier)
                        CALL EZspline_init(L2_spl,maxPnt,bcs1,ier)
                        Bhat_spl%isHermite = 1
                        L2_spl%isHermite = 1
                        CALL EZspline_setup(Bhat_spl,Bhat(ialpha,:),ier)
                        CALL EZspline_setup(L2_spl,kp1(ialpha,:),ier)
                        A2 = one/bmax; B2 = one/bmin
                        CALL DQAG(TEMfunc_proll,A2,B2,EPSABS,EPSREL,KEY,result,abserr,neval,ier,&
                                    128,NW,last,iwork,work)
                        vqqprox(ialpha,:)=-result ! minus sign necessary, assumes omega* negative
                        !PRINT *,'temproxyadvanced',result,ier
                        CALL EZspline_free(Bhat_spl,ier)
                        CALL EZspline_free(L2_spl,ier)
                     END DO
                  CASE('tem_bounce_tau')
                     DO ialpha = 1, nalpha0_
                        bmax = MAXVAL(Bhat(ialpha,:))
                        bmin = MINVAL(Bhat(ialpha,:))
                        ier = 0
                        IF (EZspline_allocated(Bhat_spl)) CALL EZspline_free(Bhat_spl,ier)
                        IF (EZspline_allocated(L2_spl)) CALL EZspline_free(L2_spl,ier)
                        CALL EZspline_init(Bhat_spl,maxPnt,bcs1,ier)
                        CALL EZspline_init(L2_spl,maxPnt,bcs1,ier)
                        Bhat_spl%isHermite = 1
                        L2_spl%isHermite = 1
                        CALL EZspline_setup(Bhat_spl,Bhat(ialpha,:),ier)
                        CALL EZspline_setup(L2_spl,kp1(ialpha,:),ier)
                        A2 = one/bmax; B2 = one/bmin
                        CALL DQAG(TEMfunc_proll_omegatau,A2,B2,EPSABS,EPSREL,KEY,result,abserr,neval,ier,&
                                    128,NW,last,iwork,work)
                        vqqprox(ialpha,:)=-result ! minus sign necessary, assumes omega* negative
                        !PRINT *,'temproxyadvanced',result,ier
                        CALL EZspline_free(Bhat_spl,ier)
                        CALL EZspline_free(L2_spl,ier)
                     END DO
                  CASE('all')
                     ! prox1
                     qqfac=(one+one/(rhsrkp*(one+(tau_s*slocav)**2)))
                     dkp  = rkp_p-rkp_cr
                     dkpfac = zero
                     WHERE(dkp > zero) dkpfac = dkp
                     vkp1fac = zero
                     WHERE(rkp1av < zero) vkp1fac = - rkp1av
                     vqqprox = (dk/rkp_p)*sqrt(tau*vkp1fac*dkpfac)*qqfac
                     txport_q_all(1,ik,:,:) = vqqprox(:,:)
                     ! prox1b
                     qqfac=(one+one/(rhsrkp*(one+(zero*slocav)**2)))
                     dkp  = rkp_p-rkp_cr
                     dkpfac = zero
                     WHERE(dkp > zero) dkpfac = dkp
                     vkp1fac = zero
                     WHERE(rkp1av < zero) vkp1fac = - rkp1av
                     vqqprox = (dk/rkp_p)*sqrt(tau*vkp1fac*dkpfac)*qqfac
                     txport_q_all(2,ik,:,:) = vqqprox(:,:)
                     !prox1d
                     qqfac=g11
                     dkp  = rkp_p-rkp_cr
                     dkpfac = zero
                     WHERE(dkp > zero) dkpfac = dkp
                     vkp1fac = zero
                     WHERE(rkp1av < zero) vkp1fac = - rkp1av
                     vqqprox = (dk/rkp_p)*sqrt(tau*vkp1fac*dkpfac)*qqfac
                     txport_q_all(3,ik,:,:) = vqqprox(:,:)
                     !'tem_overlap'
                     DO ialpha = 1, nalpha0_
                        bmax = MAXVAL(Bhat(ialpha,:))
                        bmin = MINVAL(Bhat(ialpha,:))
                        b0   = 0.5_rprec*(bmin+bmax)
                        bamp = 0.5_rprec*(bmax-bmin)
                        vqqprox(ialpha,:) = kp1(ialpha,:)*(Bhat(ialpha,:)-b0)/bamp  ! no minus sign necessary
                     END DO
                     txport_q_all(4,ik,:,:) = vqqprox(:,:)
                     !'temproxyadvanced'
                     DO ialpha = 1, nalpha0_
                        bmax = MAXVAL(Bhat(ialpha,:))
                        bmin = MINVAL(Bhat(ialpha,:))
                        ier = 0
                        IF (EZspline_allocated(Bhat_spl)) CALL EZspline_free(Bhat_spl,ier)
                        IF (EZspline_allocated(L2_spl)) CALL EZspline_free(L2_spl,ier)
                        CALL EZspline_init(Bhat_spl,maxPnt,bcs1,ier)
                        CALL EZspline_init(L2_spl,maxPnt,bcs1,ier)
                        Bhat_spl%isHermite = 1
                        L2_spl%isHermite = 1
                        CALL EZspline_setup(Bhat_spl,Bhat(ialpha,:),ier)
                        CALL EZspline_setup(L2_spl,kp1(ialpha,:),ier)
                        A2 = one/bmax; B2 = one/bmin
                        CALL DQAG(TEMfunc_proll,A2,B2,EPSABS,EPSREL,KEY,result,abserr,neval,ier,&
                                    128,NW,last,iwork,work)
                        vqqprox(ialpha,:)=-result ! minus sign necessary, assumes omega* negative
                        !PRINT *,'temproxyadvanced',result,ier
                        CALL EZspline_free(Bhat_spl,ier)
                        CALL EZspline_free(L2_spl,ier)
                     END DO
                     txport_q_all(5,ik,:,:) = vqqprox(:,:)
                     !'tem_bounce_tau'
                     DO ialpha = 1, nalpha0_
                        bmax = MAXVAL(Bhat(ialpha,:))
                        bmin = MINVAL(Bhat(ialpha,:))
                        ier = 0
                        IF (EZspline_allocated(Bhat_spl)) CALL EZspline_free(Bhat_spl,ier)
                        IF (EZspline_allocated(L2_spl)) CALL EZspline_free(L2_spl,ier)
                        CALL EZspline_init(Bhat_spl,maxPnt,bcs1,ier)
                        CALL EZspline_init(L2_spl,maxPnt,bcs1,ier)
                        Bhat_spl%isHermite = 1
                        L2_spl%isHermite = 1
                        CALL EZspline_setup(Bhat_spl,Bhat(ialpha,:),ier)
                        CALL EZspline_setup(L2_spl,kp1(ialpha,:),ier)
                        A2 = one/bmax; B2 = one/bmin
                        CALL DQAG(TEMfunc_proll_omegatau,A2,B2,EPSABS,EPSREL,KEY,result,abserr,neval,ier,&
                                    128,NW,last,iwork,work)
                        vqqprox(ialpha,:)=-result ! minus sign necessary, assumes omega* negative
                        !PRINT *,'temproxyadvanced',result,ier
                        CALL EZspline_free(Bhat_spl,ier)
                        CALL EZspline_free(L2_spl,ier)
                     END DO
                     txport_q_all(6,ik,:,:) = vqqprox(:,:)
                  CASE('gene','gene_serial')
!DEC$ IF DEFINED (GENE)
                     DO ialpha = 1, nalpha0_
   	                  ! Set some MPI Defaults
                        WRITE(num_str,'(2(A,I3.3))') '_',ik,'_',ialpha
           	            omp_level = MPI_THREAD_SINGLE
                        n_procs_sim = 1
                        ! First we write the parameters file
                        diagdir = '.'
                        geomdir = '.'
                        magn_geometry = 'gist'
                        geomfile = gist_filename
                        gene_dir_str = 'skip_parfile'
                        file_extension = '_'//TRIM(proc_string)//TRIM(num_str)
                        gene_str = '_'//TRIM(proc_string)//TRIM(num_str)
                        ch_in = 'no'
                        mype_gl = myid

                        ! RUN GENE
                        IF (lscreen .and. (ik==1)) WRITE(6,'(a)') ' ::::::::::::::::::::        LINEAR GENE RUN     ::::::::::::::::::::::'
                        IF (.not.(lscreen .and. (ik==1))) OPEN(UNIT=6,FILE='log_gene.'//TRIM(gene_str(2:256)))
                        CALL rungene(MPI_COMM_SELF,gene_dir_str,gene_str,ch_in)
                        CALL erase_stop_file
                        IF (.not.(lscreen .and. (ik==1))) CLOSE(UNIT=6)

                        ! Open EIGENVALUE FILE
                        gene_eig_r = zero
                        gene_eig_i = one
                        temp_str = 'eigenvalues_'//TRIM(proc_string)//TRIM(num_str)
                        INQUIRE(FILE=TRIM(temp_str),EXIST=file_test_eig)
                        IF (file_test_eig) THEN
                           CALL safe_open(iunit_eigen,iflag,TRIM(temp_str),'old','formatted')
                           READ(iunit_eigen,*) 
                           READ(iunit_eigen,*) gene_eig_r,gene_eig_i
                           CLOSE(iunit_eigen)
                        ENDIF
                        
                        ! Open Omega File
                        temp_str = 'omega_'//TRIM(proc_string)//TRIM(num_str)
                        INQUIRE(FILE=TRIM(temp_str),EXIST=file_test_omega)
                        IF (file_test_omega) THEN
                           CALL safe_open(iunit_eigen,iflag,TRIM(temp_str),'old','formatted')
                           READ(iunit_eigen,*) gene_kymin,gene_eig_r,gene_eig_i
                           CLOSE(iunit_eigen)
                        ENDIF

                        ! Set Turbulent Flux value
                        vqqprox(ialpha,:) = gene_eig_r
                        IF (lscreen .and. (ik==1)) WRITE(6,'(a)') ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
                     END DO
                  CASE('gene_parallel')
                     DO ialpha = 1, nalpha0_
                        kx_gene = shat*kymin*phi0_gl(ialpha)
                        ! Write a new parameters file with the proper gist filename
                        WRITE(num_str,'(2(A,I3.3))') '_',ik,'_',ialpha
                        geomfile = gist_filename
                        gene_str = '_'//TRIM(proc_string)//TRIM(num_str)

                        ! Join up with the workers
                        lscreen_gene = .false.
                        IF (lscreen .and. (ik == 1) .and. (ialpha == 1)) lscreen_gene = .true.
                        CALL stellopt_paraexe(txport_proxy,gene_str,lscreen_gene)

                        ! Open EIGENVALUE FILE
                        gene_eig_r = zero
                        gene_eig_i = one
                        temp_str = 'eigenvalues_'//TRIM(proc_string)//TRIM(num_str)
                        INQUIRE(FILE=TRIM(temp_str),EXIST=file_test_eig)
                        IF (file_test_eig) THEN
                           ier  = 0
                           CALL safe_open(iunit_eigen,iflag,TRIM(temp_str),'old','formatted')
                           READ(iunit_eigen,*) ! HEADER
                           DO
                              READ(iunit_eigen,*,IOSTAT=ier) temp1, temp2 !gene_eig_r,gene_eig_i
                              IF (ier < 0) EXIT
                              IF (temp1 > 0) gene_eig_r = gene_eig_r+temp1
                              !gene_eig_i = gene_eig_i+temp2
                           END DO
                           !PRINT *,gene_eig_r
                           CLOSE(iunit_eigen)
                        ENDIF
                        
                        ! Open Omega File
                        temp_str = 'omega_'//TRIM(proc_string)//TRIM(num_str)
                        INQUIRE(FILE=TRIM(temp_str),EXIST=file_test_omega)
                        IF (file_test_omega) THEN
                           CALL safe_open(iunit_eigen,iflag,TRIM(temp_str),'old','formatted')
                           READ(iunit_eigen,*) gene_kymin,gene_eig_r,gene_eig_i
                           CLOSE(iunit_eigen)
                        ENDIF

                        ! Set Turbulent Flux value
                        vqqprox(ialpha,:) = gene_eig_r
                     END DO
!DEC$ ENDIF
               END SELECT
               WHERE (vqqprox /= vqqprox) vqqprox = bigno
               txport_q(ik,:,:) = vqqprox(:,:)
               IF (lscreen .and. (ik==1)) WRITE(6,'(A)') '     s    alpha0   iota    shear    Bref    Lref  Q_turb'
               IF (lscreen) WRITE(6,'(7(F8.4))')  s,alpha_start_txport,ABS(iot),shat,Ba,a,SUM(SUM(vqqprox,DIM=2)/maxPnt,DIM=1)/nalpha0_
               WRITE(iunit,'(1pe14.6)') SUM(SUM(vqqprox,DIM=2)/maxPnt,DIM=1)/nalpha0_
               DO i = 1, maxPnt
                  WRITE(iunit,'(ES22.14)') vqqprox(:,i)
               END DO
               CLOSE(iunit)
               CALL FLUSH(6)
            END DO
            DEALLOCATE(phi0_gl,g11,g12,g22,Bhat,abs_jac,L2,L1,dBdt,kp1)
            DEALLOCATE(temp_arr,tout_arr)
            DEALLOCATE(rll,rllm,vrbl,sloc,rkp1av,slocav,rkpn,rkpg,rkk,rllmp,rllmpp,vv,dkp,dkpfac,vkp1fac,qqfac,vqqprox)
         CASE('spec')
      END SELECT
      IF (lscreen) WRITE(6,'(a)') ' -------------------------  TRANSPORT CALC DONE  ----------------------'
!DEC$ ENDIF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_txport
