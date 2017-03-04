cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Miscellaneous subroutines 
c   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c----------------------------------------------------------------------------
c   This subroutine performs a polynomial interpolation of the data table
c       defined by xa and ya of length(N) at the point x.  Output is y with
c       estimated error dy.  Based on modified Neville's algorithm from 
c       Numerical Recipies W. Press et.al.
c     JL 7/22/09
c----------------------------------------------------------------------------
c
      subroutine pol_int(xa,ya,N,x,y,dy)
      use penta_kind_mod
      implicit none
      integer(iknd) :: N
      real(rknd), dimension(N) :: xa, ya
      real(rknd) :: x,y,dy
      !local variables
      real(rknd), dimension(N) :: C,D, diff_vec
      integer(iknd) :: i,ns, m
      real(rknd) :: diff, difftmp, h0, hm, dh, W

      !find the nearest entry
      ns=1
      diff=dabs(x-xa(1))
      do i=1,N
        difftmp=dabs(x-xa(i))
        if ( difftmp .lt. diff ) then
          ns=i
          diff=difftmp
        endif
      enddo

      !initialize C,D tableau
      C=ya
      D=ya

      !initial approximation to y
      y=ya(ns)

      !move through column of tableau and calculate Cmi,Dmi
      ns=ns-1
      do m=1,N-1        !column
        do i=1,N-m      !row
          !components of C,D
          h0=xa(i)-x
          hm=xa(i+m)-x
          dh=h0-hm
          !check for equal xa values, which would give DIV0
          if (dh .eq. 0) then
            write(*,*) 'error in interp, 2 xa are equal'
            stop
          endif
          !define C,D for this column
          W=C(i+1)-D(i)
          C(i)=h0*W/dh
          D(i)=hm*W/dh
        enddo

        !choose which branch to take
        if (2_iknd*ns .lt. N-m) then
          dy=C(ns+1)
        else
          dy=D(ns)
          ns=ns-1
        endif
        y=y+dy
      enddo

      end subroutine pol_int
c
c----------------------------------------------------------------------------
c   This subroutine defines the classical friction coefficients and 
c    returns them as the matrix L(species1,species2) where
c    species1,2=[electrons, ion1, ion2, ...] and each L(s1,s2) is a 4x4
c    matrix of friction coefficients: [ls1s2_11, -ls1s2_12; -ls1s2_21, -ls1s2_22].
c     The total matrix lmat is thus num_species*2 square.
c   JL 7/2009
c----------------------------------------------------------------------------
c
      subroutine define_friction_coeffs(masses,charges,v_ths,Temps,dens,
     1 loglambda,num_species,lmat)
      use penta_kind_mod
      use phys_const
      implicit none
      !dummy variables
      integer(iknd) :: num_species
      real(rknd), dimension(num_species) :: masses, charges, v_ths, 
     1  Temps, dens
      real(rknd) :: loglambda
      real(rknd), dimension(num_species*2,num_species*2) :: lmat
      !local variables
      integer(iknd) :: ind1, ind2, ispec1, ispec2
      real(rknd), dimension(num_species) :: Mak_p1,
     1  Mak_p2, tau_ak,  Kak, Mak00, Mak10, Mak11, Mak01,
     2  Nak11 
      real(rknd) :: tau_coeff, ma, na, charge_a, Ta, vth_a, lab11, 
     1  lab12, lab21, lab22, Kba, Mba10, Nab01, Nab00, Nab10, Nab11,
     2  tau_ab, nama
      real(rknd), parameter :: one=1._rknd

      !coefficient on collision times
      tau_coeff=3._rknd*dsqrt(pi)*pi*eps0**2_iknd

      !This nested loop defines the upper half (including the diagonal) of
      ! the large lmat matrix

      !loop over species a (row)
      do ispec1 = 1,num_species
        !assign parameters for species a
        ma=masses(ispec1)
        Ta=Temps(ispec1)
        charge_a=charges(ispec1)
        na=dens(ispec1)
        vth_a=v_ths(ispec1)
        Kak=(v_ths/vth_a)**2_iknd !Normalized thermal velocities

        !calculate moments of collision operator Mij,Nij
        Mak_p1=-one*(one + ma/masses)
        Mak_p2=one + Kak
        Mak00=Mak_p1*Mak_p2**(-1.5_rknd)
        Mak10=1.5_rknd*Mak_p1*Mak_p2**(-2.5_rknd)
        Mak11=-1._rknd*(13._rknd/4._rknd+4._rknd*Kak+
     1    7.5_rknd*Kak**2_iknd)*Mak_p2**(-2.5_rknd)
        Mak01=Mak10
        Nak11=6.75_rknd*(Ta/Temps)*Kak*Mak_p2**(-2.5_rknd)

        !calculate collision times
        tau_ak=tau_coeff*ma**2_iknd*vth_a**3_iknd/(dens*charge_a
     1          **2_iknd*charges**2_iknd*loglambda)

        !loop over second species (columns)
        do ispec2 = ispec1,ispec1+(num_species-ispec1)

          !assign values for this species calculated above
          Nab00=-Mak00(ispec2)
           Nab10=-Mak10(ispec2)
          Nab11= Nak11(ispec2)
          tau_ab=tau_ak(ispec2)
  
          !calculate remaining values
          Kba=(vth_a/v_ths(ispec2))**2_iknd
          Mba10=-1.5_rknd*(one+masses(ispec2)/ma)*(one+Kba)**(-2.5_rknd)
          Nab01=-Mba10*(Ta*vth_a)/(Temps(ispec2)*v_ths(ispec2))

          nama=na*ma

          !define friction matrix components (4x4 matrices)
          if (ispec1 == ispec2) then  !diagonal
            lab11=nama*( sum(Mak00/tau_ak) + Nab00/tau_ab )
            lab21=nama*( sum(Mak10/tau_ak) + Nab10/tau_ab )
            lab22=nama*( sum(Mak11/tau_ak) + Nab11/tau_ab )
            lab12=nama*( sum(Mak01/tau_ak) + Nab01/tau_ab )
          else                      !off-diagonal
            lab11=nama*( Nab00/tau_ab )
            lab21=nama*( Nab10/tau_ab )
            lab22=nama*( Nab11/tau_ab )
            lab12=nama*( Nab01/tau_ab )
          endif

        !define indices for full matrix and assign components
        ind1=(ispec1 - 1)*2+1
        ind2=(ispec2 - 1)*2+1
 
        lmat(ind1,ind2)=lab11
        lmat(ind1,ind2+1)=-lab12
        lmat(ind1+1,ind2)=-lab21
        lmat(ind1+1,ind2+1)=lab22

        enddo !species b loop
      enddo  !species a loop
  
      !define other half of lmat (lower half)
      do ispec1=2,num_species
        do ispec2=1,(ispec1 - 1)

          !define indices for full matrix and assign components from transpose
          ! of upper half
          ind1=(ispec1 - 1)*2+1
          ind2=(ispec2 - 1)*2+1
    
          lmat(ind1:ind1+1,ind2:ind2+1)
     1        =transpose(lmat(ind2:ind2+1,ind1:ind1+1))
        enddo
      enddo
      end subroutine define_friction_coeffs
c
c----------------------------------------------------------------------------
c   This subroutine calculates the parallel, poloidal and toroidal (in 
c     Boozer coordinates) flows.  These are returned in vectors arranged
c     as [electron, I1, I2, ....]
c   JL 7/2009
c----------------------------------------------------------------------------
c
      subroutine calculate_flows(num_species,L1,L2,L3,N1,N2,N3,Xvec,
     1  gamma_e,gamma_i,q_e,q_i,charges,B_uprl,B_qprl,u_theta,u_zeta)
      use penta_kind_mod
      use phys_const
      use pprof_pass
      use vmec_var_pass
      implicit none
      !dummy variables
      integer(iknd) :: num_species
      real(rknd), dimension(2*num_species+1) :: Xvec
      real(rknd), dimension(num_species) :: L1,L2,L3,N1,N2,N3,
     1  B_uprl,B_qprl, u_theta, u_zeta, charges
      real(rknd) :: gamma_e, q_e
      real(rknd), dimension(num_species-1) :: gamma_i, q_i
      !local variables
      real(rknd), dimension(2,2) :: Ntmp, Ntmp_inv, Ltmp
      real(rknd), dimension(2) :: tmp1, tmp2, UandQ
      integer(iknd) :: ispec, tmp_ind
      real(rknd) :: ucoeff

      !calculate species' flows
      do ispec=1,num_species
        !define 4x4 thermal coefficient matrices and inverse of Na
        Ntmp(1,1:2)=[N1(ispec), N2(ispec)]
        Ntmp(2,1:2)=[N2(ispec), N3(ispec)]
        call mat_2by2_inverse(Ntmp,Ntmp_inv)
        Ltmp(1,1:2)=[L1(ispec), L2(ispec)]
        Ltmp(2,1:2)=[L2(ispec), L3(ispec)]

        !define parallel flows for ions and electrons
        tmp_ind=1+2*(ispec-1)
        tmp1=matmul(-Ltmp,Xvec(tmp_ind:tmp_ind+1))
        if (ispec .eq. 1) then    !electrons
          tmp2(1)=tmp1(1)+gamma_e
          tmp2(2)=tmp1(2)+q_e/elem_charge/Te
          UandQ=Bsq*matmul(Ntmp_inv,tmp2)
          B_uprl(ispec)=UandQ(1)
          B_qprl(ispec)=UandQ(2)*5._rknd/2._rknd*elem_charge*Te*ne
        else                      !ions
          tmp2(1)=tmp1(1)+gamma_i(ispec-1)
          tmp2(2)=tmp1(2)+q_i(ispec-1)/elem_charge/Ti(ispec-1)
          UandQ=Bsq*matmul(Ntmp_inv,tmp2)
          B_uprl(ispec)=UandQ(1)
          B_qprl(ispec)=UandQ(2)*5._rknd/2._rknd*elem_charge*
     1      Ti(ispec-1)*ni(ispec-1)
        endif

        !define poloidal and toroidal flows
        ucoeff=4._rknd*pi**2_iknd/vol_p
        u_theta(ispec)=chip*ucoeff*(B_uprl(ispec)/Bsq - 
     1    Xvec(tmp_ind)*bzeta/(charges(ispec)*chip*Bsq))
        u_zeta(ispec)=psip*ucoeff*(B_uprl(ispec)/Bsq +
     1    Xvec(tmp_ind)*btheta/(charges(ispec)*psip*Bsq))
      enddo
      end subroutine calculate_flows
c
c----------------------------------------------------------------------------
c   This subroutine calculates the Pfirsh-Schlueter particle and heat
c       fluxes.  These are returned in vectors arranged as [electron, I1, I2, ....]
c   JL 7/2009
c----------------------------------------------------------------------------
c
      subroutine calculate_PS_flux(charges,num_species,lmat,Xvec,U2,
     1  gamma_PS,q_PS)
      use penta_kind_mod
      use phys_const
      implicit none
      !dummy variables
      integer(iknd) :: num_species
      real(rknd) :: U2
      real(rknd),dimension(num_species*2,num_species*2) :: lmat
      real(rknd) :: Xvec(2*num_species+1)
      real(rknd),dimension(num_species) :: gamma_PS, q_PS, charges
      !local variables
      real(rknd) :: lX_sum1,lX_sum2
      integer(iknd) :: ispec1, ispec2, spec1_ind, spec2_ind

      !loop over species and calculate product of friction coefficients
      !  and thermodynamic forces, then PS fluxes
      do ispec1=1,num_species
        spec1_ind=(ispec1-1)*2+1
        !sum product over species
        lX_sum1=0._rknd
        lX_sum2=0._rknd
        do ispec2=1,num_species
          spec2_ind=(ispec2-1)*2+1
          lX_sum1=lX_sum1+sum(lmat(spec1_ind,spec2_ind:spec2_ind+1)*
     1      Xvec(spec2_ind:spec2_ind+1))/charges(ispec2)
          lX_sum2=lX_sum2+sum(lmat(spec1_ind+1,spec2_ind:spec2_ind+1)*
     1      Xvec(spec2_ind:spec2_ind+1)/charges(ispec2))
        enddo
        !calculate fluxes
        gamma_PS(ispec1)=-U2*lX_sum1/charges(ispec1)
        q_PS(ispec1)=-U2*lX_sum2/charges(ispec1)
      enddo
      end subroutine calculate_PS_flux
c
c----------------------------------------------------------------------------
c   This subroutine determines the ambipolar roots as gamma_e=sum(Z*gammai)
c     where the sum is over ion species.  The zeros are found from linear
c     interpolation and using the functions zeroin and rad_flux.  
c   The roots (er_roots(:)) and number of roots (num_roots) are passed 
c       using the module Er_roots_pass.
c   If an even number of roots is found, the first root only is chosen.
c   JL 7/2009
c----------------------------------------------------------------------------
c
      subroutine find_Er_roots(gamma_e,gamma_i,Er_test_vals,Z_ion,
     1  num_Er_test,nis,qg_all_i)
      use penta_kind_mod
      use Er_fit_pass
      use Er_roots_pass
	use penta_functions_mod
      implicit none
      !dummy variables
      integer(iknd) :: num_Er_test, nis
      real(rknd), dimension(num_Er_test) :: gamma_e, Er_test_vals
      real(rknd),dimension(nis) :: Z_ion(nis)
      !local variables
      integer(iknd) :: ispec, ie, iroot, ind, root_diff
      integer(iknd), allocatable :: I_vr(:)
      real(rknd), dimension(num_Er_test) :: I_vr_tmp, qg_all_i
      real(rknd), dimension(nis,num_Er_test) :: gamma_i
      real(rknd) :: a_test, b_test, tol_Er
      !external functions
!      real(rknd),external :: zeroin,rad_flux
    
      !sum ion fluxes times charge numbers
      qg_all_i=0._rknd
      do ispec=1,nis
        qg_all_i=qg_all_i+Z_ion(ispec)*gamma_i(ispec,:)
      enddo

      !find where the sign switches
      num_roots=0_iknd
      do ie=1_iknd,num_Er_test-1_iknd
        a_test=dsign(1._rknd,gamma_e(ie)-qg_all_i(ie))
        b_test=dsign(1._rknd,gamma_e(ie+1)-qg_all_i(ie+1))
        if (a_test .ne. b_test) then
          num_roots=num_roots+1_iknd
          I_vr_tmp(num_roots)=ie
        endif
      enddo

      !allocate variables by number of roots
      allocate(I_vr(num_roots))
      allocate(er_roots(num_roots))

      !indices of the roots
      I_vr=I_vr_tmp(1:num_roots)

      !check for zero or even number of roots
      if (num_roots .eq. 0) then
        write(*,*) 'No roots found in search range'
        stop
      elseif ( mod(num_roots,2) .eq. 0) then
        write(*,*) 'Even number of roots found, choosing first only.'
        I_vr=I_vr(1)
      endif

      !find zeros using linear interpolation
      tol_er = 1.e-10_rknd    !tolerance on zeroin function
c     num_pts=2               !number of points used in fit (around root), order of polynomial fit+1

      !allocate fit variables
      allocate(diff_qg(num_pts))
      allocate(Er_fit(num_pts))

      !loop over roots and interpolate to find fine root
      do iroot=1,num_roots
        !this index
        ind=I_vr(iroot)

c          !make sure there are not multiple roots in num_pts range
c       if ((num_roots .gt. 1) .and. (iroot .ne. num_roots)) then 
c         root_diff=I_vr(iroot+1)-ind
c         if (root_diff .le. num_pts-1) then
c           write(*,*) 'Distance between roots too small for ',
c     1       'interpolation.  Increase num_Er_pts.'
c           write(*,*) I_vr
c           stop
c         endif
c       endif

        !define difference in flux and Er over fit range
        diff_qg=gamma_e(ind:ind+1)-qg_all_i(ind:ind+1)
        Er_fit=Er_test_vals(ind:ind+1)
    
        !fit the values over range and use zeroin and rad_flux to find zero
        er_roots(iroot) = zeroin(Er_fit(1),Er_fit(num_pts)
     1    ,rad_flux,tol_er)

      enddo
      !deallocate variables
      deallocate(I_vr,diff_qg,Er_fit)
      end subroutine find_Er_roots
c
c----------------------------------------------------------------------
c   This subroutine calculates the flux and force matrices.  These
c     matrices are those which multiply the vectors of fluxes and
c     forces in the equation A*fluxes=B*forces.  Note that these
c     fluxes are the total fluxes minus the PS fluxes and do not
c     include the bootstrap currents.  The matrix arrangements are
c       the same as lmat above.
c   JL 7/2009
c----------------------------------------------------------------------
c
      subroutine define_flux_force_mat(num_s,lmat,L1,L2,L3,M1,M2,M3,
     1 N1,N2,N3,charges,dens,flux_mat,force_mat) 
      use penta_kind_mod
      use phys_const
      use vmec_var_pass
      implicit none
      !dummy variables
      integer(iknd), intent(in) :: num_s
      real(rknd), intent(in), dimension(num_s*2,num_s*2) :: lmat
      real(rknd), intent(in), dimension(num_s) :: L1,L2,L3,M1,M2,M3,N1,
     1 N2,N3,charges, dens
      real(rknd),intent(out),dimension(num_s*2,num_s*2) :: flux_mat
      real(rknd),intent(out),dimension(num_s*2,num_s*2+1) :: force_mat
      !local variables
      integer(iknd) :: ispec1, ispec2, ii1, ii2, ind1, ind2
      real(rknd), dimension(2,2) :: Ntmp, Ntmp_inv, Mtmp, Ltmp, ll_tmp,
     1  tmp1,tmp2,tmp3,tmp4,tmp5


      !define flux and force matrices
      do ispec1=1,num_s   !row (species a)
        do ispec2=1,num_s !column (species b)
          !whole matrix indices
          ind1=(ispec1 - 1)*2+1
          ind2=(ispec2 - 1)*2+1 

          !define thermal coefficient matrices
          Ntmp(1,1:2)=[N1(ispec2), N2(ispec2)]
          Ntmp(2,1:2)=[N2(ispec2), N3(ispec2)]
          call mat_2by2_inverse(Ntmp,Ntmp_inv)  !inv(N)
          Mtmp(1,1:2)=[M1(ispec2), M2(ispec2)]
          Mtmp(2,1:2)=[M2(ispec2), M3(ispec2)]
          Ltmp(1,1:2)=[L1(ispec2), L2(ispec2)]
          Ltmp(2,1:2)=[L2(ispec2), L3(ispec2)]
          !define friction coefficient matrices
          ll_tmp(1,1:2)=[lmat(ind1,ind2),lmat(ind1,ind2+1)]
          ll_tmp(2,1:2)=[lmat(ind1+1,ind2),lmat(ind1+1,ind2+1)]
     
          !perform matrix algebra
          tmp1=matmul(Mtmp,Ntmp_inv)
          tmp2=-matmul(ll_tmp,Ntmp_inv)*Bsq
          tmp3=matmul(Ntmp_inv,Ltmp)
          tmp4=matmul(Mtmp,tmp3)
          tmp5=-matmul(ll_tmp,tmp3)*Bsq
        
          !matrix assignment (loop over 2x2 matrices)
          if (ispec1==ispec2) then    !diagonal terms
            do ii1=1,2
              do ii2=1,2
                flux_mat(ind1+ii1-1,ind2+ii2-1)=tmp1(ii1,ii2)
     1            +tmp2(ii1,ii2)
                force_mat(ind1+ii1-1,ind2+ii2-1)=tmp4(ii1,ii2)
     1            +tmp5(ii1,ii2)-Ntmp(ii1,ii2)
              enddo
            enddo 
            else                      !off diagonal terms
            do ii1=1,2
              do ii2=1,2
                flux_mat(ind1+ii1-1,ind2+ii2-1)=tmp2(ii1,ii2)
                force_mat(ind1+ii1-1,ind2+ii2-1)=tmp5(ii1,ii2)
              enddo
            enddo 
          endif !diagonal check

        enddo !column
      enddo  !row

      !add last column of force_mat (XE)
      do ispec1=1,num_s
        ii1=1+(ispec1-1)*2
        ii2=2+(ispec1-1)*2

        force_mat(ii1,num_s*2+1)=dens(ispec1)
     1    *charges(ispec1)*dsqrt(Bsq)
        force_mat(ii2,num_s*2+1)=0._rknd
      enddo
      end subroutine define_flux_force_mat
c
c-----------------------------------------------------------------------
c   This subroutine inverts a 2x2 matrix
c   -mod by JL for speed
c-----------------------------------------------------------------------
c
      subroutine mat_2by2_inverse(a,a_inverse)
      use penta_kind_mod
      implicit none
      real(rknd), dimension(2,2):: a, a_inverse
      real(rknd) :: denom
      denom = a(1,1)*a(2,2) - a(1,2)*a(2,1)
      if(denom .eq. 0.) then
       write(*,'("error: singular matrix")')
       write(*,*) a(1,1),a(1,2)
       write(*,*) a(2,1),a(2,2)
       stop 25
      endif
      a_inverse(1,1) = a(2,2)
      a_inverse(2,2) = a(1,1)
      a_inverse(1,2) = -a(1,2)
      a_inverse(2,1) = -a(2,1)

      a_inverse=a_inverse/denom
      return
      end subroutine mat_2by2_inverse

c
c-----------------------------------------------------------------------
c   -This subroutine calculates the spline fit coefficients to the 
c    L,M,N* files and passes them through coeff_var_pass
c-----------------------------------------------------------------------
c
      subroutine fit_coeffs
      use penta_kind_mod
      use coeff_var_pass
      use parameter_pass
      use bspline
      implicit none
      real(rknd), dimension(:) , allocatable :: cmul, efield,
     1 log_efield

      !allocate spline arrays
      allocate(xt_cl(num_cl+kcord),xt_el(num_el+keord))
      allocate(c_spll(num_cl,num_el))
      allocate(xt_cm(num_cm+kcord),xt_em(num_em+keord))
      allocate(c_splm(num_cm,num_em))
      allocate(xt_cn(num_cn+kcord),xt_en(num_en+keord))
      allocate(c_spln(num_cn,num_en))

      !L star

      !if logarithmic interpolation is to be used, use a normalized efield
      if ( log_interp ) then
        allocate(log_efield(num_el),efield(num_el),cmul(num_cl))
        cmul=dlog(cmul_ls)
        log_efield=dlog(efield_ls)
        efield=(log_efield-minval(log_efield))
     1    /(maxval(log_efield)-minval(log_efield))
        eminl=minval(log_efield); emaxl=maxval(log_efield)
        deallocate(log_efield)
      else
        cmul=cmul_ls
        efield=efield_ls
        eminl=minval(efield); emaxl=maxval(efield)
      endif
      cminl=minval(cmul); cmaxl=maxval(cmul)

      call dbsnak(num_el,efield,keord,xt_el)    !compute 'not-a-knot' sequence
      call dbsnak(num_cl,cmul,kcord,xt_cl)      !compute 'not-a-knot' sequence
      call dbs2in(num_cl,cmul,num_el,efield,coef2d_ls, ! fit
     1            num_cl,kcord,keord,xt_cl,xt_el,c_spll)    
      deallocate(efield,cmul)

      !M star

      !if logarithmic interpolation is to be used, use a normalized efield
      if ( log_interp ) then
        allocate(log_efield(num_em),efield(num_em),cmul(num_cm))
        cmul=dlog(cmul_ms)
        log_efield=dlog(efield_ms)
        efield=(log_efield-minval(log_efield))
     1         /(maxval(log_efield)-minval(log_efield))
        eminm=minval(log_efield); emaxm=maxval(log_efield)
        deallocate(log_efield)
      else
        cmul=cmul_ms
        efield=efield_ms
        eminm=minval(efield); emaxm=maxval(efield)
      endif
      cminm=minval(cmul); cmaxm=maxval(cmul)

      call dbsnak(num_em,efield,keord,xt_em)    !compute 'not-a-knot' sequence
      call dbsnak(num_cm,cmul,kcord,xt_cm)      !compute 'not-a-knot' sequence
      call dbs2in(num_cm,cmul,num_em,efield,coef2d_ms, ! fit
     1          num_cm,kcord,keord,xt_cm,xt_em,c_splm)    
      deallocate(efield,cmul)

      !N star

      !if logarithmic interpolation is to be used
      if ( log_interp ) then
        allocate(log_efield(num_en),efield(num_en),cmul(num_cn))
        cmul=dlog(cmul_ns)
        log_efield=dlog(efield_ns)
        efield=(log_efield-minval(log_efield))
     1         /(maxval(log_efield)-minval(log_efield))
        eminn=minval(log_efield); emaxn=maxval(log_efield)
        deallocate(log_efield)
      else
        cmul=cmul_ns
        efield=efield_ns
        eminn=minval(efield); emaxn=maxval(efield)
      endif
      cminn=minval(cmul); cmaxn=maxval(cmul)

      call dbsnak(num_en,efield,keord,xt_en)    !compute 'not-a-knot' sequence
      call dbsnak(num_cn,cmul,kcord,xt_cn)      !compute 'not-a-knot' sequence
      call dbs2in(num_cn,cmul,num_en,efield,coef2d_ns, ! fit
     1            num_cn,kcord,keord,xt_cn,xt_en,c_spln)    
      deallocate(efield,cmul)
      end subroutine fit_coeffs
c
c-----------------------------------------------------------------------
c   This subroutine defines the thermal transport coefficients
c     by convolving the monoenergetic transport coefficients over
c     the local maxwellian using the subroutine energy_conv.  The
c     calculated thermal arrays are arranged as [electron, I1,I2...]
c   JL 7/2009
c-----------------------------------------------------------------------
c
      subroutine define_thermal_mats(Er,nis,v_ths,ns,masses,qs,
     1  loglambda,L1,L2,L3,M1,M2,M3,N1,N2,N3)
      use penta_kind_mod
      use phys_const
      use coeff_var_pass
      use thermal_pass
      use parameter_pass
      implicit none
      !dummy variables
      integer(iknd) :: nis
      real(rknd) :: Er, loglambda
      real(rknd), dimension(nis+1) :: v_ths,ns,qs,L1,L2,L3,
     1  M1,M2,M3,N1,N2,N3,masses, charges
      !local variables
      integer(iknd) :: ispec, mat_ind, I
      real(rknd), dimension(3) :: Lth, Mth,Nth
      integer(iknd),dimension(nis) ::  ikeep 
      logical :: logopt

      !define |Er| and loglambda
      abs_Er=dabs(Er)
      loglam=loglambda

      !allocate vectors
      allocate(nb(nis),vtb(nis),qb(nis))

      !loop over species and integrate coefficients
      do ispec = 1,nis+1
        !test species parameters
        ma=masses(ispec)
        qa=qs(ispec)
        na=ns(ispec)
        vta=v_ths(ispec)

        !field species parameters
        !ikeep=(/1:ispec-1, ispec+1:nis+1/)  !all indices not equal to ispec
        ikeep=(/(I,I=1,ispec-1), (I,I=ispec+1,nis+1)/)  !all indices not equal to ispec
        qb=qs(ikeep)
        nb=ns(ikeep)
        vtb=v_ths(ikeep)

        !integrate La
        logopt=.false.      !don't use the log of the coefficient
        call energy_conv(num_cl,num_el,xt_cl,xt_el,c_spll,logopt,Lth,
     1    eminl, emaxl, cminl, cmaxl)
        L1(ispec) = Lth(1);L2(ispec) = Lth(2);L3(ispec) = Lth(3)

        !integrate Ma
        logopt=.true.      !use the log of the coefficient
        call energy_conv(num_cm,num_em,xt_cm,xt_em,c_splm,logopt,Mth,
     1    eminm, emaxm, cminm, cmaxm)
        M1(ispec) = Mth(1);M2(ispec) = Mth(2);M3(ispec) = Mth(3)
 
        !integrate Na
        logopt=.false.      !don't use the log of the coefficient
        call energy_conv(num_cn,num_en,xt_cn,xt_en,c_spln,logopt,Nth,
     1    eminn, emaxn, cminn, cmaxn)
        N1(ispec) = Nth(1);N2(ispec) = Nth(2);N3(ispec) = Nth(3)

      enddo !species loop

      deallocate(nb,vtb,qb)
      end subroutine define_thermal_mats
c
c-----------------------------------------------------------------------
c   This subroutine convolves the monoenergetic transport coefficients over
c     the local maxwellian.
c   JL 7/2009 modified from Don's DKES_coef
c-----------------------------------------------------------------------
c
      subroutine energy_conv(num_c,num_e,xt_c,xt_e,c_spl,logopt,coeff,
     1 emin,emax,cmin,cmax)  
      use penta_kind_mod
      use phys_const
      use intfun_pass
      use thermal_pass
      use parameter_pass
	use penta_functions_mod
      implicit none
      !dummy variables
      integer(iknd), intent(in) :: num_c, num_e
      logical, intent(in) :: logopt
      real(rknd), intent(in) :: xt_c(num_c+kcord), xt_e(num_e+keord)
      real(rknd), intent(in) :: c_spl(num_c+kcord,num_e+keord)
      real(rknd), intent(in) :: emin,emax,cmin,cmax
      real(rknd), intent(out), dimension(3) :: coeff
      !local variables
      integer(iknd) :: err, jind, iK
      integer(iknd), dimension(3) :: jvals
      real(rknd) :: int_flag, nofun, coeff_tmp, dK, Ktest
      !external functions
      !real(rknd), external :: intfun
      
      allocate(xt_c_pass(num_c+kcord))
      allocate(xt_e_pass(num_e+keord))
      allocate(c_spl_pass(num_c+kcord,num_e+keord))

      !assign variables for passing
      num_c_pass=num_c
      num_e_pass=num_e
      log_coeff=logopt
      xt_c_pass=xt_c
      xt_e_pass=xt_e
      c_spl_pass=c_spl
      emin_pass=emin
      emax_pass=emax
      cmin_pass=cmin
      cmax_pass=cmax

      !loop over j=1,3 and integrate
      jvals=[1_iknd,2_iknd,3_iknd]

      do jind=1,3
       
        !test j
        jval=jvals(jind)

        !either use quanc8 or rectangular approximation
        if (use_quanc8) then         !quanc8
          call quanc8(intfun, Kmin, Kmax,epsabs, epsrel,
     1               coeff_tmp, err, nofun, int_flag)
          !check for integrator error
          if(abs(int_flag) .gt. 1.e-16_rknd) then
            write(*,*) 'warning: quanc8 error, flag = ',int_flag
            write(*,*) 'May be ok, proceeding'
          endif
        else                       !rectangular
          coeff_tmp=0._rknd
          dK=(Kmax-Kmin)/(num_K_steps-1._rknd)
          do iK=1,num_K_steps
            Ktest=Kmin+(real(iK,rknd)-1._rknd)*dK   !test K
            coeff_tmp=intfun(Ktest)+coeff_tmp       !sum(coeff(K))
          enddo
        coeff_tmp=coeff_tmp*dK                    !sum(coeff(K))*dK
        endif !integrator choice
    
        coeff(jind)=coeff_tmp*na*2._rknd*vta*ma/dsqrt(pi)
      enddo !jvals

      deallocate(xt_c_pass,xt_e_pass,c_spl_pass)

      end subroutine energy_conv
c
c-----------------------------------------------------------------------
c   This subroutine calculates the perpendicular collision frequency
c       of species a on species b.  The external function derf is 
c     used to get the double precision error function, but if this
c     is included in a library the function can be omitted.
c   JL 7/2009
c-----------------------------------------------------------------------
c
      subroutine calc_perp_coll_freq(va,xb,qa,qb,ma,nb,Kb,loglambda
     1  ,vta,Ka,nu_perp_a)
      use penta_kind_mod
      use phys_const
	use penta_functions_mod
      implicit none
      !dummy variables
      real(rknd),intent(in) :: va,xb,qa,qb,ma,nb,Kb,loglambda,vta,Ka
      real(rknd),intent(out) :: nu_perp_a
      !local variables
      real(rknd) :: va3, nu_ab, qa2, qb2, H_ab, one
      !external functions
      !real(rknd), external :: derf

      !constant
      one=1._rknd

      !define powers for speed
      va3=Ka*va*vta**2_iknd
      qa2=qa**2_iknd
      qb2=qb**2_iknd

      !calculate reference collision freq of a on b 
      nu_ab=nb*qa2*qb2*loglambda/
     1       (ma**2_iknd*va3*4._rknd*pi*eps0**2_iknd)
      !calculate potential
      H_ab=(one-one/(2._rknd*Kb))*derf(xb) + dexp(-Kb)/(xb*dsqrt(pi))
      !collision frequency
      nu_perp_a=H_ab*nu_ab
      end subroutine calc_perp_coll_freq
c
c-----------------------------------------------------------------------
c   This subroutine forms the thermodynamic force vector.
c   JL 7/2009
c-----------------------------------------------------------------------
c
      subroutine form_Xvec(Er_test,Z_ion,B_Eprl,ns,Xvec)
      use penta_kind_mod
      use phys_const
      use pprof_pass
      use vmec_var_pass
      implicit none
      !dummy variables
      integer(iknd),intent(in) :: ns
      real(rknd),intent(in) :: Er_test, B_Eprl
      real(rknd),intent(in), dimension(ns) :: Z_ion
      real(rknd),intent(out),dimension(ns*2+3) :: Xvec
        
      !   Species charge sign has been included, elem_charge has been
      !   factored out of Er term because Ta is in eV.

      !electrons
      Xvec(1) = -Te*elem_charge*(dnedr/ne + dTedr/Te + Er_test/Te)
      Xvec(2) = -elem_charge*dTedr
      !ions
      Xvec(3:(ns-1)*2+3:2) = -elem_charge*Ti*
     1  (dnidr/ni + dTidr/Ti - Z_ion*Er_test/Ti)
      Xvec(4:(ns-1)*2+4:2) = -elem_charge*dTidr
      !E_||
      Xvec(ns*2+3)=B_Eprl/dsqrt(Bsq);

      return
      end subroutine form_Xvec
c
c-----------------------------------------------------------------------
c   !Subroutine to find the inverse of a square matrix
c   !Author : Louisda16th a.k.a Ashwith J. Rego
c   !Reference : Algorithm has been well explained in:
c   !http://math.uww.edu/~mcfarlat/inverse.htm           
c   !http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
c   modified by JL 7/2009 for kind precision
c-----------------------------------------------------------------------
c
      SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
      use penta_kind_mod
      IMPLICIT NONE
      !Declarations
      INTEGER(iknd), INTENT(IN) :: n
      INTEGER(iknd), INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
      REAL(rknd), INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
      REAL(rknd), INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
    
      LOGICAL :: FLAG = .TRUE.
      INTEGER(iknd) :: i, j, k, l
      REAL(rknd) :: m
      REAL(rknd), DIMENSION(n,2*n) :: augmatrix !augmented matrix
    
    !Augment input matrix with an identity matrix
      DO i = 1, n
        DO j = 1, 2*n
            IF (j <= n ) THEN
                augmatrix(i,j) = matrix(i,j)
            ELSE IF ((i+n) == j) THEN
                augmatrix(i,j) = 1._rknd
            Else
                augmatrix(i,j) = 0._rknd
            ENDIF
        END DO
      END DO
    
      !Reduce augmented matrix to upper traingular form
      DO k =1, n-1
        IF (augmatrix(k,k) == 0._rknd) THEN
            FLAG = .FALSE.
            DO i = k+1, n
                IF (augmatrix(i,k) /= 0._rknd) THEN
                    DO j = 1,2*n
                        augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                    END DO
                    FLAG = .TRUE.
                    EXIT
                ENDIF
                IF (FLAG .EQV. .FALSE.) THEN
                    PRINT*, "Matrix is non - invertible"
                    inverse = 0._rknd
                    errorflag = -1_iknd
                    return
                ENDIF
            END DO
        ENDIF
        DO j = k+1, n           
            m = augmatrix(j,k)/augmatrix(k,k)
            DO i = k, 2*n
                augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
            END DO
        END DO
      END DO
    
      !Test for invertibility
      DO i = 1, n
        IF (augmatrix(i,i) == 0) THEN
            PRINT*, "Matrix is non - invertible"
            inverse = 0._rknd
            errorflag = -1_iknd
            return
        ENDIF
      END DO
    
      !Make diagonal elements as 1
      DO i = 1 , n
        m = augmatrix(i,i)
        DO j = i , (2 * n)              
               augmatrix(i,j) = (augmatrix(i,j) / m)
        END DO
      END DO
    
      !Reduced right side half of augmented matrix to identity matrix
      DO k = n-1, 1, -1
        DO i =1, k
        m = augmatrix(i,k+1)
            DO j = k, (2*n)
                augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
            END DO
        END DO
      END DO              
    
      !store answer
      DO i =1, n
        DO j = 1, n
            inverse(i,j) = augmatrix(i,j+n)
        END DO
      END DO
      errorflag = 0_iknd
      END SUBROUTINE FINDinv