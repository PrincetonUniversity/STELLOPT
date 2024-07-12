
      subroutine woflam(trigsu, trigsv, ifaxu, ifaxv, irho)

c  Evaluate the lambda-dependent part of the functions W1(lambda)
c  for all LAMBDA.
c  Calculates the sum over n and m of (mR+nS)/(m-nq) * dalphamn/dlambda
c      * <V||/V (mn)>/<V||/V>*lambda*(-2).
c  Here, m is the poloidal mode number and n is the toroidal mode
c  number/periods.
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: irho
      integer, dimension(13) :: ifaxu, ifaxv
      real(rprec), dimension(3*nthetah/2 + 1) :: trigsu
      real(rprec), dimension(2*nzetah) :: trigsv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1,
     1   d18 = 1.0e-18_dp, xlam = 0.96_dp
      integer, parameter :: n_lam_coarse = 97, n_lam = 137

c  as discussed below, the mesh is split with 96 + 40 intervals

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l, j, i, n, m
      real(rprec), dimension(:,:), allocatable :: a11
      complex(rprec), dimension(:,:), allocatable ::
     1   alphamn, vmn
      real(rprec) :: qn, numer, avg_vpov, denom
      real(rprec), dimension(n_lam) :: xlam_f
      real(rprec), dimension(n_lam-1) :: xlam_h
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , EXTERNAL :: sumit
C-----------------------------------------------
c
c
      allocate (a11(nthetah+2,nzetah),
     5   alphamn(-mbuse:mbuse,0:nbuse),
     8   vmn(-mbuse:mbuse,0:nbuse), stat=l)

      if (l .ne. 0) stop 'allocation error in woflam'

c-----------------------------------------------------------------------
c///////////////////////////////////////////////////////////////////////
c-----------------------------------------------------------------------

c  Meshing for lamda integral.  The lambda integral has a near singular
c  value at lambda = 1.  To evaluate, construct dalpha/dlamda on full
c  mesh.  Evaluate the rest on the full mesh.  Becasue of the (almost)
c  singularity, split the mesh into coarse (0 to 0.96) and fine (.96 to
c  1.0) components.  This was tested on both qas and qos configuratons
c  and gives results that are a percent or better on the integral and
c  ten to twenty times better for the total current.  To save storage,
c  there is only one loop on lambda.  All grids are local.

c-----------------------------------------------------------------------
c///////////////////////////////////////////////////////////////////////
c-----------------------------------------------------------------------


      if (isymm0 .ne. 0) then
         aiterm1(irho) = zero
         return
      endif

c  First contruct the meshes.  This scheme is hardwired.

      xlam_f(1) = zero
      xlam_h(1) = 0.005_dp

      do l = 2, n_lam_coarse-1
         xlam_f(l) = xlam_f(l-1) + 0.01_dp
         xlam_h(l) = xlam_h(l-1) + 0.01_dp
      enddo

      xlam_f(n_lam_coarse) = 0.96_dp               !!End coarse mesh/Begin fine mesh
      xlam_h(n_lam_coarse) = 0.9605_dp             !!Begin of fine half-mesh

      do l = n_lam_coarse+1, n_lam-1
         xlam_f(l) = xlam_f(l-1) + 0.001_dp
         xlam_h(l) = xlam_h(l-1) + 0.001_dp
      enddo

      xlam_f(n_lam) = one

c  Loop over LAMBDA = 0 to 1, calculating alpha on the full mesh, and
c  the rest of the integral on the half mesh.  The jacobian for alpha
c  puts in a b2avg but one of the b values in the denominator cancels the
c  b/bmax in alpha.  This form of the expression paralles that in the
c  hamada paper but with flux surface averages for the alpha.  This results
c  in a form that does not have the beta found in the boozer paper
c  see refereces at beginning.

      do l = 1, n_lam
         a11(:nthetah,:nzetah) =
     2      sqrt(abs(one - xlam_f(l)*bfield(:nthetah,:nzetah))+D18)
     3      *b2avg(irho)/bmax1(irho)**2/bfield(:nthetah,:nzetah)
         a11(nthetah+1,:nzetah) = zero
         a11(nthetah+2,:nzetah) = zero
         call do_fft (a11, fmn, trigsu, trigsv, ifaxu, ifaxv, nthetah,
     1      nzetah, mbuse, nbuse)
         if(l .eq. 1) then
            alphamn = fmn
            cycle                       !need to calculate to alphas and difference
         end if                         !before a term for the integral can be evaluated.


c  save alpha at lambda = 1

         if(l .eq. n_lam) alpha1mn = fmn

c  form d_alphalmn  at this point we have the value of alpha for l in fmn
c  and the previous value of alpha (l-1) in alphamn.  Store the difference in
c  fmn using vmn as a temp variable to hold alpha, then update alpha for
c  the next cycle.

         vmn = fmn
         fmn = fmn - alphamn
         alphamn = vmn

c- Now do the fft to get <exp(i(m*theta-n*zeta)V||/V> on the half mesh

         a11(:nthetah,:nzetah) =
     1      sqrt(abs(one - xlam_h(l-1)*bfield(:nthetah,:nzetah)) + D18)
     2      *(b2avg(irho)/bfield(:nthetah,:nzetah)*bmax1(irho))**2

         a11(nthetah+1,:nzetah) = zero
         a11(nthetah+2,:nzetah) = zero

         call do_fft (a11, vmn, trigsu, trigsv, ifaxu, ifaxv, nthetah,
     1      nzetah, mbuse, nbuse)

c  This needs to be divided by <V||/V>.  But this is exactly the real part of the
c  0,0 term of this transform.
         avg_vpov = real(vmn(0,0))


c- for only those harmonics that are going to be used in the sum.

         qn = periods*qsafety(irho)*zetasign
         rfmn(0,0) = zero
         do m = -mbuse, mbuse
            do n = 0, nbuse

               denom = m + n*qn
               if (n.ne.0 .or. m.ne.0) then

                  numer = denom/(denom**2 + (damp_bs*m)**2)

                  rfmn(m,n) = (m*capr(irho) + n*periods*caps(irho))*
     1             real(fmn(m,n)*vmn(m,n))*(-2*numer)
               endif
            end do
         end do
c
c- The sum in W is the SUM over all n and m (excluding (0,0)) of rfmn.
c  But only the nonegative m have been calculated and stored.  Therefore
c  the sum becomes
c    2 * SUM over m = -mbuse to mbuse and n = 1 to nbuse of rfmn
c     plus SUM over n = -mbuse to mbuse of fn0
c     minus f(0,0)

         w1(l-1) = xlam_h(l-1)*sumit(rfmn,mbuse,nbuse)/avg_vpov

c- End loop over lambda.

      end do
c  evaluate integral

      aiterm1(irho) = sum(w1)
      aiterm1(irho) = -0.75_dp*aiterm1(irho)*qsafety(irho)
     1    /ftrapped(irho) * (one + aiogar(irho)/qsafety(irho))

      deallocate (a11, alphamn, vmn, stat=l)

      end subroutine woflam
