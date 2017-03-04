      MODULE fourier
!     
!     WRITTEN 08-09-06 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: CONVERTS QUANTITIES FROM FOURIER SPACE TO REAL SPACE AND VICEVERSA
!
!       NOTE: CALL FIXARRAY must be used once before calling any of the Fourier subroutines
!       to calculate all necessary cosine/sine factors on fixed mesh of collocation angles
!       NS, NTOR, MPOL and NFP used as in the wout.file
!
      USE stel_kinds
      USE stel_constants
      USE island_params, ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i,           &
                   mpol=>mpol_i, ntor=>ntor_i, nfp=>nfp_i
      USE timer_mod
      USE descriptor_mod, ONLY: INHESSIAN, DIAGONALDONE
      IMPLICIT NONE
!
!     VARIABLE DECLARATIONS  
!
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: orthonorm
      CONTAINS
 
#if defined(SKS)
      SUBROUTINE TOIJSP_PAR(XMN, XUV, UFLAG, VFLAG, IPARITY, IHALF, nsmin, nsmax)
!
!       DESCRIPTION:  This subroutine moves a quantity X to real space by summing over its Fourier harmonics X_MN.
!       IF UFLAG = 1, the first poloidal derivative of that quantity is obtained
!       IF VFLAG = 1, the first toroidal derivative of that quantity is obtained
!       IF both flags are equal to one, it gives the joint poloidal, toroidal derivative
!       IPARITY = 0, means the quantity X (NOT any of its derivatives!) is COSINE (even); = 1, X is SINE (odd)
!       IHALF = 0, means quantity is on the radial full mesh; =1, on the half radial mesh
!
      USE stel_kinds
      USE nscalingtools, ONLY: startglobrow, endglobrow
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN):: ihalf                         ! IHALF = 0, FULL (radial) MESH); = 1, HALF MESH     
      INTEGER, INTENT(IN):: nsmin, nsmax 
      REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN):: xmn
      REAL(rprec), DIMENSION(ntheta,nzeta,ns), INTENT(OUT):: xuv
      INTEGER, INTENT(IN):: uflag, vflag                  ! UFLAG/VFLAG = order of poloidal/toroidal derivative
      INTEGER, INTENT(IN):: iparity                       ! IPARITY = 0, cosine (EVEN); = 1, sine (ODD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: m, n, jk, lk, js
      INTEGER :: in 
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: work0, work1, work2, xuv1
      REAL(dp) :: skston, skstoff
      INTEGER :: counter, nx
      INTEGER, SAVE :: PTOIJSPPASS=1
!-----------------------------------------------
      CALL second0(skston)

      ALLOCATE(work0(0:mpol,-ntor:ntor), xuv1(ntheta,nzeta),            &         
               work1(ntheta,-ntor:ntor), work2(ntheta,-ntor:ntor),      &
               stat=lk)
      IF (lk .ne. 0) STOP 'Allocation error in TOIJSP_PAR'
      in = 1 +ihalf

      NS_LOOP: DO js = nsmin, nsmax

      work0(:,:) = orthonorm(0:mpol,-ntor:ntor)*xmn(:,:,js)
      work1 = zero; work2 = zero 

      DO jk = 1, ntheta
        DO m = 0, mpol                                           ! First DO sum over poloidal modes
          IF (uflag == 0) THEN
            work1(jk,:) = work1(jk,:) + work0(m,:)*cosmu(jk, m)
            work2(jk,:) = work2(jk,:) + work0(m,:)*sinmu(jk, m)
          ELSE IF (uflag == 1) THEN                              ! Poloidal derivative requested: USE COSMUM, SINMUM
            work1(jk,:) = work1(jk,:) - work0(m,:)*sinmum(jk, m)
            work2(jk,:) = work2(jk,:) + work0(m,:)*cosmum(jk, m)
          ELSE
            STOP 'UFLAG > 1 in TOIJSP subroutine'
          ENDIF
        ENDDO
      ENDDO
 
      xuv1 = zero                                                 ! Make sure that SPATIAL variable is zeroed

      DO lk = 1, nzeta

        IF (iparity == 0) THEN                                   ! Do first N=0 mode

          IF (vflag == 0) THEN
            xuv1(:,lk) = work1(:,0) 
          ELSEIF (vflag == 1) THEN
            xuv1(:,lk) = zero
          ELSE        
            STOP 'VFLAG > 1 in TOIJSP subroutine'            
          ENDIF

        ELSE

          IF (vflag == 0) THEN
            xuv1(:,lk) = work2(:,0) 
          ELSEIF (vflag == 1) THEN
            xuv1(:,lk) = zero
          ELSE
            STOP 'VFLAG > 1 in TOUVSP subroutine' 
          ENDIF          

        ENDIF          

        DO n = 1, ntor                                            ! Then sum over N>0 and N<0 toroidal modes                                                     
          IF (iparity == 0) THEN                                  ! COSINE series

            IF (vflag == 0) THEN
              xuv1(:,lk) =  xuv1(:,lk)                          &
              + (work1(:,n) + work1(:,-n))*cosnv(lk, n)         &
              - (work2(:,n) - work2(:,-n))*sinnv(lk, n)
            ELSEIF (vflag == 1) THEN                              ! First toroidal derivative requested
              xuv1(:,lk) =  xuv1(:,lk)                          &
              - (work1(:,n) + work1(:,-n))*sinnvn(lk, n)        &
              - (work2(:,n) - work2(:,-n))*cosnvn(lk, n)
            ELSE
              STOP 'VFLAG > 1 in TOIJSP subroutine'
            ENDIF

          ELSE                                                    ! SINE series

            IF (vflag == 0) THEN
              xuv1(:,lk) =  xuv1(:,lk)                          &
              + (work1(:,n) - work1(:,-n))*sinnv(lk, n)         &
              + (work2(:,n) + work2(:,-n))*cosnv(lk, n)
            ELSEIF (vflag == 1) THEN                              ! First toroidal derivative requested
              xuv1(:,lk) =  xuv1(:,lk)                          &
              + (work1(:,n) - work1(:,-n))*cosnvn(lk, n)        &
              - (work2(:,n) + work2(:,-n))*sinnvn(lk, n)
            ELSE
              STOP 'VFLAG > 1 in TOIJSP subroutine'
            ENDIF

          ENDIF
        ENDDO
      ENDDO

      xuv(:,:,js) = xuv1(:,:)

      END DO NS_LOOP

      DEALLOCATE(work0, work1, work2, xuv1)

      CALL second0(skstoff)
      IF(DIAGONALDONE.AND.INHESSIAN) toijsp_time=toijsp_time+(skstoff-skston)
      time_toijsp = time_toijsp + (skstoff-skston)

      PTOIJSPPASS = PTOIJSPPASS+1

      END SUBROUTINE TOIJSP_PAR

      SUBROUTINE TOMNSP_PAR(XUV, XMN, IPARITY, IHALF, nsmin, nsmax)
!
!       Description: This subroutine moves a quantity X to Fourier space producing its harmonics X_MN by
!         summing over all the presribed (toroidal and poloidal) collocation points:
!           theta_j = j*pi/M, (j=0,...,M);   zeta_k = k*2*pi/(2N+1), k = 0,...., 2N
!         where M = mpol - 1 and N = ntor.
!       IPARITY = 0, means the quantity X (NOT any of its derivatives) is even; = 1, X is odd
!       IHALF = 0, means quantity is on the radial full mesh; =1, on the half radial mesh
!
      USE stel_kinds
      USE nscalingtools, ONLY: startglobrow, endglobrow
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ihalf                                        ! IHALF = 0, FULL (radial) MESH); = 1, HALF MESH
      INTEGER, INTENT(IN) :: nsmin, nsmax
      REAL(rprec), DIMENSION(ntheta,nzeta,ns), INTENT(IN):: xuv
      REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(OUT):: xmn
      INTEGER, INTENT(IN):: iparity                                       ! IPARITY = 0, cosine (EVEN); = 1; sine (ODD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: m, n, jk, lk, in, js
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: work1, work2, xuv1, xmn1
      REAL(dp), ALLOCATABLE, DIMENSION(:)   :: t1, t2, work3, work4
      REAL(dp) :: skston, skstoff
      INTEGER :: ii, jj, kk 
      INTEGER :: counter, nx
      INTEGER, SAVE :: PTOMNSPPASS=1
!-----------------------------------------------

      CALL second0(skston)

      in = 1+ihalf
      ALLOCATE(work1(0:mpol,nzeta), work2(0:mpol,nzeta),                &
               work3(0:mpol), work4(0:mpol),                            &
               xuv1(ntheta,nzeta),  xmn1(0:mpol,-ntor:ntor),            &
               t1(nzeta), t2(nzeta), stat=ii)
      IF (ii .ne. 0) STOP 'Allocation error in tomnsp_par'

      NS_LOOP: DO js = nsmin, nsmax

      xuv1(:,:) = xuv(:,:,js)
      xmn1      = 0

      DO m = 0, mpol
        t1 = 0; t2 = 0
        DO jk = 1, ntheta                                                 ! First, add over poloidal collocation angles
          t1 = t1 + xuv1(jk,:)*cosmui(jk, m)                              ! Note use of cosmuI/sinmuI; include norm. factors
          t2 = t2 + xuv1(jk,:)*sinmui(jk, m)
        ENDDO
        work1(m,:) = t1;  work2(m,:) = t2
      ENDDO


      DO n = 0, ntor  
        DO lk = 1, nzeta                                                  ! Then, add over toroidal collocation angles
          IF (iparity == 0) THEN                                          ! COSINE series              
            work3(:) = work1(:,lk)*cosnv(lk, n)
            work4(:) = work2(:,lk)*sinnv(lk, n)
            xmn1(:,n) =  xmn1(:,n) + work3(:) - work4(:)
            IF (n .NE. 0) THEN
              xmn1(1:,-n) = xmn1(1:,-n) + work3(1:) + work4(1:)
            ENDIF
          ELSE                                                            ! SINE series
            work3(:) = work1(:,lk)*sinnv(lk, n)
            work4(:) = work2(:,lk)*cosnv(lk, n)
            xmn1(:,n) =  xmn1(:,n) + work3(:) + work4(:)
            IF (n .NE. 0) THEN
              xmn1(1:,-n) = xmn1(1:,-n) - work3(1:) + work4(1:)
            ENDIF
          ENDIF
        ENDDO
      ENDDO 

      xmn1(0,-ntor:-1) = zero                                             ! Redundant: To avoid counting (0,-n) for non-zero n.

      xmn(:,:,js) = orthonorm(0:mpol,-ntor:ntor)*xmn1(:,:)

      END DO NS_LOOP

      DEALLOCATE(work1, work2, work3, work4, xuv1, xmn1, t1, t2)

      CALL second0(skstoff)
      IF(DIAGONALDONE.AND.INHESSIAN) tomnsp_time=tomnsp_time+(skstoff-skston)

      time_tomnsp = time_tomnsp + (skstoff-skston)

      PTOMNSPPASS = PTOMNSPPASS+1
      END SUBROUTINE TOMNSP_PAR
#endif

        SUBROUTINE GETORIGIN (XUV, MODE, IPARITY)
!
!       DESCRIPTION: Computes the origin value of a half-mesh quantity (stores it at index 1)
!                    from the m=mode (cos,sin(mu)-average) of the quantity at js=2
!       IPARITY = 0, means the quantity X (NOT any of its derivatives!) is COSINE (even); = 1, X is SINE (odd)
!       MODE         the poloidal mode number to be retained
        USE stel_kinds
        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER, INTENT(in) :: iparity, mode
        REAL(rprec), DIMENSION(ntheta,nzeta,ns), INTENT(INOUT) :: xuv
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER, PARAMETER :: js1=1, js2=2
        REAL(rprec), DIMENSION(0:mpol,-ntor:ntor) :: xmn
!-----------------------------------------------
        CALL tomnsp2(xuv, xmn, iparity, js2)
        
        IF (mode .gt. 0) xmn(0:mode-1,:) = 0
        IF (mode .lt. mpol) xmn(mode+1:,:) = 0

        CALL toijsp2(xmn, xuv, iparity, js1)

        END SUBROUTINE GETORIGIN

        SUBROUTINE TOIJSP(XMN, XUV, UFLAG, VFLAG, IPARITY, IHALF)
!
!       DESCRIPTION:  This subroutine moves a quantity X to real space by summing over its Fourier harmonics X_MN.
!       IF UFLAG = 1, the first poloidal derivative of that quantity is obtained
!       IF VFLAG = 1, the first toroidal derivative of that quantity is obtained
!       IF both flags are equal to one, it gives the joint poloidal, toroidal derivative
!       IPARITY = 0, means the quantity X (NOT any of its derivatives!) is COSINE (even); = 1, X is SINE (odd)
!       IHALF = 0, means quantity is on the radial full mesh; =1, on the half radial mesh
!
        USE stel_kinds
        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER, INTENT(IN):: ihalf                         ! IHALF = 0, FULL (radial) MESH); = 1, HALF MESH     
        REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns),                   &
          INTENT(IN):: xmn
        REAL(rprec), DIMENSION(ntheta,nzeta,ns),                        &
          INTENT(OUT):: xuv
        INTEGER, INTENT(IN):: uflag, vflag                  ! UFLAG/VFLAG = order of poloidal/toroidal derivative
        INTEGER, INTENT(IN):: iparity                       ! IPARITY = 0, cosine (EVEN); = 1, sine (ODD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: m, n, jk, lk, js
        INTEGER :: in 
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:):: work0,             &
          work1, work2
        REAL(rprec) :: ton, toff, skston, skstoff
!-----------------------------------------------
        CALL second0(ton)
        skston=ton

        ALLOCATE(work0(0:mpol,-ntor:ntor,ns))                    ! RS (10/03/06): Added WORK0 for n=0/m=0 normalization
        in = 1+ihalf

        work0 = xmn
        DO js = in,ns
           work0(:,:,js) = orthonorm(0:mpol,-ntor:ntor)*work0(:,:,js)
        END DO
             
        ALLOCATE(work1(ntheta,-ntor:ntor,ns),                           &
                 work2(ntheta,-ntor:ntor,ns))
        work1 = zero; work2 = zero 
        

        DO jk = 1, ntheta
          DO m = 0, mpol                                         ! First DO sum over poloidal modes
            IF (uflag == 0) THEN
              work1(jk,:,in:ns) = work1(jk,:,in:ns) +                   &
                work0(m,:,in:ns)*cosmu(jk, m)
              work2(jk,:,in:ns) = work2(jk,:,in:ns) +                   &
                work0(m,:,in:ns)*sinmu(jk, m)
            ELSE IF (uflag == 1) THEN                            ! Poloidal derivative requested: USE COSMUM, SINMUM
              work1(jk,:,in:ns) = work1(jk,:,in:ns) -                   &
                work0(m,:,in:ns)*sinmum(jk, m)
              work2(jk,:,in:ns) = work2(jk,:,in:ns) +                   &
                work0(m,:,in:ns)*cosmum(jk, m)
            ELSE
              STOP 'UFLAG > 1 in TOIJSP subroutine'
            ENDIF
          ENDDO
        ENDDO
        DEALLOCATE(work0)


        xuv = zero                                                 ! Make sure that SPATIAL variable is zeroed

        DO lk = 1, nzeta
          
          IF (iparity == 0) THEN                                   ! Do first N=0 mode
          
            IF (vflag == 0) THEN
               xuv(:,lk,in:ns) = work1(:,0,in:ns) 
            ELSEIF (vflag == 1) THEN
               xuv(:,lk,in:ns) = zero
            ELSE        
               STOP 'VFLAG > 1 in TOIJSP subroutine'            
            ENDIF
          
          ELSE
           
            IF (vflag == 0) THEN
                xuv(:,lk,in:ns) = work2(:,0,in:ns) 
            ELSEIF (vflag == 1) THEN
                xuv(:,lk,in:ns) = zero
            ELSE
               STOP 'VFLAG > 1 in TOUVSP subroutine' 
            ENDIF          
          
          ENDIF          
          
          DO n = 1, ntor                                          ! Then sum over N>0 and N<0 toroidal modes                                                     
            IF (iparity == 0) THEN                                ! COSINE series

              IF (vflag == 0) THEN
                xuv(:,lk,in:ns) =  xuv(:,lk,in:ns)                      &
                  + (work1(:,n,in:ns) + work1(:,-n,in:ns))*cosnv(lk, n) &
                  - (work2(:,n,in:ns) - work2(:,-n,in:ns))*sinnv(lk, n)
              ELSEIF (vflag == 1) THEN                            ! First toroidal derivative requested
                xuv(:,lk,in:ns) =  xuv(:,lk,in:ns)                      &
                  - (work1(:,n,in:ns) + work1(:,-n,in:ns))*sinnvn(lk, n) &
                  - (work2(:,n,in:ns) - work2(:,-n,in:ns))*cosnvn(lk, n)
              ELSE
                STOP 'VFLAG > 1 in TOIJSP subroutine'
              ENDIF

            ELSE                                                  ! SINE series

              IF (vflag == 0) THEN
                xuv(:,lk,in:ns) =  xuv(:,lk,in:ns)                       &
                  + (work1(:,n,in:ns) - work1(:,-n,in:ns))*sinnv(lk, n)  &
                  + (work2(:,n,in:ns) + work2(:,-n,in:ns))*cosnv(lk, n)
              ELSEIF (vflag == 1) THEN                            ! First toroidal derivative requested
                xuv(:,lk,in:ns) =  xuv(:,lk,in:ns)                       &
                  + (work1(:,n,in:ns) - work1(:,-n,in:ns))*cosnvn(lk, n) &
                  - (work2(:,n,in:ns) + work2(:,-n,in:ns))*sinnvn(lk, n)
              ELSE
                STOP 'VFLAG > 1 in TOIJSP subroutine'
              ENDIF
              
            ENDIF
          ENDDO
        ENDDO

        DEALLOCATE(work1, work2)

        CALL second0(toff)
        skstoff=toff
#if defined(SKS)
        IF(DIAGONALDONE.AND.INHESSIAN) toijsp_time=toijsp_time+(skstoff-skston)
#endif
        time_toijsp = time_toijsp + (toff-ton)

        END SUBROUTINE TOIJSP


        SUBROUTINE TOMNSP(XUV, XMN, IPARITY, IHALF)
!
!       Description: This subroutine moves a quantity X to Fourier space producing its harmonics X_MN by
!         summing over all the presribed (toroidal and poloidal) collocation points:
!           theta_j = j*pi/M, (j=0,...,M);   zeta_k = k*2*pi/(2N+1), k = 0,...., 2N
!         where M = mpol - 1 and N = ntor.
!       IPARITY = 0, means the quantity X (NOT any of its derivatives) is even in mu+nv; = 1, X is odd
!       IHALF = 0, means quantity is on the radial full mesh; =1, on the half radial mesh
!
        USE stel_kinds
        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER, INTENT(IN):: ihalf                         ! IHALF = 0, FULL (radial) MESH); = 1, HALF MESH
        REAL(rprec), DIMENSION(ntheta,nzeta,ns),                        &
          INTENT(IN):: xuv
        REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns),                   &
          INTENT(OUT):: xmn
        INTEGER, INTENT(IN):: iparity                          ! IPARITY = 0, cosine (EVEN); = 1; sine (ODD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: m, n, jk, lk, in, js
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:):: work1, work2,      &
          work3, work4
        REAL(rprec) :: ton, toff, skston, skstoff
!-----------------------------------------------
        CALL second0(ton)
        skston=ton

        ALLOCATE(work1(0:mpol,nzeta,ns), work2(0:mpol,nzeta,ns))
        work1 = zero; work2 = zero

        in = 1+ihalf

        DO m = 0, mpol
          DO jk = 1, ntheta                                            ! First, add over poloidal collocation angles
            work1(m,:,in:ns) = work1(m,:,in:ns) +                       &
              xuv(jk,:,in:ns)*cosmui(jk, m)                       ! Note use of cosmuI/sinmuI; include norm. factors
            work2(m,:,in:ns) = work2(m,:,in:ns) +                       &
              xuv(jk,:,in:ns)*sinmui(jk, m)
          ENDDO
        ENDDO

        xmn = zero                                                  ! Make sure that FOURIER variable is zeroed
        ALLOCATE(work3(0:mpol,nzeta,ns), work4(0:mpol,nzeta, ns))
        work3 = zero; work4 = zero

        DO n = 0, ntor  
          DO lk = 1, nzeta                                      ! Then, add over toroidal collocation angles
            IF (iparity == 0) THEN                              ! COSINE series              
              work3(:,lk,in:ns) = work1(:,lk,in:ns)*cosnv(lk, n)
              work4(:,lk,in:ns) = work2(:,lk,in:ns)*sinnv(lk, n)
              xmn(:,n,in:ns) =  xmn(:,n,in:ns)                          &
                 + work3(:,lk,in:ns) - work4(:,lk,in:ns)
              IF (n .ne. 0) THEN
                xmn(1:,-n,in:ns) = xmn(1:,-n,in:ns)                     &
                 + work3(1:,lk,in:ns) + work4(1:,lk,in:ns)
              ENDIF
            ELSE                                                ! SINE series
              work3(:,lk,in:ns) = work1(:,lk,in:ns)*sinnv(lk, n)
              work4(:,lk,in:ns) = work2(:,lk,in:ns)*cosnv(lk, n)
              xmn(:,n,in:ns) =  xmn(:,n,in:ns)                          &
                + work3(:,lk,in:ns) + work4(:,lk,in:ns)
              IF (n .ne. 0) THEN
                xmn(1:,-n,in:ns) = xmn(1:,-n,in:ns)                     &
                - work3(1:,lk,in:ns) + work4(1:,lk,in:ns)
              ENDIF
            ENDIF
          ENDDO
        ENDDO 
 
        xmn(0,-ntor:-1,:) = zero                               ! Redundant: To avoid counting (0,-n) for non-zero n.

        DO js = in,ns
           xmn(:,:,js) = orthonorm(0:mpol,-ntor:ntor)*xmn(:,:,js)
        END DO
         
        DEALLOCATE(work1, work2, work3, work4)

        CALL second0(toff)
#if defined(SKS)
        skstoff=toff
        IF(DIAGONALDONE.AND.INHESSIAN) tomnsp_time=tomnsp_time+(skstoff-skston)
#endif
        time_tomnsp = time_tomnsp + (toff-ton)

        END SUBROUTINE TOMNSP
      
!SPH 10/15/14
        SUBROUTINE TOMNSP_PEST(XUV, XMN, LMNS, IPARITY, IHALF)
!
!       Description: This subroutine moves a quantity X to Fourier space producing its harmonics X_MN by
!         summing over all the presribed (toroidal and poloidal) collocation points:
!           theta_j = j*pi/M, (j=0,...,M);   zeta_k = k*2*pi/(2N+1), k = 0,...., 2N
!         where M = mpol - 1 and N = ntor.
!       XUV is in VMEC coordinates; the XMN are harmonics in PEST coordinates, where
!              UP(s,u,v) = U + LAMBDA(s,u,v)
!       IPARITY = 0, means the quantity X (NOT any of its derivatives) is even in mu+nv; = 1, X is odd
!       IHALF = 0, means quantity is on the radial full mesh; =1, on the half radial mesh
!
        USE stel_kinds
        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER, INTENT(IN):: ihalf                            ! IHALF = 0, FULL (radial) MESH); = 1, HALF MESH
        REAL(rprec), DIMENSION(ntheta,nzeta,ns),                        &
          INTENT(IN):: xuv
        REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns),                   &
          INTENT(OUT):: xmn
        REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns),                   &
          INTENT(IN) :: lmns
        INTEGER, INTENT(IN):: iparity                          ! IPARITY = 0, cosine (EVEN); = 1; sine (ODD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: m, n, jk, lk, in, js, istat
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: work1, work2,     &
          work3, work4, lambda, dupdu, cos1l, sin1l, cosml, sinml
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: cosmuip, sinmuip
        REAL(rprec) :: ton, toff, skston, skstoff
!-----------------------------------------------
        CALL second0(ton)
        skston=ton

!
!       1. COMPUTE LAMBDA
!
        ALLOCATE (lambda(ntheta,nzeta,ns), dupdu(ntheta,nzeta,ns),      &  
                  cos1l(ntheta,nzeta,ns), cosml(ntheta,nzeta,ns),       &
                  sin1l(ntheta,nzeta,ns), sinml(ntheta,nzeta,ns),       &
                  stat=istat)
        IF (istat .NE. 0) STOP 'Allocation error #1 in TOMNSP_PEST'
        
        CALL TOIJSP(lmns, lambda, 0, 0, 1, 0)
        CALL TOIJSP(lmns, dupdu,  1, 0, 1, 0)
!       transform from up to u
        dupdu = 1 + dupdu                                              

!       INTERPOLATE TO FULL GRID IF ihalf=0
        IF (ihalf .EQ. 0) THEN
           lambda(:,:,1) = lambda(:,:,2)    !could do better, why bother?
           dupdu (:,:,1) = dupdu (:,:,2)
           DO js=2,ns-1
              lambda(:,:,js) = .5_dp*(lambda(:,:,js)+lambda(:,:,js+1))
              dupdu (:,:,js) = .5_dp*(dupdu (:,:,js)+dupdu (:,:,js+1))
           END DO
           lambda(:,:,ns) = 2*lambda(:,:,ns)-lambda(:,:,ns-1)           
           dupdu (:,:,ns) = 2*dupdu (:,:,ns)-dupdu (:,:,ns-1)           
        END IF

        cos1l = COS(lambda); sin1l = SIN(lambda)
!
!       2. DO Fourier Transform in PEST UP=u+lambda, VP=v Coordinates
!          Use recurrence relations below for EFFICIENCY
!
!          cosmuiP==cos(mu_p) = cosmui(u,m)*cos(m*lambda) - sinmui(u,m)*sin(m*lambda)
!          sinmuiP==sin(mu_p) = sinmui(u,m)*cos(m*lambda) + cosmui(u,m)*sin(m*lambda)
!          Let cos(m*lambda) == cosml(s,u,v,m)
!              sin(m*lambda) == sinml(s,u,v,m)
!          use recurrence formulae to compute these efficiently
!              cosml(m=0) = 1;  sinml(m=0) = 0
!              compute cosml(m=1), sinml(m=1)
!          then
!              cosml(m+1) = cosml(m)*cosml(1) - sinml(m)*sinml(1)
!              sinml(m+1) = sinml(m)*cosml(1) + cosml(m)*sinml(1)
!
        ALLOCATE(work1(0:mpol,nzeta,ns), work2(0:mpol,nzeta,ns),        &
                 work3(ntheta,nzeta,ns), cosmuiP(nzeta,ns),             &
                 sinmuiP(nzeta,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation error #2 in TOMNSP_PEST'
        work1 = zero; work2 = zero

        in = 1+ihalf

        DO m = 0, mpol
           IF (m .EQ. 0) THEN
              cosml = dupdu;  sinml=0
           ELSE IF (m .EQ. 1) THEN
              cosml = cos1l*dupdu; sinml = sin1l*dupdu
           ELSE
              work3 = cosml
              cosml = cosml*cos1l - sinml*sin1l
              sinml = sinml*cos1l + work3*sin1l
           END IF
           DO jk = 1, ntheta                                             ! First, add over poloidal collocation angles
              cosmuiP(:,in:ns) = cosmui(jk,m)*cosml(jk,:,in:ns)         &
                      - sinmui(jk,m)*sinml(jk,:,in:ns)
              sinmuiP(:,in:ns) = sinmui(jk,m)*cosml(jk,:,in:ns)         &
                      + cosmui(jk,m)*sinml(jk,:,in:ns)
              work1(m,:,in:ns) = work1(m,:,in:ns)                       &
                               + xuv(jk,:,in:ns)*cosmuiP(:,in:ns)        ! Note use of cosmuI/sinmuI; include norm. factors
              work2(m,:,in:ns) = work2(m,:,in:ns)                       &
                               + xuv(jk,:,in:ns)*sinmuiP(:,in:ns)
           ENDDO
        END DO

        DEALLOCATE (work3, cosmuip, sinmuip)
        xmn = zero                                                       ! Make sure that FOURIER variable is zeroed
        ALLOCATE(work3(0:mpol,nzeta,ns), work4(0:mpol,nzeta, ns),       &
                 stat=istat)
        IF (istat .NE. 0) STOP 'Allocation error #3 in TOMNSP_PEST'
        work3 = zero; work4 = zero

        DO n = 0, ntor  
          DO lk = 1, nzeta                                               ! Then, add over toroidal collocation angles
            IF (iparity == 0) THEN                                       ! COSINE series              
              work3(:,lk,in:ns) = work1(:,lk,in:ns)*cosnv(lk, n)
              work4(:,lk,in:ns) = work2(:,lk,in:ns)*sinnv(lk, n)
              xmn(:,n,in:ns) =  xmn(:,n,in:ns)                          &
                 + work3(:,lk,in:ns) - work4(:,lk,in:ns)
              IF (n .ne. 0) THEN
                xmn(1:,-n,in:ns) = xmn(1:,-n,in:ns)                     &
                 + work3(1:,lk,in:ns) + work4(1:,lk,in:ns)
              ENDIF
            ELSE                                                ! SINE series
              work3(:,lk,in:ns) = work1(:,lk,in:ns)*sinnv(lk, n)
              work4(:,lk,in:ns) = work2(:,lk,in:ns)*cosnv(lk, n)
              xmn(:,n,in:ns) =  xmn(:,n,in:ns)                          &
                + work3(:,lk,in:ns) + work4(:,lk,in:ns)
              IF (n .ne. 0) THEN
                xmn(1:,-n,in:ns) = xmn(1:,-n,in:ns)                     &
                - work3(1:,lk,in:ns) + work4(1:,lk,in:ns)
              ENDIF
            ENDIF
          ENDDO
        ENDDO 
 
        xmn(0,-ntor:-1,:) = zero                               ! Redundant: To avoid counting (0,-n) for non-zero n.

        DO js = in,ns
           xmn(:,:,js) = orthonorm(0:mpol,-ntor:ntor)*xmn(:,:,js)
        END DO
         
        DEALLOCATE(work1, work2, work3, work4, lambda, cos1l, sin1l,    &
                   sinml, cosml, dupdu)

        CALL second0(toff)
#if defined(SKS)
        skstoff=toff
        IF(DIAGONALDONE.AND.INHESSIAN) tomnsp_time=tomnsp_time+(skstoff-skston)
#endif
        time_tomnsp = time_tomnsp + (toff-ton)

        END SUBROUTINE TOMNSP_PEST


        SUBROUTINE TOMNSP2(XUV, XMN, IPARITY, JS)
!
!       Description: This subroutine moves a quantity X to Fourier space producing its harmonics XMN by
!         summing over all the prescribed (toroidal and poloidal) collocation points:
!           theta_j = j*pi/M, (j=0,...,M);   zeta_k = k*2*pi/(2N+1), k = 0,...., 2N
!         where M = mpol - 1 and N = ntor.
!       IPARITY = 0, means the quantity X (NOT any of its derivatives) is even; = 1, X is odd
!       JS:  Surface where xmn is computed
!
        USE stel_kinds
        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        REAL(rprec), DIMENSION(ntheta,nzeta,ns),                        &
          INTENT(IN):: xuv
        REAL(rprec), DIMENSION(0:mpol,-ntor:ntor),                      &
          INTENT(OUT):: xmn
        INTEGER, INTENT(IN):: iparity, js
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: m, n, jk, lk, in
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: work1, work2,       &
          work3, work4
        REAL(rprec) :: ton, toff
!-----------------------------------------------
        CALL second0(ton)

        ALLOCATE(work1(0:mpol,nzeta), work2(0:mpol,nzeta))
        work1 = zero; work2 = zero

        DO m = 0, mpol
          DO jk = 1, ntheta                                               ! First, add over poloidal collocation angles
            work1(m,:) = work1(m,:) + xuv(jk,:,js)*cosmui(jk, m)          ! Note use of cosmuI/sinmuI; include norm. factors
            work2(m,:) = work2(m,:) + xuv(jk,:,js)*sinmui(jk, m)
          ENDDO
        ENDDO

        xmn = zero                                                        ! Make sure that FOURIER variable is zeroed
        ALLOCATE(work3(0:mpol,nzeta), work4(0:mpol,nzeta))
        work3 =zero; work4 = zero

        DO n = 0, ntor  
          DO lk = 1, nzeta                                      ! Then, add over toroidal collocation angles
            IF (iparity == 0) THEN                              ! COSINE series              
              work3(:,lk) = work1(:,lk)*cosnv(lk, n)
              work4(:,lk) = work2(:,lk)*sinnv(lk, n)
              xmn(:,n) =  xmn(:,n) + work3(:,lk) - work4(:,lk)
              IF (n .ne. 0) THEN
                xmn(1:mpol,-n) = xmn(1:mpol,-n)                         &
                 + work3(1:mpol,lk) + work4(1:mpol,lk)
              ENDIF
            ELSE                                                ! SINE series
              work3(:,lk) = work1(:,lk)*sinnv(lk, n)
              work4(:,lk) = work2(:,lk)*cosnv(lk, n)
              xmn(:,n) =  xmn(:,n) + work3(:,lk) + work4(:,lk)
              IF (n .ne. 0) THEN
                xmn(1:mpol,-n) = xmn(1:mpol,-n)                         &
                 - work3(1:mpol,lk) + work4(1:mpol,lk)
              ENDIF
            ENDIF
          ENDDO
        ENDDO 
 

        xmn(0,-ntor:-1) = zero              ! Redundant: To avoid counting (0,-n) for non-zero n.
        xmn(:,:) = orthonorm(0:mpol,-ntor:ntor)*xmn(:,:)
         
        DEALLOCATE(work1, work2, work3, work4)

        CALL second0(toff)
        time_tomnsp = time_tomnsp + (toff-ton)

        END SUBROUTINE TOMNSP2

        
        SUBROUTINE TOIJSP2(XMN, XUV, IPARITY, JS)
        USE stel_kinds
        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        REAL(rprec), DIMENSION(0:mpol,-ntor:ntor),                 &
          INTENT(IN):: xmn
        REAL(rprec), DIMENSION(ntheta,nzeta,ns),                   &
          INTENT(INOUT):: xuv
        INTEGER, INTENT(IN):: iparity, js                               ! IPARITY = 0, cosine (EVEN); = 1, sine (ODD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: m, n, jk, lk
        INTEGER :: in 
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: work0,         &
          work1, work2
        REAL(rprec) :: ton, toff
!-----------------------------------------------
        CALL second0(ton)

        ALLOCATE(work0(0:mpol,-ntor:ntor))                              ! RS (10/03/06): Added WORK0 for n=0/m=0 normalization

        work0 = xmn
        work0(:,:) = orthonorm(0:mpol,-ntor:ntor)*work0(:,:)
             
        ALLOCATE(work1(ntheta,-ntor:ntor), work2(ntheta,-ntor:ntor))
        work1 = zero; work2 = zero 

        DO jk = 1, ntheta
          DO m = 0, mpol                                                ! First DO sum over poloidal modes
              work1(jk,:) = work1(jk,:) + work0(m,:)*cosmu(jk, m)
              work2(jk,:) = work2(jk,:) + work0(m,:)*sinmu(jk, m)
          ENDDO
        ENDDO
        DEALLOCATE(work0)

        xuv(:,:,js) = 0                                                 ! Make sure that SPATIAL variable is zeroed

        DO lk = 1, nzeta
          
          IF (iparity == 0) THEN                                        ! Do first N=0 mode
            xuv(:,lk,js) = work1(:,0) 
          ELSE
            xuv(:,lk,js) = work2(:,0) 
          ENDIF          
          
          DO n = 1, ntor                                          ! Then sum over N>0 and N<0 toroidal modes                                                     
            IF (iparity == 0) THEN                                ! COSINE series
                xuv(:,lk,js) =  xuv(:,lk,js)                            &
                  + (work1(:,n) + work1(:,-n))*cosnv(lk, n)             &
                  - (work2(:,n) - work2(:,-n))*sinnv(lk, n)
            ELSE                                                  ! SINE series
                xuv(:,lk,js) =  xuv(:,lk,js)                            &
                  + (work1(:,n) - work1(:,-n))*sinnv(lk, n)             &
                  + (work2(:,n) + work2(:,-n))*cosnv(lk, n)
            ENDIF
          ENDDO
        ENDDO

        DEALLOCATE(work1, work2)

        CALL second0(toff)
        time_toijsp = time_toijsp + (toff-ton)

        END SUBROUTINE TOIJSP2

      SUBROUTINE INIT_FOURIER
!
!     This subroutine computes the cosine-sine factors that will be needed when moving between
!     Fourier and real space. All normalizations are contained in the poloidal quantities 
!     used in the Fourier to real space transformation: SINMUI, COSMUI
!
!     Fourier representations are assumed to have STELLARATOR symmetry:
!         1. COSINE (iparity=0):      C(u, v) = sum_u sum_v (C_mn COS(mu + n*nfp*v))
!         2. SINE   (iparity=1):      S(u, v) = sum_u sum_v (S_mn SIN(mu + n*nfp*v))
!
!     The number of collocation points have been set initially equal to the number of modes:
!       theta_j = j*pi/M, (j=0,...,M);   zeta_k = k*2*pi/(2N+1), k = 0,...., 2N
!       where M = mpol - 1 and N = ntor.
!
        IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        REAL(rprec), PARAMETER :: nfactor = 1.41421356237310_dp          ! RS: Normalization of m=0/n=0 Fourier modes.
        INTEGER        :: jk, m, n, lk, istat, ntheta1
        REAL(rprec)    :: arg, cosi, sini
!-----------------------------------------------
!       TO CHANGE TO FULL u-GRID
!       1) nu_i = 2*nu_i
!       2) ntheta1 = ntheta
!       3) arg = twopi*(jk-1)/ntheta1
        
        ntheta1 = ntheta-1
!       EXTEND mpol, ntor FOR POST-PROCESSING OF PRESSURE B dot GRAD P
!        mpol = mpol32_i;  ntor = ntor32_i

        ALLOCATE(cosmu(ntheta,0:mpol), cosmum(ntheta,0:mpol),           &
                 cosmui(ntheta,0:mpol), sinmu(ntheta,0:mpol),           &
                 sinmum(ntheta,0:mpol), sinmui(ntheta,0:mpol),          &
                 cosnv(nzeta,0:ntor), cosnvn(nzeta,0:ntor),             &
                 sinnv(nzeta,0:ntor), sinnvn(nzeta,0:ntor), stat=istat)

        IF (istat .ne. 0) STOP 'Allocation error in fixarray'
        
        dnorm_i = one/(ntheta1*nzeta)                              ! ST:Norm for Fourier-Space -> UV-Space

        DO jk = 1, ntheta
          arg = pi*(jk-1)/ntheta1                                  ! Poloidal collocation points*m
          DO m = 0, mpol
            cosi = COS(m*arg)
            sini = SIN(m*arg)

            CALL Round_Trig(cosi, sini)

            cosmu (jk, m)  = cosi                                 ! Used in UV-SPace to Fourier-space
            cosmum(jk, m)  = m*cosi                               ! Used in UV-SPace to Fourier-space
            cosmui(jk, m) = dnorm_i*cosi                          ! Used in Fourier-Space to UV-space
            sinmu (jk, m)  = sini                                 ! Used in UV-SPace to Fourier-space
            sinmum(jk, m)  = m*sini                               ! Used in UV-SPace to Fourier-space
            sinmui(jk, m) = dnorm_i*sini                          ! Used in Fourier-Space to UV-space
            IF (jk == 1 .or. jk == ntheta) THEN
               cosmui(jk, m) = 0.5_dp*cosmui(jk, m)               ! 0, pi poloidal points are divided by 2 when integrating
               sinmui(jk, m) = 0.5_dp*sinmui(jk, m)
            ENDIF
          ENDDO
        ENDDO

!
!FOR NOW, WE ARE IN 1 FIELD PERIOD
!
        DO lk = 1, nzeta
          arg = twopi*(lk-1)/nzeta                                ! Toroidal collocation points*ABS(n*nfp_i)
          DO n = 0, ntor
            cosi = COS(n*arg)                                     ! sign of "n" taken care later through ISIGN
            sini = SIN(n*arg)
            
            CALL Round_Trig(cosi, sini)

            cosnv (lk, n)  = cosi                                 ! All normalization factors goes in poloidal quantities
            cosnvn(lk, n)  = n*nfp*cosi                         ! NFP is number of field periods
            sinnv (lk, n)  = sini
            sinnvn(lk, n)  = n*nfp*sini
          ENDDO
        ENDDO

!COMPUTE orthonorm FACTOR FOR cos(mu+nv), sin(mu+nv) BASIS (NOT cos(mu)cos(nv) basis!)
!BE SURE TO FIX THIS WHEN DOING ASYMMETRIC CASE!

        ALLOCATE (orthonorm(0:mpol,-ntor:ntor), stat=istat)
        IF (istat .ne. 0) STOP 'Allocate ORTHONORM failed in FIXARRAY!'
        orthonorm = nfactor                                       ! SQRT(2)
        orthonorm(0,0) = 1
        IF (ntheta .eq. (mpol+1)) orthonorm(mpol,:) = 1

        CALL Test_Fourier

      END SUBROUTINE INIT_FOURIER

      SUBROUTINE ROUND_TRIG(cosi, sini)
        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        REAL(dp), INTENT(inout) :: cosi, sini
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        REAL(dp) :: eps
!-----------------------------------------------
        eps = 10*EPSILON(eps)

        IF (ABS(cosi-1) .le. eps) THEN
           cosi = 1
           sini = 0
        ELSE IF (ABS(cosi+1) .le. eps) THEN
           cosi = -1
           sini = 0
        ELSE IF (ABS(sini-1) .le. eps) THEN
           cosi = 0
           sini = 1
        ELSE IF (ABS(sini+1) .le. eps) THEN
           cosi = 0
           sini = -1
        END IF

      END SUBROUTINE Round_Trig


      SUBROUTINE DEALLOC_FOURIER
      IMPLICIT NONE
      INTEGER :: istat

      DEALLOCATE(cosmu, cosmum, cosmui, sinmu, sinmum, sinmui,          &
                 cosnv, cosnvn, sinnv, sinnvn, orthonorm, stat=istat) 
      IF (istat .ne. 0) STOP 'Allocation error in dealloc_fixarray'
         
      END SUBROUTINE DEALLOC_FOURIER

      SUBROUTINE TEST_FOURIER
      USE descriptor_mod, ONLY: iam, nprocs
#if defined(SKS)
      USE timer_mod, ONLY: test_fourier_time
#endif      
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER     :: m, n, lk, istat, n1, iblock, ns_save
      REAL(rprec) :: arg, cosi, sini, testcc, testsc,                   &
                        testcs, testss, ton, toff
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: testr, testi
!-----------------------------------------------
      CALL second0(ton)

      ns_save = ns
      ns = 1
!CHECK ORTHONORMALITY OF BASIS
      ALLOCATE (testr(ntheta,nzeta,ns), testi(0:mpol,-ntor:ntor,ns),    &
                stat=istat)
      testr = 0
      testi = 0
      iblock= 0

      NTOR_LOOP: DO n = 0, ntor
      MPOL_LOOP: DO m = 0, mpol
          iblock = iblock+1
          IF (MOD(iblock-1, nprocs) .NE. iam) CYCLE
          DO lk = 1, nzeta
            testr(:,lk,ns) = cosmu(:,m)*cosnv(lk,n)
          END DO
          CALL tomnsp(testr, testi, 0,0)
          CALL toijsp(testi, testr, 0,0,0,0)
          DO lk = 1, nzeta
            testr(:,lk,ns) = ABS(testr(:,lk,ns) - cosmu(:,m)*cosnv(lk,n))
          END DO
          IF (ANY(testr(:,:,ns) .GT. 1.E-12_dp)) THEN
            PRINT *,'ORTHOGONALITY TEST FAILED FOR MODES ',m,n
            STOP
          END IF

!TEST U and V DERIVATIVES (cos(mu+nv))
          testi = 0
          testi(m,n,ns) = 1
          CALL toijsp(testi, testr,1,0,0,0)
          CALL tomnsp(testr, testi, 1,0)
          testi(m,n,ns) = testi(m,n,ns) + m
          IF (ANY(ABS(testi(:,:,ns)) .GT. 1.E-12_dp)) THEN
            PRINT *,'ORTHOGONALITY TEST FAILED FOR MODES ',m,n
            STOP 'ORTHOGONALITY TEST FAILED FOR MODES'
          END IF

          testi = 0
          testi(m,n,ns) = 1
          CALL toijsp(testi, testr,0,1,0,0)
          CALL tomnsp(testr, testi, 1,0)
          testi(m,n,ns) = testi(m,n,ns) + n*nfp
          IF (ANY(ABS(testi(:,:,ns)) .GT. 1.E-12_dp)) THEN
            PRINT *,'ORTHOGONALITY TEST FAILED FOR MODES ',m,n
            STOP
          END IF

          DO lk = 1, nzeta
            testr(:,lk,ns) = sinmu(:,m)*cosnv(lk,n)
          END DO
          CALL tomnsp(testr, testi, 1,0)
          CALL toijsp(testi, testr, 0,0,1,0)
          DO lk = 1, nzeta
            testr(:,lk,ns) = ABS(testr(:,lk,ns) - sinmu(:,m)*cosnv(lk,n))
          END DO
          IF (ANY(testr(:,:,ns) .GT. 1.E-12_dp)) THEN
            PRINT *,'ORTHOGONALITY TEST FAILED FOR MODES ',m,n
            STOP
          END IF


!TEST U and V DERIVATIVES (sin(mu+nv))
          testi = 0
          testi(m,n,ns) = 1
          CALL toijsp(testi, testr,1,0,1,0)
          CALL tomnsp(testr, testi, 0,0)
          testi(m,n,ns) = testi(m,n,ns) - m
          IF (ANY(ABS(testi(:,:,ns)) .gt. 1.E-12_dp)) THEN
            PRINT *,'ORTHOGONALITY TEST FAILED FOR MODES ',m,n
            STOP
          END IF

          testi = 0
          testi(m,n,ns) = 1
          CALL toijsp(testi, testr,0,1,1,0)
          CALL tomnsp(testr, testi, 0,0)
          testi(m,n,ns) = testi(m,n,ns) - n*nfp
          IF (ANY(ABS(testi(:,:,ns)) .GT. 1.E-12_dp)) THEN
            PRINT *,'ORTHOGONALITY TEST FAILED FOR MODES ',m,n
            STOP
          END IF

          DO lk = 1, nzeta
            testr(:,lk,ns) = cosmu(:,m)*sinnv(lk,n)
          END DO
          CALL tomnsp(testr, testi, 1,0)
          CALL toijsp(testi, testr, 0,0,1,0)
          DO lk = 1, nzeta
            testr(:,lk,ns) = ABS(testr(:,lk,ns) - cosmu(:,m)*sinnv(lk,n))
          END DO
          IF (ANY(testr(:,:,ns) .GT. 1.E-12_dp)) THEN
            PRINT *,'ORTHOGONALITY TEST FAILED FOR MODES ',m,n
            STOP
          END IF

!TEST U and V DERIVATIVES (sin(mu-nv))
          n1 = -n
          IF (m .EQ. 0) n1 = n
          testi = 0
          testi(m,n1,ns) = 1
          CALL toijsp(testi, testr,1,0,1,0)
          CALL tomnsp(testr, testi, 0,0)
          testi(m,n1,ns) = testi(m,n1,ns) - m
          IF (ANY(ABS(testi(:,:,ns)) .gt. 1.E-12_dp)) THEN
            PRINT *,'ORTHOGONALITY TEST FAILED FOR MODES ',m,n
            STOP
          END IF

          testi = 0
          testi(m,n1,ns) = 1
          CALL toijsp(testi, testr,0,1,1,0)
          CALL tomnsp(testr, testi, 0,0)
          testi(m,n1,ns) = testi(m,n1,ns) - n1*nfp
          IF (ANY(ABS(testi(:,:,ns)) .GT. 1.E-12_dp)) THEN
            PRINT *,'ORTHOGONALITY TEST FAILED FOR MODES ',m,n
            STOP
          END IF

          DO lk = 1, nzeta
            testr(:,lk,ns) = sinmu(:,m)*sinnv(lk,n)
          END DO
          CALL tomnsp(testr, testi, 0,0)
          CALL toijsp(testi, testr, 0,0,0,0)
          DO lk = 1, nzeta
            testr(:,lk,ns) = ABS(testr(:,lk,ns) - sinmu(:,m)*sinnv(lk,n))
          END DO
          IF (ANY(testr(:,:,ns) .GT. 1.E-12_dp)) THEN
            PRINT *,'ORTHOGONALITY TEST FAILED FOR MODES ',m,n
            STOP
          END IF

!TEST U and V DERIVATIVES (cos(mu-nv))
          testi = 0
          testi(m,n1,ns) = 1
          CALL toijsp(testi, testr,1,0,0,0)
          CALL tomnsp(testr, testi, 1,0)
          testi(m,n1,ns) = testi(m,n1,ns) + m
          IF (ANY(ABS(testi(:,:,ns)) .GT. 1.E-12_dp)) THEN
            PRINT *,'ORTHOGONALITY TEST FAILED FOR MODES ',m,n
            STOP
          END IF

          testi = 0
          testi(m,n1,ns) = 1
          CALL toijsp(testi, testr,0,1,0,0)
          CALL tomnsp(testr, testi, 1,0)
          testi(m,n1,ns) = testi(m,n1,ns) + n1*nfp
          IF (ANY(ABS(testi(:,:,ns)) .GT. 1.E-12_dp)) THEN
             PRINT *,'ORTHOGONALITY TEST FAILED FOR MODES ',m,n
             STOP
          END IF

         END DO MPOL_LOOP
       END DO NTOR_LOOP

       ns = ns_save
       DEALLOCATE (testr, testi)

       CALL second0(toff)
#if defined(SKS)
       test_fourier_time=toff-ton
#endif       
       IF (iam .eq. 0) THEN
          PRINT 100,'ORTHOGONALITY TEST SUCCEEDED - TIME: ', toff-ton,' s'
          WRITE (33, 100) 'ORTHOGONALITY TEST SUCCEEDED - TIME: ', toff-ton,' s'
       END IF

 100   FORMAT(1x,a,1p,e12.2,a)

       END SUBROUTINE TEST_FOURIER


      END MODULE fourier
