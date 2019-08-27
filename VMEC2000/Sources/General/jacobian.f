      SUBROUTINE jacobian_par
      USE vmec_input, ONLY: nzeta
      USE vmec_main, ONLY: ohs, nrzt, irst, nznt, iter2
      USE vmec_params, ONLY: meven, modd
      USE realspace
      USE vmec_dim, ONLY: ns, ntheta3
      USE vforces, pr12 => parmn_o, pzu12 => parmn_e, pru12 => pazmn_e, 
     &             prs => pbzmn_e, pzs => pbrmn_e, ptau => pazmn_o
      USE parallel_include_module

      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(dp), PARAMETER :: zero = 0, p5 = 0.5_dp, p25 = p5*p5
      REAL(dp) :: dphids
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, nsmin, nsmax, lnsnum
      REAL(dp) :: ltaumax, ltaumin
      REAL(dp) :: taumax, taumin
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: temp
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: minarr, maxarr
      REAL(dp) :: t1, t2, tjacon, tjacoff
C-----------------------------------------------
      CALL second0(tjacon)

      nsmin=MAX(2,tlglob); nsmax=t1rglob;
      dphids = p25
      irst = 1

      DO i = nsmin, nsmax
         pru12(:,i) = p5*(pru(:,i,meven) + pru(:,i-1,meven) +
     &                    pshalf(:,i)*(pru(:,i,modd) +
     &                                 pru(:,i-1,modd)))
         pzs(:,i)   = ohs*(pz1(:,i,meven) - pz1(:,i-1,meven) +
     &                     pshalf(:,i)*(pz1(:,i,modd) -
     &                                  pz1(:,i-1,modd)))
         ptau(:,i) = pru12(:,i)*pzs(:,i)
     &             + dphids*(pru(:,i,modd)*pz1(:,i,modd) +
     &                       pru(:,i-1,modd)*pz1(:,i-1,modd) +
     &                       (pru(:,i,meven)*pz1(:,i,modd) +
     &                        pru(:,i-1,meven)*pz1(:,i-1,modd)) /
     &                       pshalf(:,i))
      END DO

      DO i = nsmin, nsmax
         pzu12(:,i) = p5*(pzu(:,i,meven) + pzu(:,i-1,meven) +
     &                    pshalf(:,i)*(pzu(:,i,modd) +
     &                                 pzu(:,i-1,modd)))
         prs(:,i)   = ohs*(pr1(:,i,meven) - pr1(:,i-1,meven) +
     &                     pshalf(:,i)*(pr1(:,i,modd) -
     &                                  pr1(:,i-1,modd)))
         pr12(:,i)  = p5*(pr1(:,i,meven) + pr1(:,i-1,meven) +
     &                    pshalf(:,i)*(pr1(:,i,modd) +
     &                                 pr1(:,i-1,modd)))
         ptau(:,i) = ptau(:,i) - prs(:,i)*pzu12(:,i)
     &             - dphids*(pzu(:,i,modd)*pr1(:,i,modd) +
     &                       pzu(:,i-1,modd)*pr1(:,i-1,modd) +
     &                       (pzu(:,i,meven)*pr1(:,i,modd) +
     &                        pzu(:,i-1,meven)*pr1(:,i-1,modd)) /
     &                       pshalf(:,i))
      END DO

      ALLOCATE(temp(1:nznt))
      temp(:)=ptau(:,2)
      ptau(:,1)=temp(:)
      DEALLOCATE(temp)

      ltaumax=MAXVAL(ptau(:,nsmin:nsmax))
      ltaumin=MINVAL(ptau(:,nsmin:nsmax))
!      ltaumax=MAXVAL(ptau(:,tlglob:trglob))
!      ltaumin=MINVAL(ptau(:,tlglob:trglob))

      taumax=ltaumax
      taumin=ltaumin

      IF (nranks.GT.1.AND.grank.LT.nranks) THEN
         CALL second0(t1)
         CALL MPI_Allreduce(ltaumax,taumax,1,MPI_REAL8,
     &                      MPI_MAX,NS_COMM,MPI_ERR)
         CALL MPI_Allreduce(ltaumin,taumin,1,MPI_REAL8,
     &                      MPI_MIN,NS_COMM,MPI_ERR)
         CALL second0(t2)
         allreduce_time = allreduce_time + (t2-t1)
      END IF

      IF (taumax*taumin .lt. zero) THEN
         irst = 2
      END IF

      CALL second0(tjacoff)
      jacobian_time=jacobian_time+(tjacoff-tjacon)

      END SUBROUTINE jacobian_par

      SUBROUTINE jacobian
      USE vmec_main, ONLY: ohs, nrzt, irst, iter2
      USE vmec_params, ONLY: meven, modd
      USE realspace
      USE vmec_dim, ONLY: ns
      USE vforces, r12 => armn_o, ru12 => azmn_e, zu12 => armn_e, 
     &             rs => bzmn_e, zs => brmn_e, tau => azmn_o  !,z12 => blmn_e,

      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(dp), PARAMETER :: zero = 0, p5 = 0.5_dp, p25 = p5*p5
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l
      REAL(dp) :: taumax, taumin, dphids, temp(nrzt/ns)

C-----------------------------------------------
!
!     (RS, ZS)=(R, Z) SUB S, (RU12, ZU12)=(R, Z) SUB THETA(=U)
!     AND TAU=SQRT(G)/R ARE DIFFERENCED ON HALF MESH
!     NOTE: LOOPS WERE SPLIT TO ALLOW EFFICIENT MEMORY BUS USAGE
!
!     SQRT(G) = R*TAU IS COMPUTED IN BCOVAR
!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING MORE THAN ONE POINTER!
!
!
!     HERE, TAU = (Ru * Zs - Rs * Zu). THE DERIVATIVES OF SHALF = SQRT(s)
!     WERE COMPUTED EXPLICITLY AS: d(shalf)/ds = .5/shalf
!
!     NOTE: z12 IS USED IN RECONSTRUCTION PART OF CODE ONLY; COULD BE ELIMINATED...
!
!
      dphids = p25
      irst = 1

CDIR$ IVDEP
      DO l = 2,nrzt
        ru12(l) = p5*(ru(l,meven) + ru(l-1,meven) +
     &                shalf(l)*(ru(l,modd) + ru(l-1,modd)))
        zs(l)   = ohs*(z1(l,meven) - z1(l-1,meven) +
     &                 shalf(l)*(z1(l,modd) - z1(l-1,modd)))
        tau(l)  = ru12(l)*zs(l)
     &          + dphids*(ru(l,modd)*z1(l,modd) +
     &                    ru(l-1,modd)*z1(l-1,modd) +
     &                    (ru(l,meven)*z1(l,modd) +
     &                     ru(l-1,meven)*z1(l-1,modd))/shalf(l))
      ENDDO


CDIR$ IVDEP
      DO l = 2,nrzt
         zu12(l) = p5*(zu(l,meven) + zu(l-1,meven)
     &           + shalf(l)*(zu(l,modd)  + zu(l-1,modd)))
         rs(l)   = ohs*(r1(l,meven) - r1(l-1,meven)
     &           + shalf(l)*(r1(l,modd)  - r1(l-1,modd)))
         r12(l)  = p5*(r1(l,meven) + r1(l-1,meven) +
     &             shalf(l)*(r1(l,modd)  + r1(l-1,modd)))
         tau(l)  = tau(l) - rs(l)*zu12(l) -
     &             dphids*(zu(l,modd)*r1(l,modd) +
     &                     zu(l-1,modd)*r1(l-1,modd) +
     &                     (zu(l,meven)*r1(l,modd) +
     &                      zu(l-1,meven)*r1(l-1,modd))/shalf(l))
      END DO

!
!     TEST FOR SIGN CHANGE IN JACOBIAN
!
      temp(:) = tau(2:nrzt:ns)
      tau(1:nrzt:ns) = temp(:)
      taumax = MAXVAL(tau(2:nrzt))
      taumin = MINVAL(tau(2:nrzt))
      IF (taumax*taumin .lt. zero) THEN
         irst = 2
      END IF

      END SUBROUTINE jacobian
