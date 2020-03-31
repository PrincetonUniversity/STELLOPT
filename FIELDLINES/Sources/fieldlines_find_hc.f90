!-----------------------------------------------------------------------
!     Function:      fieldlines_find_hc
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/03/2013
!     Description:   This subroutine finds the homoclinic tangles given
!                    an initial guess for the x-points.
!                         fieldlines_find_hc
!                                  |
!                                E04CCA
!                                  |
!                               faxis_nag
!                                  |
!                              follow_single
!                                  |
!                              ODE Solvers
!                                  |
!                                fblin's
!
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_find_hc
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_grid
      USE fieldlines_lines
      USE fieldlines_runtime
      IMPLICIT NONE
!-----------------------------------------------------------------------
!     Local Variables                                                    
!-----------------------------------------------------------------------
      INTEGER     :: i, j, num_hc, dex
      REAL(rprec) :: r0,z0,r1,z1,rerr,zerr,alpha
      REAL(rprec) :: detq,traq,l1r,l2r,l1i,l2i
      REAL(rprec) :: q4(4),q_map(4),evecr(4),eveci(4)
      REAL(rprec) :: vec1(4),vec2(4)
      REAL(rprec) :: ident(4) = (/1,0,0,1/)
      COMPLEX(8)  :: sqrtval
!-----------------------------------------------------------------------
!     External Functions
!          E04CBF               NAG mimizier
!          E04CBK               NAG Dummy function
!-----------------------------------------------------------------------
      EXTERNAL E04FCF, monit_find_axis, faxis_nag, E04FDZ
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      alpha = 0.75
      num_hc = COUNT(r_hc > 0) ! Number of homoclines
      j=1
      r_start = -1
      z_start = -1
      phi_start = 0.0
      phi_end   = 0.0
      dex = 0
      WRITE(6,'(A)')      '===========HOMOCLINES=========='
      WRITE(6,'(A, I3)')  '   # OF HOMOCLINES  = ',num_hc
      WRITE(6,'(A, I4)')  '   LINES PER HOMOC. = ',2*num_hcp+1
      WRITE(6,'(A, F8.5)')'   LENGTH OF LINE   = ',2*delta_hc
      DO i = 1, num_hc
         WRITE(6,'(A)')    '-------------------------------'
         WRITE(6,'(A, I3.3)')         '   WORKING ON HOMOCLINE: ',i
         r0 = r_hc(i)
         z0 = z_hc(i)
         WRITE(6,'(A,F8.5,A,F8.5,A)') '      XPOINT0[R,Z]   = [',r0,',',z0,']'
         j=1
         DO
            r1=r0
            z1=z0
            CALL follow_single(r1,0.0,z1,pi2,q4)
            rerr = r0-r1
            zerr = z0-z1
            IF (sqrt(rerr*rerr+zerr*zerr) < (follow_tol*10)) EXIT
            q_map(1) = q4(1)-1.0
            q_map(2) = q4(2)
            q_map(3) = q4(3)
            q_map(4) = q4(4)-1.0
            r0   = r0 + alpha*(q_map(4) * rerr - q_map(2) * zerr) / (q_map(1)*q_map(4)-q_map(2)*q_map(3))
            z0   = z0 + alpha*(-q_map(3) * rerr + q_map(1) * zerr) / (q_map(1)*q_map(4)-q_map(2)*q_map(3))
            j= j+1
            IF (j > 1000) EXIT
         END DO
         IF (j > 1000) WRITE(6,'(A)') 'WARNING: NITER > 1000'
         r0 = r1
         z0 = z1
         WRITE(6,'(A,F8.5,A,F8.5,A)') '      XPOINTF[R,Z]   = [',r0,',',z0,']'
         CALL FLUSH(6)
         ! Now return eigenvalues of q
         ! From http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
         detq = (q4(1)*q4(4)-q4(2)*q4(3))
         traq = q4(1)+q4(4)
         sqrtval = sqrt(0.25*traq*traq-detq)
         l1r   = 0.5*traq + REAL(sqrtval)
         l1i   =          + AIMAG(sqrtval)
         l2r   = 0.5*traq - REAL(sqrtval)
         l2i   =          - AIMAG(sqrtval)
         WRITE(6,'(A,F8.5,SP,F8.5,A)') '      EIGN. VAL 1    =  ',l1r,l1i,'i '
         WRITE(6,'(A,F8.5,SP,F8.5,A)') '      EIGN. VAL 2    =  ',l2r,l2i,'i '
         ! Now return the eigenvectors
         vec1 = q4-l1r*ident
         vec2 = q4-l2r*ident
         evecr(1) = vec2(2)
         evecr(3) = vec2(4)
         evecr(2) = vec1(1)
         evecr(4) = vec1(3)
         vec1 = q4-l1i*ident
         vec2 = q4-l2i*ident
         eveci(1) = vec2(2)
         eveci(3) = vec2(4)
         eveci(2) = vec1(1)
         eveci(4) = vec1(3)
         WRITE(6,'(A,F8.5,A,F8.5,A)') '      EIGN. VEC 1    = [',evecr(1),',',evecr(3),']'
         WRITE(6,'(A,F8.5,A,F8.5,A)') '      EIGN. VEC 2    = [',evecr(2),',',evecr(4),']'
         WRITE(6,'(A,F8.5,A,F8.5,A)') '      EIGN. VEC 1i   = [',eveci(1),',',eveci(3),']'
         WRITE(6,'(A,F8.5,A,F8.5,A)') '      EIGN. VEC 2i   = [',eveci(2),',',eveci(4),']'
         CALL FLUSH(6)
         DO j = -num_hcp,num_hcp
            dex = dex + 1
            !dex = j + num_hcp + 1 + (i-1)*(2*num_hcp+1)
            r_start(dex) = r0 + evecr(1)*j*delta_hc/num_hcp
            z_start(dex) = z0 + evecr(3)*j*delta_hc/num_hcp
            phi_start(dex) = 0.0
            phi_end(dex)   = pi2*25
         END DO
         nlines = dex
         DO j = -num_hcp,num_hcp
            dex = dex + 1
            !dex = j + num_hcp + 1 + (i-1)*(2*num_hcp+1)+nlines
            r_start(dex) = r0 + evecr(2)*j*delta_hc/num_hcp
            z_start(dex) = z0 + evecr(4)*j*delta_hc/num_hcp
            phi_start(dex) = 0.0
            phi_end(dex)   = -pi2*25
         END DO
         nlines = dex
      END DO
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fieldlines_find_hc
