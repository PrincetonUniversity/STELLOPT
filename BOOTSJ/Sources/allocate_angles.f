

      subroutine allocate_angles
      use parambs
      use vmec0
      implicit none
      integer :: istat
      INTEGER :: ier_arr(16)
      
      ier_arr = 0
      IF (ALLOCATED(bfield)) DEALLOCATE(bfield); 
      ALLOCATE(bfield(nthetah,nzetah), stat = ier_arr(1))
      IF (ALLOCATED(gsqrt_b)) DEALLOCATE(gsqrt_b); 
      ALLOCATE(gsqrt_b(nthetah,nzetah), stat = ier_arr(2))
      IF (ALLOCATED(sinmi)) DEALLOCATE(sinmi); 
      ALLOCATE(sinmi(-mbuse:mbuse,nthetah), stat = ier_arr(3))
      IF (ALLOCATED(cosmi)) DEALLOCATE(cosmi); 
      ALLOCATE(cosmi(-mbuse:mbuse,nthetah), stat = ier_arr(4))
      IF (ALLOCATED(sinnj)) DEALLOCATE(sinnj); 
      ALLOCATE(sinnj(0:nbuse,nzetah), stat = ier_arr(5))
      IF (ALLOCATED(cosnj)) DEALLOCATE(cosnj); 
      ALLOCATE(cosnj(0:nbuse,nzetah), stat = ier_arr(6))
      IF (ALLOCATED(theta)) DEALLOCATE(theta); 
      ALLOCATE(theta(nthetah), stat = ier_arr(7))
      IF (ALLOCATED(zetah)) DEALLOCATE(zetah); 
      ALLOCATE(zetah(nzetah), stat = ier_arr(8))
      IF (ALLOCATED(sinmim)) DEALLOCATE(sinmim); 
      ALLOCATE(sinmim(-mbuse:mbuse,nthetahm), stat = ier_arr(9))
      IF (ALLOCATED(cosmim)) DEALLOCATE(cosmim); 
      ALLOCATE(cosmim(-mbuse:mbuse,nthetahm), stat = ier_arr(10))
      IF (ALLOCATED(sinnjm)) DEALLOCATE(sinnjm); 
      ALLOCATE(sinnjm(0:nbuse,nzetahm), stat = ier_arr(11))
      IF (ALLOCATED(cosnjm)) DEALLOCATE(cosnjm); 
      ALLOCATE(cosnjm(0:nbuse,nzetahm), stat = ier_arr(12))
      IF (ALLOCATED(b2obm)) DEALLOCATE(b2obm); 
      ALLOCATE(b2obm(nthetah,nzetah), stat = ier_arr(13))
      IF (ALLOCATED(bfieldm)) DEALLOCATE(bfieldm); 
      ALLOCATE(bfieldm(nthetahm,nzetahm), stat = ier_arr(14))
      IF (ANY(ier_arr.ne.0)) stop 'allocation error in allocate_angles'
      
      
!      allocate(bfield(nthetah,nzetah), gsqrt_b(nthetah,nzetah),
!     1  sinmi(-mbuse:mbuse,nthetah), sinnj(0:nbuse,nzetah),
!     2  cosmi(-mbuse:mbuse,nthetah), cosnj(0:nbuse,nzetah),
!     3  theta(nthetah), zetah(nzetah),
!     4  sinmim(-mbuse:mbuse,nthetahm), sinnjm(0:nbuse,nzetahm),
!     5  cosmim(-mbuse:mbuse,nthetahm), cosnjm(0:nbuse,nzetahm),
!     6  b2obm(nthetah,nzetah), bfieldm(nthetahm,nzetahm), stat=istat)

!      if (istat .ne. 0) stop 'allocation error in allocate_angles'

      end subroutine allocate_angles
