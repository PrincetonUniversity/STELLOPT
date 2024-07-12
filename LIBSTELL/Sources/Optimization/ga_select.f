      SUBROUTINE ga_select(mate,ipick)
c#######################################################################
c
c  This routine SELECTs the better of two possible parents for mating.
c
      USE ga_mod
      IMPLICIT NONE
      INTEGER :: mate, ipick, IFirst, isecond
      SAVE
c
      IF(ipick+1.gt.npopsiz) CALL ga_shuffle(ipick)
      IFirst=ipick
      isecond=ipick+1
      ipick=ipick+2
      IF(fitness(ifirst).gt.fitness(isecond)) THEN
         mate=ifirst
      ELSE
         mate=isecond
      END IF
c     WRITE(3,*)'select',ifirst,isecond,fitness(ifirst),fitness(isecond)
c
      END SUBROUTINE ga_select
