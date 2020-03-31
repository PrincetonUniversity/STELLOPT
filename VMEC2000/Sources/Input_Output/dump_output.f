      MODULE dump_output
      USE vmec_main
      USE realspace
      IMPLICIT NONE
      
      INTEGER, PRIVATE :: js, lu, lv

      CONTAINS

      SUBROUTINE dump_special

!     Dumps out any "special" information the user might want for debugging, etc

      js = (ns-1)/4 + 1

      WRITE (66,100) ns, ntheta, ntheta3, nzeta
 100  FORMAT("NS: ",i4," NU: ",i4," NU2: ",i4," NV: ",i4)
      WRITE (66, 110) js
 110  FORMAT("JS POINT: ", i4, /)

      WRITE (66, *)"          R1          Z1          RU          ZU",
     1             "          RV          ZV"
      WRITE (66, *)
      CALL WRITE_RZL(r1,z1,ru,zu,rv,zv,js)     

      END SUBROUTINE dump_special

      SUBROUTINE WRITE_RZL(r1,z1,ru,zu,rv,zv,jspt)
      REAL(rprec), INTENT(in), DIMENSION(ns,nzeta,ntheta3,0:1) ::
     1    r1, z1, ru, zu, rv, zv
      REAL(rprec)    :: factor, temp1, temp2, temp3, temp4, temp5, temp6
      INTEGER :: jspt

      factor = sqrts(jspt)

      DO lu = 1, ntheta3
         WRITE (66, 100) lu
         DO lv = 1, nzeta
            temp1 = r1(jspt,lv,lu,0) + factor*r1(jspt,lv,lu,1)
            temp2 = z1(jspt,lv,lu,0) + factor*z1(jspt,lv,lu,1)
            temp3 = ru(jspt,lv,lu,0) + factor*ru(jspt,lv,lu,1)
            temp4 = zu(jspt,lv,lu,0) + factor*zu(jspt,lv,lu,1)
            temp5 = rv(jspt,lv,lu,0) + factor*rv(jspt,lv,lu,1)
            temp6 = zv(jspt,lv,lu,0) + factor*zv(jspt,lv,lu,1)
            WRITE (66, 200) lv, temp1, temp2, temp3, 
     1                          temp4, temp5, temp6
         END DO
      END DO

 100  FORMAT ("lu: ", i4)
 200  FORMAT (i4, 1p,6e12.4)

      END SUBROUTINE WRITE_RZL

      END MODULE dump_output
